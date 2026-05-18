# /// script
# requires-python = ">=3.9"
# dependencies = ["sounddevice", "scipy", "numpy"]
# ///

"""
Real-time MVDR beamforming demo.

Device layout (48 kHz, 7 channels):
  ch 0-3  : raw microphone signals
  ch 4    : cascade signal
  ch 5    : loopback / echo reference
  ch 6    : DSP-chip processed output

Outputs (all at 16 kHz):
  out_ch0_raw.wav  – raw channel 0 (downsampled)
  out_mvdr.wav     – MVDR beamformer output
  out_dsp.wav      – DSP chip output channel (ch 6, downsampled)

Run:
  uv run run_mvdr.py
"""

import ctypes
import os
import queue
import sys
import threading
import wave

import numpy as np
import sounddevice as sd
import scipy.signal as sig

# ── DLL ───────────────────────────────────────────────────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))
dll_path = os.path.join(HERE, "btk.dll")
if not os.path.exists(dll_path):
    sys.exit(f"ERROR: {dll_path} not found. Run compile_msvc.bat first.")

btk = ctypes.CDLL(dll_path)
btk.btk_beamforming_init.restype       = None
btk.btk_beamforming_init.argtypes      = []
btk.btk_beamforming_process.restype    = ctypes.c_int
btk.btk_beamforming_process.argtypes   = [
    ctypes.POINTER(ctypes.c_short),  # input_channels[4][64]  (row-major)
    ctypes.POINTER(ctypes.c_short),  # out[64]
]
btk.btk_beamforming_set_location.restype  = None
btk.btk_beamforming_set_location.argtypes = [ctypes.c_short, ctypes.c_short]

# btk_beamforming_init() also calls set_location(0, 45) internally
btk.btk_beamforming_init()
btk.btk_beamforming_set_location(90, 30)  # 左侧 30° 仰角

# ── Constants ─────────────────────────────────────────────────────────────────
DEV_SR   = 48_000   # device sample rate
ALGO_SR  = 16_000   # algorithm sample rate (hard-coded in the C)
RATIO    = DEV_SR // ALGO_SR   # 3
D        = 64                  # algorithm frame size (samples @ 16 kHz)
DEV_BLOK = D * RATIO           # 192 samples per block @ 48 kHz

N_MIC    = 4   # channels 0-3
DSP_CH   = 6   # change to 4 or 5 if you want cascade / loopback instead
N_DEV_CH = 7

# ── Anti-alias LP filter for 3× decimation (cutoff 7.5 kHz) ──────────────────
aa_sos = sig.butter(8, 7500, btype="low", fs=DEV_SR, output="sos")
zi_mic = [sig.sosfilt_zi(aa_sos) for _ in range(N_MIC)]
zi_dsp = sig.sosfilt_zi(aa_sos)

# ── Output buffers ────────────────────────────────────────────────────────────
buf_ch0  = []   # raw ch0  @ 16 kHz
buf_mvdr = []   # MVDR     @ 16 kHz
buf_dsp  = []   # DSP ref  @ 16 kHz

# ── Threading ─────────────────────────────────────────────────────────────────
audio_q    = queue.Queue()
stop_event = threading.Event()
remainder  = np.empty((0, N_DEV_CH), dtype=np.float32)


def audio_callback(indata, frames, time_info, status):
    """Sounddevice callback – runs in a high-priority audio thread."""
    if status:
        print(f"[stream] {status}", file=sys.stderr)
    audio_q.put(indata.copy())   # indata is float32 in [-1, 1]


def processing_thread():
    global remainder, zi_mic, zi_dsp

    while not (stop_event.is_set() and audio_q.empty()):
        try:
            chunk = audio_q.get(timeout=0.3)   # shape: (frames, N_DEV_CH)
        except queue.Empty:
            continue

        if len(remainder):
            remainder = np.vstack([remainder, chunk])
        else:
            remainder = chunk

        while len(remainder) >= DEV_BLOK:
            blk = remainder[:DEV_BLOK].astype(np.float64)
            remainder = remainder[DEV_BLOK:]

            # Low-pass + decimate 4 mic channels: 48 kHz → 16 kHz
            mic_16k = np.empty((N_MIC, D), dtype=np.float64)
            for c in range(N_MIC):
                filt, zi_mic[c] = sig.sosfilt(aa_sos, blk[:, c], zi=zi_mic[c])
                mic_16k[c] = filt[::RATIO]

            # Low-pass + decimate DSP reference channel
            dsp_filt, zi_dsp = sig.sosfilt(aa_sos, blk[:, DSP_CH], zi=zi_dsp)
            dsp_16k = dsp_filt[::RATIO]

            # Scale float → int16 for the C library
            scale = 32767.0
            mic_i16 = np.ascontiguousarray(
                (mic_16k * scale).clip(-32768, 32767), dtype=np.int16
            )

            # MVDR – pass raw numpy data pointers (row-major [4][64])
            out_i16 = np.zeros(D, dtype=np.int16)
            btk.btk_beamforming_process(
                mic_i16.ctypes.data_as(ctypes.POINTER(ctypes.c_short)),
                out_i16.ctypes.data_as(ctypes.POINTER(ctypes.c_short)),
            )

            buf_ch0.extend(mic_i16[0].tolist())
            buf_mvdr.extend(out_i16.tolist())
            buf_dsp.extend(
                (dsp_16k * scale).clip(-32768, 32767).astype(np.int16).tolist()
            )


# ── Device selection ──────────────────────────────────────────────────────────
print("Input devices with >= 7 channels:")
for i, d in enumerate(sd.query_devices()):
    if d["max_input_channels"] >= N_DEV_CH:
        print(f"  [{i:2d}] {d['name']}")

raw_idx = input("\nDevice index (Enter = system default): ").strip()
dev_idx = int(raw_idx) if raw_idx else None

# ── Start ─────────────────────────────────────────────────────────────────────
t = threading.Thread(target=processing_thread, daemon=True)
t.start()

print("\nRecording… press Enter to stop.\n")
with sd.InputStream(
    device=dev_idx,
    samplerate=DEV_SR,
    channels=N_DEV_CH,
    dtype="float32",
    blocksize=DEV_BLOK,
    callback=audio_callback,
):
    input()

stop_event.set()
t.join(timeout=5)


# ── Save WAV ──────────────────────────────────────────────────────────────────
def save_wav(path, samples, sr):
    arr = np.array(samples, dtype=np.int16)
    with wave.open(path, "w") as f:
        f.setnchannels(1)
        f.setsampwidth(2)
        f.setframerate(sr)
        f.writeframes(arr.tobytes())
    dur = len(arr) / sr
    print(f"  {path}  —  {dur:.1f}s  ({len(arr)} samples)")


out_dir = HERE
print("\nSaving…")
save_wav(os.path.join(out_dir, "out_ch0_raw.wav"), buf_ch0,  ALGO_SR)
save_wav(os.path.join(out_dir, "out_mvdr.wav"),    buf_mvdr, ALGO_SR)
save_wav(os.path.join(out_dir, "out_dsp.wav"),     buf_dsp,  ALGO_SR)
print("Done.")
