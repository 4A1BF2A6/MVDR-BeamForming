@echo off
setlocal

set VCVARS="C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"

if not exist %VCVARS% (
    echo ERROR: vcvarsall.bat not found.
    exit /b 1
)

call %VCVARS% x64
if %errorlevel% neq 0 (
    echo ERROR: vcvarsall.bat failed.
    exit /b 1
)

pushd "%~dp0"

cl /nologo /O2 /LD /TC ^
   beamformerMVDR.c fft_radix.c pseudoinverse.c ^
   /Febtk.dll ^
   /link /DEF:btk.def

if %errorlevel% equ 0 (
    echo.
    echo SUCCESS: btk.dll created.
    del /Q btk.exp btk.lib 2>nul
    del /Q beamformerMVDR.obj fft_radix.obj pseudoinverse.obj 2>nul
) else (
    echo.
    echo FAILED.
)

popd
