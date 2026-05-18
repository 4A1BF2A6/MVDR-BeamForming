#ifndef __FFT_RADIX_H__
#define __FFT_RADIX_H__

typedef enum{
    fft_radix_forward  = -1,
    fft_radix_backward = +1
} fft_radix_direction;

int fft_radix2_transform(double *data,int stride, int n, fft_radix_direction sign);

#endif