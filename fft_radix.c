#include <stdio.h>
#include <math.h>

#include "fft_radix.h"

#define M_PI       3.14159265358979323846264338328      /* pi */

#define REAL(a,stride,i) ((a)[2*(stride)*(i)])
#define IMAG(a,stride,i) ((a)[2*(stride)*(i)+1])

static int fft_binary_logn(const int n){
    int ntest;
    int binary_logn = 0;
    int k = 1;

    while (k < n){
        k *= 2;
        binary_logn++;
    }

    ntest = (1 << binary_logn);

    if (n != ntest){
        return -1; /* n is not a power of 2 */
    }

    return binary_logn;
}

static int fft_radix2_bitreverse_order(double *data, const int stride, const int n, int logn){
    /* This is the Goldrader bit-reversal algorithm */
    int i;
    int j = 0;

    logn = 0; /* not needed for this algorithm */

    for (i = 0; i < n - 1; i++){
        int k = n / 2;

        if (i < j){
            const double tmp_real = REAL(data, stride, i);
            const double tmp_imag = IMAG(data, stride, i);
            REAL(data, stride, i) = REAL(data, stride, j);
            IMAG(data, stride, i) = IMAG(data, stride, j);
            REAL(data, stride, j) = tmp_real;
            IMAG(data, stride, j) = tmp_imag;
        }

        while (k <= j){
            j = j - k;
            k = k / 2;
        }

        j += k;
    }

    return 0;
}

int fft_radix2_transform(double *data,int stride, int n, fft_radix_direction sign){
    int result ;
    int dual;
    int bit; 
    int logn = 0;
    int status;

    if (n == 1) /* identity operation */{
        return 0 ;
    }

    /* make sure that n is a power of 2 */
    result = fft_binary_logn(n) ;

    if (result == -1){
        printf("n is not a power of 2\n");
    }
    else{
        logn = result;
    }

    /* apply fft recursion */
    dual = n / 2;

    for (bit = 0; bit < logn; bit++){
        double w_real = 1.0;
        double w_imag = 0.0;

        const double theta = 2.0 * ((int)sign) * M_PI / ((double)(2 * dual));

        const double s = sin(theta);
        const double t = sin(theta / 2.0);
        const double s2 = 2.0 * t * t;

        int a, b;

        for (b = 0; b < dual; b++){
            for (a = 0; a < n; a += 2 * dual){
                const int i = b + a;
                const int j = b + a + dual;

                const double t1_real = REAL(data, stride, i) + REAL(data, stride, j);
                const double t1_imag = IMAG(data, stride, i) + IMAG(data, stride, j);
                const double t2_real = REAL(data, stride, i) - REAL(data, stride, j);
                const double t2_imag = IMAG(data, stride, i) - IMAG(data, stride, j);

                REAL(data, stride, i) = t1_real;
                IMAG(data, stride, i) = t1_imag;
                REAL(data, stride, j) = w_real * t2_real - w_imag * t2_imag;
                IMAG(data, stride, j) = w_real * t2_imag + w_imag * t2_real;
            }

            /* trignometric recurrence for w-> exp(i theta) w */
            const double tmp_real = w_real - s * w_imag - s2 * w_real;
            const double tmp_imag = w_imag + s * w_real - s2 * w_imag;
            w_real = tmp_real;
            w_imag = tmp_imag;
        }
        dual /= 2;
    }

    /* bit reverse the ordering of output data for decimation in frequency algorithm */
    status = fft_radix2_bitreverse_order(data, stride, n, logn);

    return 0;
}

