#include <stdio.h>
#include <stdbool.h>
#include <math.h>

int i4_min(int i1, int i2){
    if (i1 < i2){
        return i1;
    }
    else{
        return i2;
    }
}

int i4_max(int i1, int i2){
    if (i2 < i1){
        return i1;
    }
    else{
        return i2;
    }
}

float r4_max ( float x, float y ){
    if (y < x){
        return x;
    }
    else{
        return y;
    }
}

float r4_abs(float x){
    float value;

    if (0.0 <= x){
        value = x;
    }
    else{
        value = -x;
    }
    return value;
}

void p_complex_mul(float a_Re,float a_Lm, float b_Re, float b_Lm, float *c_Re, float *c_Lm){
    *c_Re = a_Re * b_Re - a_Lm * b_Lm;
    *c_Lm = a_Re * b_Lm + a_Lm * b_Re;
}

float p_complex_abs(float z_Re, float z_Lm){
	return hypot(z_Re, z_Lm);
}

void p_complex_div(float a_Re,float a_Lm, float b_Re, float b_Lm, float *c_Re, float *c_Lm){
    double s = 1.0 / p_complex_abs(b_Re,b_Lm);

	double sbr = s * b_Re;
	double sbi = s * b_Lm;

	*c_Re = (a_Re * sbr + a_Lm * sbi) * s;
	*c_Lm = (a_Lm * sbr - a_Re * sbi) * s;
}

//****************************************************************************80
//
//  Purpose:
//
//    CSWAP interchanges two vectors.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    Jack Dongarra
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, complex <float> CX[], one of the vectors to swap.
//
//    Input, int INCX, the increment between successive entries of CX.
//
//    Input/output, complex <float> CY[], one of the vectors to swap.
//
//    Input, int INCY, the increment between successive elements of CY.
//
void cswap ( int n, float *cx_Re, float *cx_Lm, int incx, float* cy_Re, float* cy_Lm, int incy ){

    float ctemp_Re,ctemp_Lm;
    int i;
    int ix;
    int iy;

    if (n <= 0){
        return;
    }

    if (incx == 1 && incy == 1){
        for (i = 0; i < n; i++){
            ctemp_Re = cx_Re[i];
            ctemp_Lm = cx_Lm[i];

            cx_Re[i] = cy_Re[i];
            cx_Lm[i] = cy_Lm[i];

            cy_Re[i] = ctemp_Re;
            cy_Lm[i] = ctemp_Lm;
        }
    }
    else{
        if (0 <= incx){
            ix = 0;
        }
        else{
            ix = (-n + 1) * incx;
        }

        if (0 <= incy){
            iy = 0;
        }
        else{
            iy = (-n + 1) * incy;
        }

        for (i = 0; i < n; i++){
            ctemp_Re = cx_Re[ix];
            ctemp_Lm = cx_Lm[ix];

            cx_Re[ix] = cy_Re[iy];
            cx_Lm[ix] = cy_Lm[iy];

            cy_Re[iy] = ctemp_Re;
            cy_Lm[iy] = ctemp_Lm;

            ix = ix + incx;
            iy = iy + incy;
        }
    }

    return;
}

//****************************************************************************80
//
//  Purpose:
//
//    SCNRM2 returns the euclidean norm of a vector.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    SCNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
//            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, complex <float> X[], the vector.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, float SCNRM2, the norm of the vector.
//

float scnrm2(int n, float *x_Re, float *x_Lm, int incx){
    int i;
    int ix;
    float scale;
    float ssq;
    float temp;
    float value;

    if (n < 1 || incx < 1){
        value = 0.0;
    }
    else{
        scale = 0.0;
        ssq = 1.0;
        ix = 0;

        for (i = 0; i < n; i++){
            if (x_Re[ix] != 0.0){
                temp = r4_abs(x_Re[ix]);
                if (scale < temp){
                    ssq = 1.0 + ssq * pow(scale / temp, 2);
                    scale = temp;
                }
                else{
                    ssq = ssq + pow(temp / scale, 2);
                }
            }

            if (x_Lm[ix] != 0.0){
                temp = r4_abs(x_Lm[ix]);
                if (scale < temp){
                    ssq = 1.0 + ssq * pow(scale / temp, 2);
                    scale = temp;
                }
                else{
                    ssq = ssq + pow(temp / scale, 2);
                }
            }
            ix = ix + incx;
        }
        value = scale * sqrt(ssq);
        //printf("%f,%f ",x_Re[ix],x_Lm[ix]);
    }
    return value;
}

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIGN returns the sign of an R4.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose sign is desired.
//
//    Output, float R4_SIGN, the sign of X.
//
float r4_sign ( float x ){
    float value;

    if (x < 0.0){
        value = -1.0;
    }
    else{
        value = 1.0;
    }
    return value;
}

//****************************************************************************80
//
//  Purpose:
//
//    CABS1 returns the L1 norm of a number.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L1 norm of a complex number is the sum of the absolute values
//    of the real and imaginary components.
//
//    CABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )
//
//  Modified:
//
//    31 March 2007
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <float> Z, the number whose norm is desired.
//
//    Output, float CABS1, the L1 norm of Z.
//
float cabs1(float z_Re,float z_Lm){
    float value;
    value = r4_abs(z_Re) + r4_abs(z_Lm);
    return value;
}

//****************************************************************************80
//
//  Purpose:
//
//    CABS2 returns the L2 norm of a number.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L2 norm of a complex number is the square root of the sum 
//    of the squares of the real and imaginary components.
//
//    CABS2 ( Z ) = sqrt ( ( real ( Z ) )**2 + ( imaginary ( Z ) )**2 )
//
//  Modified:
//
//    10 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <float> Z, the number whose norm is desired.
//
//    Output, float CABS2, the L2 norm of Z.
//
float cabs2(float z_Re, float z_Lm){
    float value;
    value = sqrt(pow(z_Re, 2) + pow(z_Lm, 2));
    return value;
}

//****************************************************************************80
//
//  Purpose:
//
//    CSIGN2 is a transfer-of-sign function.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L2 norm is used.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> Z1, Z2, the arguments.
//
//    Output, complex <float> CSIGN2,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
void csign2(float z1_Re, float z1_Lm, float z2_Re, float z2_Lm, float *value_Re, float *value_Lm){
    if(cabs2(z2_Re,z2_Lm) == 0.0){
        *value_Re = 0.0;
        *value_Lm = 0.0;
    }
    else{
        *value_Re = cabs2(z1_Re,z1_Lm) * (z2_Re / cabs2(z2_Re,z2_Lm));
        *value_Lm = cabs2(z1_Re,z1_Lm) * (z2_Lm / cabs2(z2_Re,z2_Lm));
    }
    //printf("%g,%g,%g,%g,%g,%g\n",z1_Re,z1_Lm,z2_Re,z2_Lm,*value_Re,*value_Lm);
}

//****************************************************************************80
//
//  Purpose:
//
//    CSCAL scales a vector by a constant.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    Jack Dongarra
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, complex <float> CA, the multiplier.
//
//    Input/output, complex <float> CX[], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of CX.
//
void cscal(int n, float ca_Re, float ca_Lm, float *cx_Re, float *cx_Lm, int incx){
    int i;

    if (n <= 0 || incx <= 0){
        return;
    }

    if (incx == 1){
        for (i = 0; i < n; i++){
            p_complex_mul(ca_Re,ca_Lm,cx_Re[i],cx_Lm[i],&cx_Re[i],&cx_Lm[i]); 
            //printf("%g,%g ",ca_Re,ca_Lm);
        }
    }
    else{
        for (i = 0; i < n; i++){
            p_complex_mul(ca_Re,ca_Lm,cx_Re[i*incx],cx_Lm[i*incx],&cx_Re[i*incx],&cx_Lm[i*incx]);
        }
    }
    //printf("\n");
    return;
}

//****************************************************************************80
//
//  Purpose:
//
//    CDOTC forms the conjugated dot product of two vectors.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//  Modified:
//
//    10 April 2006
//
//  Author:
//
//    Jack Dongarra
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, complex <float> CX[], the first vector.
//
//    Input, int INCX, the increment between successive entries in CX.
//
//    Input, complex <float> CY[], the second vector.
//
//    Input, int INCY, the increment between successive entries in CY.
//
//    Output, complex <float> CDOTC, the conjugated dot product of
//    the corresponding entries of CX and CY.
//
void cdotc(int n, float *cx_Re, float *cx_Lm, int incx, 
          float *cy_Re, float *cy_Lm, int incy, float *value_Re, float *value_Lm){
    int i;
    int ix;
    int iy;

    *value_Re = 0.0;
    *value_Lm = 0.0;

    float tmp_Re,tmp_Lm;

    if (n <= 0){
        return;
    }

    if (incx == 1 && incy == 1){
        for (i = 0; i < n; i++){            
            p_complex_mul(cx_Re[i],-cx_Lm[i],cy_Re[i],cy_Lm[i],&tmp_Re,&tmp_Lm); //conj cx[i]
            *value_Re = *value_Re + tmp_Re;
            *value_Lm = *value_Lm + tmp_Lm;
            //printf("%f,%f,%f,%f\n",tmp_Re,tmp_Lm,*value_Re,*value_Lm);
        }
    }
    else{
        if (0 <= incx){
            ix = 0;
        }
        else{
            ix = (-n + 1) * incx;
        }

        if (0 <= incy){
            iy = 0;
        }
        else{
            iy = (-n + 1) * incy;
        }

        for(i = 0; i < n; i++){
            p_complex_mul(cx_Re[ix],-cx_Lm[ix],cy_Re[iy],cy_Lm[iy],&tmp_Re,&tmp_Lm); //conj cx[ix]
            *value_Re = *value_Re + tmp_Re;
            *value_Lm = *value_Lm + tmp_Lm;
            ix = ix + incx;
            iy = iy + incy;
        }
    }
    return;
}

//****************************************************************************80
//
//  Purpose:
//
//    CAXPY computes a constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//  Modified:
//
//    10 April 2006
//
//  Author:
//
//    Jack Dongarra
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in CX and CY.
//
//    Input, complex <float> CA, the multiplier of CX.
//
//    Input, complex <float> CX[], the first vector.
//
//    Input, int INCX, the increment between successive entries of CX.
//
//    Input/output, complex <float> CY[], the second vector.
//    On output, CY(*) has been replaced by CY(*) + CA * CX(*).
//
//    Input, int INCY, the increment between successive entries of CY.
//
void caxpy(int n, float ca_Re, float ca_Lm, float* cx_Re, float* cx_Lm,
           int incx, float* cy_Re, float* cy_Lm, int incy){
    int i;
    int ix;
    int iy;

    float tmp_Re;
    float tmp_Lm;

    if(n <= 0){
        return;
    }

    if(cabs1(ca_Re,ca_Lm) == 0.0){
        return;
    }

    if(incx != 1 || incy != 1){
        if (0 <= incx){
            ix = 0;
        }
        else{
            ix = (-n + 1) * incx;
        }

        if (0 <= incy){
            iy = 0;
        }
        else{
            iy = (-n + 1) * incy;
        }

        for (i = 0; i < n; i++){
            p_complex_mul(ca_Re,ca_Lm,cx_Re[ix],cx_Lm[ix],&tmp_Re,&tmp_Lm);
            cy_Re[iy] = cy_Re[iy] + tmp_Re;
            cy_Lm[iy] = cy_Lm[iy] + tmp_Lm;
            ix = ix + incx;
            iy = iy + incy;
        }
    }
    else{
        for (i = 0; i < n; i++){
            p_complex_mul(ca_Re,ca_Lm,cx_Re[i],cx_Lm[i],&tmp_Re,&tmp_Lm);
            cy_Re[i] = cy_Re[i] + tmp_Re;
            cy_Lm[i] = cy_Lm[i] + tmp_Lm;
        }
    }

    return;
}

//****************************************************************************80
//
//  Purpose:
//
//    SROTG constructs a float Givens plane rotation.
//
//  Discussion:
//
//    Given values A and B, this routine computes
//
//    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
//          = sign ( B ) if abs ( A ) <= abs ( B );
//
//    R     = SIGMA * ( A * A + B * B );
//
//    C = A / R if R is not 0
//      = 1     if R is 0;
//
//    S = B / R if R is not 0,
//        0     if R is 0.
//
//    The computed numbers then satisfy the equation
//
//    (  C  S ) ( A ) = ( R )
//    ( -S  C ) ( B ) = ( 0 )
//
//    The routine also computes
//
//    Z = S     if abs ( A ) > abs ( B ),
//      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
//      = 1     if C is 0.
//
//    The single value Z encodes C and S, and hence the rotation:
//
//    If Z = 1, set C = 0 and S = 1;
//    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
//    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    Jack Dongarra
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input/output, float *SA, *SB,  On input, SA and SB are the values
//    A and B.  On output, SA is overwritten with R, and SB is
//    overwritten with Z.
//
//    Output, float *C, *S, the cosine and sine of the Givens rotation.
//
void srotg(float *sa, float *sb, float *c, float *s){
    float r;
    float roe;
    float scale;
    float z;

    if (r4_abs(*sb) < r4_abs(*sa)){
        roe = *sa;
    }
    else{
        roe = *sb;
    }

    scale = r4_abs(*sa) + r4_abs(*sb);

    if (scale == 0.0){
        *c = 1.0;
        *s = 0.0;
        r = 0.0;
    }
    else{
        r = scale * sqrt((*sa / scale) * (*sa / scale) + (*sb / scale) * (*sb / scale));
        r = r4_sign(roe) * r;
        *c = *sa / r;
        *s = *sb / r;
        //printf("%f,%f\n",*c,*s);
    }

    if (0.0 < r4_abs(*c) && r4_abs(*c) <= *s){
        z = (double)1.0 / (double)*c;
        //printf("%f,%f\n",*c,z);
    }
    else{
        z = *s;
        //printf("%f\n",z);
    }

    *sa = r;
    *sb = z;

    return;
}


//****************************************************************************80
//
//  Purpose:
//
//    CSROT applies a plane rotation.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The cosine C and sine S are real and the vectors CX and CY are complex.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    Jack Dongarra
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, complex <float> CX[], one of the vectors to be rotated.
//
//    Input, int INCX, the increment between successive entries of CX.
//
//    Input/output, complex <float> CY[], one of the vectors to be rotated.
//
//    Input, int INCY, the increment between successive elements of CY.
//
//    Input, float C, S, parameters (presumably the cosine and sine of
//    some angle) that define a plane rotation.
//
void csrot(int n, float *cx_Re, float *cx_Lm, int incx, float *cy_Re, float *cy_Lm,
           int incy, float c, float s){
    
    float ctemp_Re,ctemp_Lm;
    int i;
    int ix;
    int iy;

    if (n <= 0){
        return;
    }

    if (incx == 1 && incy == 1){
        for (i = 0; i < n; i++){
            ctemp_Re = c * cx_Re[i] + s * cy_Re[i];
            ctemp_Lm = c * cx_Lm[i] + s * cy_Lm[i];

            cy_Re[i] = c * cy_Re[i] - s * cx_Re[i];
            cy_Lm[i] = c * cy_Lm[i] - s * cx_Lm[i];

            cx_Re[i] = ctemp_Re;
            cx_Lm[i] = ctemp_Lm;
        }
    }
    else{
        if (0 <= incx){
            ix = 0;
        }
        else{
            ix = (-n + 1) * incx;
        }

        if (0 <= incy){
            iy = 0;
        }
        else{
            iy = (-n + 1) * incy;
        }

        for (i = 0; i < n; i++){
            ctemp_Re = c * cx_Re[ix] + s * cy_Re[iy];
            ctemp_Lm = c * cx_Lm[ix] + s * cy_Lm[iy];

            cy_Re[iy] = c * cy_Re[iy] - s * cx_Re[ix];
            cy_Lm[iy] = c * cy_Lm[iy] - s * cx_Lm[ix];

            cx_Re[ix] = ctemp_Re;
            cx_Lm[ix] = ctemp_Lm;
            
            ix = ix + incx;
            iy = iy + incy;
        }
    }
    return;
}

//****************************************************************************80
//
//  Purpose:
//
//    CSVDC applies the singular value decompostion to an N by P matrix.
//
//  Discussion:
//
//    The routine reduces a complex <float> N by P matrix X, by unitary transformations
//    U and V, to diagonal form.
//
//    The diagonal elements, S(I), are the singular values of Z.  The
//    columns of U are the corresponding left singular vectors,
//    and the columns of V the right singular vectors.
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <float> X[LDX*P]; on input, the matrix whose singular value
//    decomposition is to be computed.  X is destroyed on output.
//
//    Input, int LDX, the leading dimension of X.  N <= LDX.
//
//    Input, int N, the number of rows of the matrix.
//
//    Input, int P, the number of columns of the matrix X.
//
//    Output, complex <float> S[MM], where MM = min ( N + 1, P ), the first min ( N, P )
//    entries of S contain the singular values of X arranged in descending
//    order of magnitude.
//
//    Output, complex <float> E[MM], where MM = min ( N + 1, P ),
//    ordinarily contains zeros on output.  However, see the discussion
//    of INFO for exceptions.
//
//    Output, complex <float> U[LDU*K].  If JOBA == 1 then K == n; if JOBA >= 2,
//    then K == min ( N, P ).  U contains the matrix of left singular vectors.
//    U is not referenced if JOBA == 0.  If N <= P or if JOBA > 2,
//    then U may be identified with X in the subroutine call.
//
//    Input, int LDU, the leading dimension of U.  N <= LDU.
//
//    Output, complex <float> V[LDV*P], if requested, the matrix of right singular
//    vectors.  If P <= N, V may be identified with X in the subroutine call.
//
//    Input, int LDV, the leading dimension of V.  P <= LDV.
//
//    Input, int JOB, controls the computation of the singular vectors.
//    It has the decimal expansion AB meaning:
//    A =  0, do not compute the left singular vectors.
//    A =  1, return the N left singular vectors in U.
//    A >= 2, returns the first min ( N, P ) left singular vectors in U.
//    B =  0, do not compute the right singular vectors.
//    B =  1, return the right singular vectors in V.
//
//    Output, int CSVDC, the value of INFO.  The singular values and their
//    corresponding singular vectors are correct for entries,
//    S(INFO+1), S(INFO+2), ..., S(M).  Here M = min ( N, P ).  Thus if
//    INFO == 0, all the singular values and their vectors are correct.
//    In any event, the matrix
//      B = hermitian(U)*X*V
//    is the bidiagonal matrix with the elements of S on its diagonal
//    and the elements of E on its super-diagonal.  Hermitian(U)
//    is the conjugate-transpose of U.  Thus the singular values of X
//    and B are the same.
//

int csvdc(float *x_Re, float *x_Lm, int ldx, int n, int p,
          float *s_Re, float *s_Lm, float *e_Re, float *e_Lm,
          float *u_Re, float *u_Lm, int ldu,
          float *v_Re, float *v_Lm, int ldv, int job){

    float b;
    float c;
    float cs;
    float el;
    float emm1;
    float f;
    float g;
    int i;
    int info;
    int iter;
    int j;
    int jobu;
    int k;
    int kase;
    int kk;
    int l;
    int ll;
    int lls;
    int lp1;
    int ls;
    int lu;
    int m;
    int maxit = 30;
    int mm;
    int mm1;
    int mp1;
    int nct;
    int nctp1;
    int ncu;
    int nrt;
    int nrtp1;

    float r_Re;
    float r_Lm;

    float scale;
    float shift;
    float sl;
    float sm;
    float smm1;
    float sn;

    float t_Re;
    float t_Lm;

    float t1;
    float test;
    bool wantu;
    bool wantv;

    float ztest;

    float work_Re[8];  /* fixed-size; MSVC does not support VLA (n <= N_CHANNEL = 4) */
    float work_Lm[8];

    float tmp_Re;
    float tmp_Lm;
    //
    //  Determine what is to be computed.
    //
    wantu = false;
    wantv = false;
    jobu = (job % 100) / 10;

    if (1 < jobu){
        ncu = i4_min(n, p);
    }
    else{
        ncu = n;
    }

    if (jobu != 0){
        wantu = true;
    }

    if ((job % 10) != 0){
        wantv = true;
    }
    //
    //  Reduce X to bidiagonal form, storing the diagonal elements
    //  in S and the super-diagonal elements in E.
    //
    info = 0;
    nct = i4_min(n - 1, p);
    nrt = i4_max(0, i4_min(p - 2, n));
    lu = i4_max(nct, nrt);

    #if 0
    for(int i = 0; i < 4*4; i++){
      printf("%g,%g ",x_Re[i],x_Lm[i]);
    }
    printf("\n");
    getchar();
    #endif

    for(l = 1; l <= lu; l++){
        lp1 = l + 1;
        //
        //  Compute the transformation for the L-th column and
        //  place the L-th diagonal in S(L).
        //
        if(l <= nct){
            s_Re[l - 1] = scnrm2(n - l + 1, x_Re + l - 1 + (l - 1) * ldx, x_Lm + l - 1 + (l - 1) * ldx, 1);
            s_Lm[l - 1] = 0.0;

            if(cabs1(s_Re[l - 1],s_Lm[l - 1]) != 0.0){
                if(cabs1(x_Re[l - 1 + (l - 1) * ldx],x_Lm[l - 1 + (l - 1) * ldx]) != 0.0){
                    csign2(s_Re[l - 1], s_Lm[l - 1], x_Re[l - 1 + (l - 1) * ldx], x_Lm[l - 1 + (l - 1) * ldx],&s_Re[l - 1], &s_Lm[l - 1]);
                }
                p_complex_div(1.0,0.0,s_Re[l - 1],s_Lm[l - 1],&t_Re,&t_Lm);
                //printf("%g,%g,%g,%g\n",s_Re[l-1],s_Lm[l-1],t_Re,t_Lm);
                cscal(n - l + 1, t_Re, t_Lm, x_Re + l - 1 + (l - 1) * ldx, x_Lm + l - 1 + (l - 1) * ldx, 1);
                x_Re[l - 1 + (l - 1) * ldx] = 1.0 + x_Re[l - 1 + (l - 1) * ldx];
                x_Lm[l - 1 + (l - 1) * ldx] = 0.0 + x_Lm[l - 1 + (l - 1) * ldx];
                #if 0
                for(int i = 0; i < 4*4; i++){
                    printf("%g,%g ",x_Re[i],x_Lm[i]);
                }
                printf("\n");
                getchar();
                #endif                
            }

            s_Re[l - 1] = -s_Re[l - 1];
            s_Lm[l - 1] = -s_Lm[l - 1];
        }

        for(j = lp1; j <= p; j++){
            if(l <= nct){
                if (cabs1(s_Re[l - 1],s_Lm[l - 1]) != 0.0){
                    cdotc(n - l + 1, x_Re + l - 1 + (l - 1) * ldx, x_Lm + l - 1 + (l - 1) * ldx, 1, 
                          x_Re + l - 1 + (j - 1) * ldx, x_Lm + l - 1 + (j - 1) * ldx, 1,
                          &t_Re, &t_Lm);
                    p_complex_div(-t_Re,-t_Lm,x_Re[l - 1 + (l - 1) * ldx],x_Lm[l - 1 + (l - 1) * ldx],&t_Re,&t_Lm);
                    //printf("%g,%g ",t_Re,t_Lm);
                    caxpy(n - l + 1, t_Re, t_Lm, x_Re + l - 1 + (l - 1) * ldx, x_Lm + l - 1 + (l - 1) * ldx, 1, 
                          x_Re + l - 1 + (j - 1) * ldx, x_Lm + l - 1 + (j - 1) * ldx, 1);
                }
            }
            //
            //  Place the L-th row of X into E for the
            //  subsequent calculation of the row transformation.
            //
            e_Re[j - 1] = x_Re[l - 1 + (j - 1) * ldx];
            e_Lm[j - 1] = -x_Lm[l - 1 + (j - 1) * ldx];
        }
        //
        //  Place the transformation in U for subsequent back multiplication.
        //
        if (wantu && l <= nct){
            for (i = l; i <= n; i++){
                u_Re[i - 1 + (l - 1) * ldu] = x_Re[i - 1 + (l - 1) * ldx];
                u_Lm[i - 1 + (l - 1) * ldu] = x_Lm[i - 1 + (l - 1) * ldx];
            }
        }

        if (l <= nrt){
            //
            //  Compute the L-th row transformation and place the
            //  L-th super-diagonal in E(L).
            //
            e_Re[l - 1] = scnrm2(p - l, e_Re + lp1 - 1, e_Lm + lp1 - 1, 1);
            e_Lm[l - 1] = 0.0;

            if (cabs1(e_Re[l - 1],e_Lm[l - 1]) != 0.0){
                if (cabs1(e_Re[lp1 - 1],e_Lm[lp1 - 1]) != 0.0){
                    csign2(e_Re[l - 1],e_Lm[l - 1],e_Re[lp1 - 1],e_Lm[lp1 - 1],&e_Re[l - 1],&e_Lm[l - 1]);
                }

                p_complex_div(1.0,0.0,e_Re[l - 1],e_Lm[l - 1],&t_Re,&t_Lm);
                cscal(p - l, t_Re, t_Lm, e_Re + lp1 - 1, e_Lm + lp1 -1,1);
                
                e_Re[lp1 - 1] = 1.0 + e_Re[lp1 - 1];
                e_Lm[lp1 - 1] = 0.0 + e_Lm[lp1 - 1];
            }

            e_Re[l - 1] = - e_Re[l - 1];
            e_Lm[l - 1] = e_Lm[l - 1];
            //printf("%g,%g ",e_Re[l-1],e_Lm[l-1]);

            //
            //  Apply the transformation.
            //
            if (lp1 <= n && cabs1(e_Re[l - 1],e_Lm[l - 1]) != 0.0){
                for (j = lp1; j <= n; j++){
                    work_Re[j - 1] = 0.0;
                    work_Lm[j - 1] = 0.0;
                }

                for (j = lp1; j <= p; j++){
                    caxpy(n - l, e_Re[j - 1], e_Lm[j - 1], x_Re + lp1 - 1 + (j - 1) * ldx, x_Lm + lp1 - 1 + (j - 1) * ldx, 1, 
                         work_Re + lp1 - 1, work_Lm + lp1 - 1, 1);
                }

                for (j = lp1; j <= p; j++){
                    p_complex_div(-e_Re[j - 1],-e_Lm[j - 1],e_Re[lp1 - 1],e_Lm[lp1 - 1],&tmp_Re,&tmp_Lm);

                    caxpy(n - l, tmp_Re, -tmp_Lm, work_Re + lp1 - 1,work_Lm + lp1 - 1, 1, 
                          x_Re + lp1 - 1 + (j - 1) * ldx, x_Lm + lp1 - 1 + (j - 1) * ldx, 1);
                }
            }
            //
            //  Place the transformation in V for subsequent back multiplication.
            //
            if (wantv){
                for (i = lp1; i <= p; i++){
                    v_Re[i - 1 + (l - 1) * ldv] = e_Re[i - 1];
                    v_Lm[i - 1 + (l - 1) * ldv] = e_Lm[i - 1];
                    //printf("%g,%g ",v_Re[i - 1 + (l - 1) * ldv],v_Lm[i - 1 + (l - 1) * ldv]);
                }
            }
        }
    }
    #if 0
    printf("\n");
    getchar();
    #endif
    //
    //  Set up the final bidiagonal matrix of order M.
    //
    m = i4_min(p, n + 1);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;

    if (nct < p){
        s_Re[nctp1 - 1] = x_Re[nctp1 - 1 + (nctp1 - 1) * ldx];
        s_Lm[nctp1 - 1] = x_Lm[nctp1 - 1 + (nctp1 - 1) * ldx];
    }

    if (n < m){
        s_Re[m - 1] = 0.0;
        s_Lm[m - 1] = 0.0;
    }

    if (nrtp1 < m){
        e_Re[nrtp1 - 1] = x_Re[nrtp1 - 1 + (m - 1) * ldx];
        e_Lm[nrtp1 - 1] = x_Lm[nrtp1 - 1 + (m - 1) * ldx];
    }

    e_Re[m - 1] = 0.0;
    e_Lm[m - 1] = 0.0;
    //
    //  If required, generate U.
    //
    if (wantu){
        for (j = nctp1; j <= ncu; j++){
            for (i = 1; i <= n; i++){
                u_Re[i - 1 + (j - 1) * ldu] = 0.0;
                u_Lm[i - 1 + (j - 1) * ldu] = 0.0;
            }
            u_Re[j - 1 + (j - 1) * ldu] = 1.0;
            u_Lm[j - 1 + (j - 1) * ldu] = 0.0;
        }
        for (ll = 1; ll <= nct; ll++){
            l = nct - ll + 1;

            if (cabs1(s_Re[l - 1],s_Lm[l - 1]) != 0.0){
                lp1 = l + 1;
                
                for (j = l + 1; j <= ncu; j++){
                    cdotc(n - l + 1, u_Re + l - 1 + (l - 1) * ldu, u_Lm + l - 1 + (l - 1) * ldu, 1, 
                          u_Re + l - 1 + (j - 1) * ldu, u_Lm + l - 1 + (j - 1) * ldu, 1, &t_Re, &t_Lm);
                    //printf("%g,%g ",t_Re,t_Lm);      
                    p_complex_div(-t_Re,-t_Lm,u_Re[l - 1 + (l - 1) * ldu],u_Lm[l - 1 + (l - 1) * ldu],&t_Re,&t_Lm);
                    //printf("%f,%f ",t_Re,t_Lm);
                    caxpy(n - l + 1, t_Re, t_Lm, u_Re + l - 1 + (l - 1) * ldu, u_Lm + l - 1 + (l - 1) * ldu,1,
                          u_Re + l - 1 + (j - 1) * ldu, u_Lm + l - 1 + (j - 1) * ldu, 1);
                }

                cscal(n - l + 1, -1.0, 0.0, u_Re + l - 1 + (l - 1) * ldu, u_Lm + l - 1 + (l - 1) * ldu, 1);

                u_Re[l - 1 + (l - 1) * ldu] = 1.0 + u_Re[l - 1 + (l - 1) * ldu];
                u_Lm[l - 1 + (l - 1) * ldu] = 0.0 + u_Lm[l - 1 + (l - 1) * ldu];

                for (i = 1; i <= l - 1; i++){
                    u_Re[i - 1 + (l - 1) * ldu] = 0.0;
                    u_Lm[i - 1 + (l - 1) * ldu] = 0.0;
                }
            }
            else{
                for (i = 1; i <= n; i++){
                    u_Re[i - 1 + (l - 1) * ldu] = 0.0;
                    u_Lm[i - 1 + (l - 1) * ldu] = 0.0;
                }
                u_Re[l - 1 + (l - 1) * ldu] = 1.0;
                u_Lm[l - 1 + (l - 1) * ldu] = 0.0;
            }
        }
    }
    #if 0
    for(int i = 0; i < 4*4; i++){
      printf("%f,%f ",u_Re[i],u_Lm[i]);
    }
    #endif
    //printf("\n");
    //getchar();
    //
    //  If it is required, generate V.
    //
    if (wantv){
        for (ll = 1; ll <= p; ll++){
            l = p - ll + 1;
            lp1 = l + 1;

            if (l <= nrt){
                if (cabs1(e_Re[l - 1],e_Lm[l - 1]) != 0.0){
                    for (j = lp1; j <= p; j++){

                        cdotc(p - l, v_Re + lp1 - 1 + (l - 1) * ldv,v_Lm + lp1 - 1 + (l - 1) * ldv, 1, 
                              v_Re + lp1 - 1 + (j - 1) * ldv, v_Lm + lp1 - 1 + (j - 1) * ldv, 1, &t_Re, &t_Lm);
                        
                        p_complex_div(-t_Re,-t_Lm,v_Re[lp1 - 1 + (l - 1) * ldv],v_Lm[lp1 - 1 + (l - 1) * ldv],&t_Re,&t_Lm);
                        //printf("%f,%f ",t_Re,t_Lm); 

                        caxpy(p - l, t_Re, t_Lm, v_Re + lp1 - 1 + (l - 1) * ldv, v_Lm + lp1 - 1 + (l - 1) * ldv, 1, 
                                v_Re + lp1 - 1 + (j - 1) * ldv, v_Lm + lp1 - 1 + (j - 1) * ldv, 1);
                    }
                }
            }
            for (i = 1; i <= p; i++){
                v_Re[i - 1 + (l - 1) * ldv] = 0.0;
                v_Lm[i - 1 + (l - 1) * ldv] = 0.0;
            }
            v_Re[l - 1 + (l - 1) * ldv] = 1.0;
            v_Lm[l - 1 + (l - 1) * ldv] = 0.0;
        }
    }

    #if 0
    for(int i = 0; i < 4*4; i++){
      printf("%f,%f ",v_Re[i],v_Lm[i]);
    }
    printf("\n");
    getchar();    
    #endif

    //
    //  Transform S and E so that they are real.
    //
    for (i = 1; i <= m; i++){
        if (cabs1(s_Re[i - 1],s_Lm[i - 1]) != 0.0){
            t_Re = p_complex_abs(s_Re[i - 1],s_Lm[i - 1]);
            t_Lm = 0.0;
            p_complex_div(s_Re[i - 1],s_Lm[i - 1],t_Re,t_Lm,&r_Re,&r_Lm);
            s_Re[i - 1] = t_Re;
            s_Lm[i - 1] = t_Lm;

            if (i < m){
                p_complex_div(e_Re[i - 1],e_Lm[i - 1],r_Re,r_Lm,&e_Re[i - 1],&e_Lm[i - 1]);
            }

            if (wantu){
                cscal(n, r_Re, r_Lm, u_Re + 0 + (i - 1) * ldu, u_Lm + 0 + (i - 1) * ldu, 1);
            }
        }

        if (i == m){
            break;
        }

        if (cabs1(e_Re[i - 1],e_Lm[i - 1]) != 0.0){
            t_Re = p_complex_abs(e_Re[i - 1],e_Lm[i - 1]);
            t_Lm = 0.0;
            p_complex_div(t_Re,t_Lm,e_Re[i - 1],e_Lm[i - 1],&r_Re,&r_Lm);
            e_Re[i - 1] = t_Re;
            e_Lm[i - 1] = t_Lm;
            p_complex_mul(s_Re[i],s_Lm[i],r_Re,r_Lm,&s_Re[i],&s_Lm[i]);

            if (wantv){
                cscal(p, r_Re, r_Lm, v_Re + 0 + i * ldv, v_Lm + 0 + i * ldv, 1);
            }
        }
    }

    #if 0
    for(int i = 0; i < 4+4; i++){
      printf("%f,%f ",e_Re[i],e_Lm[i]);
    }
    printf("\n");
    getchar();    
    #endif

    //
    //  Main iteration loop for the singular values.
    //
    mm = m;
    iter = 0;

    for (;;){
        //
        //  Quit if all the singular values have been found.
        //
        if (m == 0){
            break;
        }
        //
        //  If too many iterations have been performed, set flag and return.
        //
        if (maxit <= iter){
            info = m;
            break;
        }
        //
        //  This section of the program inspects for negligible elements in S and E.
        //
        //  On completion, the variables KASE and L are set as follows.
        //
        //  KASE = 1     if S(M) and E(L-1) are negligible and L < M
        //  KASE = 2     if S(L) is negligible and L < M
        //  KASE = 3     if E(L-1) is negligible, L < M, and
        //               S(L), ..., S(M) are not negligible (QR step).
        //  KASE = 4     if E(M-1) is negligible (convergence).
        //
        for (ll = 1; ll <= m; ll++){
            l = m - ll;

            if (l == 0){
                break;
            }

            test = p_complex_abs(s_Re[l - 1],s_Lm[l - 1]) + p_complex_abs(s_Re[l],s_Lm[l]);
            ztest = test + p_complex_abs(e_Re[l - 1],e_Lm[l - 1]);

            if (ztest == test){
                e_Re[l - 1] = 0.0;
                e_Lm[l - 1] = 0.0;
                break;
            }
        }

        if (l == m - 1){
            kase = 4;
        }
        else{
            lp1 = l + 1;
            mp1 = m + 1;

            for (lls = lp1; lls <= mp1; lls++){
                ls = m - lls + lp1;

                if (ls == l){
                    break;
                }

                test = 0.0;

                if (ls != m){
                    test = test + p_complex_abs(e_Re[ls - 1],e_Lm[ls - 1]);
                }

                if (ls != l + 1){
                    test = test + p_complex_abs(e_Re[ls - 2],e_Lm[ls - 2]);
                }

                ztest = test + p_complex_abs(s_Re[ls - 1],s_Lm[ls - 1]);

                if (ztest == test){
                    s_Re[ls - 1] = 0.0;
                    s_Lm[ls - 1] = 0.0;
                    break;
                }
            }
            if (ls == l){
                kase = 3;
            }
            else if (ls == m){
                kase = 1;
            }
            else{
                kase = 2;
                l = ls;
            }
        }

        l = l + 1;
        //
        //  Deflate negligible S(M).
        //
        if (kase == 1){
            mm1 = m - 1;
            f = e_Re[m - 2];

            e_Re[m - 2] = 0.0;
            e_Lm[m - 2] = 0.0;

            for (kk = 1; kk <= mm1; kk++){
                k = mm1 - kk + l;
                t1 = s_Re[k - 1];
                srotg(&t1, &f, &cs, &sn);
                
                s_Re[k - 1] = t1;
                s_Lm[k - 1] = 0.0;

                if (k != l){
                    f = -sn * e_Re[k - 2];
                    e_Re[k - 2] = cs * e_Re[k - 2];
                    e_Lm[k - 2] = cs * e_Lm[k - 2];
                }

                if (wantv){
                    csrot(p, v_Re + 0 + (k - 1) * ldv,v_Lm + 0 + (k - 1) * ldv, 1, 
                            v_Re + 0 + (m - 1) * ldv,v_Lm + 0 + (m - 1) * ldv, 1, cs, sn);
                }
            }
        }
        //
        //  Split at negligible S(L).
        //
        else if (kase == 2){
            f = e_Re[l - 2];

            e_Re[l - 2] = 0.0;
            e_Lm[l - 2] = 0.0;

            for (k = l; k <= m; k++){
                t1 = s_Re[k - 1];
                srotg(&t1, &f, &cs, &sn);
                s_Re[k - 1] = t1;
                s_Lm[k - 1] = 0.0;

                f = -sn * e_Re[k - 1];

                e_Re[k - 1] = cs * e_Re[k - 1];
                e_Lm[k - 1] = cs * e_Lm[k - 1];

                if (wantu){
                    csrot(n, u_Re + 0 + (k - 1) * ldu, u_Lm + 0 + (k - 1) * ldu, 1, 
                         u_Re + 0 + (l - 2) * ldu, u_Lm + 0 + (l - 2) * ldu, 1, cs, sn);
                }
            }
        }
        //
        //  Perform one QR step.
        //
        else if (kase == 3){
            //
            //  Calculate the shift.
            //
            scale = r4_max(p_complex_abs(s_Re[m - 1],s_Lm[m - 1]),
                           r4_max(p_complex_abs(s_Re[m - 2],s_Lm[m - 2]),
                                  r4_max(p_complex_abs(e_Re[m - 2],e_Lm[m - 2]),
                                         r4_max(p_complex_abs(s_Re[l - 1],s_Lm[l - 1]), p_complex_abs(e_Re[l - 1],e_Lm[l - 1])))));

            sm = s_Re[m - 1] / scale;
            smm1 = s_Re[m - 2] / scale;
            emm1 = e_Re[m - 2] / scale;
            sl = s_Re[l - 1] / scale;
            el = e_Re[l - 1] / scale;
            b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
            c = (sm * emm1) * (sm * emm1);
            shift = 0.0;

            if (b != 0.0 || c != 0.0){
                shift = sqrt(b * b + c);
                if (b < 0.0){
                    shift = -shift;
                }
                shift = c / (b + shift);
            }

            f = (sl + sm) * (sl - sm) + shift;
            g = sl * el;

            //
            //  Chase zeros.
            //
            mm1 = m - 1;

            for (k = l; k <= mm1; k++){
                //printf("%f,%f,%f,%f ",f,g,cs,sn);
                srotg(&f, &g, &cs, &sn);

                if (k != l){
                    e_Re[k - 2] = f;
                    e_Lm[k - 2] = 0.0;
                }

                f = cs * s_Re[k - 1] + sn * e_Re[k - 1];
                e_Re[k - 1] = cs * e_Re[k - 1] - sn * s_Re[k - 1];
                e_Lm[k - 1] = cs * e_Lm[k - 1] - sn * s_Lm[k - 1];
                g = sn * s_Re[k];
                s_Re[k] = cs * s_Re[k];
                s_Lm[k] = cs * s_Lm[k];

                if (wantv){
                    csrot(p, v_Re + 0 + (k - 1) * ldv, v_Lm + 0 + (k - 1) * ldv, 1, 
                    v_Re + 0 + k * ldv, v_Lm + 0 + k * ldv, 1, cs, sn);
                }
                
                srotg(&f, &g, &cs, &sn);

                s_Re[k - 1] = f;
                s_Lm[k - 1] = 0.0;

                f = cs * e_Re[k - 1] + sn * s_Re[k];
                s_Re[k] = -sn * e_Re[k - 1] + cs * s_Re[k];
                s_Lm[k] = -sn * e_Lm[k - 1] + cs * s_Lm[k];

                g = sn * e_Re[k];

                e_Re[k] = cs * e_Re[k];
                e_Lm[k] = cs * e_Lm[k];

                if (wantu && k < n){
                    csrot(n, u_Re + 0 + (k - 1) * ldu, u_Lm + 0 + (k - 1) * ldu, 1, 
                          u_Re + 0 + k * ldu, u_Lm + 0 + k * ldu, 1, cs, sn);
                }
            }
            e_Re[m - 2] = f;
            e_Lm[m - 2] = 0.0;
            iter = iter + 1;

        }
        //
        //  Convergence.
        //
        else if (kase == 4){
            //
            //  Make the singular value positive.
            //
            if (s_Re[l - 1] < 0.0){
                s_Re[l - 1] = -s_Re[l - 1];
                s_Lm[l - 1] = -s_Lm[l - 1];

                if (wantv){
                    cscal(p, -1.0, 0.0, v_Re + 0 + (l - 1) * ldv, v_Lm + 0 + (l - 1) * ldv, 1);
                }
            }
            //
            //  Order the singular values.
            //
            while (l != mm){
                if (s_Re[l] <= s_Re[l - 1]){
                    break;
                }

                t_Re = s_Re[l - 1];
                t_Lm = s_Lm[l - 1];

                s_Re[l - 1] = s_Re[l];
                s_Lm[l - 1] = s_Lm[l];

                s_Re[l] = t_Re;
                s_Lm[l] = t_Lm;

                if (wantv && l < p){
                    cswap(p, v_Re + 0 + (l - 1) * ldv, v_Lm + 0 + (l - 1) * ldv, 1, 
                          v_Re + 0 + l * ldv, v_Lm + 0 + l * ldv, 1);
                }

                if (wantu && l < n){
                    cswap(n, u_Re + 0 + (l - 1) * ldu, u_Lm + 0 + (l - 1) * ldu, 1, 
                          u_Re + 0 + l * ldu, u_Lm + 0 + l * ldu, 1);
                }
                l = l + 1;
            }
            iter = 0;
            m = m - 1;
        }
    }

    return info;
}

/**
   @brief calculate the Moore-Penrose pseudoinverse matrix

   @param gsl_matrix_complex *A[in] an input matrix. A[M][N] where M > N
   @param gsl_matrix_complex *invA[out] an pseudoinverse matrix of A
   @param float dThreshold[in]

   @note if M>N, invA * A = E. Otherwise, A * invA = E.
 */

#define N_CHANNEL 4

bool pseudoinverse(double (*A_Re)[N_CHANNEL], double (*A_Lm)[N_CHANNEL],
                  double (*invA_Re)[N_CHANNEL], double (*invA_Lm)[N_CHANNEL], float dThreshold)
{

    float a_Re[N_CHANNEL * N_CHANNEL];
    float a_Lm[N_CHANNEL * N_CHANNEL];

    float s_Re[N_CHANNEL + N_CHANNEL];
    float s_Lm[N_CHANNEL + N_CHANNEL];

    float e_Re[N_CHANNEL + N_CHANNEL];
    float e_Lm[N_CHANNEL + N_CHANNEL];

    float u_Re[N_CHANNEL * N_CHANNEL];
    float u_Lm[N_CHANNEL * N_CHANNEL];

    float v_Re[N_CHANNEL * N_CHANNEL];
    float v_Lm[N_CHANNEL * N_CHANNEL];

    float tmp_Re,tmp_Lm;

    int info;

    bool ret = true;

    for (int i = 0; i < N_CHANNEL; i++){
        for (int j = 0; j < N_CHANNEL; j++){
            a_Re[i + j * N_CHANNEL] = (float)A_Re[i][j];
            a_Lm[i + j * N_CHANNEL] = (float)A_Lm[i][j];
        }
    }

    info = csvdc(a_Re, a_Lm, N_CHANNEL, N_CHANNEL, N_CHANNEL, s_Re, s_Lm, e_Re, e_Lm, u_Re, u_Lm, N_CHANNEL, v_Re, v_Lm, N_CHANNEL, 11);
    if (info != 0){
        printf("CSVDC returned nonzero INFO\n"); 
        ret = false; //return(true);
    }

    #if 0
    for(int i = 0; i < 4*4; i++){
      printf("%f,%f ",u_Re[i],u_Lm[i]);
    }
    printf("\n");
    getchar();
    #endif

    for (int k = 0; k < N_CHANNEL; k++){
        if (p_complex_abs(s_Re[k],s_Lm[k]) < dThreshold){
            s_Re[k] = 0.0;
            s_Lm[k] = 0.0;
            fprintf(stderr, "pseudoinverse: s[%d] = 0 because of %e < %e\n", k, p_complex_abs(s_Re[k],s_Lm[k]), dThreshold);
            ret = false;
        }
        else{
            p_complex_div(1.0,0.0,s_Re[k],s_Lm[k],&s_Re[k],&s_Lm[k]);
        }
    }

    for (int i = 0; i < N_CHANNEL; i++){
        for (int j = 0; j < N_CHANNEL; j++){
            float x_Re = 0.0;
            float x_Lm = 0.0;

            for (int k = 0; k < N_CHANNEL /*minN*/; k++){

                p_complex_mul(v_Re[j + k * N_CHANNEL],v_Lm[j + k * N_CHANNEL],s_Re[k],s_Lm[k],&tmp_Re,&tmp_Lm);
                p_complex_mul(tmp_Re,tmp_Lm,u_Re[i + k * N_CHANNEL],-u_Lm[i + k * N_CHANNEL], &tmp_Re,&tmp_Lm);

                x_Re = x_Re + tmp_Re;
                x_Lm = x_Lm + tmp_Lm;
            }
            invA_Re[j][i] = x_Re;
            invA_Lm[j][i] = x_Lm;
        }
    }
    return (ret);
}