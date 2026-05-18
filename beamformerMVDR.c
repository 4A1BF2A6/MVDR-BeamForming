#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "M128-m4-r1.h"
#include "fft_radix.h"
#include "mvdr_weights.h"

// M_ : 256 number of subband
// m_ : 4 filter length
// R_ : 2
#define M_ 128
#define m_ 4
#define R_ 2
#define D_ 64

#define N_CHANNEL 4
#define NC 1

#define M_PI 3.14159265358979323846
#define SOUNDSPEED 343740.0
#define SAMPLERATE 16000

typedef struct real_buffer_float{
	float data_float[m_*R_][M_];
	int zeros;

} *RealBuffer_Float;

RealBuffer_Float real_buffer_float_create(){

	RealBuffer_Float real_buf_f = (RealBuffer_Float)malloc(sizeof(struct real_buffer_float));

	int i,j;
	for(i = 0; i < m_*R_; i++){
		for(j = 0; j < M_; j++){
			real_buf_f->data_float[i][j] = 0.0f;
		}
	}

	real_buf_f->zeros = m_*R_ - 1;

	return real_buf_f;
}

int real_buffer_float_insert(RealBuffer_Float real_buf_f,float *data){
	real_buf_f->zeros = (real_buf_f->zeros + 1) % (m_*R_);
	memcpy(real_buf_f->data_float[real_buf_f->zeros],data,sizeof(float)*M_);
}

int real_buffer_float_insert_reverse(RealBuffer_Float real_buf_f,float *data){
	real_buf_f->zeros = (real_buf_f->zeros + 1) % (m_*R_);
	for(int i = 0; i < M_; i++){
		real_buf_f->data_float[real_buf_f->zeros][i] = data[M_-1-i];
	}
}

float real_buffer_float_get(RealBuffer_Float real_buf_f,int timeX, int binX){
	int head = (real_buf_f->zeros + m_*R_ - timeX) % (m_*R_); //最新的数据
	return real_buf_f->data_float[head][binX];
}

typedef struct real_buffer_double{
	double data_double[m_*R_][M_];
	int zeros;

} *RealBuffer_Double;

RealBuffer_Double real_buffer_double_create(){

	RealBuffer_Double real_buf_d = (RealBuffer_Double)malloc(sizeof(struct real_buffer_double));

	int i,j;
	for(i = 0; i < m_*R_; i++){
		for(j = 0; j < M_; j++){
			real_buf_d->data_double[i][j] = 0.0f;
		}
	}

	real_buf_d->zeros = m_*R_ - 1;

	return real_buf_d;
}

int real_buffer_double_insert(RealBuffer_Double real_buf_d,double *data){
	real_buf_d->zeros = (real_buf_d->zeros + 1) % (m_*R_);
	memcpy(real_buf_d->data_double[real_buf_d->zeros],data,sizeof(double)*M_);
}

int real_buffer_double_insert_reverse(RealBuffer_Double real_buf_d,double *data){
	real_buf_d->zeros = (real_buf_d->zeros + 1) % (m_*R_);
	for(int i = 0; i < M_; i++){
		real_buf_d->data_double[real_buf_d->zeros][i] = data[M_-1-i];
	}
}

double real_buffer_double_get(RealBuffer_Double real_buf_d,int timeX, int binX){
	int head = (real_buf_d->zeros + m_*R_ - timeX) % (m_*R_); //最新的数据
	return real_buf_d->data_double[head][binX];
}

/**********************************AnalysisBank*********************************/
typedef struct _OverSampledDFTAnalysisBank_{
	RealBuffer_Float analysis_buf;
	float audio_data_overlap_float[M_];
	double polyphase_output_[M_*2];

}*OverSampledDFTAnalysisBank;


OverSampledDFTAnalysisBank overSampled_dtf_analysisBank_created(){
	OverSampledDFTAnalysisBank analysisBank = (OverSampledDFTAnalysisBank)malloc(sizeof(struct _OverSampledDFTAnalysisBank_));
	analysisBank->analysis_buf = real_buffer_float_create();
	for(int i = 0; i < D_; i++){
		analysisBank->audio_data_overlap_float[i] = 0.0f;
	}
	return analysisBank;
}

int overSampled_dtf_analysisBank_process(OverSampledDFTAnalysisBank analysisBank,short *audio_data,double *output_data){
	for(int i = 0; i < D_; i++){
		analysisBank->audio_data_overlap_float[i+D_] = (float)audio_data[i];
	}

	#if 0
	for(int i = 0; i < 128; i++){
		printf("%d ",audio_data[i]);
	}
	printf("\n");
	getchar();
	#endif

	#if 0
	for(int i = 0; i < 256; i++){
		printf("%g ",analysisBank->audio_data_overlap_float[i]);
	}
	printf("\n");
	getchar();
	#endif

	real_buffer_float_insert_reverse(analysisBank->analysis_buf,analysisBank->audio_data_overlap_float);

	for(int i = 0; i < D_; i++){
		analysisBank->audio_data_overlap_float[i] = analysisBank->audio_data_overlap_float[i+D_];
	}

	// calculate outputs of polyphase filters
	// M_ : 256 number of subband
	// m_ : 4 filter length
	// R_ : 2
	for (int m = 0; m < M_; m++) {
		double sum  = 0.0;
		for (int k = 0; k < m_; k++){
			//从1024个滤波器中选择对应的系数 0+256*0 0+256*1 0+256*2 0+256*3 每次取4个滤波器系数
			//R_ * k 0 2 4 6 第0,2,4,6队列，按列取  
			sum += h_M128_m4_r1[m+M_*k] * real_buffer_float_get(analysisBank->analysis_buf,R_ * k, m);//长度为4的滤波器输出的结果
			//printf("%g ",sum);
		}

		//polyphase_output_ 256*2 double
		//256点IDFT 先准备好复数数据
		analysisBank->polyphase_output_[2*m]   = sum;
		analysisBank->polyphase_output_[2*m+1] = 0.0;
	}
	//printf("\n");
	//getchar();

	//gsl_fft_complex_radix2_backward(analysisBank->polyphase_output_, /* stride= */ 1, M_);
	fft_radix2_transform(analysisBank->polyphase_output_,1, M_, fft_radix_backward);

	#if 0
	for(int i = 0; i < 256; i++){
		printf("%g,%g ",analysisBank->polyphase_output_[2*i],analysisBank->polyphase_output_[2*i+1]);
	}
	printf("\n");
	getchar();
	#endif	
	memcpy(output_data,analysisBank->polyphase_output_,sizeof(double)*M_*2);

	return 0;
}

/**********************************SynthesisBank********************************/
typedef struct _OverSampledDFTSynthesisBank_{
	RealBuffer_Double synthesis_buf;
	double synthesis_overlap[2][M_];
	int count;

}*OverSampledDFTSynthesisBank;


OverSampledDFTSynthesisBank overSampled_dtf_synthesisBank_created(){
	OverSampledDFTSynthesisBank synthesisBank = (OverSampledDFTSynthesisBank)malloc(sizeof(struct _OverSampledDFTSynthesisBank_));
	synthesisBank->synthesis_buf = real_buffer_double_create();
	for(int i = 0; i < M_; i++){
		synthesisBank->synthesis_overlap[0][i] = 0.0f;
		synthesisBank->synthesis_overlap[1][i] = 0.0f;
	}
	synthesisBank->count = 0;

	return synthesisBank;
}

int overSampled_dtf_synthesisBank_process(OverSampledDFTSynthesisBank synthesisBank,double *input_data,short *output_data){
	double polyphase_input_real_[M_];
	double synthesis_sum[M_];
	float  synthesis_output[D_];

	//gsl_fft_complex_radix2_forward(input_data, /* stride= */ 1, M_);
	fft_radix2_transform(input_data,1, M_, fft_radix_forward);
	
	for(int i = 0; i < M_; i++){
		polyphase_input_real_[i] = input_data[2*i];
	}

	real_buffer_double_insert(synthesisBank->synthesis_buf,polyphase_input_real_);

	if(synthesisBank->count >= (R_*m_-1)){
		// M_ 256 m_ 4
		// R_ * k 0 2 4 6
		// buffer_里面已经准备好了256*8的数据了
		// polyphase(M_ - m - 1, k)   255+256*0 255+256*1 255+256*2 255+256*3=1023 从最后面开始取 
		for (int m = 0; m < M_; m++) {
			double sum  = 0.0;
			for (int k = 0; k < m_; k++){
				sum += g_M128_m4_r1[M_-m-1+M_*k] * real_buffer_double_get(synthesisBank->synthesis_buf,R_ * k, m);
			}
			synthesis_sum[m] = sum;
			//printf("%g ",sum);
		}
		//printf("\n");
		//getchar();
		memcpy(synthesisBank->synthesis_overlap[1],synthesisBank->synthesis_overlap[0],sizeof(double)*M_);
		memcpy(synthesisBank->synthesis_overlap[0],synthesis_sum,sizeof(double)*M_);

		for(int i = 0; i < D_; i++){
			synthesis_output[i] = 0.0f;
		}

		for (int sampX = 0; sampX < R_; sampX++){
			for (int d = 0; d < D_; d++){
				// (D_ - d - 1) 127,126,125....
				// R_ - sampX - 1   1,0
				// d + sampX * D_   128,129,130,...,255    0,1,2,...,127
				synthesis_output[D_-d -1] = synthesis_output[D_-d-1] + synthesisBank->synthesis_overlap[R_-sampX-1][d+sampX*D_];
			}
		}
		
		#if 0
		for(int i = 0; i < 128; i++){
			printf("%g ",synthesis_output[i]);
		}
		printf("\n");
		getchar();
		#endif

		for(int i = 0; i < D_; i++){
			output_data[i] = synthesis_output[i];
		}

		return 1;
	}
	else{
		synthesisBank->count++;
		return 0;
	}
}

#if 0
float micPos[4][3] = {
	{-32.25,0.0,0.0},
	{0.0,-32.25,0.0},
	{32.25,0.0,0.0},
	{0.0,32.25,0.0}
};
#endif

#if 1
float micPos[4][3] = {
	{-33, 0, 0.000},
	{-11, 0, 0.000},
	{ 11, 0, 0.000},
	{ 33, 0, 0.000},
};
#endif

void calcDelaysPolar2(float azimuth, float elevation, float micPos[][3], float *delays,int N){
	float c_x = - sin( elevation ) * cos( azimuth );
	float c_y = - sin( elevation ) * sin( azimuth );
	float c_z = - cos( elevation );

	for(int i = 0; i < N; i++){
		delays[i] = (c_x * micPos[i][0] + c_y * micPos[i][1] + c_z * micPos[i][2]) / SOUNDSPEED;
		//printf("%f\n",delay[i]);
	}
}

void complex_polar(double r, double theta, double *Re,double *Lm){
	*Re = r*cos(theta);
	*Lm = r*sin(theta);
}

/**********************************MVDR Releated*******************************/
#define MVDR_DBL_EPSILON        2.2204460492503131e-16
#define MVDR_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define MVDR_SQRT_DBL_EPSILON   1.4901161193847656e-08

struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

/* Chebyshev expansion for f(t) = sinc((t+1)/2), -1 < t < 1
 */
static double sinc_data[17] = {
  1.133648177811747875422,
 -0.532677564732557348781,
 -0.068293048346633177859,
  0.033403684226353715020,
  0.001485679893925747818,
 -0.000734421305768455295,
 -0.000016837282388837229,
  0.000008359950146618018,
  0.000000117382095601192,
 -0.000000058413665922724,
 -0.000000000554763755743,
  0.000000000276434190426,
  0.000000000001895374892,
 -0.000000000000945237101,
 -0.000000000000004900690,
  0.000000000000002445383,
  0.000000000000000009925
};

static cheb_series sinc_cs = {
  sinc_data,
  16,
  -1, 1,
  10
};

/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (sin(x)/x - 1)/(x*x)
 */
static double sin_data[12] = {
  -0.3295190160663511504173,
   0.0025374284671667991990,
   0.0006261928782647355874,
  -4.6495547521854042157541e-06,
  -5.6917531549379706526677e-07,
   3.7283335140973803627866e-09,
   3.0267376484747473727186e-10,
  -1.7400875016436622322022e-12,
  -1.0554678305790849834462e-13,
   5.3701981409132410797062e-16,
   2.5984137983099020336115e-17,
  -1.1821555255364833468288e-19
};

static cheb_series sin_cs = {
  sin_data,
  11,
  -1, 1,
  11
};

/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (2(cos(x) - 1)/(x^2) + 1) / x^2
 */
static double cos_data[11] = {
  0.165391825637921473505668118136,
 -0.00084852883845000173671196530195,
 -0.000210086507222940730213625768083,
  1.16582269619760204299639757584e-6,
  1.43319375856259870334412701165e-7,
 -7.4770883429007141617951330184e-10,
 -6.0969994944584252706997438007e-11,
  2.90748249201909353949854872638e-13,
  1.77126739876261435667156490461e-14,
 -7.6896421502815579078577263149e-17,
 -3.7363121133079412079201377318e-18
};
static cheb_series cos_cs = {
  cos_data,
  10,
  -1, 1,
  10
};

struct sf_result_struct {
	double val;
	double err;
};
typedef struct sf_result_struct sf_result;

static inline int cheb_eval_e(const cheb_series * cs, const double x, sf_result * result){
	int j;
	double d = 0.0;
	double dd = 0.0;

	double y = (2.0 * x - cs->a - cs->b) / (cs->b - cs->a);
	double y2 = 2.0 * y;

	double e = 0.0;

	for (j = cs->order; j >= 1; j--)
	{
		double temp = d;
		d = y2 * d - dd + cs->c[j];
		e += fabs(y2 * temp) + fabs(dd) + fabs(cs->c[j]);
		dd = temp;
	}

	{
		double temp = d;
		d = y * d - dd + 0.5 * cs->c[0];
		e += fabs(y * temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
	}

	result->val = d;
	result->err = MVDR_DBL_EPSILON * e + fabs(cs->c[cs->order]);

	return 0;
}

#define MVDR_SIGN(x)    ((x) >= 0.0 ? 1 : -1)
#define MVDR_IS_ODD(n)  ((n) & 1)

int sf_sin_e(double x, sf_result * result){

	const double P1 = 7.85398125648498535156e-1;
	const double P2 = 3.77489470793079817668e-8;
	const double P3 = 2.69515142907905952645e-15;

	const double sgn_x = MVDR_SIGN(x);
	const double abs_x = fabs(x);

	if (abs_x < MVDR_ROOT4_DBL_EPSILON){
		const double x2 = x * x;
		result->val = x * (1.0 - x2 / 6.0);
		result->err = fabs(x * x2 * x2 / 100.0);
		return 0;
	}
	else{
		double sgn_result = sgn_x;
		double y = floor(abs_x / (0.25 * M_PI));
		int octant = y - ldexp(floor(ldexp(y, -3)), 3);
		int stat_cs;
		double z;

		if (MVDR_IS_ODD(octant)){
			octant += 1;
			octant &= 07;
			y += 1.0;
		}

		if (octant > 3){
			octant -= 4;
			sgn_result = -sgn_result;
		}

		z = ((abs_x - y * P1) - y * P2) - y * P3;

		if (octant == 0){
			sf_result sin_cs_result;
			const double t = 8.0 * fabs(z) / M_PI - 1.0;
			stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
			result->val = z * (1.0 + z * z * sin_cs_result.val);
		}
		else{
			/* octant == 2 */
			sf_result cos_cs_result;
			const double t = 8.0 * fabs(z) / M_PI - 1.0;
			stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
			result->val = 1.0 - 0.5 * z * z * (1.0 - z * z * cos_cs_result.val);
		}

		result->val *= sgn_result;

		if (abs_x > 1.0 / MVDR_DBL_EPSILON){
			result->err = fabs(result->val);
		}
		else if (abs_x > 100.0 / MVDR_SQRT_DBL_EPSILON){
			result->err = 2.0 * abs_x * MVDR_DBL_EPSILON * fabs(result->val);
		}
		else if (abs_x > 0.1 / MVDR_SQRT_DBL_EPSILON){
			result->err = 2.0 * MVDR_SQRT_DBL_EPSILON * fabs(result->val);
		}
		else{
			result->err = 2.0 * MVDR_DBL_EPSILON * fabs(result->val);
		}

		return stat_cs;
	}
}

static int sf_sinc(double x,sf_result *result){
	const double ax = fabs(x);
	if(ax < 0.8){
		/* Do not go to the limit of the fit since
		* there is a zero there and the Chebyshev
		* accuracy will go to zero.
		*/
		return cheb_eval_e(&sinc_cs, 2.0*ax-1.0, result);
	}
    else if(ax < 100.0){
		/* Small arguments are no problem.
		* We trust the library sin() to
		* roughly machine precision.
		*/
		result->val = sin(M_PI * ax)/(M_PI * ax);
		result->err = 2.0 * MVDR_DBL_EPSILON * fabs(result->val);
		return 0;
	}
    else{
		/* Large arguments must be handled separately.
		*/
		const double r = M_PI*ax;
		sf_result s;
		int stat_s = sf_sin_e(r, &s);
		result->val = s.val/r;
		result->err = s.err/r + 2.0 * MVDR_DBL_EPSILON * fabs(result->val);
		return stat_s;
    }
}

/* return |z| */
double complex_abs(double z_Re, double z_Lm){
	return hypot(z_Re, z_Lm);
}

/* z=a/b */
static void complex_div(double a_Re, double a_Lm, double b_Re, double b_Lm, double *z_Re, double *z_Lm){                             
	double s = 1.0 / complex_abs(b_Re,b_Lm);

	double sbr = s * b_Re;
	double sbi = s * b_Lm;

	*z_Re = (a_Re * sbr + a_Lm * sbi) * s;
	*z_Lm = (a_Lm * sbr - a_Re * sbi) * s;

}

void matrix_complex_identity(double (*data_Re)[N_CHANNEL], double (*data_Lm)[N_CHANNEL]){
	for (int i = 0; i < N_CHANNEL; i++){
		for (int j = 0; j < N_CHANNEL; j++){
			if(i == j){
				data_Re[i][j] = 1.0;
				data_Lm[i][j] = 0.0;
			}
			else{
				data_Re[i][j] = 0.0;
				data_Lm[i][j] = 0.0;
			}
		}
	}
}

//These functions compute the matrix-vector product and sum: y = \alpha A x + \beta y
static void complex_zgemv(double alpha_Re,double alpha_Lm,
						  double (*A_Re)[N_CHANNEL],double (*A_Lm)[N_CHANNEL],
						  double *X_Re,double *X_Lm,
						  double beta_Re,double beta_Lm,
						  double *Y_Re,double *Y_Lm){
	
	int lenX = N_CHANNEL;
	int lenY = N_CHANNEL;

	if ((alpha_Re == 0.0 && alpha_Lm == 0.0) && (beta_Re == 1.0 && beta_Lm == 0.0))
		return;

	 /* form  y := beta*y */
	if (beta_Re == 0.0 && beta_Lm == 0.0) {
		for (int i = 0; i < lenY; i++) {
			Y_Re[i] = 0.0;
			Y_Lm[i] = 0.0;
		}
  	}
	else if (!(beta_Re == 1.0 && beta_Lm == 0.0)){
		for (int i = 0; i < lenY; i++){
			Y_Re[i] = Y_Re[i] * beta_Re - Y_Lm[i] * beta_Lm;
			Y_Lm[i] = Y_Re[i] * beta_Lm + Y_Lm[i] * beta_Re;
		}
	}

	if (alpha_Re == 0.0 && alpha_Lm == 0.0)
    	return;

	/* form  y := alpha*A^H*x + y */
	for (int j = 0; j < lenX; j++){
		double tmpR = alpha_Re * X_Re[j] - alpha_Lm * X_Lm[j];
		double tmpI = alpha_Re * X_Lm[j] + alpha_Lm * X_Re[j];

		for (int i = 0; i < lenY; i++){
			Y_Re[i] += A_Re[j][i] * tmpR - (-A_Lm[j][i]) * tmpI;
			Y_Lm[i] += A_Re[j][i] * tmpI + (-A_Lm[j][i]) * tmpR;
		}
	}
}

/**********************************MVDR Beamformer*******************************/
typedef struct _BeamformerMVDR_{
	double wq_Re[M_][N_CHANNEL];
	double wq_Lm[M_][N_CHANNEL];

	double snapshots_Re[M_][N_CHANNEL];
	double snapshots_Lm[M_][N_CHANNEL];

	double R_Re[M_/2+1][N_CHANNEL][N_CHANNEL];
	double R_Lm[M_/2+1][N_CHANNEL][N_CHANNEL];

	double invR_Re[M_/2+1][N_CHANNEL][N_CHANNEL];
	double invR_Lm[M_/2+1][N_CHANNEL][N_CHANNEL];

	double wmvdr_Re[M_/2+1][N_CHANNEL];
	double wmvdr_Lm[M_/2+1][N_CHANNEL];

	float diagonal_weights[M_/2+1];


}*BeamformerMVDR;

BeamformerMVDR beamformerMVDR_created(){
	BeamformerMVDR beamformerMVDR = (BeamformerMVDR)malloc(sizeof(struct _BeamformerMVDR_));
	
	return beamformerMVDR;
}

void beamformerMVDR_calcMainlobe(BeamformerMVDR beamformerMVDR,float *delays){
	int fftLen2 = M_ / 2;
	double ReTmp,LmTmp;

    // calculate weights of a direct component.
	for(int chanX = 0; chanX < N_CHANNEL; chanX++){
		complex_polar(1.0,0.0,&ReTmp,&LmTmp);
		beamformerMVDR->wq_Re[0][chanX] = ReTmp/N_CHANNEL;
		beamformerMVDR->wq_Lm[0][chanX] = LmTmp/N_CHANNEL;
		//printf("%g:%g\n",beamformerMVDR->wq_Re[0][chanX],beamformerMVDR->wq_Lm[0][chanX]);
	}

	// calculate weights from FFT bin 1 to fftLen - 1 by using the property of the symmetry.
	for(int fbinX = 1; fbinX < fftLen2; fbinX++) {
		for(int chanX = 0; chanX < N_CHANNEL; chanX++) {
			double val = -2.0 * M_PI * fbinX * delays[chanX] * SAMPLERATE / M_;

			complex_polar(1.0,val,&ReTmp,&LmTmp);
			beamformerMVDR->wq_Re[fbinX][chanX] = ReTmp/N_CHANNEL;
			beamformerMVDR->wq_Lm[fbinX][chanX] = LmTmp/N_CHANNEL;
			//printf("%g,%g ",beamformerMVDR->wq_Re[fbinX][chanX],beamformerMVDR->wq_Lm[fbinX][chanX]);

			complex_polar(1.0,-val,&ReTmp,&LmTmp);
			beamformerMVDR->wq_Re[M_-fbinX][chanX] = ReTmp/N_CHANNEL;
			beamformerMVDR->wq_Lm[M_-fbinX][chanX] = LmTmp/N_CHANNEL;
			//printf("%g,%g ",beamformerMVDR->wq_Re[M_-fbinX][chanX],beamformerMVDR->wq_Lm[M_-fbinX][chanX]);
		}
	}

	// for exp(-j*pi)
	for(int chanX = 0; chanX < N_CHANNEL; chanX++) {
		double val = -M_PI * SAMPLERATE * delays[chanX];
		complex_polar(1.0,val,&ReTmp,&LmTmp);
		beamformerMVDR->wq_Re[fftLen2][chanX] = ReTmp/N_CHANNEL;
		beamformerMVDR->wq_Lm[fftLen2][chanX] = LmTmp/N_CHANNEL;
		//printf("%g,%g\n",beamformerMVDR->wq_Re[fftLen2][chanX],beamformerMVDR->wq_Lm[fftLen2][chanX]);
	}
}

void beamformerMVDR_set_diffuse_noise_model(BeamformerMVDR beamformerMVDR, const float(*micPositions)[3], float samplerate, float sspeed){
	double dm[N_CHANNEL][N_CHANNEL];

	// calculate the distance matrix.
	for(int m = 0; m < N_CHANNEL; m++){
		for(int n = 0; n < m; n++){
			double dx = micPositions[m][0] - micPositions[n][0];
			double dy = micPositions[m][1] - micPositions[n][1];
			double dz = micPositions[m][2] - micPositions[m][2];

			dm[m][n] = sqrt(dx * dx + dy * dy + dz * dz);
		}
	}

	for(int fbinX = 0; fbinX <= M_/2; fbinX++){
		double omega_d_c = 2.0 * samplerate * fbinX / (M_ * sspeed);

		for(int m = 0; m < N_CHANNEL; m++){
			for (int n = 0; n < m; n++){
				sf_result result;
				sf_sinc(omega_d_c * dm[m][n], &result);
				
				beamformerMVDR->R_Re[fbinX][m][n] = result.val; //Gamma_mn
				beamformerMVDR->R_Lm[fbinX][m][n] = 0.0;
			}
		}

		for(int m = 0; m < N_CHANNEL; m++){
			beamformerMVDR->R_Re[fbinX][m][m] = 1.0;
			beamformerMVDR->R_Lm[fbinX][m][m] = 0.0;
		}
		
		for(int m = 0; m < N_CHANNEL; m++){
			for(int n = m + 1; n < N_CHANNEL; n++){
				beamformerMVDR->R_Re[fbinX][m][n] = beamformerMVDR->R_Re[fbinX][n][m];
				beamformerMVDR->R_Lm[fbinX][m][n] = beamformerMVDR->R_Lm[fbinX][n][m];
			}
		}
	}

	#if 0
	for(int i = 0; i < M_/2+1; i++){
		for(int j = 0; j < N_CHANNEL; j++){
			for(int k = 0; k < N_CHANNEL; k++){
				printf("%g,%g ",beamformerMVDR->R_Re[i][j][k],beamformerMVDR->R_Lm[i][j][k]);
			}
		}
	}
	printf("\n");
	getchar();	
	#endif
}

void beamformerMVDR_divide_nondiagonal_elements(BeamformerMVDR beamformerMVDR,int fbinX, float mu){
	for (int chanX = 0; chanX < N_CHANNEL; chanX++){
		for (int chanY = 0; chanY < N_CHANNEL; chanY++){
			if (chanX != chanY){
				complex_div(beamformerMVDR->R_Re[fbinX][chanX][chanY], beamformerMVDR->R_Lm[fbinX][chanX][chanY], 
				1.0+mu, 0.0, &beamformerMVDR->R_Re[fbinX][chanX][chanY],&beamformerMVDR->R_Lm[fbinX][chanX][chanY]);
			}
		}
	}
}

void beamformerMVDR_divide_all_nondiagonal_elements(BeamformerMVDR beamformerMVDR,float mu){
	for(int fbinX = 0; fbinX <= M_/2; fbinX++){
		beamformerMVDR_divide_nondiagonal_elements(beamformerMVDR,fbinX, mu);
	}

	#if 0
	for(int i = 0; i < M_/2+1; i++){
		for(int j = 0; j < N_CHANNEL; j++){
			for(int k = 0; k < N_CHANNEL; k++){
				printf("%g,%g ",beamformerMVDR->R_Re[i][j][k],beamformerMVDR->R_Lm[i][j][k]);
			}
		}
	}
	printf("\n");
	getchar();	
	#endif
}

extern bool pseudoinverse(double (*A_Re)[N_CHANNEL], double (*A_Lm)[N_CHANNEL],
                  double (*invA_Re)[N_CHANNEL], double (*invA_Lm)[N_CHANNEL], float dThreshold);

void beamformerMVDR_calc_mvdr_weights(BeamformerMVDR beamformerMVDR,float samplerate, float dThreshold, int calcInverseMatrix){

	double tmpH_Re[N_CHANNEL];
	double tmpH_Lm[N_CHANNEL]; 

	double val1_Re = 1.0;
	double val1_Lm = 0.0;

	double val0_Re = 0.0;
	double val0_Lm = 0.0;

	double Lambda_Re;
	double Lambda_Lm;

	for(int chanX=0; chanX < N_CHANNEL; chanX++){
		beamformerMVDR->wmvdr_Re[0][chanX] = val1_Re;
		beamformerMVDR->wmvdr_Lm[0][chanX] = val1_Lm;
	}

	for(int fbinX=1; fbinX<=M_/2; fbinX++){
		double norm_Re;
		double norm_Lm;

		//calculate the inverse matrix of the coherence matrix
		if (1 == calcInverseMatrix){
			bool ret = pseudoinverse(beamformerMVDR->R_Re[fbinX],beamformerMVDR->R_Lm[fbinX], beamformerMVDR->invR_Re[fbinX],beamformerMVDR->invR_Lm[fbinX], dThreshold);
			#if 0
			for(int i = 0; i < N_CHANNEL; i++){
				for(int j = 0; j < N_CHANNEL; j++){
					printf("%g,%g ",beamformerMVDR->invR_Re[fbinX][i][j],beamformerMVDR->invR_Lm[fbinX][i][j]);
				}
			}
			printf("\n");
			//getchar();
			#endif
			if (false == ret){
				matrix_complex_identity(beamformerMVDR->invR_Re[fbinX],beamformerMVDR->invR_Lm[fbinX]);
			}
		}

		complex_zgemv(val1_Re, val1_Lm, beamformerMVDR->invR_Re[fbinX], beamformerMVDR->invR_Lm[fbinX], 
						beamformerMVDR->wq_Re[fbinX], beamformerMVDR->wq_Lm[fbinX], val0_Re, val0_Lm ,tmpH_Re, tmpH_Lm); // tmpH = invR^H * d

		#if 0
		for(int i = 0; i < N_CHANNEL; i++){
			printf("%f,%f ",tmpH_Re[i],tmpH_Lm[i]);
		}
		printf("\n");
		getchar();		
		#endif
		Lambda_Re = 0.0f;
		Lambda_Lm = 0.0f;
		for (int i = 0; i < N_CHANNEL; i++){
			Lambda_Re += tmpH_Re[i] * beamformerMVDR->wq_Re[fbinX][i] + tmpH_Lm[i] * beamformerMVDR->wq_Lm[fbinX][i];
			Lambda_Lm += tmpH_Re[i] * beamformerMVDR->wq_Lm[fbinX][i] - tmpH_Lm[i] * beamformerMVDR->wq_Re[fbinX][i];
		}

		norm_Re = Lambda_Re * N_CHANNEL;
		norm_Lm = Lambda_Lm * N_CHANNEL;

		for (int chanX = 0; chanX < N_CHANNEL; chanX++){
			complex_div(tmpH_Re[chanX],tmpH_Lm[chanX],norm_Re,norm_Lm,
						&beamformerMVDR->wmvdr_Re[fbinX][chanX],&beamformerMVDR->wmvdr_Lm[fbinX][chanX]);			
		}
	}

	#if 0
	for(int i = 0; i < 256/2+1; i++)
	for(int j = 0; j < 4; j++){
		printf("%f,%f ",beamformerMVDR->wmvdr_Re[i][j],beamformerMVDR->wmvdr_Lm[i][j]);
	}
	printf("\n");
	//getchar();
	#endif
}

void beamformerMVDR_process(BeamformerMVDR beamformerMVDR,double analysis_data[][M_*2],double *output_data){	
	for(int i = 0; i < M_; i++){
		for(int j = 0; j < N_CHANNEL; j++){
			beamformerMVDR->snapshots_Re[i][j] = analysis_data[j][2*i];
			beamformerMVDR->snapshots_Lm[i][j] = analysis_data[j][2*i+1];
		}
	}

	//原始数据
	#if 0
	for(int i = 0; i < M_; i++){
		printf("%g,%g ",analysis_data[0][2*i],analysis_data[0][2*i+1]);
		//printf("%g,%g ",analysis_data[1][2*i],analysis_data[1][2*i+1]);
		//printf("%g,%g ",analysis_data[2][2*i],analysis_data[2][2*i+1]);
		//printf("%g,%g ",analysis_data[3][2*i],analysis_data[3][2*i+1]);
	}
	printf("\n");
	getchar();
	#endif

	int fftLen2 = M_/2;
	//(a+bi)(c+di)=(ac-bd)+(bc+ad)i
	//conj (ac+bd) + (bc-ad)i
	output_data[0] = 0.0f;
	output_data[1] = 0.0f;

	// calculate a direct component.
	for(int i = 0; i < N_CHANNEL; i++){
		output_data[0] += beamformerMVDR->wmvdr_Re[0][i]*beamformerMVDR->snapshots_Re[0][i] + beamformerMVDR->wmvdr_Lm[0][i]*beamformerMVDR->snapshots_Lm[0][i];
		output_data[1] += beamformerMVDR->wmvdr_Re[0][i]*beamformerMVDR->snapshots_Lm[0][i] - beamformerMVDR->wmvdr_Lm[0][i]*beamformerMVDR->snapshots_Re[0][i];
		//printf("%g,%g\n",beamformerMVDR->wq_Re[0][i],beamformerMVDR->wq_Lm[0][i]);
		//printf("%g,%g\n",beamformerMVDR->snapshots_Re[0][i],beamformerMVDR->snapshots_Lm[0][i]);
	}
	//printf("O:%g,%g\n",output_data[0],output_data[1]);
	//getchar();

	// calculate outputs from bin 1 to fftLen - 1 by using the property of the symmetry.
	for (int fbinX = 1; fbinX <= fftLen2; fbinX++) {
		double val_Re = 0.0f;
		double val_Lm = 0.0f;
		for(int i = 0; i < N_CHANNEL; i++){
			val_Re += beamformerMVDR->wmvdr_Re[fbinX][i]*beamformerMVDR->snapshots_Re[fbinX][i] + beamformerMVDR->wmvdr_Lm[fbinX][i]*beamformerMVDR->snapshots_Lm[fbinX][i];
			val_Lm += beamformerMVDR->wmvdr_Re[fbinX][i]*beamformerMVDR->snapshots_Lm[fbinX][i] - beamformerMVDR->wmvdr_Lm[fbinX][i]*beamformerMVDR->snapshots_Re[fbinX][i];
		}

		if( fbinX < fftLen2 ){
			output_data[2*fbinX] = val_Re;
			output_data[2*fbinX+1] = val_Lm;
			
			output_data[2*(M_-fbinX)] = val_Re;
			output_data[2*(M_-fbinX)+1] = -val_Lm;
		}
		else{
			output_data[2*fftLen2] = val_Re;
			output_data[2*fftLen2+1] = val_Lm;
		}
	}

	// 数据
	#if 0
	for(int i = 0; i < M_; i++){
		printf("%g,%g ",beamformerMVDR->snapshots_Re[i][0],beamformerMVDR->snapshots_Lm[i][0]);
		//printf("%g,%g ",beamformerMVDR->snapshots_Re[i][1],beamformerMVDR->snapshots_Lm[i][1]);
		//printf("%g,%g ",beamformerMVDR->snapshots_Re[i][2],beamformerMVDR->snapshots_Lm[i][2]);
		//printf("%g,%g ",beamformerMVDR->snapshots_Re[i][3],beamformerMVDR->snapshots_Lm[i][3]);	
	}
	printf("\n");
	getchar();	
	#endif

	// beam输出
	#if 0
	for(int i = 0; i < M_; i++){
		printf("%g,%g ",output_data[2*i],output_data[2*i+1]);
	}
	printf("\n");
	//getchar();
	#endif
}


/**********************************ZelinskiPostFilter********************************/
typedef struct _ZelinskiPostFilter_{
	BeamformerMVDR beamformerMVDR;
	int count;
	double alpha;
	int type;
	double prevCSDs_Re[M_][N_CHANNEL*N_CHANNEL];
	double prevCSDs_Lm[M_][N_CHANNEL*N_CHANNEL];
	double wp1_Re[M_];
	double wp1_Lm[M_];

}*ZelinskiPostFilter;

typedef enum {
  TYPE_ZELINSKI1_REAL = 0x01,
  TYPE_ZELINSKI1_ABS  = 0x02,
  TYPE_APAB = 0x04,
  TYPE_ZELINSKI2 = 0x08,
  NO_USE_POST_FILTER = 0x00
} PostfilterType;

ZelinskiPostFilter zelinskiPostFilter_created(BeamformerMVDR beamformerMVDR,double alpha,int type){
	ZelinskiPostFilter zelinskiPostFilter = (ZelinskiPostFilter)malloc(sizeof(struct _ZelinskiPostFilter_));
	zelinskiPostFilter->beamformerMVDR = beamformerMVDR;
	zelinskiPostFilter->count = 0;
	zelinskiPostFilter->alpha = alpha;
	zelinskiPostFilter->type = type;

	for(int i = 0; i < M_; i++){
		for(int j = 0; j < N_CHANNEL*N_CHANNEL; j++){
			zelinskiPostFilter->prevCSDs_Re[i][j] = 0.0f;
			zelinskiPostFilter->prevCSDs_Lm[i][j] = 0.0f;
		}
	}

	return zelinskiPostFilter;
}

//this function compensates a time delay for the signal at each channel.
static void time_alignment_(double *wq_Re,double *wq_Lm,double *snapshots_Re,double *snapshots_Lm,
						    double *output_Re,double *output_Lm){

  double dsf_Re,dsf_Lm;
  double xsf_Re,xsf_Lm;
  double ysf_Re,ysf_Lm;

  for(int i = 0;i < N_CHANNEL; i++){
    dsf_Re = wq_Re[i];
	dsf_Lm = -wq_Lm[i];
	
	xsf_Re = snapshots_Re[i];
	xsf_Lm = snapshots_Lm[i];
	
	ysf_Re = dsf_Re*xsf_Re - dsf_Lm*xsf_Lm;
	ysf_Lm = dsf_Re*xsf_Lm + dsf_Lm*xsf_Re;

	output_Re[i] = ysf_Re;
	output_Lm[i] = ysf_Lm;
  }
}

//calculate cross spectral density (CSD).
static void calc_CSD_(double xi_Re, double xi_Lm, double xj_Re, double xj_Lm, 
					  double prevPhi_Re, double prevPhi_Lm, double *out_data_Re,double *out_data_Lm,double alpha){

	double xjA_Re = xj_Re;
	double xjA_Lm = -xj_Lm;														

	double xixjA_Re = xi_Re * xjA_Re - xi_Lm * xjA_Lm;
	double xixjA_Lm = xi_Re * xjA_Lm + xi_Lm * xjA_Re;

	if(alpha > 0.0){
		*out_data_Re = prevPhi_Re*alpha + xixjA_Re*(1.0-alpha);
		*out_data_Lm = prevPhi_Lm*alpha + xixjA_Lm*(1.0-alpha);
	}
	else{
		*out_data_Re = xixjA_Re;
		*out_data_Lm = xixjA_Lm;
	}
}

#define SPECTRAL_FLOOR 0.0001 // to avoid  very small transfer function components

static double ZelinskiFilter_f(double *wq_Re,double *wq_Lm,double *snapshots_Re,double *snapshots_Lm,
								double *prevCSDf_Re,double *prevCSDf_Lm,
							    double alpha, int pfType){

	double timeAlignedSignalf_Re[N_CHANNEL];
	double timeAlignedSignalf_Lm[N_CHANNEL];
	double sum_Re;
	double sum_Lm;

	time_alignment_(wq_Re,wq_Lm,snapshots_Re,snapshots_Lm,timeAlignedSignalf_Re,timeAlignedSignalf_Lm);
	
	#if 0
	for(int i = 0; i < N_CHANNEL; i++){
		printf("%g,%g ",timeAlignedSignalf_Re[i],timeAlignedSignalf_Lm[i]);
	}
	#endif

	sum_Re = 0.0f;
	sum_Lm = 0.0f;

	int idx;
	double xi_Re,xi_Lm;
	double xj_Re,xj_Lm;
	double prevPhi_Re,prevPhi_Lm;
	double estCSD_Re,estCSD_Lm;

	for(int i = 0;i < N_CHANNEL-1; i++){
		for(int j = i+1;j < N_CHANNEL;j++){
			idx = i * N_CHANNEL + j;
			xi_Re = timeAlignedSignalf_Re[i];
			xi_Lm = timeAlignedSignalf_Lm[i];
			xj_Re = timeAlignedSignalf_Re[j];
			xj_Lm = timeAlignedSignalf_Lm[j];
			prevPhi_Re = prevCSDf_Re[idx];
			prevPhi_Lm = prevCSDf_Lm[idx];
			calc_CSD_(xi_Re, xi_Lm, xj_Re, xj_Lm, prevPhi_Re, prevPhi_Lm, &estCSD_Re, &estCSD_Lm, alpha);
			#if 0
			printf("%g,%g ",estCSD_Re,estCSD_Lm);
			#endif
			sum_Re = sum_Re + estCSD_Re;
			sum_Lm = sum_Lm + estCSD_Lm;
			prevCSDf_Re[idx] = estCSD_Re;
			prevCSDf_Lm[idx] = estCSD_Lm;
		}
	}
	
	double numerator=0.0;

	// the real operation is performed 
    if(TYPE_ZELINSKI1_REAL & pfType){
		numerator = sum_Re;
		if( numerator < 0.0 )
			numerator = 0.0;
			//fprintf(stderr,"Zelinski post-filter  1 is used\n");
	}
	else{
		numerator = sqrt(sum_Re*sum_Re+sum_Lm*sum_Lm);
	}

	double estPSD, prevPSD;
	double denominator = 0.0f;
	double tmpC_Re;
	double tmpC_Lm;

	for(int i = 0; i < N_CHANNEL; i++){
		idx = i * N_CHANNEL + i;
		xi_Re = timeAlignedSignalf_Re[i];
		xi_Lm = timeAlignedSignalf_Lm[i];

		if( alpha > 0.0 ){
			prevPSD = prevCSDf_Re[idx];
			estPSD = alpha * prevPSD + (1.0-alpha) * (xi_Re*xi_Re + xi_Lm*xi_Lm);
		}
		else{
			estPSD = xi_Re*xi_Re + xi_Lm*xi_Lm;
		}
		denominator += estPSD;
		tmpC_Re = estPSD;
		tmpC_Lm = 0.0f;
		prevCSDf_Re[idx] = tmpC_Re;
		prevCSDf_Lm[idx] = tmpC_Lm;
	}

	double W1f = 1.0;
	W1f = (numerator / denominator) * (2.0 / (N_CHANNEL - 1.0));
    // to avoid artificial amplification,
    if(W1f >= 1.0) W1f =  1.0;
    if(W1f < SPECTRAL_FLOOR) W1f = SPECTRAL_FLOOR;

 	#if 0
	printf("%g\n",W1f);
	getchar();
	#endif

	return W1f;
}

void zelinskiPostFilter_process(ZelinskiPostFilter zelinskiPostFilter,double* beam_input_data,double *output_data){	
	double alpha;
	int type;
	double r;
	double wf_Re,wf_Lm;
	

	if(zelinskiPostFilter->count == 0){
		type = NO_USE_POST_FILTER;
	}
	else{
		type = zelinskiPostFilter->type;
	}

	if(zelinskiPostFilter->count > 1){
		alpha =  zelinskiPostFilter->alpha;
	}
	else{
		alpha = 0.0;
	}

	if(zelinskiPostFilter->count < 2){
		zelinskiPostFilter->count++;
	}

	int fftLen2 = M_/2;

	for(int fbinX=0;fbinX<=fftLen2;fbinX++){
		r = ZelinskiFilter_f(zelinskiPostFilter->beamformerMVDR->wq_Re[fbinX],zelinskiPostFilter->beamformerMVDR->wq_Lm[fbinX],
						 	 zelinskiPostFilter->beamformerMVDR->snapshots_Re[fbinX],zelinskiPostFilter->beamformerMVDR->snapshots_Lm[fbinX],
						 	 zelinskiPostFilter->prevCSDs_Re[fbinX],zelinskiPostFilter->prevCSDs_Lm[fbinX],
							 alpha,type);

		wf_Re = r * cos (0);
		wf_Lm = r * sin (0);
		
		//printf("%g,%g,%g ",r,wf_Re,wf_Lm);

		zelinskiPostFilter->wp1_Re[fbinX] = wf_Re;
		zelinskiPostFilter->wp1_Lm[fbinX] = wf_Lm;

		if( fbinX > 0 && fbinX < fftLen2 ){  //substitute a conjugate component
			zelinskiPostFilter->wp1_Re[M_ - fbinX] = wf_Re;
			zelinskiPostFilter->wp1_Lm[M_ - fbinX] = -wf_Lm;
		}
	}
	#if 0
	printf("\n");
    getchar();
	#endif

	memcpy(output_data,beam_input_data,512*sizeof(double));

	if(NO_USE_POST_FILTER == type){// CSD is just updated.
    	return;
	}
	
	double outf_Re,outf_Lm;
	

	for(int fbinX = 0;fbinX <= fftLen2; fbinX++){
		wf_Re = zelinskiPostFilter->wp1_Re[fbinX];
		wf_Lm = zelinskiPostFilter->wp1_Lm[fbinX];

		outf_Re = wf_Re*output_data[2*fbinX] - wf_Lm*output_data[2*fbinX+1];
		outf_Lm = wf_Re*output_data[2*fbinX+1] + wf_Lm*output_data[2*fbinX];

		//printf("%g,%g ",outf_Re,outf_Lm);

		output_data[2*fbinX]   = outf_Re;
		output_data[2*fbinX+1] = outf_Lm;

		if(fbinX > 0 && fbinX < fftLen2){// substitute a conjugate component
			output_data[2*(M_ - fbinX)] = outf_Re;
			output_data[2*(M_ - fbinX)+1] = -outf_Lm;
		}
	}

	#if 0
	printf("\n");
	getchar();
	#endif

	#if 0
	for(int i = 0; i < 256; i++){
		printf("%g,%g ",output_data[2*i],output_data[2*i+1]);
	}
	printf("\n");
	getchar();
	#endif
}

OverSampledDFTAnalysisBank analysisBank1;
OverSampledDFTAnalysisBank analysisBank2;
OverSampledDFTAnalysisBank analysisBank3;
OverSampledDFTAnalysisBank analysisBank4;
OverSampledDFTSynthesisBank synthesisBank;

BeamformerMVDR beamformerMVDR;
ZelinskiPostFilter zelinskiPostFilter;

void btk_beamforming_set_number(BeamformerMVDR beamformerMVDR,int theta_index,int phi_index){
	memcpy(&beamformerMVDR->wmvdr_Re[0][0],wmvdr_Re[3*theta_index+phi_index],sizeof(double)*65*4);
	memcpy(&beamformerMVDR->wmvdr_Lm[0][0],wmvdr_Lm[3*theta_index+phi_index],sizeof(double)*65*4);
	//memcpy(&beamformerMVDR->wmvdr_Re[0][0],wmvdr_Re[0],sizeof(double)*65*4);
	//memcpy(&beamformerMVDR->wmvdr_Lm[0][0],wmvdr_Lm[0],sizeof(double)*65*4);
}

void btk_beamforming_set_location(short theta,short phi){
	double theta_double = M_PI * theta / 180;
	double phi_dobule = M_PI * (90-phi) /180; 
	float delays[4];
	calcDelaysPolar2(theta_double, phi_dobule, micPos, delays, N_CHANNEL);
	beamformerMVDR_calcMainlobe(beamformerMVDR,delays);
	beamformerMVDR_set_diffuse_noise_model(beamformerMVDR,micPos,SAMPLERATE,SOUNDSPEED);
	beamformerMVDR_divide_all_nondiagonal_elements(beamformerMVDR,0.1);
	beamformerMVDR_calc_mvdr_weights(beamformerMVDR,SAMPLERATE,1.0E-8,1);	
}

int btk_beamforming_process(short input_channels[][D_],short *out){

	double analysis_output[N_CHANNEL][M_*2];
	double beam_output[M_*2];
	double pf_output[M_*2];

	overSampled_dtf_analysisBank_process(analysisBank1,input_channels[0],analysis_output[0]);
	overSampled_dtf_analysisBank_process(analysisBank2,input_channels[1],analysis_output[1]);
	overSampled_dtf_analysisBank_process(analysisBank3,input_channels[2],analysis_output[2]);
	overSampled_dtf_analysisBank_process(analysisBank4,input_channels[3],analysis_output[3]);

	beamformerMVDR_process(beamformerMVDR,analysis_output,beam_output);
	//zelinskiPostFilter_process(zelinskiPostFilter,beam_output,pf_output);

	int ret = overSampled_dtf_synthesisBank_process(synthesisBank,beam_output,out);

	return ret;
}	

int btk_beamforming_process_doa(short input_channels[][D_],short *out,int theta_index,int phi_index){

	double analysis_output[N_CHANNEL][M_*2];
	double beam_output[M_*2];
	double pf_output[M_*2];
	overSampled_dtf_analysisBank_process(analysisBank1,input_channels[0],analysis_output[0]);
	overSampled_dtf_analysisBank_process(analysisBank2,input_channels[1],analysis_output[1]);
	overSampled_dtf_analysisBank_process(analysisBank3,input_channels[2],analysis_output[2]);
	overSampled_dtf_analysisBank_process(analysisBank4,input_channels[3],analysis_output[3]);

	btk_beamforming_set_number(beamformerMVDR,theta_index,phi_index);
	beamformerMVDR_process(beamformerMVDR,analysis_output,beam_output);
	zelinskiPostFilter_process(zelinskiPostFilter,beam_output,pf_output);

	int ret = overSampled_dtf_synthesisBank_process(synthesisBank,pf_output,out);

	return ret;
}	

void btk_beamforming_init(){
	analysisBank1  = overSampled_dtf_analysisBank_created();
	analysisBank2  = overSampled_dtf_analysisBank_created();
	analysisBank3  = overSampled_dtf_analysisBank_created();
	analysisBank4  = overSampled_dtf_analysisBank_created();

	synthesisBank = overSampled_dtf_synthesisBank_created();

	beamformerMVDR = beamformerMVDR_created();
	
	zelinskiPostFilter = zelinskiPostFilter_created(beamformerMVDR,0.7,2);

	short zero_data[N_CHANNEL][D_];
	short beam_out[D_];
	memset(&zero_data[0][0],0x00,sizeof(short)*D_*N_CHANNEL);
	for(int i = 0; i < 7; i++){
		btk_beamforming_process(zero_data,beam_out);
	}

	//btk_beamforming_set_number(0);
	btk_beamforming_set_location(0,45);
}


#if 0
int main(){
	FILE *input1 = fopen("input1.wav","rb");
	FILE *input2 = fopen("input2.wav","rb");
	FILE *input3 = fopen("input3.wav","rb");
	FILE *input4 = fopen("input4.wav","rb");

	FILE *output = fopen("output.pcm","wb");

	OverSampledDFTAnalysisBank analysisBank1  = overSampled_dtf_analysisBank_created();
	OverSampledDFTAnalysisBank analysisBank2  = overSampled_dtf_analysisBank_created();
	OverSampledDFTAnalysisBank analysisBank3  = overSampled_dtf_analysisBank_created();
	OverSampledDFTAnalysisBank analysisBank4  = overSampled_dtf_analysisBank_created();

	OverSampledDFTSynthesisBank synthesisBank = overSampled_dtf_synthesisBank_created();

	BeamformerMVDR  beamformerMVDR = beamformerMVDR_created();

	ZelinskiPostFilter zelinskiPostFilter = zelinskiPostFilter_created(beamformerMVDR,0.3,2);

	float delays[4];
	calcDelaysPolar2(0.0f, 3.1415926/4, micPos, delays, N_CHANNEL);
	beamformerMVDR_calcMainlobe(beamformerMVDR,delays);
	beamformerMVDR_set_diffuse_noise_model(beamformerMVDR,micPos,SAMPLERATE,SOUNDSPEED);
	beamformerMVDR_divide_all_nondiagonal_elements(beamformerMVDR,0.1);
	beamformerMVDR_calc_mvdr_weights(beamformerMVDR,SAMPLERATE,1.0E-8,1);

	short audio_data1[D_];
	short audio_data2[D_];
	short audio_data3[D_];
	short audio_data4[D_];

	int ret;
	fseek(input1, 44, SEEK_SET);
	fseek(input2, 44, SEEK_SET);
	fseek(input3, 44, SEEK_SET);
	fseek(input4, 44, SEEK_SET);
		
	double analysis_output[N_CHANNEL][M_*2];
	double beam_output[M_*2];
	double pf_output[M_*2];
	short  final_output[D_];

	int count = 0;

	while(1){
		ret = fread(audio_data1,D_*2,1,input1);
		if(ret <= 0){
			return -1;
		}

		ret = fread(audio_data2,D_*2,1,input2);
		ret = fread(audio_data3,D_*2,1,input3);
		ret = fread(audio_data4,D_*2,1,input4);

		overSampled_dtf_analysisBank_process(analysisBank1,audio_data1,analysis_output[0]);
		overSampled_dtf_analysisBank_process(analysisBank2,audio_data2,analysis_output[1]);
		overSampled_dtf_analysisBank_process(analysisBank3,audio_data3,analysis_output[2]);
		overSampled_dtf_analysisBank_process(analysisBank4,audio_data4,analysis_output[3]);

		beamformerMVDR_process(beamformerMVDR,analysis_output,beam_output);

		zelinskiPostFilter_process(zelinskiPostFilter,beam_output,pf_output);

		ret = overSampled_dtf_synthesisBank_process(synthesisBank,pf_output,final_output);
		if(ret == 1){
			fwrite(final_output,D_*2,1,output);
		}
	}
}
#endif