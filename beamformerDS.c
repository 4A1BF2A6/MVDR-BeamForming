#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "beamformBTK.h"
#include "M128-m4-r1.h"
#include "fft_radix.h"

//#define MIC_ORG

//#define MIC_C980_PRO

#define MIC_JUPITER

// M_ : 256 number of subband
// m_ : 4 filter length
// R_ : 2
#define M_ 128
#define m_ 4
#define R_ 2
#define D_ BTK_D

#define N_CHANNEL 4
#define NC 1

#define M_PI 3.14159265358979323846
#define SOUNDSPEED 343740.0
#define SAMPLERATE 16000

typedef struct real_buffer_float{
	float data_float[m_*R_][M_];
	int zeros;

} *RealBuffer_Float;


#define MAX_VAL 32767.0
#define MIN_VAL -32768.0

#define FLOAT2SHORT_SAT(x) ((x>MAX_VAL)?MAX_VAL:(x<MIN_VAL)?MIN_VAL:x)

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

int overSampled_dtf_analysisBank_pre_process(OverSampledDFTAnalysisBank analysisBank,short *audio_data){
	#if 0
	for(int i = 0; i < D_; i++){
		printf("%d ",audio_data[i]);
	}
	printf("\n");
	getchar();
	#endif

	for(int i = 0; i < D_; i++){
		analysisBank->audio_data_overlap_float[i+D_] = (float)audio_data[i];
	}

	real_buffer_float_insert_reverse(analysisBank->analysis_buf,analysisBank->audio_data_overlap_float);
	#if 0
	for(int i = 0; i < 256; i++){
		printf("%g ",analysisBank->audio_data_overlap_float[i]);
	}
	printf("\n");
	getchar();
	#endif

	for(int i = 0; i < D_; i++){
		analysisBank->audio_data_overlap_float[i] = analysisBank->audio_data_overlap_float[i+D_];
	}
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
				synthesis_output[D_-d -1] = FLOAT2SHORT_SAT(synthesis_output[D_-d-1] + synthesisBank->synthesis_overlap[R_-sampX-1][d+sampX*D_]);
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

#ifdef MIC_ORG
float micPos[4][3] = {
	{-32.25,0.0,0.0},
	{0.0,-32.25,0.0},
	{32.25,0.0,0.0},
	{0.0,32.25,0.0}
};
#elif defined(MIC_C980_PRO)
float micPos[4][3] = {
	{-33, 0, 0.000},
	{-11, 0, 0.000},
	{ 11, 0, 0.000},
	{ 33, 0, 0.000},
};
#elif defined(MIC_JUPITER)
float micPos[4][3] = {
	{-24, 0, 0.000},
	{-8, 0, 0.000},
	{ 8, 0, 0.000},
	{ 24, 0, 0.000},
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

/**********************************D&S Beamformer*******************************/
typedef struct _BeamformerDS_{
	double wq_Re[M_][N_CHANNEL];
	double wq_Lm[M_][N_CHANNEL];

	double snapshots_Re[M_][N_CHANNEL];
	double snapshots_Lm[M_][N_CHANNEL];

}*BeamformerDS;

BeamformerDS beamformerDS_created(){
	BeamformerDS beamformerDS = (BeamformerDS)malloc(sizeof(struct _BeamformerDS_));
	
	return beamformerDS;
}

void beamformerDS_calcMainlobe(BeamformerDS beamformerDS,float *delays){
	int fftLen2 = M_ / 2;
	double ReTmp,LmTmp;

    // calculate weights of a direct component.
	for(int chanX = 0; chanX < N_CHANNEL; chanX++){
		complex_polar(1.0,0.0,&ReTmp,&LmTmp);
		beamformerDS->wq_Re[0][chanX] = ReTmp/N_CHANNEL;
		beamformerDS->wq_Lm[0][chanX] = LmTmp/N_CHANNEL;
		//printf("%g:%g\n",beamformerDS->wq_Re[0][chanX],beamformerDS->wq_Lm[0][chanX]);
	}

	// calculate weights from FFT bin 1 to fftLen - 1 by using the property of the symmetry.
	for(int fbinX = 1; fbinX < fftLen2; fbinX++) {
		for(int chanX = 0; chanX < N_CHANNEL; chanX++) {
			double val = -2.0 * M_PI * fbinX * delays[chanX] * SAMPLERATE / M_;

			complex_polar(1.0,val,&ReTmp,&LmTmp);
			beamformerDS->wq_Re[fbinX][chanX] = ReTmp/N_CHANNEL;
			beamformerDS->wq_Lm[fbinX][chanX] = LmTmp/N_CHANNEL;
			//printf("%g,%g ",beamformerDS->wq_Re[fbinX][chanX],beamformerDS->wq_Lm[fbinX][chanX]);

			complex_polar(1.0,-val,&ReTmp,&LmTmp);
			beamformerDS->wq_Re[M_-fbinX][chanX] = ReTmp/N_CHANNEL;
			beamformerDS->wq_Lm[M_-fbinX][chanX] = LmTmp/N_CHANNEL;
			//printf("%g,%g ",beamformerDS->wq_Re[M_-fbinX][chanX],beamformerDS->wq_Lm[M_-fbinX][chanX]);
		}
	}

	// for exp(-j*pi)
	for(int chanX = 0; chanX < N_CHANNEL; chanX++) {
		double val = -M_PI * SAMPLERATE * delays[chanX];
		complex_polar(1.0,val,&ReTmp,&LmTmp);
		beamformerDS->wq_Re[fftLen2][chanX] = ReTmp/N_CHANNEL;
		beamformerDS->wq_Lm[fftLen2][chanX] = LmTmp/N_CHANNEL;	
		//printf("%g,%g\n",beamformerDS->wq_Re[fftLen2][chanX],beamformerDS->wq_Lm[fftLen2][chanX]);
	}
}

void beamformerDS_process(BeamformerDS beamformerDS,double analysis_data[][M_*2],double *output_data){	
	for(int i = 0; i < M_; i++){
		for(int j = 0; j < N_CHANNEL; j++){
			beamformerDS->snapshots_Re[i][j] = analysis_data[j][2*i];
			beamformerDS->snapshots_Lm[i][j] = analysis_data[j][2*i+1];
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
		output_data[0] += beamformerDS->wq_Re[0][i]*beamformerDS->snapshots_Re[0][i] + beamformerDS->wq_Lm[0][i]*beamformerDS->snapshots_Lm[0][i];
		output_data[1] += beamformerDS->wq_Re[0][i]*beamformerDS->snapshots_Lm[0][i] - beamformerDS->wq_Lm[0][i]*beamformerDS->snapshots_Re[0][i];
		//printf("%g,%g\n",beamformerDS->wq_Re[0][i],beamformerDS->wq_Lm[0][i]);
		//printf("%g,%g\n",beamformerDS->snapshots_Re[0][i],beamformerDS->snapshots_Lm[0][i]);
	}
	//printf("O:%g,%g\n",output_data[0],output_data[1]);
	//getchar();

	// calculate outputs from bin 1 to fftLen - 1 by using the property of the symmetry.
	for (int fbinX = 1; fbinX <= fftLen2; fbinX++) {
		double val_Re = 0.0f;
		double val_Lm = 0.0f;
		for(int i = 0; i < N_CHANNEL; i++){
			val_Re += beamformerDS->wq_Re[fbinX][i]*beamformerDS->snapshots_Re[fbinX][i] + beamformerDS->wq_Lm[fbinX][i]*beamformerDS->snapshots_Lm[fbinX][i];
			val_Lm += beamformerDS->wq_Re[fbinX][i]*beamformerDS->snapshots_Lm[fbinX][i] - beamformerDS->wq_Lm[fbinX][i]*beamformerDS->snapshots_Re[fbinX][i];
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

	// 滤波器系数
	#if 0
	for(int i = 0; i < M_; i++){
		printf("%g,%g ",beamformerDS->wq_Re[i][0],beamformerDS->wq_Lm[i][0]);
		printf("%g,%g ",beamformerDS->wq_Re[i][1],beamformerDS->wq_Lm[i][1]);
		printf("%g,%g ",beamformerDS->wq_Re[i][2],beamformerDS->wq_Lm[i][2]);
		printf("%g,%g ",beamformerDS->wq_Re[i][3],beamformerDS->wq_Lm[i][3]);
	}
	printf("\n");
	getchar();
	#endif

	// 数据
	#if 0
	for(int i = 0; i < M_; i++){
		printf("%g,%g ",beamformerDS->snapshots_Re[i][0],beamformerDS->snapshots_Lm[i][0]);
		//printf("%g,%g ",beamformerDS->snapshots_Re[i][1],beamformerDS->snapshots_Lm[i][1]);
		//printf("%g,%g ",beamformerDS->snapshots_Re[i][2],beamformerDS->snapshots_Lm[i][2]);
		//printf("%g,%g ",beamformerDS->snapshots_Re[i][3],beamformerDS->snapshots_Lm[i][3]);	
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
	getchar();
	#endif
}

/**********************************ZelinskiPostFilter********************************/
typedef struct _ZelinskiPostFilter_{
	BeamformerDS beamformerDS;
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

ZelinskiPostFilter zelinskiPostFilter_created(BeamformerDS beamformerDS,double alpha,int type){
	ZelinskiPostFilter zelinskiPostFilter = (ZelinskiPostFilter)malloc(sizeof(struct _ZelinskiPostFilter_));
	zelinskiPostFilter->beamformerDS = beamformerDS;
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
		r = ZelinskiFilter_f(zelinskiPostFilter->beamformerDS->wq_Re[fbinX],zelinskiPostFilter->beamformerDS->wq_Lm[fbinX],
						 	zelinskiPostFilter->beamformerDS->snapshots_Re[fbinX],zelinskiPostFilter->beamformerDS->snapshots_Lm[fbinX],
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

BeamformerDS beamformerDS;
ZelinskiPostFilter zelinskiPostFilter;


void btk_beamforming_set_location(short theta,short phi){
	double theta_double = M_PI * theta / 180;
	double phi_dobule = M_PI * phi /180; 
	float delays[4];
	calcDelaysPolar2(theta_double, phi_dobule, micPos, delays, N_CHANNEL);
	beamformerDS_calcMainlobe(beamformerDS,delays);
}

void btk_beamforming_set_location_num(int doa_num){
	if(doa_num == 0){
		return;
	}
	double theta_double = M_PI * (doa_num-1)*36 / 180;
	double phi_dobule = M_PI * 45 /180;
	float delays[4];
	calcDelaysPolar2(theta_double, phi_dobule, micPos, delays, N_CHANNEL);
	beamformerDS_calcMainlobe(beamformerDS,delays);
}

int btk_beamforming_process(short input_channels[][D_],short *out){

	double analysis_output[N_CHANNEL][M_*2];
	double beam_output[M_*2];
	double pf_output[M_*2];

	overSampled_dtf_analysisBank_process(analysisBank1,input_channels[0],analysis_output[0]);
	overSampled_dtf_analysisBank_process(analysisBank2,input_channels[1],analysis_output[1]);
	overSampled_dtf_analysisBank_process(analysisBank3,input_channels[2],analysis_output[2]);
	overSampled_dtf_analysisBank_process(analysisBank4,input_channels[3],analysis_output[3]);

	beamformerDS_process(beamformerDS,analysis_output,beam_output);
	//zelinskiPostFilter_process(zelinskiPostFilter,beam_output,pf_output);

	int ret = overSampled_dtf_synthesisBank_process(synthesisBank,beam_output,out);

	return ret;
}

void btk_beamforming_init(){
	analysisBank1  = overSampled_dtf_analysisBank_created();
	analysisBank2  = overSampled_dtf_analysisBank_created();
	analysisBank3  = overSampled_dtf_analysisBank_created();
	analysisBank4  = overSampled_dtf_analysisBank_created();

	synthesisBank = overSampled_dtf_synthesisBank_created();

	beamformerDS = beamformerDS_created();
	
	zelinskiPostFilter = zelinskiPostFilter_created(beamformerDS,0.7,2);

	short zero_data[N_CHANNEL][D_];
	short beam_out[D_];
	memset(&zero_data[0][0],0x00,sizeof(short)*D_*N_CHANNEL);
	for(int i = 0; i < 7; i++){
		btk_beamforming_process(zero_data,beam_out);
	}

	btk_beamforming_set_location(90,45);
}


#if 0
void main(){
	FILE *input1 = fopen("input1.pcm","rb");
	FILE *input2 = fopen("input2.pcm","rb");
	FILE *input3 = fopen("input3.pcm","rb");
	FILE *input4 = fopen("input4.pcm","rb");

	FILE *output = fopen("output.pcm","wb");

	//printf("%s,%d %p %p %p %p\n",__FUNCTION__,__LINE__,input1,input2,input3,input4);
#if 0
	OverSampledDFTAnalysisBank analysisBank1  = overSampled_dtf_analysisBank_created();
	OverSampledDFTAnalysisBank analysisBank2  = overSampled_dtf_analysisBank_created();
	OverSampledDFTAnalysisBank analysisBank3  = overSampled_dtf_analysisBank_created();
	OverSampledDFTAnalysisBank analysisBank4  = overSampled_dtf_analysisBank_created();

	OverSampledDFTSynthesisBank synthesisBank = overSampled_dtf_synthesisBank_created();

	BeamformerDS  beamformerDS = beamformerDS_created();

	ZelinskiPostFilter zelinskiPostFilter = zelinskiPostFilter_created(beamformerDS,0.7,2);

	float delays[4];
	calcDelaysPolar2(0.0f, 3.1415926/4, micPos, delays, N_CHANNEL);
	beamformerDS_calcMainlobe(beamformerDS,delays);
#else
	btk_beamforming_init();
#endif

	int ret;
	// fseek(input1, 44, SEEK_SET);
	// fseek(input2, 44, SEEK_SET);
	// fseek(input3, 44, SEEK_SET);
	// fseek(input4, 44, SEEK_SET);

	double analysis_output[4][512];
	double beam_output[512];
	double pf_output[512];
	short  final_output[64];
	short beamforming_in[4][64] = {{0}};

	int count = 0;

	while(1){
		ret = fread(beamforming_in[0],64*2,1,input1);
		if(ret <= 0){
			return;
		}
		fread(beamforming_in[1],64*2,1,input2);
		fread(beamforming_in[2],64*2,1,input3);
		fread(beamforming_in[3],64*2,1,input4);
#if 0
		overSampled_dtf_analysisBank_process(analysisBank1,beamforming_in[0],analysis_output[0]);
		overSampled_dtf_analysisBank_process(analysisBank2,beamforming_in[1],analysis_output[1]);
		overSampled_dtf_analysisBank_process(analysisBank3,beamforming_in[2],analysis_output[2]);
		overSampled_dtf_analysisBank_process(analysisBank4,beamforming_in[3],analysis_output[3]);

		beamformerDS_process(beamformerDS,analysis_output,beam_output);

		//zelinskiPostFilter_process(zelinskiPostFilter,beam_output,pf_output);

		ret = overSampled_dtf_synthesisBank_process(synthesisBank,beam_output,final_output);
#else
		btk_beamforming_process(beamforming_in, final_output);
#endif
		if(ret == 1){
			fwrite(final_output,64*2,1,output);
		}
	}
}
#endif