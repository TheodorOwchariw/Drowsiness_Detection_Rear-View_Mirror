/**
    @file: dsp_lib.c

  	@brief: This library includes basic functions like FFT and windowing functions etc.
*/

/* ===========================================================================
** Copyright (C) 2018-2019 Infineon Technologies AG
** All rights reserved.
** ===========================================================================
**
** ===========================================================================
** This document contains proprietary information of Infineon Technologies AG.
** Passing on and copying of this document, and communication of its contents
** is not permitted without Infineon's prior written authorization.
** ===========================================================================
*/

/*
==============================================================================
   1. INCLUDE FILES
==============================================================================
 */

#include "dsp_lib.h"

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
 */

#define FAST_TIME_FFT_WIN_LEN			(256U)		/**< Max Window length for FFT windowing over fast time */

#define SLOW_TIME_FFT_WIN_LEN			(16U)		/**< Max Window length for FFT windowing over slow time */

/*
==============================================================================
   3. DATA
==============================================================================
 */

const arm_cfft_instance_f32 *fast_time_twiddle_table;
const arm_cfft_instance_f32 *slow_time_twiddle_table;

// BlackmanHarris Window without 3dB compensation
static float FastTime_Window[] = {
		0.000060,0.000069,0.000095,0.000138,0.000201,0.000283,0.000386,0.000511,0.000662,0.000839,0.001046,0.001285,0.001560,0.001874,0.002230,0.002634,
		0.003089,0.003600,0.004172,0.004811,0.005521,0.006309,0.007180,0.008141,0.009199,0.010361,0.011633,0.013022,0.014537,0.016185,0.017974,0.019912,
		0.022007,0.024267,0.026702,0.029319,0.032127,0.035136,0.038353,0.041787,0.045448,0.049343,0.053482,0.057872,0.062523,0.067441,0.072636,0.078115,
		0.083885,0.089953,0.096327,0.103011,0.110014,0.117340,0.124993,0.132979,0.141302,0.149965,0.158971,0.168321,0.178018,0.188062,0.198453,0.209191,
		0.220273,0.231697,0.243461,0.255559,0.267988,0.280741,0.293812,0.307192,0.320874,0.334847,0.349102,0.363627,0.378409,0.393436,0.408693,0.424165,
		0.439837,0.455692,0.471712,0.487879,0.504173,0.520575,0.537064,0.553619,0.570218,0.586838,0.603456,0.620049,0.636593,0.653063,0.669435,0.685683,
		0.701782,0.717708,0.733434,0.748935,0.764185,0.779160,0.793834,0.808181,0.822177,0.835799,0.849021,0.861821,0.874175,0.886062,0.897460,0.908348,
		0.918706,0.928515,0.937756,0.946413,0.954469,0.961909,0.968718,0.974885,0.980396,0.985242,0.989413,0.992902,0.995700,0.997804,0.999209,0.999912,
		0.999912,0.999209,0.997804,0.995700,0.992902,0.989413,0.985242,0.980396,0.974885,0.968718,0.961909,0.954469,0.946413,0.937756,0.928515,0.918706,
		0.908348,0.897460,0.886062,0.874175,0.861821,0.849021,0.835799,0.822177,0.808181,0.793834,0.779160,0.764185,0.748935,0.733434,0.717708,0.701782,
		0.685683,0.669435,0.653063,0.636593,0.620049,0.603456,0.586838,0.570218,0.553619,0.537064,0.520575,0.504173,0.487879,0.471712,0.455692,0.439837,
		0.424165,0.408693,0.393436,0.378409,0.363627,0.349102,0.334847,0.320874,0.307192,0.293812,0.280741,0.267988,0.255559,0.243461,0.231697,0.220273,
		0.209191,0.198453,0.188062,0.178018,0.168321,0.158971,0.149965,0.141302,0.132979,0.124993,0.117340,0.110014,0.103011,0.096327,0.089953,0.083885,
		0.078115,0.072636,0.067441,0.062523,0.057872,0.053482,0.049343,0.045448,0.041787,0.038353,0.035136,0.032127,0.029319,0.026702,0.024267,0.022007,
		0.019912,0.017974,0.016185,0.014537,0.013022,0.011633,0.010361,0.009199,0.008141,0.007180,0.006309,0.005521,0.004811,0.004172,0.003600,0.003089,
		0.002634,0.002230,0.001874,0.001560,0.001285,0.001046,0.000839,0.000662,0.000511,0.000386,0.000283,0.000201,0.000138,0.000095,0.000069,0.000060 };

// Chebychev Window with 60dB gain
static float SlowTime_Window[] = {
		0.024082,0.078911,0.182342,0.337575,0.531631,0.733976,0.903277,1.000000,1.000000,0.903277,0.733976,0.531631,0.337575,0.182342,0.078911,0.024082 };

/*
==============================================================================
   4. EXPORTED FUNCTIONS
==============================================================================
 */

void fft_init(void)
{
	/* Fast Time FFT */
#if (RANGE_FFT_SIZE == 1024U)
	fast_time_twiddle_table = &arm_cfft_sR_f32_len1024;		// as FFT_SIZE = 1024 points
#elif (RANGE_FFT_SIZE == 512U)
	fast_time_twiddle_table = &arm_cfft_sR_f32_len512;		// as FFT_SIZE = 512 points
#elif (RANGE_FFT_SIZE == 256U)
	fast_time_twiddle_table = &arm_cfft_sR_f32_len256;		// as FFT_SIZE = 256 points
#elif (RANGE_FFT_SIZE == 128U)
	fast_time_twiddle_table = &arm_cfft_sR_f32_len128;		// as FFT_SIZE = 128 points
#else
	fast_time_twiddle_table = &arm_cfft_sR_f32_len64;		// as FFT_SIZE = 64 points
#endif

	/* Slow Time FFT */
#if (DOPPLER_FFT_SIZE == 64U)
	slow_time_twiddle_table = &arm_cfft_sR_f32_len64;		// as FFT_SIZE = 64 points
#else
	slow_time_twiddle_table = &arm_cfft_sR_f32_len32;		// as FFT_SIZE = 32 points
#endif
}

//============================================================================

void compute_fft_signal(float* i_data, float* q_data, uint16_t Nd, uint16_t Nf, float if_scale,
						FFT_Input_t fft_type, FFT_Direction_t fft_direction,
						float* i_mean, float* q_mean, float* complex_fft_signal)
{
	uint32_t idx;
	uint32_t dec_idx;   // decimation index of FFT window in case selected data points are less than Window length
	float *fft_window;
	const arm_cfft_instance_f32 *cfft_twiddle_table;


	if (fft_direction == FFT_FAST_TIME)
	{
		fft_window = FastTime_Window;

		cfft_twiddle_table = fast_time_twiddle_table;

		dec_idx = FAST_TIME_FFT_WIN_LEN / Nd;
	}
	else /* fft_direction == FFT_SLOW_TIME */
	{
		fft_window = SlowTime_Window;

		cfft_twiddle_table = slow_time_twiddle_table;

		dec_idx = SLOW_TIME_FFT_WIN_LEN / Nd;
	}

	/* Find the mean in i_data */
	arm_mean_f32(i_data, Nd, i_mean);

	/* Find the mean in q_data */
	arm_mean_f32(q_data, Nd, q_mean);

	if (fft_type == FFT_INPUT_REAL_I)
	{
		/* Interleaved (re = I & im = 0) samples as {re[0], im[0], re[1], im[1], ...} */
		for (idx = 0; idx < Nd; idx++)
		{
			complex_fft_signal[2 * idx + 0] = 2.0 * (float)(i_data[idx] - *i_mean) * if_scale * fft_window[idx*dec_idx];	// additional scaling by 2 for real input FFT
			complex_fft_signal[2 * idx + 1] = 0;
		}
	}
	else if (fft_type == FFT_INPUT_REAL_Q)
	{
		/* Interleaved (re = Q & im = 0) samples as {re[0], im[0], re[1], im[1], ...} */
		for (idx = 0; idx < Nd; idx++)
		{
			complex_fft_signal[2 * idx + 0] = 2.0 * (float)(q_data[idx] - *q_mean) * if_scale * fft_window[idx*dec_idx];	// additional scaling by 2 for real input FFT
			complex_fft_signal[2 * idx + 1] = 0;
		}
	}
	else
	{
		/* Interleaved (re = I & im = Q) samples as {re[0], im[0], re[1], im[1], ...} */
		for (idx = 0; idx < Nd; idx++)
		{
			complex_fft_signal[2 * idx + 0] = (float)(i_data[idx] - *i_mean) * if_scale * fft_window[idx*dec_idx];
			complex_fft_signal[2 * idx + 1] = (float)(q_data[idx] - *q_mean) * if_scale * fft_window[idx*dec_idx];
		}
	}

	/* Zero Padding */
	for (uint16_t idx = Nd; idx < Nf; idx++)
	{
		complex_fft_signal[2 * idx + 0] = 0;
		complex_fft_signal[2 * idx + 1] = 0;
	}

	/* Processing the floating-point complex FFT. */
	arm_cfft_f32(cfft_twiddle_table, complex_fft_signal, 0, 1);
}

//============================================================================

void compute_fft_spectrum(float* fft_input_signal, uint32_t Nf, float* fft_output_spectrum)
{
	/* Convert to real magnitude data */
	arm_cmplx_mag_f32(fft_input_signal, fft_output_spectrum, Nf);
}
//============================================================================

target_angle_data compute_angle(float if1_i, float if1_q, float if2_i, float if2_q, double d_old, int16_t angle_offset_deg, float wave_length_ant_spacing_ratio)
{

	//float wave_length_ant_spacing_ratio = (TGT_WAVE_LENGTH_MM / TGT_ANTENNA_SPACING_MM);

	target_angle_data temp;

	float rx1_ang, rx2_ang;

	float d_phi;

	float delta_angle;

	float target_angle;

	rx1_ang = get_phase(if1_i, if1_q); //- (double)(0.13 * PI);

	rx2_ang = get_phase(if2_i, if2_q);

	d_phi = (rx1_ang - rx2_ang);

	if (d_phi <= 0)
	{
		d_phi += 2*PI;
	}
	d_phi -= PI;

	if ((uint32_t)d_old == IGNORE_NAN)
	{
		target_angle = 0;
	}
	else if (d_phi > d_old + 0.9* PI || d_phi < d_old - 0.9* PI)
	{
	   d_phi = d_old;
	}

	/* Arcus sinus (-PI/2 to PI/2), input= -1..1 */
	target_angle = asin(d_phi * wave_length_ant_spacing_ratio / (2*PI));

	target_angle = target_angle * 180 / PI;	// Angle (-90...90°)

	target_angle = target_angle + (double)((int32_t) angle_offset_deg + ANGLE_QUANTIZATION * 0.5);

	delta_angle  = fmodf(target_angle , (double)ANGLE_QUANTIZATION);

	target_angle -= delta_angle;

	temp.d_phi = d_phi;

	temp.target_angle = target_angle;

	return temp;
}

//============================================================================

double get_phase(float real, float imag)
{
	double phi;

	/* Phase angle (0 to 2Pi) */
	if((real > 0) && (imag >= 0))		// 1st quadrant
	{
		phi = atan((double)imag / (double)real);
	}
	else if((real < 0) && (imag >= 0))	// 2nd quadrant
	{
		phi = atan((double)imag / (double)real) + PI;
	}
	else if((real < 0) && (imag <= 0)) 	// 3rd quadrant
	{
		phi = atan((double)imag / (double)real) + PI;
	}
	else if((real > 0) && (imag <= 0)) 	// 4th quadrant
	{
		phi = atan((double)imag / (double)real) + 2*PI;
	}
	else if((real == 0) && (imag > 0))
	{
		phi = PI/2;
	}
	else if((real == 0) && (imag < 0))
	{
		phi = 3*PI/2;
	}
	else
	{
		phi = 0;	// Indeterminate
	}

	return(phi);
}

/* --- End of File -------------------------------------------------------- */
