/**
    @file: config_error.h

    @brief: Compilation errors check for configuration file values defined by user.
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

#ifndef CONFIG_ERROR_H_
#define CONFIG_ERROR_H_

/* Enable C linkage if header is included in C++ files */
#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*
==============================================================================
   1. INCLUDE FILES
==============================================================================
 */

#include "config.h"
#include "position2go_bsp.h"

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
 */

//===================== General settings error check =======================//

/* Check valid range of Modulation */
#ifdef FW_MODULATION_TYPE
	#if (FW_MODULATION_TYPE < 0) || (FW_MODULATION_TYPE > 1)
		#error "FW_MODULATION_TYPE in config.h file is out of range [0 , 1]"
	#endif
#endif

/* Check the number of chirps in Doppler Modulation */
#if (FW_MODULATION_TYPE != 1) && (NUM_OF_CHIRPS != 1)
	#error "NUM_OF_CHIRPS in config.h must be 1 for Doppler modulation mode"
#endif

/* Check valid range of Number of chirps */
#ifdef NUM_OF_CHIRPS
	#if (NUM_OF_CHIRPS <= 0) || (NUM_OF_CHIRPS > 64)
		#error "NUM_OF_CHIRPS in config.h file is out of range [1 , 64]"
	#endif
#endif

/* Check valid range of samples per chirp */
#ifdef SAMPLES_PER_CHIRP
	#if (SAMPLES_PER_CHIRP < BSP_MIN_NUM_SAMPLES_PER_CHIRP) || (SAMPLES_PER_CHIRP > BSP_MAX_NUM_SAMPLES_PER_CHIRP)
		#error "SAMPLES_PER_CHIRP in config.h file is out of range [32 , 256]"
	#endif
#endif

/* Check valid range of frame period */
#ifdef FRAME_PERIOD_MSEC
	#if (FRAME_PERIOD_MSEC < 1)
		#error "FRAME_PERIOD_MS in config.h file is set below the minimum value of 1 ms."
	#endif
#endif

//==================== BGT / PGA settings error check ======================//

/* Check valid range of duty cycle mode enable/disable */
#ifdef DUTY_CYCLE_ENABLE
	#if (DUTY_CYCLE_ENABLE < 0) || (DUTY_CYCLE_ENABLE > 1)
		#error "DUTY_CYCLE_ENABLE in config.h file is out of range [0 , 1]"
	#endif
#endif

/* Check valid range of BGT TX power */
#ifdef BGT_TX_POWER
	#if (BGT_TX_POWER < 0) || (BGT_TX_POWER > BSP_MAX_TX_POWER_LEVEL)
		#error "BGT_TX_POWER in config.h file is out of range [0 , 7]"
	#endif
#endif

/* Check valid range of LNA gain enable/disable */
#ifdef LNA_GAIN_ENABLE
	#if (LNA_GAIN_ENABLE < 0) || (LNA_GAIN_ENABLE > 1)
		#error "LNA_GAIN_ENABLE in config.h file is out of range [0 , 1]"
	#endif
#endif

/* Check valid range of PGA gain */
#ifdef PGA_GAIN
	#if (PGA_GAIN < 0) || (PGA_GAIN > 7)
		#error "PGA_GAIN in config.h file is out of range [0 , 7]"
	#endif
#endif

/* Check signals parts formats of sampled radar data */
#ifdef RADAR_SIGNAL_PART
	#if defined(DOPPLER_SUPPORTED) || defined(FMCW_SUPPORTED)
		#if (RADAR_SIGNAL_PART != 2)
			#error "RADAR_SIGNAL_PART for I or Q only is not yet supported"
		#endif
	#else
		#if (RADAR_SIGNAL_PART < 0) || (RADAR_SIGNAL_PART > 2)
			#error "RADAR_SIGNAL_PART in config.h file is out of range [0 , 2]"
		#endif
	#endif
#endif

//======================== DSP settings error check ========================//

/* Check valid range of FFT length with zero padding */
#ifdef RANGE_FFT_SIZE
	#if (RANGE_FFT_SIZE < 64) || (RANGE_FFT_SIZE > 512)
		#error "RANGE_FFT_SIZE in config.h file is out of range [64 , 512]"
	#endif
#endif

/* Check valid Doppler of FFT length with zero padding */
#ifdef DOPPLER_FFT_SIZE
	#if (DOPPLER_FFT_SIZE < 32) || (DOPPLER_FFT_SIZE > 64)
		#error "DOPPLER_FFT_SIZE in config.h file is out of range [32 , 64]"
	#endif
#endif

/* Check valid range of FFT input type */
#ifdef FFT_INPUT_TYPE
	#if (FFT_INPUT_TYPE < 0) || (FFT_INPUT_TYPE > 2)
		#error "FFT_INPUT_TYPE in config.h file is out of range [0 , 2]"
	#endif
#endif

/* Check valid range of maximum number of targets to be detected */
#ifdef MAX_NUM_OF_TARGETS
	#if (MAX_NUM_OF_TARGETS < 1) || (MAX_NUM_OF_TARGETS > 5)
		#error "NUM_OF_TARGETS in config.h file is out of range [1 , 5]"
	#endif
#endif

/* Check valid range of RX antenna selection */
#ifdef RX_ANTENNA_SELECTION
	#if (RX_ANTENNA_SELECTION < 1) || (RX_ANTENNA_SELECTION > 2)
		#error "RX_ANTENNA_SELECTION in config.h file is out of range [1 , 2]"
	#endif
#endif

/* Check valid range of Algorithm processing time */
#ifdef ALGO_PROCESS_TIME_USEC
	#if (ALGO_PROCESS_TIME_USEC < 0)
		#error "ALGO_PROCESS_TIME_USEC in config.h file is set below the minimum value of 0 ms."
	#endif
#endif

//=============== Tracking and Filter settings error check =================//

/* Check valid range of tracking enable/disable */
#ifdef TRACKING_ENABLE
	#if (TRACKING_ENABLE < 0) || (TRACKING_ENABLE > 1)
		#error "TRACKING_ENABLE in config.h file is out of range [0 , 1]"
	#endif
#endif

/* Check valid range of number of detected targets to be tracked */
#ifdef CURRENT_NUM_OF_TRACKS
	#if (CURRENT_NUM_OF_TRACKS < 1) || (CURRENT_NUM_OF_TRACKS > 5)
		#error "CURRENT_NUM_OF_TRACKS in config.h file is out of range [1 , 5]"
	#endif
#endif

/* Check valid range of moving average window length */
#ifdef MVG_AVG_LEN
	#if (MVG_AVG_LEN < 1) || (MVG_AVG_LEN > 10)
		#error "MVG_AVG_LEN in config.h file is out of range [1 , 10]"
	#endif
#endif

/* Check valid range of length of the median filter for angle */
#ifdef MAX_MEDIAN_FILTER_LEN
	#if (MAX_MEDIAN_FILTER_LEN < 1) || ((MAX_MEDIAN_FILTER_LEN % 2) == 0)
		#error "MAX_MEDIAN_FILTER_LEN in config.h file is out of range [1 , ..] or is not odd number"
	#endif
#endif

/* Check valid range of length of the median filter for angle */
#ifdef MEDIAN_FILTER_LEN
	#if (MEDIAN_FILTER_LEN < 1) || ((MEDIAN_FILTER_LEN % 2) == 0)
		#error "MEDIAN_FILTER_LEN in config.h file is out of range [1 , ..] or is not odd number"
	#endif

	#if !defined(MAX_MEDIAN_FILTER_LEN)
		#error "MAX_MEDIAN_FILTER_LEN not defined"
	#endif

	#if (MEDIAN_FILTER_LEN > MAX_MEDIAN_FILTER_LEN)
		#error "MEDIAN_FILTER_LEN in config.h file is out of range [1 , ..]"
	#endif
#endif

/* Check valid range of number of degrees to be quantized */
#ifdef ANGLE_QUANTIZATION
	#if (ANGLE_QUANTIZATION < 1)
		#error "ANGLE_QUANTIZATION in config.h file is out of range [1 , ..]"
	#endif
#endif

/* Check valid range of MTI filter enable/disable */
#ifdef MTI_FILTER_ENABLE
	#if (MTI_FILTER_ENABLE < 0) || (MTI_FILTER_ENABLE > 1)
		#error "MTI_FILTER_ENABLE in config.h file is out of range [0 , 1]"
	#endif
#endif

/* Check valid range of MTI filter weight */
#ifdef MTI_FILTER_LEN
	#if (MTI_FILTER_LEN < 1) || (MTI_FILTER_LEN > BSP_MAX_MTI_FILTER_LEN)
		#error "MTI_FILTER_LEN in config.h file is out of range [1 , 1000]"
	#endif
#endif

//======================= FMCW settings error check ========================//

/* Check enabled supported algorithms */
#if ((defined(FMCW_SUPPORTED) && (FW_MODULATION_TYPE == 0)) || \
    (defined(DOPPLER_SUPPORTED) && (FW_MODULATION_TYPE == 1)))
	#error "Please ensure to enable the appropriate algorithm define under your project config.h file"
#endif

/* Check valid range of Bandwidth */
#ifdef BANDWIDTH_MHZ
	#if (BANDWIDTH_MHZ < 1) || (BANDWIDTH_MHZ > 200)
		#error "BANDWIDTH_MHZ in config.h file is out of range [1 , 200]"
	#endif
#endif

/* Check valid range of Chirp-time */
#ifdef CHIRP_TIME_USEC
	#if (CHIRP_TIME_USEC < BSP_MIN_CHIRP_TIME_USEC) || (CHIRP_TIME_USEC > BSP_MAX_CHIRP_TIME_USEC)
		#error "CHIRP_TIME_USEC in config.h file is out of range [50 , 3000]"
	#endif
	#if (CHIRP_TIME_USEC < 300) && (SAMPLES_PER_CHIRP > 256)
		#error "Sampling rate not supported with current chirp time less than 300us and SAMPLES_PER_CHIRP greater than 256"
	#endif
	#if (CHIRP_TIME_USEC < 150) && (SAMPLES_PER_CHIRP > 128)
		#error "Sampling rate not supported with current chirp time less than 200us and SAMPLES_PER_CHIRP greater than 128"
	#endif
	#if (CHIRP_TIME_USEC < 75) && (SAMPLES_PER_CHIRP > 64)
		#error "Sampling rate not supported with current chirp time less than 100us and SAMPLES_PER_CHIRP greater than 64"
	#endif
#endif

/* Check valid range of Minimum range to be filtered out by FMCW algorithm */
#if defined (MINIMUM_RANGE_CM) && defined (MAXIMUM_RANGE_CM)
	#if (MINIMUM_RANGE_CM < 0)
		#error "MINIMUM_RANGE_CM in config.h file should be non-negative value"
	#endif
	#if (MINIMUM_RANGE_CM > MAXIMUM_RANGE_CM)
		#error "MINIMUM_RANGE_CM in config.h file should not be greater than MAXIMUM_RANGE_CM"
	#endif
#endif

/* Check valid range of Maximum range to be filtered out by FMCW algorithm */
#if defined (MAXIMUM_RANGE_CM) && defined (MINIMUM_RANGE_CM)
	#if (MAXIMUM_RANGE_CM < 0)
		#error "MAXIMUM_RANGE_CM in config.h file should be non-negative value"
	#endif
	#if (MAXIMUM_RANGE_CM < MINIMUM_RANGE_CM)
		#error "MAXIMUM_RANGE_CM in config.h file should not be less than MINIMUM_RANGE_CM"
	#endif
#endif

/* Check valid range of FMCW FFT detection threshold */
#ifdef RANGE_DETECTION_THRESHOLD
	#if (RANGE_DETECTION_THRESHOLD < 0)
		#error "RANGE_DETECTION_THRESHOLD in config.h file should be non-negative value"
	#endif
#endif

//===================== Doppler settings error check =======================//

/* Check valid range of Doppler Sampling frequency in Hz */
#ifdef DOPPLER_SAMPLING_FREQ_HZ
	#if (DOPPLER_SAMPLING_FREQ_HZ < 0)
		#error "DOPPLER_SAMPLING_FREQ_HZ in config.h file should be non-negative value"
	#endif
	#if (DOPPLER_SAMPLING_FREQ_HZ > 100000)
		#error "DOPPLER_SAMPLING_FREQ_HZ in config.h file is set too high, it should be less than 100kHz"
	#endif
#endif

/* Check valid range of Minimum speed to be filtered out by Doppler algorithm */
#if defined (MINIMUM_SPEED_KMH) && defined (MAXIMUM_SPEED_KMH)
	#if (MINIMUM_SPEED_KMH < 0)
		#error "MINIMUM_SPEED_KMH in config.h file should be non-negative value"
	#endif
	#if (MINIMUM_SPEED_KMH > MAXIMUM_SPEED_KMH)
		#error "MINIMUM_SPEED_KMH in config.h file should not be greater than MAXIMUM_SPEED_KMH"
	#endif
#endif

/* Check valid range of Maximum speed to be filtered out by Doppler algorithm */
#if defined (MAXIMUM_SPEED_KMH) && defined (MINIMUM_SPEED_KMH)
	#if (MAXIMUM_SPEED_KMH < 0)
		#error "MAXIMUM_SPEED_KMH in config.h file should be non-negative value"
	#endif
	#if (MAXIMUM_SPEED_KMH < MINIMUM_SPEED_KMH)
		#error "MAXIMUM_SPEED_KMH in config.h file should not be less than MINIMUM_SPEED_KMH"
	#endif
#endif

/* Check valid range of FFT spectrum threshold to detect a target in Doppler */
#ifdef SPEED_DETECTION_THRESHOLD
	#if (SPEED_DETECTION_THRESHOLD < 0)
		#error "SPEED_DETECTION_THRESHOLD in config.h file should be non-negative value"
	#endif
#endif

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* CONFIG_ERROR_H_ */

/* --- End of File -------------------------------------------------------- */
