/**
    @file: range_doppler.h

    @brief: This file contains RangeDoppler algorithm with Monopulse angle calculation.
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

#ifndef RANGE_DOPPLER_H_
#define RANGE_DOPPLER_H_

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

#include <stdint.h>
#include <stdbool.h>
#include "config.h"
#include "dsp_lib.h"
#include "tracking.h"

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
*/

#define  MAX_MVG_AVG_LEN		(10U)		/**<  Maximum length of moving average filter */

/*
==============================================================================
   3. FUNCTION PROTOTYPES
==============================================================================
*/

/**
 * \brief  Top-level wrapper function for FMCW algorithm.
 *
 * Internally, this function do the following sub tasks;
 * 	- Subtract calibration data from raw ADC IQ
 * 	- Calls standard FMCW core functions to compute targets based on range
 */
void range_doppler_do(acq_buf_obj *p_acq_buf, const algo_settings_t *cp_algo_settings, const device_settings_t *cp_dev_settings, algo_result_t *p_algo_result);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RANGE_DOPPLER_H_ */

/* --- End of File -------------------------------------------------------- */
