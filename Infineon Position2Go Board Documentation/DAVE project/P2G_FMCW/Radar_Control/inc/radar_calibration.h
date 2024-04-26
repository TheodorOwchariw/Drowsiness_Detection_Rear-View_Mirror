/**
  @file radar_calibration.h
 
  @brief This file defines the driver API for Calibration features of
         Infineon radar sensors.
 
  This file defines types and functions that are used to set/get the
  calibration data in EEPROM Flash memory to compensate TX to RX leakage
  and some offset compensations.
 
  @note This file is not meant for being included directly. It will be
        included by the main driver API header 'ifxRadar.h' from the device
        specific driver folder.
*/
 
/* ===========================================================================
** Copyright (C) 2018-2019 Infineon Technologies AG
** All rights reserved.
** ===========================================================================
**
** ===========================================================================
** This document contains proprietary information of Infineon Technologies AG.
** Passing on and copying of this document, and communication of its contents
** is not permitted without Infineon's prior written authorisation.
** ===========================================================================
*/

#ifndef RADAR_CALIBRATION_H_
#define RADAR_CALIBRATION_H_

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
#include "position2go_bsp.h"

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
 */

#define 	DEFAULT_CALIBRATION_VALUE		(0U)		/**<  Default calibration value to be written in to the Flash */

#define		MAX_CALIB_SAMPLES_PER_CHANNEL	(256U)		/**< Number of data samples stored in the calibration EEPROM per ADC channel */

#define		ADC_CALIB_LEN_BYTES				(BSP_NUM_OF_ADC_CHANNELS * MAX_CALIB_SAMPLES_PER_CHANNEL * 2U)

/*
==============================================================================
   3. TYPES
==============================================================================
 */

/**
 * \brief Data structure for Algorithm related Calibrations values.
 * @{
 */
typedef struct
{
	uint16_t 	distance_offset_cm;		/**< Distance offset in units of cm */
	int16_t		angle_offset_deg;		/**< Angle offset in units of degrees */
} 	Algo_Calibrations_t;
/** @} */

/*
==============================================================================
   4. DATA
==============================================================================
 */

/* Data memory section in the SRAM containing the current calibration contents of EEPROM */
extern uint16_t    adc_i1_calibration[MAX_CALIB_SAMPLES_PER_CHANNEL * BSP_NUM_OF_ADC_CHANNELS];
extern uint16_t*   adc_q1_calibration;
extern uint16_t*   adc_i2_calibration;
extern uint16_t*   adc_q2_calibration;

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* RADAR_CALIBRATION_H_ */

/* --- End of File -------------------------------------------------------- */
