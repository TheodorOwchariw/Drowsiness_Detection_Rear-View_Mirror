/**
    @file: application.h

    @brief: This file includes the function declarations application layer.
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

#ifndef APPLICATION_H_
#define APPLICATION_H_

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

#include "DAVE.h"
#include "version.h"
#include "config.h"
#include "config_error.h"
#include "radar_control.h"

/* Host Communication Library includes */
#include "EndpointCalibration.h"
#include "EndpointTargetDetection.h"
#include "EndpointRadarADCXMC.h"
#include "EndpointRadarBase.h"
#include "EndpointRadarDoppler.h"
#include "EndpointRadarFMCW.h"
#include "EndpointRadarP2G.h"
#include "EndpointRadarIndustrial.h"
#include "Protocol.h"

/* Algorithms' includes */
#ifdef FMCW_SUPPORTED
#include "range_doppler.h"
#endif /* FMCW_SUPPORTED */

#ifdef DOPPLER_SUPPORTED
#include "doppler.h"
#endif /* DOPPLER_SUPPORTED */
  
/*
==============================================================================
   2. TYPES
==============================================================================
 */

/**
 * \brief Defines the main radar application states.
 * @{
 */
typedef enum
{
   APP_IDLE                        = 1U, /* The default application state */
   APP_APPLY_DEVICE_SETTINGS,
   APP_RADAR_DATA_ACQUISITION_COMPLETED, /* Data acquisition completed and ADC raw data is available in the buffer */
   APP_CHECK_FOR_SETTINGS_UPDATE,        /* Check if a new hardware settings are required before to start acquisition for the next frame*/
} Radar_App_State_t;
/** @} */

/*
==============================================================================
   3. FUNCTION PROTOTYPES
==============================================================================
*/

/**
 * \brief  This function initializes all firmware and algorithm parameters, starts the
 *         radar data acquisition and initiates the communication with Radar GUI.
 *
 * \param[in]  None
 *
 * \return None
 */
void app_init(void);

/**
 * \brief  This function registers algorithm processing function in application layer
 *
 * \param[in] algo_processor  Processing algorithm need to be assigned to this
 *                            function pointer of signature "void *funcptr(void)".
 *
 * \return None
 */
void app_register_algo_process(algorithm algo_processor);

/**
 * \brief  This function is the main application process, called from the main program
 *         loop and it manages the communication with the Radar GUI.
 *
 * \param[in]  None
 *
 * \return None
 */
void app_process(void);

/**
 * \brief  This function is used to de-initialize and free all used
 *         resources by the application.
 *
 * \param[in]  None
 *
 * \return None
 */
void app_destructor(void);

/**
 * \brief  This function is used for software reset.
 *
 * \param[in]  None
 *
 * \return None
 */
void system_reset(void);

/**
 * \brief  The main radar application process.
 *
 * \param[in]  None
 *
 * \return None
 */
void radar_app_process(void);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* APPLICATION_H_ */

/* --- End of File -------------------------------------------------------- */
