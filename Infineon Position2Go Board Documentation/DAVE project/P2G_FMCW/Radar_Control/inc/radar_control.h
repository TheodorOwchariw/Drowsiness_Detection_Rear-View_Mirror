/**
    @file: radar_control.h

  	@brief: Common API definition valid for all radar operations related to MCU, BGT, PLL
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

#ifndef RADAR_CONTROL_H_
#define RADAR_CONTROL_H_

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

#include "acq_buffer.h"
#include "datastore.h"
  
/*
==============================================================================
   2. DEFINITIONS
==============================================================================
*/

/**
 * \brief Timing constants, and Key parameters of the radar can be restricted to these values.
 * @{
 */

#define 	PLL_LOCK_TIME_USEC		(50U)		/**< This is the minimum time in microseconds (usec) to program the PLL */

#define 	BGT_VCO_LOCK_TIME_USEC	(1000U) 	/**< This is the minimum time in microseconds (usec) to program after BGT and PLL
													 power-up in duty cycle mode. After 1500us BGT Q2 divider settled to zero */

/** @} */

/*
==============================================================================
   3. TYPES
==============================================================================
*/

/**
 * \brief Defines radar function callback ID.
 * @{
 */
typedef enum
{
  RADAR_ACQUISITION_STARTED_CB_ID   = 1U,   /**< Start of data acquisition callback ID */
  RADAR_ACQUISITION_DONE_CB_ID      = 2U    /**< End of data acquisition callback ID   */
} Radar_Callback_ID_t;
/** @} */

/**
 * \brief Defines radar data acquisition modes: Manual, Automatic or Test mode.
 * @{
 */
typedef enum
{
  RADAR_MANUAL_TRIGGER    = 0U,   /**< Manual triggering data acquisition */
  RADAR_AUTO_TRIGGER      = 1U,   /**< Automatic triggering data acquisition */
  RADAR_TEST_MODE         = 2U    /**< Test mode, not used yet */
} Radar_Data_Acq_Mode_t;
/** @} */

typedef void(*algorithm)(acq_buf_obj *p_acq_buf, const algo_settings_t *cp_algo_settings, const device_settings_t *cp_dev_settings, algo_result_t *p_algo_result);

/**
 * \brief Defines radar hardware state settings
 * @{
 */
typedef enum
{
 	  INIT,
	  ENTER_POWERSAVING,
	  LEAVE_POWERSAVING,
	  POWER_DOWN,
	  UPDATE,
	  START_MANUAL_ACQ,
	  START_AUTO_ACQ,
	  STOP_ACQ,
} hw_state_setting_t;
/** @} */

typedef void* Radar_Handle_t;

/** @} */

/*
==============================================================================
   4. DATA
==============================================================================
 */

extern volatile uint32_t   g_adc_sampling_completed;   /**< Global variable to inform about the sampling done status, used by algorithm to start processing */

extern volatile uint32_t   g_do_processing;            /**< Global variable to enable/disable post processing after data acquisition, valid range [0 - 1] */

/*
==============================================================================
   5. FUNCTION PROTOTYPES
==============================================================================
 */

/**
 * \brief	This function registers user callback functions in BSP layer.
 *          It delegates user callback to the BSP low layer driver.
 *
 * \param[in] 	cb_id		 Radar callback ID \ref Radar_Callback_ID_t
 * \param[in] 	*pcallback   A pointer to user callback function
 *
 * \return One of the following status can be return:
 *   	   - \ref RADAR_STATUS_SUCCESS    Success
 *   	   - \ref RADAR_STATUS_FAIL       Fail
 *   	   - \ref RADAR_STATUS_NA         Feature not available
 */
uint16_t radar_register_callback(Radar_Callback_ID_t cb_id, void (*pcallback)(void *parameter));

/**
 * \brief	This function unregisters user callback functions from BSP layer.
 *
 * \param[in] 	cb_id		Radar callback ID \ref Radar_Callback_ID_t
 *
 * \return One of the following status can be return:
 *   	   - \ref RADAR_STATUS_SUCCESS    Success
 *   	   - \ref RADAR_STATUS_FAIL       Fail
 *   	   - \ref RADAR_STATUS_NA         Feature not available
 */
uint16_t radar_unregister_callback(Radar_Callback_ID_t cb_id);

/**
 * \brief	This function applies radar hardware settings.
 *
 * \param[in] 	device           A handle to the radar sensor device.
 * \param[in] 	state		     Radar hardware state \ref hw_state_setting_t
 * \param[in] 	*p_hw_settings   A pointer to hardware configuration structure \ref device_settings_t
 *
 * \return One of the following status can be return:
 *   	   - \ref RADAR_STATUS_SUCCESS    Success
 *   	   - \ref RADAR_STATUS_FAIL       Fail
 *   	   - \ref RADAR_STATUS_NA         Feature not available
 */
uint32_t radar_apply_hw_settings(Radar_Handle_t device, hw_state_setting_t state, device_settings_t *p_hw_settings);

/**
 * \brief	This function returns the radar data acquisition mode.
 *
 * \return One of the following modes can be return \ref Radar_Data_Acq_Mode_t:
 *   	   - \ref RADAR_MANUAL_TRIGGER    Manual trigger
 *   	   - \ref RADAR_AUTO_TRIGGER      Automatic trigger
 *   	   - \ref RADAR_TEST_MODE         Test mode
 */
Radar_Data_Acq_Mode_t radar_get_data_acq_mode(void);

/**
 * \brief	This function sets the radar data acquisition mode: Manual, Automatic or Test mode.
 *
 * \param[in] 	new_mode		Radar data acquisition mode \ref Radar_Data_Acq_Mode_t
 */
void radar_set_data_acq_mode(Radar_Data_Acq_Mode_t new_mode);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RADAR_CONTROL_H_ */

/* --- End of File -------------------------------------------------------- */
