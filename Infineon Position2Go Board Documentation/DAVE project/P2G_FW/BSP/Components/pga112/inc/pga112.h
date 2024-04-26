/**
    @file: pga112.h

  	@brief: This file includes the constants declarations and functions prototypes
    for PGA112xxx amplifier driver.
    This file needs to be included to make the full driver APIs available.
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

#ifndef PGA112_H_
#define PGA112_H_

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
#include "bsp_common.h"

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
 */

/**
 * @brief The maximum number of supported gain by PGA112 device.
 */
#define PGA112_MAX_NUMBER_SUPPORTED_GAIN       	(7U)        /**< binary gain 128 */

/**
 * @brief PGA112 default configuration value.
 */
#define PGA112_BASE_CONF			(PGA112_CMD_WRITE | PGA112_BINARY_GAIN_32 | PGA112_CH1)

/*
==============================================================================
   3. TYPES
==============================================================================
 */

/**
 * \brief SPI commands, used to control PGA112 device.
 * @{
 */
typedef enum
{
  PGA112_CMD_READ           = 0x6A00U,     /**< Read command */
  PGA112_CMD_WRITE          = 0x2A00U,     /**< Write command */
  PGA112_CMD_NOP_WRITE      = 0x0000U,     /**< No OPeration command */
  PGA112_CMD_SDN_DIS_WRITE  = 0xE100U,     /**< Exit Shutdown mode */
  PGA112_CMD_SDN_EN_WRITE   = 0xE1F1U      /**< Enter Shutdown mode */
} Pga112_Command_t;
/** @} */

/**
 * @brief PGA112 binary gain levels.
 */
typedef enum
{
  PGA112_BINARY_GAIN_1     = 0x0000U,     /**< binary gain 1 */
  PGA112_BINARY_GAIN_2     = 0x0010U,     /**< binary gain 2 */
  PGA112_BINARY_GAIN_4     = 0x0020U,     /**< binary gain 4 */
  PGA112_BINARY_GAIN_8     = 0x0030U,     /**< binary gain 8 */
  PGA112_BINARY_GAIN_16    = 0x0040U,     /**< binary gain 16 */
  PGA112_BINARY_GAIN_32    = 0x0050U,     /**< binary gain 32 */
  PGA112_BINARY_GAIN_64    = 0x0060U,     /**< binary gain 64 */
  PGA112_BINARY_GAIN_128   = 0x0070U      /**< binary gain 128 */
} Pga112_Binary_Gain_t;
/** @} */

/**
 * @brief PGA112 used channels for RX gain.
 */
typedef enum
{
  PGA112_CH0_VCAL          = 0x0000U,     /**< Input MUX channel 0 and VCAL input */
  PGA112_CH1               = 0x0001U      /**< Input MUX channel 1 */
} Pga112_Channel_t;
/** @} */

/*
==============================================================================
   4. FUNCTION PROTOTYPES
==============================================================================
 */

/**
 * \brief This function configures and initializes PGA112 device.
 *
 * \param[in]  gain_level	Unsigned 8-bit integer, whose values are matched to the Power levels defined by enum \ref PGA_Binary_Gain_t.
 *				There are eight binary gain levels [0 - 7] that can be set with PGA112.
 *
 * \return One of the following status can be return:
 *   	   - \ref PGA_STATUS_SUCCESS	Success
 *   	   - \ref PGA_STATUS_FAIL       Fail
 *   	   - \ref PGA_STATUS_NA         Feature not available
 */
Pga_Status_t pga112_init(uint16_t gain_level);

/**
 * \brief This function sets and updates the power level of PGA112 at the Rx side.
 *
 * \param[in]  gain_level	Unsigned 16-bit integer, whose values are matched to the Power levels defined by enum \ref PGA_Binary_Gain_t.
 *				There are eight binary gain levels [0 - 7] that can be set with PGA112.
 *
 * \return One of the following status can be return:
 *   	   - \ref PGA_STATUS_SUCCESS    Success
 *   	   - \ref PGA_STATUS_FAIL       Fail
 *   	   - \ref PGA_STATUS_NA         Feature not available
 */
Pga_Status_t pga112_set_gain(uint16_t gain_level);

/**
 * \brief This function returns the current power level of PGA112 at the Rx side.
 *
 * \param[out]  *gain_level	Pointer to Unsigned 16-bit integer, whose values are within the range of [1 - 8].
 *
 * \return One of the following status can be return:
 *   	   - \ref PGA_STATUS_SUCCESS    Success
 *   	   - \ref PGA_STATUS_FAIL       Fail
 *   	   - \ref PGA_STATUS_NA         Feature not available
 */
Pga_Status_t pga112_get_gain(uint16_t* gain_level);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* PGA112_H_ */

/* --- End of File -------------------------------------------------------- */
