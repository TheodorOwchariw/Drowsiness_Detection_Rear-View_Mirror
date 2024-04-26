/**
    @file: position2go_pga.h

    @brief: This file includes the constants declarations and functions prototypes
            for PGA112xxx amplifier driver. This file needs to be included to 
            make the full driver APIs available.
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

#ifndef POSITION2GO_PGA_H_
#define POSITION2GO_PGA_H_

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

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
 */

/**
 * @brief Define the PGA maximum gain level.
 */
#define MAX_PGA_GAIN_LEVEL     PGA112_MAX_NUMBER_SUPPORTED_GAIN

/*
==============================================================================
   3. FUNCTION PROTOTYPES
==============================================================================
 */

/**
 * \brief This function is used to configure and initial  PGA112 device.
 *
 * \param[in]  gain_level	Unsigned 8-bit integer, whose values are matched to the Power levels defined by enum \ref PGA_Binary_Gain_t.
 *							There are eight binary gain levels [0 - 7] that can be set with PGA112.
 */
Pga_Status_t pga_init(uint16_t gain_level);

/**
 * \brief This function is used to set the power level of PGA112 at the Rx side.
 *
 * \param[in]  gain_level	Unsigned 16-bit integer, whose values are matched to the Power levels defined by enum \ref PGA_Binary_Gain_t.
 *							There are eight binary gain levels [0 - 7] that can be set with PGA112.
 */
Pga_Status_t pga_set_gain(uint16_t gain_level);

/**
 * \brief This function is used to get the current power level of PGA112 at the Rx side.
 *
 * \param[out]  *gain_level	Pointer to Unsigned 16-bit integer, whose values are within the range of [1 - 8].
 */
void pga_get_gain(uint16_t* gain_level);

/**
 * \brief This function enables the LDO for PGA.
 */
void pga_ldo_enable(void);

/**
 * \brief This function disables the LDO for PGA.
 */
void pga_ldo_disable(void);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* POSITION2GO_PGA_H_ */

/* --- End of File -------------------------------------------------------- */
