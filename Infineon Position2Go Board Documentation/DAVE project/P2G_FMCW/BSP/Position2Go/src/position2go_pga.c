/**
    @file: position2go_pga.c

    @brief: This file provides driver for PGA112xxx Zero-Drift, Programmable Gain
            Amplifier. It contains a set of functions to configure and use PGA112 device.
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

#include "bsp_common.h"
#include "position2go_pga.h"
#include "pga112.h"

/*
==============================================================================
   2. EXPORTED FUNCTIONS
==============================================================================
 */

Pga_Status_t pga_init(uint16_t gain_level)
{
  /* Initialize the PGA112 device */
  return pga112_init(gain_level);
}

//============================================================================

void pga_get_gain(uint16_t* gain_level)
{
  pga112_get_gain(gain_level);
}

//============================================================================

Pga_Status_t pga_set_gain(uint16_t gain_level)
{
  /* Updates the PGA112 gain */
  return pga112_set_gain(gain_level);
}

//============================================================================

void pga_ldo_enable(void)
{
  DIGITAL_IO_SetOutputHigh(&DIGITAL_IO_SPI_M_CS_PGA);
  DIGITAL_IO_SetOutputHigh(&DIGITAL_IO_PGA_LDO_ENABLE);
}

//============================================================================

void pga_ldo_disable(void)
{
  DIGITAL_IO_SetOutputLow(&DIGITAL_IO_SPI_DATA_PGA);
  DIGITAL_IO_SetOutputLow(&DIGITAL_IO_SPI_M_CS_PGA);
  DIGITAL_IO_SetOutputLow(&DIGITAL_IO_PGA_LDO_ENABLE);
}

/* --- End of File -------------------------------------------------------- */
