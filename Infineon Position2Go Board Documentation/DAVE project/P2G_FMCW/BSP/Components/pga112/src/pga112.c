/**
    @file: pga112.c

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

#include "pga112.h"

/*
==============================================================================
   2. DATA
==============================================================================
 */

static volatile  uint16_t pga112_global_gain = 0; /**< Current PGA112 gain */

static volatile  uint16_t pga112_global_config = PGA112_BASE_CONF; /**< Current PGA112 configuration */

static SPI_DIGITAL_IOs_t  pga112_spi_pins; /**< SPI GPIOs used to configure PGA112 device */

/*
==============================================================================
   3. LOCAL FUNCTION PROTOTYPES
==============================================================================
 */

/**
 * \brief This function returns the PGA112 gain value from the selected index.
 *
 * \param[in]  gain_idx  index of the PGA112 gain value (from 0 to 7)
 *
 * \return PGA112 gain value
 */
static Pga112_Binary_Gain_t pga112_get_binary_gain(uint16_t gain_idx);

/*
==============================================================================
   4. EXPORTED FUNCTIONS
==============================================================================
 */

Pga_Status_t pga112_init(uint16_t gain_level)
{
  /* Set the digital IOs used by SPI interface */
  pga112_spi_pins.spi_device = PGA_DEVICE;
  pga112_spi_pins.gpio_cs_pin = (DIGITAL_IO_t*)&DIGITAL_IO_SPI_M_CS_PGA;
  pga112_spi_pins.gpio_data_pin = (DIGITAL_IO_t*)&DIGITAL_IO_SPI_DATA_PGA;
  pga112_spi_pins.gpio_clk_pin = (DIGITAL_IO_t*)&DIGITAL_IO_SPI_M_CLK;

  return (pga112_set_gain(gain_level));
}

//============================================================================

Pga_Status_t pga112_get_gain(uint16_t* gain_level)
{
  *gain_level = pga112_global_gain;

  return (PGA_STATUS_SUCCESS);
}

//============================================================================

Pga_Status_t pga112_set_gain(uint16_t gain_level)
{
  Pga112_Binary_Gain_t  gain_val;
  uint16_t shuffled_data = 0;

  if(gain_level > PGA112_MAX_NUMBER_SUPPORTED_GAIN)
  {
    return (PGA_STATUS_FAIL);
  }

  /* Get PGA112 gain value from the selected index */
  pga112_global_gain = gain_level;
  gain_val = pga112_get_binary_gain(gain_level);

  pga112_global_config = (PGA112_CMD_WRITE | gain_val | PGA112_CH1);

  /* Lower byte goes to upper byte position */
  shuffled_data  = ((pga112_global_config & 0x00FF) << 8);

  /* Upper byte goes to lower byte position */
  shuffled_data |= ((pga112_global_config & 0xFF00) >> 8);

  /* Send command to PGA112 */
  spi_transmit_data((uint8_t *)&shuffled_data, 2, (void*) &pga112_spi_pins);

  return (PGA_STATUS_SUCCESS);
}

/*
==============================================================================
   5. LOCAL FUNCTIONS
==============================================================================
 */

static Pga112_Binary_Gain_t pga112_get_binary_gain(uint16_t gain_idx)
{
	Pga112_Binary_Gain_t gain_val;

  switch (gain_idx)
  {
  case 0U:
    gain_val = PGA112_BINARY_GAIN_1;
    break;

  case 1U:
    gain_val = PGA112_BINARY_GAIN_2;
    break;

  case 2U:
    gain_val = PGA112_BINARY_GAIN_4;
    break;

  case 3U:
    gain_val = PGA112_BINARY_GAIN_8;
    break;

  case 5U:
    gain_val = PGA112_BINARY_GAIN_32;
    break;

  case 6U:
    gain_val = PGA112_BINARY_GAIN_64;
    break;

  case 7U:
    gain_val = PGA112_BINARY_GAIN_128;
    break;

  default:
  case 4U:
    gain_val = PGA112_BINARY_GAIN_16;
    break;
  }

  return (gain_val);
}

/* --- End of File -------------------------------------------------------- */
