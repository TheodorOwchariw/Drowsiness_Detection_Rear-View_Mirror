/**
    @file: position2go_bgt.c

    @brief: This file implements the function declarations for Infineon's BGT
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

#include "position2go_bgt.h"
#include "bgt24mtr1x.h"
#include "config.h"

/*
==============================================================================
   2. DATA
==============================================================================
 */

static SPI_DIGITAL_IOs_t bgt24mtr1x_spi_pins; /**< SPI GPIOs used to configure BGT24MTR1x device */

/*
==============================================================================
   3. EXPORTED FUNCTIONS
==============================================================================
 */

void bgt_init(uint32_t lna_gain, uint32_t power_level)
{
  /* Set the digital IOs used by SPI interface */
  bgt24mtr1x_spi_pins.spi_device    = BGT_DEVICE;
  bgt24mtr1x_spi_pins.gpio_cs_pin   = (DIGITAL_IO_t*)&DIGITAL_IO_SPI_M_CS_BGT24;
  bgt24mtr1x_spi_pins.gpio_data_pin = (DIGITAL_IO_t*)&DIGITAL_IO_SPI_M_DATA;
  bgt24mtr1x_spi_pins.gpio_clk_pin  = (DIGITAL_IO_t*)&DIGITAL_IO_SPI_M_CLK;
  
  bgt24mtr1x_init((Bgt24mtr1x_LNAgain_t)lna_gain, (Bgt24mtr1x_Power_t)power_level);
  
  bgt24mtr1x_ana_vref_tx();
}

//============================================================================

void bgt_start_tx(void)
{
  bgt24mtr1x_start_tx();
}

//============================================================================

void bgt_stop_tx(void)
{
  bgt24mtr1x_stop_tx();
}

//============================================================================

void bgt_power_up(void)
{
  /* Before turning on BGT, we should keep CE pin high.
   * CE pin is active low, so it should keep high until it is activated.
   */
  DIGITAL_IO_SetOutputHigh(&DIGITAL_IO_SPI_M_CS_BGT24);
  DIGITAL_IO_SetOutputLow(&DIGITAL_IO_BGT_POWER_ENABLE);
}

//============================================================================

void bgt_power_down(void)
{
  /* After turning off BGT, we should keep SPI's signals low
   * to avoid offset voltage at BGT's VCC.
   *
   * If they are above 0.3V which is bias voltage inside BGT, then inside BGT turns on.
   * It makes offset voltage at BGT's VCC.
   */
  DIGITAL_IO_SetOutputHigh(&DIGITAL_IO_BGT_POWER_ENABLE);
  DIGITAL_IO_SetOutputLow(&DIGITAL_IO_SPI_M_CS_BGT24);
  DIGITAL_IO_SetOutputLow(&DIGITAL_IO_SPI_M_DATA);
}

//============================================================================

void bgt_set_tx_power(uint8_t power_level)
{
  bgt24mtr1x_set_tx_power((Bgt24mtr1x_Power_t)power_level);
}

//============================================================================

uint8_t bgt_get_tx_power(void)
{
  return bgt24mtr1x_get_tx_power();
}

//============================================================================

void bgt_lna_gain_enable(void)
{
  bgt24mtr1x_lna_gain_enable();
}

//============================================================================

void bgt_lna_gain_disable(void)
{
  bgt24mtr1x_lna_gain_disable();
}

//============================================================================

uint8_t bgt_lna_gain_is_enable(void)
{
  return bgt24mtr1x_lna_gain_is_enable();
}

//============================================================================

void bgt_set_config(uint16_t config_val)
{
  bgt24mtr1x_set_config(config_val);
}

//============================================================================

uint16_t bgt_get_config(void)
{
  return bgt24mtr1x_get_config();
}

//============================================================================

void bgt_ana_temp(void)
{
  bgt24mtr1x_ana_temp();
}

//============================================================================

void bgt_ana_vout_tx(void)
{
  bgt24mtr1x_ana_vout_tx();
}

//============================================================================

void bgt_ana_vref_tx(void)
{
  bgt24mtr1x_ana_vref_tx();
}

//============================================================================

uint16_t bgt_get_ana_config(void)
{
  return bgt24mtr1x_get_ana_config();
}

//============================================================================

void bgt_lowest_power_with_q2_disable(void)
{
  bgt24mtr1x_set_config((uint16_t)BGT24MTR1X_POWER_CONF);
}

//============================================================================

void bgt_ldo_enable(void)
{
  DIGITAL_IO_SetOutputHigh(&DIGITAL_IO_BGT_LDO_ENABLE);
}

//============================================================================

void bgt_ldo_disable(void)
{
  DIGITAL_IO_SetOutputLow(&DIGITAL_IO_BGT_LDO_ENABLE);
}

/* --- End of File -------------------------------------------------------- */
