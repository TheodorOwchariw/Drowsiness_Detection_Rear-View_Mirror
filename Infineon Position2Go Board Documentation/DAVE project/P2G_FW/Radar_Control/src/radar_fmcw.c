/**
    @file: radar_fmcw.c

    @brief: This file contains the driver API for FMCW operation of
            Infineon radar sensors.
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

#include "radar_api.h"
#include "radar_control.h"
#include "ds_radar_interface.h"
#include "position2go_bsp.h"

/*
==============================================================================
   2. EXPORTED FUNCTIONS
==============================================================================
 */

uint16_t radar_set_fmcw_configuration(Radar_Handle_t device, 
                                      const Fmcw_Configuration_t* configuration)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  float tmp_bandwidth_mhz;
  
  /* Configuration change is only allowed while automatic trigger is not working */
  if (radar_get_data_acq_mode() != RADAR_MANUAL_TRIGGER)
  {
    return (RADAR_ERR_BUSY);
  }
  
  /* Check if configuration is valid */
  /* ------------------------------- */
  if ((configuration->lower_frequency_kHz < BSP_MIN_RF_FREQUENCY_KHZ) ||
      (configuration->lower_frequency_kHz > BSP_MAX_RF_FREQUENCY_KHZ) ||
        (configuration->upper_frequency_kHz < BSP_MIN_RF_FREQUENCY_KHZ) ||
          (configuration->upper_frequency_kHz > BSP_MAX_RF_FREQUENCY_KHZ) ||
            (configuration->lower_frequency_kHz > configuration->upper_frequency_kHz))
  {
    return (RADAR_ERR_FREQUENCY_OUT_OF_RANGE);
  }
  
  tmp_bandwidth_mhz = ((float)(configuration->upper_frequency_kHz - configuration->lower_frequency_kHz)) / 1000.0f;
  
  if ((tmp_bandwidth_mhz < BSP_MIN_BANDWIDTH_MHZ) || (tmp_bandwidth_mhz > BSP_MAX_BANDWIDTH_MHZ))
  {
    return (RADAR_ERR_FREQUENCY_OUT_OF_RANGE);
  }
  
  if (configuration->direction != RADAR_DIR_UPCHIRP_ONLY)
  {
    return (RADAR_ERR_UNSUPPORTED_DIRECTION);
  }
  
  if (configuration->tx_power > BSP_MAX_TX_POWER_LEVEL)
  {
    return (RADAR_ERR_POWER_OUT_OF_RANGE);
  }
  
  ds_ep_store_fmcw_config(configuration);
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_fmcw_configuration(Radar_Handle_t device, 
                                      Fmcw_Configuration_t* configuration)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  /* Get the current FMCW configuration from data store */
  ds_ep_fetch_fmcw_config(configuration);
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_bandwidth_per_second(Radar_Handle_t device, 
                                        uint32_t* bandwidth_per_second_MHz_s)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  uint32_t remainder;
  uint32_t divisor;
  Frame_Format_t frame_format;
  
  /* Get the current FMCW configuration from data store */
  ds_ep_fetch_frame_format(&frame_format);
  
  /* Unit is MHz/s. */
  /*
  * Actual calculation "uBandwidthPerSecond = (uUpperFMCWFrequency - uLowerFMCWFrequency) * 1000000 / chirp_duration"
  * is split up to avoid overflow of 32 bit words
  */
  remainder = (uint32_t)((ds_device_get_settings()->pll_upper_frequency_kHz - ds_device_get_settings()->pll_lower_frequency_kHz) *
              (BSP_REFERENCE_OSC_FREQ_HZ / 1000000));
  
  uint32_t sample_rate_divider = BSP_REFERENCE_OSC_FREQ_HZ / ds_device_get_settings()->adc_sampling_freq_Hz;
  
  divisor = frame_format.num_samples_per_chirp * sample_rate_divider;
  
  *bandwidth_per_second_MHz_s = remainder / divisor;
  
  remainder -= *bandwidth_per_second_MHz_s * divisor;
  
  *bandwidth_per_second_MHz_s *= 1000;
  
  *bandwidth_per_second_MHz_s += (remainder * 1000) / divisor;
  
  return (RADAR_ERR_OK);
}

/* --- End of File -------------------------------------------------------- */
