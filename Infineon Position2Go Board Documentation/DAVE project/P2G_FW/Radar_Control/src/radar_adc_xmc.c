/**
    @file: radar_adc_xmc.c

    @brief: This file contains the basic driver API to configure the ADC that's
            built in XMC microcontrollers used in Infineon radar sensors.
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

#include "ds_radar_interface.h"
#include "radar_api.h"
#include "radar_control.h"
#include "position2go_bsp.h"

/*
==============================================================================
   2. EXPORTED FUNCTIONS
==============================================================================
 */

uint16_t radar_set_adc_configuration(Radar_Handle_t device, const Adc_Configuration_t* configuration)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  uint32_t chirp_time_usec;
  
  /* Configuration change is only allowed while automatic trigger is not working */
  if (radar_get_data_acq_mode() != RADAR_MANUAL_TRIGGER)
  {
    return (RADAR_ERR_BUSY);
  }
  
  /* Check if configuration is valid */
  /* ------------------------------- */
  if ((configuration->samplerate_Hz == 0) || (configuration->samplerate_Hz > BSP_MAX_ADC_FREQ_HZ))
  {
    return (RADAR_ERR_SAMPLERATE_OUT_OF_RANGE);
  }
  
  if (configuration->resolution != BSP_ADC_RESOLUTION)
  {
    return (RADAR_ERR_UNSUPPORTED_RESOLUTION);
  }
  
  /* Apply ADC configuration */
  /* ----------------------- */
  //sample_rate_divider = ((REFERENCE_OSC_FREQ_HZ << 1) / configuration->samplerate_Hz + 1) >> 1;
  
  if (ds_device_get_settings()->pll_modulation_mode == MODULATION_FMCW)
  {
    chirp_time_usec = (uint32_t) (((float)ds_device_get_settings()->num_samples_per_chirp * 1000.0f ) / configuration->samplerate_Hz * 1000.0f);
    
    if (radar_set_chirp_time(device, chirp_time_usec) != RADAR_ERR_OK)
    {
      return (RADAR_ERR_SAMPLERATE_OUT_OF_RANGE);
    }
  }
  else
  {
    /* Apply ADC configuration */
    /* ----------------------- */
    if (radar_set_sampling_freq(device, MODULATION_DOPPLER, ds_device_get_settings()->adc_sampling_freq_Hz) != RADAR_ERR_OK)
    {
      return (RADAR_ERR_SAMPLERATE_OUT_OF_RANGE);
    }
  }
  
  /* Remember configuration (adjust frequency to rounded value) */
  /* ---------------------------------------------------------- */
  /* Copy the new adc configuration to the data store*/
  ds_ep_store_adc_config(configuration);
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_adc_configuration(Radar_Handle_t device, Adc_Configuration_t* configuration)
{
  XMC_UNUSED_ARG(device);  /* Suppress compiler warnings of unused variable */
  
  if(configuration != NULL)
  {
    ds_ep_fetch_adc_config(configuration);
    return (RADAR_ERR_OK);
  }
  return (RADAR_ERR_FEATURE_UNAVAILABLE);
}

/* --- End of File -------------------------------------------------------- */
