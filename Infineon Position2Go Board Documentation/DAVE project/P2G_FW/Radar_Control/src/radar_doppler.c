/**
    @file: radar_doppler.c

  	@brief: This file contains the driver API for Doppler operation of
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

uint16_t radar_set_doppler_configuration(Radar_Handle_t device, const Doppler_Configuration_t* const configuration)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  /* Configuration change is only allowed while automatic trigger is not working */
  if (radar_get_data_acq_mode() != RADAR_MANUAL_TRIGGER)
  {
    return (RADAR_ERR_BUSY);
  }
  
  /* Check if configuration is valid */
  /* ------------------------------- */
  if((configuration->frequency_kHz < BSP_MIN_RF_FREQUENCY_KHZ) ||
     (configuration->frequency_kHz > BSP_MAX_RF_FREQUENCY_KHZ))
  {
    return (RADAR_ERR_FREQUENCY_OUT_OF_RANGE);
  }
  
  if (configuration->tx_power > BSP_MAX_TX_POWER_LEVEL)
  {
    return (RADAR_ERR_POWER_OUT_OF_RANGE);
  }
  
  /* Copy the new Doppler configuration to the data store*/
  ds_ep_store_doppler_config(configuration);
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_doppler_configuration(Radar_Handle_t device, Doppler_Configuration_t* configuration)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  /* Get the current Doppler configuration from data store */
  ds_ep_fetch_doppler_config(configuration);
  return (RADAR_ERR_OK);
}

/* --- End of File -------------------------------------------------------- */
