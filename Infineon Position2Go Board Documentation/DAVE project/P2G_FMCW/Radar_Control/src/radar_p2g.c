/**
  @file radar_p2g.c

  @brief This file contains the driver API for Position2Go Board of
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
** is not permitted without Infineon's prior written authorisation.
** ===========================================================================
*/

/*
==============================================================================
   1. INCLUDE FILES
==============================================================================
 */

#include "radar_api.h"
#include "ds_radar_interface.h"

/*
==============================================================================
   2. EXPORTED FUNCTIONS
==============================================================================
 */

uint16_t radar_set_gain_level(Radar_Handle_t device, uint16_t gain_level)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  if (gain_level > MAX_PGA_GAIN_LEVEL)
  {
    return (RADAR_ERR_UNSUPPORTED_PGA_GAIN);
  }
  
  /* Update the PGA Gain level */
  ds_device_get_settings()->pga_rx_gain_level = gain_level;
  ds_device_get_settings()->isGainlevelUpdated = 1;
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_gain_level(Radar_Handle_t device, uint16_t* gain_level)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  /* Get current PGA Gain level from driver */
  *gain_level = ds_device_get_settings()->pga_rx_gain_level;
  
  return (RADAR_ERR_OK);
}

/* --- End of File -------------------------------------------------------- */
