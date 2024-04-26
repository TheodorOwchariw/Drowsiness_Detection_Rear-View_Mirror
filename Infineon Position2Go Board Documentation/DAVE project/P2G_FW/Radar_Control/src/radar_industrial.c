/**
    @file: radar_industrial.c

    @brief: This file contains the driver API for Industrial operation of
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
#include "ds_radar_interface.h"

/*
==============================================================================
   2. EXPORTED FUNCTIONS
==============================================================================
 */

uint16_t radar_set_duty_cycle(Radar_Handle_t device, const uint8_t flag)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  if (flag != ds_device_get_settings()->power_duty_cycle_enable_flag)
  {
    ds_device_get_settings()->is_duty_cycle_enable_updated = 1;
    ds_device_get_settings()->power_duty_cycle_enable_flag = flag;
  }
  
  return (RADAR_ERR_OK);
}

//============================================================================

void radar_get_duty_cycle(Radar_Handle_t device, uint8_t* flag)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  *flag = ds_device_get_settings()->power_duty_cycle_enable_flag;
}

//============================================================================

void radar_disable_lna_gain(Radar_Handle_t device)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  ds_device_get_settings()->bgt_rx_lna_gain_enable_flag = false;
}

//============================================================================

void radar_enable_lna_gain(Radar_Handle_t device)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  ds_device_get_settings()->bgt_rx_lna_gain_enable_flag = true;
}

//============================================================================

uint8_t radar_get_lna_gain_enable_status(Radar_Handle_t device)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  return (ds_device_get_settings()->bgt_rx_lna_gain_enable_flag);
}

/* --- End of File -------------------------------------------------------- */
