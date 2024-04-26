/**
   @file radar_target_detection.c

   @brief This file contains the driver API for Target Detection of
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

#include "ds_radar_interface.h"
#include "radar_api.h"

/*
==============================================================================
   2. EXPORTED FUNCTIONS
==============================================================================
 */

uint16_t radar_set_dsp_settings(Radar_Handle_t device, const DSP_Settings_t* configuration)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  algo_settings_t new_settings;
  memset(&new_settings, 0, sizeof(algo_settings_t));
  
  /* Fetch the complete existing data, to update only partially the new data. */
  ds_algo_fetch_settings(&new_settings);
  
#ifdef FMCW_SUPPORTED
  new_settings.tracking_enable 		  	 = configuration->enable_tracking;
  new_settings.mti_filter_enable	 	 = configuration->enable_mti_filter;
  new_settings.num_of_tracks			 = configuration->num_of_tracks;
  new_settings.mti_filter_len			 = configuration->mti_filter_length;
  new_settings.median_filter_len		 = configuration->median_filter_length;
  new_settings.mvg_avg_len				 = configuration->range_mvg_avg_length;
  new_settings.min_distance_cm			 = configuration->min_range_cm;
  new_settings.max_distance_cm			 = configuration->max_range_cm;
  new_settings.range_detection_threshold = configuration->range_threshold;
#endif
  
#if defined(DOPPLER_SUPPORTED) || defined(FMCW_SUPPORTED)
  new_settings.min_speed_kmh = configuration->min_speed_kmh;
  new_settings.max_speed_kmh = configuration->max_speed_kmh;
  new_settings.speed_detection_threshold = configuration->speed_threshold;
#endif
  
  /* Check for new settings */
  if (ds_algo_check_settings(&new_settings) == 0)
  {
    /* Apply new settings */
    ds_algo_store_settings(&new_settings);
    return (RADAR_ERR_OK);
  }
  else
  {
    return (RADAR_ERR_PARAMETER_OUT_OF_RANGE);
  }
}

//============================================================================

uint16_t radar_get_dsp_settings(Radar_Handle_t device, DSP_Settings_t* configuration)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  algo_settings_t *p_ds_config = ds_algo_get_settings();
  
  configuration->enable_tracking      = p_ds_config->tracking_enable;
  configuration->enable_mti_filter    = p_ds_config->mti_filter_enable;
  configuration->num_of_tracks        = p_ds_config->num_of_tracks;
  configuration->mti_filter_length    = p_ds_config->mti_filter_len;
  configuration->median_filter_length = p_ds_config->median_filter_len;
  configuration->range_mvg_avg_length = p_ds_config->mvg_avg_len;
  configuration->min_range_cm         = p_ds_config->min_distance_cm;
  configuration->max_range_cm         = p_ds_config->max_distance_cm;
  configuration->range_threshold      = p_ds_config->range_detection_threshold;
  
  configuration->min_speed_kmh        = p_ds_config->min_speed_kmh;
  configuration->max_speed_kmh        = p_ds_config->max_speed_kmh;
  configuration->speed_threshold      = p_ds_config->speed_detection_threshold;
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_target_info(Radar_Handle_t device, Target_Info_t* target_info, uint8_t* target_count)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */

  ds_ep_fetch_target_list(target_info, target_count);
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_range_detection_threshold(Radar_Handle_t device, uint16_t *threshold)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */

  ds_ep_fetch_range_detection_threshold(threshold);
  return (RADAR_ERR_OK);
}

/* --- End of File -------------------------------------------------------- */
