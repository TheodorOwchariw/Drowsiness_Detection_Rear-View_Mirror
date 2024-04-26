/**
    @file  ds_radar_interface.c

  	@brief This file includes the function to access data structure from data store.
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
#include "version.h"

/*
==============================================================================
   2. DATA
==============================================================================
 */

/* ***************** Get data structures from datastore ******************* */

extern float result_range_detection_threshold;

/* Global constant structure to store the device info, this only holds compile time defines, so no need to store it in RAM! */
const Device_Info_t cg_device_info_struct = {
  .description           = DEMO_DESCRIPTION,
  .major_version_hw      = HW_VERSION_MAJOR,
  .minor_version_hw      = HW_VERSION_MINOR,
  .num_tx_antennas       = BSP_NUM_TX_ANTENNAS,
  .num_rx_antennas       = BSP_NUM_RX_ANTENNAS,
  .num_temp_sensors      = BSP_NUM_TEMP_SENSORS,
  .max_tx_power          = BSP_MAX_TX_POWER_LEVEL,
  .min_rf_frequency_kHz  = BSP_MIN_RF_FREQUENCY_KHZ,
  .max_rf_frequency_kHz  = BSP_MAX_RF_FREQUENCY_KHZ,
  .interleaved_rx        = 0,
  .data_format 	         = RADAR_RX_DATA_COMPLEX
};

/*
==============================================================================
   3. EXPORTED FUNCTIONS
==============================================================================
 */

void ds_ep_fetch_adc_config(Adc_Configuration_t *ptr)
{
	device_settings_t *p_dev   = ds_device_get_settings();
	ptr->samplerate_Hz 		   = p_dev->adc_sampling_freq_Hz;
	ptr->resolution    		   = p_dev->adc_resolution;
	ptr->use_post_calibration  = p_dev->adc_use_post_calibration;
}

//============================================================================

uint16_t ds_ep_store_adc_config(const Adc_Configuration_t *p_config)
{
  uint16_t status = RADAR_ERR_BUSY;

  if(p_config != NULL)
  {
	  device_settings_t *p_dev		   = ds_device_get_settings();
	  p_dev->adc_sampling_freq_Hz 	   = p_config->samplerate_Hz;
	  p_dev->adc_resolution    		   = p_config->resolution;
	  p_dev->adc_use_post_calibration  = p_config->use_post_calibration;

	  status = RADAR_ERR_OK;
  }
  return (status);
}

//============================================================================

void ds_ep_fetch_doppler_config(Doppler_Configuration_t *config)
{
	device_settings_t *p_dev  = ds_device_get_settings();
    config->frequency_kHz 	  = p_dev->pll_frequency_kHz;
    config->tx_power 		  = p_dev->bgt_tx_power_level;
}

//============================================================================

uint16_t ds_ep_store_doppler_config(const Doppler_Configuration_t *config)
{
  uint16_t status = RADAR_ERR_BUSY;

  if(config != NULL)
  {
	  device_settings_t *p_dev		   = ds_device_get_settings();
	  p_dev->pll_frequency_kHz 		   = config->frequency_kHz;
	  p_dev->bgt_tx_power_level 	   = config->tx_power;
	  p_dev->isUpdated_doppler_config  = 1;

	  status = RADAR_ERR_OK;
  }
  return (status);
}

//============================================================================

void ds_ep_fetch_fmcw_config(Fmcw_Configuration_t *config)
{
  device_settings_t *p_dev	   		= ds_device_get_settings();
  config->lower_frequency_kHz  		= p_dev->pll_lower_frequency_kHz;
  config->upper_frequency_kHz  		= p_dev->pll_upper_frequency_kHz;
  config->direction 		   		= RADAR_DIR_UPCHIRP_ONLY;
  config->tx_power 			   		= p_dev->bgt_tx_power_level;
  config->up_chirp_duration_nsec 	= p_dev->pll_chirp_time_usec * 1000U;
  config->down_chirp_duration_nsec  = (PLL_RAMP_DOWM_TIME_USEC) * 1000U;
  config->chirp_to_chirp_delay_nsec = (PLL_STEADY_STATE_USEC)* 1000U;
}

//============================================================================

uint16_t ds_ep_store_fmcw_config(const Fmcw_Configuration_t *config)
{
  uint16_t status = RADAR_ERR_BUSY;

  if(config != NULL)
  {
	  if(config->direction != RADAR_DIR_UPCHIRP_ONLY)
	  {
		  status = RADAR_ERR_UNSUPPORTED_DIRECTION;
	  }
	  else
	  {
		  device_settings_t *p_dev		  = ds_device_get_settings();
		  p_dev->pll_lower_frequency_kHz  = config->lower_frequency_kHz;
		  p_dev->pll_upper_frequency_kHz  = config->upper_frequency_kHz;
		  p_dev->bgt_tx_power_level 	  = config->tx_power;
		  p_dev->isUpdated_fmcw_config 	  = 1;

		  status = RADAR_ERR_OK;
	  }
  }
  return (status);
}

//============================================================================

const Device_Info_t *ds_ep_get_device_info(void)
{
  return (&cg_device_info_struct);
}

//============================================================================

void ds_ep_fetch_frame_format(Frame_Format_t *p_frame_format)
{
	device_settings_t *p_dev			  = ds_device_get_settings();
	p_frame_format->num_samples_per_chirp = p_dev->num_samples_per_chirp;
	p_frame_format->num_chirps_per_frame  = p_dev->pll_num_of_chirps_per_frame;
	p_frame_format->rx_mask 		      = p_dev->rx_antenna_mask;
	p_frame_format->signal_part 	      = RADAR_SIGNAL_I_AND_Q;
}

//============================================================================

uint16_t ds_ep_store_frame_format(const Frame_Format_t *p_new_frame_format)
{
  uint16_t status = RADAR_ERR_BUSY;

  if(p_new_frame_format != NULL)
  {
	  if(p_new_frame_format->signal_part != RADAR_SIGNAL_I_AND_Q)
	  {
		  status = RADAR_ERR_UNAVAILABLE_SIGNAL_PART;
	  }
	  else
	  {
		  device_settings_t *p_dev				= ds_device_get_settings();
		  p_dev->num_samples_per_chirp 	        = p_new_frame_format->num_samples_per_chirp;
		  p_dev->pll_num_of_chirps_per_frame 	= p_new_frame_format->num_chirps_per_frame;
		  p_dev->rx_antenna_mask	            = p_new_frame_format->rx_mask;

		  status = RADAR_ERR_OK;
	  }
  }
  return (status);
}

//============================================================================

void ds_ep_fetch_target_list(Target_Info_t* target_info_ptr, uint8_t* target_count)
{
  uint8_t num_targets = 0;

  if (ds_device_get_settings()->pll_modulation_mode == MODULATION_FMCW)
  {
	ds_target_list_t *p_target_list = ds_get_target_list();
    for (uint8_t idx = 0; idx < p_target_list->max_num_of_targets; idx++)
    {
      if (p_target_list->elems[idx].target_id != 0)
      {
        target_info_ptr[num_targets].target_id 	   = idx + 1;
    	target_info_ptr[num_targets].level 		   = p_target_list->elems[idx].strength;
        target_info_ptr[num_targets].radius 	   = p_target_list->elems[idx].range;
        target_info_ptr[num_targets].radial_speed  = p_target_list->elems[idx].speed * 3.6f; // m/s => km/h
        target_info_ptr[num_targets].azimuth 	   = p_target_list->elems[idx].angle;

        num_targets++;
      }
    }
  }

  *target_count = num_targets;
}

//============================================================================

void ds_ep_fetch_range_detection_threshold(uint16_t* threshold)
{
  uint16_t temp = (uint16_t) result_range_detection_threshold;

  *threshold = temp;
}

void ds_ep_fetch_fft_config(FFT_configuration_t *p_fft_config)
{
	p_fft_config->fft_window_size 	= RANGE_FFT_SIZE;
	p_fft_config->mean_removed 		= 0;
	p_fft_config->window_type 		= WINDOW_BLACKMANHARRIS;
	p_fft_config->scaling_factor	= 1.0F;
}

/* --- End of File -------------------------------------------------------- */
