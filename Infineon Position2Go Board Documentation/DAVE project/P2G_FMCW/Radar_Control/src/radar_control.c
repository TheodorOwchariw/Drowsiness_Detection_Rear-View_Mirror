/**
    @file: radar_control.c

    @brief: This file contains common API definition valid for all radar 
            operations related to MCU, BGT, PLL.
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
   2. DATA
==============================================================================
 */

volatile  uint32_t    g_adc_sampling_completed = false;

volatile  uint32_t    g_do_processing = true;

Radar_Data_Acq_Mode_t  s_data_acq_mode = RADAR_MANUAL_TRIGGER;

extern Radar_Handle_t  h_radar_device;

/*
==============================================================================
    3. EXPORTED FUNCTIONS
==============================================================================
 */

uint16_t radar_init(Radar_Handle_t device)
{
  uint16_t status = RADAR_ERR_BUSY;
  uint16_t num_of_bytes;
  void* calibration_data_ptr = NULL;
  
  /* PLL parameters setup */
  /* -------------------- */
  float temp_pll_lower_freq_MHz = ds_device_get_settings()->pll_lower_frequency_kHz / 1000.0f;
  float temp_pll_upper_freq_MHz = ds_device_get_settings()->pll_upper_frequency_kHz / 1000.0f;
  
  /* Set PLL upper and lower frequencies for FMCW modulation */
  pll_set_upper_lower_frequency(temp_pll_lower_freq_MHz, temp_pll_upper_freq_MHz);
  
#if FW_MODULATION_TYPE == 1U
  
  ds_device_get_settings()->pll_modulation_mode = MODULATION_FMCW;
  status = radar_set_pll_frequency(device, temp_pll_upper_freq_MHz);
  
#else
  
  ds_device_get_settings()->pll_modulation_mode = MODULATION_DOPPLER;
  status = radar_set_pll_frequency(device, BSP_DOPPLER_BASE_FREQ_MHZ);
  
#endif
  
  /* BGT, PLL, PGA, USB and DMA Setup */
  /* -------------------------------- */
  
  acq_buf_obj *p_acq_buf = ds_get_active_acq_buf();
  
  bsp_dma_set_destination_addr_from_acq_buf_obj(p_acq_buf, 0);
  
  status |= (uint16_t) bsp_init();
  
  /* Load calibration data & Timers setup */
  /* ------------------------------------ */
  radar_init_calibration(device);
  
  /* Read ADC calibration data */
  status |= radar_read_calibration(device, CALIBRATION_TARGET_FLASH, CALIBRATION_DATA_ADC, &calibration_data_ptr, &num_of_bytes);
  
  /* Read algorithm calibration data */
  status |= radar_read_calibration(device, CALIBRATION_TARGET_FLASH, CALIBRATION_DATA_ALGO, &calibration_data_ptr, &num_of_bytes);
  
  status |= radar_set_sampling_freq(device, ds_device_get_settings()->pll_modulation_mode, ds_device_get_settings()->adc_sampling_freq_Hz);
  
  return (status);
}

//============================================================================

uint16_t radar_start_acquisition(Radar_Handle_t device)
{

  uint16_t status = RADAR_ERR_BUSY;
  
  if (ds_device_get_settings()->frame_period_usec > 0)
  {
    status = radar_set_automatic_frame_trigger(device, ds_device_get_settings()->frame_period_usec);
  }
  
  return (status);
}

//============================================================================

uint16_t radar_stop_acquisition(Radar_Handle_t device)
{

  uint16_t status = RADAR_ERR_BUSY;
  
  while(g_adc_sampling_completed != true);
  
  status = radar_set_automatic_frame_trigger(device , 0);
  
  return (status);
}

//============================================================================

uint16_t radar_set_bandwidth(Radar_Handle_t device, const float bandwidth_MHz)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  if ((bandwidth_MHz > 0) && (bandwidth_MHz <= BSP_MAX_BANDWIDTH_MHZ))
  {
    ds_device_get_settings()->pll_bandwidth_MHz = bandwidth_MHz;
    
    pll_set_bandwidth(bandwidth_MHz);
    pll_set_update_config_flag(true);
    
    return (RADAR_ERR_OK);
  }
  else
  {
    return (RADAR_ERR_PARAMETER_OUT_OF_RANGE);
  }
}

//============================================================================

void radar_get_bandwidth(Radar_Handle_t device, float* bandwidth_MHz)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */

  *bandwidth_MHz = ds_device_get_settings()->pll_bandwidth_MHz;
}

//============================================================================
uint16_t radar_set_pll_frequency(Radar_Handle_t device, float freq_MHz)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  if ((freq_MHz >= BSP_MIN_RF_FREQUENCY_KHZ/1000.0f) && (freq_MHz <= BSP_MAX_RF_FREQUENCY_KHZ/1000.0f))
  {
    pll_set_frequency(freq_MHz);
    pll_set_update_config_flag(true);
    
    return (RADAR_ERR_OK);
  }
  else
  {
    return (RADAR_ERR_PARAMETER_OUT_OF_RANGE);
  }
}

//============================================================================

void radar_get_pll_frequency(Radar_Handle_t device, float* freq_MHz)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  *freq_MHz = pll_get_frequency();
}

//============================================================================

static uint16_t radar_check_chirp_time(const uint32_t chirp_time_usec)
{
  if ((ds_device_get_settings()->pll_modulation_mode == MODULATION_FMCW) &&
      (chirp_time_usec >= BSP_MIN_CHIRP_TIME_USEC) && (chirp_time_usec <= BSP_MAX_CHIRP_TIME_USEC))
  {
    if ((chirp_time_usec < 300 && ds_device_get_settings()->num_samples_per_chirp > 256) ||
        (chirp_time_usec < 150 && ds_device_get_settings()->num_samples_per_chirp > 128) ||
        (chirp_time_usec < 75  && ds_device_get_settings()->num_samples_per_chirp > 64))
    {
      return (RADAR_ERR_PARAMETER_OUT_OF_RANGE);
    }
  }
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_set_chirp_time(Radar_Handle_t device, const uint32_t chirp_time_usec)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  uint16_t retValue = radar_check_chirp_time(chirp_time_usec);
  
  if(retValue)
  {
    return retValue;
  }
  
  /* ---- chirp time was verified being o.k. */
  ds_device_get_settings()->pll_chirp_time_usec = chirp_time_usec;
  
  /* Set the current PLL chirp-time, in units of usec */
  pll_set_chirp_time(ds_device_get_settings()->pll_chirp_time_usec);
   
  bsp_timer_stop_clear(&TIMER_ADC_TRIG);
  
  uint32_t sampling_period_usec = (uint32_t)(((float)ds_device_get_settings()->pll_chirp_time_usec * 100.0f) / (float)ds_device_get_settings()->num_samples_per_chirp);
  
  if (TIMER_SetTimeInterval(&TIMER_ADC_TRIG, sampling_period_usec) != TIMER_STATUS_SUCCESS)
  {
    return (RADAR_ERR_PARAMETER_OUT_OF_RANGE);
  }
  
  /* Update the PLL Configuration */
  pll_set_update_config_flag(true);
  
  return (RADAR_ERR_OK);
}

//============================================================================

void radar_get_chirp_time(Radar_Handle_t device, uint32_t* chirp_time_usec)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  *chirp_time_usec = ds_device_get_settings()->pll_chirp_time_usec;
}

//============================================================================

void radar_get_frame_period(Radar_Handle_t device, uint32_t* time_msec)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  *time_msec = (ds_device_get_settings()->frame_period_usec / 1000);
}

//============================================================================

uint16_t radar_set_sampling_freq(Radar_Handle_t device, Modulation_Type_t modulation_type, const float freq_Hz)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  TIMER_STATUS_t status = TIMER_STATUS_FAILURE;
  uint32_t sampling_period_usec = 0;
  
  if (((uint32_t)freq_Hz <= 0) || ((uint32_t)freq_Hz > BSP_MAX_ADC_FREQ_HZ))
  {
    return (RADAR_ERR_BUSY);
  }
  
  bsp_timer_stop_clear(&TIMER_ADC_TRIG);
  
  if (modulation_type == MODULATION_FMCW)
  {
    sampling_period_usec = (uint32_t)(((float) ds_device_get_settings()->pll_chirp_time_usec * 100.0f) / (float)ds_device_get_settings()->num_samples_per_chirp);
  }
  else
  {
    sampling_period_usec = (uint32_t)(1000.0f * 1000.0f * 100.0f / freq_Hz);
  }
  
  status = TIMER_SetTimeInterval(&TIMER_ADC_TRIG, sampling_period_usec);
  
  if (status == TIMER_STATUS_SUCCESS)
  {
    return (RADAR_ERR_OK);
  }
  else
  {
    return (RADAR_ERR_BUSY);
  }
}

//============================================================================

void radar_get_sampling_freq(Radar_Handle_t device, float* freq_Hz)
{
  XMC_UNUSED_ARG(device); ; /* Suppress compiler warnings of unused variable */
  
  *freq_Hz = (1000.0f * 1000.0f) * 100.0f / ds_device_get_settings()->adc_sampling_freq_Hz;
}

//============================================================================

uint16_t radar_register_callback(Radar_Callback_ID_t cb_id, void (*pcallback)(void *parameter))
{
  uint16_t status = RADAR_ERR_OK;
  
  /* Check if the callback function is valid */
  if(pcallback == NULL)
  {
    return(RADAR_ERR_BUSY);
  }
  
  /* Register user callback function based on callback ID */
  switch (cb_id)
  {
  case RADAR_ACQUISITION_STARTED_CB_ID:
    /* Callback executed by BSP layer at the beginning of data acquisition */
    bsp_register_callback(BSP_ACQUISITION_STARTED_CB_ID, pcallback);
    break;
    
  case RADAR_ACQUISITION_DONE_CB_ID:
    /* Callback executed by BSP layer at the end of data acquisition process */
    bsp_register_callback(BSP_ACQUISITION_DONE_CB_ID, pcallback);
    break;
    
  default:
    /* Return error status */
    status =  RADAR_ERR_BUSY;
    break;
  }
  
  return (status);
}

//============================================================================

uint16_t radar_unregister_callback(Radar_Callback_ID_t cb_id)
{
  uint16_t status = RADAR_ERR_OK;
  
  /* Unregister user callback function based on callback ID */
  switch (cb_id)
  {
  case RADAR_ACQUISITION_STARTED_CB_ID:
    bsp_unregister_callback(BSP_ACQUISITION_STARTED_CB_ID);
    break;
    
  case RADAR_ACQUISITION_DONE_CB_ID:
    bsp_unregister_callback(BSP_ACQUISITION_DONE_CB_ID);
    break;
    
  default:
    /* Return error status */
    status = RADAR_ERR_BUSY;
    break;
  }
  
  return (status);
}

//============================================================================

static uint32_t radar_apply_new_hw_setting(Radar_Handle_t device, device_settings_t *p_new)
{
  bsp_components_power_down();
  
  device_settings_t *p_curr = ds_device_get_shadow_settings();
  
  /* Set BGT RX LNA gain */
  if (p_new->bgt_rx_lna_gain_enable_flag == true)
  {
    bgt_lna_gain_enable();
  }
  else
  {
    bgt_lna_gain_disable();
  }
  
  if (p_curr->pll_chirp_time_usec != p_new->pll_chirp_time_usec || (p_curr->num_samples_per_chirp != p_new->num_samples_per_chirp))
  {
    if (ds_device_get_settings()->pll_modulation_mode == MODULATION_FMCW)
	{
      pll_set_chirp_time(p_new->pll_chirp_time_usec);
      pll_set_update_config_flag(true);
    
      uint32_t sampling_period_usec = (uint32_t)(((float) p_new->pll_chirp_time_usec * 100.0f) / (float)p_new->num_samples_per_chirp);

      if (TIMER_SetTimeInterval(&TIMER_ADC_TRIG, sampling_period_usec) != TIMER_STATUS_SUCCESS)
      {
        return (RADAR_ERR_PARAMETER_OUT_OF_RANGE);
      }
	}
  }
  
  if (p_curr->pll_num_of_chirps_per_frame != p_new->pll_num_of_chirps_per_frame)
  {
    /* Update the number of chirps count per frame */
    bsp_set_num_chirps_per_frame(p_new->pll_num_of_chirps_per_frame);
    pll_set_update_config_flag(true);
  }
  
  if (p_curr->num_samples_per_chirp != p_new->num_samples_per_chirp)
  {
    bsp_set_num_samples_per_dma_transfer(p_new->num_samples_per_chirp);
  }
  
  if (p_new->isUpdated_doppler_config)
  {
    float freq_MHz;
    
    /* Apply RF configuration */
    /* ---------------------- */
    bgt_set_tx_power(p_new->bgt_tx_power_level);
    
    freq_MHz = (float)(p_new->pll_frequency_kHz / 1000.0f); // to bring into the MHz
    
    if (ds_device_get_settings()->pll_modulation_mode == MODULATION_DOPPLER)
    {
      radar_set_pll_frequency(device, freq_MHz);
    }
    
    p_new->isUpdated_doppler_config = 0;
  }
  
  if(p_new->isUpdated_fmcw_config)
  {
    float bandwidth_mhz;
    float pll_base_frequency_mhz;

    /* Apply RF configuration */
    /* ---------------------- */
    bandwidth_mhz = ((float)(p_new->pll_upper_frequency_kHz - p_new->pll_lower_frequency_kHz)) / 1000.0f;
    
    pll_base_frequency_mhz = (float)(p_new->pll_upper_frequency_kHz) / 1000.0f;
    
    // This is redundant, but we leave it for now, needs to be tested!
    ds_device_get_settings()->pll_bandwidth_MHz = bandwidth_mhz;
    
    pll_set_frequency(pll_base_frequency_mhz);
    pll_set_bandwidth(bandwidth_mhz);
    pll_set_upper_lower_frequency(p_new->pll_lower_frequency_kHz / 1000.0f, p_new->pll_upper_frequency_kHz / 1000.0f);
    pll_set_update_config_flag(true);
    
    bgt_set_tx_power(p_new->bgt_tx_power_level);
    
    p_new->isUpdated_fmcw_config = 0;
  }

  /* Set the Radar duty cycle */
  if(p_new->is_duty_cycle_enable_updated)
  {
    bsp_set_duty_cycle_enable_flag(p_new->power_duty_cycle_enable_flag);
    p_new->is_duty_cycle_enable_updated = 0;
  }

  if ((p_curr->pll_num_of_chirps_per_frame != p_new->pll_num_of_chirps_per_frame) || (p_curr->num_samples_per_chirp != p_new->num_samples_per_chirp))
  {
    /* Recreate acq-buffer object! */
    /* The buffer is statically allocated to the max supported buffer size, so it is save to re-use it!
     * checks should have been applied earlier! */
    acq_buf_obj *p_acq_buf = ds_get_active_acq_buf();
    
    uint8_t *p_temp_buf = p_acq_buf->p_acq_buf;
    uint32_t reserved_size = p_acq_buf->internals.size_of_acq_buffer;

    acq_buf_obj acq_buf = create_acq_data_buffer_obj(p_temp_buf, reserved_size,
                                                     p_new->pll_num_of_chirps_per_frame,
                                                     p_new->num_samples_per_chirp, 2, 2, 12);  // 2 --> 2 I channels and two Q-channels, from two RX-antennas, 12-bit ADC resolution
    
    ds_set_active_acq_buf(acq_buf);
  }
  
  bsp_components_power_up();
  
  return 0; // confirm that everything was o.k.
}

//============================================================================

uint32_t radar_apply_hw_settings(Radar_Handle_t device, hw_state_setting_t state, device_settings_t *p_hw_settings)
{
  uint32_t retValue = RADAR_ERR_OK;
  
  switch(state)
  {
  case INIT:
    retValue = radar_apply_new_hw_setting(device, p_hw_settings);
    break;
    
  case ENTER_POWERSAVING:
    //bsp_... dear down
    break;
    
  case LEAVE_POWERSAVING:
    break;
    
  case POWER_DOWN:
    bgt_stop_tx();
    break;
    
  case UPDATE:
    if(ds_is_device_settings_updated())
    {
      /* stop the timer for ADC acquisition, while updating the values! */
      bsp_timer_stop_clear(&TIMER_ADC_TRIG);
      
      retValue = radar_apply_new_hw_setting(device, p_hw_settings);
    }
    break;
    
  case START_MANUAL_ACQ:
    {
      acq_buf_obj *p_acq_buf = ds_get_active_acq_buf();
      
      bsp_dma_set_destination_addr_from_acq_buf_obj(p_acq_buf, 0);
      
      bsp_reset_frame_counter();
      
      g_adc_sampling_completed = false;
    }
    break;
    
  case START_AUTO_ACQ:
    break;
    
  case STOP_ACQ:
    break;
    
  default:
    break;
  }
  
  return retValue;
}

//============================================================================

Radar_Data_Acq_Mode_t radar_get_data_acq_mode(void)
{
  return s_data_acq_mode;
}

//============================================================================

void radar_set_data_acq_mode(Radar_Data_Acq_Mode_t new_mode)
{
  s_data_acq_mode = new_mode;
}

/* --- End of File -------------------------------------------------------- */
