/**
    @file: radar_base.c

    @brief: This file contains the basic driver API for Infineon radar sensors.
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

#include "version.h"
#include "radar_api.h"
#include "radar_control.h"
#include "ds_radar_interface.h"
#include "position2go_bsp.h"

/*
==============================================================================
   2. DATA
==============================================================================
 */

size_t frame_size = 0;

/* Global Data Buffer used to store the raw IQ data samples for FMCW chirps,
   first (SAMPLES_PER_CHIRP * 2) samples are for I&Q data for RX1,
   followed by (SAMPLES_PER_CHIRP * 2) samples for I&Q data for RX2 */
uint32_t data_fifo[BSP_MAX_ADC_BUFFER_SIZE_BYTES/sizeof(uint32_t)];  
                                   					 	 	 	 	 	

static const Driver_Version_t driver_version_s =
{
  /*.uMajor    =*/ FW_VERSION_MAJOR,
  /*.uMinor    =*/ FW_VERSION_MINOR,
  /*.uRevision =*/ FW_VERSION_REVISION
};

/*
==============================================================================
   3. EXPORTED FUNCTIONS
==============================================================================
 */

const Driver_Version_t* radar_get_driver_version(void)
{
  return (&driver_version_s);
}

//============================================================================

uint8_t radar_get_number_of_available_devices(void)
{
  return bsp_get_number_of_available_radar_devices();
}

//============================================================================

uint16_t radar_open_device(int32_t device_number, Radar_Handle_t* device)
{
  XMC_UNUSED_ARG(device_number); /* Suppress compiler warnings of unused variable */
  
  Frame_Format_t frame_format;
  uint16_t retValue = RADAR_ERR_OK;
  int32_t status = 0;
  uint16_t num_of_bytes;
  void* calibration_data_ptr = NULL;
  acq_buf_obj acq_buf = { .status = not_initialized } ;
    
  /* Initialize the driver instance with NULL, until it was successfully initialized */
  *device = NULL;
  
  /* Initialize driver data structure */
  /* -------------------------------- */
  status = ds_init();
  status |= radar_read_calibration(device, CALIBRATION_TARGET_FLASH, CALIBRATION_DATA_ADC, &calibration_data_ptr, &num_of_bytes);
  if(status > 0)
  {
    return (RADAR_ERR_UNSUPPORTED_FRAME_FORMAT);
  }
  
  //========================================================================
  //------------------------- Frame Format configurations ------------------
  //========================================================================
  ds_ep_fetch_frame_format(&frame_format);
  
  //=========================================================================
  //------------------------- General configurations ------------------------
  //=========================================================================
  radar_set_data_acq_mode(RADAR_MANUAL_TRIGGER);
  
  /* If the frame format contains a 0, this makes no sense. */
  if ((frame_format.rx_mask == 0) ||
      (frame_format.num_chirps_per_frame  == 0) ||
        (frame_format.num_samples_per_chirp == 0) ||
          (frame_format.num_chirps_per_frame  > BSP_MAX_NUM_CHIRPS_PER_FRAME) ||
            (frame_format.num_samples_per_chirp > BSP_MAX_NUM_SAMPLES_PER_CHIRP))
  {
    return (RADAR_ERR_UNSUPPORTED_FRAME_FORMAT);
  }
  
  /* Check frame memory limits */
  /* ------------------------- */
  frame_size = frame_format.num_chirps_per_frame * frame_format.num_samples_per_chirp *
    ds_ep_get_device_info()->num_rx_antennas * 2U * sizeof(uint16_t); // complex IQ samples per antenna
  
  if (frame_size > BSP_MAX_ADC_BUFFER_SIZE_BYTES)
  {
    return (RADAR_ERR_UNSUPPORTED_FRAME_FORMAT);
  }
  
  /* Initialize sensor hardware */
  /* -------------------------- */
  acq_buf = create_acq_data_buffer_obj((uint8_t*) data_fifo, BSP_MAX_ADC_BUFFER_SIZE_BYTES,
                                       frame_format.num_chirps_per_frame,
                                       frame_format.num_samples_per_chirp, 2, 2, 12);
  
  ds_set_active_acq_buf(acq_buf);
  
  ds_device_get_settings()->num_samples_per_chirp = frame_format.num_samples_per_chirp;
  
  /* Set the number of samples to captured during each chirp
  This value will be used in BSP to configure the DMA block size and destination address */
  bsp_set_num_samples_per_dma_transfer(ds_device_get_settings()->num_samples_per_chirp);
  
  ds_device_get_settings()->pll_num_of_chirps_per_frame = frame_format.num_chirps_per_frame;
  
  /* Set the number of chirps count per frame */
  bsp_set_num_chirps_per_frame(frame_format.num_chirps_per_frame);
  
  /* BGT, PLL, DMA, Timers & Calibration initialization */
  /* -------------------------------------------------- */
  retValue = radar_init(device);
  
  g_adc_sampling_completed = false;
  
  return (retValue);
}

//============================================================================

void radar_close_device(Radar_Handle_t device)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  bgt_stop_tx();
}

//============================================================================

uint16_t radar_get_device_info(Radar_Handle_t device, Device_Info_t* device_info)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */

  /* Get device information:
      - Sensor description
      - Hardware Version information
      - Sensor features
      - Data format
  */
  *device_info = *ds_ep_get_device_info();
   
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_temperature(Radar_Handle_t device, uint8_t temp_sensor, int32_t* temperature_001C)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  /* Check temperature sensor, only single sensor supported */
  if (temp_sensor == 0)
  {
    *temperature_001C = (int32_t) (bsp_bgt_get_temp_value() * 1000.0f);	// (expected value is actual temperature in C° x 1000)
    
    return (RADAR_ERR_OK);
  }
  else
  {
    return (RADAR_ERR_SENSOR_DOES_NOT_EXIST);
  }
}

//============================================================================

uint16_t radar_get_tx_power(Radar_Handle_t device, uint8_t tx_antenna, int32_t* tx_power_001dBm)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  /* Check number of TX antenna */
  if (tx_antenna == 0)
  {
    *tx_power_001dBm = (int32_t)(bsp_bgt_get_txpower_value() * 1000.0f); // return TX power
    
    return (RADAR_ERR_OK);
  }
  
  return (RADAR_ERR_ANTENNA_DOES_NOT_EXIST);
}

//============================================================================

uint16_t radar_get_chirp_duration(Radar_Handle_t device, uint32_t* chirp_duration_nsec)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  /* Calculate chirp timing from current settings */
  /* -------------------------------------------- */
  /* Unit is ns. (factor 10^9 is distributed to denominator and denominator for fixed point precision) */
  
  *chirp_duration_nsec = (ds_device_get_settings()->pll_chirp_time_usec * 1000U);
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_min_frame_interval(Radar_Handle_t device, uint32_t* min_frame_interval_usec)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  uint32_t bgt_duty_cycle_delay_usec;
  
  if(ds_device_get_settings()->power_duty_cycle_enable_flag == true)
  {
    bgt_duty_cycle_delay_usec = 3000U; /* approximate value calculated from bgt power up sequence 100 + 400 + 100 + update configuration time + PLL_LOCK_TIME_us + BGT_VCO_LOCK_TIME_us */
  }
  else
  {
    bgt_duty_cycle_delay_usec = 1000U; /* approximate delay to accommodate processing time */
  }
  
#if defined (FMCW_SUPPORTED) && defined (DOPPLER_SUPPORTED)
  
  *min_frame_interval_usec = ds_device_get_settings()->pll_num_of_chirps_per_frame * (ds_device_get_settings()->pll_chirp_time_usec + PLL_RAMP_DOWM_TIME_USEC + PLL_STEADY_STATE_USEC) + bgt_duty_cycle_delay_usec +
    ((SAMPLES_PER_CHIRP * 1000U / DOPPLER_SAMPLING_FREQ_HZ) * 1000U) + 2 * ALGO_PROCESS_TIME_USEC; // units in us
  
#elif defined (FMCW_SUPPORTED)
  *min_frame_interval_usec = ds_device_get_settings()->pll_num_of_chirps_per_frame * (ds_device_get_settings()->pll_chirp_time_usec + PLL_RAMP_DOWM_TIME_USEC + PLL_STEADY_STATE_USEC) + bgt_duty_cycle_delay_usec + ALGO_PROCESS_TIME_USEC; // units in us
  
#elif defined (DOPPLER_SUPPORTED)
  *min_frame_interval_usec = ds_device_get_settings()->pll_num_of_chirps_per_frame * ((ds_device_get_settings()->num_samples_per_chirp * 1000U * 1000U) / DOPPLER_SAMPLING_FREQ_HZ) + bgt_duty_cycle_delay_usec + ALGO_PROCESS_TIME_USEC; // units in us
  
#else
  
  if(ds_device_get_settings()->pll_modulation_mode == MODULATION_FMCW)
  {
    *min_frame_interval_usec = ds_device_get_settings()->pll_num_of_chirps_per_frame * (ds_device_get_settings()->pll_chirp_time_usec + PLL_RAMP_DOWM_TIME_USEC + PLL_STEADY_STATE_USEC) + bgt_duty_cycle_delay_usec; // units in us
  }
  else /* Doppler modulation */
  {
    /* min_frame_interval_usec = Sampling period x No. of chirps x No. of samples per chirp x 1000000 (for microseconds)
    Sampling period =  1 / DOPPLER_SAMPLING_FREQ_Hz */
    *min_frame_interval_usec = ds_device_get_settings()->pll_num_of_chirps_per_frame * ((ds_device_get_settings()->num_samples_per_chirp * 1000U * 1000U) / DOPPLER_SAMPLING_FREQ_HZ) + bgt_duty_cycle_delay_usec; // units in us
  }
  
#endif
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_set_frame_format(Radar_Handle_t device, const Frame_Format_t* frame_format)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  /* Configuration change is only allowed while automatic trigger is not working */
  if (radar_get_data_acq_mode() != RADAR_MANUAL_TRIGGER)
  {
    return (RADAR_ERR_BUSY);
  }
  
  /* Check if configuration is valid */
  /* ------------------------------- */
  
  /* This sensor can only acquire real signal data */
  if (frame_format->signal_part != RADAR_SIGNAL_I_AND_Q)
  {
    return (RADAR_ERR_UNAVAILABLE_SIGNAL_PART);
  }
  /* This sensor support RX1, RX2 and both at the same time antennas enabled */
  if (frame_format->rx_mask > BSP_RX_ANTENNA_MASK )
  {
    return (RADAR_ERR_ANTENNA_DOES_NOT_EXIST);
  }
  /* LMX2491 supports only a limited number of ramps before auto turn off. */
  if ((frame_format->num_chirps_per_frame == 0) || (frame_format->num_chirps_per_frame > BSP_MAX_NUM_CHIRPS_PER_FRAME))
  {
    return (RADAR_ERR_UNSUPPORTED_FRAME_FORMAT);
  }
  
  /* BGT can store only a limited number of samples per chirp */
  if ((frame_format->num_samples_per_chirp == 0) || (frame_format->num_samples_per_chirp > BSP_MAX_NUM_SAMPLES_PER_CHIRP))
  {
    return (RADAR_ERR_UNSUPPORTED_FRAME_FORMAT);
  }
  
  /* For Doppler modulation mode the number of chirp must be 1 */
  if ((frame_format->num_chirps_per_frame != 1) && (ds_device_get_settings()->pll_modulation_mode == MODULATION_DOPPLER))
  {
    return (RADAR_ERR_UNSUPPORTED_FRAME_FORMAT);
  }
  
  /* Check if samples per chirp are in [32,64,128,256] range */
  if ((frame_format->num_samples_per_chirp != 256U) &&
      (frame_format->num_samples_per_chirp != 128U) &&
        (frame_format->num_samples_per_chirp != 64U)  &&
          (frame_format->num_samples_per_chirp != 32U))
  {
    return (RADAR_ERR_SAMPLERATE_OUT_OF_RANGE);
  }
  
  uint32_t temp_samplerate_Hz = (uint32_t) (((float)frame_format->num_samples_per_chirp * 1000000.0f ) / ds_device_get_settings()->pll_chirp_time_usec );
  radar_set_sampling_freq(device, MODULATION_FMCW, temp_samplerate_Hz);

  /* Check frame memory limits */
  /* ------------------------- */
  frame_size = frame_format->num_chirps_per_frame * frame_format->num_samples_per_chirp * 2U *
    ds_ep_get_device_info()->num_rx_antennas * sizeof(uint16_t); // complex IQ samples per antenna
  
  if (frame_size > BSP_MAX_ADC_BUFFER_SIZE_BYTES)
  {
    return (RADAR_ERR_UNSUPPORTED_FRAME_FORMAT);
  }
  
  /* Store new frame format in store */
  ds_ep_store_frame_format(frame_format);
  
  g_adc_sampling_completed = false;
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_frame_format(Radar_Handle_t device, Frame_Format_t* p_frame_format)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  /* Return current frame format of driver object */
  ds_ep_fetch_frame_format(p_frame_format);
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_get_frame(Radar_Handle_t device, Frame_Info_t* frame_info, uint8_t wait_for_data)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  Frame_Format_t frame_format;
  
  ds_ep_fetch_frame_format(&frame_format);
  
  acq_buf_obj *p_acq_buf = ds_get_active_acq_buf();
  
  if (radar_get_data_acq_mode() != RADAR_AUTO_TRIGGER)	// without automaticFrameTrigger
  {
    bsp_trigger_new_frame();
  }
  else // with automaticFrameTrigger
  {
    if ((!wait_for_data) && (!g_adc_sampling_completed))
    {
      return (RADAR_ERR_TIMEOUT);
    }
  }
  
  while (g_adc_sampling_completed == false); // block GetFrame if data is not ready
  
  /* Data acquired */
  /* ------------- */
  //memcpy(transport_buffer, p_acq_buf->p_acq_buf, frame_size);	// copy ADC buffer to be send to avoid data corruption by Automatic trigger
  
  frame_info->sample_data           = p_acq_buf->p_acq_buf;
  frame_info->num_rx_antennas       = ds_ep_get_device_info()->num_rx_antennas;
  frame_info->num_chirps            = frame_format.num_chirps_per_frame;
  frame_info->rx_mask               = frame_format.rx_mask;
  frame_info->adc_resolution        = ds_device_get_settings()->adc_resolution;
  frame_info->interleaved_rx        = 0;
  frame_info->frame_number          = bsp_get_frame_counter();
  frame_info->data_format           = RADAR_RX_DATA_COMPLEX;
  frame_info->temperature_001C      = (int32_t) (bsp_bgt_get_temp_value() * 1000);	// (expected value is actual temperature in C° x 1000)
  frame_info->num_samples_per_chirp = frame_format.num_samples_per_chirp;

  g_adc_sampling_completed          = false;
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_set_automatic_frame_trigger(Radar_Handle_t device, uint32_t frame_interval_usec)
{
  /* Stop trigger timer first  */
  /* ------------------------- */
  bsp_timer_stop_clear(&TIMER_FRAME_TRIG);
  
  radar_set_data_acq_mode(RADAR_MANUAL_TRIGGER);
  
  /* Restart trigger timer */
  /* --------------------- */
  if (frame_interval_usec > 0)
  {
    uint32_t min_frame_interval_usec;
    
    radar_get_min_frame_interval(device, &min_frame_interval_usec);
    
    if (frame_interval_usec < min_frame_interval_usec)
    {
      return (RADAR_ERR_UNSUPPORTED_FRAME_INTERVAL);
    }
    
    /* Update time period of Automatic Trigger */
    if (TIMER_SetTimeInterval(&TIMER_FRAME_TRIG, frame_interval_usec * 100U) != TIMER_STATUS_SUCCESS)
    {
      return (RADAR_ERR_UNSUPPORTED_FRAME_INTERVAL);
    }
    
    ds_device_get_settings()->frame_period_usec = frame_interval_usec;
    
    /* Reset the current frame counter */
    bsp_reset_frame_counter();
    
    g_adc_sampling_completed = false;
    
    radar_set_data_acq_mode(RADAR_AUTO_TRIGGER);
    
    bsp_timer_start(&TIMER_FRAME_TRIG);
  }
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_test_antennas(Radar_Handle_t device, uint8_t tx_antenna_mask, uint8_t rx_antenna_mask, uint32_t frequency_kHz, uint8_t tx_power)
{
  /* Suppress compiler warnings of unused variables */
  XMC_UNUSED_ARG(device);
  XMC_UNUSED_ARG(tx_antenna_mask);
  XMC_UNUSED_ARG(rx_antenna_mask);
  XMC_UNUSED_ARG(frequency_kHz);
  XMC_UNUSED_ARG(tx_power);
  
  return (RADAR_ERR_FEATURE_UNAVAILABLE);
}

//============================================================================

void radar_get_fft_configuration(Radar_Handle_t device, FFT_configuration_t *p_fft_config)
{
	XMC_UNUSED_ARG(device);

	ds_ep_fetch_fft_config(p_fft_config);
}

/* --- End of File -------------------------------------------------------- */
