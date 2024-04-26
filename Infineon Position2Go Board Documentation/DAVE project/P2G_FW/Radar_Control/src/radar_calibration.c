/**
    @file: radar_calibration.c

  	@brief: This file contains the driver API Calibration features of
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
#include "radar_calibration.h"
#include "radar_control.h"
#include "position2go_bsp.h"
#include "ds_radar_interface.h"

/*
==============================================================================
   2. DATA
==============================================================================
 */

uint16_t   adc_i1_calibration[BSP_NUM_OF_ADC_CHANNELS * MAX_CALIB_SAMPLES_PER_CHANNEL];
uint16_t*  adc_q1_calibration = &adc_i1_calibration[MAX_CALIB_SAMPLES_PER_CHANNEL * 1];
uint16_t*  adc_i2_calibration = &adc_i1_calibration[MAX_CALIB_SAMPLES_PER_CHANNEL * 2];
uint16_t*  adc_q2_calibration = &adc_i1_calibration[MAX_CALIB_SAMPLES_PER_CHANNEL * 3];

Algo_Calibrations_t  algo_calibration_s = {0,0};

extern Radar_Handle_t h_radar_device;

/*
==============================================================================
   3. LOCAL FUNCTION PROTOTYPES
==============================================================================
 */

/**
 * \brief  This function saves the calibration data to the SRAM or Flash based on the target memory defined by argument \ref Calibration_Target_t.
 *
 * \param[in]	target	Target of the calibration defined by \ref Calibration_Target_t.
 */
static void save_adc_calibration(Calibration_Target_t target);

/**
 * \brief  This function reads the calibration data from the SRAM or Flash based on the target memory defined by argument \ref Calibration_Target_t.
 *
 * \param[in]	target	Target of the calibration defined by \ref Calibration_Target_t.
 *
 * \return		Number of bytes occupied by calibration samples for all ADC channels
 */
static uint32_t read_adc_calibration(Calibration_Target_t target);

/**
 * \brief  This function erase calibration data from the SRAM or Flash based on the target memory defined by argument \ref Calibration_Target_t.
 *
 * \param[in]	target	Target of the calibration defined by \ref Calibration_Target_t.
 */
static void clear_adc_calibration(Calibration_Target_t target);

/**
 * \brief  This function saves the algorithm related calibrations in the  SRAM/Flash memory based on target argument type.
 *
 * \param[in]	target					Target of the calibration defined by \ref Calibration_Target_t.
 * \param[in]	*algo_calibration_ptr	A const pointer to the structure defined by \ref Algo_Calibrations_t.
 */
static void save_algo_calibration(Calibration_Target_t target, const Algo_Calibrations_t* algo_calibration_ptr);

/**
 * \brief  This function reads the algorithm related calibrations from SRAM/Flash memory based on target argument type.
 *
 * \param[in]	target	Target memory of the calibration defined by \ref Calibration_Target_t.
 *
 * \return		Number of bytes used in algo calibration structure
 */
static uint32_t read_algo_calibration(Calibration_Target_t target);

/**
 * \brief  This function clears the algorithm related calibrations from SRAM/Flash memory based on target argument type.
 *
 * \param[in]	target	Target memory of the calibration defined by \ref Calibration_Target_t.
 */
static void clear_algo_calibration(Calibration_Target_t target);

/**
 * \brief  This function writes the calibration data to the Flash.
 */
static void write_calibration_flash_data(void);

/**
 * \brief  This is a helper function to sets the angle offset compensation to get zero angle for target in front of the radar.
 *
 * \param[in]	offset_deg	Signed integer 16-bit value between [-30, +30] in unit of degrees.
 */
static void set_angle_offset(int16_t offset_deg);

/**
 * \brief  This is a helper function to sets the range offset compensation.
 *
 * \param[in]	offset_cm	Unsigned int16 signed value between [0, 100] in units of cm.
 */
static void set_range_offset(uint16_t offset_cm);

/*
==============================================================================
   4. EXPORTED FUNCTIONS
==============================================================================
 */

void radar_init_calibration(Radar_Handle_t device)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */

  memset(adc_i1_calibration, DEFAULT_CALIBRATION_VALUE, ADC_CALIB_LEN_BYTES);
}

//============================================================================

uint16_t radar_save_calibration(Radar_Handle_t device, Calibration_Target_t target, Calibration_Data_t calib_data_type, const void* calibration_data_ptr)
{
  if (calib_data_type == CALIBRATION_DATA_ADC)
  {
    uint16_t status  = 1;
    Frame_Info_t frame_info;
    
    g_do_processing = false;
    
    while(status)
    {
      status = radar_get_frame(device, &frame_info, true);
    }
    
    save_adc_calibration(target);
    g_do_processing = true;
  }
  else
  {
    
    save_algo_calibration(target, (const Algo_Calibrations_t*)calibration_data_ptr);
  }
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_read_calibration(Radar_Handle_t device, Calibration_Target_t target, Calibration_Data_t calib_data_type, void** data_ptr, uint16_t* num_of_bytes)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  if (calib_data_type == CALIBRATION_DATA_ADC)
  {
    *num_of_bytes = read_adc_calibration(target);
    *data_ptr = adc_i1_calibration;
  }
  else
  {
    *num_of_bytes = read_algo_calibration(target);
    *data_ptr = &algo_calibration_s;
  }
  
  return (RADAR_ERR_OK);
}

//============================================================================

uint16_t radar_clear_calibration(Radar_Handle_t device, Calibration_Target_t target, Calibration_Data_t calib_data_type)
{
  XMC_UNUSED_ARG(device); /* Suppress compiler warnings of unused variable */
  
  if (calib_data_type == CALIBRATION_DATA_ADC)
  {
    clear_adc_calibration(target);
  }
  else
  {
    clear_algo_calibration(target);
  }
  
  return (RADAR_ERR_OK);
}

/*
==============================================================================
   5. LOCAL FUNCTIONS
==============================================================================
 */

static void save_adc_calibration(Calibration_Target_t target)
{
  uint32_t cal_idx = 0;
  uint32_t interp_factor = MAX_CALIB_SAMPLES_PER_CHANNEL / ds_device_get_settings()->num_samples_per_chirp;	// interpolation factor

  acq_buf_obj *p_acq_buf = ds_get_active_acq_buf();
  
  uint16_t *adc_temp[4];
  adc_temp[0] = (uint16_t *)get_buffer_address_by_chirp(p_acq_buf, 0, 0);
  adc_temp[1] = (uint16_t *)get_buffer_address_by_chirp(p_acq_buf, 1, 0);
  adc_temp[2] = (uint16_t *)get_buffer_address_by_chirp(p_acq_buf, 2, 0);
  adc_temp[3] = (uint16_t *)get_buffer_address_by_chirp(p_acq_buf, 3, 0);

  /* Fill the calibration buffer with raw ADC data */
  for (uint32_t i = 0; i < ds_device_get_settings()->num_samples_per_chirp; i++)
  {
    adc_i1_calibration[cal_idx] = (adc_temp[0])[i];
    adc_q1_calibration[cal_idx] = (adc_temp[1])[i];
    adc_i2_calibration[cal_idx] = (adc_temp[2])[i];
    adc_q2_calibration[cal_idx] = (adc_temp[3])[i];
    
    cal_idx += interp_factor;
  }
  
  /* Permanently write the calibration data in to the Flash */
  if (target == CALIBRATION_TARGET_FLASH)
  {
    /* Read the flashed algorithm parameters */
    bsp_flash_read_data(ADC_CALIB_LEN_BYTES, (uint8_t*)&algo_calibration_s, sizeof(algo_calibration_s));
    
    write_calibration_flash_data();
  }
}

//============================================================================

static void clear_adc_calibration(Calibration_Target_t target)
{
  Radar_Handle_t device = NULL;
  
  /* Reset SRAM buffer to a default value */
  radar_init_calibration(device);
  
  if (target == CALIBRATION_TARGET_FLASH)
  {
    bsp_flash_read_data(ADC_CALIB_LEN_BYTES, (uint8_t*)&algo_calibration_s, sizeof(algo_calibration_s));
    
    write_calibration_flash_data();
  }
}

//============================================================================

static uint32_t read_adc_calibration(Calibration_Target_t target)
{
	Radar_Handle_t device = NULL;

	if (target == CALIBRATION_TARGET_FLASH)
	{
		if (bsp_flash_is_empty())
		{
			/* Flash memory is empty, set the initial data to default value */
			radar_init_calibration(device);
		}
		else
		{
			/* Read calibration stored in flash */
			bsp_flash_read_data(0U, (uint8_t*)adc_i1_calibration, ADC_CALIB_LEN_BYTES);
		}
	}

	return (ADC_CALIB_LEN_BYTES);
}

//============================================================================

static void save_algo_calibration(Calibration_Target_t target, const Algo_Calibrations_t*  algo_calibration_ptr)
{
	algo_calibration_s = *algo_calibration_ptr;

	if (target == CALIBRATION_TARGET_FLASH)
	{
		read_adc_calibration(CALIBRATION_TARGET_FLASH);

		write_calibration_flash_data();
	}

	set_angle_offset(algo_calibration_s.angle_offset_deg);
	set_range_offset(algo_calibration_s.distance_offset_cm);
}

//============================================================================

static uint32_t read_algo_calibration(Calibration_Target_t target)
{
	size_t algo_size = sizeof(algo_calibration_s);

	if (target == CALIBRATION_TARGET_FLASH)
	{
		if (bsp_flash_is_empty())
		{
			memset(&algo_calibration_s, DEFAULT_CALIBRATION_VALUE, algo_size);
		}
		else
		{
			bsp_flash_read_data(ADC_CALIB_LEN_BYTES, (uint8_t*)&algo_calibration_s, algo_size);
		}
	}

	set_angle_offset(algo_calibration_s.angle_offset_deg);
	set_range_offset(algo_calibration_s.distance_offset_cm);

	return (algo_size);
}

//============================================================================

static void clear_algo_calibration(Calibration_Target_t target)
{
	memset(&algo_calibration_s, 0, sizeof(algo_calibration_s));

	if (target == CALIBRATION_TARGET_FLASH)
	{
		read_adc_calibration(CALIBRATION_TARGET_FLASH);

		write_calibration_flash_data();
	}

	set_angle_offset(algo_calibration_s.angle_offset_deg);
	set_range_offset(algo_calibration_s.distance_offset_cm);
}

//============================================================================

static void write_calibration_flash_data(void)
{
	bsp_flash_write_data(0U, (uint8_t*)adc_i1_calibration, ADC_CALIB_LEN_BYTES);

	bsp_flash_write_data(ADC_CALIB_LEN_BYTES, (uint8_t*)&algo_calibration_s, sizeof(algo_calibration_s));

	bsp_flash_flush_data();
}

//============================================================================

static void set_angle_offset(int16_t offset)
{
	if ((offset > BSP_MIN_ANGLE_OFFSET_DEG) || (offset < BSP_MAX_ANGLE_OFFSET_DEG))
	{
		ds_algo_get_settings()->angle_offset_deg = offset;
	}
}

//============================================================================

static void set_range_offset(uint16_t offset)
{
	if (offset < BSP_MAX_RANGE_OFFSET_CM)
	{
		ds_algo_get_settings()->range_offset_cm = offset;
	}
}

/* --- End of File -------------------------------------------------------- */
