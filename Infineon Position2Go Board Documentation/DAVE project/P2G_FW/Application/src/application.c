/**
    @file: application.c

    @brief: This file includes the functions for the Position2Go
            Firmware Application Layer.
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

#include "application.h"
#include "datastore.h"

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
 */

/* Define the number of supported endpoints */
#define NUMBER_SUPPORTED_ENDPOINTS         (6U)

/*
==============================================================================
   3. DATA
==============================================================================
 */

/* Firmware version declaration */
Firmware_Information_t firmware_information =
{
	"Position2Go BGT24MTR12 Demo Radar Board",
	FW_VERSION_MAJOR,
	FW_VERSION_MINOR,
	FW_VERSION_REVISION
};

/* Radar handler declaration */
Radar_Handle_t h_radar_device;

/* A global buffer to hold all endpoints used by Host Communication library */
Endpoint_t endpoint_list[NUMBER_SUPPORTED_ENDPOINTS];

/* Global variable to store the state if USB initialization, used with Radar GUI */
volatile uint8_t g_usb_init = false;

/* Radar application demo global state */
Radar_App_State_t radar_app_state = APP_IDLE;

/*
==============================================================================
   4. USER CALLBACK FUNCTION PROTOTYPES
==============================================================================
 */

/* This callback function will be executed by radar control when completing data acquisition.
 *
 */
void acq_completed_cb(void);

/* This callback function...
 *
 */
void algo_completed_cb(void);

/* This callback function will be executed by radar control when starting data acquisition process,
 * it is used to get information about the beginning of data acquisition.
 */
void app_data_acq_started_cb(void *parameter);

/* This callback function will be executed by radar control at the end of data acquisition,
 * it is used to inform upper layers that data sampling is done and data is ready in user buffer.
 */
void app_data_acq_done_cb(void *parameter);

/* This function will be executed at the end of data acquisition.
 * This function pointer should be set to the assigned algorithm before starting data processing.
 */
static algorithm algo_processing_function = NULL;

/*
==============================================================================
   5. EXPORTED FUNCTIONS
==============================================================================
 */

void app_init(void)
{
  /* Initialize USB */
  if(bsp_usb_init(BSP_USB_TIMEOUT) == 0)
  {
    g_usb_init = true;
  }

  /* Open the available radar device */
  if(radar_open_device(0, &h_radar_device) != RADAR_ERR_OK)
  {
    XMC_DEBUG("Opening Radar device failed\n");
    while(1U);
  }

  /* Register user callback functions */
  radar_register_callback(RADAR_ACQUISITION_STARTED_CB_ID, app_data_acq_started_cb);
  radar_register_callback(RADAR_ACQUISITION_DONE_CB_ID, app_data_acq_done_cb);

  /* Register all endpoints to the communication protocol */
  Endpoint_t ep_list[NUMBER_SUPPORTED_ENDPOINTS] =
  {
    EP_RADAR_BASE_ENDPOINT_LIST_ENTRY(h_radar_device),
    EP_RADAR_ADCXMC_ENDPOINT_LIST_ENTRY(h_radar_device),
#if (FW_MODULATION_TYPE == DOPPLER_MODULATION) /* Doppler modulation */
    EP_RADAR_DOPPLER_ENDPOINT_LIST_ENTRY(h_radar_device),
#else
    EP_RADAR_FMCW_ENDPOINT_LIST_ENTRY(h_radar_device),
#endif
    EP_RADAR_INDUSTRIAL_ENDPOINT_LIST_ENTRY(h_radar_device),
    EP_RADAR_P2G_ENDPOINT_LIST_ENTRY(h_radar_device),
    EP_CALIBRATION_ENDPOINT_LIST_ENTRY(h_radar_device)
  };

  memcpy(endpoint_list, ep_list, NUMBER_SUPPORTED_ENDPOINTS * sizeof(Endpoint_t));

  /* Initialize the internal protocol state and the communication interface */
  protocol_init(endpoint_list, sizeof(endpoint_list)/sizeof(endpoint_list[0]), SYSTIMER_GetTime, BSP_HOST_TIMEOUT_USEC, system_reset);

  /* Start data acquisition at the defined frame interval */
  if (radar_start_acquisition(h_radar_device) != RADAR_ERR_OK)
  {
    XMC_DEBUG("Starting Radar data acquisition failed\n");
    while(1U);
  }
}

//===========================================================================

void app_register_algo_process(algorithm algo_processor)
{
  /* Register algorithm processing function */
  algo_processing_function = algo_processor;
}

//===========================================================================

void app_process(void)
{
  if (g_usb_init == true)
  {
    /* Run communication protocol */
    protocol_run();
  }
  else
  {
    /* Initialize the USB low layer driver and start USB data communication */
    if(bsp_usb_init(BSP_USB_TIMEOUT) == 0)
    {
      g_usb_init = true;
    }
  }

  /* Main radar application process */
  radar_app_process();
}

//===========================================================================

void radar_app_process(void)
{
  switch(radar_app_state)
  {
  case APP_IDLE:
    break;

  case APP_RADAR_DATA_ACQUISITION_COMPLETED:
    #if defined(FMCW_SUPPORTED) || defined(DOPPLER_SUPPORTED)
      if (g_do_processing)
      {
        /* Start the processing of the algorithm assigned in the function pointer 'algo_processor' */
    	if (algo_processing_function != NULL)
    	{
    	  algo_processing_function();
    	}
      }
    #endif

    /* Check if a new hardware settings are required */
    radar_app_state = APP_CHECK_FOR_SETTINGS_UPDATE;
    break;

   case APP_CHECK_FOR_SETTINGS_UPDATE:
	/* Check if a new hardware settings are required */
    if(ds_is_device_settings_updated())
	{
	  radar_app_state = APP_APPLY_DEVICE_SETTINGS;
	}
	else
	{
	  radar_app_state = APP_IDLE;
	}
	break;
  case APP_APPLY_DEVICE_SETTINGS:
    if(ds_is_device_settings_updated())
	{
	  device_settings_t *p_dev_settings = ds_device_get_settings();

	  if(radar_apply_hw_settings(h_radar_device, UPDATE, p_dev_settings) == RADAR_ERR_OK)
	  {
	    ds_device_settings_copy_to_shadow();
	  }
	}
    radar_app_state = APP_IDLE;
    break;

  default:
    break;
  }
}

//===========================================================================

void app_destructor(void)
{
  /* Unregister user callback functions */
  radar_unregister_callback(RADAR_ACQUISITION_STARTED_CB_ID);
  radar_unregister_callback(RADAR_ACQUISITION_DONE_CB_ID);
}

/*
==============================================================================
   6. USER CALLBACK FUNCTIONS
==============================================================================
 */

/* Here an example of user callback implementation */

/* This callback function will be executed by radar control when starting data acquisition process
 *  it can be used to get information about the beginning of data acquisition
 */
void app_data_acq_started_cb(void *parameter)
{
  /* Turn On LED at the beginning of data acquisition */
  bsp_led_on();
}

//===========================================================================

/* This callback function will be executed by radar control at the end of data acquisition,
 * it is used to inform upper layers that data sampling is done and data is ready in user buffer
 */
void app_data_acq_done_cb(void *parameter)
{
  /* Turn Off LED at the end of data acquisition */
  bsp_led_off();

  /* Set the main application state at the end of data acquisition */
  g_adc_sampling_completed = true;
  radar_app_state = APP_RADAR_DATA_ACQUISITION_COMPLETED;
  acq_completed_cb();
}

//===========================================================================

/* This function yet to be filled in case of software reset */
void system_reset(void)
{

}

/* --- End of File -------------------------------------------------------- */
