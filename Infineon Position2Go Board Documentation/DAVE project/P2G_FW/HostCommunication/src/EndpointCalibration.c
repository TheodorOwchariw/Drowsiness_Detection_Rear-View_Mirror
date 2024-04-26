/**
 * \file EndpointCalibration.c
 *
 * \brief ...
 *
 * ...
 *
 */

/* ===========================================================================
** Copyright (C) 2018 Infineon Technologies AG
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

#include "EndpointCalibration.h"
#include "EndpointRadarChangeCodes.h"
#include "PayloadHelper.h"

/*
==============================================================================
   2. LOCAL DEFINITIONS
==============================================================================
*/

/**
 * \internal
 * \defgroup EndpointCalibrationCodes
 *
 * \brief These are the command codes to identify the payload type.
 *
 * Each payload message of the supported endpoint starts with one of these
 * commend codes.
 *
 * @{
 */

#define  MSG_SET_ADC_FLASH_CALIBRATION		0x01	/**< A message to save ADC calibration data from EEPROM */

#define  MSG_GET_ADC_FLASH_CALIBRATION     	0x02	/**< A message to read ADC calibration data from EEPROM */

#define  MSG_CLEAR_ADC_FLASH_CALIBRATION    0x03	/**< A message to delete ADC calibration data from EEPROM */

#define  MSG_SET_ADC_SRAM_CALIBRATION      	0x04	/**< A message to save ADC calibration data in SRAM */

#define  MSG_GET_ADC_SRAM_CALIBRATION      	0x05	/**< A message to read ADC calibration data in SRAM */

#define  MSG_CLEAR_ADC_SRAM_CALIBRATION     0x06	/**< A message to delete ADC calibration data from SRAM */

#define  MSG_SET_ALGO_FLASH_CALIBRATION		0x07  	/**< A message to store the Algo calibration in the Flash memory. */

#define  MSG_GET_ALGO_FLASH_CALIBRATION		0x08  	/**< A message to read current Algo calibration data stored in the Flash memory. */

#define  MSG_CLEAR_ALGO_FLASH_CALIBRATION	0x09  	/**< A message to clear the Algo calibration data from Flash memory. */

#define  MSG_SET_ALGO_SRAM_CALIBRATION		0x0A  	/**< A message to store the Algo calibration data in the SRAM memory. */

#define  MSG_GET_ALGO_SRAM_CALIBRATION		0x0B  	/**< A message to read current Algo calibration data stored in the SRAM memory. */

#define  MSG_CLEAR_ALGO_SRAM_CALIBRATION	0x0C  	/**< A message to clear the Algo calibration data from SRAM memory. */

/** @} */

/*
==============================================================================
   3. LOCAL FUNCTION PROTOTYPES
==============================================================================
*/

static uint16_t send_adc_calibration_data(uint8_t endpoint, Radar_Handle_t radar_driver, Calibration_Target_t target);

static uint16_t send_algo_calibration_data(uint8_t endpoint, Radar_Handle_t radar_driver, Calibration_Target_t target);

/*
==============================================================================
  4. LOCAL FUNCTIONS
==============================================================================
*/

static uint16_t send_adc_calibration_data(uint8_t endpoint, Radar_Handle_t radar_driver, Calibration_Target_t target)
{
	uint16_t error_code = 0;

	uint16_t num_of_bytes;

	void* calibration_data_ptr = NULL;

	//----------------------------------------------------------------------

	error_code = radar_read_calibration(radar_driver, target, CALIBRATION_DATA_ADC, &calibration_data_ptr, &num_of_bytes);

	if (error_code == RADAR_ERR_OK)
	{
		uint8_t message_header[3];

		/* setup message header */
		if (target == CALIBRATION_TARGET_FLASH)
		{
			wr_payload_u8(message_header, 0, MSG_GET_ADC_FLASH_CALIBRATION);
		}
		else
		{
			wr_payload_u8(message_header, 0, MSG_GET_ADC_SRAM_CALIBRATION);
		}

		wr_payload_u16(message_header, 1, num_of_bytes);

		/* send message header */
		protocol_send_header(endpoint, sizeof(message_header) + num_of_bytes);

		protocol_send_payload(message_header, sizeof(message_header));

		protocol_send_payload((uint8_t*)calibration_data_ptr, num_of_bytes);

		protocol_send_tail();
	}

	return error_code;
}

//==============================================================================

static uint16_t send_algo_calibration_data(uint8_t endpoint, Radar_Handle_t radar_driver, Calibration_Target_t target)
{
	uint16_t error_code = 0;

	uint16_t num_of_bytes;

	void* calibration_data_ptr;

	//----------------------------------------------------------------------

	error_code = radar_read_calibration(radar_driver, target, CALIBRATION_DATA_ALGO, &calibration_data_ptr, &num_of_bytes);

	if (error_code == RADAR_ERR_OK)
	{
		uint8_t message_header[5];

		/* setup message header */

		if (target == CALIBRATION_TARGET_FLASH)
		{
			wr_payload_u8(message_header, 0, MSG_GET_ALGO_FLASH_CALIBRATION);
		}
		else
		{
			wr_payload_u8(message_header, 0, MSG_GET_ALGO_SRAM_CALIBRATION);
		}

		wr_payload_u16(message_header, 1, ((Algo_Calibrations_t*)calibration_data_ptr)->distance_offset_cm);
		wr_payload_i16(message_header, 3, ((Algo_Calibrations_t*)calibration_data_ptr)->angle_offset_deg);

		/* send message header */
		protocol_send_header(endpoint, sizeof(message_header));

		protocol_send_payload(message_header, sizeof(message_header));

		protocol_send_tail();
	}

	return error_code;
}

/*
==============================================================================
   5. EXPORTED FUNCTIONS
==============================================================================
 */

uint16_t ep_calibration_handle_message(uint8_t endpoint, uint8_t* message_data, uint16_t num_bytes, void* context)
{
	Radar_Handle_t radar_driver = (Radar_Handle_t) context;

	uint16_t error_code;

	switch (message_data[0])
	{

	case MSG_GET_ADC_FLASH_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = send_adc_calibration_data(endpoint, radar_driver, CALIBRATION_TARGET_FLASH);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_GET_ADC_SRAM_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = send_adc_calibration_data(endpoint, radar_driver, CALIBRATION_TARGET_SRAM);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_SET_ADC_FLASH_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = radar_save_calibration(radar_driver, CALIBRATION_TARGET_FLASH, CALIBRATION_DATA_ADC, NULL);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_SET_ADC_SRAM_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = radar_save_calibration(radar_driver, CALIBRATION_TARGET_SRAM,CALIBRATION_DATA_ADC, NULL);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_CLEAR_ADC_FLASH_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = radar_clear_calibration(radar_driver, CALIBRATION_TARGET_FLASH, CALIBRATION_DATA_ADC);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_CLEAR_ADC_SRAM_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = radar_clear_calibration(radar_driver, CALIBRATION_TARGET_SRAM, CALIBRATION_DATA_ADC);

			return error_code;
		}

		break;

		//===============================================================

	case MSG_GET_ALGO_FLASH_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = send_algo_calibration_data(endpoint, radar_driver, CALIBRATION_TARGET_FLASH);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_GET_ALGO_SRAM_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = send_algo_calibration_data(endpoint, radar_driver, CALIBRATION_TARGET_SRAM);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_SET_ALGO_FLASH_CALIBRATION:

		if (num_bytes == 5)
		{
			Algo_Calibrations_t algo_calibration_settings;

			/* extract parameters from message */

			algo_calibration_settings.distance_offset_cm = rd_payload_u16(message_data, 1);
			algo_calibration_settings.angle_offset_deg   = rd_payload_i16(message_data, 3);

			error_code = radar_save_calibration(radar_driver, CALIBRATION_TARGET_FLASH, CALIBRATION_DATA_ALGO, &algo_calibration_settings);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_SET_ALGO_SRAM_CALIBRATION:

		if (num_bytes == 5)
		{
			Algo_Calibrations_t algo_calibration_settings;

			/* extract parameters from message */
			algo_calibration_settings.distance_offset_cm = rd_payload_u16(message_data, 1);
			algo_calibration_settings.angle_offset_deg   = rd_payload_i16(message_data, 3);

			error_code = radar_save_calibration(radar_driver, CALIBRATION_TARGET_SRAM, CALIBRATION_DATA_ALGO, &algo_calibration_settings);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_CLEAR_ALGO_FLASH_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = radar_clear_calibration(radar_driver, CALIBRATION_TARGET_FLASH, CALIBRATION_DATA_ALGO);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	case MSG_CLEAR_ALGO_SRAM_CALIBRATION:

		if (num_bytes == 1)
		{
			error_code = radar_clear_calibration(radar_driver, CALIBRATION_TARGET_SRAM, CALIBRATION_DATA_ALGO);

			return error_code;
		}

		break;

		//---------------------------------------------------------------

	default:

		break;
	}

	return PROTOCOL_STATUS_INVALID_PAYLOAD;
}

//==============================================================================

void ep_calibration_handle_change (uint8_t endpoint, void* context, uint32_t what)
{
	Radar_Handle_t radar_driver = (Radar_Handle_t)context;

	(void) radar_driver; /* Suppress compiler warnings of unused variable */

	switch (what)
	{
	case MSG_SET_ALGO_FLASH_CALIBRATION:
	case MSG_SET_ALGO_SRAM_CALIBRATION:
		break;
	default:
		break;
	}
}

/* --- End of File -------------------------------------------------------- */
