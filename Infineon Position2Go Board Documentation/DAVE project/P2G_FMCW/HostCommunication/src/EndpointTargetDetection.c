/**
 * \file EndpointTargetDetection.c
 *
 * \brief ...
 *
 * ...
 *
 */

/* ===========================================================================
 ** Copyright (C) 2016-2018 Infineon Technologies AG
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

#include "EndpointTargetDetection.h"
#include "EndpointRadarChangeCodes.h"
#include "PayloadHelper.h"

/*
==============================================================================
   2. LOCAL DEFINITIONS
==============================================================================
 */

/**
 * \internal
 * \defgroup EndpointTargetDetectionCodes
 *
 * \brief These are the command codes to identify the payload type.
 *
 * Each payload message of the supported endpoint starts with one of these
 * commend codes.
 *
 * @{
 */
#define  MSG_GET_DSP_SETTINGS		    0x00 	/**< A message to retrieve DSP settings */

#define  MSG_SET_DSP_SETTINGS           0x01 	/**< A message to set DSP settings */

#define  MSG_GET_TARGETS                0x02	/**< A message to get targets information */

#define  MSG_GET_RANGE_THRESHOLD		0x03    /**< A message to get range detection threshold */

/** @} */

Target_Info_t  target_list[MAX_NUM_OF_TARGETS];		/**< Defined global structured array to avoid stack over-flow issues */

/*
==============================================================================
   5. LOCAL FUNCTION PROTOTYPES
==============================================================================
 */

static uint16_t send_dsp_settings (uint8_t endpoint, Radar_Handle_t radar_driver);

static uint16_t send_targets_info (uint8_t endpoint, Radar_Handle_t radar_driver);

static uint16_t send_range_threshold (uint8_t endpoint, Radar_Handle_t radar_driver);

/*
==============================================================================
  6. LOCAL FUNCTIONS
==============================================================================
 */

static uint16_t send_dsp_settings (uint8_t endpoint, Radar_Handle_t radar_driver)
{
	uint16_t error_code = 0;

	DSP_Settings_t  dsp_settings;

	error_code = radar_get_dsp_settings(radar_driver, &dsp_settings);

	if (error_code == RADAR_ERR_OK)
	{
		uint8_t message[27];

		/* compile message */
		wr_payload_u8 (message,  0, MSG_GET_DSP_SETTINGS);

		wr_payload_u8 (message,  1, dsp_settings.range_mvg_avg_length);
		wr_payload_u16(message,  2, dsp_settings.min_range_cm);
		wr_payload_u16(message,  4, dsp_settings.max_range_cm);
		wr_payload_u16(message,  6, dsp_settings.min_speed_kmh);
		wr_payload_u16(message,  8, dsp_settings.max_speed_kmh);
		wr_payload_u16(message, 10, dsp_settings.min_angle_degree);
		wr_payload_u16(message, 12, dsp_settings.max_angle_degree);
		wr_payload_u16(message, 14, dsp_settings.range_threshold);
		wr_payload_u16(message, 16, dsp_settings.speed_threshold);

		wr_payload_u16(message, 18, 0); // see MMWSW-603
		wr_payload_u8 (message, 20, dsp_settings.enable_tracking);
		wr_payload_u8 (message, 21, dsp_settings.num_of_tracks);

		wr_payload_u8 (message, 22, dsp_settings.median_filter_length);
		wr_payload_u8 (message, 23, dsp_settings.enable_mti_filter);

		wr_payload_u16(message, 24, dsp_settings.mti_filter_length);
		wr_payload_u8 (message, 26, 0); // see MMWSW-603

		/* send message */
		protocol_send_header(endpoint, sizeof(message));
		protocol_send_payload(message, sizeof(message));
		protocol_send_tail();
	}

	return error_code;
}

//==============================================================================

static uint16_t send_targets_info (uint8_t endpoint, Radar_Handle_t radar_driver)
{
	uint16_t i;

	uint16_t error_code = 0;

	uint8_t num_of_targets = 0;

	uint8_t message = MSG_GET_TARGETS;

	for (i = 0; i < MAX_NUM_OF_TARGETS; i++)		// initialize target list
	{
		target_list[i].target_id = 0;
		target_list[i].level = 0.f;
		target_list[i].radius = 0.f;
		target_list[i].azimuth = 0.f;
		target_list[i].elevation = 0.f;
		target_list[i].radial_speed = 0.f;
		target_list[i].azimuth_speed = 0.f;
		target_list[i].elevation_speed = 0.f;
	}

	error_code = radar_get_target_info(radar_driver, target_list, &num_of_targets);

	if (error_code == RADAR_ERR_OK)
	{
		/* send message */
		protocol_send_header(endpoint, num_of_targets * sizeof(Target_Info_t) + 1 + 1);

		protocol_send_payload (&message, 1);

		protocol_send_payload ((const uint8_t*) &num_of_targets, 1);

		for (i = 0; i < num_of_targets; i++)
		{
			protocol_send_payload ((uint8_t*) &target_list[i], sizeof(Target_Info_t));
		}

		protocol_send_tail();
	}

	return error_code;
}

//==============================================================================

static uint16_t send_range_threshold (uint8_t endpoint, Radar_Handle_t radar_driver)
{
	uint16_t error_code = 0;

	uint16_t threshold = 0;

	error_code = radar_get_range_detection_threshold(radar_driver, &threshold);

	if (error_code == RADAR_ERR_OK)
	{
		uint8_t message[3];

		/* compile message */
		wr_payload_u8 (message, 0, MSG_GET_RANGE_THRESHOLD);
		wr_payload_u16(message, 1, threshold);

		/* send message */
		protocol_send_header(endpoint, sizeof(message));
		protocol_send_payload(message, sizeof(message));
		protocol_send_tail();
	}

	return error_code;
}

/*
==============================================================================
   7. EXPORTED FUNCTIONS
==============================================================================
 */

uint16_t ep_target_detection_handle_message (uint8_t endpoint, uint8_t* message_data, uint16_t num_bytes, void* context)
{
	Radar_Handle_t radar_driver = (Radar_Handle_t) context;

	switch (message_data[0])
	{
	case MSG_GET_DSP_SETTINGS:

		if (num_bytes == 1)
		{
			return send_dsp_settings(endpoint, radar_driver);
		}

		break;

	case MSG_SET_DSP_SETTINGS:

		if (num_bytes >= 18)
		{
			DSP_Settings_t dsp_settings;

			/* extract parameters from message */
			dsp_settings.range_mvg_avg_length = rd_payload_u8(message_data, 1);
			dsp_settings.min_range_cm     = rd_payload_u16(message_data, 2);
			dsp_settings.max_range_cm     = rd_payload_u16(message_data, 4);
			dsp_settings.min_speed_kmh    = rd_payload_u16(message_data, 6);
			dsp_settings.max_speed_kmh    = rd_payload_u16(message_data, 8);
			dsp_settings.min_angle_degree = rd_payload_u16(message_data, 10);
			dsp_settings.max_angle_degree = rd_payload_u16(message_data, 12);
			dsp_settings.range_threshold  = rd_payload_u16(message_data, 14);
			dsp_settings.speed_threshold  = rd_payload_u16(message_data, 16);

			// skip uint16 entry for index 18
			dsp_settings.enable_tracking  = rd_payload_u8 (message_data, 20);
			dsp_settings.num_of_tracks    = rd_payload_u8 (message_data, 21);

			dsp_settings.median_filter_length = rd_payload_u8 (message_data, 22);
			dsp_settings.enable_mti_filter    = rd_payload_u8 (message_data, 23);

			dsp_settings.mti_filter_length = rd_payload_u16(message_data, 24);
			// skip uint8 entry for index 26

			/* apply DSP Settings */
			return radar_set_dsp_settings(radar_driver, &dsp_settings);
		}

		break;

	case MSG_GET_TARGETS:

		if (num_bytes == 1)
		{
			return send_targets_info (endpoint, radar_driver);
		}

		break;

	case MSG_GET_RANGE_THRESHOLD:

		if (num_bytes == 1)
		{
			return send_range_threshold (endpoint, radar_driver);
		}

		break;

	default:
		break;
	}

	return PROTOCOL_STATUS_INVALID_PAYLOAD;

}

//==============================================================================

void ep_target_detection_handle_change (uint8_t endpoint, void* context, uint32_t what)
{
	Radar_Handle_t radar_driver = (Radar_Handle_t)context;

	switch (what)
	{
	case EP_RADAR_CHNG_FRAME_FORMAT:
		send_targets_info(endpoint, radar_driver);
		break;

	default:
		break;
	}
}

/* --- End of File -------------------------------------------------------- */
