/**
 * \file EndpointRadarFMCW.c
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
** is not permitted without Infineon's prior written authorisation.
** ===========================================================================
*/

/*
==============================================================================
   1. INCLUDE FILES
==============================================================================
*/
#include "EndpointRadarFMCW.h"
#include "EndpointRadarChangeCodes.h"
#include "PayloadHelper.h"

/*
==============================================================================
   2. LOCAL DEFINITIONS
==============================================================================
*/

/**
 * \defgroup MessageTypes
 *
 * \brief This end point knows these message types.
 *
 * The first payload byte of every message is one of these codes.
 *
 * @{
 */
#define MSG_GET_CONFIGURATION     0x00  /**< A message to query the FMCW
                                             configuration. */
#define MSG_SET_CONFIGURATION     0x01  /**< A message containing the FMCW
                                             configuration. */
#define MSG_GET_BW_PER_SECOND     0x02  /**< A message to query the bandwidth
                                             per second. */
#define MSG_SET_BW_PER_SECOND     0x02  /**< A message containing the
                                             bandwidth per second. */

#define MSG_GET_ADV_CONFIGURATION 0x03  /**< A message to query the FMCW
                                             configuration. */
#define MSG_SET_ADV_CONFIGURATION 0x04  /**< A message containing the FMCW
                                             configuration. */


/** @} */


static uint32_t g_host_supports_advanced_configuration = 0;


/*
==============================================================================
   3. LOCAL FUNCTION PROTOTYPES
==============================================================================
*/

static uint16_t send_fmcw_configuration(uint8_t endpoint,
                                        Radar_Handle_t radar_driver);

static uint16_t send_bandwidth_per_second(uint8_t endpoint,
                                          Radar_Handle_t radar_driver);

/*
==============================================================================
   4. LOCAL FUNCTIONS
==============================================================================
*/

static uint16_t send_fmcw_configuration(uint8_t endpoint,
                                        Radar_Handle_t radar_driver)
{
    uint16_t error_code;
    Fmcw_Configuration_t fmcw_configuration;

    /* read FMCW configuration from device */
    error_code = radar_get_fmcw_configuration(radar_driver,
                                              &fmcw_configuration);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[11 + 12];
        size_t len = 11;

        /* compile message */
        wr_payload_u8 (message,   0, MSG_SET_CONFIGURATION);
        wr_payload_u32(message,   1, fmcw_configuration.lower_frequency_kHz);
        wr_payload_u32(message,   5, fmcw_configuration.upper_frequency_kHz);
        wr_payload_u8 (message,   9, fmcw_configuration.direction);
        wr_payload_u8 (message,  10, fmcw_configuration.tx_power);

        if( g_host_supports_advanced_configuration ) {
			wr_payload_u32 (message, 11, fmcw_configuration.up_chirp_duration_nsec);
			wr_payload_u32 (message, 15, fmcw_configuration.down_chirp_duration_nsec);
			wr_payload_u32 (message, 19, fmcw_configuration.chirp_to_chirp_delay_nsec);
			len += 4*sizeof(uint32_t);
        }
        /* send message */
        protocol_send_header(endpoint, len);
        protocol_send_payload(message, len);
        protocol_send_tail();
    }
    return error_code;
}

static uint16_t send_bandwidth_per_second(uint8_t endpoint,
                                          Radar_Handle_t radar_driver)
{
    uint16_t error_code;
    uint32_t bandwidth_per_second;

    /* read bandwidth per second from device */
    error_code = radar_get_bandwidth_per_second(radar_driver,
                                                &bandwidth_per_second);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[5];

        /* compile message */
        wr_payload_u8 (message, 0, MSG_SET_BW_PER_SECOND);
        wr_payload_u32(message, 1, bandwidth_per_second);

        /* send message */
        protocol_send_header(endpoint, sizeof(message));
        protocol_send_payload(message, sizeof(message));
        protocol_send_tail();
    }
    return error_code;
}

/*
==============================================================================
   5. EXPORTED FUNCTIONS
==============================================================================
*/

uint16_t ep_radar_fmcw_handle_message(uint8_t endpoint, uint8_t* message_data,
                                      uint16_t num_bytes, void* context)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;

    switch (message_data[0])
    {
    case MSG_GET_CONFIGURATION:
        if (num_bytes == 1)
        {
        	// a local variable does not work, because the updates are also sent on updates!
        	g_host_supports_advanced_configuration = 0;
            return send_fmcw_configuration(endpoint, radar_driver);
        }
        break;
    case MSG_GET_ADV_CONFIGURATION:
		if (num_bytes == 1)
		{
			g_host_supports_advanced_configuration = 1;
			return send_fmcw_configuration(endpoint, radar_driver);
		}
		break;

    case MSG_SET_CONFIGURATION:
        if (num_bytes == 11)
        {
            unsigned error_code;
            Fmcw_Configuration_t config;

            /* extract parameters from message */
            config.lower_frequency_kHz = rd_payload_u32(message_data,  1);
            config.upper_frequency_kHz = rd_payload_u32(message_data,  5);
            config.direction           = rd_payload_u8 (message_data,  9);
            config.tx_power            = rd_payload_u8 (message_data, 10);

            /* apply FMCW configuration */
            error_code = radar_set_fmcw_configuration(radar_driver, &config);

            /* inform other endpoints of the change */
            protocol_broadcast_change(radar_driver,
                                      EP_RADAR_CHNG_FMCW_SETTINGS);

            return error_code;
        }
        break;

    case MSG_GET_BW_PER_SECOND:
        if (num_bytes == 1)
        {
            return send_bandwidth_per_second(endpoint, radar_driver);
        }
        break;

    default:
        break;
    }
    return PROTOCOL_STATUS_INVALID_PAYLOAD;
}

void ep_radar_fmcw_handle_change(uint8_t endpoint, void* context,
                                 uint32_t what)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;

    switch (what)
    {
    case EP_RADAR_CHNG_FRAME_FORMAT:
        send_bandwidth_per_second(endpoint, radar_driver);
        break;

    case EP_RADAR_CHNG_FMCW_SETTINGS:
        send_fmcw_configuration(endpoint, radar_driver);
        send_bandwidth_per_second(endpoint, radar_driver);
        break;

    case EP_RADAR_CHNG_ADC_SETTINGS:
        send_bandwidth_per_second(endpoint, radar_driver);
        break;

    default:
        break;
    }
}

/* --- End of File -------------------------------------------------------- */
