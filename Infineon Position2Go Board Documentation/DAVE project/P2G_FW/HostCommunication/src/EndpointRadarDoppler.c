/**
 * \file EndpointRadarDoppler.c
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

#include "EndpointRadarDoppler.h"
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
#define MSG_GET_CONFIGURATION     0x00 /**< A message to query the Doppler
                                            radar configuration. */
#define MSG_SET_CONFIGURATION     0x01 /**< A message containing the Doppler
                                            radar configuration. */
/** @} */

/*
==============================================================================
   3. LOCAL FUNCTION PROTOTYPES
==============================================================================
*/

static uint16_t send_doppler_configuration(uint8_t endpoint,
                                           Radar_Handle_t radar_driver);

/*
==============================================================================
   4. LOCAL FUNCTIONS
==============================================================================
*/

static uint16_t send_doppler_configuration(uint8_t endpoint,
                                           Radar_Handle_t radar_driver)
{
    uint16_t error_code;
    Doppler_Configuration_t doppler_configuration;

    /* read FMCW configuration from device */
    error_code = radar_get_doppler_configuration(radar_driver,
                                                 &doppler_configuration);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[6];

        /* compile message */
        wr_payload_u8 (message, 0, MSG_SET_CONFIGURATION);
        wr_payload_u32(message, 1, doppler_configuration.frequency_kHz);
        wr_payload_u8 (message, 5, doppler_configuration.tx_power);

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

uint16_t ep_radar_doppler_handle_message(uint8_t endpoint,
                                         uint8_t* message_data,
                                         uint16_t num_bytes,
                                         void* context)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;

    switch (message_data[0])
    {
    case MSG_GET_CONFIGURATION:
        if (num_bytes == 1)
        {
            return send_doppler_configuration(endpoint, radar_driver);
        }
        break;

    case MSG_SET_CONFIGURATION:
        if (num_bytes == 6)
        {
            unsigned error_code;
            Doppler_Configuration_t configuration;

            /* extract parameters from message */
            configuration.frequency_kHz = rd_payload_u32(message_data, 1);
            configuration.tx_power      = rd_payload_u8 (message_data, 5);

            /* apply FMCW configuration */
            error_code = radar_set_doppler_configuration(radar_driver,
                                                         &configuration);

            /* inform other endpoints of the change */
            protocol_broadcast_change(radar_driver,
                                      EP_RADAR_CHNG_DOPPLER_SETTINGS);

            return error_code;
        }
        break;

    default:
        break;
    }
    return PROTOCOL_STATUS_INVALID_PAYLOAD;
}

void ep_radar_doppler_handle_change(uint8_t endpoint, void* context,
                                    uint32_t what)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;

    if (what == EP_RADAR_CHNG_DOPPLER_SETTINGS)
        send_doppler_configuration(endpoint, radar_driver);
}

/* --- End of File -------------------------------------------------------- */
