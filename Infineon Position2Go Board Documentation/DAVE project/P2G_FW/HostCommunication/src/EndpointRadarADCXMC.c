/**
 * \file EndpointRadarADCXMC.c
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
#include "EndpointRadarADCXMC.h"
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
#define MSG_GET_CONFIGURATION     0x00 /**< A message to query the ADC
                                            configuration. */
#define MSG_SET_CONFIGURATION     0x01 /**< A message containing the ADC
                                            configuration. */
/** @} */

/*
==============================================================================
   3. LOCAL FUNCTION PROTOTYPES
==============================================================================
*/

uint16_t send_adc_configuration(uint8_t endpoint,
                                       Radar_Handle_t radar_driver);


uint8_t adc_endpoint = 0;

/*
==============================================================================
  4. LOCAL FUNCTIONS
==============================================================================
*/

uint16_t send_adc_configuration(uint8_t endpoint,
                                       Radar_Handle_t radar_driver)
{
    uint16_t error_code;
    Adc_Configuration_t adc_configuration;

    adc_endpoint = endpoint;

    /* read ADC configuration from device */
    error_code = radar_get_adc_configuration(radar_driver,
                                             &adc_configuration);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[7];

        /* compile message */
        wr_payload_u8 (message, 0, MSG_SET_CONFIGURATION);
        wr_payload_u32(message, 1, adc_configuration.samplerate_Hz);
        wr_payload_u8 (message, 5, adc_configuration.resolution);
        wr_payload_u8 (message, 6, adc_configuration.use_post_calibration);

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

uint16_t ep_radar_adcxmc_handle_message(uint8_t endpoint,
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
            return send_adc_configuration(endpoint, radar_driver);
        }
        break;

    case MSG_SET_CONFIGURATION:
        if (num_bytes == 7)
        {
            unsigned error_code;
            Adc_Configuration_t config;

            /* extract parameters from message */
            config.samplerate_Hz        = rd_payload_u32(message_data, 1);
            config.resolution           = rd_payload_u8 (message_data, 5);
            config.use_post_calibration = rd_payload_u8 (message_data, 6);

            /* apply ADC configuration */
            error_code = radar_set_adc_configuration(radar_driver, &config);

            /* inform other endpoints of the change */
            protocol_broadcast_change(radar_driver,
                                      EP_RADAR_CHNG_ADC_SETTINGS);

            return error_code;
        }
        break;

    default:
        break;
    }
    return PROTOCOL_STATUS_INVALID_PAYLOAD;
}

void ep_radar_adcxmc_handle_change(uint8_t endpoint, void* context,
                                   uint32_t what)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;

    if (what == EP_RADAR_CHNG_ADC_SETTINGS)
        send_adc_configuration(endpoint, radar_driver);
}

/* --- End of File -------------------------------------------------------- */
