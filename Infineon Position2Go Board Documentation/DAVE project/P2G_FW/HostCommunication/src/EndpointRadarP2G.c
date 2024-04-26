/**
 * \file EndpointRadarP2G.c
 *
 * \brief Extension of radar endpoint for Position2Go demo board functionality
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
** is not permitted without Infineon's prior written authorisation.
** ===========================================================================
*/

/*
==============================================================================
   1. INCLUDE FILES
==============================================================================
*/
#include "EndpointRadarP2G.h"
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
 
#define MSG_GET_PGA_LEVEL 			0x00  	/**< A message to query current PGA gain level. */
											   
#define MSG_SET_PGA_LEVEL			0x01  	/**< A message containing PGA gain level. */											  
											   
/** @} */

/*
==============================================================================
   3. LOCAL FUNCTION PROTOTYPES
==============================================================================
*/

static uint16_t send_pga_level(Radar_Handle_t device, uint8_t endpoint);

/*
==============================================================================
   4. LOCAL FUNCTIONS
==============================================================================
*/

static uint16_t send_pga_level(Radar_Handle_t device, uint8_t endpoint)
{
    uint16_t error_code;
    
	uint16_t gain_level;

    /* get current PGA Gain level from driver */
	error_code = radar_get_gain_level(device, &gain_level);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[3];

        /* compile message */
        wr_payload_u8 (message, 0, MSG_SET_PGA_LEVEL);
        wr_payload_u16(message, 1, gain_level);

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

uint16_t ep_radar_p2g_handle_message(uint8_t endpoint,
                                     uint8_t* message_data,
                                       uint16_t num_bytes,
                                       void* context)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;
	
    (void) radar_driver; /* Suppress compiler warnings of unused variable */

	uint16_t error_code;

    switch (message_data[0])
    {
		
    case MSG_GET_PGA_LEVEL:
	
        if (num_bytes == 1)
        {
			error_code = send_pga_level(radar_driver, endpoint);
			
            return error_code;
        }
        break;
	
    case MSG_SET_PGA_LEVEL:
	
        if (num_bytes == 3)
        {    		
            uint16_t  gain_level;

            /* extract parameters from message */
            gain_level = rd_payload_u16(message_data, 1);
            
            /* apply PGA Gain Level */
            error_code = radar_set_gain_level(radar_driver, gain_level);

            return error_code;
        }
        break;
	
    default:
	
		break;
    }

    return PROTOCOL_STATUS_INVALID_PAYLOAD;
}

void ep_radar_p2g_handle_change(uint8_t endpoint, void* context,
                                uint32_t what)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;

    (void) radar_driver; /* Suppress compiler warnings of unused variable */

    switch (what)
    {
    default:
        break;
    }
}

/* --- End of File -------------------------------------------------------- */
