/**
 * \file EndpointRadarIndustrial.c
 *
 * \brief Radar endpoint containing functionality pertinent to Industrial 
 *        radar demo kits, e.g. Distance2Go, Position2Go etc.
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
#include "EndpointRadarIndustrial.h"
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
 
#define MSG_DISABLE_DUTY_CYCLE 			0x00  	/**< A message to disable duty cycle. */
											   
#define MSG_ENABLE_DUTY_CYCLE			0x01  	/**< A message to enable duty cycle. */	

#define MSG_IS_ENABLE_DUTY_CYCLE		0x02  	/**< A message to check the status of duty cycle. */

#define MSG_DISABLE_BGT_LNA 			0x03  	/**< A message to disable BGT Rx LNA gain. */
											   
#define MSG_ENABLE_BGT_LNA				0x04  	/**< A message to enable BGT Rx LNA gain. */

#define MSG_IS_ENABLE_BGT_LNA			0x05  	/**< A message to check the status of BGT Rx LNA gain. */
											   
/** @} */

/*
==============================================================================
   3. EXPORTED FUNCTIONS
==============================================================================
*/

uint16_t ep_radar_industrial_handle_message(uint8_t endpoint, uint8_t* message_data,
											uint16_t num_bytes,	void* context)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;

	uint16_t error_code = 0;

    switch (message_data[0])
    {
		
    case MSG_DISABLE_DUTY_CYCLE:
	
        if (num_bytes == 1)
        {
			error_code = radar_set_duty_cycle(radar_driver, 0);
			
            return error_code;
        }
        break;
		
	case MSG_ENABLE_DUTY_CYCLE:
	
        if (num_bytes == 1)
        {
			error_code = radar_set_duty_cycle(radar_driver, 1);
			
            return error_code;
        }
        break;

	case MSG_IS_ENABLE_DUTY_CYCLE:

		if (num_bytes == 1)
		{
			uint8_t message[2];

			uint8_t flag;

			radar_get_duty_cycle(radar_driver, &flag);

			/* compile message */
			wr_payload_u8(message, 0, MSG_IS_ENABLE_DUTY_CYCLE);
			wr_payload_u8(message, 1, flag);

			/* send message */
			protocol_send_header(endpoint, sizeof(message));
			protocol_send_payload(message, sizeof(message));
			protocol_send_tail();

			return RADAR_ERR_OK;
		}

		break;

    case MSG_DISABLE_BGT_LNA:
	
        if (num_bytes == 1)
        {    		
        	radar_disable_lna_gain(radar_driver);
			
            return RADAR_ERR_OK;
        }
        break;
		
	case MSG_ENABLE_BGT_LNA:
	
        if (num_bytes == 1)
        {    		
        	radar_enable_lna_gain(radar_driver);

            return RADAR_ERR_OK;
        }
        break;

	case MSG_IS_ENABLE_BGT_LNA:

		if (num_bytes == 1)
		{
			uint8_t message[2];
			uint8_t flag;

			flag = radar_get_lna_gain_enable_status(radar_driver);
			/* compile message */
			wr_payload_u8(message, 0, MSG_IS_ENABLE_BGT_LNA);
			wr_payload_u8(message, 1, flag);

			/* send message */
			protocol_send_header(endpoint, sizeof(message));
			protocol_send_payload(message, sizeof(message));
			protocol_send_tail();

			return RADAR_ERR_OK;
		}

		break;

    default:
	
		break;
    }

    return PROTOCOL_STATUS_INVALID_PAYLOAD;
}

void ep_radar_industrial_handle_change(uint8_t endpoint, void* context,
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
