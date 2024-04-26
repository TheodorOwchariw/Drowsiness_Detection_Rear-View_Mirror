/**
 * \file EndpointRadarIndustrial.h
 *
 * \brief Radar endpoint containing functionality pertinent to Industrial Radar Demo kits like D2G and P2G etc.
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

#ifndef HOSTCOMMUNICATION_ENDPOINTRADARINDUSTRIAL_H_
#define HOSTCOMMUNICATION_ENDPOINTRADARINDUSTRIAL_H_

/*
==============================================================================
   1. INCLUDE FILES
==============================================================================
*/
#include "Protocol.h"
#include "ifxRadar.h"

/* Enable C linkage if header is included in C++ files */
#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
*/

#define EP_RADAR_INDUSTRIAL_ENDPOINT_LIST_ENTRY(context) \
{ \
    /*.endpoint_type    = */ 0x52494E44, /* ASCII code 'RIND' */ \
                                         /* (= Radar Industrial) */ \
    /*.endpoint_version = */ 1, \
    /*.handle_message   = */ ep_radar_industrial_handle_message, \
    /*.context          = */ context, \
    /*.handle_change    = */ ep_radar_industrial_handle_change \
}

/*
==============================================================================
   5. FUNCTION PROTOTYPES AND INLINE FUNCTIONS
==============================================================================
*/

uint16_t ep_radar_industrial_handle_message(uint8_t endpoint, uint8_t* message_data,
											uint16_t num_bytes,	void* context);

void ep_radar_industrial_handle_change(uint8_t endpoint, void* context,
                                  	   uint32_t what);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* HOSTCOMMUNICATION_ENDPOINTRADARINDUSTRIAL_H_ */

/* --- End of File -------------------------------------------------------- */
