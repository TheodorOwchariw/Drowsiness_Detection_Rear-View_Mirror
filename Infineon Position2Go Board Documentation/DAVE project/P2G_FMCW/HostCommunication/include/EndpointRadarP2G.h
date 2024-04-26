/**
 * \file EndpointRadarP2G.h
 *
 * \brief Extension of radar endpoint for Position2Go demo board
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

#ifndef HOSTCOMMUNICATION_ENDPOINTRADARP2G_H_
#define HOSTCOMMUNICATION_ENDPOINTRADARP2G_H_

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

#define EP_RADAR_P2G_ENDPOINT_LIST_ENTRY(context) \
{ \
    /*.endpoint_type    = */ 0x52503247, /* ASCII code 'RP2G' */ \
                                         /* (= Radar Extension P2G) */ \
    /*.endpoint_version = */ 1, \
    /*.handle_message   = */ ep_radar_p2g_handle_message, \
    /*.context          = */ context, \
    /*.handle_change    = */ ep_radar_p2g_handle_change \
}

/*
==============================================================================
   3. FUNCTION PROTOTYPES AND INLINE FUNCTIONS
==============================================================================
*/

uint16_t ep_radar_p2g_handle_message(uint8_t endpoint,
                                     uint8_t* message_data,
                                     uint16_t num_bytes, void* context);

void ep_radar_p2g_handle_change(uint8_t endpoint, void* context,
                                uint32_t what);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* HOSTCOMMUNICATION_ENDPOINTRADARP2G_H_ */

/* --- End of File -------------------------------------------------------- */
