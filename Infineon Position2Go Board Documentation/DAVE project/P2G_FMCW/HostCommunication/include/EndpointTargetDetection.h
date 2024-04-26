/**
 * \file EndpointTargetDetection.h
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

#ifndef HOSTCOMMUNICATION_ENDPOINTTARGETDETECTION_H_
#define HOSTCOMMUNICATION_ENDPOINTTARGETDETECTION_H_

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
#define EP_TARGET_DETECTION_ENDPOINT_LIST_ENTRY(context) \
{ \
	/*.uEndpointType    = */ 0x52544443, /* ASCII code 'RTDC' */ \
	/*.uEndpointVersion = */ 1, \
	/*.handleMessage    = */ ep_target_detection_handle_message, \
	/*.context         = */ context, \
	/*.handleChange     = */ ep_target_detection_handle_change \
}

/*
==============================================================================
   3. FUNCTION PROTOTYPES AND INLINE FUNCTIONS
==============================================================================
*/

uint16_t ep_target_detection_handle_message(uint8_t endpoint,
                                         uint8_t* message_data,
                                         uint16_t num_bytes, void* context);

void ep_target_detection_handle_change(uint8_t endpoint, void* context,
                                    uint32_t what);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* HOSTCOMMUNICATION_ENDPOINTTARGETDETECTION_H_ */

/* --- End of File -------------------------------------------------------- */
