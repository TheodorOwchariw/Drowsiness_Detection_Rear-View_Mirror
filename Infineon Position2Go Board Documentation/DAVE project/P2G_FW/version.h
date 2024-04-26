/**
    @file: version.h

  	@brief: This file includes the versions for Firmware revision.
*/

/* ===========================================================================
** Copyright (C) 2018-2019 Infineon Technologies AG
** All rights reserved.
** ===========================================================================
**
** ===========================================================================
** This document contains proprietary information of Infineon Technologies AG.
** Passing on and copying of this document, and communication of its contents
** is not permitted without Infineon's prior written authorization.
** ===========================================================================
*/

#ifndef VERSION_H_
#define VERSION_H_

/* Enable C linkage if header is included in C++ files */
#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */
 	
/*
==============================================================================
   1. DEFINITIONS
==============================================================================
*/

#define DEMO_DESCRIPTION            "Position2Go BGT24MTR12 Demo Radar Board"


#define FW_VERSION_MAJOR  			(1U)			/**< Firmware major version */
#define FW_VERSION_MINOR  			(0U)			/**< Firmware minor version */
#define FW_VERSION_REVISION 		(3U)			/**< Firmware patch version */


/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* VERSION_H_ */

/* --- End of File -------------------------------------------------------- */
