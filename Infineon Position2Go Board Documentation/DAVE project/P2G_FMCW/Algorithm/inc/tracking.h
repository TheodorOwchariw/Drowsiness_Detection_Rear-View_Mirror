/**
    @file: tracking.h

    @brief: This file contains open loop tracking with data association.
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

#ifndef TRACKING_H_
#define TRACKING_H_

/* Enable C linkage if header is included in C++ files */
#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*
==============================================================================
   1. INCLUDE FILES
==============================================================================
*/

#include <stdint.h>
#include <stdbool.h>
#include "config.h"
#include "datastore.h"
  
/*
==============================================================================
   2. DEFINITIONS
==============================================================================
*/

/*
==============================================================================
   3. TYPES
==============================================================================
*/
  
/**
 * \brief Data structure for Median filtering.
 * @{
 */
typedef struct
{
  uint8_t  is_full;
  uint32_t median_filter_len;
  float    buffer[MAX_MEDIAN_FILTER_LEN];
} Median_Filtering_t;

/**
 * \brief Data structure for current measurements used in data association.
 * @{
 */
typedef struct
{
  uint16_t is_associated;
  float    strength;
  float    range;
  float    speed;
  float    rx1_angle_arg_re;
  float    rx1_angle_arg_im;
  float    rx2_angle_arg_re;
  float    rx2_angle_arg_im;
} Measurement_elem_t;

/**
 * \brief Data structure for Track parameters used in tracking.
 * @{
 */
typedef struct
{
  uint8_t  track_id;
  uint8_t  is_alived;
  uint16_t speed_count;
  uint16_t range_change_flag;
  uint16_t lifetime_counter;
  uint32_t measurement_counter;
  float    strength;
  float    range;
  float    speed;
  float    angle;
  float    rx1_angle_arg_re[NUM_OF_CHIRPS];
  float    rx1_angle_arg_im[NUM_OF_CHIRPS];
  float    rx2_angle_arg_re[NUM_OF_CHIRPS];
  float    rx2_angle_arg_im[NUM_OF_CHIRPS];
  float    d_phi;
} Tracking_Params_t;

/**
 * \brief Data structure for Tracking List.
 * @{
 */
typedef struct
{
  uint32_t num_of_tracks;
  uint32_t max_num_of_tracks;
  Tracking_Params_t elems[CURRENT_NUM_OF_TRACKS];
} tracking_list_t;

/*
==============================================================================
   4. EXPORTED FUNCTION PROTOTYPES
==============================================================================
*/

/**
 * \brief  Assigns new measurements to the tracks and initializes new tracks.
 */

void data_association(tracking_list_t *p_tracks, Measurement_elem_t* target_measurements,
                      uint16_t num_of_targets, const algo_settings_t *cp_algo_settings, uint32_t frame_period_usec, uint32_t num_of_chirps);

/**
 * \brief  Resets the Median filter for all tracks in case Median filter length changes.
 */
void median_filter_init(uint32_t median_filter_len);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* TRACKING_H_ */

/* --- End of File -------------------------------------------------------- */
