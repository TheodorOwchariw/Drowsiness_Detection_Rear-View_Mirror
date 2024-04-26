/**
    \file   tracking.c

  	\brief  Open loop tracking with data association
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


/*
==============================================================================
   1. INCLUDE FILES
==============================================================================
 */

#include <stdlib.h>
#include "tracking.h"
#include "dsp_lib.h"
#include "datastore.h"

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
 */

#define  LIFE_TIME_COUNT	(20U)			/**< Number of frames after which track is killed */

#define  GHOST_LIFE_TIME	(5U)			/**< Number of frames after which ghost track is killed */

#define  DELTA_PX_CM		(120.0f)		/**< Size of prediction window, 3 bins */

/*
==============================================================================
   3. DATA
==============================================================================
 */

static Median_Filtering_t median_angle_arr[CURRENT_NUM_OF_TRACKS];

static float rx_angle_fft[2*DOPPLER_FFT_SIZE];

static float rx_angle_fft_spectrum[DOPPLER_FFT_SIZE];

/*
==============================================================================
   4. LOCAL FUNCTIONS
==============================================================================
 */

static uint8_t get_next_free_trackID(tracking_list_t *p_tracks);

static void clear_track_elem(Tracking_Params_t* track_ptr);

static uint32_t assign_track(tracking_list_t *p_tracks, Measurement_elem_t* target_measurement, uint8_t id, uint16_t angle_offset_deg, float wave_length_ant_spacing_ratio, float min_angle);

static void update_track(Tracking_Params_t* track_ptr, Measurement_elem_t* target_measurement, uint8_t track_id,
		float px_track_predict, float range_detection_threshold, int16_t angle_offset_deg, uint32_t num_of_chirps,
		float wave_length_ant_spacing_ratio);

static int compare_float(const void* a, const void* b);

static float median_filtering(Median_Filtering_t *track_median_arr, float new_input);

//===========================================================================

void data_association(tracking_list_t *p_tracks, Measurement_elem_t* target_measurements, uint16_t num_of_targets,
		const algo_settings_t *cp_algo_settings, uint32_t frame_period_usec, uint32_t num_of_chirps )
{

	float px_track_predict = 0;
	float px_measured = 0;
	float px_min, px_max;
	float curr_neighbor, last_neighbor;

	//========= CASE 1: all tracks are empty ===========
	if(p_tracks->num_of_tracks == 0 && num_of_targets != 0)
	{
		// clear all tracks:
		for (uint8_t i = 0; i < p_tracks->max_num_of_tracks; i++)
		{
			memset(&p_tracks->elems[i], 0 ,sizeof(p_tracks->elems[i]) );
		}

		for (uint8_t j = 0; j < num_of_targets; j++)
		{
			if(target_measurements[j].speed != 0)		// only moving targets
			{
				uint8_t next_free_id = get_next_free_trackID( p_tracks );

				(void) assign_track(p_tracks, &target_measurements[j], next_free_id, cp_algo_settings->angle_offset_deg,
						cp_algo_settings->wave_length_ant_spacing_ratio, cp_algo_settings->min_angle_for_track_assignment);

				if (p_tracks->num_of_tracks >= p_tracks->max_num_of_tracks)
				{
					break;
				}
			}
		}
	}

	//========= CASE 2: update tracks already assigned ===========
	else
	{
		int8_t tgt_idx;

		for (uint8_t i = 0; i < p_tracks->max_num_of_tracks; i++)
		{
			if (p_tracks->elems[i].is_alived == 0)
			{
				continue;
			}

			px_track_predict = p_tracks->elems[i].range - (p_tracks->elems[i].speed * ((frame_period_usec/1000) * 0.1f) );	// range is in cm, thus speed is cm/sec

			px_max = px_track_predict + DELTA_PX_CM;

			px_min = px_track_predict - DELTA_PX_CM;

			if (px_min < 0)
			{
				px_min = 0;
			}

			//=========== SNN loop =============

			tgt_idx = -1;

			curr_neighbor = 0;

			last_neighbor = 10000;

			for (uint8_t j = 0; j < num_of_targets; j++)
			{
				if(target_measurements[j].is_associated == 0)
				{
					px_measured = target_measurements[j].range;

					if((px_min <= px_measured) && (px_measured <= px_max))
					{
						curr_neighbor = fabsf(px_measured - px_track_predict);

						if (curr_neighbor < last_neighbor)
						{
							last_neighbor = px_measured;

							tgt_idx = j;
						}
					}
				}
			}

			if (tgt_idx == -1)
			{
				p_tracks->elems[i].lifetime_counter += 1;

				p_tracks->elems[i].range = 0.5 * (p_tracks->elems[i].range + px_track_predict);
			}
			else
			{
				p_tracks->elems[i].lifetime_counter /= 2;

				p_tracks->elems[i].measurement_counter += 1;

				update_track(&p_tracks->elems[i], &target_measurements[tgt_idx], i, px_track_predict * 0.5,
								cp_algo_settings->range_detection_threshold , cp_algo_settings->angle_offset_deg, num_of_chirps,
								cp_algo_settings->wave_length_ant_spacing_ratio);

				target_measurements[tgt_idx].is_associated = 1;
			}
		}

		//========= CASE 3: Kill expired tracks ===========

		int32_t ghost_target_diff = 0;

		for (uint32_t j = 0; j < p_tracks->max_num_of_tracks; j++)
		{
			ghost_target_diff = p_tracks->elems[j].lifetime_counter - p_tracks->elems[j].measurement_counter - GHOST_LIFE_TIME;

			if (p_tracks->elems[j].lifetime_counter > LIFE_TIME_COUNT || ghost_target_diff > 0  ||
					p_tracks->elems[j].speed_count > 2*LIFE_TIME_COUNT )
			{
				clear_track_elem(&p_tracks->elems[j]);  // init track data
				p_tracks->num_of_tracks -= 1;

				median_angle_arr[j].is_full = 0;
			}
		}

		//========= CASE 4: assign new tracks  ===========

		if (p_tracks->num_of_tracks < p_tracks->max_num_of_tracks)
		{
			for (uint8_t j = 0; j < num_of_targets; j++)
			{
				if(target_measurements[j].is_associated == 0 && target_measurements[j].speed != 0)
				{
					uint8_t next_free_id = get_next_free_trackID( p_tracks );

					(void) assign_track(p_tracks, &target_measurements[j], next_free_id, cp_algo_settings->angle_offset_deg,
							cp_algo_settings->wave_length_ant_spacing_ratio, cp_algo_settings->min_angle_for_track_assignment );

					if (p_tracks->num_of_tracks >= p_tracks->max_num_of_tracks)
					{
						break;
					}
				}
			}
		}
	}

}

//===========================================================================

static void clear_track_elem(Tracking_Params_t* track_ptr)
{
	memset(track_ptr, 0, sizeof(Tracking_Params_t));
	track_ptr->d_phi = IGNORE_NAN;
}

//===========================================================================

static void update_track(Tracking_Params_t* track_ptr, Measurement_elem_t* target_measurement, uint8_t track_id,
		float px_track_predict, float range_detection_threshold, int16_t angle_offset_deg, uint32_t num_of_chirps,
		float wave_length_ant_spacing_ratio)
{

	float last_range = track_ptr->range;

	float i1_mean, q1_mean;

	float i2_mean, q2_mean;

	float max_value;

	uint32_t idx;

	float alpha_range ;

	if (fabs(target_measurement->speed) <= 0.3)
	{
		alpha_range = 0.25;
	}
	else if (fabs(target_measurement->speed) > 0.3 && fabs(target_measurement->speed) < 0.5)
	{
		alpha_range = 0.50;
	}
	else if (fabs(target_measurement->speed) >= 0.5 && fabs(target_measurement->speed) < 1.0)
	{
		alpha_range = 0.8;
	}
	else if (fabs(target_measurement->speed) >= 1.0)
	{
		alpha_range = 0.99;
	}
	else
	{
		alpha_range = 1.0/16.0;
	}

	//------------------------- update range -----------------------------------

	track_ptr->range = (1.0f - alpha_range) * track_ptr->range +
			alpha_range * (target_measurement->range * 0.5f + px_track_predict);

	if (track_ptr->range_change_flag == 0)
	{
		if ((track_ptr->range < last_range - 40) || (track_ptr->range > last_range + 40))
		{
			track_ptr->range_change_flag = 1;
		}
	}

	//------------------------- update speed -----------------------------------

	track_ptr->speed = target_measurement->speed;

	if (track_ptr->speed == 0)
	{
		track_ptr->speed_count += 1;
	}
	else
	{
		track_ptr->speed_count = 0;
	}

	//------------------------- update angle -----------------------------------

	track_ptr->strength = target_measurement->strength;

	if (track_ptr->strength > 1.0 * range_detection_threshold)
	{
		if (track_ptr->speed_count == 0)
		{
			track_ptr->rx1_angle_arg_re[0] = target_measurement->rx1_angle_arg_re;

			track_ptr->rx1_angle_arg_im[0] = target_measurement->rx1_angle_arg_im;

			track_ptr->rx2_angle_arg_re[0] = target_measurement->rx2_angle_arg_re;

			track_ptr->rx2_angle_arg_im[0] = target_measurement->rx2_angle_arg_im;

			target_angle_data T;

			T = compute_angle(track_ptr->rx1_angle_arg_re[0], track_ptr->rx1_angle_arg_im[0],
					track_ptr->rx2_angle_arg_re[0], track_ptr->rx2_angle_arg_im[0],track_ptr->d_phi, angle_offset_deg, wave_length_ant_spacing_ratio);

			track_ptr->angle = T.target_angle;
			track_ptr->d_phi = T.d_phi;

			track_ptr->angle = median_filtering(&median_angle_arr[track_id], track_ptr->angle);
		}
		else if ((track_ptr->speed_count % (num_of_chirps + 1)) != 0)
		{

			idx = (track_ptr->speed_count % (num_of_chirps + 1)) - 1;

			track_ptr->rx1_angle_arg_re[idx] = target_measurement->rx1_angle_arg_re;

			track_ptr->rx1_angle_arg_im[idx] = target_measurement->rx1_angle_arg_im;

			track_ptr->rx2_angle_arg_re[idx] = target_measurement->rx2_angle_arg_re;

			track_ptr->rx2_angle_arg_im[idx] = target_measurement->rx2_angle_arg_im;

			target_angle_data T;

			T = compute_angle(track_ptr->rx1_angle_arg_re[idx], track_ptr->rx1_angle_arg_im[idx],
					track_ptr->rx2_angle_arg_re[idx], track_ptr->rx2_angle_arg_im[idx],track_ptr->d_phi, angle_offset_deg, wave_length_ant_spacing_ratio);

			track_ptr->angle = T.target_angle;
			track_ptr->d_phi = T.d_phi;

			track_ptr->angle = median_filtering(&median_angle_arr[track_id], track_ptr->angle);
		}
		else
		{
			//---------------------------------------- Rx1 ------------------------------------------

			compute_fft_signal(track_ptr->rx1_angle_arg_re, track_ptr->rx1_angle_arg_im, num_of_chirps, DOPPLER_FFT_SIZE, 1.0,
					FFT_INPUT_COMPLEX, FFT_SLOW_TIME, &i1_mean, &q1_mean, rx_angle_fft);

			compute_fft_spectrum(rx_angle_fft, DOPPLER_FFT_SIZE, rx_angle_fft_spectrum);

			arm_max_f32(rx_angle_fft_spectrum, (uint32_t)DOPPLER_FFT_SIZE, &max_value, &idx);

			if (idx != 0)
			{
				i1_mean = rx_angle_fft[2*idx + 0];
				q1_mean = rx_angle_fft[2*idx + 1];
			}

			//---------------------------------------- Rx2 ------------------------------------------

			compute_fft_signal(track_ptr->rx2_angle_arg_re, track_ptr->rx2_angle_arg_im, num_of_chirps, DOPPLER_FFT_SIZE, 1.0,
					FFT_INPUT_COMPLEX, FFT_SLOW_TIME, &i2_mean, &q2_mean, rx_angle_fft);

			compute_fft_spectrum(rx_angle_fft, DOPPLER_FFT_SIZE, rx_angle_fft_spectrum);

			arm_max_f32(rx_angle_fft_spectrum, (uint32_t)DOPPLER_FFT_SIZE, &max_value, &idx);

			if (idx != 0)
			{
				i2_mean = rx_angle_fft[2*idx + 0];
				q2_mean = rx_angle_fft[2*idx + 1];
			}

			target_angle_data T;

			T = compute_angle(i1_mean, q1_mean, i2_mean, q2_mean,track_ptr->d_phi, angle_offset_deg, wave_length_ant_spacing_ratio);

			track_ptr->angle = T.target_angle;
			track_ptr->d_phi = T.d_phi;

			track_ptr->angle = median_filtering(&median_angle_arr[track_id], track_ptr->angle);
		}
	}
}

//===========================================================================

static uint32_t assign_track(tracking_list_t *p_tracks, Measurement_elem_t* target_measurement, uint8_t id, uint16_t angle_offset_deg, float wave_length_ant_spacing_ratio, float min_angle)
{

	uint32_t retValue = 0;

	target_angle_data T;
	T = compute_angle(target_measurement->rx1_angle_arg_re, target_measurement->rx1_angle_arg_im,
			target_measurement->rx2_angle_arg_re, target_measurement->rx2_angle_arg_im, IGNORE_NAN, angle_offset_deg, wave_length_ant_spacing_ratio);

	if( fabs(T.target_angle) < min_angle)
	{

		Tracking_Params_t *p_track = &p_tracks->elems[id];

		p_track->track_id = id;
		p_track->is_alived = 1;

		p_track->measurement_counter = 1;

		p_track->range_change_flag = 1;
		p_track->speed_count = 1;

		p_track->strength = target_measurement->strength;
		p_track->range    = target_measurement->range;
		p_track->speed    = target_measurement->speed;

		p_track->rx1_angle_arg_re[0] = target_measurement->rx1_angle_arg_re;
		p_track->rx1_angle_arg_im[0] = target_measurement->rx1_angle_arg_im;
		p_track->rx2_angle_arg_re[0] = target_measurement->rx2_angle_arg_re;
		p_track->rx2_angle_arg_im[0] = target_measurement->rx2_angle_arg_im;

		p_track->angle  = T.target_angle;
		p_track->d_phi  = T.d_phi;

		// calculate the medium filter for the track angle
		median_angle_arr[id].is_full = 0;
		median_filtering(&median_angle_arr[id], p_track->angle);

		target_measurement->is_associated = 1;

		// Inc the number of valid tracks:
		p_tracks->num_of_tracks += 1;
		retValue = 1;
	}
	else
	{
		clear_track_elem(&p_tracks->elems[id]);
		retValue = 0;
	}

	return retValue;
}

//===========================================================================

static uint8_t get_next_free_trackID(tracking_list_t *p_tracks )
{
	uint8_t i = 0;

	while(i < p_tracks->max_num_of_tracks)
	{
		if(p_tracks->elems[i].is_alived == 0)
		{
			break;
		}
		i++;
	}
	return i;
}

//===========================================================================

static int compare_float(const void *a, const void *b)
{
	int retval = 0;

	float a_f = *(float*)a;
	float b_f = *(float*)b;

	if (a_f > b_f)
	{
		retval = 1;
	}
	else if (a_f < b_f)
	{
		retval = -1;
	}

	return retval;
}

//===========================================================================

float median_filtering(Median_Filtering_t *track_median_arr, float new_input)
{

	if( track_median_arr->median_filter_len > MAX_MEDIAN_FILTER_LEN)
		track_median_arr->median_filter_len = MAX_MEDIAN_FILTER_LEN;

	if (track_median_arr->is_full == 0)
	{
		for (uint32_t j = 0; j < track_median_arr->median_filter_len; j++)
		{
			track_median_arr->buffer[j] = new_input;
		}

		track_median_arr->is_full = 1;

		return new_input;
	}
	else
	{
		float sorting_arr[MAX_MEDIAN_FILTER_LEN];

		uint32_t len = MAX_MEDIAN_FILTER_LEN-1;
		if( len > track_median_arr->median_filter_len-1)
			len = track_median_arr->median_filter_len-1;

		for (uint32_t j = 0; j < len; j++)
		{

			sorting_arr[j] = track_median_arr->buffer[j+1];		// shift the array left in order to add new value

			track_median_arr->buffer[j] = track_median_arr->buffer[j+1];

		}

		track_median_arr->buffer[track_median_arr->median_filter_len-1] = new_input;

		sorting_arr[track_median_arr->median_filter_len-1] = new_input;

		qsort(sorting_arr, track_median_arr->median_filter_len, sizeof(float), compare_float);

		return sorting_arr[track_median_arr->median_filter_len / 2];
	}
}

//===========================================================================

void median_filter_init( uint32_t median_filter_len )
{
	for (uint32_t j = 0; j < sizeof(median_angle_arr)/sizeof(median_angle_arr[0]); j++)
	{
		median_angle_arr[j].is_full = 0;
		median_angle_arr[j].median_filter_len = MAX_MEDIAN_FILTER_LEN;
		if( median_filter_len < MAX_MEDIAN_FILTER_LEN)
			median_angle_arr[j].median_filter_len = median_filter_len;
		// The buffer values of the structure is not cleared. The buffer values are initalized, when then filter is used the first time!
	}

}

/* --- End of File -------------------------------------------------------- */
