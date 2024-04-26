/**
 * \file EndpointRadarBase.c
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
#include "EndpointRadarBase.h"
#include "EndpointRadarChangeCodes.h"
#include "PayloadHelper.h"
#include "EndpointRadarADCXMC.h"

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
#define MSG_FRAME_DATA            0x00 /**< A message containing radar frame
                                            data. */
#define MSG_GET_FRAME_DATA        0x01 /**< A message to retrieve frame data.
                                            */
#define MSG_SET_AUTOMATIC_TRIGGER 0x02 /**< A message to configure the
                                            automatic frame trigger. */
#define MSG_ENABLE_TEST_MODE      0x03 /**< A message to configure antenna
                                            test mode. */
#define MSG_GET_DRIVER_VERSION    0x20 /**< A message to query device
                                            information. */
#define MSG_SET_DRIVER_VERSION    0x21 /**< A message containing device
                                            information. */
#define MSG_GET_DEVICE_INFO       0x22 /**< A message to query device
                                            information. */
#define MSG_SET_DEVICE_INFO       0x23 /**< A message containing device
                                            information. */
#define MSG_GET_TEMPRATURE        0x30 /**< A message to query a measured
                                            temperature value. */
#define MSG_SET_TEMPRATURE        0x31 /**< A message containing a measured
                                            temperature value. */
#define MSG_GET_TX_POWER          0x32 /**< A message to query a measured
                                            TX power value. */
#define MSG_SET_TX_POWER          0x33 /**< A message containing a measured
                                            TX power value. */
#define MSG_GET_CHRIP_DURATION    0x34 /**< A message to query the current
                                            chirp duration. */
#define MSG_SET_CHRIP_DURATION    0x35 /**< A message containing the current
                                            chirp duration. */
#define MSG_GET_MIN_INTERVAL      0x36 /**< A message to query the minimum
                                            frame interval. */
#define MSG_SET_MIN_INTERVAL      0x37 /**< A message containing the minimum
                                            frame interval. */
#define MSG_GET_FRAME_FORMAT      0x40 /**< A message to query the frame
                                            format. */
#define MSG_SET_FRAME_FORMAT      0x41 /**< A message containing frame format
                                            parameters. */

#define MSG_GET_FFT_INFO		  0x42 /**< A message containing the FFT configuration
											parameters. */
#define MSG_SET_FFT_INFO		  0x43 /**< A message containing the FFT configuration
											parameters. */
/** @} */

/*
==============================================================================
   3. LOCAL FUNCTION PROTOTYPES
==============================================================================
*/

static void send_driver_version(uint8_t endpoint);

static uint16_t send_frame_data(uint8_t endpoint, Radar_Handle_t radar_driver,
                                uint8_t wait);

static uint16_t send_device_info(uint8_t endpoint,
                                 Radar_Handle_t radar_driver);

static uint16_t send_temperature(uint8_t endpoint,
                                 Radar_Handle_t radar_driver,
                                 uint8_t temp_sensor);

static uint16_t send_tx_power(uint8_t endpoint, Radar_Handle_t radar_driver,
                              uint8_t tx_antenna);

static uint16_t send_chirp_duration(uint8_t endpoint,
                                    Radar_Handle_t radar_driver);

static uint16_t send_min_frame_interval(uint8_t endpoint,
                                        Radar_Handle_t radar_driver);

static uint16_t send_frame_format(uint8_t endpoint,
                                 Radar_Handle_t radar_driver);

static uint16_t send_fft_config(uint8_t endpoint,
                                 Radar_Handle_t radar_driver);

/*
==============================================================================
  4. LOCAL FUNCTIONS
==============================================================================
*/

static void send_driver_version(uint8_t endpoint)
{
    const Driver_Version_t* driver_version = radar_get_driver_version();

    uint8_t message[4];

    /* compile message */
    wr_payload_u8(message, 0, MSG_SET_DRIVER_VERSION);
    wr_payload_u8(message, 1, driver_version->major);
    wr_payload_u8(message, 2, driver_version->minor);
    wr_payload_u8(message, 3, driver_version->revision);

    /* send message */
    protocol_send_header(endpoint, sizeof(message));
    protocol_send_payload(message, sizeof(message));
    protocol_send_tail();
}



static size_t packing_helper_12Bit( uint8_t *p_dest, uint16_t *p_src, size_t numSamples)
{
	/* numBytes needs to be devided by 2 */
	if( (numSamples & 1) != 0)
		return -1;

	size_t   processed_output_bytes = 0;

	while (numSamples >= 2)
	{
		/* pack two samples into one 32 bit word */
		uint32_t packed_word;
		packed_word  = *p_src++;
		packed_word |= *p_src++ << 12;

		/* write packed word into buffer */
		wr_payload_u32(p_dest, 0, packed_word);
		p_dest += 3;

		numSamples -= 2;
		processed_output_bytes += 3;
	}

	return processed_output_bytes;
}



static uint16_t send_frame_data(uint8_t endpoint, Radar_Handle_t radar_driver,
                                uint8_t wait)
{
    uint16_t error_code;

    Frame_Info_t frame_info;

    /* get frame data */
    error_code = radar_get_frame(radar_driver, &frame_info, wait);

    if (error_code == RADAR_ERR_OK)
    {

    	/* Global buffer to have a copy of fresh ADC data to avoid frame corruption at Host */
    	static uint32_t transport_buffer[BSP_MAX_ADC_BUFFER_SIZE_BYTES/sizeof(uint32_t)];

        /* send message containing frame data */
        uint32_t num_samples;
        uint32_t data_size;
        uint8_t message_header[18];
        uint32_t sample_start = 0;

        uint32_t num_used_ant = 2;

        if (frame_info.adc_resolution == 12)
        {
        	/* Selection of Antenna to output is currently only supported for 12 Bit ADC resolution! */
        	if( frame_info.rx_mask == 1 )
			{
				sample_start = 0;
				num_used_ant = 1;

			}
			else if( frame_info.rx_mask == 2 )
			{
				num_used_ant = 1;
				sample_start = frame_info.num_chirps *
						frame_info.num_samples_per_chirp *
						num_used_ant *
						(frame_info.data_format == RADAR_RX_DATA_REAL ? 1 : 2);
			}
        }
        /* setup message header */
        wr_payload_u8 (message_header,  0, MSG_FRAME_DATA);
        wr_payload_u32(message_header,  1, frame_info.frame_number);
        wr_payload_u32(message_header,  5, frame_info.num_chirps);
        wr_payload_u8 (message_header,  9, num_used_ant);
        wr_payload_u32(message_header, 10, frame_info.num_samples_per_chirp);
        wr_payload_u8 (message_header, 14, frame_info.rx_mask);
        wr_payload_u8 (message_header, 15, frame_info.data_format);
        wr_payload_u8 (message_header, 16, frame_info.adc_resolution);
        wr_payload_u8 (message_header, 17, frame_info.interleaved_rx);

        /* calculate size of data */
        /* calculate total number of bits and divide by 8 to get number of
           bytes. Always round upwards */
        num_samples = frame_info.num_chirps *
                      frame_info.num_samples_per_chirp *
                      num_used_ant *
                      (frame_info.data_format == RADAR_RX_DATA_REAL ? 1 : 2);

        data_size = num_samples * frame_info.adc_resolution;
        data_size = (data_size >> 3) + ((data_size & 0x07) ? 1 : 0);

        /* send message header */
        protocol_send_header(endpoint, data_size + sizeof(message_header));
        protocol_send_payload(message_header, sizeof(message_header));

        if (frame_info.adc_resolution == 12)
        {
			uint8_t *destPtr = (uint8_t *)&transport_buffer[0];
			uint8_t *srcPtr  = (uint8_t *)frame_info.sample_data;

			uint32_t channel_inc = frame_info.num_samples_per_chirp*sizeof(uint16_t);
			uint32_t chirp_inc   = channel_inc * frame_info.num_rx_antennas * (frame_info.data_format == RADAR_RX_DATA_REAL ? 1 : 2);

			uint32_t processed_out_bytes = 0;

			for(uint32_t chirp = 0; chirp < frame_info.num_chirps; chirp++)
			{
				for(uint32_t ant_idx= 0; ant_idx < frame_info.num_rx_antennas; ant_idx++)  // max Ant mask!
				{
					uint8_t mask = (1 << ant_idx);
					if( (mask & frame_info.rx_mask) != 0 )
					{
						uint32_t num_chan_per_ant = (frame_info.data_format == RADAR_RX_DATA_REAL ? 1 : 2);
						for(uint32_t ch_idx = 0; ch_idx < num_chan_per_ant; ch_idx++)
						{
							uint32_t in_offset = (chirp_inc * chirp) + (channel_inc*(ant_idx*num_chan_per_ant + ch_idx ));
							size_t out_offset = packing_helper_12Bit(destPtr, (uint16_t *)&srcPtr[in_offset], frame_info.num_samples_per_chirp );
			
							destPtr += out_offset;
							processed_out_bytes += out_offset;
						}
					}
				}
			}
			/* now send the packed data at once */
			protocol_send_payload((uint8_t*)transport_buffer, data_size);
        }

        else
        {
             uint16_t sample_bit_mask = (1 << frame_info.adc_resolution) - 1;
 			 uint16_t sample_bit_stream = 0;
			 uint16_t write_bit_position = 0;

			 for (unsigned i = 0; i < num_samples; ++i)
			 {
				const uint16_t* unpacked_data = &((uint16_t*)frame_info.sample_data)[sample_start];
				 sample_bit_stream |= (unpacked_data[i] & sample_bit_mask) <<
									   write_bit_position;
				 write_bit_position += frame_info.adc_resolution;

				 protocol_send_payload((uint8_t*)&sample_bit_stream,
									   write_bit_position >> 3);
				 sample_bit_stream >>= write_bit_position & 0xF8;
				 write_bit_position &= 0x07;
			 }
			 if (write_bit_position != 0)
				 protocol_send_payload((uint8_t*)&sample_bit_stream, 1);
        }


        protocol_send_tail();

        /* send another message, containing the temperature */
        {
            // Send the temperature, only when the value is changing!
        	static int32_t previous_temp = 0;
            if( frame_info.temperature_001C != previous_temp )
            {

				uint8_t message[6];

				/* compile message */
				wr_payload_u8 (message, 0, MSG_SET_TEMPRATURE);
				wr_payload_u8 (message, 1, 0);
				wr_payload_i32(message, 2, frame_info.temperature_001C);

				/* send message */
				protocol_send_header(endpoint, sizeof(message));
				protocol_send_payload(message, sizeof(message));
				protocol_send_tail();
				previous_temp = frame_info.temperature_001C;
            }
        }
    }
    return error_code;
}

static uint16_t send_device_info(uint8_t endpoint,
                                 Radar_Handle_t radar_driver)
{
    uint16_t error_code;
    Device_Info_t device_info;

    /* read device info from device */
    error_code = radar_get_device_info(radar_driver, &device_info);

    if (error_code == RADAR_ERR_OK)
    {
        uint16_t string_length;
        uint8_t message[17];

        /* compile message */
        wr_payload_u8 (message,  0, MSG_SET_DEVICE_INFO);
        wr_payload_u32(message,  1, device_info.min_rf_frequency_kHz);
        wr_payload_u32(message,  5, device_info.max_rf_frequency_kHz);
        wr_payload_u8 (message,  9, device_info.num_tx_antennas);
        wr_payload_u8 (message, 10, device_info.num_rx_antennas);
        wr_payload_u8 (message, 11, device_info.max_tx_power);
        wr_payload_u8 (message, 12, device_info.num_temp_sensors);
        wr_payload_u8 (message, 13, device_info.major_version_hw);
        wr_payload_u8 (message, 14, device_info.minor_version_hw);
        wr_payload_u8 (message, 15, device_info.interleaved_rx);
        wr_payload_u8 (message, 16, device_info.data_format);

        /* count string length (add 1 for the terminating 0) */
        string_length = 0;
        while (device_info.description[string_length] != 0)
            ++string_length;
        ++string_length;

        /* send message */
        protocol_send_header(endpoint, sizeof(message) + string_length);
        protocol_send_payload(message, sizeof(message));
        protocol_send_payload((const uint8_t*)device_info.description,
                              string_length);
        protocol_send_tail();
    }
    return error_code;
}

static uint16_t send_temperature(uint8_t endpoint,
                                 Radar_Handle_t radar_driver,
                                 uint8_t temp_sensor)
{
    uint16_t error_code;
    int32_t temperature;

    /* read temperature value */
    error_code = radar_get_temperature(radar_driver, temp_sensor,
                                       &temperature);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[6];

        /* compile message */
        wr_payload_u8 (message, 0, MSG_SET_TEMPRATURE);
        wr_payload_u8 (message, 1, temp_sensor);
        wr_payload_i32(message, 2, temperature);

        /* send message */
        protocol_send_header(endpoint, sizeof(message));
        protocol_send_payload(message, sizeof(message));
        protocol_send_tail();
    }
    return error_code;
}

static uint16_t send_tx_power(uint8_t endpoint, Radar_Handle_t radar_driver,
                              uint8_t tx_antenna)
{
    uint16_t error_code;
    int32_t tx_power;

    /* read TX power value value */
    error_code =  radar_get_tx_power(radar_driver, tx_antenna, &tx_power);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[6];

        /* compile message */
        wr_payload_u8 (message, 0, MSG_SET_TX_POWER);
        wr_payload_u8 (message, 1, tx_antenna);
        wr_payload_i32(message, 2, tx_power);

        /* send message */
        protocol_send_header(endpoint, sizeof(message));
        protocol_send_payload(message, sizeof(message));
        protocol_send_tail();
    }
    return error_code;
}

static uint16_t send_chirp_duration(uint8_t endpoint,
                                    Radar_Handle_t radar_driver)
{
    uint16_t error_code;
    uint32_t chirp_duration;

    /* read chirp duration */
    error_code = radar_get_chirp_duration(radar_driver, &chirp_duration);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[5];

        /* compile message */
        wr_payload_u8 (message, 0, MSG_SET_CHRIP_DURATION);
        wr_payload_u32(message, 1, chirp_duration);

        /* send message */
        protocol_send_header(endpoint, sizeof(message));
        protocol_send_payload(message, sizeof(message));
        protocol_send_tail();
    }
    return error_code;
}

static uint16_t send_min_frame_interval(uint8_t endpoint,
                                        Radar_Handle_t radar_driver)
{
    uint16_t error_code;
    uint32_t min_frame_interval;

    /* read minimum frame interval */
    error_code = radar_get_min_frame_interval(radar_driver,
                                              &min_frame_interval);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[5];

        /* compile message */
        wr_payload_u8 (message, 0, MSG_SET_MIN_INTERVAL);
        wr_payload_u32(message, 1, min_frame_interval);

        /* send message */
        protocol_send_header(endpoint, sizeof(message));
        protocol_send_payload(message, sizeof(message));
        protocol_send_tail();
    }
    return error_code;
}

static uint16_t send_frame_format(uint8_t endpoint,
                                 Radar_Handle_t radar_driver)
{
    uint16_t error_code;
    Frame_Format_t frame_format;

    /* get current frame format from radar driver */
    error_code = radar_get_frame_format(radar_driver, &frame_format);

    if (error_code == RADAR_ERR_OK)
    {
        uint8_t message[11];

        /* compile message */
        wr_payload_u8 (message,  0, MSG_SET_FRAME_FORMAT);
        wr_payload_u32(message,  1, frame_format.num_samples_per_chirp);
        wr_payload_u32(message,  5, frame_format.num_chirps_per_frame);
        wr_payload_u8 (message,  9, frame_format.rx_mask);
        wr_payload_u8 (message, 10, frame_format.signal_part);

        /* send message */
        protocol_send_header(endpoint, sizeof(message));
        protocol_send_payload(message, sizeof(message));
        protocol_send_tail();
    }
    return error_code;
}

static uint16_t send_fft_config(uint8_t endpoint,
                                 Radar_Handle_t radar_driver)
{

	FFT_configuration_t fft_conf;
	radar_get_fft_configuration(radar_driver,&fft_conf);

	union {
			uint32_t  u32_value;
			float	  f32_value;
		} convert_helper;
	convert_helper.f32_value = fft_conf.scaling_factor;

    uint8_t message[11];

	/* compile message */
	wr_payload_u8 (message,  0, MSG_SET_FFT_INFO);
	wr_payload_u32(message,  1, fft_conf.fft_window_size);
	wr_payload_u8(message,   5, (uint8_t)fft_conf.window_type);
	wr_payload_u8 (message,  6, fft_conf.mean_removed);
	wr_payload_u32(message,  7, convert_helper.u32_value);

	/* send message */
	protocol_send_header(endpoint, sizeof(message));
	protocol_send_payload(message, sizeof(message));
	protocol_send_tail();

	return RADAR_ERR_OK;
}
/*
==============================================================================
   5. EXPORTED FUNCTIONS
==============================================================================
*/

uint16_t ep_radar_base_handle_message(uint8_t endpoint, uint8_t* message_data,
                                      uint16_t num_bytes, void* context)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;

    switch (message_data[0])
    {
    case MSG_GET_FRAME_DATA:
        if (num_bytes == 2)
        {
            uint8_t wait;

            /* extract parameters from message */
            wait = rd_payload_u8(message_data, 1);

            return send_frame_data(endpoint, radar_driver, wait);
        }
        break;

    case MSG_SET_AUTOMATIC_TRIGGER:
        if (num_bytes == 5)
        {
            unsigned error_code;
            uint32_t frame_interval;

            /* extract parameters from message */
            frame_interval = rd_payload_u32(message_data, 1);

            /* apply new settings */
            error_code = radar_set_automatic_frame_trigger(radar_driver,
                                                           frame_interval);

            return error_code;
        }
        break;

    case MSG_ENABLE_TEST_MODE:
        if (num_bytes == 8)
        {
            unsigned error_code;
            uint8_t tx_mask;
            uint8_t rx_mask;
            uint32_t frequency;
            uint8_t tx_power;

            /* extract parameters from message */
            tx_mask   = rd_payload_u8 (message_data, 1);
            rx_mask   = rd_payload_u8 (message_data, 2);
            frequency = rd_payload_u32(message_data, 3);
            tx_power  = rd_payload_u8 (message_data, 7);

            /* apply new settings */
            error_code = radar_test_antennas(radar_driver, tx_mask, rx_mask,
                                             frequency, tx_power);

            return error_code;
        }
        break;

    case MSG_GET_DRIVER_VERSION:
        if (num_bytes == 1)
        {
            /* send driver version info back to host */
            send_driver_version(endpoint);

            return RADAR_ERR_OK;
        }
        break;

    case MSG_GET_DEVICE_INFO:
        if (num_bytes == 1)
        {
            return send_device_info(endpoint, radar_driver);
        }
        break;

    case MSG_GET_TEMPRATURE:
    if (num_bytes == 2)
        {
            uint8_t temp_sensor;

            /* extract parameters from message */
            temp_sensor = rd_payload_u8(message_data, 1);

            /* send temperature value back to host */
            return send_temperature(endpoint, radar_driver, temp_sensor);
        }
        break;

    case MSG_GET_TX_POWER:
    if (num_bytes == 2)
        {
            uint8_t tx_antenna;

            /* extract parameters from message */
            tx_antenna = rd_payload_u8(message_data, 1);

            return send_tx_power(endpoint, radar_driver, tx_antenna);
        }
        break;

    case MSG_GET_CHRIP_DURATION:
    if (num_bytes == 1)
        {
            return send_chirp_duration(endpoint, radar_driver);
        }
        break;

    case MSG_GET_MIN_INTERVAL:
    if (num_bytes == 1)
        {
            return send_min_frame_interval(endpoint, radar_driver);
        }
        break;

    case MSG_GET_FRAME_FORMAT:
        if (num_bytes == 1)
        {
            return send_frame_format(endpoint, radar_driver);
        }
        break;

    case MSG_SET_FRAME_FORMAT:
        if (num_bytes == 11)
        {
            unsigned error_code;
            Frame_Format_t format;

            /* extract parameters from message */
            format.num_samples_per_chirp = rd_payload_u32(message_data,  1);
            format.num_chirps_per_frame  = rd_payload_u32(message_data,  5);
            format.rx_mask               = rd_payload_u8 (message_data,  9);
            format.signal_part           = rd_payload_u8 (message_data, 10);

            /* apply new frame format */
            error_code = radar_set_frame_format(radar_driver, &format);

            /* inform other endpoints of the change */
            protocol_broadcast_change(radar_driver,
                                      EP_RADAR_CHNG_FRAME_FORMAT);

            return error_code;
        }
        break;
    case MSG_GET_FFT_INFO:
		if (num_bytes == 1)
		{
			return send_fft_config(endpoint, radar_driver);
		}
		break;

    default:
        break;
    }
    return PROTOCOL_STATUS_INVALID_PAYLOAD;
}

void ep_radar_base_handle_change(uint8_t endpoint, void* context,
                                 uint32_t what)
{
    Radar_Handle_t radar_driver = (Radar_Handle_t)context;

    switch (what)
    {
    case EP_RADAR_CHNG_FRAME_FORMAT:
        send_frame_format(endpoint, radar_driver);
        send_chirp_duration(endpoint, radar_driver);
        send_min_frame_interval(endpoint, radar_driver);

        {
        	extern uint8_t adc_endpoint;

        	if( adc_endpoint != 0)
        		send_adc_configuration(adc_endpoint,radar_driver); // change
        }
        break;

    case EP_RADAR_CHNG_FMCW_SETTINGS:
    case EP_RADAR_CHNG_DOPPLER_SETTINGS:
    case EP_RADAR_CHNG_TX_MODE:
        {
            uint8_t tx_antenna = 0;
            uint16_t error_code = RADAR_ERR_OK;

            while (error_code == RADAR_ERR_OK)
            {
                error_code = send_tx_power(endpoint, radar_driver,
                                           tx_antenna++);
            }
        }
        break;

    case EP_RADAR_CHNG_ADC_SETTINGS:
        send_chirp_duration(endpoint, radar_driver);
        send_min_frame_interval(endpoint, radar_driver);
        break;

    default:
        break;
    }
}

/* --- End of File -------------------------------------------------------- */
