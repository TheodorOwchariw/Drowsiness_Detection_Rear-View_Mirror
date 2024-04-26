/**
    @file: bsp_common.h

    @brief: This file defines the API of the SPI communication protocol.
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

#ifndef BSP_COMMON_H_
#define BSP_COMMON_H_

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

#include "DAVE.h"

/*
==============================================================================
   2. DEFINITIONS
==============================================================================
*/

#define DISABLED	        (0U)
#define ENABLED		        (1U)

#define DOPPLER_MODULATION	(0U)
#define FMCW_MODULATION		(1U)

/*
==============================================================================
   3. TYPES
==============================================================================
 */

/**
 * \brief Select SPI device for read/write operation.
 * @{
 */
typedef enum
{
  BGT_DEVICE               = 0U,   /**< Select BGT for SPI read/write operation */
  LMX_DEVICE               = 1U,   /**< Select LMX for SPI read/write operation */
  PGA_DEVICE               = 2U    /**< Select PGA for SPI read/write operation */
} SPI_Device_Type_t;
/** @} */

/**
 * @brief Global structure contains the SPI Digital IOs used to configure PGA, LMX and BGT devices.
 */
typedef struct
{
  SPI_Device_Type_t spi_device;    /**< SPI device type, refer to SPI_Device_Type_t */
  DIGITAL_IO_t* gpio_cs_pin;       /**< Data structure for SPI_CS GPIO pin */
  DIGITAL_IO_t* gpio_data_pin;     /**< Data structure for SPI_DATA GPIO pin */
  DIGITAL_IO_t* gpio_clk_pin;      /**< Data structure for SPI_CLK GPIO pin */
} SPI_DIGITAL_IOs_t;

/**
 * \brief Returned status of PGA functions.
 * @{
 */
typedef enum
{
  PGA_STATUS_SUCCESS  = 0U,        /**< Success */
  PGA_STATUS_FAIL     = 1U,        /**< Fail */
  PGA_STATUS_NA       = 2U         /**< Feature not available */
} Pga_Status_t;
/** @} */

/*
==============================================================================
   4. FUNCTION PROTOTYPES
==============================================================================
 */

/**
 * \brief  This function is used to transmit the data to the BGT, PGA and LMX devices through SPI protocol.
 *
 * \param[in]	*data_ptr		Pointer to unsigned 8-bit memory, which is to be transmitted via SPI protocol
 * \param[in]	num_of_bytes	Unsigned 8-bit value, representing the number of data bytes to be transmitted via SPI protocol
 * \param[in]	*peripheral_ptr	Void pointer pointing to the target SPI peripheral
 *
 */
void spi_transmit_data(uint8_t* data_ptr, uint8_t num_of_bytes, void* peripheral_ptr);

/**
 * \brief  This function reads out LMX PLL registers.
 *
 * \param[in]  address  16-bit address value, pointing to the first register to start reading out
 * \param[in]  num_of_bytes  Number of consecutive bytes to be read out starting form the given address
 * \param[in]  *peripheral_ptr	Void pointer pointing to the target SPI peripheral
 * \param[out] *out_data_ptr  Pointer of type uint8_t to a memory location where output bytes to be written
 *
 */
void spi_read_register(uint16_t address, uint8_t num_of_bytes, void* peripheral_ptr, uint8_t* out_data_ptr);

/* --- Close open blocks -------------------------------------------------- */

/* Disable C linkage for C++ files */
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* BSP_COMMON_H_ */

/* --- End of File -------------------------------------------------------- */
