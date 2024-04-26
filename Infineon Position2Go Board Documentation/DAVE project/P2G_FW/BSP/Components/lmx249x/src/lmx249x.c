/**
    @file: lmx249x.c

    @brief: This file provides low level driver for LMX2491RTW low noise fractional-N PLL
            It contains a set of functions to configure and use LMX2491RTW PLL.
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

#include "lmx249x.h"

/*
==============================================================================
   2. EXPORTED FUNCTIONS
==============================================================================
 */

void lmx249x_init(Lmx249x_Object_t* pThis, const Lmx249x_Hardware_Setup_t* pSetup,
                  sendSPIFunction sendSPI, void* pDataForSendSPI)
{
  uint8_t uSPIData[17];
  
  /* Set Internal Parameters */
  pThis->dExternalDivideFactor = (1.0 / pSetup->uExternalDivider); /* inverse RF-divider within the BGT chip */
  pThis->dPFDCycleTime = 1.0f / (pSetup->dReferenceFreq * (pSetup->eReferenceDoubler ? 2.0 : 1.0) / pSetup->uReferenceDivider );
  pThis->uReg58 = 0;
  pThis->sendSPI = sendSPI;
  pThis->pDataForSendSPI = pDataForSendSPI;
  
  /* Do a Power on Reset and set all PLL registers to the default value */
  uSPIData[0] = 0;
  uSPIData[1] = 2;
  uSPIData[2] = (1 << 2);
  pThis->sendSPI(uSPIData, 3, pThis->pDataForSendSPI);
  
  /* Setup Routing and Drive for Pins TRIG1, TRIG2, MOD and MUXout */
  uSPIData[LMX249x_REG_IDX(25, 15, 35)] = 0x41 | /* these bits must be set */
    ((pSetup->eTrig1PinFunction   >> 2) & 0x08) | /* TRIG1_MUX[5] shifted to bit 3 */
      ((pSetup->eTrig2PinFunction   >> 1) & 0x10) | /* TRIG2_MUX[5] shifted to bit 4 */
        ((pSetup->eMUXoutPinFunction  >> 0) & 0x20) | /* MUXout_MUX[5] stays at bit 5 */
          ((pSetup->eModPinFunction     << 2) & 0x80);  /* MOD_MUX[5] shifted to at bit 7 */
  
  uSPIData[LMX249x_REG_IDX(25, 15, 36)] = ((pSetup->eTrig1PinFunction & 0x1F) << 3) |
    pSetup->eTrig1PinDriveMode;
  
  uSPIData[LMX249x_REG_IDX(25, 15, 37)] = ((pSetup->eTrig2PinFunction & 0x1F) << 3) |
    pSetup->eTrig2PinDriveMode;
  
  uSPIData[LMX249x_REG_IDX(25, 15, 38)] = ((pSetup->eModPinFunction & 0x1F) << 3) |
    pSetup->eModPinDriveMode;
  
  uSPIData[LMX249x_REG_IDX(25, 15, 39)] = ((pSetup->eMUXoutPinFunction & 0x1F) << 3) |
    pSetup->eMUXoutPinDriveMode;
  
  /* Setup Reference Scaling */
  uSPIData[LMX249x_REG_IDX(25, 15, 25)] = (pSetup->uReferenceDivider >> 0) & 0xFF;
  uSPIData[LMX249x_REG_IDX(25, 15, 26)] = (pSetup->uReferenceDivider >> 8) & 0xFF;
  uSPIData[LMX249x_REG_IDX(25, 15, 27)] =  pSetup->eReferenceDoubler |
    (pSetup->eOscInMode << 2);
  
  /* Setup Charge Pump configuration */
  uSPIData[LMX249x_REG_IDX(25, 15, 27)] |= (pSetup->eChargePumpPulseWidth << 3);
  uSPIData[LMX249x_REG_IDX(25, 15, 28)] = ((pSetup->eChargePumpPolarity) << 5) |
    (pSetup->uChargePumpCurrent & 0x1F);
  uSPIData[LMX249x_REG_IDX(25, 15, 29)] =  (pSetup->uChargePumpCurrentFS & 0x1F);
  
  /* Setup speed up settings */
  uSPIData[LMX249x_REG_IDX(25, 15, 27)] |= (pSetup->eCycleSlipReduction << 5);
  uSPIData[LMX249x_REG_IDX(25, 15, 32)]  = (pSetup->uFastLockTimer & 0xFF);
  uSPIData[LMX249x_REG_IDX(25, 15, 29)] |= (pSetup->uFastLockTimer >> 3) & 0xE0;
  
  /* Setup lock detection */
  uSPIData[LMX249x_REG_IDX(25, 15, 30)] =  pSetup->uChargePumpThresholdLo & 0x3F;
  uSPIData[LMX249x_REG_IDX(25, 15, 31)] =  pSetup->uChargePumpThresholdHi & 0x3F;
  
  uSPIData[LMX249x_REG_IDX(25, 15, 33)] =  pSetup->uLockDetectNumGoodEdge;
  uSPIData[LMX249x_REG_IDX(25, 15, 34)] = (pSetup->uLockDetectNumBadEdge & 0x1F) |
    (pSetup->eLockDetectWindow << 5);
  
  uSPIData[0] = 0;
  uSPIData[1] = 25 + 15 - 1;
  pThis->sendSPI(uSPIData, 15 + 2, pThis->pDataForSendSPI);
}

//============================================================================

void lmx249x_set_power_state(Lmx249x_Object_t* pThis, Lmx249x_Power_State_t eState)
{
  uint8_t uSPIData[3];
  uSPIData[0] = 0;
  uSPIData[1] = 2;
  uSPIData[2] = eState;
  pThis->sendSPI(uSPIData, 3, pThis->pDataForSendSPI);
}

//============================================================================

Lmx249x_Error_Code_t lmx249x_set_frequency(Lmx249x_Object_t* pThis, double dBaseFrequency,
                                           Lmx249x_Frac_Order_t eFracOrder,
                                           Lmx249x_Frac_Dither_t eDitherMode)
{
  uint8_t uSPIData[11]; /* Data buffer that will be passed to the SPI interface */
  double dRelFrequency;
  uint32_t iFactorN;
  int32_t iFracDenominator;
  int32_t iFracNumarator;
  
  /* Disable ramp (just in case a ramp is currently in progress) */
  uSPIData[0] = 0;
  uSPIData[1] = 58;
  uSPIData[2] = 0;
  pThis->sendSPI(uSPIData, 3, pThis->pDataForSendSPI);
  
  /* Setup frequency */
  dRelFrequency = dBaseFrequency * pThis->dExternalDivideFactor * pThis->dPFDCycleTime;  /* divider ratio between PLL-RF-in and PFD-frequency */
  
  iFactorN = (int32_t)dRelFrequency;
  dRelFrequency -= iFactorN;
  iFracDenominator = 1<<24;
  iFracNumarator = (uint32_t)(dRelFrequency * iFracDenominator + 0.5);
  
  uSPIData[LMX249x_REG_IDX(16, 9, 16)] = (int8_t) ((iFactorN >>  0) & 0xFF);
  uSPIData[LMX249x_REG_IDX(16, 9, 17)] = (int8_t) ((iFactorN >>  8) & 0xFF);
  uSPIData[LMX249x_REG_IDX(16, 9, 18)] = (int8_t) ((iFactorN >> 16) & 0x03) |
    (eFracOrder << 4)        |
      (eDitherMode << 2);
  
  uSPIData[LMX249x_REG_IDX(16, 9, 19)] = (iFracNumarator >>  0) & 0xFF;
  uSPIData[LMX249x_REG_IDX(16, 9, 20)] = (iFracNumarator >>  8) & 0xFF;
  uSPIData[LMX249x_REG_IDX(16, 9, 21)] = (iFracNumarator >> 16) & 0xFF;
  
  iFracDenominator -= 1;
  uSPIData[LMX249x_REG_IDX(16, 9, 22)] = (iFracDenominator >>  0) & 0xFF;
  uSPIData[LMX249x_REG_IDX(16, 9, 23)] = (iFracDenominator >>  8) & 0xFF;
  uSPIData[LMX249x_REG_IDX(16, 9, 24)] = (iFracDenominator >> 16) & 0xFF;
  
  /* Send register sequence to PLL */
  uSPIData[0] = 0;
  uSPIData[1] = 16 + 9 - 1;
  pThis->sendSPI(uSPIData, 9 + 2, pThis->pDataForSendSPI);
  
  return (LMX249x_ERROR_CODE_OK);
}

//============================================================================

Lmx249x_Error_Code_t lmx249x_configure_ramps(Lmx249x_Object_t* pThis,
                                             const Lmx249x_Ramp_Global_t* pGlobalSettings,
                                             const Lmx249x_Ramp_Section_t* pRampSections,
                                             uint8_t uNumSections)
{
  uint8_t uSPIData[29];  /* Data buffer that will be passed to the SPI interface */
  const double dFrequencyToNFactor = pThis->dExternalDivideFactor * pThis->dPFDCycleTime;
  const double dFracDenominator = (double)(1<<24);
  uint8_t uComparator0Enable = 0;
  uint8_t uComparator1Enable = 0;
  int8_t uIdx;
  uint16_t uBaseRegister;
  const uint8_t uNumRegs = 7;
  const Lmx249x_Ramp_Section_t* pThisSection;
  uint32_t uRampLength;
  uint8_t uDelayFlag = 0;
  int32_t uCounterInc;
  uint8_t uCompEnMask;
  double dRelFrequency;
  uint32_t iFactorN;
  int32_t iFracNumarator;
  int64_t iRampComp0;
  int64_t iRampComp1;
  int64_t iRampLimitLow;
  int64_t iRampLimitHigh;
  int64_t iFSKDev;
  
  /* Check if the number of range is in a valid range */
  if ((uNumSections < 1) || (uNumSections > 8))
    return (LMX249x_ERROR_CODE_INVALID_NUMBER_OF_RAMPS);
  
  /* Disable ramp (just in case a ramp is currently in progress) */
  uSPIData[0] = 0;
  uSPIData[1] = 58;
  uSPIData[2] = 0;
  pThis->sendSPI(uSPIData, 3, pThis->pDataForSendSPI);
  
  /* Setup up the ramp sections */
  for (uIdx = 0; uIdx < uNumSections; ++uIdx)
  {
    /* Setup register buffer */
    uBaseRegister = (86 + uIdx * 7);
    pThisSection = &pRampSections[uIdx];
    
    /* Convert given ramp parameters to counter values */
    uRampLength =  (int32_t) ((pThisSection->dTramp) / pThis->dPFDCycleTime);
    if (uRampLength > 0xFFFF)
    {
      /* If the ramp is too long, divide ramp length by two and set the delay flag which doubles the ramp time */
      uRampLength >>= 1;
      uDelayFlag = 0x80;
    }
    
    /* If the transition to the next sections is triggered by the length of this section, the length must be
    * at least 1, otherwise the counter seems to wrap around and the section will be longer than expected.
    */
    if ((uRampLength == 0) && (pThisSection->eNextTrig == LMX249x_RAMP_NEXT_TRIG_RAMPX_LEN))
      uRampLength = 1;
    
    uCounterInc =  (int32_t) (pThisSection->dFreqShift * dFrequencyToNFactor * dFracDenominator / (double)uRampLength);
    
    /* Set ramp_increment */
    uSPIData[LMX249x_REG_IDX(0, uNumRegs, 0)] = (uint8_t) ((uCounterInc >>  0) & 0xFF);
    uSPIData[LMX249x_REG_IDX(0, uNumRegs, 1)] = (uint8_t) ((uCounterInc >>  8) & 0xFF);
    uSPIData[LMX249x_REG_IDX(0, uNumRegs, 2)] = (uint8_t) ((uCounterInc >> 16) & 0xFF);
    
    /* Set flags and increment */
    uSPIData[LMX249x_REG_IDX(0, uNumRegs, 3)] = uDelayFlag |
      (uint8_t) ((pThisSection->eFastlock) << 6) |
        (uint8_t) ((uCounterInc >> 24) & 0x3F);
    
    /* Set ramp_length */
    uSPIData[LMX249x_REG_IDX(0, uNumRegs, 4)] = ((uRampLength >> 0) & 0xFF);
    uSPIData[LMX249x_REG_IDX(0, uNumRegs, 5)] = ((uRampLength >> 8) & 0xFF);
    
    /* Set flags */
    uSPIData[LMX249x_REG_IDX(0, uNumRegs, 6)] = ((int8_t) ((pThisSection->uNext &0x07)  << 5)) |
      ((int8_t)  (pThisSection->eNextTrig << 3)) |
        ((int8_t)  (pThisSection->eReset    << 2)) |
          ((int8_t)  (pThisSection->eFlag   << 0));
    
    /* Set comparator enable bit */
    uCompEnMask = 1 << uIdx;
    uComparator0Enable |= (pThisSection->eComparators & LMX249x_RAMP_USE_COMPARATOR_1) ? uCompEnMask : 0;
    uComparator1Enable |= (pThisSection->eComparators & LMX249x_RAMP_USE_COMPARATOR_1) ? uCompEnMask : 0;
    
    /* Write the register data to the chip, write highest address (the one of the register written first) at the
    * beginning of the data buffer.
    */
    uBaseRegister += uNumRegs - 1;
    uSPIData[0] = (uint8_t)((uBaseRegister >> 8) & 0xFF);
    uSPIData[1] = (uint8_t)((uBaseRegister >> 0) & 0xFF);
    pThis->sendSPI(uSPIData, uNumRegs + 2, pThis->pDataForSendSPI);
  }
  
  /* Wet comparator enable bits */
  uSPIData[LMX249x_REG_IDX(58, 27, 64)] = uComparator0Enable;
  uSPIData[LMX249x_REG_IDX(58, 27, 69)] = uComparator1Enable;
  
  /* Setup base frequency */
  dRelFrequency = pGlobalSettings->dBaseFrequency * dFrequencyToNFactor;  /* divider ratio between PLL-RF-in and PFD-frequency */
  
  iFactorN = (int32_t)dRelFrequency;
  dRelFrequency -= iFactorN;
  iFracNumarator = (uint32_t) (dRelFrequency * dFracDenominator);
  
  uSPIData[LMX249x_REG_IDX(16, 9, 16)] = (int8_t) ((iFactorN >>  0) & 0xFF);
  uSPIData[LMX249x_REG_IDX(16, 9, 17)] = (int8_t) ((iFactorN >>  8) & 0xFF);
  
  uSPIData[LMX249x_REG_IDX(16, 9, 18)] = (int8_t) (((iFactorN >> 16) & 0x03) | (pGlobalSettings->eFracOrder << 4) | (pGlobalSettings->eDitherMode << 2));
  
  uSPIData[LMX249x_REG_IDX(16, 9, 19)] = (iFracNumarator >>  0) & 0xFF;
  uSPIData[LMX249x_REG_IDX(16, 9, 20)] = (iFracNumarator >>  8) & 0xFF;
  uSPIData[LMX249x_REG_IDX(16, 9, 21)] = (iFracNumarator >> 16) & 0xFF;
  
  /* Fractional denominator is always 0xFFFFFF when ramp generator is active */
  uSPIData[LMX249x_REG_IDX(16, 9, 22)] = 0xFF;
  uSPIData[LMX249x_REG_IDX(16, 9, 23)] = 0xFF;
  uSPIData[LMX249x_REG_IDX(16, 9, 24)] = 0xFF;
  
  /*Send register sequence to PLL */
  uSPIData[0] = 0;
  uSPIData[1] = 16 + 9 - 1;
  pThis->sendSPI(uSPIData, 9 + 2, pThis->pDataForSendSPI);
  
  /* Setup ramp range and comparator values */
  /* The formula in the LMX249x data sheet for the following settings is very misleading. The used formular was
  * found by debugging. The essence here is to specify the limits as the difference (limFreq - baseFreq), while
  * limFreq is the N factor with fractional part, but baseFreq must be the integer N factor (no rounding, just truncating).
  */
  iRampComp0     = (int64_t) ((pGlobalSettings->dComp0Freq    * dFrequencyToNFactor - iFactorN) * dFracDenominator);
  iRampComp1     = (int64_t) ((pGlobalSettings->dComp1Freq    * dFrequencyToNFactor - iFactorN) * dFracDenominator);
  iRampLimitLow  = (int64_t) ((pGlobalSettings->dMinFrequency * dFrequencyToNFactor - iFactorN) * dFracDenominator);
  iRampLimitHigh = (int64_t) ((pGlobalSettings->dMaxFrequency * dFrequencyToNFactor - iFactorN) * dFracDenominator);
  
  uSPIData[LMX249x_REG_IDX(58, 27, 60)] = (iRampComp0 >>  0) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 61)] = (iRampComp0 >>  8) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 62)] = (iRampComp0 >> 16) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 63)] = (iRampComp0 >> 24) & 0xFF;
  
  uSPIData[LMX249x_REG_IDX(58, 27, 65)] = (iRampComp1 >>  0) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 66)] = (iRampComp1 >>  8) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 67)] = (iRampComp1 >> 16) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 68)] = (iRampComp1 >> 24) & 0xFF;
  
  uSPIData[LMX249x_REG_IDX(58, 27, 75)] = (iRampLimitLow >>  0) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 76)] = (iRampLimitLow >>  8) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 77)] = (iRampLimitLow >> 16) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 78)] = (iRampLimitLow >> 24) & 0xFF;
  
  uSPIData[LMX249x_REG_IDX(58, 27, 79)] = (iRampLimitHigh >>  0) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 80)] = (iRampLimitHigh >>  8) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 81)] = (iRampLimitHigh >> 16) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 82)] = (iRampLimitHigh >> 24) & 0xFF;
  
  uSPIData[LMX249x_REG_IDX(58, 27, 70)] = ((iRampComp0    & 0x100000000) ? 1 << 0 : 0) |
    ((iRampComp1    & 0x100000000) ? 1 << 1 : 0) |
      ((iRampLimitLow   & 0x100000000) ? 1 << 3 : 0) |
        ((iRampLimitHigh  & 0x100000000) ? 1 << 4 : 0);
  
  /* Define FSK deviation */
  iFSKDev = (int64_t) ((pGlobalSettings->dDeviationFrequency - pGlobalSettings->dBaseFrequency)
                       * dFrequencyToNFactor * (1 << 24));
  
  uSPIData[LMX249x_REG_IDX(58, 27, 71)] = (iFSKDev >>  0) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 72)] = (iFSKDev >>  8) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 73)] = (iFSKDev >> 16) & 0xFF;
  uSPIData[LMX249x_REG_IDX(58, 27, 74)] = (iFSKDev >> 24) & 0xFF;
  
  uSPIData[LMX249x_REG_IDX(58, 27, 70)] |= ((iFSKDev    & 0x100000000) ? 1 << 2 : 0)
    | (pGlobalSettings->eDevTrigger << 5);
  
  /* Define Ramp Trigger Sources and Modulation type */
  pThis->uReg58 = (pGlobalSettings->eRampClock) /* 0 << 1) */  | /* internal or external clock */
    (pGlobalSettings->eModulation) /* 0 << 2) */ | /* phase or frequency modulation */
      (pGlobalSettings->eTriggerA << 4);          /* Trigger A source */
  uSPIData[LMX249x_REG_IDX(58, 27, 58)] = pThis->uReg58;
  
  uSPIData[LMX249x_REG_IDX(58, 27, 59)] = (pGlobalSettings->eTriggerB)   | /* Trigger B source */
    (pGlobalSettings->eTriggerC << 4);   /* Trigger C source */
  
  /* Setup ramp counter */
  uSPIData[LMX249x_REG_IDX(58, 27, 83)] = (pGlobalSettings->uNumRamps & 0xFF);
  uSPIData[LMX249x_REG_IDX(58, 27, 84)] = ((pGlobalSettings->uNumRamps >> 8)& 0x1F)  |
    (pGlobalSettings->eAutoOff << 5)     |
      (pGlobalSettings->eRampCountTrigger << 6);
  
  /* Send register sequence to PLL */
  uSPIData[0] = 0;
  uSPIData[1] = 58 + 27 - 1;
  pThis->sendSPI(uSPIData, 27 + 2, pThis->pDataForSendSPI);
  
  return (LMX249x_ERROR_CODE_OK);
}

//============================================================================

void lmx249x_enable_ramps(Lmx249x_Object_t* pThis, uint8_t bEnable)
{
  uint8_t uSPIData[4];
  uSPIData[0] = 0;              /* High address byte */
  uSPIData[1] = 58;             /* Low address byte */
  uSPIData[2] = pThis->uReg58;  /* Configuration of register 58 */
  
  /* enable ramps */
  if (bEnable != 0)
    uSPIData[2] |= 1;
  else
    uSPIData[2] &= 0xFE;
  
  pThis->sendSPI(uSPIData, 3, pThis->pDataForSendSPI);
}

//============================================================================

double lmx249x_get_real_frequency_shift(Lmx249x_Object_t* pThis, double FreqShift_MHz, uint32_t FreqShift_time_usec)
{
  double RealFreqShift;
  const double FracDenominator = (double)(1<<24);
  uint32_t FreqShiftSteps;
  double FreqShiftPerStep;
  
  FreqShiftSteps =  (uint32_t) (FreqShift_time_usec / pThis->dPFDCycleTime);
  
  FreqShiftPerStep = (FreqShift_MHz * pThis->dExternalDivideFactor) / FreqShiftSteps;
  
  FreqShiftPerStep = ((uint32_t)((FreqShiftPerStep * FracDenominator) * pThis->dPFDCycleTime) ) / FracDenominator / pThis->dPFDCycleTime;
  
  RealFreqShift = FreqShiftPerStep * FreqShiftSteps /  pThis->dExternalDivideFactor;
  
  return (RealFreqShift);
}

//============================================================================

void lmx249x_trigger_ramp(void)
{
  /* Wait until PLL is ready to generate the next chirp */
  while (DIGITAL_IO_GetInput(&DIGITAL_IO_PLL_MOD) != 0);
  
  /* Start the chirp */
  DIGITAL_IO_SetOutputHigh(&DIGITAL_IO_PLL_TRIG2);
}

//============================================================================

void lmx249x_release_ramp_trigger(void)
{
  DIGITAL_IO_SetOutputLow(&DIGITAL_IO_PLL_TRIG2);
}

//============================================================================

void lmx249x_power_up(void)
{
  DIGITAL_IO_SetOutputHigh(&DIGITAL_IO_PLL_CE);
}

//============================================================================

void lmx249x_power_down(void)
{
  DIGITAL_IO_SetOutputLow(&DIGITAL_IO_PLL_CE);
}

/* --- End of File -------------------------------------------------------- */
