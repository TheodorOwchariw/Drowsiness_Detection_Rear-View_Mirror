===============================================================================
Folder Contents
===============================================================================

- P2G_FW.hex
  *  Precompiled FW only project that provides raw Radar ADC data based on Doppler
     or FMCW modulation (configurable).
  *  Use this project to develop your own algorithm.

- P2G_FMCW.hex
  *  Precompiled FW including FW described above (P2G_FW.hex) and in addition 
     adds implementation of FMCW algorithm to measure distance, speed and angle
     including tracking of multiple objects.
  *  Use this project in case you are interested in outdoor people tracking use case.
