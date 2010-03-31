/*! \file
  Header file for CUDA device management routines
*/

#ifndef DEVICE_MANAGEMENT_CUDA_H
#define DEVICE_MANAGEMENT_CUDA_H


#if defined(__cplusplus)
extern "C" {
#endif

  extern int gEnableGPUTiming;

  //! Acquires device designated by FREESURFER_CUDA_DEVICE environment variable
  void AcquireCUDADevice( void );
  
  //! Prints out all the timers
  void PrintGPUtimers( void );


#if defined(__cplusplus)
};
#endif


#endif
