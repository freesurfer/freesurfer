/*! \file
  File to hold routines which take care of CUDA device management
*/

#include <iostream>
#include <cstdlib>
using namespace std;

#include "cudacheck.h"

#include "devicemanagement.h"

#include "mriconvolve_cuda.hpp"
#include "mrimean_cuda.hpp"
#include "mrivol2vol_cuda.hpp"

#ifdef GCAMORPH_ON_GPU
#include "gcamorphenergy.hpp"
#include "gcamorphgpu.hpp"
#include "gcamorphtermgpu.hpp"
#include "gcamorphcpu.hpp"
#endif

#include "mrilabels_cuda.hpp"

// ==============================================


void AcquireCUDADevice( void ) {
  /*!
    Looks up the environment variable "FREESURFER_CUDA_DEVICE"
    attempts to acquire that GPU.
    On success, it will print out the name of the acquired GPU
  */
  char *devString;
  int iDevice;

  cout << "Acquiring CUDA device" << endl;

  // Get the environment variable
  devString = getenv( "FREESURFER_CUDA_DEVICE" );
  if( devString == NULL ) {
    cout << "Using default device" << endl;
  } else {
    iDevice = atoi( devString );
    cout << "Device " << iDevice << " requested" << endl;
    CUDA_SAFE_CALL( cudaSetDevice( iDevice ) );
  }

  // Acquire the device
  int *d_tmp;
  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_tmp, sizeof(int) ) );
  CUDA_SAFE_CALL( cudaFree( d_tmp ) );

  // Verify and print out device name
  CUDA_SAFE_CALL( cudaGetDevice( &iDevice ) );

  cudaDeviceProp devProp;
  
  CUDA_SAFE_CALL( cudaGetDeviceProperties( &devProp, iDevice ) );
  cout << "CUDA device: " << devProp.name << endl;
  
}





// ==============================================


void PrintGPUtimers( void ) {

  GPU::Algorithms::MRIconvolve::ShowTimings();
  GPU::Algorithms::MRImean::ShowTimings();
  GPU::Algorithms::MRIvol2vol::ShowTimings();
#if GCAMORPH_ON_GPU
  GPU::Classes::GCAmorphGPU::ShowTimings();
  GPU::Algorithms::GCAmorphEnergy::ShowTimings();
  GPU::Algorithms::GCAmorphTerm::ShowTimings();
  Freesurfer::GCAmorphCPU::ShowTimings();
#endif
  GPU::Algorithms::MRIlabels::ShowTimings();
}
