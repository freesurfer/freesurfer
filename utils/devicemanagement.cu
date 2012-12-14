/*! \file
  File to hold routines which take care of CUDA device management
*/

#include <iostream>
#include <cstdlib>
using namespace std;

#include "nvcc_version.h"

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


void AcquireCUDADevice( void )
{
  /*!
    Looks up the environment variable "FREESURFER_CUDA_DEVICE"
    attempts to acquire that GPU.
    On success, it will print out the name of the acquired GPU
  */
  char *devString;
  int iDevice;

  cout << nvcc_version << endl;
  int driverVersion=0, runtimeVersion=0;
  if( cudaSuccess == cudaDriverGetVersion( &driverVersion ) )
  {
    cudaRuntimeGetVersion( &runtimeVersion );
    cout << "Driver : "
         << driverVersion/1000 << "." << driverVersion%1000
         << endl;
    cout << "Runtime : "
         << runtimeVersion/1000 << "." << runtimeVersion%1000
         << endl;
    cout << endl;
  }

  cout << "Acquiring CUDA device" << endl;

  // Get the environment variable
  devString = getenv( "FREESURFER_CUDA_DEVICE" );
  if( devString == NULL )
  {
    cout << "Using default device" << endl;
  }
  else
  {
    iDevice = atoi( devString );
    cout << "Device " << iDevice << " requested" << endl;
    if( cudaSuccess != cudaSetDevice( iDevice ) )
    {
      cerr << "ERROR: Unable to set CUDA device " << iDevice << endl;
      exit( EXIT_FAILURE );
    }
  }

  // Acquire the device
  int *d_tmp;
  if( cudaSuccess != cudaMalloc( (void**)&d_tmp, sizeof(int) ) )
  {
    cerr << "ERROR: Unable to acquire a CUDA device!" << endl;
    exit( EXIT_FAILURE );
  }
  CUDA_SAFE_CALL( cudaFree( d_tmp ) );

  // Verify and print out device name
  CUDA_SAFE_CALL( cudaGetDevice( &iDevice ) );

  cudaDeviceProp devProp;

  CUDA_SAFE_CALL( cudaGetDeviceProperties( &devProp, iDevice ) );
  cout << "CUDA device: " << devProp.name << endl;

#ifdef GCAMORPH_ON_GPU
  if( devProp.major < 2 )
  {
    cerr << "Must have compute capability 2 for GCAMORPH_ON_GPU"
         << endl;
    exit( EXIT_FAILURE );
  }
#endif

}





// ==============================================


void PrintGPUtimers( void )
{

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
