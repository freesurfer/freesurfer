/**
 * @file cudadetect.c
 * @brief Detects whether cuda library is present in the system.
 *
 * This program checks whether cuda libraries are present in the system.
 * If it finds the shared library, it queries for the available total RAM.
 * If the program succeeds it returns 0. otherwise exits with 1.
 * If devices are found, then device info is printed.
 */
/*
 * Original Author: Krish Subramaniam
 * CVS Revision Info:
 * $Id: cudadetect.cpp,v 1.9 2011/03/26 19:43:19 nicks Exp $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 * Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and
 * proprietary rights in and to this software and related documentation and
 * any modifications thereto.  Any use, reproduction, disclosure, or
 * distribution of this software and related documentation without an
 * express license agreement from NVIDIA Corporation is strictly prohibited.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

#ifdef FS_CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>
static int dump();
// Defines for GPU Architecture types 
// (using the SM version to determine the # of cores per SM
static int nGpuArchCoresPerSM[] = { -1, 8, 32 };
// end of GPU Architecture definitions

#endif /* FS_CUDA */

#ifdef __linux
const char* cudalibname = "libcudart.so";
#endif

#ifdef __APPLE__
const char* cudalibname = "libcudart.dylib";
#endif

/* we make use of dlopen and dlclose calls of Linux and MacOSX to check whether
 * the libraries exist. Then we make sure one of the basic calls in that
 * library works.
 * This is a small sanity test to make sure the library isn't corrupted.
 *
 * CUDA is only supported in Linux, MacOSX and Windows. Freesurfer doesn't
 * have Windows support.
 * So, only Linux and MacOSX libraries are queried. both support dlopen(),
 * so there is no further preprocessor conditionals
 */

int main(int argc, char **argv)
{
  void* handle;
  char* error;
  void (*pGetDeviceCount)(int *);
  int number ;

  printf("Detecting CUDA... ");
  handle = dlopen(cudalibname, RTLD_LAZY);
  if (!handle)
  {
    printf("Could not open CUDA lib. *** CUDA not detected! ***\n");
    exit(1);
  }

  /* pGetDeviceCount is the function handle for the function cudaGetDeviceCount
   * in the library*/
  pGetDeviceCount = (void (*)(int *))dlsym(handle, "cudaGetDeviceCount");
  if ((error = dlerror()) != NULL)
  {
    fprintf (stderr, "%s\n", error);
    fprintf (stderr, "libcudart doesn't seem to have cudaGetDeviceCount "
             "function. It's probably corrupted.\n");
    exit(1);
  }
  (*pGetDeviceCount)(&number);
  dlclose(handle);

  if (number==0)
  {
    printf("*** No CUDA enabled device(s) detected! ***\n");
    return 1;
  }

#ifdef FS_CUDA
  // print device info (if cuda is installed)
  return dump();
#else
  printf("%d CUDA enabled device(s) detected.\n",number);
  return 0;
#endif /* FS_CUDA */
}


#ifdef FS_CUDA
/////////////////////////////////////////////////////////////////////////////
// Device info dump (this code is taken from the SDK's deviceQuery example)
/////////////////////////////////////////////////////////////////////////////
static int dump()
{
  int deviceCount = 0;

  // This function call returns 0 if there are no CUDA capable devices.
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0)
  {
    printf("There is no device supporting CUDA\n");
    return (1);
  }

  int dev = 0;
  for (dev = 0; dev < deviceCount; ++dev)
  {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);

    if (dev == 0)
    {
      // This function call returns 9999 for both major & minor fields,
      // if no CUDA capable devices are present
      if (deviceProp.major == 9999 && deviceProp.minor == 9999)
      {
        printf("There is no device supporting CUDA.\n");
        return (1);
      }
      else if (deviceCount == 1)
        printf("There is 1 device supporting CUDA:\n");
      else
        printf("There are %d devices supporting CUDA:\n", deviceCount);
    }
    printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

    int driverVersion = 0, runtimeVersion = 0;
    cudaDriverGetVersion(&driverVersion);
    printf("  CUDA Driver Version:                           %d.%d\n",
           driverVersion/1000, driverVersion%100);
    cudaRuntimeGetVersion(&runtimeVersion);
    printf("  CUDA Runtime Version:                          %d.%d\n",
           runtimeVersion/1000, runtimeVersion%100);

    printf("  CUDA Capability Major revision number:         %d\n",
           deviceProp.major);
    printf("  CUDA Capability Minor revision number:         %d\n",
           deviceProp.minor);

    printf("  Total amount of global memory:                 %u bytes\n",
           (unsigned int)deviceProp.totalGlobalMem);

    printf("  Number of multiprocessors:                     %d\n",
           deviceProp.multiProcessorCount);
    printf("  Number of cores:                               %d\n",
           nGpuArchCoresPerSM[deviceProp.major] *
           deviceProp.multiProcessorCount);

    printf("  Total amount of constant memory:               %u bytes\n",
           (unsigned int)deviceProp.totalConstMem);
    printf("  Total amount of shared memory per block:       %u bytes\n",
           (unsigned int)deviceProp.sharedMemPerBlock);
    printf("  Total number of registers available per block: %d\n",
           deviceProp.regsPerBlock);
    printf("  Warp size:                                     %d\n",
           deviceProp.warpSize);
    printf("  Maximum number of threads per block:           %d\n",
           deviceProp.maxThreadsPerBlock);
    printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
           deviceProp.maxThreadsDim[0],
           deviceProp.maxThreadsDim[1],
           deviceProp.maxThreadsDim[2]);
    printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
           deviceProp.maxGridSize[0],
           deviceProp.maxGridSize[1],
           deviceProp.maxGridSize[2]);
    printf("  Maximum memory pitch:                          %u bytes\n",
           (unsigned int)deviceProp.memPitch);
    printf("  Texture alignment:                             %u bytes\n",
           (unsigned int)deviceProp.textureAlignment);
    printf("  Clock rate:                                    %.2f GHz\n",
           deviceProp.clockRate * 1e-6f);

    printf("  Concurrent copy and execution:                 %s\n",
           deviceProp.deviceOverlap ? "Yes" : "No");

    printf("  Run time limit on kernels:                     %s\n",
           deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
    printf("  Integrated:                                    %s\n",
           deviceProp.integrated ? "Yes" : "No");
    printf("  Support host page-locked memory mapping:       %s\n",
           deviceProp.canMapHostMemory ? "Yes" : "No");
    printf("  Compute mode:                                  %s\n",
           deviceProp.computeMode == cudaComputeModeDefault ?
           "Default (multiple host threads can use this device "
           "simultaneously)" :
           deviceProp.computeMode == cudaComputeModeExclusive ?
           "Exclusive (only one host thread at a time can use this device)" :
           deviceProp.computeMode == cudaComputeModeProhibited ?
           "Prohibited (no host thread can use this device)" :
           "Unknown");
  }

  return (0);
}
#endif /* FS_CUDA */
