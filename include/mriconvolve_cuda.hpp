/**
 * @file  mriconvolve_cuda.hpp
 * @brief Holds MRI convolution functions for the GPU.
 *
 * Implements MRI convolutions on the GPU. These routines will be hooked
 * into the higher level routines in the main library.
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.4 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef MRI_CONVOLVE_CUDA_HPP
#define MRI_CONVOLVE_CUDA_HPP

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>


#include "chronometer.hpp"

#include "mriframegpu.hpp"

namespace GPU
{

namespace Algorithms
{

//! Class to contain MRI convolution algorithms
/*!
  This is a class designed to encapsulate all of
  the GPU convolution algorithms.
  It maintains internal timers for various
  operations, and when the destructor is called
  it can print these to stdout.
  A static instance of this class is called
  by the C routines in the Freesurfer suite.
*/
class MRIconvolve
{
public:
  // ---------------------------------------
  //! Constructor with stream (also default)
  MRIconvolve( const cudaStream_t s = 0 ) : stream(s),
    d_kernel(NULL),
    kernelAllocSize(0),
    h_workspace(NULL),
    workSize(0) {}

  //! Destructor
  ~MRIconvolve( void );

  //! Routine to print timers to std::cout
  static void ShowTimings( void );

  // ---------------------------------------

  //! Dispatch for 1D convolution of unknown type
  void Convolve1D( const MRI* src, MRI* dst,
                   const int srcFrame, const int dstFrame,
                   const float *kernel,
                   const unsigned int kernelLength,
                   const int axis );

  //! Dispatch routine with transfers
  template<typename T, typename U>
  void Dispatch1D( const MRI* src, MRI* dst,
                   const int srcFrame, const int dstFrame,
                   const int axis ) const;

  //! Runs the 1D convolution kernel on the GPU
  template<typename T, typename U>
  void RunGPU1D( const GPU::Classes::MRIframeGPU<T> &src,
                 GPU::Classes::MRIframeGPU<U> &dst,
                 const int axis ) const;


  // ---------------------------------------

  //! Dispatch for 3D convolution of unknown MRI with single kernel
  void Convolve3D( const MRI* src, MRI* dst,
                   const int srcFrame, const int dstFrame,
                   const float *kernel,
                   const unsigned int kernelLength );

  //! Wrapper for 3D convolution with single kernel
  template<typename T>
  void Conv3DSingleKernelDispatch( const MRI* src,
                                   MRI* dst,
                                   const int srcFrame,
                                   const int dstFrame ) const;

  //! Dispatch for 3D convolution of unknown MRI with separate kernels
  void Convolve3D( const MRI* src,
                   MRI* dst,
                   const int srcFrame,
                   const int dstFrame,
                   const float* xKernel,
                   const unsigned int xL,
                   const float* yKernel,
                   const unsigned int yL,
                   const float* zKernel,
                   const unsigned int zL );

  //! Wrapper for a 3D convolution with different kernels
  template<typename T>
  void Convolve3DMultiKernelDispatch( const MRI* src,
                                      MRI* dst,
                                      const int srcFrame,
                                      const int dstFrame,
                                      const float* xKernel,
                                      const unsigned int xL,
                                      const float* yKernel,
                                      const unsigned int yL,
                                      const float* zKernel,
                                      const unsigned int zL );


  // ---------------------------------------

  //! Gets the given convolution kernel onto the GPU.
  void BindKernel( const float *kernel, const unsigned int nVals );

  //! Unbinds the kernel texture
  void UnbindKernel( void );


  // ==============================================================
private:
  //! Stream to use for operations
  cudaStream_t stream;

  //! Device memory to hold convolution kernel
  float *d_kernel;
  //! Currently allocated convolution kernel size in floats (may not be fully used)
  unsigned int kernelAllocSize;

  //! Private pinned memory workspace
  mutable char* h_workspace;
  //! Size of private workspace
  mutable size_t workSize;

  static SciGPU::Utilities::Chronometer tMem;
  static SciGPU::Utilities::Chronometer tHostMem;
  static SciGPU::Utilities::Chronometer tSend, tRecv;
  static SciGPU::Utilities::Chronometer tTextureKernel;
  static SciGPU::Utilities::Chronometer tCompute1d, tTotal1d;
  static SciGPU::Utilities::Chronometer tCompute3d, tTotal3d;

  // ----------------------------------------------
  // Suppress copying
  MRIconvolve( const MRIconvolve& src ) : stream(0),
    d_kernel(NULL),
    kernelAllocSize(0),
    h_workspace(NULL),
    workSize(0)
  {
    std::cerr << __FUNCTION__
              << ": Please do not copy"
              << std::endl;
    abort();
  }

  MRIconvolve& operator=( const MRIconvolve& src )
  {
    std::cerr << __FUNCTION__
              << ": Please do not copy"
              << std::endl;
    abort();
  }


  // ----------------------------------------------

  //! Allocates convolution kernel array on the device
  void AllocateKernel( const unsigned int nFloats );

  //! Releases device array associated with convolution kernel
  void ReleaseKernel( void );

  //! Ensures internal pinned memory buffer is at least of size nBytes
  void AllocateWorkspace( const size_t nBytes ) const;

  //! Releases internal pinned memory buffer
  void ReleaseWorkspace( void ) const;


  //! Dispatch wrapper for 1D convolutions based on dst type
  template<typename T>
  void DispatchWrap1D( const MRI* src, MRI* dst,
                       const int srcFrame, const int dstFrame,
                       const int axis ) const;

};


}
}


#endif
