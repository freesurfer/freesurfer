/**
 * @file  mrimean_cuda.hpp
 * @brief Holds MRI mean function for the GPU (Class definition)
 *
 * Implements MRI mean function on the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.2 $
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

#ifndef MRI_MEAN_CUDA_HPP
#define MRI_MEAN_CUDA_HPP

#include "mriframegpu.hpp"

namespace GPU
{
namespace Algorithms
{


//! Wrapper class for the MRI mean algorithm
class MRImean
{
public:
  //! Constructor with stream (also default constructor)
  MRImean( const cudaStream_t s = 0 ) : stream(s),
    h_workspace(NULL),
    workSize(0) {}

  //! Destructor
  ~MRImean( void );

  // Show timings on std.out
  static void ShowTimings( void );

  // -------------------------------------

  //! Dispatch for data on the CPU of unknown type
  void DoMean( const MRI* src, MRI* dst,
               const unsigned int wSize,
               const unsigned int srcFrame = 0,
               const unsigned int dstFrame = 0 ) const;


  //! Templated dispatch for known data types
  template<typename T, typename U>
  void MeanDispatch( const MRI* src, MRI* dst,
                     const unsigned int wSize,
                     const int srcFrame, const int dstFrame ) const;

  //! Runs the mean filtering kernel
  template<typename T, typename U>
  void RunGPU( const GPU::Classes::MRIframeGPU<T> &src,
               GPU::Classes::MRIframeGPU<U> &dst,
               const unsigned int wSize ) const;

  // #####################################
private:
  //! Stream which should be used for this instance
  cudaStream_t stream;
  //! Private pinned memory workspace
  mutable char* h_workspace;
  //! Size of private workspace
  mutable size_t workSize;

  static SciGPU::Utilities::Chronometer tMem, tHostMem;
  static SciGPU::Utilities::Chronometer tSend, tRecv, tCompute;
  static SciGPU::Utilities::Chronometer tTotal;

  //! Wrapper function
  template<typename T>
  void DispatchWrap( const MRI* src, MRI* dst,
                     const unsigned int wSize,
                     const int srcFrame, const int dstFrame ) const;

  //! Ensures internal pinned memory buffer is at least of size nBytes
  void Allocate( const size_t nBytes ) const;

  //! Releases internal pinned memory buffer
  void Release( void ) const;

  // =========================
  // Prevent copying

  MRImean( const MRImean& src ) : stream(0),
    h_workspace(NULL),
    workSize(0)
  {
    std::cerr << __FUNCTION__
              << ": Please don't copy these objects"
              << std::endl;
    exit( EXIT_FAILURE );
  }

  MRImean& operator=( const MRImean& src )
  {
    std::cerr << __FUNCTION__
              << ": Please don't copy these objects"
              << std::endl;
    exit( EXIT_FAILURE );
  }


};
}
}





#endif
