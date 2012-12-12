/**
 * @file  volumegpucompare.hpp
 * @brief Holds routines to compare VolumeGPU types
 *
 * Holds routines (inside a class) to perform comparisons
 * between instances of VolumeGPUs.
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.6 $
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

#include <cstdlib>

#include <iostream>

#include <thrust/device_new_allocator.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <thrust/reduce.h>


#include "cudacheck.h"
#include "cudatypeutils.hpp"

#include "volumegpu.hpp"

#ifndef VOLUME_GPU_COMPARE_CUDA_H
#define VOLUME_GPU_COMPARE_CUDA_H



namespace GPU
{
namespace Algorithms
{

const unsigned int kAbsDiffBlockSize = 16;
const unsigned int kErrL2BlockSize = 16;


//! Kernel to compute absolute differences between all voxels
template<typename T>
__global__
void ComputeAbsDiffs( const GPU::Classes::VolumeArgGPU<T> a,
                      const GPU::Classes::VolumeArgGPU<T> b,
                      T *diffs )
{

  const unsigned int ix = threadIdx.x + ( blockIdx.x * blockDim.x );
  const unsigned int iy = threadIdx.y + ( blockIdx.y * blockDim.y );


  // Check if our voxel is in range
  if( !a.InVolume( ix, iy, 0 ) )
  {
    return;
  }

  // Loop over each z slice
  for( unsigned int iz = 0; iz < a.dims.z; iz++ )
  {
    // Compute the 1D index of the current location
    const unsigned int iLoc = a.Index1D( ix, iy, iz );

    // Get the difference
    // Convert to float to avoid 'unsigned' trouble
    float diff = static_cast<float>(a(ix,iy,iz)) - b(ix,iy,iz);
    diffs[iLoc] = fabs( diff );
  }
}


//! Kernel for the error in the L2 norm
template<typename T>
__global__
void ErrL2Compute( const GPU::Classes::VolumeArgGPU<T> cmp,
                   const GPU::Classes::VolumeArgGPU<T> ref,
                   double *errors,
                   double *refL2s )
{


  const unsigned int ix = threadIdx.x + ( blockIdx.x * blockDim.x );
  const unsigned int iy = threadIdx.y + ( blockIdx.y * blockDim.y );


  // Check if our voxel is in range
  if( !ref.InVolume( ix, iy, 0 ) )
  {
    return;
  }

  // Loop over each z slice
  for( unsigned int iz = 0; iz < ref.dims.z; iz++ )
  {
    // Compute the 1D index of the current location
    const unsigned int iLoc = ref.Index1D( ix, iy, iz );

    double tmp = cmp(ix,iy,iz) - ref(ix,iy,iz);
    errors[iLoc] = tmp*tmp;
    refL2s[iLoc] = ref(ix,iy,iz) * ref(ix,iy,iz);
  }

}

class VolumeGPUCompare
{

public:

  //! Finds voxels with maximum difference
  template<typename T>
  void MaxDiff( const GPU::Classes::VolumeGPU<T>& a,
                const GPU::Classes::VolumeGPU<T>& b,
                T& maxDiff,
                dim3& maxLoc ) const
  {
    /*!
      Routine to compare two volumes, and find
      the voxel with the maximum difference.
    */
    const dim3 aDims = a.GetDims();
    const dim3 bDims = b.GetDims();

    // Verify dimensions
    if( aDims != bDims )
    {
      std::cerr << __FUNCTION__
                << ": Dimension mismatch"
                << std::endl;
      abort();
    }


    // Allocate 'difference' array
    thrust::device_ptr<T> d_absDiffs;

    const unsigned int nVoxels = aDims.x * aDims.y * aDims.z;

    d_absDiffs = thrust::device_new<T>( nVoxels );


    // Run the difference kernel
    dim3 grid, threads;
    threads.x = threads.y = kAbsDiffBlockSize;
    threads.z = 1;

    grid = a.CoverBlocks( kAbsDiffBlockSize );
    grid.z = 1;

    ComputeAbsDiffs<T>
    <<<grid,threads>>>
    ( a, b, thrust::raw_pointer_cast( d_absDiffs ) );
    CUDA_CHECK_ERROR( "ComputeAbsDiffs kernel failed!\n" );


    // Extract the maximum and its location
    thrust::device_ptr<T> d_maxLoc;
    d_maxLoc = thrust::max_element( d_absDiffs, d_absDiffs + nVoxels );

    maxDiff = *d_maxLoc;
    maxLoc = a.Index3D( d_maxLoc - d_absDiffs );

    // Release 'difference' array
    thrust::device_delete( d_absDiffs );

  }


  //! Computes error in the L2 norm
  template<typename T>
  double ErrL2Norm( const GPU::Classes::VolumeGPU<T>& cmp,
                    const GPU::Classes::VolumeGPU<T>& ref ) const
  {
    /*!
      Routine to compute the error in the L2 norm between
      two volumes
    */
    const dim3 cmpDims = cmp.GetDims();
    const dim3 refDims = ref.GetDims();

    // Verify dimensions
    if( refDims != cmpDims )
    {
      std::cerr << __FUNCTION__
                << ": Dimension mismatch"
                << std::endl;
      abort();;
    }

    // Compute number of voxels
    const unsigned int nVoxels = refDims.x * refDims.y * refDims.z;

    // Allocate thrust arrays
    thrust::device_ptr<double> d_err;
    thrust::device_ptr<double> d_reference;
    d_err = thrust::device_new<double>( nVoxels );
    d_reference = thrust::device_new<double>( nVoxels );

    // Run the kernel
    dim3 grid, threads;
    threads.x = threads.y = kErrL2BlockSize;
    threads.z = 1;

    grid = ref.CoverBlocks( kErrL2BlockSize );
    grid.z = 1;

    ErrL2Compute<T><<<grid,threads>>>
    ( cmp, ref,
      thrust::raw_pointer_cast( d_err ),
      thrust::raw_pointer_cast( d_reference ) );
    CUDA_CHECK_ERROR( "ErrL2Compute kernel failed!\n" );

    // Extract sums
    double totErr = thrust::reduce( d_err, d_err+nVoxels );
    double totRef = thrust::reduce( d_reference, d_reference+nVoxels );

    // Release thrust arrays
    thrust::device_delete( d_err );
    thrust::device_delete( d_reference );

    return( sqrt( totErr / totRef ) );

  }


};



}
}



#endif
