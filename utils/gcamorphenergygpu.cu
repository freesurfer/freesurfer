/**
 * @file  gcamorphenergygpu.cu
 * @brief Holds routines to compute GCAmorph energies on the GPU
 *
 * This file holds a variety of routines which compute energies for
 * mri_ca_register
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.44 $
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
 */

#ifdef GCAMORPH_ON_GPU

#include <memory>

#include <thrust/device_new_allocator.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

#include "macros.h"
#include "cma.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"

#include "gcamorphenergy.hpp"

// Stolen from gcamorph.c
#define MIN_STD 2
#define MIN_VAR (MIN_STD*MIN_STD)


//! Texture reference for an unsigned char mri
texture<unsigned char, 3, cudaReadModeNormalizedFloat> dt_mri_uchar;

//! Texture for vx in SmoothnessEnergy
texture<float, 3, cudaReadModeElementType> dt_smooth_vx;
//! Texture for vy in SmoothnessEnergy
texture<float, 3, cudaReadModeElementType> dt_smooth_vy;
//! Texture for vz in SmoothnessEnergy
texture<float, 3, cudaReadModeElementType> dt_smooth_vz;
//! Texure for 'invalid' in SmoothnessEnergy
texture<char, 3, cudaReadModeElementType> dt_smooth_invalid;


// ==============================================================

namespace GPU
{
namespace Algorithms
{


const unsigned int kGCAmorphLLEkernelSize = 16;
const unsigned int kGCAmorphJacobEnergyKernelSize = 16;
const unsigned int kGCAmorphLabelEnergyKernelSize = 16;
const unsigned int kGCAmorphSmoothEnergyKernelSize = 16;

//! Templated texture fetch
template<typename T>
__device__ float FetchMRIVoxel( const float3 r )
{
  /*!
  Does look ups into the textures.
  Recall that the textures are configured to return
  normalised floats (except for the float texture!)
  in order to enable linear filtering.
  This in turn requires a conversion from the normalised
  value back to the true value.
  Unspecialised version writes junk
  */
  return( 100000 );
}

template<>
__device__ float FetchMRIVoxel<unsigned char>( const float3 r )
{

  float texVal;
  texVal = tex3D( dt_mri_uchar, r.x+0.5f, r.y+0.5f, r.z+0.5f );
  texVal *= UCHAR_MAX;

  return( texVal );
}


__device__
float GCAmorphDist( const float mean, const float variance,
                    const float val )
{
  float v;
  v = val - mean;
  v = v*v;
  v /= variance;

  return( v );
}

__device__
bool IsUnknown( const int label )
{
  bool res;

  res = (label==Unknown);
  res = res || (label==255);
  res = res || (label==Bright_Unknown);
  res = res || (label==Dark_Unknown);

  return(res);
}


//! Kernel to compute whether a voxel should be included in LLE calculation
__global__
void ComputeGoodLLE( const GPU::Classes::VolumeArgGPU<char> invalid,
                     const GPU::Classes::VolumeArgGPU<int> label,
                     const GPU::Classes::VolumeArgGPU<int> status,
                     GPU::Classes::VolumeArgGPU<char> good )
{

  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;
  /*!
  Computes whether each voxel should be included in the LLE
  calculation.
  The 'border' loops should really use shared memory, but
  right now, it's pointless worrying about this.
  Getting data onto the GPU takes two orders of magnitude
  more time
  */

  // Loop over z slices
  for( unsigned int iz = 0; iz< good.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( good.InVolume(ix,iy,iz) )
    {

      // Start by flagging as 'bad'
      good(ix,iy,iz) = 0;

      // Is it valid?
      if( invalid(ix,iy,iz) == GCAM_POSITION_INVALID )
      {
        // If not, go to next z slice
        continue;
      }

      // What's the status?
      if( status(ix,iy,iz) &
          (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD) )
      {
        // Go to next z slice
        continue;
      }

      // Don't use unknowns unless they border known
      if( IS_UNKNOWN(label(ix,iy,iz)) )
      {
        unsigned int diffLabels = 0;
        const int myLabel = label(ix,iy,iz);

        for( unsigned int k=max(0,iz-1);
             k<=min(invalid.dims.z-1,iz+1);
             k++ )
        {
          for( unsigned int j=max(0,iy-1);
               j<=min(invalid.dims.y-1,iy+1);
               j++ )
          {
            for( unsigned int i=max(0,ix-1);
                 i<=min(invalid.dims.x-1,ix+1);
                 i++ )
            {
              if( label(i,j,k) != myLabel )
              {
                diffLabels++;
              }
            }
          }
        }

        if( diffLabels == 0 )
        {
          // Go to next z slice
          continue;
        }
      }

      // If we get to here, it's OK
      good(ix,iy,iz) = 1;
    }
  }
}



template<typename T>
__global__
void ComputeLLE( const GPU::Classes::VolumeArgGPU<float> rx,
                 const GPU::Classes::VolumeArgGPU<float> ry,
                 const GPU::Classes::VolumeArgGPU<float> rz,
                 const GPU::Classes::VolumeArgGPU<char> good,
                 const GPU::Classes::VolumeArgGPU<float> mean,
                 const GPU::Classes::VolumeArgGPU<float> variance,
                 float* energies )
{

  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  float myEnergy;


  // Loop over z slices
  for( unsigned int iz = 0; iz < rx.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( rx.InVolume(ix,iy,iz) )
    {

      const unsigned int iLoc = rx.Index1D( ix, iy, iz );

      // See if we want to do this pixel
      if( good(ix,iy,iz) == 0 )
      {
        energies[iLoc] = 0;
        continue;
      }

      float3 r = make_float3( rx( ix, iy, iz ),
                              ry( ix, iy, iz ),
                              rz( ix, iy, iz ) );

      // Get the MRI value
      // For some reason, bounds checking causes lots of zeros
      float mriVal = FetchMRIVoxel<T>( r );

      // Compute contribution to the energy
      if( variance(ix,iy,iz) >= 0 )
      {
        // We have a valid variance
        myEnergy = GCAmorphDist( mean(ix,iy,iz),
                                 variance(ix,iy,iz),
                                 mriVal );
        myEnergy += log( variance(ix,iy,iz) );
      }
      else
      {
        myEnergy = mriVal*mriVal / MIN_VAR;
      }

      energies[iLoc] = myEnergy;
    }
  }
}

// --------------------------------------------------

//! Implements FZERO macro
__device__ bool NearZero( const float f )
{

  bool res = false;
  if( fabsf(f) < 0.0000001f )
  {
    res = true;
  }

  return( res );
}


//! Calculates delta based on exponent
__device__ float JacobDelta( const float exponent )
{
  /*!
  Calculates the value of delta for ComputeJacobEnergy.
  This states
  \f[
  \delta = \log \left ( 1 + e^{x} \right )
  \f]
  Uses Taylor expansions to avoid hideous accuracy losses.
  Note \f$2^{12} \approx e^{8.3}\f$ is used to assess
  changeover points for single precision.
  */
  float delta;

  const float kExpLimit = 8;
  const float kMaxExp = 200;

  if( exponent > kMaxExp )
  {
    delta = 0;
  }
  else if( exponent > kExpLimit )
  {
    // Can drop the 1 in the logarithm
    delta = exponent;
  }
  else if ( exponent < -kExpLimit )
  {
    // Taylor expand ln( 1 + y )
    float y = exp( exponent );

    delta = y * ( 1 + ( y * ( -0.5f ) ) );
  }
  else
  {
    // Do full calculation
    delta = log( 1 + exp( exponent ) );
  }

  return( delta );
}

//! Kernel to implement loops of gcamComputeJacobianEnergy
__global__
void ComputeJacobEnergy( const GPU::Classes::VolumeArgGPU<char> invalid,
                         const GPU::Classes::VolumeArgGPU<float> origArea1,
                         const GPU::Classes::VolumeArgGPU<float> origArea2,
                         const GPU::Classes::VolumeArgGPU<float> area1,
                         const GPU::Classes::VolumeArgGPU<float> area2,
                         const float exp_k,
                         const float thick,
                         float *energies )
{

  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  float ratio, exponent, delta;

  // Loop over z slices
  for( unsigned int iz = 0; iz < invalid.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( invalid.InVolume(ix,iy,iz) )
    {

      const unsigned int iLoc = invalid.Index1D( ix, iy, iz );

      // Check for validity
      if( invalid(ix,iy,iz) != GCAM_VALID )
      {
        energies[iLoc] = 0;
        continue;
      }

      float myEnergy = 0;

      // Look at area1
      if( !NearZero( origArea1(ix,iy,iz) ) )
      {
        ratio = area1(ix,iy,iz) / origArea1(ix,iy,iz);
        exponent = -exp_k * ratio;

        delta = JacobDelta( exponent );

        myEnergy += ( delta * thick );
      }

      if( !NearZero( origArea2(ix,iy,iz) ) )
      {
        ratio = area2(ix,iy,iz) / origArea2(ix,iy,iz);
        exponent = -exp_k * ratio;

        delta = JacobDelta( exponent );

        myEnergy += ( delta * thick );
      }

      energies[iLoc] = myEnergy;

    }

  }
}


// --------------------------------------------------

//! Kernel to mirror gcamLabelEnergy
__global__
void ComputeLabelEnergy( const GPU::Classes::VolumeArgGPU<char> invalid,
                         const GPU::Classes::VolumeArgGPU<int> status,
                         const GPU::Classes::VolumeArgGPU<float> labelDist,
                         float *energies )
{

  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  for( unsigned int iz = 0; iz < invalid.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( invalid.InVolume(ix,iy,iz) )
    {

      const unsigned int iLoc = invalid.Index1D( ix, iy, iz );

      if( (invalid(ix,iy,iz) == GCAM_POSITION_INVALID) ||
          ( (status(ix,iy,iz) & GCAM_LABEL_NODE) == 0 ) )
      {
        energies[iLoc] = 0;
      }
      else
      {
        energies[iLoc] = fabsf( labelDist(ix,iy,iz) );
      }
    }
  }
}

// --------------------------------------------------

//! Kernel to perform \f$c = a - b\f$ on volumes
template<typename T>
__global__
void SubtractVolumes( const GPU::Classes::VolumeArgGPU<T> a,
                      const GPU::Classes::VolumeArgGPU<T> b,
                      GPU::Classes::VolumeArgGPU<T> c )
{

  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  for( unsigned int iz = 0; iz < a.dims.z; iz++ )
  {
    if( a.InVolume(ix,iy,iz) )
    {

      c(ix,iy,iz) = a(ix,iy,iz) - b(ix,iy,iz);
    }
  }

}


__device__ float Fetchvx( const int ix,
                          const int iy,
                          const int iz )
{
  return( tex3D( dt_smooth_vx, ix+0.5f, iy+0.5f, iz+0.5f ) );
}

__device__ float Fetchvy( const int ix,
                          const int iy,
                          const int iz )
{
  return( tex3D( dt_smooth_vy, ix+0.5f, iy+0.5f, iz+0.5f ) );
}

__device__ float Fetchvz( const int ix,
                          const int iy,
                          const int iz )
{
  return( tex3D( dt_smooth_vz, ix+0.5f, iy+0.5f, iz+0.5f ) );
}


//! Kernel to compute the smoothness energy
__global__
void SmoothnessKernel( const GPU::Classes::VolumeArgGPU<char> invalid,
                       float *energies )
{
  /*!
  The main computation of gcamSmoothnessEnergy.
  Note that this kernel accesses 'invalid' through both
  the argument and a texture.
  The argument is there so we can use the 'InVolume' method.
  The texture is used for boundary conditions which will
  exactly match gcamSmoothnessEnergy.
  It would be preferable to change these BC, and dump
  the texture.
  However, that means changing the CPU code.
  */

  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  for( unsigned int iz = 0; iz < invalid.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( invalid.InVolume(ix,iy,iz) )
    {

      float node_sse = 0;
      unsigned int num = 0;

      // Only calculate if we're not invalid
      if( invalid(ix,iy,iz) != GCAM_POSITION_INVALID )
      {


        // Re-order loop nest in gcamSmoothnessEnergy
        for( int zk=-1; zk<=1; zk++ )
        {
          int zn = iz + zk;

          for( int yk=-1; yk<=1; yk++ )
          {
            int yn = iy + yk;

            for( int xk=-1; xk<=1; xk++ )
            {
              int xn = ix + xk;

              // Don't include self
              if( (!xk) && (!yk) && (!zk) )
              {
                continue;
              }

              // Don't include invalid neighbours

              // Use texture to get boundary conditions
              const char myInvalid = tex3D( dt_smooth_invalid,
                                            xn+0.5f, yn+0.5f, zn+0.5f );
              if( myInvalid == GCAM_POSITION_INVALID )
              {
                continue;
              }


              float dx = Fetchvx(xn,yn,zn) - Fetchvx(ix,iy,iz);
              float dy = Fetchvy(xn,yn,zn) - Fetchvy(ix,iy,iz);
              float dz = Fetchvz(xn,yn,zn) - Fetchvz(ix,iy,iz);

              node_sse += (dx*dx) + (dy*dy) + (dz*dz);
              num++;
            }
          }
        }

        // Take average
        if( num > 0 )
        {
          node_sse /= num;
        }
      }

      // Assign value
      const unsigned int iLoc = invalid.Index1D( ix, iy, iz );
      energies[iLoc] = node_sse;
    }
  }

}


// ##############################################################

void GCAmorphEnergy::ShowTimings( void )
{
#ifdef CUDA_SHOW_TIMINGS
  std::cout << "==================================" << std::endl;
  std::cout << "GPU GCAmorph energy timers" << std::endl;
  std::cout << "--------------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
  std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
  std::cout << "Timings may not be accurate" << std::endl;
#endif
  std::cout << std::endl;

  std::cout << "LLEdispatch:" << std::endl;
  std::cout << "Total         : " << GCAmorphEnergy::tLLEdispatch
            << std::endl;
  std::cout << std::endl;

  std::cout << "LogLikelihoodEnergy:" << std::endl;
  std::cout << "    Find good : " << GCAmorphEnergy::tLLEgood
            << std::endl;
  std::cout << "      Compute : " << GCAmorphEnergy::tLLEcompute
            << std::endl;
  std::cout << "Total         : " << GCAmorphEnergy::tLLEtot
            << std::endl;
  std::cout << std::endl;

  std::cout << "JacobianEnergy:" << std::endl;
  std::cout << "      Compute : " << GCAmorphEnergy::tJacobCompute
            << std::endl;
  std::cout << "Total         : " << GCAmorphEnergy::tJacobTot
            << std::endl;
  std::cout << std::endl;

  std::cout << "LabelEnergy:" << std::endl;
  std::cout << "      Compute : " << GCAmorphEnergy::tLabelCompute
            << std::endl;
  std::cout << "Total         : " << GCAmorphEnergy::tLabelTot
            << std::endl;
  std::cout << std::endl;

  std::cout << "SmoothnessEnergy:" << std::endl;
  std::cout << "  Subtraction : " << GCAmorphEnergy::tSmoothSubtract
            << std::endl;
  std::cout << "      Compute : " << GCAmorphEnergy::tSmoothCompute
            << std::endl;
  std::cout << "Total         : " << GCAmorphEnergy::tSmoothTot
            << std::endl;
  std::cout << std::endl;

  std::cout << "==================================" << std::endl;
#endif
}

// --------------------------------------------------


template<typename T>
float GCAmorphEnergy::LogLikelihoodEnergy( const GPU::Classes::GCAmorphGPU& gcam,
    const GPU::Classes::MRIframeGPU<T>& mri ) const
{
  /*!
  This the the host side function for
  gcamLogLikelihoodEnergy on the GPU.
  Note that a GCAmorphGPU implicitly only has one input for
  each location.
  This means that each covariance is just a variance,
  and negative values flag
  */

  GCAmorphEnergy::tLLEtot.Start();


  // Make sure the GCAM is sane
  gcam.CheckIntegrity();

  const dim3 gcamDims = gcam.d_rx.GetDims();
  const unsigned int nVoxels = gcamDims.x * gcamDims.y * gcamDims.z;
  // Create a 'flag' array
  GPU::Classes::VolumeGPU<char> d_good;
  d_good.Allocate( gcamDims );

  // Allocate thrust arrays
  thrust::device_ptr<float> d_energies;
  d_energies = thrust::device_new<float>( nVoxels );

  // Get the MRI into a texture
  std::auto_ptr<GPU::Classes::CTfactory> mri_texture( this->BindMRI( mri ) );


  // Run the computation
  dim3 grid, threads;
  threads.x = threads.y = kGCAmorphLLEkernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphLLEkernelSize );
  grid.z = 1;

  GCAmorphEnergy::tLLEgood.Start();
  ComputeGoodLLE<<<grid,threads,0,this->stream>>>( gcam.d_invalid,
      gcam.d_label,
      gcam.d_status,
      d_good );
  CUDA_CHECK_ERROR( "ComputeGood kernel failed!\n" );
  GCAmorphEnergy::tLLEgood.Stop();


  GCAmorphEnergy::tLLEcompute.Start();
  ComputeLLE<T><<<grid,threads,0,this->stream>>>
  ( gcam.d_rx, gcam.d_ry, gcam.d_rz,
    d_good,
    gcam.d_mean, gcam.d_variance,
    thrust::raw_pointer_cast( d_energies ) );
  CUDA_CHECK_ERROR( "ComputeLLE kernel failed!\n" );
  GCAmorphEnergy::tLLEcompute.Stop();

  CUDA_SAFE_CALL( cudaStreamSynchronize( this->stream ) );

  // Get the sum of the energies
  float energy = thrust::reduce( d_energies, d_energies+nVoxels );

  // Release thrust arrays
  thrust::device_delete( d_energies );

  GCAmorphEnergy::tLLEtot.Stop();

  return( energy );
}



template<typename T>
float GCAmorphEnergy::LLEdispatch( const GCA_MORPH *gcam,
                                   const MRI* mri ) const
{

  float energy;

  GCAmorphEnergy::tLLEdispatch.Start();

  GPU::Classes::GCAmorphGPU myGCAM;
  myGCAM.SendAll( gcam );

  GPU::Classes::MRIframeGPU<T> myMRI;
  myMRI.Allocate( mri );
  myMRI.Send( mri, 0 );

  energy = this->LogLikelihoodEnergy( myGCAM, myMRI );

  GCAmorphEnergy::tLLEdispatch.Stop();

  return( energy );

}

// --------------------------------------------------


float GCAmorphEnergy::ComputeJacobianEnergy( const GPU::Classes::GCAmorphGPU& gcam,
    const float mriThick ) const
{

  GCAmorphEnergy::tJacobTot.Start();

  // Make sure the GCAM is sane
  gcam.CheckIntegrity();

  const dim3 gcamDims = gcam.d_rx.GetDims();
  const unsigned int nVoxels = gcamDims.x * gcamDims.y * gcamDims.z;

  // Allocate thrust arrays
  thrust::device_ptr<float> d_energies;
  d_energies = thrust::device_new<float>( nVoxels );


  // Run the computation
  dim3 grid, threads;
  threads.x = threads.y = kGCAmorphJacobEnergyKernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphJacobEnergyKernelSize );
  grid.z = 1;

  GCAmorphEnergy::tJacobCompute.Start();
  ComputeJacobEnergy<<<grid,threads>>>
  ( gcam.d_invalid,
    gcam.d_origArea1, gcam.d_origArea2,
    gcam.d_area1, gcam.d_area2,
    gcam.exp_k, mriThick,
    thrust::raw_pointer_cast( d_energies ) );
  CUDA_CHECK_ERROR( "ComputeJacobEnergy kernel failed!\n" );
  GCAmorphEnergy::tJacobCompute.Stop();

  // Get the sum of the energies
  float jEnergy = thrust::reduce( d_energies, d_energies+nVoxels );

  // Release thrust arrays
  thrust::device_delete( d_energies );

  GCAmorphEnergy::tJacobTot.Stop();

  return( jEnergy );
}

// --------------------------------------------------

//! Implementation of gcamLabelEnergy for the GPU
float GCAmorphEnergy::LabelEnergy( const GPU::Classes::GCAmorphGPU& gcam ) const
{

  GCAmorphEnergy::tLabelTot.Start();

  // Make sure GCAM is sane
  gcam.CheckIntegrity();

  const dim3 gcamDims = gcam.d_rx.GetDims();
  const unsigned int nVoxels = gcamDims.x * gcamDims.y * gcamDims.z;

  // Allocate thrust arrays
  thrust::device_ptr<float> d_energies;
  d_energies = thrust::device_new<float>( nVoxels );

  // Run the computation
  dim3 grid, threads;
  threads.x = threads.y = kGCAmorphLabelEnergyKernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphLabelEnergyKernelSize );
  grid.z = 1;

  GCAmorphEnergy::tLabelCompute.Start();
  ComputeLabelEnergy<<<grid,threads>>>
  ( gcam.d_invalid,
    gcam.d_status,
    gcam.d_labelDist,
    thrust::raw_pointer_cast( d_energies ) );
  CUDA_CHECK_ERROR( "ComputeLabelEnergy kernel failed!\n" );
  GCAmorphEnergy::tLabelCompute.Stop();

  // Get the sum
  float lEnergy = thrust::reduce( d_energies, d_energies+nVoxels );

  // Release thrust arrays
  thrust::device_delete( d_energies );

  GCAmorphEnergy::tLabelTot.Stop();

  return( lEnergy );
}


// --------------------------------------------------

//! Implementation of gcamSmoothnessEnergy for the GPU
float GCAmorphEnergy::SmoothnessEnergy( GPU::Classes::GCAmorphGPU& gcam ) const
{
  /*!
  Computes the smoothness energy of the given gcam.
  Mirrors the gcamSmoothnessEnergy routine in gcamorph.c.
  The only reason gcam is not declared 'const' is because
  we have to get the 'invalid' field into a texture,
  and hence needs its cudaArray.
  We use a texture for this field to get easy boundary
  conditions - although ideally the BC would be changed.
  However, that means changing the CPU code too.
  */

  GCAmorphEnergy::tSmoothTot.Start();

  // Make sure GCAM is sane
  gcam.CheckIntegrity();

  const dim3 gcamDims = gcam.d_rx.GetDims();
  const unsigned int nVoxels = gcamDims.x * gcamDims.y * gcamDims.z;

  // Allocate thrust arrays
  thrust::device_ptr<float> d_energies;
  d_energies = thrust::device_new<float>( nVoxels );

  // Compute vx, vy, and vz (see gcamSmoothnessEnergy)
  GPU::Classes::VolumeGPU<float> vx, vy, vz;

  vx.Allocate( gcamDims );
  vy.Allocate( gcamDims );
  vz.Allocate( gcamDims );

  dim3 grid, threads;
  threads.x = threads.y = kGCAmorphSmoothEnergyKernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphSmoothEnergyKernelSize );
  grid.z = 1;

  GCAmorphEnergy::tSmoothSubtract.Start();
  SubtractVolumes<float>
  <<<grid,threads>>>
  ( gcam.d_rx, gcam.d_origx, vx );
  CUDA_CHECK_ERROR( "SubtractVolumes failed for x!" );

  SubtractVolumes<float>
  <<<grid,threads>>>
  ( gcam.d_ry, gcam.d_origy, vy );
  CUDA_CHECK_ERROR( "SubtractVolumes failed for y!" );

  SubtractVolumes<float>
  <<<grid,threads>>>
  ( gcam.d_rz, gcam.d_origz, vz );
  CUDA_CHECK_ERROR( "SubtractVolumes failed for z!" );
  GCAmorphEnergy::tSmoothSubtract.Stop();

  // Get vx, vy and vz into CUDA arrays
  GPU::Classes::CTfactory vxArray( vx, dt_smooth_vx );
  GPU::Classes::CTfactory vyArray( vy, dt_smooth_vy );
  GPU::Classes::CTfactory vzArray( vz, dt_smooth_vz );

  GPU::Classes::CTfactory invalidArray( gcam.d_invalid,
                                        dt_smooth_invalid );

  // Run the main kernel
  GCAmorphEnergy::tSmoothCompute.Start();
  SmoothnessKernel<<<grid,threads>>>
  ( gcam.d_invalid,
    thrust::raw_pointer_cast( d_energies ) );
  CUDA_CHECK_ERROR( "SmoothnessKernel failed!" );
  GCAmorphEnergy::tSmoothCompute.Stop();

  // Get the total
  float smoothEnergy = thrust::reduce( d_energies, d_energies+nVoxels );

  // Release thrust arrays
  thrust::device_delete( d_energies );


  GCAmorphEnergy::tSmoothTot.Stop();

  return( smoothEnergy );
}


// --------------------------------------------------

//! Implementation of gcamComputeSSE for the GPU
template<typename T>
float GCAmorphEnergy::ComputeSSE( GPU::Classes::GCAmorphGPU& gcam,
                                  const GPU::Classes::MRIframeGPU<T>& mri,
                                  GCA_MORPH_PARMS *parms ) const
{
  /*!
  Calls all the necessary things on the GPU to match
  gcamComputeSSE.
  For now, just a bunch of tests, to see if we have
  all the required routines in place
  */

  float l_sse, j_sse, s_sse;
  float label_sse;

  l_sse = j_sse = s_sse = 0;
  label_sse = 0;

  if( !DZERO(parms->l_area_intensity) )
  {
    std::cerr << "gcamCreateNodeLookupTable not on GPU yet!"
              << std::endl;
    exit( EXIT_FAILURE );
  }

  // Compute the metric properties
  int invalid;
  gcam.ComputeMetricProperties( invalid );


  // Get the log likelihood energy
  if( !DZERO(parms->l_log_likelihood) || !DZERO(parms->l_likelihood) )
  {
    l_sse = MAX(parms->l_log_likelihood, parms->l_likelihood) *
            this->LogLikelihoodEnergy( gcam, mri );
  }

  if( !DZERO(parms->l_multiscale) )
  {
    std::cerr << "gcamMultiscaleEnergy not on GPU yet!"
              << std::endl;
    exit( EXIT_FAILURE );
  }

  if( !DZERO(parms->l_dtrans) )
  {
    std::cerr << "gcamDistanceTransformEnergy not on GPU yet!"
              << std::endl;
    exit( EXIT_FAILURE );
  }


  if( !DZERO(parms->l_label) )
  {
    label_sse = this->LabelEnergy( gcam );
  }

  if( !DZERO(parms->l_binary) )
  {
    std::cerr << "gcamBinaryEnergy not on GPU yet!" << std::endl;
    exit( EXIT_FAILURE );
  }

  if( !DZERO(parms->l_area_intensity) )
  {
    std::cerr << "gcamAreaIntensityEnergy not on GPU yet!"
              << std::endl;
    exit( EXIT_FAILURE );
  }

  if( !DZERO(parms->l_map) )
  {
    std::cerr << "gcamMapEnergy not on GPU yet!" << std::endl;
    exit( EXIT_FAILURE );
  }


  if( !DZERO(parms->l_expansion) )
  {
    std::cerr << "gcamExpansionEnergy not on GPU yet!" << std::endl;
    exit( EXIT_FAILURE );
  }


  if( !DZERO(parms->l_distance) )
  {
    std::cerr << "gcamDistanceEnergy not on GPU yet!" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Compute the Jacobian energy
  if( !DZERO(parms->l_jacobian) )
  {
    j_sse = parms->l_jacobian *
            this->ComputeJacobianEnergy( gcam, mri.GetThickness() );
  }

  if( !DZERO(parms->l_area) )
  {
    std::cerr << "gcamAreaEnergy not on GPU yet!" << std::endl;
    exit( EXIT_FAILURE );
  }

  if( !DZERO(parms->l_area_smoothness) )
  {
    std::cerr << "gcamAreaEnergy not on GPU yet!" << std::endl;
    exit( EXIT_FAILURE );
  }

  if( !DZERO(parms->l_smoothness) )
  {
    s_sse = parms->l_smoothness * this->SmoothnessEnergy( gcam );
  }

  // Note extra 'l'
  if( !DZERO(parms->l_lsmoothness) )
  {
    std::cerr << "gcamLSmoothnessEnergy not on GPU yet!" << std::endl;
    exit( EXIT_FAILURE );
  }

  if( !DZERO(parms->l_spring) )
  {
    std::cerr << "gcamSpringEnergy not on GPU yet!" << std::endl;
    exit( EXIT_FAILURE );
  }

  float sse;

  sse = l_sse + j_sse + s_sse;
  sse += label_sse;

  return( sse );
}



// --------------------------------------------------

#if 0
template<typename T>
float GCAmorphEnergy::ComputeRMS( GPU::Classes::GCAmorphGPU& gcam,
                                  const GPU::Classes::MRIframeGPU<T>& mri,
                                  GCA_MORPH_PARMS *parms ) const
{

  float sse = this->ComputeSSE( gcam, mri, parms );

  const dim3 dims = gcam.d_rx.GetDims();
  float nVoxels = dims.x;
  nVoxels *= dims.y;
  nVoxels *= dims.z;

  float rms = sqrtf( sse/nVoxels );

  return( rms );
}
#endif

template<typename T>
float GCAmorphEnergy::RMSdispatch( GPU::Classes::GCAmorphGPU& gcam,
                                   const MRI *mri,
                                   GCA_MORPH_PARMS *parms ) const
{


  GPU::Classes::MRIframeGPU<T> mriGPU;

  mriGPU.Allocate( mri );
  mriGPU.Send( mri, 0 );

  float rms = this->ComputeRMS( gcam, mriGPU, parms );

  return( rms );
}


// --------------------------------------------------------



template<typename T>
GPU::Classes::CTfactory* GCAmorphEnergy::BindMRI( const GPU::Classes::MRIframeGPU<T>& mri ) const
{
  std::cerr << __PRETTY_FUNCTION__
            << ": Unrecognised MRI type" << std::endl;
  abort();
}

// Texture binding specialisations

template<>
GPU::Classes::CTfactory* GCAmorphEnergy::BindMRI<unsigned char>( const GPU::Classes::MRIframeGPU<unsigned char>& mri ) const
{

  GPU::Classes::CTfactory* tmp = NULL;
  tmp = new GPU::Classes::CTfactory( mri, dt_mri_uchar,
                                     cudaFilterModeLinear );

  return( tmp );
}


// --------------------------------------------------------

// Declare static members

//! Timer for LLEdispatch
SciGPU::Utilities::Chronometer GCAmorphEnergy::tLLEdispatch;

//! Timer for log likelihood energy
SciGPU::Utilities::Chronometer GCAmorphEnergy::tLLEtot;
//! Timer for LLE 'good' assesments
SciGPU::Utilities::Chronometer GCAmorphEnergy::tLLEgood;
//! Timer for LLE calculation
SciGPU::Utilities::Chronometer GCAmorphEnergy::tLLEcompute;

//! Timer for Jacobian Energy
SciGPU::Utilities::Chronometer GCAmorphEnergy::tJacobTot;
//! Timer for Jacobian Energy kernel
SciGPU::Utilities::Chronometer GCAmorphEnergy::tJacobCompute;

//! Mutable for Label Energy
SciGPU::Utilities::Chronometer GCAmorphEnergy::tLabelTot;
//! Timer for Label Energy kernel
SciGPU::Utilities::Chronometer GCAmorphEnergy::tLabelCompute;

//! Timer for smoothness energy
SciGPU::Utilities::Chronometer GCAmorphEnergy::tSmoothTot;
//! Timer for subtractions in smoothness energy
SciGPU::Utilities::Chronometer GCAmorphEnergy::tSmoothSubtract;
//! Timer for smoothness computation itself
SciGPU::Utilities::Chronometer GCAmorphEnergy::tSmoothCompute;





}
}










static GPU::Algorithms::GCAmorphEnergy myEnergy;


//! Wrapper around GPU class for LLE
float gcamLogLikelihoodEnergyGPU( const GCA_MORPH *gcam,
                                  const MRI* mri )
{

  float energy;

  switch( mri->type )
  {

  case MRI_UCHAR:
    energy = myEnergy.LLEdispatch<unsigned char>( gcam, mri );
    break;


  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised MRI type" << std::endl;
    exit( EXIT_FAILURE );
  }

  return( energy );

}


//! Wrapper around GPU class for JacobianEnergy
float gcamJacobianEnergyGPU( const GCA_MORPH *gcam,
                             const MRI* mri )
{

  float energy;

  const float thick = ( mri ? mri->thick : 1.0 );

  GPU::Classes::GCAmorphGPU myGCAM;
  myGCAM.SendAll( gcam );

  energy = myEnergy.ComputeJacobianEnergy( myGCAM, thick );

  return( energy );
}


//! Wrapper around GPU class for the LabelEnergy
float gcamLabelEnergyGPU( const GCA_MORPH *gcam )
{

  float energy;

  GPU::Classes::GCAmorphGPU myGCAM;
  myGCAM.SendAll( gcam );

  energy = myEnergy.LabelEnergy( myGCAM );

  return( energy );
}


//! Wrapper around GPU class for the SmoothnessEnergy
float gcamSmoothnessEnergyGPU( const GCA_MORPH *gcam )
{

  float energy;

  GPU::Classes::GCAmorphGPU myGCAM;
  myGCAM.SendAll( gcam );

  energy = myEnergy.SmoothnessEnergy( myGCAM );

  return( energy );
}



//! Wrapper for GCAMcomputeRMS
float gcamComputeRMSonGPU( GCA_MORPH *gcam,
                           const MRI* mri,
                           GCA_MORPH_PARMS *parms )
{

  float rms;

  GPU::Classes::GCAmorphGPU myGCAM;
  myGCAM.SendAll( gcam );


  switch( mri->type )
  {

  case MRI_UCHAR:
    rms = myEnergy.RMSdispatch<unsigned char>( myGCAM, mri, parms );
    break;

  case MRI_FLOAT:
    rms = myEnergy.RMSdispatch<float>( myGCAM, mri, parms );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised MRI type" << std::endl;
    exit( EXIT_FAILURE );
  }


  myGCAM.RecvAll( gcam );

  return( rms );
}


#endif
