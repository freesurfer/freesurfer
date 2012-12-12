/**
 * @file  gcamorphtermgpu.cu
 * @brief Holds routines to compute GCAmorph terms on the GPU
 *
 * This file holds a variety of routines which compute terms for
 * mri_ca_register
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.36 $
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
#include "voxlist.h"
#include "mrinorm.h"
#include "gca.h"
#include "error.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"

#include "cudatypeutils.hpp"

#include "gcamorphtermgpu.hpp"

// Stolen from gcamorph.c
#define MIN_STD 2
#define MIN_VAR (MIN_STD*MIN_STD)


//! Texure for 'invalid' in Smoothness
texture<char, 3, cudaReadModeElementType> dt_smooth_invalid;

//! Texture reference for unsigned char mri
texture<unsigned char, 3, cudaReadModeNormalizedFloat> dt_mri_uchar;
//! Texture reference for unsigned char smoothed mri
texture<unsigned char, 3, cudaReadModeNormalizedFloat> dt_mri_smooth_uchar;


// ==============================================================


namespace GPU
{
namespace Algorithms
{

const unsigned int kGCAmorphSmoothTermKernelSize = 16;

const unsigned int kGCAmorphJacobTermKernelSize = 16;

const unsigned int kGCAmorphLLTermKernelSize = 16;

const unsigned int kGCAmorphLabelTermCopyDeltasKernelSize = 16;
const unsigned int kGCAmorphLabelTermPostAntConsistencyKernelSize = 16;
const unsigned int kGCAmorphLabelTermFinalKernelSize = 16;

// ##############################################################

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



//! Kernel to compute the smoothness term
__global__
void SmoothnessTermKernel( const GPU::Classes::VolumeArgGPU<float> vx,
                           const GPU::Classes::VolumeArgGPU<float> vy,
                           const GPU::Classes::VolumeArgGPU<float> vz,
                           GPU::Classes::VolumeArgGPU<float> dx,
                           GPU::Classes::VolumeArgGPU<float> dy,
                           GPU::Classes::VolumeArgGPU<float> dz,
                           const float l_smoothness )
{
  /*!
  Main computation of the smoothness term.
  */

  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  for( unsigned int iz = 0; iz < dx.dims.z; iz++ )
  {
    // Only compute if ix, iy & iz are inside the bounding box
    if( dx.InVolume(ix,iy,iz) )
    {

      // Only calculate if we're not invalid
      const char myValid = tex3D( dt_smooth_invalid,
                                  ix+0.5f, iy+0.5f, iz+0.5f );
      if( myValid != GCAM_POSITION_INVALID )
      {

        const float myvx = vx(ix,iy,iz);
        const float myvy = vy(ix,iy,iz);
        const float myvz = vz(ix,iy,iz);

        float ddx, ddy, ddz;
        ddx = ddy = ddz = 0;

        unsigned int num = 0;

        // Re-order loop nest in gcamSmoothnessEnergy
        for( int zk=-1; zk<=1; zk++ )
        {
          int zn = iz + zk;
          zn = min( static_cast<int>(dx.dims.z-1), zn );
          zn = max( 0, zn );

          for( int yk=-1; yk<=1; yk++ )
          {
            int yn = iy + yk;
            yn = min( static_cast<int>(dx.dims.y-1), yn );
            yn = max( 0, yn );

            for( int xk=-1; xk<=1; xk++ )
            {
              int xn = ix + xk;
              xn = min( static_cast<int>(dx.dims.x-1), xn );
              xn = max( 0, xn );

              // Don't include self
              if( (!xk) && (!yk) && (!zk) )
              {
                continue;
              }

              // Don't include invalid neighbours

              // Use texture to get boundary conditions
              const char neighbourInvalid = tex3D( dt_smooth_invalid,
                                                   xn+0.5f,
                                                   yn+0.5f,
                                                   zn+0.5f );
              if( neighbourInvalid == GCAM_POSITION_INVALID )
              {
                continue;
              }


              ddx += ( vx(xn,yn,zn) - myvx );
              ddy += ( vy(xn,yn,zn) - myvy );
              ddz += ( vz(xn,yn,zn) - myvz );

              num++;

            }
          }
        }

        if( num > 0 )
        {
          // Rescale
          ddx *= l_smoothness / num;
          ddy *= l_smoothness / num;
          ddz *= l_smoothness / num;

          // Update the GCAmorph
          dx(ix,iy,iz) += ddx;
          dy(ix,iy,iz) += ddy;
          dz(ix,iy,iz) += ddz;
        }

      } // End if( node valid )


    } // End if( inVolume )

  } // End loop over z
}



void GCAmorphTerm::Smoothness( GPU::Classes::GCAmorphGPU& gcam,
                               const float l_smoothness ) const
{
  /*!
  Computes ther smoothness term of the given gcam, following
  the routine gcamSmoothnessTerm in gcamorph.c.
  */

  if( DZERO(l_smoothness) )
  {
    return;
  }

  GCAmorphTerm::tSmoothTot.Start();

  // Make sure GCAM is sane
  gcam.CheckIntegrity();

  const dim3 gcamDims = gcam.d_rx.GetDims();

  // Compute vx, vy, and vz (see gcamSmoothnessEnergy)
  GPU::Classes::VolumeGPU<float> vx, vy, vz;

  vx.Allocate( gcamDims );
  vy.Allocate( gcamDims );
  vz.Allocate( gcamDims );

  dim3 grid, threads;
  threads.x = threads.y = kGCAmorphSmoothTermKernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphSmoothTermKernelSize );
  grid.z = 1;

  GCAmorphTerm::tSmoothSubtract.Start();
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
  GCAmorphTerm::tSmoothSubtract.Stop();


  // Also have to get the 'invalid' field to its texture
  GPU::Classes::CTfactory( gcam.d_invalid, dt_smooth_invalid );

  // Run the main kernel
  GCAmorphTerm::tSmoothCompute.Start();
  SmoothnessTermKernel<<<grid,threads>>>
  ( vx, vy, vz, gcam.d_dx, gcam.d_dy, gcam.d_dz, l_smoothness );
  CUDA_CHECK_ERROR( "SmoothnessTermKernelFailed!" );
  GCAmorphTerm::tSmoothCompute.Stop();


  GCAmorphTerm::tSmoothTot.Stop();

}

// ##############################################################

__global__
void NormKernel( const GPU::Classes::VolumeArgGPU<float> dx,
                 const GPU::Classes::VolumeArgGPU<float> dy,
                 const GPU::Classes::VolumeArgGPU<float> dz,
                 float* norms )
{

  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  float myNorm;

  // Loop over z slices
  for( unsigned int iz = 0; iz < dx.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( dx.InVolume(ix,iy,iz) )
    {

      const unsigned int iLoc = dx.Index1D( ix, iy, iz );

      myNorm = ( dx(ix,iy,iz) * dx(ix,iy,iz) )
               + ( dy(ix,iy,iz) * dy(ix,iy,iz) )
               + ( dz(ix,iy,iz) * dz(ix,iy,iz) );

      myNorm = sqrtf( myNorm );

      norms[iLoc] = myNorm;
    }
  }
}



__device__
void JacobianTermAtNode( const GPU::Classes::VolumeArgGPU<float>& x,
                         const GPU::Classes::VolumeArgGPU<float>& y,
                         const GPU::Classes::VolumeArgGPU<float>& z,
                         const GPU::Classes::VolumeArgGPU<float>& area1,
                         const GPU::Classes::VolumeArgGPU<float>& area2,
                         const GPU::Classes::VolumeArgGPU<float>& origArea1,
                         const GPU::Classes::VolumeArgGPU<float>& origArea2,
                         const GPU::Classes::VolumeArgGPU<char>& invalid,
                         const int i, const int j, const int k,
                         const float l_jacobian, const float exp_k,
                         float& pdx, float& pdy, float& pdz )
{
  /*!
  This is a direct copy of the corresponding CPU routine
  gcamJacobianTermAtNode.
  It's almost certain that lots of data could be re-used in
  shared memory
  */
  const unsigned int kAreaNeighbours = 8;
  const float kMaxExp = 200;

  const dim3 dims = x.dims;

  // Zero results
  pdx = pdy = pdz = 0;

  float3 vgrad = make_float3( 0, 0, 0 );

  int invert;
  dim3 nLoc, niLoc, njLoc, nkLoc;
  float orig_area, area;

  for( unsigned int n=0; n<kAreaNeighbours; n++ )
  {
    invert = 1;

    switch( n )
    {
    default:
    case 0: // Central node
      if( (i+1 >= dims.x) || (j+1 >= dims.y) || (k+1 >= dims.z) )
      {
        continue; // Go to next iteration of for loop
      }
      nLoc = make_uint3(i,j,k);
      niLoc = make_uint3(i+1,j,k);
      njLoc = make_uint3(i,j+1,k);
      nkLoc = make_uint3(i,j,k+1);
      break;

    case 1: // i-1 node
      if( (i == 0) || (j+1 >= dims.y) || (k+1 >= dims.z) )
      {
        continue;
      }
      nLoc = make_uint3(i-1,j,k);
      niLoc = make_uint3(i,j,k);
      njLoc = make_uint3(i-1,j+1,k);
      nkLoc = make_uint3(i-1,j,k+1);
      break;

    case 2: // j-1 node
      if( (i+1 >= dims.x) || (j == 0) || (k+1 >= dims.z) )
      {
        continue;
      }
      nLoc = make_uint3(i,j-1,k);
      niLoc = make_uint3(i+1,j-1,k);
      njLoc = make_uint3(i,j,k);
      nkLoc = make_uint3(i,j-1,k+1);
      break;

    case 3: // k-1 node
      if( (i+1 >= dims.x) || (j+1 >= dims.y) || (k == 0) )
      {
        continue;
      }
      nLoc = make_uint3(i,j,k-1);
      niLoc = make_uint3(i+1,j,k-1);
      njLoc = make_uint3(i,j+1,k-1);
      nkLoc = make_uint3(i,j,k);
      break;

    case 4: // Left-handed central node
      if( (i==0) || (j==0) || (k==0) )
      {
        continue; // Go to next iteration of for loop
      }
      invert = -1;
      nLoc = make_uint3(i,j,k);
      niLoc = make_uint3(i-1,j,k);
      njLoc = make_uint3(i,j-1,k);
      nkLoc = make_uint3(i,j,k-1);
      break;

    case 5: // i+1 node
      if( (i+1>= dims.x) || (j==0) || (k==0) )
      {
        continue;
      }
      invert = -1;
      nLoc = make_uint3(i+1,j,k);
      niLoc = make_uint3(i,j,k);
      njLoc = make_uint3(i+1,j-1,k);
      nkLoc = make_uint3(i+1,j,k-1);
      break;

    case 6: // j+1 node
      if( (i==0) || (j+1>=dims.y) || (k==0) )
      {
        continue;
      }
      invert = -1;
      nLoc = make_uint3(i,j+1,k);
      niLoc = make_uint3(i-1,j+1,k);
      njLoc = make_uint3(i,j,k);
      nkLoc = make_uint3(i,j+1,k-1);
      break;

    case 7: // k+1 node
      if( (i==0) || (j==0) || (k+1>=dims.z) )
      {
        continue;
      }
      invert = -1;
      nLoc = make_uint3(i,j,k+1);
      niLoc = make_uint3(i-1,j,k+1);
      njLoc = make_uint3(i,j-1,k+1);
      nkLoc = make_uint3(i,j,k);
      break;
    }

    if( invert > 0 )
    {
      // Right handed co-ordinate system
      orig_area = origArea1( nLoc.x, nLoc.y, nLoc.z );
      area = area1( nLoc.x, nLoc.y, nLoc.z );
    }
    else
    {
      // Left handed co-ordinate system
      orig_area = origArea2( nLoc.x, nLoc.y, nLoc.z );
      area = area2( nLoc.x, nLoc.y, nLoc.z );
    }

    // Check for area (volume...?) close to zero
    if( fabs(orig_area) < 1e-6f )
    {
      continue;
    }

    // Check for validity
    if( (invalid(nLoc.x,nLoc.y,nLoc.z) == GCAM_POSITION_INVALID) ||
        (invalid(niLoc.x,niLoc.y,niLoc.z) == GCAM_POSITION_INVALID) ||
        (invalid(njLoc.x,njLoc.y,njLoc.z) == GCAM_POSITION_INVALID) ||
        (invalid(nkLoc.x,nkLoc.y,nkLoc.z) == GCAM_POSITION_INVALID) )
    {
      continue;
    }

    // We've got somewhere
    const float3 vi = make_float3( x(niLoc.x,niLoc.y,niLoc.z) - x(nLoc.x,nLoc.y,nLoc.z),
                                   y(niLoc.x,niLoc.y,niLoc.z) - y(nLoc.x,nLoc.y,nLoc.z),
                                   z(niLoc.x,niLoc.y,niLoc.z) - z(nLoc.x,nLoc.y,nLoc.z) );

    const float3 vj = make_float3( x(njLoc.x,njLoc.y,njLoc.z) - x(nLoc.x,nLoc.y,nLoc.z),
                                   y(njLoc.x,njLoc.y,njLoc.z) - y(nLoc.x,nLoc.y,nLoc.z),
                                   z(njLoc.x,njLoc.y,njLoc.z) - z(nLoc.x,nLoc.y,nLoc.z) );

    const float3 vk = make_float3( x(nkLoc.x,nkLoc.y,nkLoc.z) - x(nLoc.x,nLoc.y,nLoc.z),
                                   y(nkLoc.x,nkLoc.y,nkLoc.z) - y(nLoc.x,nLoc.y,nLoc.z),
                                   z(nkLoc.x,nkLoc.y,nkLoc.z) - z(nLoc.x,nLoc.y,nLoc.z) );

    float ratio = area / orig_area;

    float exponent = exp_k * ratio;
    if( exponent > kMaxExp )
    {
      exponent = kMaxExp;
    }

    // Might want to consider taylor series here....
    float delta = ( invert * exp_k / orig_area ) * ( 1 / ( 1 + exp(exponent) ) );

    // Compute and accumulate cross products
    float3 vjxk, vkxi, vixj, vtmp;
    switch( n )
    {
    default:
    case 4:
    case 0: // Central node
      vjxk = cross( vj, vk );
      vkxi = cross( vk, vi );
      vixj = cross( vi, vj );
      vtmp = vixj + vjxk + vkxi;
      vtmp *= -delta;
      break;

    case 5: // i+1
    case 1: // i-1
      vjxk = cross( vj, vk );
      vtmp = vjxk * delta;
      break;

    case 6: // j+1
    case 2: // j-1
      vkxi = cross( vk, vi );
      vtmp = delta * vkxi;
      break;

    case 7: // k+1
    case 3: // k-1
      vixj = cross( vi, vj );
      vtmp = delta * vixj;
      break;
    }

    vgrad += vtmp;
  } // End of for loop

  // Set the results
  pdx = l_jacobian * vgrad.x;
  pdy = l_jacobian * vgrad.y;
  pdz = l_jacobian * vgrad.z;
}


__global__
void JacobianTermKernel( const GPU::Classes::VolumeArgGPU<float> x,
                         const GPU::Classes::VolumeArgGPU<float> y,
                         const GPU::Classes::VolumeArgGPU<float> z,
                         const GPU::Classes::VolumeArgGPU<float> area1,
                         const GPU::Classes::VolumeArgGPU<float> area2,
                         const GPU::Classes::VolumeArgGPU<float> origArea1,
                         const GPU::Classes::VolumeArgGPU<float> origArea2,
                         const GPU::Classes::VolumeArgGPU<char> invalid,
                         const float l_jacobian, const float exp_k,
                         const float maxNorm, const float jacScale,
                         GPU::Classes::VolumeArgGPU<float> dx,
                         GPU::Classes::VolumeArgGPU<float> dy,
                         GPU::Classes::VolumeArgGPU<float> dz )
{
  /*!
  Does the main work of gcamJacobianTerm
  */
  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  for( unsigned int iz = 0; iz < dx.dims.z; iz++ )
  {
    // Only compute if ix, iy & iz are inside the bounding box
    if( dx.InVolume(ix,iy,iz) )
    {

      // Check for validity
      if( invalid(ix,iy,iz) == GCAM_POSITION_INVALID )
      {
        continue;
      }

      float ndx, ndy, ndz;

      JacobianTermAtNode( x, y, z,
                          area1, area2,
                          origArea1, origArea2,
                          invalid,
                          ix, iy, iz,
                          l_jacobian, exp_k,
                          ndx, ndy, ndz );

      float norm = (ndx*ndx) + (ndy*ndy) + (ndz*ndz);
      norm = sqrtf( norm );
      if( norm > maxNorm*jacScale && maxNorm > 0 && norm > 1 )
      {
        ndx *= maxNorm / norm;
        ndy *= maxNorm / norm;
        ndz *= maxNorm / norm;
      }

      dx(ix,iy,iz) += ndx;
      dy(ix,iy,iz) += ndy;
      dz(ix,iy,iz) += ndz;
    }
  }
}





void GCAmorphTerm::Jacobian( GPU::Classes::GCAmorphGPU& gcam,
                             const float l_jacobian,
                             const float jac_scale ) const
{


  if( DZERO(l_jacobian) )
  {
    return;
  }

  GCAmorphTerm::tJacobTot.Start();

  const dim3 gcamDims = gcam.d_rx.GetDims();
  const unsigned int nVoxels = gcamDims.x * gcamDims.y * gcamDims.z;

  dim3 grid, threads;


  // Allocate Thrust arrays
  thrust::device_ptr<float> d_norms;
  d_norms = thrust::device_new<float>( nVoxels );

  // Compute the norms
  GCAmorphTerm::tJacobMaxNorm.Start();
  threads.x = threads.y = kGCAmorphJacobTermKernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphJacobTermKernelSize );
  grid.z = 1;

  NormKernel<<<grid,threads>>>( gcam.d_dx,
                                gcam.d_dy,
                                gcam.d_dz,
                                thrust::raw_pointer_cast( d_norms ) );
  CUDA_CHECK_ERROR( "NormKernel failed!\n" );

  float maxNorm = *thrust::max_element( d_norms, d_norms+nVoxels );
  GCAmorphTerm::tJacobMaxNorm.Stop();

  // Main computation
  GCAmorphTerm::tJacobCompute.Start();
  JacobianTermKernel<<<grid,threads>>>( gcam.d_rx, gcam.d_ry, gcam.d_rz,
                                        gcam.d_area1, gcam.d_area2,
                                        gcam.d_origArea1, gcam.d_origArea2,
                                        gcam.d_invalid,
                                        l_jacobian, gcam.exp_k,
                                        maxNorm, jac_scale,
                                        gcam.d_dx, gcam.d_dy, gcam.d_dz );
  CUDA_CHECK_ERROR( "JacobianTermKernel failed!" );
  GCAmorphTerm::tJacobCompute.Stop();


  // Release Thrust arrays
  thrust::device_delete( d_norms );


  GCAmorphTerm::tJacobTot.Stop();
}


// ##############################################################

template<typename T>
__device__ float FetchMRIvoxel( const float x,
                                const float y,
                                const float z )
{
  return( 10000 );
}

template<>
__device__ float FetchMRIvoxel<unsigned char>( const float x,
    const float y,
    const float z )
{
  float texVal;
  texVal = tex3D( dt_mri_uchar, x+0.5f, y+0.5f, z+0.5f );
  texVal *= UCHAR_MAX;

  return( texVal );
}

// ---

template<typename T>
__device__ float FetchMRIsmoothVoxel( const float x,
                                      const float y,
                                      const float z )
{
  return( 50000 );
}

template<>
__device__ float FetchMRIsmoothVoxel<unsigned char>( const float x,
    const float y,
    const float z )
{
  float texVal;
  texVal = tex3D( dt_mri_smooth_uchar, x+0.5f, y+0.5f, z+0.5f );
  texVal *= UCHAR_MAX;

  return( texVal );
}

// ---

template<typename T>
GPU::Classes::CTfactory* GCAmorphTerm::BindMRI( const GPU::Classes::MRIframeGPU<T>& mri ) const
{
  std::cerr << __PRETTY_FUNCTION__
            << ": Unrecognised MRI type" << std::endl;
  abort();
  return( NULL );
}

template<>
GPU::Classes::CTfactory* GCAmorphTerm::BindMRI<unsigned char>( const GPU::Classes::MRIframeGPU<unsigned char>& mri ) const
{
  GPU::Classes::CTfactory *tmp ;
  tmp = new GPU::Classes::CTfactory( mri, dt_mri_uchar,
                                     cudaFilterModeLinear );

  return( tmp );
}


// ---

template<typename T>
GPU::Classes::CTfactory* GCAmorphTerm::BindMRIsmooth( const GPU::Classes::MRIframeGPU<T>& mri ) const
{
  std::cerr << __PRETTY_FUNCTION__
            << ": Unrecognised MRI type" << std::endl;
  abort();
  return( NULL );
}

template<>
GPU::Classes::CTfactory* GCAmorphTerm::BindMRIsmooth<unsigned char>( const GPU::Classes::MRIframeGPU<unsigned char>& mri ) const
{
  GPU::Classes::CTfactory *tmp ;
  tmp = new GPU::Classes::CTfactory( mri, dt_mri_smooth_uchar,
                                     cudaFilterModeLinear );

  return( tmp );
}





// ---

template<typename U>
__device__
void SmoothGradient( const float x, const float y, const float z,
                     const float3 sizes,
                     float& dx, float& dy, float& dz )
{

  float xp1, xm1, yp1, ym1, zp1, zm1;

  xp1 = FetchMRIsmoothVoxel<U>( x+1, y, z );
  xm1 = FetchMRIsmoothVoxel<U>( x-1, y, z );

  yp1 = FetchMRIsmoothVoxel<U>( x, y+1, z );
  ym1 = FetchMRIsmoothVoxel<U>( x, y-1, z );

  zp1 = FetchMRIsmoothVoxel<U>( x, y, z+1 );
  zm1 = FetchMRIsmoothVoxel<U>( x, y, z-1 );

  dx = (xp1-xm1) / ( 2 * sizes.x );
  dy = (yp1-ym1) / ( 2 * sizes.y );
  dz = (zp1-zm1) / ( 2 * sizes.z );
}


template<typename T, typename U>
__global__
void LogLikelihoodKernel( const GPU::Classes::VolumeArgGPU<char> invalid,
                          const GPU::Classes::VolumeArgGPU<int> label,
                          const GPU::Classes::VolumeArgGPU<int> status,
                          const GPU::Classes::VolumeArgGPU<float> rx,
                          const GPU::Classes::VolumeArgGPU<float> ry,
                          const GPU::Classes::VolumeArgGPU<float> rz,
                          const GPU::Classes::VolumeArgGPU<float> mean,
                          const GPU::Classes::VolumeArgGPU<float> variance,
                          const float3 mriSizes,
                          const float l_log_likelihood,
                          GPU::Classes::VolumeArgGPU<float> dx,
                          GPU::Classes::VolumeArgGPU<float> dy,
                          GPU::Classes::VolumeArgGPU<float> dz )
{

  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  // Loop over z slices
  for( unsigned int iz = 0; iz < rx.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( rx.InVolume(ix,iy,iz) )
    {

      // Check validity
      if( invalid(ix,iy,iz) == GCAM_POSITION_INVALID )
      {
        continue;
      }

      // Check status
      if( status(ix,iy,iz) &
          (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD) )
      {
        continue;
      }

      // Only use unknowns if they border knowns
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

      // This ends the checks

      float mriVal = FetchMRIvoxel<T>( rx(ix,iy,iz),
                                       ry(ix,iy,iz),
                                       rz(ix,iy,iz) );

      float vMean, invVariance;

      // This check is equivalent to "if( gcamn->gc )"
      if( variance(ix,iy,iz) >= 0 )
      {
        vMean = mean(ix,iy,iz);
        invVariance = 1/variance(ix,iy,iz);
      }
      else
      {
        vMean = 0;
        invVariance = 1.0f/MIN_VAR;
      }

      float mydx, mydy, mydz;

      SmoothGradient<U>( rx(ix,iy,iz), ry(ix,iy,iz), rz(ix,iy,iz),
                         mriSizes,
                         mydx, mydy, mydz );

      float norm = sqrtf( mydx*mydx + mydy*mydy + mydz*mydz );
      if( !FZERO(norm) )
      {
        mydx /= norm;
        mydy /= norm;
        mydz /= norm;
      }

      vMean -= mriVal;
#define MAX_ERROR 1000.0f
      if( fabs( vMean ) > MAX_ERROR )
      {
        vMean = copysign( MAX_ERROR, vMean );
      }

      vMean *= invVariance;

      if( IS_UNKNOWN( label(ix,iy,iz) ) )
      {
        if( FZERO(mriVal) )
        {
          if( vMean >= 0.5f )
          {
            vMean = 0.5f;
          }
        }
        else
        {
          if( vMean > 2 )
          {
            vMean = 2;
          }
        }
      }


      dx(ix,iy,iz) += l_log_likelihood * mydx * vMean;
      dy(ix,iy,iz) += l_log_likelihood * mydy * vMean;
      dz(ix,iy,iz) += l_log_likelihood * mydz * vMean;

    }
  }
}





template<typename T, typename U>
void GCAmorphTerm::LogLikelihood( GPU::Classes::GCAmorphGPU& gcam,
                                  const GPU::Classes::MRIframeGPU<T>& mri,
                                  const GPU::Classes::MRIframeGPU<U>& mri_smooth,
                                  double l_log_likelihood ) const
{
  /*!
  Computes the Log Likelihood term on the GPU.
  Assumes that both MRI data structures are already in their
  cudaArrays.
  */
  if( DZERO(l_log_likelihood) )
  {
    return;
  }

  GCAmorphTerm::tLogLikelihoodTot.Start();

  // Get the MRI textures set up (assumes MRIs already in cudaArrays)
  std::auto_ptr<GPU::Classes::CTfactory> mriCT( this->BindMRI( mri ) );
  std::auto_ptr<GPU::Classes::CTfactory> mri_smoothCT( this->BindMRIsmooth( mri_smooth ) );

  // Run the computation
  dim3 threads, grid;
  threads.x = threads.y = kGCAmorphLLTermKernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphLLTermKernelSize );
  grid.z = 1;

  GCAmorphTerm::tLogLikelihoodCompute.Start();
  LogLikelihoodKernel<T,U><<<grid,threads>>>( gcam.d_invalid,
      gcam.d_label,
      gcam.d_status,
      gcam.d_rx,
      gcam.d_ry,
      gcam.d_rz,
      gcam.d_mean,
      gcam.d_variance,
      mri_smooth.GetSizes(),
      l_log_likelihood,
      gcam.d_dx,
      gcam.d_dy,
      gcam.d_dz );
  CUDA_CHECK_ERROR( "LogLikelihoodKernel failed!" );
  GCAmorphTerm::tLogLikelihoodCompute.Stop();


  GCAmorphTerm::tLogLikelihoodTot.Stop();
}

// --

template<typename T, typename U>
void GCAmorphTerm::LLtermDispatch( GCA_MORPH *gcam,
                                   const MRI* mri,
                                   const MRI* mri_smooth,
                                   double l_log_likelihood ) const
{
  GPU::Classes::GCAmorphGPU myGCAM;
  GPU::Classes::MRIframeGPU<T> myMRI;
  GPU::Classes::MRIframeGPU<U> myMRIsmooth;

  // Set up the GCAmorph
  myGCAM.SendAll( gcam );

  // Handle the MRIs
  myMRI.Allocate( mri );
  myMRI.Send( mri, 0 );

  myMRIsmooth.Allocate( mri_smooth );
  myMRIsmooth.Send( mri_smooth, 0 );


  // Run the computation
  this->LogLikelihood( myGCAM, myMRI, myMRIsmooth, l_log_likelihood );

  // Retrieve results
  myGCAM.RecvAll( gcam );
}

// --

template<typename T>
void GCAmorphTerm::LLTmrismoothDispatch( GCA_MORPH *gcam,
    const MRI*  mri,
    const MRI* mri_smooth,
    double l_log_likelihood ) const
{

  switch( mri_smooth->type )
  {

  case MRI_UCHAR:
    this->LLtermDispatch<T,unsigned char>( gcam, mri, mri_smooth,
                                           l_log_likelihood );
    break;

  case MRI_FLOAT:
    this->LLtermDispatch<T,float>( gcam, mri, mri_smooth,
                                   l_log_likelihood );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised type for mri_smooth "
              << mri_smooth->type << std::endl;
    exit( EXIT_FAILURE );
  }

}

// --

void GCAmorphTerm::LLTDispatch( GCA_MORPH *gcam,
                                const MRI*  mri,
                                const MRI* mri_smooth,
                                double l_log_likelihood ) const
{

  switch( mri->type )
  {

  case MRI_UCHAR:
    this->LLTmrismoothDispatch<unsigned char>( gcam, mri, mri_smooth,
        l_log_likelihood );
    break;

  case MRI_FLOAT:
    this->LLTmrismoothDispatch<float>( gcam, mri, mri_smooth,
                                       l_log_likelihood );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised type for mri "
              << mri->type << std::endl;
    exit( EXIT_FAILURE );
  }

}


// ##############################################################

const int kMaxMLEdist = 1;

__global__
void LabelCopyDeltasKernel( const GPU::Classes::VolumeArgGPU<char> invalid,
                            const GPU::Classes::VolumeArgGPU<int> status,
                            GPU::Classes::VolumeArgGPU<float> dy,
                            GPU::Classes::VolumeArgGPU<float> labelDist,
                            const GPU::Classes::MRIframeOnGPU<float> mriDist ,
                            const float l_label )
{
  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  // Loop over z slices
  for( unsigned int iz = 0; iz < invalid.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( invalid.InVolume(ix,iy,iz) )
    {

      // Validity check
      if( invalid(ix,iy,iz) ||
          ( (status(ix,iy,iz) & GCAM_LABEL_NODE) == 0 ) )
      {
        continue;
      }

      float localdy = mriDist(ix,iy,iz);
      labelDist(ix,iy,iz) = localdy;

      if( fabs(localdy) > kMaxMLEdist )
      {
        localdy *= kMaxMLEdist / fabs(localdy);
      }

      dy(ix,iy,iz) += l_label * localdy;

    }
  }

}


void GCAmorphTerm::LabelCopyDeltas( GPU::Classes::GCAmorphGPU& gcam,
                                    const GPU::Classes::MRIframeGPU<float>& mri_dist,
                                    const float l_label ) const
{

  GCAmorphTerm::tLabelCopyDeltas.Start();

  // Run the kernel
  dim3 grid, threads;
  threads.x = threads.y = kGCAmorphLabelTermCopyDeltasKernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphLabelTermCopyDeltasKernelSize );
  grid.z = 1;

  LabelCopyDeltasKernel<<<grid,threads>>>( gcam.d_invalid,
      gcam.d_status,
      gcam.d_dy,
      gcam.d_labelDist,
      mri_dist,
      l_label );
  CUDA_CHECK_ERROR( "LabelCopyDeltasKernel failed!\n" );

  GCAmorphTerm::tLabelCopyDeltas.Stop();

}



// -------------------------------------------------------------------

__global__
void LabelPostAntConsistKernel( const GPU::Classes::VolumeArgGPU<char> invalid,
                                const GPU::Classes::VolumeArgGPU<float> dy,
                                GPU::Classes::VolumeArgGPU<int> status,
                                GPU::Classes::MRIframeOnGPU<float> mriDist,
                                int *nRemoved )
{
  /*!
  This kernel is awful. If you look, you'll see that each voxel
  depends on those on either side in both the x and z directions.
  It gives different results to the CPU version
  */
  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  // Loop over z slices
  for( unsigned int iz = 0; iz < invalid.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( invalid.InVolume(ix,iy,iz) )
    {

      // Validity check
      if( invalid(ix,iy,iz) ||
          ( (status(ix,iy,iz) & GCAM_LABEL_NODE) == 0 ) )
      {
        continue;
      }

      /*
         If signs of nodes to left and right are different,
         don't trust this one
      */
      if( (ix<invalid.dims.x-1) && (ix>0) )
      {

        if( ( (status(ix-1,iy,iz) & GCAM_LABEL_NODE)==0 ) ||
            ( (status(ix+1,iy,iz) & GCAM_LABEL_NODE)==0 ) )
        {
          // Only if they are both label nodes
          continue;
        }

        if( ( dy(ix-1,iy,iz) * dy(ix+1,iy,iz) > 0 ) &&
            ( dy(ix-1,iy,iz) * dy(ix,iy,iz) < 0 ) )
        {
          status(ix,iy,iz) = GCAM_USE_LIKELIHOOD;
          mriDist(ix,iy,iz) = 0;
          atomicAdd( nRemoved, 1 );
          continue;
        }

      }

      // Repeat in z direction
      if( (iz<invalid.dims.z-1) && (iz>0) )
      {

        if( ( (status(ix,iy,iz-1) & GCAM_LABEL_NODE)==0 ) ||
            ( (status(ix,iy,iz+1) & GCAM_LABEL_NODE)==0 ) )
        {
          continue;
        }

        if( ( dy(ix,iy,iz+1)*dy(ix,iy,iz-1) > 0 ) &&
            ( dy(ix,iy,iz+1)*dy(ix,iy,iz) < 0 ) )
        {
          status(ix,iy,iz) = GCAM_USE_LIKELIHOOD;
          mriDist(ix,iy,iz) = 0;
          atomicAdd( nRemoved, 1 );
        }
      }

    } // if inVolume
  } // iz Loop


}

int GCAmorphTerm::LabelPostAntConsistency( GPU::Classes::GCAmorphGPU& gcam,
    GPU::Classes::MRIframeGPU<float>& mri_dist ) const
{

  int nRemoved;
  int *d_nRemoved;

  GCAmorphTerm::tLabelPostAntConsistency.Start();

  // Allocate memory
  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_nRemoved, sizeof(int) ) );

  // Zero
  nRemoved = 0;
  CUDA_SAFE_CALL( cudaMemcpy( d_nRemoved, &nRemoved, sizeof(int),
                              cudaMemcpyHostToDevice ) );

  // Run the kernel
  dim3 threads, grid;
  threads.x = threads.y = kGCAmorphLabelTermPostAntConsistencyKernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphLabelTermPostAntConsistencyKernelSize );
  grid.z = 1;
  LabelPostAntConsistKernel<<<grid,threads>>>( gcam.d_invalid,
      gcam.d_dy,
      gcam.d_status,
      mri_dist,
      d_nRemoved );
  CUDA_CHECK_ERROR( "LabelPostAndConsistKernel failed!\n" );

  // Retrieve result
  CUDA_SAFE_CALL( cudaMemcpy( &nRemoved, d_nRemoved, sizeof(int),
                              cudaMemcpyDeviceToHost ) );

  // Release memory
  CUDA_SAFE_CALL( cudaFree( d_nRemoved ) );

  GCAmorphTerm::tLabelPostAntConsistency.Stop();

  return( nRemoved );
}

// -------------------------------------------------------------------

__global__
void LabelFinalKernel( const GPU::Classes::VolumeArgGPU<char> invalid,
                       const GPU::Classes::VolumeArgGPU<int> status,
                       const GPU::Classes::MRIframeOnGPU<float> mriDist,
                       GPU::Classes::VolumeArgGPU<float> dy,
                       GPU::Classes::VolumeArgGPU<float> labelDist,
                       const float l_label,
                       int *num )
{
  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  // Loop over z slices
  for( unsigned int iz = 0; iz < invalid.dims.z; iz++ )
  {

    // Only compute if ix, iy & iz are inside the bounding box
    if( invalid.InVolume(ix,iy,iz) )
    {

      if( invalid(ix,iy,iz) ||
          ( (status(ix,iy,iz) & GCAM_LABEL_NODE) == 0 ) )
      {
        continue;
      }

      dy(ix,iy,iz) = l_label * mriDist(ix,iy,iz);

      if( fabs( dy(ix,iy,iz)/l_label ) >= 1 )
      {
        atomicAdd( num, 1 );
      }

      labelDist(ix,iy,iz) = dy(ix,iy,iz);
    }
  }
}


int GCAmorphTerm::LabelFinalUpdate( GPU::Classes::GCAmorphGPU& gcam,
                                    const GPU::Classes::MRIframeGPU<float>& mri_dist,
                                    const float l_label ) const
{

  int num;
  int *d_num;

  GCAmorphTerm::tLabelFinal.Start();

  // Allocate memory for d_num
  CUDA_SAFE_CALL( cudaMalloc( (void**)&d_num, sizeof(int) ) );

  // Zero it
  num = 0;
  CUDA_SAFE_CALL( cudaMemcpy( d_num, &num, sizeof(int),
                              cudaMemcpyHostToDevice ) );


  // Run the kernel
  dim3 threads, grid;
  threads.x = threads.y = kGCAmorphLabelTermFinalKernelSize;
  threads.z = 1;

  grid = gcam.d_rx.CoverBlocks( kGCAmorphLabelTermFinalKernelSize );
  grid.z = 1;

  LabelFinalKernel<<<grid,threads>>>( gcam.d_invalid,
                                      gcam.d_status,
                                      mri_dist,
                                      gcam.d_dy,
                                      gcam.d_labelDist,
                                      l_label,
                                      d_num );
  CUDA_CHECK_ERROR( "LabelFinalKernel failed!\n" );



  // Retrieve d_num
  CUDA_SAFE_CALL( cudaMemcpy( &num, d_num, sizeof(int),
                              cudaMemcpyDeviceToHost ) );

  // Release d_num
  CUDA_SAFE_CALL( cudaFree( d_num ) );

  GCAmorphTerm::tLabelFinal.Stop();

  return( num );
}


// -------------------------------------------------------------------

int GCAmorphTerm::RemoveLabelOutliers( Freesurfer::GCAmorphCPU& gcam,
                                       MRI *mri_dist,
                                       const int whalf,
                                       const double thresh ) const
{

  int         nremoved, nremoved_total, n, i, vox_to_examine, niters;
  MRI         *mri_std, *mri_ctrl, *mri_tmp, *mri_ctrl_tmp;
  VOXEL_LIST  *vl;
  float       diff, val0, oval;
  int         del, xv, yv, zv, xo, yo, zo, x, y, z;
  double      max_change;

  GCAmorphTerm::tRemoveOutliers.Start();

  // Get hold of dimensions etc.
  gcam.CheckIntegrity();

  unsigned int width, height, depth;

  gcam.rx.GetDims( width, height, depth );
  const int nx = width;
  const int ny = height;
  const int nz = depth;

  Freesurfer::VolumeArgCPU<int> status( gcam.status );
  Freesurfer::VolumeArgCPU<float> dy( gcam.dy );


  mri_ctrl = MRIcloneDifferentType( mri_dist, MRI_UCHAR );

  for (z = 0 ; z < nz ; z++)
  {
    for (y = 0 ; y < ny ; y++)
    {
      for (x = 0 ; x < nx ; x++)
      {
        if( status(x,y,z) & GCAM_LABEL_NODE)
        {
          MRIsetVoxVal( mri_ctrl, x, y, z, 0, CONTROL_MARKED );
        }
      }
    }
  }

  niters = 100 ;
  for( nremoved_total = i = 0; i < niters; i++)
  {

    mri_std = MRIstdNonzero(mri_dist, NULL, mri_dist, whalf*2+1) ;
    vl = VLSTcreate(mri_std, .001, 1e10, NULL, 0, 0) ;
    VLSTsort(vl, vl) ;

    vox_to_examine = static_cast<int>( ceil( vl->nvox*.05 ) );
    for( nremoved = n = 0; n < vox_to_examine; n++ )
    {
      x = vl->xi[n] ;
      y = vl->yi[n] ;
      z = vl->zi[n];

      if( (status(x,y,z) & GCAM_LABEL_NODE) == 0 )
      {
        continue;
      }

      val0 = MRIgetVoxVal(mri_dist, x, y, z, 0) ;
      del = 0 ;

      for (xo = -whalf ; xo <= whalf && !del ; xo++)
      {
        xv = mri_dist->xi[x+xo] ;
        for (yo = -whalf ; yo <= whalf && !del; yo++)
        {
          yv = mri_dist->yi[y+yo] ;
          for (zo = -whalf ; zo <= whalf && !del ; zo++)
          {
            zv = mri_dist->zi[z+zo] ;
            oval = MRIgetVoxVal(mri_dist, xv, yv, zv, 0);

            if (!FZERO(oval))
            {
              diff = fabs(oval - val0) ;
              if (diff > thresh)   /* shouldn't ever be this big */
              {
                if (fabs(val0) > fabs(oval) && (val0*oval < 0))
                {
                  /*if their signs are different and difference is big */
                  del = 1 ;
                }
              }
            }
          }
        }
      }

      if( del )
      {
        nremoved++ ;
        MRIFvox(mri_dist, x, y, z) = 0 ;
        MRIsetVoxVal(mri_ctrl, x, y, z, 0, CONTROL_NONE) ;
        dy(x,y,z) = 0;

        // Inferior is y+1
        if( status(x,y+1,z) & GCAM_LABEL_NODE )
        {
          nremoved++ ;
          //  gcamn_inf->status = GCAM_USE_LIKELIHOOD ;
          MRIsetVoxVal(mri_dist, x, y+1, z, 0, 0) ;
          MRIsetVoxVal(mri_ctrl, x, y+1, z, 0, CONTROL_NONE) ;
          dy(x,y+1,z) = 0 ;
        }

        // Superior is y-1
        if( status(x,y-1,z) & GCAM_LABEL_NODE )
        {
          nremoved++ ;
          //          gcamn_sup->status = GCAM_USE_LIKELIHOOD ;
          dy(x,y-1,z) = 0 ;
          MRIsetVoxVal(mri_dist, x, y-1, z, 0, 0) ;
          MRIsetVoxVal(mri_ctrl, x, y-1, z, 0, CONTROL_NONE) ;
        }
      }
    }
#if 0
    if( true )
    {
      char fname[STRLEN] ;
      sprintf( fname, "gpu_dist_after%d.mgz",i );
      MRIwrite(mri_dist, fname) ;
    }
#endif

    nremoved_total += nremoved ;
    MRIfree(&mri_std) ;
    if( nremoved == 0 )
    {
      // nothing happened
      break ;
    }
    VLSTfree(&vl); // keep it for last iteration
  }




  /* now use soap bubble smoothing to estimate
  label offset of deleted locations */
  mri_tmp = NULL;
  mri_ctrl_tmp = NULL;

  for (i = 0 ; i < 100 ; i++)
  {
    max_change = 0.0 ;
    mri_tmp = MRIcopy(mri_dist, mri_tmp);
    mri_ctrl_tmp = MRIcopy( mri_ctrl, mri_ctrl_tmp );

    for( z = 0; z < nz; z++ )
    {
      for( y = 0; y < ny; y++ )
      {
        for( x = 0; x < nx; x++ )
        {
          int    xi, yi, zi, xk, yk, zk, num ;
          double mean ;

          if( (MRIgetVoxVal(mri_ctrl, x, y, z, 0) == CONTROL_MARKED) ||
              ((status(x,y,z) & GCAM_LABEL_NODE) == 0) )
          {
            continue;
          }

          for (xk = -1, num = 0, mean = 0.0 ; xk <= 1 ; xk++)
          {
            xi = mri_ctrl->xi[x+xk] ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_ctrl->yi[y+yk] ;
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                zi = mri_ctrl->zi[z+zk];


                if( (MRIgetVoxVal(mri_ctrl, xi, yi, zi, 0) == CONTROL_MARKED) ||
                    (MRIgetVoxVal(mri_ctrl, xi, yi, zi, 0) == CONTROL_NBR) )
                {
                  mean += MRIgetVoxVal(mri_dist, xi, yi, zi, 0) ;
                  num++ ;
                }
              }
            }
          }
          if (num > 0)
          {
            float val ;
            val = MRIgetVoxVal(mri_tmp, x, y, z, 0) ;
            mean /= num ;
            if (fabs(val-mean) > max_change)
            {
              max_change = fabs(val-mean) ;
            }
            MRIsetVoxVal(mri_tmp, x, y, z, 0, mean) ;
            MRIsetVoxVal(mri_ctrl_tmp, x, y, z, 0, CONTROL_TMP) ;
            status(x,y,z) = (GCAM_IGNORE_LIKELIHOOD | GCAM_LABEL_NODE) ;
          }
        }
      }
    }

    MRIcopy(mri_tmp, mri_dist) ;
    MRIcopy( mri_ctrl_tmp, mri_ctrl );

#if 0
    if( true )
    {
      char fname[STRLEN] ;
      sprintf( fname, "gpu_dist_after%d.mgz",i );
      MRIwrite(mri_dist, fname) ;
    }
#endif

    MRIreplaceValuesOnly(mri_ctrl, mri_ctrl, CONTROL_TMP, CONTROL_NBR) ;
    if( max_change < 0.05 )
    {
      break;
    }

  }


  MRIfree( &mri_ctrl );
  MRIfree( &mri_tmp );
  MRIfree( &mri_ctrl_tmp );

  GCAmorphTerm::tRemoveOutliers.Stop();

  return( nremoved );
}


int GCAmorphTerm::RemoveLabelOutliersDispatch( GPU::Classes::GCAmorphGPU& gcam,
    MRI *mri_dist,
    const int whalf,
    const double thresh ) const
{
  Freesurfer::GCAmorphCPU myGCAM;

  int nRemoved;

  myGCAM.AllocateFromTemplate( gcam );
  myGCAM.GetFromGPU( gcam );

  nRemoved = this->RemoveLabelOutliers( myGCAM, mri_dist, whalf, thresh );

  myGCAM.PutOnGPU( gcam );

  return( nRemoved );
}



// -------------------------------------------------------------------


#define MAX_MLE_DIST 1

void GCAmorphTerm::LabelMainLoop( Freesurfer::GCAmorphCPU& gcam,
                                  const MRI *mri,
                                  MRI *mri_dist,
                                  const double l_label,
                                  const double label_dist ) const
{

  int x, y, z;
  int xn, yn, zn;
  int best_label, sup_wm, sup_ven, wm_label;
  int n;
  double dy;
  GCA_NODE *gcan;
  GC1D *wm_gc;
  float vals[MAX_GCA_INPUTS], yi, yk, min_dist;


  GCAmorphTerm::tLabelMainLoop.Start();

  // Get hold of dimensions etc.
  gcam.CheckIntegrity();

  unsigned int width, height, depth;

  gcam.rx.GetDims( width, height, depth );
  const int nx = width;
  const int ny = height;
  const int nz = depth;

  const Freesurfer::VolumeArgCPU<char> invalid( gcam.invalid );
  const Freesurfer::VolumeArgCPU<int> label( gcam.label );
  Freesurfer::VolumeArgCPU<int> status( gcam.status );
  const Freesurfer::VolumeArgCPU<float> rx( gcam.rx );
  const Freesurfer::VolumeArgCPU<float> ry( gcam.ry );
  const Freesurfer::VolumeArgCPU<float> rz( gcam.rz );
  Freesurfer::VolumeArgCPU<float> gcamdy( gcam.dy );
  Freesurfer::VolumeArgCPU<float> labelDist( gcam.labelDist );


  for( z=0; z < nz; z++ )
  {
    for( y=0; y < ny; y++ )
    {
      for( x=0; x < nx; x++ )
      {

        const char& currInvalid( invalid(x,y,z) );

        // Skip invalid nodes
        if( currInvalid == GCAM_POSITION_INVALID )
        {
          continue;
        }

        // Can't do top or bottom layers
        if( (y == ny-1) || (y == 0) )
        {
          continue;
        }

        const float& rxCurr( rx(x,y,z) );
        const float& ryCurr( ry(x,y,z) );
        const float& rzCurr( rz(x,y,z) );

        /*
          only process nodes which are hippocampus
          superior to something else,
          or white matter inferior to hippocampus
        */
        if( y == 0 || y == ny-1 ||
            ryCurr == 0 || ryCurr == mri->height-1 )
        {
          continue;
        }

        const int& labelCurr( label(x,y,z) );
        const int& labelInf( label(x,y+1,z) );
        const int& labelSup( label(x,y-1,z) );

        if( !IS_HIPPO(labelCurr) && !IS_WM(labelCurr) )
        {
          continue;
        }

        if ( !IS_WM(labelCurr) )     /* only do white matter for now */
        {
          continue;
        }

        if (
          ((IS_HIPPO(labelCurr) && IS_WM(labelInf)) ||
           (IS_WM(labelCurr) && IS_HIPPO(labelSup))) == 0)
        {
          continue ;  /* only hippo above wm, or wm below hippo */
        }

        if( IS_HIPPO(labelCurr) )
        {
          load_vals( mri,
                     rxCurr, ryCurr+1, rzCurr,
                     vals, 1 ); // Last argument is ninputs == 1
        }
        else
        {
          load_vals( mri,
                     rxCurr, ryCurr, rzCurr,
                     vals, 1); // Last argument is ninputs == 1
        }

        if( GCApriorToNode( gcam.gca,
                            x, y, z,
                            &xn, &yn, &zn ) != NO_ERROR )
        {
          continue ;
        }

        if( (IS_HIPPO(labelCurr) && labelCurr == Left_Hippocampus) ||
            (IS_WM(labelCurr) &&
             labelCurr == Left_Cerebral_White_Matter))
        {
          wm_label = Left_Cerebral_White_Matter ;
        }
        else
        {
          wm_label = Right_Cerebral_White_Matter ;
        }

        wm_gc = GCAfindPriorGC( gcam.gca, x, y, z, wm_label );

        if (wm_gc == NULL)
        {
          continue;
        }

        gcan = GCAbuildRegionalGCAN( gcam.gca, xn, yn, zn, 3 );

        // ventral DC is indistinguishible from temporal wm pretty much
        for (n = 0 ; n < gcan->nlabels; n++)
        {
          if ((gcan->labels[n] == Left_VentralDC ||
               gcan->labels[n] == Right_VentralDC ||
               gcan->labels[n] == Brain_Stem) &&
              gcan->gcs[n].means[0] > 90)
          {
            gcan->labels[n] = wm_label ;
          }
        }

        dy = 0;
        min_dist = label_dist+1 ;
        sup_ven = sup_wm = 0;  /* if can't find any wm superior, then
              must be partial volume and don't trust */
#define SAMPLE_DIST 0.1
        for (yk = -label_dist ; yk <= label_dist ; yk += SAMPLE_DIST)
        {
          yi = ryCurr+yk ;   /* sample inferiorly */
          if ((yi >= (mri->height-1)) || (yi <= 0))
          {
            break ;
          }

          load_vals(mri, rxCurr, yi, rzCurr, vals, 1 );

          best_label = GCAmaxLikelihoodLabel(gcan, vals, 1, NULL) ;
          if( yk < 0 )
          {
            if( IS_CSF(best_label) )
            {
              sup_ven++ ;
            }
            else if (sup_ven < 3/SAMPLE_DIST)
            {
              sup_ven = 0 ;
            }
          }


          if( IS_CSF(best_label) &&
              sup_ven > 2/SAMPLE_DIST && yk < 0)
          {
            /* shouldn't have to go
               through CSF to get to
               wm superiorly */
            min_dist = label_dist+1 ;
          }

          if( best_label != labelCurr )
          {
            continue ;
          }

          if (yk < 0 && IS_WM(best_label))
          {
            sup_wm = 1 ;
          }

          if (fabs(yk) < fabs(min_dist))
          {
            if( GCAmorphTerm::IsTemporalWM( mri, gcan,
                                            rxCurr, yi, rzCurr ) )
            {
              min_dist = yk ;
            }
          }
        }

        /* if inferior to lateral ventricle (as opposed to
           temporal horn) and can't find any
           wm above then it means the wm is partial-volumed
           and don't trust estimate */
        if (sup_ven && sup_wm == 0)
        {
          min_dist = label_dist+1 ;
        }

        if (min_dist > label_dist)  /* couldn't find any labels that match */
        {
          double log_p, max_log_p ;

          /* wm may be partial volumed - look in smaller
             nbhd for most likely location of wm */
          min_dist = 0 ;
          max_log_p = -1e20 ;
          for (yk = -label_dist/3 ; yk <= label_dist/3 ; yk += SAMPLE_DIST)
          {
            yi = ryCurr+yk ;
            if ((yi >= (mri->height-1)) || (yi <= 0))
            {
              break ;
            }

            load_vals(mri, rxCurr, yi, rzCurr, vals, 1 );
            log_p = GCAcomputeConditionalLogDensity
                    (wm_gc, vals, 1, wm_label);
            if (log_p > max_log_p)
            {
              max_log_p = log_p ;
              min_dist = yk ;
            }
          }
        }
        else   /* adjust estimated position to be at most likely value */
        {
          double log_p, max_log_p ;
          double  ykmin, ykmax ;

          /* wm may be partial volumed - look in smaller
             nbhd for most likely location of wm */
          max_log_p = -1e20 ;
          ykmin = min_dist-1 ;
          ykmax = min_dist+1 ;
          for (yk = ykmin ; yk <= ykmax ; yk += SAMPLE_DIST)
          {
            yi = ryCurr+yk ;   /* sample inferiorly */
            if ((yi >= (mri->height-1)) || (yi <= 0))
            {
              break ;
            }

            load_vals(mri, rxCurr, yi, rzCurr, vals, 1 );
            log_p =
              GCAcomputeConditionalLogDensity
              (wm_gc, vals, 1, wm_label);
            if (log_p > max_log_p)
            {
              max_log_p = log_p ;
              min_dist = yk ;
            }
          }
        }


        labelDist(x,y,z) = MRIFvox(mri_dist, x, y, z) = min_dist ;
#if 1
        if (fabs(min_dist) > MAX_MLE_DIST)
        {
          min_dist = MAX_MLE_DIST * min_dist / fabs(min_dist) ;
        }
#endif
        dy = min_dist ;
        if( !FZERO(min_dist) )
        {
          status(x,y,z) = (GCAM_IGNORE_LIKELIHOOD | GCAM_LABEL_NODE) ;

          if( IS_WM(labelInf ) &&
              ((status(x,y+1,z) & GCAM_LABEL_NODE)==0))
          {
            status(x,y+1,z) = (GCAM_IGNORE_LIKELIHOOD | GCAM_LABEL_NODE);
            gcamdy(x,y+1,z) += (l_label)*dy ;
            labelDist(x,y+1,z) =
              MRIFvox(mri_dist, x, y+1, z) = labelDist(x,y,z);

          }

          if( IS_HIPPO(labelSup) &&
              ((status(x,y-1,z) & GCAM_LABEL_NODE)==0))
          {
            status(x,y-1,z) = (GCAM_IGNORE_LIKELIHOOD | GCAM_LABEL_NODE) ;
            gcamdy(x,y-1,z) += (l_label)*dy ;
            labelDist(x,y-1,z) =
              MRIFvox(mri_dist, x, y-1, z) = labelDist(x,y,z);

          }
        }
        GCAfreeRegionalGCAN(&gcan) ;
      }
    }
  }


  GCAmorphTerm::tLabelMainLoop.Stop();
}


// --------------------------------------

#define MAX_TEMPORAL_WM 3

int GCAmorphTerm::IsTemporalWM( const MRI *mri,
                                const GCA_NODE *gcan,
                                float xf, float yf, float zf )
{
  /*!
  Re-implementation of is_temporarl_wm, but with ninputs set
  to 1.
  */
  int  yk, label, nwhite ;
  float vals[MAX_GCA_INPUTS], yi ;


  nwhite = 0 ;
  for( yk = 1; yk < 2*MAX_TEMPORAL_WM ; yk++)
  {
    yi = yf-yk ;
    if (yi < 0)
    {
      break;
    }

    load_vals( mri, xf, yi, zf, vals, 1 );
    label = GCAmaxLikelihoodLabel( gcan, vals, 1, NULL );
    if( IS_WM(label) || IS_THALAMUS(label) )
    {
      nwhite++ ;
    }
  }

  for( yk = -1; yk > -2*MAX_TEMPORAL_WM; yk-- )
  {
    yi = yf-yk ;
    if (yi < 0 || yi >= mri->height)
    {
      break ;
    }

    load_vals( mri, xf, yi, zf, vals, 1 );
    label = GCAmaxLikelihoodLabel( gcan, vals, 1, NULL );
    if( IS_WM(label) || IS_THALAMUS(label) )
    {
      nwhite++;
    }
    else
    {
      /*
        if moving inferiorly and find non-wm voxel,
        then stop counting - should be hippo
      */
      break;

    }
  }
  /* must be main body of white matter - too much of it */
  return( nwhite <= MAX_TEMPORAL_WM );
}




void GCAmorphTerm::LabelMainLoopDispatch( GPU::Classes::GCAmorphGPU& gcam,
    const MRI *mri,
    MRI *mri_dist,
    const double l_label,
    const double label_dist ) const
{

  Freesurfer::GCAmorphCPU myGCAM;

  myGCAM.AllocateFromTemplate( gcam );
  myGCAM.GetFromGPU( gcam );

  this->LabelMainLoop( myGCAM, mri, mri_dist, l_label, label_dist );

  myGCAM.PutOnGPU( gcam );

}


// ##############################################################

void GCAmorphTerm::ShowTimings( void )
{

#ifdef CUDA_SHOW_TIMINGS
  std::cout << "==================================" << std::endl;
  std::cout << "GPU GCAmorph term timers" << std::endl;
  std::cout << "------------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
  std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
  std::cout << "Timings may not be accurate" << std::endl;
#endif
  std::cout << std::endl;

  std::cout << "Smoothness:" << std::endl;
  std::cout << "  Subtraction : " << GCAmorphTerm::tSmoothSubtract
            << std::endl;
  std::cout << "      Compute : " << GCAmorphTerm::tSmoothCompute
            << std::endl;
  std::cout << "Total         : " << GCAmorphTerm::tSmoothTot
            << std::endl;
  std::cout << std::endl;


  std::cout << "Jacobian:" << std::endl;
  std::cout << "    Max. Norm : " << GCAmorphTerm::tJacobMaxNorm
            << std::endl;
  std::cout << "      Compute : " << GCAmorphTerm::tJacobCompute
            << std::endl;
  std::cout << "Total         : " << GCAmorphTerm::tJacobTot
            << std::endl;
  std::cout << std::endl;


  std::cout << "Log Likelihood:" << std::endl;
  std::cout << "      Compute : " << GCAmorphTerm::tLogLikelihoodCompute
            << std::endl;
  std::cout << "Total         : " << GCAmorphTerm::tLogLikelihoodTot
            << std::endl;
  std::cout << std::endl;

  std::cout << "Label Term:" << std::endl;
  std::cout << "        Main Loop:" << GCAmorphTerm::tLabelMainLoop
            << std::endl;
  std::cout << "  Remove Outliers:" << GCAmorphTerm::tRemoveOutliers
            << std::endl;
  std::cout << "      Copy Deltas:" << GCAmorphTerm::tLabelCopyDeltas
            << std::endl;
  std::cout << "    Post/Ant Cons:" << GCAmorphTerm::tLabelPostAntConsistency
            << std::endl;
  std::cout << "     Final Update:" << GCAmorphTerm::tLabelFinal
            << std::endl;
  std::cout << "Total            :" << GCAmorphTerm::tLabelTot
            << std::endl;

  std::cout << "==================================" << std::endl;
#endif
}


// ##############################################################

// Declare static members

SciGPU::Utilities::Chronometer GCAmorphTerm::tSmoothTot;
SciGPU::Utilities::Chronometer GCAmorphTerm::tSmoothSubtract;
SciGPU::Utilities::Chronometer GCAmorphTerm::tSmoothCompute;


SciGPU::Utilities::Chronometer GCAmorphTerm::tJacobTot;
SciGPU::Utilities::Chronometer GCAmorphTerm::tJacobMaxNorm;
SciGPU::Utilities::Chronometer GCAmorphTerm::tJacobCompute;

SciGPU::Utilities::Chronometer GCAmorphTerm::tLogLikelihoodTot;
SciGPU::Utilities::Chronometer GCAmorphTerm::tLogLikelihoodCompute;

SciGPU::Utilities::Chronometer GCAmorphTerm::tLabelMainLoop;
SciGPU::Utilities::Chronometer GCAmorphTerm::tRemoveOutliers;
SciGPU::Utilities::Chronometer GCAmorphTerm::tLabelCopyDeltas;
SciGPU::Utilities::Chronometer GCAmorphTerm::tLabelPostAntConsistency;
SciGPU::Utilities::Chronometer GCAmorphTerm::tLabelFinal;

SciGPU::Utilities::Chronometer GCAmorphTerm::tLabelTot;
}
}





// --------------------------------------------
static GPU::Algorithms::GCAmorphTerm myTerms;


//! Wrapper around GPU class for smoothness term
void gcamSmoothnessTermGPU( GCA_MORPH *gcam,
                            const float l_smoothness )
{

  GPU::Classes::GCAmorphGPU myGCAM;

  myGCAM.SendAll( gcam );

  myTerms.Smoothness( myGCAM, l_smoothness );

  myGCAM.RecvAll( gcam );

}

//! Wrapper around GPU class for jacobian term
void gcamJacobianTermGPU( GCA_MORPH *gcam,
                          const float l_jacobian,
                          const float jac_scale )
{
  GPU::Classes::GCAmorphGPU myGCAM;

  myGCAM.SendAll( gcam );

  myTerms.Jacobian( myGCAM, l_jacobian, jac_scale );

  myGCAM.RecvAll( gcam );
}

//! Wrapper around GPU class for LogLikelihood term
void gcamLogLikelihoodTermGPU( GCA_MORPH *gcam,
                               const MRI *mri,
                               const MRI *mri_smooth,
                               double l_log_likelihood )
{

  myTerms.LLTDispatch( gcam, mri, mri_smooth, l_log_likelihood );
}



//! Wrapper around GPU class for LabelTermFinal
int gcamLabelTermFinalUpdateGPU( GCA_MORPH *gcam,
                                 const MRI* mri_dist,
                                 const double l_label )
{

  GPU::Classes::GCAmorphGPU myGCAM;
  GPU::Classes::MRIframeGPU<float> mriDist;

  int num;

  myGCAM.SendAll( gcam );

  mriDist.Allocate( mri_dist );
  mriDist.Send( mri_dist, 0 );

  num = myTerms.LabelFinalUpdate( myGCAM, mriDist, l_label );

  myGCAM.RecvAll( gcam );

  return( num );
}

//! Wrapper around GPU class for LabelPostAntConsistency
int gcamLabelTermPostAntConsistencyGPU( GCA_MORPH *gcam,
                                        MRI* mri_dist )
{

  GPU::Classes::GCAmorphGPU myGCAM;
  GPU::Classes::MRIframeGPU<float> mriDist;

  int num;

  myGCAM.SendAll( gcam );

  mriDist.Allocate( mri_dist );
  mriDist.Send( mri_dist, 0 );

  num = myTerms.LabelPostAntConsistency( myGCAM, mriDist );

  myGCAM.RecvAll( gcam );
  mriDist.Recv( mri_dist, 0 );

  return( num );
}


//! Wrapper around GPU class for LabelCopyDeltas
void gcamLabelTermCopyDeltasGPU( GCA_MORPH *gcam,
                                 const MRI* mri_dist,
                                 const double l_label )
{

  GPU::Classes::GCAmorphGPU myGCAM;
  GPU::Classes::MRIframeGPU<float> mriDist;

  myGCAM.SendAll( gcam );
  mriDist.Allocate( mri_dist );
  mriDist.Send( mri_dist, 0 );

  myTerms.LabelCopyDeltas( myGCAM, mriDist, l_label );

  myGCAM.RecvAll( gcam );
}


//! Wrapper around GPU class for RemoveLablOutliers
int gcamRemoveLabelOutliersGPU( GCA_MORPH *gcam,
                                MRI* mri_dist,
                                const int whalf,
                                const double thresh )
{

  int nremoved;

  GPU::Classes::GCAmorphGPU myGCAM;

  myGCAM.SendAll( gcam );

  nremoved = myTerms.RemoveLabelOutliersDispatch( myGCAM,
             mri_dist,
             whalf,
             thresh );

  myGCAM.RecvAll( gcam );

  return( nremoved );
}


//! Wrapper around GPU class for Label Main Loop
void gcamLabelTermMainLoopGPU( GCA_MORPH *gcam, const MRI *mri,
                               MRI *mri_dist,
                               const double l_label, const double label_dist )
{

  GPU::Classes::GCAmorphGPU myGCAM;

  myGCAM.SendAll( gcam );

  myTerms.LabelMainLoopDispatch( myGCAM,
                                 mri, mri_dist,
                                 l_label, label_dist );

  myGCAM.RecvAll( gcam );
}


//! Wrapper around GPU class for whole of LabelTerm
int gcamLabelTermGPU( GCA_MORPH *gcam, const MRI *mri,
                      double l_label, double label_dist )
{

  int retVal;

  if( DZERO(l_label) )
  {
    return( NO_ERROR );
  }

  switch( mri->type )
  {
  case MRI_UCHAR:
    retVal = myTerms.LabelTermDispatch<unsigned char>( gcam, mri,
             l_label, label_dist );
    break;

  case MRI_FLOAT:
    retVal = myTerms.LabelTermDispatch<float>( gcam, mri,
             l_label, label_dist );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised destination MRI type " << mri->type
              << std::endl;
    exit( EXIT_FAILURE );
  }

  return( retVal );
}


#endif
