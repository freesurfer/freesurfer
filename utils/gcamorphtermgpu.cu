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
 *    $Author: rge21 $
 *    $Date: 2010/06/15 19:31:01 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2008,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <thrust/device_new_allocator.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

#include "macros.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"

#include "cudatypeutils.hpp"

#include "gcamorphtermgpu.hpp"


//! Texure for 'invalid' in Smoothness
texture<char, 3, cudaReadModeElementType> dt_smooth_invalid;


// ==============================================================


namespace GPU {
  namespace Algorithms {

    const unsigned int kGCAmorphSmoothTermKernelSize = 16;

    const unsigned int kGCAmorphJacobTermKernelSize = 16;


    // ##############################################################

    //! Kernel to perform \f$c = a - b\f$ on volumes
    template<typename T>
    __global__
    void SubtractVolumes( const GPU::Classes::VolumeArgGPU<T> a,
			  const GPU::Classes::VolumeArgGPU<T> b,
			  GPU::Classes::VolumeArgGPU<T> c ) {

      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;

      for( unsigned int iz = 0; iz < a.dims.z; iz++ ) {
	if( a.InVolume(ix,iy,iz) ) {

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
			       const float l_smoothness ) {
      /*!
	Main computation of the smoothness term.
      */

      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;

      for( unsigned int iz = 0; iz < dx.dims.z; iz++ ) {
	// Only compute if ix, iy & iz are inside the bounding box
	if( dx.InVolume(ix,iy,iz) ) {

	  // Only calculate if we're not invalid
	  const char myValid = tex3D( dt_smooth_invalid,
				      ix+0.5f, iy+0.5f, iz+0.5f );
	  if( myValid != GCAM_POSITION_INVALID ) {
	    
	    const float myvx = vx(ix,iy,iz);
	    const float myvy = vy(ix,iy,iz);
	    const float myvz = vz(ix,iy,iz);

	    float ddx, ddy, ddz;
	    ddx = ddy = ddz = 0;

	    unsigned int num = 0;

	    // Re-order loop nest in gcamSmoothnessEnergy
	    for( int zk=-1; zk<=1; zk++ ) {
	      int zn = iz + zk;
	      zn = min( static_cast<int>(dx.dims.z-1), zn );
	      zn = max( 0, zn );
	      
	      for( int yk=-1; yk<=1; yk++ ) {
		int yn = iy + yk;
		yn = min( static_cast<int>(dx.dims.y-1), yn );
		yn = max( 0, yn );
		
		for( int xk=-1; xk<=1; xk++ ) {
		  int xn = ix + xk;
		  xn = min( static_cast<int>(dx.dims.x-1), xn );
		  xn = max( 0, xn );
		
		  // Don't include self
		  if( (!xk) && (!yk) && (!zk) ) {
		    continue;
		  }
		  
		  // Don't include invalid neighbours
		  
		  // Use texture to get boundary conditions
		  const char neighbourInvalid = tex3D( dt_smooth_invalid,
						       xn+0.5f,
						       yn+0.5f,
						       zn+0.5f );
		  if( neighbourInvalid == GCAM_POSITION_INVALID ) {
		    continue;
		  }


		  ddx += ( vx(xn,yn,zn) - myvx );
		  ddy += ( vy(xn,yn,zn) - myvy );
		  ddz += ( vz(xn,yn,zn) - myvz );

		  num++;

		}
	      }
	    }

	    if( num > 0 ) {
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
				   const float l_smoothness ) const {
      /*!
	Computes ther smoothness term of the given gcam, following
	the routine gcamSmoothnessTerm in gcamorph.c.
      */

      if( DZERO(l_smoothness) ) {
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
      dt_smooth_invalid.normalized = false;
      dt_smooth_invalid.addressMode[0] = cudaAddressModeClamp;
      dt_smooth_invalid.addressMode[1] = cudaAddressModeClamp;
      dt_smooth_invalid.addressMode[2] = cudaAddressModeClamp;
      dt_smooth_invalid.filterMode = cudaFilterModePoint;
      
      gcam.d_invalid.AllocateArray();
      gcam.d_invalid.SendArray();
      CUDA_SAFE_CALL( cudaBindTextureToArray( dt_smooth_invalid,
					      gcam.d_invalid.GetArray() ) );
      

      // Run the main kernel
      GCAmorphTerm::tSmoothCompute.Start();
      SmoothnessTermKernel<<<grid,threads>>>
	( vx, vy, vz, gcam.d_dx, gcam.d_dy, gcam.d_dz, l_smoothness );
      CUDA_CHECK_ERROR( "SmoothnessTermKernelFailed!" );
      GCAmorphTerm::tSmoothCompute.Stop();



      // Unbind textures
      CUDA_SAFE_CALL( cudaUnbindTexture( dt_smooth_invalid ) );


      GCAmorphTerm::tSmoothTot.Stop();

    }

    // ##############################################################

    __global__
    void NormKernel( const GPU::Classes::VolumeArgGPU<float> dx,
		     const GPU::Classes::VolumeArgGPU<float> dy,
		     const GPU::Classes::VolumeArgGPU<float> dz,
		     float* norms ) {

      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;

      float myNorm;

      // Loop over z slices
      for( unsigned int iz = 0; iz < dx.dims.z; iz++ ) {

	// Only compute if ix, iy & iz are inside the bounding box
	if( dx.InVolume(ix,iy,iz) ) {

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
			     float& pdx, float& pdy, float& pdz ) {
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

      for( unsigned int n=0; n<kAreaNeighbours; n++ ) {
	invert = 1;
	
	switch( n ) {
	default:
	case 0: // Central node
	  if( (i+1 >= dims.x) || (j+1 >= dims.y) || (k+1 >= dims.z) ) {
	    continue; // Go to next iteration of for loop
	  }
	  nLoc = make_uint3(i,j,k);
	  niLoc = make_uint3(i+1,j,k);
	  njLoc = make_uint3(i,j+1,k);
	  nkLoc = make_uint3(i,j,k+1);
	  break;

	case 1: // i-1 node
	  if( (i == 0) || (j+1 >= dims.y) || (k+1 >= dims.z) ) {
	    continue;
	  }
	  nLoc = make_uint3(i-1,j,k);
	  niLoc = make_uint3(i,j,k);
	  njLoc = make_uint3(i-1,j+1,k);
	  nkLoc = make_uint3(i-1,j,k+1);
	  break;

	case 2: // j-1 node
	  if( (i+1 >= dims.x) || (j == 0) || (k+1 >= dims.z) ) {
	    continue;
	  }
	  nLoc = make_uint3(i,j-1,k);
	  niLoc = make_uint3(i+1,j-1,k);
	  njLoc = make_uint3(i,j,k);
	  nkLoc = make_uint3(i,j-1,k+1);
	  break;

	case 3: // k-1 node
	  if( (i+1 >= dims.x) || (j+1 >= dims.y) || (k == 0) ) {
	    continue;
	  }
	  nLoc = make_uint3(i,j,k-1);
	  niLoc = make_uint3(i+1,j,k-1);
	  njLoc = make_uint3(i,j+1,k-1);
	  nkLoc = make_uint3(i,j,k);
	  break;

	case 4: // Left-handed central node
	  if( (i==0) || (j==0) || (k==0) ) {
	    continue; // Go to next iteration of for loop
	  }
	  invert = -1;
	  nLoc = make_uint3(i,j,k);
	  niLoc = make_uint3(i-1,j,k);
	  njLoc = make_uint3(i,j-1,k);
	  nkLoc = make_uint3(i,j,k-1);
	  break;

	case 5: // i+1 node
	  if( (i+1>= dims.x) || (j==0) || (k==0) ) {
	    continue;
	  }
	  invert = -1;
	  nLoc = make_uint3(i+1,j,k);
	  niLoc = make_uint3(i,j,k);
	  njLoc = make_uint3(i+1,j-1,k);
	  nkLoc = make_uint3(i+1,j,k-1);
	  break;

	case 6: // j+1 node
	  if( (i==0) || (j+1>=dims.y) || (k==0) ) {
	    continue;
	  }
	  invert = -1;
	  nLoc = make_uint3(i,j+1,k);
	  niLoc = make_uint3(i-1,j+1,k);
	  njLoc = make_uint3(i,j,k);
	  nkLoc = make_uint3(i,j+1,k-1);
	  break;

	case 7: // k+1 node
	  if( (i==0) || (j==0) || (k+1>=dims.z) ) {
	    continue;
	  }
	  invert = -1;
	  nLoc = make_uint3(i,j,k+1);
	  niLoc = make_uint3(i-1,j,k+1);
	  njLoc = make_uint3(i,j-1,k+1);
	  nkLoc = make_uint3(i,j,k);
	  break;
	}

	if( invert > 0 ) {
	  // Right handed co-ordinate system
	  orig_area = origArea1( nLoc.x, nLoc.y, nLoc.z );
	  area = area1( nLoc.x, nLoc.y, nLoc.z );
	} else {
	  // Left handed co-ordinate system
	  orig_area = origArea2( nLoc.x, nLoc.y, nLoc.z );
	  area = area2( nLoc.x, nLoc.y, nLoc.z );
	}

	// Check for area (volume...?) close to zero
	if( fabs(orig_area) < 1e-6f ) {
	  continue;
	}

	// Check for validity
	if( (invalid(nLoc.x,nLoc.y,nLoc.z) == GCAM_POSITION_INVALID) ||
	    (invalid(niLoc.x,niLoc.y,niLoc.z) == GCAM_POSITION_INVALID) ||
	    (invalid(njLoc.x,njLoc.y,njLoc.z) == GCAM_POSITION_INVALID) ||
	    (invalid(nkLoc.x,nkLoc.y,nkLoc.z) == GCAM_POSITION_INVALID) ) {
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
	if( exponent > kMaxExp ) {
	  exponent = kMaxExp;
	}

	// Might want to consider taylor series here....
	float delta = ( invert * exp_k / orig_area ) * ( 1 / ( 1 + exp(exponent) ) );

	// Compute and accumulate cross products
	float3 vjxk, vkxi, vixj, vtmp;
	switch( n ) {
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
			     GPU::Classes::VolumeArgGPU<float> dz ) {
      /*!
	Does the main work of gcamJacobianTerm
      */
      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;

      for( unsigned int iz = 0; iz < dx.dims.z; iz++ ) {
	// Only compute if ix, iy & iz are inside the bounding box
	if( dx.InVolume(ix,iy,iz) ) {
	  
	  // Check for validity
	  if( invalid(ix,iy,iz) == GCAM_POSITION_INVALID ) {
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
	  if( norm > maxNorm*jacScale && maxNorm > 0 && norm > 1 ) {
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
				 const float jac_scale ) const {

      GCAmorphTerm::tJacobTot.Start();

      if( DZERO(l_jacobian) ) {
	return;
      }

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
    
    void GCAmorphTerm::ShowTimings( void ) {

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

  }
}





// --------------------------------------------
static GPU::Algorithms::GCAmorphTerm myTerms;


//! Wrapper around GPU class for smoothness term
void gcamSmoothnessTermGPU( GCA_MORPH *gcam,
			    const float l_smoothness ) {
  
  GPU::Classes::GCAmorphGPU myGCAM;
  
  myGCAM.SendAll( gcam );
  
  myTerms.Smoothness( myGCAM, l_smoothness );
  
  myGCAM.RecvAll( gcam );
  
}

//! Wrapper around GPU class for jacobian term
void gcamJacobianTermGPU( GCA_MORPH *gcam,
			  const float l_jacobian,
			  const float jac_scale ) {
  GPU::Classes::GCAmorphGPU myGCAM;
  
  myGCAM.SendAll( gcam );
  
  myTerms.Jacobian( myGCAM, l_jacobian, jac_scale );
  
  myGCAM.RecvAll( gcam );
}
