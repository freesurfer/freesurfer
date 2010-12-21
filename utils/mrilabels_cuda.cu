/**
 * @file  mrilabels_cuda.cu
 * @brief Holds various MRI 'label' routines for the GPU
 *
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/12/21 20:03:44 $
 *    $Revision: 1.1 $
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
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */
#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>

#include "mri.h"

#include "chronometer.hpp"
#include "cudacheck.h"
#include "mriframegpu.hpp"

#include "mrilabels_cuda.hpp"

// ==============================================


//! Texture for unsigned char source
texture<unsigned char, 3, cudaReadModeElementType> dt_src_uchar;


namespace GPU {
  namespace Algorithms {

    //! Helper function to get texture values
    __device__
    unsigned char FetchSrcVoxel( const int ix,
				 const int iy,
				 const int iz ) {
      unsigned char texVal;

      texVal = tex3D( dt_src_uchar, ix+0.5f, iy+0.5f, iz+0.5f );

      return( texVal );
    }


    
    //! GPU kernel for MarkLabelBorderVoxel
    template<bool sixConnect>
    __global__
    void MarkLabelBorderVoxelKernel( GPU::Classes::MRIframeOnGPU<unsigned char> dst,
				     const int label,
				     const int mark ) {
      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;
      
      // Loop over z
      for( unsigned int iz=0; iz<dst.dims.z; iz++ ) {

	// Ensure we're in range
	if( dst.InVolume(ix,iy,iz) ) {
	  const int this_label = FetchSrcVoxel( ix, iy, iz );
	  int border = 0;

	  // Loop over local volume
	  for( int xk=-1; xk<=1 && !border; xk++ ) {
	    for( int yk=-1; yk<=1 && !border; yk++ ) {
	      for( int zk=-1; zk<=1; zk++ ) {

		if( sixConnect && (abs(xk)+abs(yk)+abs(zk) != 1) ) {
		  continue;
		}

		const int that_label = FetchSrcVoxel( ix+xk, iy+yk, iz+zk );
		if( ((this_label == label) && (that_label != label)) ||
		    ((this_label != label) && (that_label == label)) ) {
		  border = 1 ;
		  break ;
		}
	      }
	    }
	  }
	  if( border ) {
	    dst(ix,iy,iz) = mark;
	  } else {
	    dst(ix,iy,iz) = 0;
	  }
	}
      }
    }
    

    // =============================================

    void MRIlabels::MarkLabelBorderVoxels( const GPU::Classes::MRIframeGPU<unsigned char>& src,
					   GPU::Classes::MRIframeGPU<unsigned char>& dst,
					   const int label,
					   const int mark,
					   const int sixConnect ) const {
      cudaArray* srcArr = NULL;

      MRIlabels::tMarkLabelBorderVoxelsTot.Start();

      // Verify we have output allocated
      dst.Allocate( src );

      // Send the input to a cuda Array
      srcArr = src.CreateArray();

      // Bind the array to the texture
      dt_src_uchar.normalized = false;
      dt_src_uchar.addressMode[0] = cudaAddressModeClamp;
      dt_src_uchar.addressMode[1] = cudaAddressModeClamp;
      dt_src_uchar.addressMode[2] = cudaAddressModeClamp;
      dt_src_uchar.filterMode = cudaFilterModePoint;

      CUDA_SAFE_CALL( cudaBindTextureToArray( dt_src_uchar, srcArr ) );

      // Run the kernel
      const unsigned int kKernelSize = 16;

      dim3 grid, threads;
      threads.x = threads.y = kKernelSize;
      threads.z = 1;

      grid = src.CoverBlocks( kKernelSize );
      grid.z = 1;

      if( sixConnect ) {
	MarkLabelBorderVoxelKernel<true><<<grid,threads>>>( dst, label, mark );
      } else {
	MarkLabelBorderVoxelKernel<false><<<grid,threads>>>( dst, label, mark );
      }
      CUDA_CHECK_ERROR( "MarkLabelBorderVoxelKernel failed!\n" );


      // Release the texture
      CUDA_SAFE_CALL( cudaUnbindTexture( dt_src_uchar ) );

      // Release the cudaArray
      CUDA_SAFE_CALL( cudaFreeArray( srcArr ) );

      MRIlabels::tMarkLabelBorderVoxelsTot.Stop();
    }



    // ==============================================================

    // Declare statics
    SciGPU::Utilities::Chronometer MRIlabels::tMarkLabelBorderVoxelsTot;


    void MRIlabels::ShowTimings( void ) {
      std::cout << "=============================================" << std::endl;
      std::cout << "GPU MRI Label timers" << std::endl;
      std::cout << "--------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
      std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
      std::cout << "Timings may not be accurate" << std::endl;
#endif
      std::cout << std::endl;

      std::cout << "MarkLabelBorderVoxels" << std::endl;
      std::cout << "Total : " << MRIlabels::tMarkLabelBorderVoxelsTot << std::endl;
      std::cout << "=============================================" << std::endl;

    }

  }
}

// ===============================================

static GPU::Algorithms::MRIlabels myLabels;

//! Wrapper for MRImarkLabelBorderVoxels
void MRImarkLabelBorderVoxelsGPU( const MRI* mri_src,
				  MRI* mri_dst,
				  int label,
				  int mark,
				  int six_connected ) {
  GPU::Classes::MRIframeGPU<unsigned char> srcGPU;
  GPU::Classes::MRIframeGPU<unsigned char> dstGPU;

  srcGPU.Allocate( mri_src );
  srcGPU.VerifyMRI( mri_src );

  srcGPU.SendFrame( mri_src, 0 );

  myLabels.MarkLabelBorderVoxels( srcGPU, dstGPU, label, mark, six_connected );

  dstGPU.RecvFrame( mri_dst, 0 );
}
