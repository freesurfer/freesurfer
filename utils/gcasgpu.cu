/**
 * @file  gcasgpu.cu
 * @brief Holds GCAS class for the GPU
 *
 * 
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/05 14:41:29 $
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
 *
 */


#include "cudacheck.h"

#include "gcasgpu.hpp"



// ====================================================

namespace GPU {

  namespace Classes {


    // ##################################################
    
    // Memory management

    void GCASonemeanGPU::Allocate( const unsigned int n ) {
      if( this->nSamplesAlloc < n ) {
	this->Release();

	CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_x),
				    n*sizeof(int) ) );
	CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_y),
				    n*sizeof(int) ) );
	CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_z),
				    n*sizeof(int) ) );

	CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_means),
				    n*sizeof(float) ) );
	CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_priors),
				    n*sizeof(float) ) );
	CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_covars),
				    n*sizeof(float) ) );

	this->nSamplesAlloc = n;
      }

      this->nSamples = n;
    }

    void GCASonemeanGPU::Release( void ) {
      if( this->nSamplesAlloc != 0 ) {
	CUDA_SAFE_CALL( cudaFree( d_x ) );
	CUDA_SAFE_CALL( cudaFree( d_y ) );
	CUDA_SAFE_CALL( cudaFree( d_z ) );
	CUDA_SAFE_CALL( cudaFree( d_means ) );
	CUDA_SAFE_CALL( cudaFree( d_priors ) );
	CUDA_SAFE_CALL( cudaFree( d_covars ) );
	this->nSamples = 0;
	this->nSamplesAlloc = 0;
      }
    }

  }
}
