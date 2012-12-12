/**
 * @file  gcapriorgpu.hpp
 * @brief Class to hold a volume of GCA priors in linear memory on the GPU
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.3 $
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

#include "cudacheck.h"

#include "gcapriorgpu.hpp"


namespace GPU
{
namespace Classes
{

void GCApriorGPU::Allocate( const long long nxDim,
                            const long long nyDim,
                            const long long nzDim,
                            const size_t num4D )
{
  // Get rid of old allocation
  this->Release();

  // Copy sizes
  this->xDim = nxDim;
  this->yDim = nyDim;
  this->zDim = nzDim;

  this->n4D = num4D;

  // Do the allocation
  const size_t nVoxels = this->xDim * this->yDim * this->zDim;

  // The offset array
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_offsets4D),
                              (nVoxels+1)*sizeof(size_t) ) );

  // Space for maxLabels
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_maxLabels),
                              nVoxels*sizeof(short) ) );

  // Space for the labels
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_labels),
                              this->n4D*sizeof(unsigned short) ) );

  // Space for the priors
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_priors),
                              this->n4D*sizeof(float) ) );

  // Space for the total_training
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_totTraining),
                              nVoxels*sizeof(int) ) );
}


// --------------------

void GCApriorGPU::Release( void )
{

  if( this->xDim != 0 )
  {
    // Release offset array
    CUDA_SAFE_CALL( cudaFree( this->d_offsets4D ) );
    this->d_offsets4D = NULL;

    // Release 3D arrays
    CUDA_SAFE_CALL( cudaFree( this->d_maxLabels ) );
    this->d_maxLabels = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_totTraining ) );
    this->d_totTraining = NULL;

    // Release 4D arrays
    CUDA_SAFE_CALL( cudaFree( this->d_labels ) );
    this->d_labels = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_priors ) );
    this->d_priors = NULL;

    // Zero sizes
    this->xDim = 0;
    this->yDim = 0;
    this->zDim = 0;
    this->n4D = 0;

  }
}


// --------------------

void GCApriorGPU::Send( const Freesurfer::GCAlinearPrior& src )
{

  // Allocate memory
  this->Allocate( src.xDim, src.yDim, src.zDim, src.n4D );

  const size_t nVoxels = this->xDim * this->yDim * this->zDim;

  // Copy offsets array
  CUDA_SAFE_CALL( cudaMemcpy( this->d_offsets4D,
                              &src.offsets4D.front(),
                              (nVoxels+1)*sizeof(size_t),
                              cudaMemcpyHostToDevice ) );

  // Copy 3D arrays
  CUDA_SAFE_CALL( cudaMemcpy( this->d_maxLabels,
                              &src.maxLabels.front(),
                              nVoxels*sizeof(short),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_totTraining,
                              &src.totTraining.front(),
                              nVoxels*sizeof(int),
                              cudaMemcpyHostToDevice ) );

  // Copy 4D arrays
  CUDA_SAFE_CALL( cudaMemcpy( this->d_labels,
                              &src.labels.front(),
                              this->n4D*sizeof(unsigned short),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_priors,
                              &src.priors.front(),
                              this->n4D*sizeof(float),
                              cudaMemcpyHostToDevice ) );
}

}
}
