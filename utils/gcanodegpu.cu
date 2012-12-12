/**
 * @file  gcanodegpu.cu
 * @brief Class to hold a volume of GCA nodes in linear memory on the GPU
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.4 $
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

#include "gcanodegpu.hpp"



namespace GPU
{
namespace Classes
{


// ----------------------------------------

void GCAnodeGPU::Allocate( const int nxDim, const int nyDim, const int nzDim,
                           const size_t num4D, const size_t num5D, const size_t num6D )
{

  // Get rid of old allocation
  this->Release();

  // Copy sizes
  this->xDim = nxDim;
  this->yDim = nyDim;
  this->zDim = nzDim;

  this->n4D = num4D;
  this->n5D = num5D;
  this->n6D = num6D;

  // Do the allocations
  const size_t nVoxels = this->xDim * this->yDim * this->zDim;

  // Offset arrays
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_offsets4D),
                              (nVoxels+1)*sizeof(size_t) ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_offsets5D),
                              (this->n4D+1)*sizeof(size_t) ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_offsets6D),
                              (this->n5D+1)*sizeof(size_t) ) );

  // 3D arrays
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_nodeMaxLabels),
                              nVoxels*sizeof(int) ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_nodeTotalTraining),
                              nVoxels*sizeof(int) ) );

  // 4D arrays
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_nodeLabels),
                              this->n4D * sizeof(unsigned short) ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_means),
                              this->n4D * sizeof(float) ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_variances),
                              this->n4D * sizeof(float) ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_nJustPriors),
                              this->n4D * sizeof(short) ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_nTraining),
                              this->n4D * sizeof(int) ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_regularised),
                              this->n4D * sizeof(char) ) );

  // 6D arrays
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_gc1dDirecLabelPriors),
                              this->n6D * sizeof(float) ) );
  CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_gc1dDirecLabels),
                              this->n6D * sizeof(unsigned short) ) );
}


// -----------

void GCAnodeGPU::Release( void )
{
  if( this->xDim != 0 )
  {
    // Release offset arrays
    CUDA_SAFE_CALL( cudaFree( this->d_offsets4D ) );
    this->d_offsets4D = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_offsets5D ) );
    this->d_offsets5D = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_offsets6D ) );
    this->d_offsets6D = NULL;

    // Release 3D arrays
    CUDA_SAFE_CALL( cudaFree( this->d_nodeMaxLabels ) );
    this->d_nodeMaxLabels = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_nodeTotalTraining ) );
    this->d_nodeTotalTraining = NULL;

    // Release 4D arrays
    CUDA_SAFE_CALL( cudaFree( this->d_nodeLabels ) );
    this->d_nodeLabels = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_means ) );
    this->d_means = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_variances ) );
    this->d_variances = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_nJustPriors ) );
    this->d_nJustPriors = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_nTraining ) );
    this->d_nTraining = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_regularised ) );
    this->d_regularised = NULL;

    // Release 6D arrays
    CUDA_SAFE_CALL( cudaFree( this->d_gc1dDirecLabelPriors ) );
    this->d_gc1dDirecLabelPriors = NULL;
    CUDA_SAFE_CALL( cudaFree( this->d_gc1dDirecLabels ) );
    this->d_gc1dDirecLabels = NULL;

    // Zero sizes
    this->xDim = 0;
    this->yDim = 0;
    this->zDim = 0;
    this->n4D = 0;
    this->n5D = 0;
    this->n6D = 0;
  }
}


// ----------------------------------------

void GCAnodeGPU::Send( const Freesurfer::GCAlinearNode& src )
{

  // Allocate memory
  this->Allocate( src.xDim, src.yDim, src.zDim,
                  src.n4D, src.n5D, src.n6D );

  const size_t nVoxels = this->xDim * this->yDim * this->zDim;

  // Copy offset arrays
  CUDA_SAFE_CALL( cudaMemcpy( this->d_offsets4D,
                              &src.offsets4D.front(),
                              (nVoxels+1)*sizeof(size_t),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_offsets5D,
                              &src.offsets5D.front(),
                              (this->n4D+1)*sizeof(size_t),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_offsets6D,
                              &src.offsets6D.front(),
                              (this->n5D+1)*sizeof(size_t),
                              cudaMemcpyHostToDevice ) );

  // Copy 3D arrays
  CUDA_SAFE_CALL( cudaMemcpy( this->d_nodeMaxLabels,
                              &src.nodeMaxLabels.front(),
                              nVoxels*sizeof(int),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_nodeTotalTraining,
                              &src.nodeTotalTraining.front(),
                              nVoxels*sizeof(int),
                              cudaMemcpyHostToDevice ) );

  // Copy 4D arrays
  CUDA_SAFE_CALL( cudaMemcpy( this->d_nodeLabels,
                              &src.nodeLabels.front(),
                              this->n4D*sizeof(unsigned short),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_means,
                              &src.means.front(),
                              this->n4D*sizeof(float),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_variances,
                              &src.variances.front(),
                              this->n4D*sizeof(float),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_nJustPriors,
                              &src.nJustPriors.front(),
                              this->n4D*sizeof(short),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_nTraining,
                              &src.nTraining.front(),
                              this->n4D*sizeof(int),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_regularised,
                              &src.regularised.front(),
                              this->n4D*sizeof(char),
                              cudaMemcpyHostToDevice ) );

  // Copy 6D arrays
  CUDA_SAFE_CALL( cudaMemcpy( this->d_gc1dDirecLabelPriors,
                              &src.gc1dDirecLabelPriors.front(),
                              this->n6D*sizeof(float),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_gc1dDirecLabels,
                              &src.gc1dDirecLabels.front(),
                              this->n6D*sizeof(unsigned short),
                              cudaMemcpyHostToDevice ) );
}

}
}
