/**
 * @file  gcasgpu.cu
 * @brief Holds GCAS class for the GPU
 *
 *
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.7 $
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

using namespace std;

#include "cudacheck.h"

#include "gcasgpu.hpp"



// ====================================================

namespace GPU
{

namespace Classes
{


void GCASampleGPU::SendGPU( GCA *gca,
                            const GCA_SAMPLE *gcaSample,
                            const MRI *mri_inputs,
                            const int nSamples )
{
  /*!
  Sends the given GCAS to the GPU, doing
  various other things whose use is currently
  obscure.
  */

  // Verify the number of inputs
  if( gca->ninputs != 1 )
  {
    cerr << __FUNCTION__
         << ": Can only have one input"
         << endl;
    exit( EXIT_FAILURE );
  }

  // Allocate device memory
  this->Allocate( nSamples );

  // Identity transform matrix
  MATRIX *identity = MatrixIdentity( 4, NULL );
  TRANSFORM *identityTransform = TransformAlloc( LINEAR_VOX_TO_VOX, NULL );
  static_cast<LTA*>(identityTransform->xform)->xforms[0].m_L = identity;

  // Allocate some arrays
  int* myx = new int[nSamples];
  int* myy = new int[nSamples];
  int* myz = new int[nSamples];

  float* covars = new float[nSamples];
  float* priors = new float[nSamples];
  float* means = new float[nSamples];

  for( int i=0; i<nSamples; i++ )
  {
    // Copy code from GCAcomputeLogSampleProbability() in gca.c
    int xp = gcaSample[i].xp;
    int yp = gcaSample[i].yp;
    int zp = gcaSample[i].zp;

    if( GCApriorToSourceVoxel( gca, mri_inputs, identityTransform,
                               xp, yp, zp,
                               &myx[i], &myy[i], &myz[i] ) != NO_ERROR )
    {
      cerr << __FUNCTION__ << ": Failed with i=" << i << endl;
      exit( EXIT_FAILURE );
    }

    // These lines require the check for ninputs==1
    covars[i] = gcaSample[i].covars[0];
    means[i] = gcaSample[i].means[0];
    priors[i] = gcaSample[i].prior;

  }

  // Send to the GPU
  CUDA_SAFE_CALL( cudaMemcpy( this->d_x, myx,
                              nSamples*sizeof(int),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_y, myy,
                              nSamples*sizeof(int),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_z, myz,
                              nSamples*sizeof(int),
                              cudaMemcpyHostToDevice ) );

  CUDA_SAFE_CALL( cudaMemcpy( this->d_means, means,
                              nSamples*sizeof(float),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_covars, covars,
                              nSamples*sizeof(float),
                              cudaMemcpyHostToDevice ) );
  CUDA_SAFE_CALL( cudaMemcpy( this->d_priors, priors,
                              nSamples*sizeof(float),
                              cudaMemcpyHostToDevice ) );

  // Release memory
  delete[] myx;
  delete[] myy;
  delete[] myz;
  delete[] covars;
  delete[] priors;
  delete[] means;

  // Following should also free the identity matrix
  TransformFree( &identityTransform );
}

// ##################################################

// Memory management

void GCASampleGPU::Allocate( const unsigned int n )
{
  if( this->nSamplesAlloc < n )
  {
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

void GCASampleGPU::Release( void )
{
  if( this->nSamplesAlloc != 0 )
  {
    cudaFree( d_x );
    cudaFree( d_y );
    cudaFree( d_z );
    cudaFree( d_means );
    cudaFree( d_priors );
    cudaFree( d_covars );
    this->nSamples = 0;
    this->nSamplesAlloc = 0;
  }
}

}
}
