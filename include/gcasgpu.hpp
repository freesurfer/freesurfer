/**
 * @file  gcasgpu.hpp
 * @brief Holds GCAS class for the GPU
 *
 *
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.8 $
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

#ifndef GCAS_GPU_CUDA_H
#define GCAS_GPU_CUDA_H

#include <cstdlib>
#include <iostream>

#include "error.h"
#include "mri.h"
#include "matrix.h"
#include "gca.h"

// ====================================================

namespace GPU
{

namespace Classes
{

// Forward declaration on class passed to kernel
class GCASonGPU;

//! Class to hold GCAS with one mean per location
class GCASampleGPU
{
public:
  // -----------------------------
  //! Default constructor
  GCASampleGPU( void ) : nSamples(0), nSamplesAlloc(0),
    d_x(NULL), d_y(NULL), d_z(NULL),
    d_means(NULL),
    d_priors(NULL),
    d_covars(NULL) {}

  // Destructor
  ~GCASampleGPU( void )
  {
    this->Release();
  }

  // -----------------------------

  //! Sends GCAS info to the GPU
  void SendGPU( GCA *gca,
                const GCA_SAMPLE *gcaSample,
                const MRI *mri_inputs,
                const int nSamples );

  // Declase related 'kernel' class a friend for private access
  friend class GCASonGPU;

  // =================================
private:
  //! Number of samples in this GCAS
  unsigned int nSamples;
  //! Space allocated for samples
  unsigned int nSamplesAlloc;
  //! List of x positions
  int *d_x;
  //! List of y positions
  int *d_y;
  //! List of z positions
  int *d_z;

  //! List of means
  float *d_means;
  //! List of priors
  float *d_priors;
  //! List of covariances
  float *d_covars;

  // -----------------------------------
  // Memory management

  //! Allocates space for n samples
  void Allocate( const unsigned int n );
  //! Releases GPU memory
  void Release( void );

  // ------------------------------------
  // Inhibit copying

  //! Copy constructor aborts
  GCASampleGPU( const GCASampleGPU& src ) : nSamples(0),
    nSamplesAlloc(0),
    d_x(NULL),
    d_y(NULL),
    d_z(NULL),
    d_means(NULL),
    d_priors(NULL),
    d_covars(NULL)
  {
    std::cerr << __FUNCTION__ << ": Please don't copy"
              << std::endl;
    abort();
  }

  //! Assignment operator aborts
  GCASampleGPU& operator=( const GCASampleGPU& src )
  {
    std::cerr << __FUNCTION__ << ": Please don't copy"
              << std::endl;
    abort();
  }
};


class GCASonGPU
{
public:
  //! Number of samples
  unsigned int nSamples;

  //! List of x positions
  int *x;
  //! List of y positions
  int *y;
  //! List of z positions
  int *z;

  //! List of means
  float *means;
  //! List of priors
  float *priors;
  //! List of covariances
  float *covars;

  //! Construct from GCASGPU
  GCASonGPU( const GCASampleGPU& src ) : nSamples(src.nSamples),
    x(src.d_x),
    y(src.d_y),
    z(src.d_z),
    means(src.d_means),
    priors(src.d_priors),
    covars(src.d_covars) {}

  //! Accessor for locations
  __device__ float3 GetLocation( const unsigned int i ) const
  {
    float3 r;
    r.x = this->x[i];
    r.y = this->y[i];
    r.z = this->z[i];

    return( r );
  }

private:
  // ================================================

  //! Default constructor should not be used
  GCASonGPU( void ) : nSamples(0),
    x(NULL), y(NULL), z(NULL),
    means(NULL), priors(NULL), covars(NULL)
  {
    std::cerr << __FUNCTION__ << ": Please don't use"
              << std::endl;
    abort();
  }

  //! Assignment operator should not be used
  GCASonGPU& operator=( const GCASonGPU& src )
  {
    std::cerr << __FUNCTION__ << ": Please don't use"
              << std::endl;
    abort();
  }


};


}

}




#endif
