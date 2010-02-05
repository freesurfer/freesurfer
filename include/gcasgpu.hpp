/**
 * @file  gcasgpu.hpp
 * @brief Holds GCAS class for the GPU
 *
 * 
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/05 17:25:28 $
 *    $Revision: 1.3 $
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

#ifndef GCAS_GPU_CUDA_H
#define GCAS_GPU_CUDA_H

#include <cstdlib>
#include <iostream>

#include "error.h"
#include "mri.h"
#include "matrix.h"
#include "gca.h"

// ====================================================

namespace GPU {

  namespace Classes {

    // Forward declaration on class passed to kernel
    class GCASonGPUonemean;

    //! Class to hold GCAS with one mean per location
    class GCASonemeanGPU {
    public:
      // -----------------------------
      //! Default constructor
      GCASonemeanGPU( void ) : nSamples(0), nSamplesAlloc(0),
			       d_x(NULL), d_y(NULL), d_z(NULL),
			       d_means(NULL),
			       d_priors(NULL),
			       d_covars(NULL) {}

      // Destructor
      ~GCASonemeanGPU( void ) {
	this->Release();
      }

      // -----------------------------

      //! Sends GCAS info to the GPU 
      void SendGPU( GCA *gca,
		    GCA_SAMPLE *gcaSample,
		    MRI *mri_inputs,
		    const int nSamples );

      // Declase related 'kernel' class a friend for private access
      friend class GCASonGPUonemean;

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
      GCASonemeanGPU( const GCASonemeanGPU& src ) : nSamples(0),
						    nSamplesAlloc(0),
						    d_x(NULL),
						    d_y(NULL),
						    d_z(NULL),
						    d_means(NULL),
						    d_priors(NULL),
						    d_covars(NULL) {
	std::cerr << __FUNCTION__ << ": Please don't copy"
		  << std::endl;
	exit( EXIT_FAILURE );
      }

      //! Assignment operator aborts
      GCASonemeanGPU& operator=( const GCASonemeanGPU& src ) {
	std::cerr << __FUNCTION__ << ": Please don't copy"
		  << std::endl;
	exit( EXIT_FAILURE );
      }
    };


    class GCASonGPUonemean {
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
      
      //! Construct from GCASonemeanGPU
      GCASonGPUonemean( const GCASonemeanGPU& src ) : nSamples(src.nSamples),
						      x(src.d_x),
						      y(src.d_y),
						      z(src.d_z),
						      means(src.d_means),
						      priors(src.d_priors),
						      covars(src.d_covars) {}

    private:
      // ================================================

      //! Default constructor should not be used
      GCASonGPUonemean( void ) : nSamples(0),
				 x(NULL), y(NULL), z(NULL),
				 means(NULL), priors(NULL), covars(NULL) {
	std::cerr << __FUNCTION__ << ": Please don't use"
		  << std::endl;
	exit( EXIT_FAILURE );
      }

      //! Assignment operator should not be used
      GCASonGPUonemean& operator=( const GCASonGPUonemean& src ) {
	std::cerr << __FUNCTION__ << ": Please don't use"
		  << std::endl;
	exit( EXIT_FAILURE );
      }
      
      
    };


  }

}




#endif
