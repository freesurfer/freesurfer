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
 *    $Date: 2010/02/05 14:41:32 $
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

#ifndef GCAS_GPU_CUDA_H
#define GCAS_GPU_CUDA_H




// ====================================================

namespace GPU {

  namespace Classes {


    //! Class to hold GCAS with one mean per location

    class GCASonemeanGPU {
    public:

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
    };



  }

}




#endif
