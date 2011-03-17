/**
 * @file  mrivol2vol_cuda.hpp
 * @brief Holds MRI vol2vol functions for the GPU. This is the class definition
 *
 * Implements MRI vol2vol on the GPU. These routines will be hooked
 * into the higher level routines in the main library.
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2011/03/17 17:18:55 $
 *    $Revision: 1.2 $
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


#ifndef MRI_VOL2VOL_CUDA_HPP
#define MRI_VOL2VOL_CUDA_HPP

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>


#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "affinegpu.hpp"

#include "ctfactory.hpp"

namespace GPU {

  namespace Algorithms {

    //! Class to perform volume-to-volume transforms on MRI data
    class MRIvol2vol {
    public:
      // ---------------------------------------

      MRIvol2vol( void ) {} ;

      ~MRIvol2vol( void ) {} ;

      static void ShowTimings( void );

      // ---------------------------------------

      void transform( const MRI* src, MRI* targ,
		      const MATRIX* transformMatrix,
		      const int InterpMode,
		      const float param );

      //! Dispatch function for the transformation kernel
      template<typename T, typename U>
      void transformGPU( const GPU::Classes::MRIframeGPU<T> &src,
			 GPU::Classes::MRIframeGPU<U> &dst,
			 const GPU::Classes::AffineTransformation& transform,
			 const int InterpMode,
			 const cudaStream_t myStream = 0 );
	

      // ------------------------------------------------------
    private:

      // Timers
      static SciGPU::Utilities::Chronometer tVol2VolMem, tVol2VolMemHost;
      static SciGPU::Utilities::Chronometer tVol2VolMRISendFrame, tVol2VolMRIRecvFrame;
      static SciGPU::Utilities::Chronometer tVol2VolCompute;
      static SciGPU::Utilities::Chronometer tVol2VolTotal;

      // --------------------------------------
      // Suppress copying

      MRIvol2vol( const MRIvol2vol& src ) {
	std::cerr << __FUNCTION__ << ": Please do not copy"
		  << std::endl;
	exit( EXIT_FAILURE );
      }

      MRIvol2vol& operator=( const MRIvol2vol& src ) {
	std::cerr << __FUNCTION__ << ": Please do not copy"
		  << std::endl;
	exit( EXIT_FAILURE );
      }

      // --------------------------------------------


      //! Templated texture binding wrapper
      template<typename T>
      GPU::Classes::CTfactory* BindSourceTexture( const GPU::Classes::MRIframeGPU<T>& src,
                                                  const int InterpMode ) {
	/*!
	  Will bind the CUDA array of the source frame to the appropriate
	  texture, with correct filtering set
	  Unspecialised version aborts
	*/
	std::cerr << __PRETTY_FUNCTION__ << ": Unrecognised type" << std::endl;
        abort();
        return(NULL);
      }

      //! Templated texture unbinding
      template<typename T>
      void UnbindSourceTexture( void ) {
	/*!
	  Will unbind the appropriate texture from the source frame.
	  Unspecialised version aborts
	*/
	std::cerr << __PRETTY_FUNCTION__ << ": Unrecognised type" << std::endl;
        abort();
      }

      // --------------------------------------------

      //! Routine to use GPU for vol2vol on all MRI frames
      template<typename T, typename U>
      void AllFramesGPU( const MRI* src, MRI* targ,
			 const MATRIX* transformMatrix,
			 const int InterpMode,
			 const float param );

      //! Dispatch routine to add templating
      template<typename T>
      void AllFramesDstDispatch( const MRI* src, MRI* targ,
				 const MATRIX* transformMatrix,
				 const int InterpMode,
				 const float param );
    };

  }

}



#endif
