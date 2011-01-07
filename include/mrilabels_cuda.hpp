/**
 * @file  mrilabels_cuda.hpp
 * @brief Holds various MRI label routines for the GPU
 *
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2011/01/07 20:45:30 $
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
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef MRI_LABELS_CUDA_HPP
#define MRI_LABELS_CUDA_HPP

#include "mri.h"
#include "mriframegpu.hpp"


namespace GPU {
  namespace Algorithms {

    //! Class to perform various label things on MRIs
    class MRIlabels {
    public:

      //! Constructor
      MRIlabels( void ) {};

      //! Destructor
      ~MRIlabels( void ) {};

      static void ShowTimings( void );

      // --------------------------------------

      void MarkLabelBorderVoxels( const GPU::Classes::MRIframeGPU<unsigned char>& src,
				  GPU::Classes::MRIframeGPU<unsigned char>& dst,
				  const int label,
				  const int mark,
				  const int sixConnect ) const;

      float VoxInLabelWithPartialVolume( const GPU::Classes::MRIframeGPU<unsigned char>& mri,
					 const GPU::Classes::MRIframeGPU<unsigned char>& mri_vals,
					 const int label,
					 GPU::Classes::MRIframeGPU<float>& mri_mixing_coeff,
					 GPU::Classes::MRIframeGPU<unsigned char>& mri_nbr_labels ) const;
					 


      // --------------------------------------
    private:

      // Timers
      static SciGPU::Utilities::Chronometer tMarkLabelBorderVoxelsTot;

      static SciGPU::Utilities::Chronometer tVoxInLabelPartVolumeTot;
      static SciGPU::Utilities::Chronometer tVoxInLabelPartVolumeCompute;

      // ---------

      //! Suppress copy constructor
      MRIlabels( const MRIlabels& src ) {
	std::cerr << __FUNCTION__
		  << ": Please do not copy" << std::endl;
	abort();
      }

      //! Suppress assignment operator
      MRIlabels& operator=( const MRIlabels& src ) {
	std::cerr << __FUNCTION__
		  << ": Please do not copy" << std::endl;
	abort();
      }

    };

  }
}


#endif
