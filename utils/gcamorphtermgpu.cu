/**
 * @file  gcamorphtermgpu.cu
 * @brief Holds routines to compute GCAmorph terms on the GPU
 *
 * This file holds a variety of routines which compute terms for
 * mri_ca_register
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/05/14 16:37:03 $
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

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"

#include "gcamorphtermgpu.hpp"



// ==============================================================


namespace GPU {
  namespace Algorithms {


    // ##############################################################
    void GCAmorphTerm::ShowTimings( void ) {

#ifdef CUDA_SHOW_TIMINGS
      std::cout << "==================================" << std::endl;
      std::cout << "GPU GCAmorph term timers" << std::endl;
      std::cout << "------------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
      std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
      std::cout << "Timings may not be accurate" << std::endl;
#endif
      std::cout << std::endl;

     std::cout << "==================================" << std::endl;
#endif
    }
  }
}
