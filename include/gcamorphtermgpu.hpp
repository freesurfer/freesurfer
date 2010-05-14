/**
 * @file  gcamorphtermgpu.hpp
 * @brief Holds routines to compute GCAmorph terms on the GPU (header)
 *
 * This file holds a variety of routines which compute terms for
 * mri_ca_register
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/05/14 16:37:04 $
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


#ifndef GCA_MORPH_TERM_GPU_HPP
#define GCA_MORPH_TERM_GPU_HPP

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"

namespace GPU {
  namespace Algorithms {

    //! Class to hold GCAMorph 'term' computations
    class GCAmorphTerm {

    public:

      //! Constructor (stream defaults to zero)
      GCAmorphTerm( const cudaStream_t s = 0 ) : stream(s) {};

      //! Destructor
      ~GCAmorphTerm( void ) {};

      //! Routine to print timing information to std.out
      static void ShowTimings( void );

      // ######################################################
    private:

      //! Stream to use for operations
      cudaStream_t stream;

    };

  }
}

#endif
