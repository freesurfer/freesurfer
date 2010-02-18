/**
 * @file  gcamorphgpu.hpp
 * @brief Holds GCA morph data on the GPU
 *
 * 
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/18 18:42:38 $
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
 *
 */


#ifndef GCA_MORPH_GPU_H
#define GCA_MORPH_GPU_H

#include "volumegpu.hpp"


// ==================================================================

namespace GPU {

  namespace Classes {

    //! Class to hold a GCA morph on the GPU
    class GCAmorphGPU {
    public:
      //! Matches x, y and z in GCAmorph
      VolumeArgGPU<float3> d_r;
      //! Matches invalid flag in GCAmorph
      VolumeArgGPU<unsigned char> d_invalid;
      //! Matches area field in GCAmorph
      VolumeArgGPU<float> d_area;
      //! Matches area1 field in GCAmorph
      VolumeArgGPU<float> d_area1;
      //! Matches area2 field in GCAmorph
      VolumeArgGPU<float> d_area2;

      // -----------------------------------------
      

    private:

    };

  }
}



#endif
