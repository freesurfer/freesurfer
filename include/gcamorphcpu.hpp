/**
 * @file  gcamorphcpu.hpp
 * @brief Holds GCA morph data on the CPU
 *
 * Mirrors GCAmorphGPU on the CPU. Yes this is messy.....
 * 
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/09/28 19:40:22 $
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

#ifdef GCAMORPH_ON_GPU

#ifndef GCA_MORPH_GPU_HPP
#define GCA_MORPH_GPU_HPP

#include <cuda_runtime.h>
#include "gcamorphgpu.hpp"

#include "volumecpu.hpp"

// =================================

namespace Freesurfer {

  class GCAmorphCPU {
  public:
    // Member variables
    VolumeCPU<float> rx, ry, rz;
    VolumeCPU<float> origx, origy, origz;
    VolumeCPU<float> dx, dy, dz;
    VolumeCPU<float> odx, ody, odz;

    VolumeCPU<float> origArea, origArea1, origArea2;
    VolumeCPU<float> area, area1, area2;

    VolumeCPU<char> invalid;
    VolumeCPU<int> label;
    VolumeCPU<int> status;
    VolumeCPU<float> labelDist;

    VolumeCPU<float> mean, variance;

    double exp_k;
    int neg;

    // -----------------------------------
    //! Default constructor
    GCAmorphCPU( void ) : rx(), ry(), rz(),
			  origx(), origy(), origz(),
			  dx(), dy(), dz(),
			  odx(), ody(), odz(),
			  origArea(), origArea1(), origArea2(),
			  area(), area1(), area2(),
			  invalid(),
			  label(),
			  status(),
			  labelDist(),
			  mean(),
			  variance(),
			  exp_k(0),
			  neg(0) {};

    //! Destructor
    ~GCAmorphCPU( void ) {};


    // ----------------------------------

    //! Routine to verify sizes
    void CheckIntegrity( void ) const;


    // ----------------------------------
    
    void AllocateAll( const unsigned int nx,
		      const unsigned int ny,
		      const unsigned int nz );

    // ----------------------------------

    //! Get data from the GPU
    void GetFromGPU( const GPU::Classes::GCAmorphGPU& src );
    
    //! Return data to the GPU
    void PutOnGPU( GPU::Classes::GCAmorphGPU& dst ) const;
    
    
  };

}

#endif




#endif
