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
 *    $Date: 2010/03/05 18:05:50 $
 *    $Revision: 1.10 $
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

#include "gcamorph.h"

#include "volumegpu.hpp"


// ==================================================================

namespace GPU {

  namespace Classes {

    //! Class to hold a GCA morph on the GPU
    class GCAmorphGPU {
      /*!
	Class to hold a GCA morph on the GPU.
	This version of the class only supports
	one 'input' in the GC1D structure
      */
    public:
      //! Matches x, y and z in GCAmorph
      VolumeGPU<float3> d_r;
      //! Matches invalid flag in GCAmorph
      VolumeGPU<unsigned char> d_invalid;
      //! Matches orig_area field in GCAmorph
      VolumeGPU<float> d_origArea;
      //! Matches area field in GCAmorph
      VolumeGPU<float> d_area;
      //! Matches area1 field in GCAmorph
      VolumeGPU<float> d_area1;
      //! Matches area2 field in GCAmorph
      VolumeGPU<float> d_area2;
      //! Matches label field in GCAMorph
      VolumeGPU<int> d_label;
      //! Matches status field in GCAMorph
      VolumeGPU<int> d_status;
      //! Matches the 'means' field of the GC1D
      VolumeGPU<float> d_mean;
      //! Matches the 'covars' field of the GC1D (a variance with only one mean). A negative value indicates that no value is stored for this or corresponding d_mean
      VolumeGPU<float> d_variance;

      // -----------------------------------------
      // Constructors & Destructor

      //! Default constructor
      GCAmorphGPU( void ) : d_r(),
			    d_invalid(),
			    d_area(),
			    d_origArea(),
			    d_area1(),
			    d_area2(),
			    d_label(),
			    d_status(),
			    d_mean(),
			    d_variance() {};

      //! Destructor
      ~GCAmorphGPU( void ) {};

      // -------------------------------------------

      //! Checks integrity of members
      void CheckIntegrity( void ) const;

      // -------------------------------------------
      // Memory management
      
      //! Allocates GPU memory for volume of given size
      void AllocateAll( const dim3& dims );

      //! Releases all the GPU memory
      void ReleaseAll( void );

      // -------------------------------------------
      // Transfer routines

      //! Sends all data to the GPU
      void SendAll( const GCAM* src );

      //! Receives all data from the GPU
      void RecvAll( GCAM* dst ) const;

      // -------------------------------------------
      // Computations

      //! Computes the properties of the metric
      void ComputeMetricProperties( int& invalid, int& neg );

    private:

    };

  }
}

#endif
