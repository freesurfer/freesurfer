/**
 * @file  gcamregisterlevel_cuda.cu
 * @brief Implementation of GCAMregisterLevel for the GPU
 *
 * Reference:
  * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
  * Structures in the Human Brain", Fischl et al.
  * (2002) Neuron, 33:341-355.
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/12/03 18:48:45 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2010,
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

#ifndef GCAM_REGISTER_LEVEL_CUDA_HPP
#define GCAM_REGISTER_LEVEL_CUDA_HPP


#ifdef GCAMORPH_ON_GPU

#include "macros.h"
#include "error.h"

#include "gcamorph.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"
#include "gcamorphenergy.hpp"

#include "gcamcomputegradient_cuda.hpp"
#include "gcamfots_cuda.hpp"


// ========================================================================

template<typename T, typename U>
int RegisterLevel( GPU::Classes::GCAmorphGPU& gcam,
		   const GPU::Classes::MRIframeGPU<T>& mri,
		   const GPU::Classes::MRIframeGPU<U>& mri_smooth,
		   GCA_MORPH_PARMS *parms );

#endif

#endif
