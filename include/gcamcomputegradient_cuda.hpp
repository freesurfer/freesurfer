/**
 * @file  gcamcomputegradient_cuda.hpp
 * @brief Header file for gcamComputeGradient for the GPU
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
 *    $Date: 2010/11/30 16:25:59 $
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

#ifndef GCAM_COMPUTE_GRADIENT_CUDA_HPP
#define GCAM_COMPUTE_GRADIENT_CUDA_HPP

#ifdef GCAMORPH_ON_GPU

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"



template<typename T, typename U>
int ComputeGradient( GPU::Classes::GCAmorphGPU& gcam,
		     const GPU::Classes::MRIframeGPU<T>& mri,
		     const GPU::Classes::MRIframeGPU<U>& mri_smooth,
		     GCA_MORPH_PARMS *parms );


#endif

#endif
