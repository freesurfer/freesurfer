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
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.2 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
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
