/**
 * @file  gcamfots_cuda.hpp
 * @brief C++ interface to gcamFindOptimalTimeStep for the GPU
 *
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

#ifndef GCA_MORPH_FOTS_HPP
#define GCA_MORPH_FOTS_HPP

#include "gcamorphgpu.hpp"
#include "mriframegpu.hpp"

template<typename T>
float FindOptimalTimestep( GPU::Classes::GCAmorphGPU& gcam,
                           GCA_MORPH_PARMS *parms,
                           const GPU::Classes::MRIframeGPU<T>& mri );

#endif
