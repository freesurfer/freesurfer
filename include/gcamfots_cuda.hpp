/**
 * @file  gcamfots_cuda.hpp
 * @brief C++ interface to gcamFindOptimalTimeStep for the GPU
 *
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/11/30 16:25:59 $
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
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
