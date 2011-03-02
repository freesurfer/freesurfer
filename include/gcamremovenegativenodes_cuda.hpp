/**
 * @file  gcamremovenegativenodes_cuda.hpp
 * @brief Implementation of gcamRemoveNegativeNodes for the GPU
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
 *    $Date: 2011/03/02 18:36:57 $
 *    $Revision: 1.2 $
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


#ifndef GCAM_REMOVE_NEGATIVE_NODES_HPP
#define GCAM_REMOVE_NEGATIVE_NODES_HPP

#include "gcamorphgpu.hpp"
#include "mriframegpu.hpp"

#include "gcamfots_cuda.hpp"
#include "gcamorphenergy.hpp"
#include "gcamcomputegradient_cuda.hpp"


const int kMaxNegIter = 100;

//! Implement gcamRemoveNegativeNodes
template<typename T>
void RemoveNegativeNodes( GPU::Classes::GCAmorphGPU& gcam,
                          const GPU::Classes::MRIframeGPU<T>& mri,
                          GCA_MORPH_PARMS *parms ) {
  /*!
    Implementation of gcamRemoveNegativeNodes for the GPU
    pipeline
  */
  GCA_MORPH_PARMS saved_parms = *parms;
  double min_dt, orig_dt = parms->orig_dt, rms, last_rms, pct_change;
  int old_neg, new_neg, i;
  GPU::Algorithms::GCAmorphEnergy gcamEnergy;

  if( gcam.neg <=0 ) {
    return;
  }

  // Try simple removal
  gcam.RemoveSingularities();
  if( gcam.neg <= 0 ) {
    return;
  }

  parms->noneg = 0 ;
  parms->l_distance
    = parms->l_log_likelihood
    = parms->l_binary
    = parms->l_multiscale
    = parms->l_spring
    = parms->l_area
    = parms->l_smoothness
    = parms->l_label = 0;
  parms->navgs = 0;
  parms->l_area = 0.0;
  parms->l_jacobian = 1;
  parms->dt = 0.1;
  parms->tol = 0.01;

  last_rms = rms = gcamEnergy.ComputeRMS( gcam, mri, parms );
  new_neg = gcam.neg;
  i = 0;

  do {
    old_neg = new_neg;

    ComputeGradient( gcam, mri, mri, parms );

    min_dt = FindOptimalTimestep( gcam, parms, mri );

    parms->dt = min_dt;

    gcam.ApplyGradient( parms );

    last_rms = rms;
    rms = gcamEnergy.ComputeRMS( gcam, mri, parms );

    parms->dt = orig_dt;
    new_neg = gcam.neg;
    pct_change = 100.0*(last_rms-rms)/(last_rms);

    printf( "iter %d, dt=%2.6f: new neg %d, old_neg %d, delta %d, rms=%2.3f\n",
            ++i, min_dt, new_neg, old_neg, old_neg-new_neg, rms, pct_change );
  }
  while( (new_neg>0) && (pct_change > parms->tol) && (i < kMaxNegIter) );

  *parms = saved_parms;
}


#endif
