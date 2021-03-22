/*
 *    $Author$
 *    $Date$
 *    $Revision$
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <iostream>

#include "fs_vnl/fs_lbfgs_observer.h"

fs_lbfgs_observer::fs_lbfgs_observer()
{
  mNumberOfOptimalUpdates = 0;

  mStepFunction = NULL;
  mStepFunctionParms = NULL;
  mUserCallbackFunction = NULL;
}

fs_lbfgs_observer::~fs_lbfgs_observer() {}

void fs_lbfgs_observer::update(double bestF, vnl_vector< double > *bestX)
{
  if (hasStepFunction() || hasUserCallbackFunction()) {
    const int n = bestX->size();

    // legacy one indexing
    float currentX[n + 1];
    copyVnlToFloat(bestX, currentX, n);

    if (hasStepFunction()) {
      (*mStepFunction)(mNumberOfOptimalUpdates, static_cast< float >(bestF), mStepFunctionParms, currentX);
    }

    if (hasUserCallbackFunction()) {
      (*mUserCallbackFunction)(currentX);
    }
  }

  mNumberOfOptimalUpdates++;
}

bool fs_lbfgs_observer::hasStepFunction() { return (mStepFunction != NULL); }

bool fs_lbfgs_observer::hasUserCallbackFunction() { return (mUserCallbackFunction != NULL); }

int fs_lbfgs_observer::getNumberOfOptimalUpdates() { return mNumberOfOptimalUpdates; }

void fs_lbfgs_observer::setStepFunction(void (*stepFunction)(int itno, float sse, void *parms, float *p), void *parms)
{
  mStepFunction = stepFunction;
  mStepFunctionParms = parms;
}

void fs_lbfgs_observer::setUserCallbackFunction(void (*userCallbackFunction)(float[]))
{
  mUserCallbackFunction = userCallbackFunction;
}

void fs_lbfgs_observer::copyVnlToFloat(const vnl_vector< double > *input, float *output, const int n)
{
  for (int i = 0; i < n; i++) {
    // legacy one indexing
    output[i + 1] = static_cast< float >((*input)(i));
  }
}
