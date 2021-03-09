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

#include "fs_vnl/fs_cost_function.h"

// TODO:
#include <iostream>

fs_cost_function::fs_cost_function(float (*function)(float[])) { SetFunction(function); }

fs_cost_function::fs_cost_function(float (*function)(float[]),
                                   void (*functionGradient)(float[], float[]),
                                   int numberOfUnknowns)
{
  SetFunction(function);
  SetFunctionGradient(functionGradient);
  set_number_of_unknowns(numberOfUnknowns);
}

fs_cost_function::~fs_cost_function()
{
  mFunction = NULL;
  mFunctionGradient = NULL;
}

void fs_cost_function::SetFunction(float (*function)(float[])) { mFunction = function; }

void fs_cost_function::SetFunctionGradient(void (*functionGradient)(float[], float[]))
{
  mFunctionGradient = functionGradient;
}

double fs_cost_function::f(const vnl_vector< double > &x)
{
  int numberOfParameters = x.size();

  float p[numberOfParameters + 1];

  // TODO: can't do x.copy_out( p ) because x is
  // of type double and p is float
  copyFromVNLToFloat(p, x, numberOfParameters);

  double value = static_cast< double >((*mFunction)(p));

  return value;
}

void fs_cost_function::gradf(vnl_vector< double > const &x, vnl_vector< double > &gradient)
{
  int numberOfParameters = static_cast< int >(x.size());

  float p[numberOfParameters + 1];
  float result[numberOfParameters + 1];

  copyFromVNLToFloat(p, x, numberOfParameters);

  (*mFunctionGradient)(p, result);

  for (int i = 0; i < numberOfParameters; i++) {
    gradient(i) = static_cast< double >(result[i + 1]);
  }
}

void fs_cost_function::copyFromVNLToFloat(float *floatVector,
                                          const vnl_vector< double > vnlVector,
                                          int numberOfParameters)
{
  for (int i = 0; i < numberOfParameters; i++) {
    // TODO: index starting at 1--legacy indexing from NR that the old
    //       functions passed in will have
    floatVector[i + 1] = static_cast< float >(vnlVector(i));
  }
}
