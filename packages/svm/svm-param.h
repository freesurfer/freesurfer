/*
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


////SVM-LIB////////////////////////////////////////////////////////////////
//
// Name: SvmParam
//
// This file defines all the parameters fro the svm optimization.
//
//  Polina Golland polina@ai.mit.edu
//
///////////////////////////////////////////////////////////////////////////



#ifndef __SVM_PARAM_H__
#define __SVM_PARAM_H__

#include "svm-kernel.h"

class SvmParam {
public:

  // SVM parameters

  double C;                              // soft margin constant
  static const double dC;
  static const char* nameC;

  Kernel kernel;                         // kernel
  static const Kernel dKernel;
  static const char* nameKernel;


  // Programming parameters

  int verbose;                           // verbose mode
  static const int dVerbose;
  static const char* nameVerbose ;

  double alphaEpsilon;                    // threshold to decide that
  static const double dAlphaEpsilon;      // alpha is zero
  static const char* nameAlphaEpsilon;

  double classEpsilon;                  // threshold for classification
  static const double dClassEpsilon;    // comparison to 0 and +/-1.
  static const char* nameClassEpsilon;


  // Interior point - specific optimization parameters

  int maxIterations;                     // max. number of iterations
  static const int dMaxIterations;
  static const char* nameMaxIterations;

  double sigDig;                          // the primal and the dual must
  static const double dSigDig;            // agree to that many digits
  static const char* nameSigDig;

  double optEpsilon;                     // initialization parameter
  static const double dOptEpsilon;
  static const char* nameOptEpsilon;


  // Constructors
  SvmParam();


  // I/O
  void parse(const char* const* argv, int argc);

  static std::ostream& printUsage(std::ostream& s);
};


#endif // __SVM_PARAM_H__




