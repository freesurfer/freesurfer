/**
 * @file  svm-param.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:17 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2007,
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




