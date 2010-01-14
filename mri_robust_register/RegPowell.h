/**
 * @file RegPowell.h
 * @brief A class to compute a robust registration using Powell
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/01/14 19:41:04 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2008-2009
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

//
// written by Martin Reuter
// Apr. 22th ,2009
//

#ifndef RegPowell_H
#define RegPowell_H

#ifdef __cplusplus
extern "C"
{
#endif
#include "matrix.h"
#include "mri.h"
#ifdef __cplusplus
}
#endif

#include <utility>
#include <string>
#include <vector>
#include "Registration.h"


class RegPowell : public Registration
{
public:
  RegPowell():Registration()
  {};
  RegPowell(MRI * s, MRI *t):Registration(s,t)
  {};
  virtual ~RegPowell()
  {};

  virtual void computeIterativeRegistration( int n,double epsit,MRI * mriS=NULL, MRI* mriT=NULL, const vnl_matrix < double > &Minit = vnl_matrix<double>(), double iscaleinit = 1.0);

protected:

  static float costFunction(float p[] );
  static RegPowell* tocurrent;
  static MRI * scf;
  static MRI * tcf;
  static int pcount;
  static vnl_matrix_fixed < double , 4, 4 > mh1;
  static vnl_matrix_fixed < double , 4, 4 > mh2;
  static int icount;

};


#endif

