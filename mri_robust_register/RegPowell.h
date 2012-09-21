/**
 * @file RegPowell.h
 * @brief A class to compute a registration using Powell
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2012/09/21 23:05:15 $
 *    $Revision: 1.10 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

/** \class RegPowell
 * \brief Class for registration using Powell method
 */
class RegPowell: public Registration
{
public:
  RegPowell() :
      Registration()
  {
  }
  ;
  virtual ~RegPowell()
  {
  }
  ;

//  void computeIterativeRegistration(int n,double epsit){Registration::computeIterativeRegistration(n,epsit);};

  virtual void computeIterativeRegistration(int n, double epsit, MRI * mriS,
      MRI* mriT, const vnl_matrix<double> &Minit, double iscaleinit);
  //! The static cost function for the Powell minimizer
  static double costFunction(const vnl_vector<double> & p);

  MRI * getHalfWayGeom()
  {
    return tcf;
  }
  ;

protected:

  virtual void setTransformation(bool is2d);
  static RegPowell* tocurrent;
  static MRI * scf;
  static MRI * tcf;
  static int pcount;
  static vnl_matrix_fixed<double, 4, 4> mh1;
  static vnl_matrix_fixed<double, 4, 4> mh2;
  static int icount;
  static int subsamp;
  static bool is2d;

};

#endif

