/**
 * @brief A class to compute a registration using Powell
 *
 */

/*
 * Original Author: Martin Reuter
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

//
// written by Martin Reuter
// Apr. 22th ,2009
//
#ifndef RegPowell_H
#define RegPowell_H

#include "matrix.h"
#include "mri.h"

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
      Registration(),xtol(1e-5),ftol(1e-5)
  {
  }

  virtual ~RegPowell()
  {
  }

  //! Set the tolerance for Powell minimizer (both x and f tollerance)
  void setTolerance(double tol) {xtol=tol;ftol=tol;}

  //! The Powell way of doing iterative registration
  virtual void computeIterativeRegistrationFull(int n, double epsit, MRI * mriS,
      MRI* mriT, const vnl_matrix<double> &Minit, double iscaleinit);
      
  //! The static cost function for the Powell minimizer
  static double costFunction(const vnl_vector<double> & p);

  //! Return the half way geometry (here use target for now)
  MRI * getHalfWayGeom()
  {
    return tcf;
  }

  //! Get Name of Registration class
  virtual std::string getClassName() {return "RegPowell";}

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
  double xtol;
  double ftol;

};

#endif

