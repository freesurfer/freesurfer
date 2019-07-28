/**
 * This is taken from vnl_powell.  This version is created so that the initial
 * direction, xi, can be set.
 */

// This is core/vnl/algo/vnl_powell.h
#ifndef fs_powell_h_
#define fs_powell_h_
#ifdef VCL_NEEDS_PRAGMA_INTERFACE
#pragma interface
#endif
//:
// \file
// \brief Powell minimizer.
// \author awf@robots.ox.ac.uk
// \date   05 Dec 00

#define export // obsolete feature "export template" used in these header files
#include <vnl/vnl_nonlinear_minimizer.h>
#undef export
#include "fs_vnl/fs_cost_function.h"

//: The ever-popular Powell minimizer.
// Derivative-free method which may be faster if your
// function is expensive to compute and many-dimensional.
// Implemented from scratch from NR.
class fs_powell : public vnl_nonlinear_minimizer
{
public:

  //: Initialize a powell with the given cost function
  fs_powell(fs_cost_function* functor)
      : functor_(functor), linmin_xtol_(1e-4), initial_step_(1.0)
  {}

  //: Run minimization, place result in x.
  ReturnCodes minimize(vnl_vector<double>& x, vnl_matrix<double>* xi);

  //: Set tolerance on line search parameter step
  //  Default value is 0.0001
  void set_linmin_xtol(double tol)
  {
    linmin_xtol_ = tol;
  }

  //: Set initial step when bracketting minima along a line
  //  Default value is 1.0
  void set_initial_step(double step)
  {
    initial_step_ = step;
  }

protected:
  fs_cost_function* functor_;

  friend class fs_powell_1dfun;
  void pub_report_eval(double e)
  {
    report_eval(e);
  }

  //: Tolerance on line search parameter step
  double linmin_xtol_;

  //: Initial step when bracketting minima along a line
  double initial_step_;
};

#endif // fs_powell_h_
