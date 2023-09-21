/**
 * @brief This is taken from vnl_powell
 *
 * This version is created so that the initial
 * direction, xi, can be set in the minimize method.
 */
/*
 * Original Author: Dennis Jen
 *    $Author$
 *    $Date$f
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
 * General inquiries and bugs: freesurfer@nmr.mgh.harvard.edu
 *
 */

// This is core/vnl/algo/vnl_powell.cxx
#ifdef VCL_NEEDS_PRAGMA_INTERFACE
#pragma implementation
#endif
#include "fs_vnl/fs_powell.h"

#include <stdio.h>   // printf
#include <stdlib.h>  // calloc and free

#include <itkVersion.h>

#if ITK_VERSION_MAJOR >= 5
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vcl_legacy_aliases.h>
#else
#include <vcl_cassert.h>
#endif

#include <vnl/algo/vnl_brent.h>
#include <vnl/vnl_math.h>

#ifdef DEBUG
#include <vcl_iostream.h>
#include <vnl/vnl_matlab_print.h>
#endif

class fs_powell_1dfun : public vnl_cost_function
{
 public:
  fs_powell *powell_;
  vnl_cost_function *f_;
  unsigned int n_;
  vnl_vector< double > x0_;
  vnl_vector< double > dx_;
  vnl_vector< double > tmpx_;
  fs_powell_1dfun(int n, vnl_cost_function *f, fs_powell *p)
      : vnl_cost_function(1), powell_(p), f_(f), n_(n), x0_(n), dx_(n), tmpx_(n)
  {
  }

  void init(vnl_vector< double > const &x0, vnl_vector< double > const &dx)
  {
    x0_ = x0;
    dx_ = dx;
    assert(x0.size() == n_);
    assert(dx.size() == n_);
  }

  double f(const vnl_vector< double > &x)
  {
    uninit(x[0], tmpx_);
    double e = f_->f(tmpx_);
    powell_->pub_report_eval(e);
    return e;
  }

  void uninit(double lambda, vnl_vector< double > &out)
  {
    for (unsigned int i = 0; i < n_; ++i) out[i] = x0_[i] + lambda * dx_[i];
  }
};

/**
 * This minimizer is modified from the original vnl version to take an initial
 * direction.
 *
 * ftol is set by minimizer.set_f_tolerance(ftol)
 * linmin_xtol_ is set by minimizer.set_linmin_xtol(linmin_xtol_)
 * xtol is set by minimizer.set_x_tolerance(xtol)
 * verbose_ is set by minimizer.set_verbose(1);
 * maxfev is set by  minimizer.set_max_function_evals(MaxIterations);
 *
 * Note: each iteration is a loop through a 1D minimization for each
 * parameter.
 *
 * @param p Parameters to be minimized.  The minimized parameters will be
 * returned in p as well.
 * @param xi Initial direction.  If left out, then it will be set to the
 * identity initially.  The direction will be returned--so be sure
 * to deallocate it.
 *

 */
void (*powell_iteration_func)(float *p, int nparms) = NULL;

vnl_nonlinear_minimizer::ReturnCodes fs_powell::minimize(vnl_vector< double > &p, vnl_matrix< double > *xi)
// double p[], double **xi, int n
{
  // verbose_ = true;
  int n = p.size();
  fs_powell_1dfun f1d(n, functor_, this);

  // This line below was replaced by the if block
  //  vnl_matrix<double> xi(n,n, vnl_matrix_identity);
  if (xi == NULL) {
    xi = new vnl_matrix< double >(n, n, vnl_matrix_identity);
  }

  vnl_vector< double > ptt(n);
  vnl_vector< double > xit(n);
  double fret = functor_->f(p);
  report_eval(fret);
  vnl_vector< double > pt = p;
  printf("fs_powell::minimize\n");
  printf("  nparams %d\n", n);
  printf("  maxfev %ld\n", maxfev);
  printf("  ftol   %lf\n", (double)ftol);
  printf("  linmin_xtol_   %lf\n", (double)linmin_xtol_);

  while (num_iterations_ < unsigned(maxfev)) {
    double fp = fret;
    int ibig = 0;
    double del = 0.0;
    printf("  powell nthiter %d: fret = %f\n", (int)num_iterations_, fret);
    if (powell_iteration_func) {
      float *float_parms = (float *)calloc(n + 1, sizeof(float));

      // indexing p at 1 because of NR legacy
      for (int i = 0; i < n; i++) {
        float_parms[i + 1] = static_cast< float >(p(i));
      }
      (*powell_iteration_func)(float_parms, n);
      free(float_parms);
    }
    for (int i = 0; i < n; i++) {
      if (verbose_) printf("    param %3d  nthiter %2d\n", (int)i, (int)num_iterations_);

      // xit = ith column of xi
      for (int j = 0; j < n; ++j) xit[j] = (*xi)[j][i];
      double fptt = fret;

      // 1D minimization along xi
      f1d.init(p, xit);
      vnl_brent brent(&f1d);
      double ax = 0.0;
      double xx = initial_step_;
      double bx = 0;
      if (verbose_) printf("  bracketing\n");
      brent.bracket_minimum(&ax, &xx, &bx);
      if (verbose_) printf("  bracket ax=%g xx=%g bx=%g\n", ax, xx, bx);
      fret = brent.minimize_given_bounds(bx, xx, ax, linmin_xtol_, &xx);
      if (verbose_) printf("  min xx=%g, fret=%g\n", xx, fret);
      f1d.uninit(xx, p);
      // Now p is minimizer along xi

      if (vcl_fabs(fptt - fret) > del) {
        del = vcl_fabs(fptt - fret);
        ibig = i;
      }
      if (verbose_) printf("  i = %d  fp=%f fret=%f\n", (int)i, fp, fret);
    }
    if (verbose_) {
      printf("  niters %d  fp=%f fret=%f\n", (int)num_iterations_, fp, fret);
      printf("  check: %lf <= %lf\n", 2.0 * vcl_fabs(fp - fret), ftol * (vcl_fabs(fp) + vcl_fabs(fret)));
    }

    if (2.0 * vcl_fabs(fp - fret) <= ftol * (vcl_fabs(fp) + vcl_fabs(fret))) {
#ifdef DEBUG
      vnl_matlab_print(vcl_cerr, xi, "xi");
#endif
      if (verbose_) printf("  converged\n");
      return CONVERGED_FTOL;
    }

    if (num_iterations_ == unsigned(maxfev)) {
      return FAILED_USER_REQUEST;
    }

    for (int j = 0; j < n; ++j) {
      ptt[j] = 2.0 * p[j] - pt[j];
      xit[j] = p[j] - pt[j];
      pt[j] = p[j];
    }

    double fptt = functor_->f(ptt);
    report_eval(fret);
    if (fptt < fp) {
      double t = 2.0 * (fp - 2.0 * fret + fptt) * vnl_math_sqr(fp - fret - del) - del * vnl_math_sqr(fp - fptt);
      if (t < 0.0) {
        f1d.init(p, xit);
        vnl_brent brent(&f1d);
        double ax = 0.0;
        double xx = 1.0;
        double bx = 0.0;
        brent.bracket_minimum(&ax, &xx, &bx);
        fret = brent.minimize_given_bounds(bx, xx, ax, linmin_xtol_, &xx);
        f1d.uninit(xx, p);

        for (int j = 0; j < n; j++) {
          (*xi)[j][ibig] = (*xi)[j][n - 1];
          (*xi)[j][n - 1] = xit[j];
        }
      }
    }
    report_iter();
  }
  if (verbose_) printf("  niters %d\n", (int)num_iterations_);
  return FAILED_USER_REQUEST;
}
