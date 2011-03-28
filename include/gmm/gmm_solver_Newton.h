// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/**@file gmm_solver_Newton.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @author  Michel Fournie <fournie@mip.ups-tlse.fr>
   @date January 24, 2006.
*/
#ifndef GMM_SOLVERS_NEWTON_H__
#define GMM_SOLVERS_NEWTON_H__ 

namespace gmm {

#include "gmm_kernel.h"

  /* ***************************************************************** */
  /*     Line search definition                                        */
  /* ***************************************************************** */

  struct abstract_newton_line_search {
    double conv_alpha, conv_r;
    size_t it, itmax, glob_it;
    virtual void init_search(double r, size_t git) = 0;
    virtual double next_try(void) = 0;
    virtual bool is_converged(double) = 0;
    virtual double converged_value(void) { return conv_alpha; };
    virtual double converged_residual(void) { return conv_r; };
    virtual ~abstract_newton_line_search() {}
  };


  struct simplest_newton_line_search : public abstract_newton_line_search {
    double alpha, alpha_mult, first_res, alpha_max_ratio, alpha_min;
    virtual void init_search(double r, size_t git) {
      glob_it = git;
      conv_alpha = alpha = double(1); conv_r = first_res = r; it = 0;
    }
    virtual double next_try(void)
    { conv_alpha = alpha; alpha *= alpha_mult; ++it; return conv_alpha; }
    virtual bool is_converged(double r) {
      conv_r = r;
      return ((it <= 1 && r < first_res)
	      || (r <= first_res * alpha_max_ratio)
	      || (conv_alpha <= alpha_min)
	      || it >= itmax);
    }
    simplest_newton_line_search
    (size_t imax = size_t(-1), double a_max_ratio = 6.0/5.0,
     double a_min = 1.0/1000.0, double a_mult = 3.0/5.0)
      : alpha_mult(a_mult), alpha_max_ratio(a_max_ratio), alpha_min(a_min)
      { itmax = imax; }
  };

  struct default_newton_line_search : public abstract_newton_line_search {
    double alpha, alpha_mult, first_res, alpha_max_ratio;
    double alpha_min, prev_res, alpha_max_augment;
    virtual void init_search(double r, size_t git) {
      glob_it = git;
      conv_alpha = alpha = double(1);
      prev_res = conv_r = first_res = r; it = 0;
    }
    virtual double next_try(void)
    { conv_alpha = alpha; alpha *= alpha_mult; ++it; return conv_alpha; }
    virtual bool is_converged(double r) {
      if (glob_it == 0 || (r < first_res / double(2))
	  || (conv_alpha <= alpha_min && r < first_res * alpha_max_augment)
	  || it >= itmax)
	{ conv_r = r; return true; }
      if (it > 1 && r > prev_res && prev_res < alpha_max_ratio * first_res)
	return true;
      prev_res = conv_r = r;
      return false;
    }
    default_newton_line_search
    (size_t imax = size_t(-1),
     double a_max_ratio = 5.0/3.0,
     double a_min = 1.0/1000.0, double a_mult = 3.0/5.0, double a_augm = 2.0)
      : alpha_mult(a_mult), alpha_max_ratio(a_max_ratio),
	alpha_min(a_min), alpha_max_augment(a_augm) { itmax = imax; }
  };


  struct systematic_newton_line_search : public abstract_newton_line_search {
    double alpha, alpha_mult, first_res;
    double alpha_min, prev_res;
    bool first;
    virtual void init_search(double r, size_t git) {
      glob_it = git;
      conv_alpha = alpha = double(1);
      prev_res = conv_r = first_res = r; it = 0; first = true;
    }
    virtual double next_try(void)
    { double a = alpha; alpha *= alpha_mult; ++it; return a; }
    virtual bool is_converged(double r) {
      // cout << "a = " << alpha / alpha_mult << " r = " << r << endl;
      if (r < conv_r || first)
	{ conv_r = r; conv_alpha = alpha / alpha_mult; first = false; }
      if ((alpha <= alpha_min*alpha_mult) || it >= itmax) return true;
      return false;
    }
    systematic_newton_line_search
    (size_t imax = size_t(-1),
     double a_min = 1.0/10000.0, double a_mult = 3.0/5.0)
      : alpha_mult(a_mult), alpha_min(a_min)  { itmax = imax; }
  };

}


#endif
