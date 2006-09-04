#ifdef VCL_NEEDS_PRAGMA_INTERFACE
#pragma implementation
#endif
//:
// \file
//
// \author Andrew W. Fitzgibbon, Oxford RRG
// \date   22 Aug 99
//
//-----------------------------------------------------------------------------

#include "fs_vnl/fs_lbfgs.h"
#include <vcl_cmath.h>
#include <vcl_iostream.h>
#include <vcl_iomanip.h> // for setw (replaces cout.form())
#include <vnl/algo/vnl_netlib.h> // lbfgs_()

extern "C" {
#include <lbfgs.h> // from netlib, for lb3_ data struct
}


//: Default constructor.
// memory is set to 5, line_search_accuracy to 0.9.
// Calls init_parameters
fs_lbfgs::fs_lbfgs():
  f_(0)
{
  init_parameters();
}


//: Constructor. f is the cost function to be minimized.
// Calls init_parameters
fs_lbfgs::fs_lbfgs(vnl_cost_function& f):
  f_(&f)
{
  init_parameters();
}


//: Called by constructors.
// Memory is set to 5, line_search_accuracy to 0.9, default_step_length to 1.
void fs_lbfgs::init_parameters()
{
  memory = 5;
  line_search_accuracy = 0.9;
  default_step_length = 1.0;
  step_function_ = NULL;
  step_function_parms_ = NULL;
  mUserCallbackFunction = NULL;
}


void fs_lbfgs::set_step_function
( void ( *step_function )
  ( int itno, float sse, void *parms, float *p ), void *parms)
{
  this->step_function_ = step_function;
  this->step_function_parms_ = parms;
}


void fs_lbfgs::set_user_callback_function
( void (*userCallbackFunction)(float []) )
{
  mUserCallbackFunction = userCallbackFunction;
}


bool fs_lbfgs::minimize(vnl_vector<double>& x)
{
  /* Local variables */
  /*     The driver for fs_lbfgs must always declare LB2 as EXTERNAL */

  int n = f_->get_number_of_unknowns();
  int m = memory; // The number of basis vectors to remember.

  int iprint[2] = {1, 0};
  vnl_vector<double> g(n);

  // Workspace
  vnl_vector<double> diag(n);

  vnl_vector<double> w(n * (2*m+1)+2*m);

  if (verbose_)
    vcl_cerr << "ERROR fs_lbfgs: n = "
             << n
             <<", memory = "
             << m
             <<", Workspace = "
             << w.size()
             << "[ "
             << ( w.size() / 128.0 / 1024.0)
             <<" MB], ErrorScale = "
             << f_->reported_error(1)
             <<", xnorm = "
             << x.magnitude()
             << vcl_endl;

  bool we_trace = (verbose_ && !trace);

  if (we_trace)
    vcl_cerr << "ERROR fs_lbfgs: ";

  double best_f = 0;
  vnl_vector<double> best_x;

  bool ok;
  this->num_evaluations_ = 0;
  this->num_iterations_ = 0;
  int iflag = 0;
  while (true) {

    // TODO: addition
    if( mUserCallbackFunction != NULL ) {

      float *currentX = new float[n+1];
      for(int i=0; i<n; i++) {
        // TODO: legacy one indexing
        currentX[ i+1 ] = static_cast< float >( x(i) );
      }
      ( *mUserCallbackFunction )( currentX );
      delete[] currentX;
    }

    // We do not wish to provide the diagonal
    // matrices Hk0, and therefore set DIAGCO to FALSE.
    logical diagco = false;

    // Set these every iter in case user changes them to bail out
    double eps = gtol; // Gradient tolerance
    double local_xtol = 1e-16;
    lb3_.gtol = line_search_accuracy; // set to 0.1 for huge
                                      //problems or cheap functions
    lb3_.stpawf = default_step_length;

    // Call function
    double f;
    f_->compute(x, &f, &g);
    if (this->num_evaluations_ == 0) {
      this->start_error_ = f;
      best_f = f;
    } else if (f < best_f) {
      best_x = x;
      best_f = f;
    }

#define print_(i,a,b,c,d) vcl_cerr<<vcl_setw(6)<<i<<' '<<vcl_setw(20)<<a<<' '\
           <<vcl_setw(20)<<b<<' '<<vcl_setw(20)<<c<<' '<<vcl_setw(20)<<d<<'\n'

    if (check_derivatives_)
      {
        vcl_cerr << "ERROR fs_lbfgs: f = "
                 << f_->reported_error(f)
                 << ", computing FD gradient\n";
        vnl_vector<double> fdg = f_->fdgradf(x);
        if (verbose_)
          {
            int l = n;
            int limit = 100;
            int limit_tail = 10;
            if (l > limit + limit_tail) {
              vcl_cerr << " [ Showing only first "
                       <<limit
                       << " components ]\n";
              l = limit;
            }
            print_("i","x","g","fdg","dg");
            print_("-","-","-","---","--");
            for (int i = 0; i < l; ++i)
              print_(i, x[i], g[i], fdg[i], g[i]-fdg[i]);
            if (n > limit) {
              vcl_cerr << "   ...\n";
              for (int i = n - limit_tail; i < n; ++i)
                print_(i, x[i], g[i], fdg[i], g[i]-fdg[i]);
            }
          }
        vcl_cerr << "   ERROR = "
                 << (fdg - g).squared_magnitude() / vcl_sqrt(double(n))
                 << "\n";
      }

    iprint[0] = trace ? 1 : -1; // -1 no o/p, 0 start and end, 1 every iter.
    iprint[1] = 0; // 1 prints X and G

    lbfgs_(&n, &m, x.data_block(),
           &f, g.data_block(),
           &diagco, diag.data_block(),
           iprint, &eps, &local_xtol, w.data_block(), &iflag);

    ++this->num_iterations_;

    if (we_trace)
      vcl_cerr << iflag << ":" << f_->reported_error(f) << " ";

    if (iflag == 0) {
      // Successful return
      this->end_error_ = f;
      ok = true;

      // TODO: this is an DJ addition for when the
      // best_x is not initialized yet
      if( best_x.size() != 0 ) {
        x = best_x;
      }
      break;
    }

    if (iflag < 0) {
      // Eeek.
      vcl_cerr << "fs_lbfgs: ** ERROR ** (iflag < 0)\n";
      ok = false;

      // TODO: this is an DJ addition for when
      // the best_x is not initialized yet
      if( best_x.size() != 0 ) {
        this->end_error_ = f;
        x = best_x;
      }

      break;
    }

    if (++this->num_evaluations_ > get_max_function_evals()) {
      failure_code_ = FAILED_TOO_MANY_ITERATIONS;
      ok = false;

      // TODO: this is an DJ addition for when the
      // best_x is not initialized yet
      if( best_x.size() != 0 ) {
        this->end_error_ = f;
        x = best_x;
      }

      break;
    }

    // step function addition
    if (step_function_ != NULL) {

      // TODO: this could potentially be very expensive
      // TODO: legacy one indexing
      float *currentX = new float[n+1];
      for(int i=0; i<n; i++) {
        // TODO: legacy one indexing
        currentX[ i+1 ] = static_cast<float>(x(i));
      }
      (*step_function_)(this->num_iterations_, static_cast<float>(f),
                        step_function_parms_, currentX);
      for(int i=0; i<n; i++) {
        // TODO: legacy one indexing
        x(i) = static_cast<double>(currentX[ i+1 ]);
      }
      delete []currentX;
    }
  }
  if (we_trace) vcl_cerr << "done\n";

  return ok;
}

