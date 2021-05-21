/**
 * @brief Wrappers for core math routines from open-sources: VNL and CEPHES.
 */
/*
 * Original Author:  Dennis Jen and Silvester Czanner
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


#ifndef NUMERICS_H_
#define NUMERICS_H_

#include "matrix.h"


  MATRIX* MatrixSqrt(MATRIX * m, MATRIX * sqrtm);

#define SPLINE_USE_QUADRATIC         0
#define SPLINE_USE_FIRST_DERIVATIVE  1
#define SPLINE_USE_SECOND_DERIVATIVE 2

  float OpenBetaIncomplete(float a, float b, float x);

  float OpenGammaIncomplete(float a, float x);

  float OpenRan1( long *seed );

  // this is only called from the matrix class
  int OpenSvdcmp( MATRIX *a, VECTOR *w, MATRIX *v );

  int OpenQRdecomposition(const MATRIX *iMatrix, MATRIX *oQ, MATRIX *oR);

  float OpenMatrixDeterminant( MATRIX *matrix );

  int OpenLUMatrixInverse( MATRIX *matrix, MATRIX *inverse );

  void OpenPowell( float p[], float **ioInitialDirection, int n, float ftol,
                   int *iter, float *fret,
                   float (*func)(float []) );
  int OpenPowell2( float p[], float **ioInitialDirection, int n, float ftol,
		   float linmintol, int maxiters, int *iter, float *fret,
		   float (*func)(float []) );

  int OpenEigenSystem( float *iaData, int icData, float *oEigenValues,
                       float *oEigenVectors );

  int OpenNonSymmetricEigenSystem( float *iaData,
                                   int icData,
                                   float *oEigenValues,
                                   float *oEigenVectors );

  int OpenDFPMin( float p[], int n, float ftol, int *iter, float *fret,
                   float(*func)(float []), void (*dfunc)(float [], float []),
                   // addition
                   void (*step_func)(int itno,
                                     float sse,
                                     void *parms,
                                     float *p),
                   void *parms,
                   void (*user_callback_function)(float[]) );

  void OpenSpline( float x[], float y[], int n, float yp1, float ypn,
                   float y2[] );

  void OpenSplint( float xa[], float ya[], float y2a[], int n, float x,
                   float *y );

  float *SplineCubicSet ( int n,
                          float t[],
                          float y[],
                          int ibcbeg,
                          float ybcbeg,
                          int ibcend,
                          float ybcend );

  float SplineCubicValue ( int n, float t[], float tval, float y[],
                           float ypp[], float *ypval, float *yppval );

  float *d3_np_fs ( int n, float a[], float b[] );

  float *vector(long nl, long nh);
  float **matrix(long nrl, long nrh, long ncl, long nch);
  void free_vector(float *v, long nl, long nh);
  void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);



#ifdef INFINITY
# define SC_POSINF INFINITY
# define SC_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define SC_POSINF HUGE_VAL
# define SC_NEGINF (-HUGE_VAL)
#endif


  //########################### RNG #####################

  typedef struct
  {
    unsigned int i;
    unsigned int j;
    unsigned int n;
    unsigned int skip;
    unsigned int carry;
    unsigned long int u[24];
  }
  sc_status_t;

  typedef struct
  {
    const char *name;
    unsigned long int min;
    unsigned long int max;
    int size;
  }
  sc_rng_type;

  typedef struct
  {
    const sc_rng_type *type;
    void  *status;
  }
  sc_rng;

  static const sc_rng_type intern_rng_type =
    {"rng_generator",
     0x00ffffffUL,
     0,
     sizeof(sc_status_t)
    };


  //##################### FUNCTION PROTOTYPES ########################

  int sc_linalg_cholesky_decomp(MATRIX *U);


  //##########
  sc_rng           *sc_rng_alloc(const sc_rng_type *T);
  void              sc_rng_free(sc_rng *r);

  void              sc_rng_set(sc_rng *r, unsigned long int seed);
  unsigned long int sc_rng_get (const sc_rng * R);

  //##########

  double sc_ran_flat(const sc_rng *r, const double a, const double b);
  double sc_ran_gaussian (const sc_rng *r, const double sigma);
  double sc_ran_fdist (const sc_rng * r, const double nu1, const double nu2);

  double sc_ran_gamma(const sc_rng * r, const double a, const double b);
  double sc_ran_chisq(const sc_rng * r, const double nu);
  double sc_ran_tdist(const sc_rng * r, const double nu);

  double sc_ran_exponential(const sc_rng * r, const double mu);
  double sc_ran_binomial_pdf( unsigned int k, double p, unsigned int n);

  double sc_cdf_flat_Q(double x, double a, double b);
  double sc_cdf_flat_Qinv(double Q, double a, double b);

  double sc_cdf_fdist_Q(double x, double nu1, double nu2);
  double sc_cdf_fdist_Qinv(double Q, double nu1, double nu2);

  double sc_cdf_tdist_Q(double x, double nu);
  double sc_cdf_tdist_Qinv(double Q, double nu);

  double sc_cdf_gaussian_Q(double x, double nu);
  double sc_cdf_gaussian_Qinv(double Q, double nu);

  double sc_cdf_chisq_Q(double x, double nu);
  double sc_cdf_chisq_Qinv(double Q, double nu);

#endif /*NUMERICS_H_*/
