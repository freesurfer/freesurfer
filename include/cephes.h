/**
 * @file  cephes.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.5 $
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


#if !defined(__CEPHES_H)
#define __CEPHES_H

#ifdef __cplusplus
extern "C"
{
#endif

  extern double bdtrc ( int k, int n, double p );
  extern double bdtr ( int k, int n, double p );
  extern double bdtri ( int k, int n, double y );
  extern double btdtr ( double a, double b, double x );
  extern double chdtrc ( double df, double x );
  extern double chdtr ( double df, double x );
  extern double chdtri ( double df, double y );
  extern double fdtrc ( double ia, double ib, double x );
  extern double fdtr ( double ia, double ib, double x );
  extern double fdtri ( double ia, double ib, double y );
  extern double gamma ( double x );
  extern double lgam ( double x );
  extern double gdtr ( double a, double b, double x );
  extern double gdtrc ( double a, double b, double x );
  extern double igamc ( double a, double x );
  extern double igam ( double a, double x );
  extern double igami ( double a, double y0 );
  extern double incbet ( double aa, double bb, double xx );
  extern double incbi ( double aa, double bb, double yy0 );
  extern int mtherr ( char *name, int code );
  extern double nbdtrc ( int k, int n, double p );
  extern double nbdtr ( int k, int n, double p );
  extern double nbdtri ( int k, int n, double p );
  extern double ndtr ( double a );
  extern double erfc ( double a );
  extern double erf ( double x );
  extern double ndtri ( double y0 );
  extern double pdtrc ( int k, double m );
  extern double pdtr ( int k, double m );
  extern double pdtri ( int k, double y );
  extern double stdtr ( double k, double t );
  extern double stdtri ( double k, double p );
  extern double log1p ( double x );
  extern double expm1 ( double x );
  extern double cos1m ( double x );
  extern double polevl ( double x, void *P, int n );
  extern double p1evl ( double x, void *P, int n );

#ifdef __cplusplus
}
#endif

#endif
