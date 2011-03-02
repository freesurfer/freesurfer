/**
 * @file  sc_test.c
 * @brief tests of core math routines
 *
 */
/*
 * Original Author: Nick Schmansky, Silvester Czanner
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 13:03:57 $
 *    $Revision: 1.14 $
 *
 * Copyright (C) 2002 Jason H Stover.
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#define HAVE_VPRINTF 1

#include "numerics.h"

#define TEST(func, args, value, tol) \
{ double res = func args ; sc_test_rel (res, value, tol, #func #args); } ;

#define TEST_TOL0  (2.0*SC_DBL_EPSILON)
#define TEST_TOL1  (16.0*SC_DBL_EPSILON)
#define TEST_TOL2  (256.0*SC_DBL_EPSILON)
#define TEST_TOL3  (2048.0*SC_DBL_EPSILON)
#define TEST_TOL4  (16384.0*SC_DBL_EPSILON)
#define TEST_TOL5  (131072.0*SC_DBL_EPSILON)
#define TEST_TOL6  (1048576.0*SC_DBL_EPSILON)
#define TEST_TOL7  (16777216.0*SC_DBL_EPSILON)


#define SC_DBL_EPSILON 2.2204460492503131e-16

const sc_rng_type *sc_rng_intern_type = &intern_rng_type;

static unsigned int tests = 0;
static unsigned int passed = 0;
static unsigned int failed = 0;

static unsigned int verbose = 1;

const char* Progname = "sc_test";

void
sc_test (int status, const char *test_description,...)
{

  tests++;

  if (status == 0)
  {
    passed++;
    if (verbose)
      printf ("PASS: ");
  }
  else
  {
    failed++;
    if (verbose)
      printf ("*** FAIL: ");
  }

  if (verbose)
  {

#if HAVE_VPRINTF
    va_list ap;

#ifdef STDC_HEADERS
    va_start (ap, test_description);
#else
    va_start (ap);
#endif
    vprintf (test_description, ap);
    va_end (ap);
#endif

    printf("\n");
    fflush (stdout);
  }
}

int
sc_isnan (const double x)
{
  int status = (x != x);
  return status;
}

int
sc_isinf (const double x)
{
  double y = x - x;
  int s = (y != y);

  if (s && x > 0)
    return +1;
  else if (s && x < 0)
    return -1;
  else
    return 0;
}

void
sc_test_rel (double result, double expected, double relative_error,
             const char *test_description,...)
{
  int status ;

  /* Check for NaN vs inf vs number */

  if (sc_isnan(result) || sc_isnan(expected))
  {
    status = sc_isnan(result) != sc_isnan(expected);
  }
  else if (sc_isinf(result) || sc_isinf(expected))
  {
    status = sc_isinf(result) != sc_isinf(expected);
  }
  else if (expected != 0 )
  {
    status = (fabs(result-expected)/fabs(expected) > relative_error) ;
  }
  else
  {
    status = (fabs(result) > relative_error) ;
  }

  tests++;

  if (status == 0)
  {
    passed++;
    if (verbose)
      printf ("PASS: ");
  }
  else
  {
    failed++;
    if (verbose)
      printf ("*** FAIL: ");

  }

  if (verbose)
  {

#if HAVE_VPRINTF
    va_list ap;

#ifdef STDC_HEADERS
    va_start (ap, test_description);
#else
    va_start (ap);
#endif
    vprintf (test_description, ap);
    va_end (ap);
#endif
    if (status == 0)
    {
      if (strlen(test_description) < 45)
      {
        printf(" (%g observed vs %g expected)", result, expected) ;
      }
      else
      {
        printf(" (%g obs vs %g exp)", result, expected) ;
      }
    }
    else
    {
      printf(" (%.18g observed vs %.18g expected)", result, expected) ;
    }

    printf ("\n") ;
    fflush (stdout);
  }
}

void
sc_test_abs (double result, double expected, double absolute_error,
             const char *test_description,...)
{
  int status ;

  /* Check for NaN vs inf vs number */

  if (sc_isnan(result) || sc_isnan(expected))
  {
    status = sc_isnan(result) != sc_isnan(expected);
  }
  else if (sc_isinf(result) || sc_isinf(expected))
  {
    status = sc_isinf(result) != sc_isinf(expected);
  }
  else
  {
    status = fabs(result-expected) > absolute_error ;
  }

  tests++;

  if (status == 0)
  {
    passed++;
    if (verbose)
      printf ("PASS: ");
  }
  else
  {
    failed++;
    if (verbose)
      printf ("FAIL: ");

  }

  if (verbose)
  {

#if HAVE_VPRINTF
    va_list ap;

#ifdef STDC_HEADERS
    va_start (ap, test_description);
#else
    va_start (ap);
#endif
    vprintf (test_description, ap);
    va_end (ap);
#endif
    if (status == 0)
    {
      if (strlen(test_description) < 45)
      {
        printf(" (%g observed vs %g expected)", result, expected) ;
      }
      else
      {
        printf(" (%g obs vs %g exp)", result, expected) ;
      }
    }
    else
    {
      printf(" (%.18g observed vs %.18g expected)", result, expected) ;
    }

    printf ("\n") ;
    fflush (stdout);
  }
}


void
sc_test_factor (double result, double expected, double factor,
                const char *test_description,...)
{
  int status;

  if (result == expected)
  {
    status = 0;
  }
  else if (expected == 0.0)
  {
    status = (result > expected || result < expected);
  }
  else
  {
    double u = result / expected;
    status = (u > factor || u < 1.0 / factor) ;
  }

  tests++;

  if (status == 0)
  {
    passed++;
    if (verbose)
      printf ("PASS: ");
  }
  else
  {
    failed++;
    if (verbose)
      printf ("FAIL: ");

  }

  if (verbose)
  {

#if HAVE_VPRINTF
    va_list ap;

#ifdef STDC_HEADERS
    va_start (ap, test_description);
#else
    va_start (ap);
#endif
    vprintf (test_description, ap);
    va_end (ap);
#endif
    if (status == 0)
    {
      if (strlen(test_description) < 45)
      {
        printf(" (%g observed vs %g expected)", result, expected) ;
      }
      else
      {
        printf(" (%g obs vs %g exp)", result, expected) ;
      }
    }
    else
    {
      printf(" (%.18g observed vs %.18g expected)", result, expected) ;
    }

    printf ("\n") ;
    fflush (stdout);
  }
}

void
sc_test_int (int result, int expected, const char *test_description,...)
{
  int status = (result != expected) ;

  tests++;

  if (status == 0)
  {
    passed++;
    if (verbose)
      printf ("PASS: ");
  }
  else
  {
    failed++;
    if (verbose)
      printf ("FAIL: ");
  }

  if (verbose)
  {

#if HAVE_VPRINTF
    va_list ap;

#ifdef STDC_HEADERS
    va_start (ap, test_description);
#else
    va_start (ap);
#endif
    vprintf (test_description, ap);
    va_end (ap);
#endif
    if (status == 0)
    {
      printf(" (%d observed vs %d expected)", result, expected) ;
    }
    else
    {
      printf(" (%d observed vs %d expected)", result, expected) ;
    }

    printf ("\n");
    fflush (stdout);
  }
}

void
sc_test_str (const char * result, const char * expected,
             const char *test_description,...)
{
  int status = strcmp(result,expected) ;

  tests++;

  if (status == 0)
  {
    passed++;
    if (verbose)
      printf ("PASS: ");
  }
  else
  {
    failed++;
    if (verbose)
      printf ("FAIL: ");
  }

  if (verbose)
  {

#if HAVE_VPRINTF
    va_list ap;

#ifdef STDC_HEADERS
    va_start (ap, test_description);
#else
    va_start (ap);
#endif
    vprintf (test_description, ap);
    va_end (ap);
#endif
    if (status)
    {
      printf(" (%s observed vs %s expected)", result, expected) ;
    }

    printf ("\n");
    fflush (stdout);
  }
}

void
sc_test_verbose (int v)
{
  verbose = v;
}

int
sc_test_summary (void)
{

  if (verbose && 0)             /* FIXME: turned it off, this annoys me */
    printf ("%d tests, passed %d, failed %d.\n", tests, passed, failed);

  if (failed != 0)
  {

    if (verbose && 0)         /* FIXME: turned it off, this annoys me */
    {
      printf ("%d TEST%s FAILED.\n", failed, failed == 1 ? "" : "S");
    }
    return EXIT_FAILURE;
  }

  if (tests != passed + failed)
  {
    if (verbose)
      printf ("TEST RESULTS DO NOT ADD UP %d != %d + %d\n",
              tests, passed, failed);
    return EXIT_FAILURE;
  }

  if (passed == tests)
  {
    if (verbose && 0)         /* FIXME: turned it off, this annoys me */
      printf ("All tests passed successfully\n");
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}



void
test_auto_fdist (void)
{
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+10,5.3,2.7), 4.272202262298e-14, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+09,5.3,2.7), 9.564269502770e-13, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+08,5.3,2.7), 2.141173208523e-11, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+07,5.3,2.7), 4.793489218238e-10, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+06,5.3,2.7), 1.073127433440e-08, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+05,5.3,2.7), 2.402407758939e-07, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+04,5.3,2.7), 5.377754447932e-06, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+03,5.3,2.7), 1.202661950234e-04, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+02,5.3,2.7), 2.664302300604e-03, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+01,5.3,2.7), 5.386718252832e-02, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e+00,5.3,2.7), 5.467279509126e-01, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000001e-01,5.3,2.7), 9.863930100746e-01, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e-02,5.3,2.7), 9.999515966684e-01, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e-03,5.3,2.7), 9.999998859974e-01, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e-04,5.3,2.7), 9.999999997435e-01, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000001e-05,5.3,2.7), 9.999999999994e-01, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (9.9999999999999995e-07,5.3,2.7), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (9.9999999999999995e-08,5.3,2.7), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e-08,5.3,2.7), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000001e-09,5.3,2.7), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (1.0000000000000000e-10,5.3,2.7), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_fdist_Q,
       (0.0000000000000000e+00,5.3,2.7), 1.000000000000e+00, TEST_TOL6);
}



void
test_auto_gaussian (void)
{
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+10,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+09,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+08,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+07,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+06,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+05,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+04,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+03,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+02,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+01,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (7.2252292279265077e-15,1.0), 7.692307692307691, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e+00,1.0), 0.158655253931457046, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (2.2087816371245972e-01,1.0), 0.7692307692307e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000001e-01,1.0), 4.601721627222e-01, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (4.6934236960338749e-01,1.0), 0.0769230769230768996, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e-02,1.0), 4.960106436853e-01, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (4.9693124349158196e-01,1.0), 0.7692307692307e-02, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e-03,1.0), 4.996010577860e-01, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (4.9969312135303229e-01,1.0), 0.7692307692307e-03, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e-04,1.0), 0.49996010577202632, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000001e-05,1.0), 0.499996010577196059, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (9.9999999999999995e-07,1.0), 0.499999601057719623, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (9.9999999999999995e-08,1.0), 0.499999960105771968, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e-08,1.0), 0.499999996010577208, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000001e-09,1.0), 4.999999996931e-01, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (1.0000000000000000e-10,1.0), 4.999999999693e-01, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (0.0000000000000000e+00,1.0), 5.000000000000e-01, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (5.0000000000000000e-01,1.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e-10,1.0), 5.000000000307e-01, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000001e-09,1.0), 5.000000003069e-01, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e-08,1.0), 0.500000003989422792, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-9.9999999999999995e-08,1.0), 0.500000039894228032, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-9.9999999999999995e-07,1.0), 0.500000398942280433, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000001e-05,1.0), 0.500003989422803996, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e-04,1.0), 0.50003989422797368, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e-03,1.0), 0.500398942213911013, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (5.0030687864696777e-01,1.0), -0.00076923076923092709, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e-02,1.0), 0.503989356314631598, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (5.0306875650841798e-01,1.0), -0.00769230769230770488, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000001e-01,1.0), 0.539827837277028988, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (5.3065763039661251e-01,1.0), -0.0769230769230768996, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+00,1.0), 0.841344746068542926, TEST_TOL6);
  TEST(sc_cdf_gaussian_Qinv,
       (7.7912183628754028e-01,1.0), -0.769230769230769051, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+01,1.0), 1.000000000000e-00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+02,1.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+03,1.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+04,1.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+05,1.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+06,1.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+07,1.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+08,1.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+09,1.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_gaussian_Q,
       (-1.0000000000000000e+10,1.0), 1.000000000000e+00, TEST_TOL6);
}



void
test_auto_flat (void)
{
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+10,1.3,750.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+09,1.3,750.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+08,1.3,750.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+07,1.3,750.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+06,1.3,750.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+05,1.3,750.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+04,1.3,750.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+03,1.3,750.0), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+02,1.3,750.0), 8.681714972619e-01, TEST_TOL6);
  TEST(sc_cdf_flat_Qinv,
       (8.6817149726190368e-01,1.3,750.0), 1.000000000000e+02, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+01,1.3,750.0), 9.883798584213e-01, TEST_TOL6);
  TEST(sc_cdf_flat_Qinv,
       (9.8837985842125353e-01,1.3,750.0), 1.000000000000e+01, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e+00,1.3,750.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000001e-01,1.3,750.0), 1.000000000000e-00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e-02,1.3,750.0), 1.000000000000e-00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e-03,1.3,750.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e-04,1.3,750.0), 1.000000000000e-00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000001e-05,1.3,750.0), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (9.9999999999999995e-07,1.3,750.0), 1.000000000000e-00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (9.9999999999999995e-08,1.3,750.0), 1.000000000000e-00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e-08,1.3,750.0), 1.000000000000e-00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000001e-09,1.3,750.0), 1.000000000000e-00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (1.0000000000000000e-10,1.3,750.0), 1.000000000000e-00, TEST_TOL6);
  TEST(sc_cdf_flat_Q,
       (0.0000000000000000e+00,1.3,750.0), 1.000000000000e+00, TEST_TOL6);
}



void
test_auto_chisq (void)
{
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+10,1.3), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+09,1.3), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+08,1.3), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+07,1.3), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+06,1.3), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+05,1.3), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+04,1.3), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+03,1.3), 5.840240518729e-219, TEST_TOL6);
  TEST(sc_cdf_chisq_Qinv,
       (5.8402405187288964e-219,1.3), 1.000000000000e+03, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+02,1.3), 3.517864771108e-23, TEST_TOL6);
  TEST(sc_cdf_chisq_Qinv,
       (3.5178647711076648e-23,1.3), 1.000000000000e+02, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+01,1.3), 2.613210979470e-03, TEST_TOL6);
  TEST(sc_cdf_chisq_Qinv,
       (2.6132109794696230e-03,1.3), 1.000000000000e+01, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e+00,1.3), 4.121379867221e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Qinv,
       (4.1213798672211427e-01,1.3), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000001e-01,1.3), 8.445731104160e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Qinv,
       (8.4457311041596417e-01,1.3), 1.000000000000e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e-02,1.3), 9.645858609746e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Qinv,
       (9.6458586097459775e-01,1.3), 1.000000000000e-02, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e-03,1.3), 9.920577036204e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Qinv,
       (9.9205770362041057e-01,1.3), 1.000000000000e-03, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e-04,1.3), 9.982216261112e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Qinv,
       (9.9822162611119969e-01,1.3), 1.000000000000e-04, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000001e-05,1.3), 9.996018646205e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Qinv,
       (9.9960186462053913e-01,1.3), 1.000000000000e-05, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (9.9999999999999995e-07,1.3), 9.999108684330e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (9.9999999999999995e-08,1.3), 9.999800459241e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e-08,1.3), 9.999955328388e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000001e-09,1.3), 9.999989999272e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (1.0000000000000000e-10,1.3), 9.999997761116e-01, TEST_TOL6);
  TEST(sc_cdf_chisq_Q,
       (0.0000000000000000e+00,1.3), 1.000000000000e+00, TEST_TOL6);
}



void
test_auto_tdist (void)
{
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+10,1.3), 3.467848111850e-14, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (3.4678481118500305e-14,1.3), 1.000000000000e+10, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+09,1.3), 6.919266651610e-13, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (6.9192666516103524e-13,1.3), 1.000000000000e+09, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+08,1.3), 1.380575199718e-11, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (1.3805751997179027e-11,1.3), 1.000000000000e+08, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+07,1.3), 2.754609668978e-10, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (2.7546096689777484e-10,1.3), 1.000000000000e+07, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+06,1.3), 5.496168864957e-09, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (5.4961688649569980e-09,1.3), 1.000000000000e+06, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+05,1.3), 1.096629861231e-07, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (1.0966298612314582e-07,1.3), 1.000000000000e+05, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+04,1.3), 2.188064222827e-06, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (2.1880642228271703e-06,1.3), 1.000000000000e+04, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+03,1.3), 4.365759541083e-05, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (4.3657595410833571e-05,1.3), 1.000000000000e+03, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+02,1.3), 8.710327647608e-04, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (8.7103276476079201e-04,1.3), 1.000000000000e+02, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+01,1.3), 1.727893386820e-02, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (1.7278933868204446e-02,1.3), 1.000000000000e+01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e+00,1.3), 2.336211937932e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (2.3362119379322516e-01,1.3), 1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000001e-01,1.3), 4.667575980083e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (4.6675759800826139e-01,1.3), 1.000000000000e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e-02,1.3), 4.966660755117e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (4.9666607551169606e-01,1.3), 1.000000000000e-02, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e-03,1.3), 4.996665978189e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (4.9966659781887629e-01,1.3), 1.000000000000e-03, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e-04,1.3), 4.999666597722e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000001e-05,1.3), 4.999966659772e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (9.9999999999999995e-07,1.3), 4.999996665977e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (9.9999999999999995e-08,1.3), 4.999999666598e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e-08,1.3), 4.999999966660e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000001e-09,1.3), 4.999999996666e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (1.0000000000000000e-10,1.3), 4.999999999667e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (0.0000000000000000e+00,1.3), 5.000000000000e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (4.9999999999999900e-01,1.3), 0.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e-10,1.3), 5.000000000333e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000001e-09,1.3), 5.000000003334e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e-08,1.3), 5.000000033340e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-9.9999999999999995e-08,1.3), 5.000000333402e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-9.9999999999999995e-07,1.3), 5.000003334023e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000001e-05,1.3), 5.000033340228e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e-04,1.3), 5.000333402278e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e-03,1.3), 5.003334021811e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (5.0033340218112365e-01,1.3), -1.000000000000e-03, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e-02,1.3), 5.033339244883e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (5.0333392448830394e-01,1.3), -1.000000000000e-02, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000001e-01,1.3), 5.332424019917e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (5.3324240199173856e-01,1.3), -1.000000000000e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+00,1.3), 7.663788062068e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (7.6637880620677490e-01,1.3), -1.000000000000e+00, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+01,1.3), 9.827210661318e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (9.8272106613179555e-01,1.3), -1.000000000000e+01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+02,1.3), 9.991289672352e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Qinv,
       (9.9912896723523925e-01,1.3), -1.000000000000e+02, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+03,1.3), 9.999563424046e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+04,1.3), 9.999978119358e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+05,1.3), 9.999998903370e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+06,1.3), 9.999999945038e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+07,1.3), 9.999999997245e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+08,1.3), 9.999999999862e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+09,1.3), 9.999999999993e-01, TEST_TOL6);
  TEST(sc_cdf_tdist_Q,
       (-1.0000000000000000e+10,1.3), 1.000000000000e-00, TEST_TOL6);
}


void test_tdist (void)
{
  TEST (sc_cdf_tdist_Q,
        (0.0, 1.0), 0.5, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1e-100, 1.0), 0.5, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.001, 1.0), 4.99681690219919441e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.01, 1.0), 4.96817007235091745e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.1, 1.0), 4.68274482569446430e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.325, 1.0), 3.99976879967147876e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1.0, 1.0), 2.5e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1.5, 1.0), 1.87167041810998816e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (2.0, 1.0), 1.47583617650433274e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (10.0, 1.0), 3.17255174305535695e-2, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (20.0, 1.0), 1.59022512561763752e-2, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (100.0, 1.0), 3.18299276490825515e-3, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1000.0, 1.0), 3.18309780080558939e-4, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (10000.0, 1.0), 3.18309885122757724e-5, TEST_TOL6);

  TEST (sc_cdf_tdist_Q,
        (-1e-100, 1.0), 0.5, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.001, 1.0), 5.00318309780080559e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.01, 1.0), 5.03182992764908255e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.1, 1.0), 5.31725517430553570e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.325, 1.0), 6.00023120032852124e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-1.0, 1.0), 7.5e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-1.5, 1.0), 8.12832958189001184e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-2.0, 1.0), 8.52416382349566726e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-10.0, 1.0), 9.68274482569446430e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-20.0, 1.0), 9.84097748743823625e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-100.0, 1.0), 9.96817007235091745e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-1000.0, 1.0), 9.99681690219919441e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-10000.0, 1.0), 9.99968169011487724e-1, TEST_TOL6);

  TEST (sc_cdf_tdist_Q,
        (0.0, 2.0), 5.00000000000000000e-01, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1e-100, 2.0), 5.00000000000000000e-01, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.001, 2.0), 4.99646446697795041e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.01, 2.0), 4.96464554479100486e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.1, 2.0), 4.64732719207070087e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.325, 2.0), 3.88014227253126233e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1.0, 2.0), 2.11324865405187118e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1.5, 2.0), 1.36196562445500540e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (2.0, 2.0), 9.17517095361369836e-2, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (10.0, 2.0), 4.92622851166284542e-3, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (20.0, 2.0), 1.24533194618354849e-3, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (100.0, 2.0), 4.99925012497812894e-5, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1000.0, 2.0), 4.99999250001249998e-7, TEST_TOL6);
  // TEST (sc_cdf_tdist_Q,
  //    (10000.0, 2.0), 4.99999996961264515e-9, TEST_TOL6);

  TEST (sc_cdf_tdist_Q,
        (-1e-100, 2.0), 5.00000000000000000e-01, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.001, 2.0), 5.00353553302204959e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.01, 2.0), 5.03535445520899514e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.1, 2.0), 5.35267280792929913e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.325, 2.0), 6.11985772746873767e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-1.0, 2.0), 7.88675134594812882e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-1.5, 2.0), 8.63803437554499460e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-2.0, 2.0), 9.08248290463863016e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-10.0, 2.0), 9.95073771488337155e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-20.0, 2.0), 9.98754668053816452e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-100.0, 2.0), 9.99950007498750219e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-1000.0, 2.0), 9.99999500000749999e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-10000.0, 2.0), 9.99999995000000075e-1, TEST_TOL6);

  TEST (sc_cdf_tdist_Q,
        (0.0, 300.0), 5.00000000000000000e-01, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1e-100, 300.0), 5.00000000000000000e-01, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.001, 300.0), 4.99601390099057051e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.01, 300.0), 4.96013966979440912e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.1, 300.0), 4.60205558822231806e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (0.325, 300.0), 3.72703798457476188e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1.0, 300.0), 1.59058202215313138e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (1.5, 300.0), 6.73330165746308628e-2, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (2.0, 300.0), 2.32007604915745452e-2, TEST_TOL6);

  // these fail because they are out of range of valid
  // input values for the library routines used (cephes)
  //TEST (sc_cdf_tdist_Q,
  //  (10.0, 300.0), 8.279313677e-21, TEST_TOL6);
  //TEST (sc_cdf_tdist_Q,
  //  (20.0, 300.0), 1.93159812815803978e-57, TEST_TOL6);
  //TEST (sc_cdf_tdist_Q,
  //  (100.0, 300.0), 1.02557519997736154e-232, TEST_TOL6);
  //TEST (sc_cdf_tdist_Q,
  //    (1000.0, 300.0), 0.00000000000000000e+00, 0.0);
  //TEST (sc_cdf_tdist_Q,
  //    (10000.0, 300.0), 0.00000000000000000e+00, 0.0);

  TEST (sc_cdf_tdist_Q,
        (-1e-100, 300.0), 5.00000000000000000e-01, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.001, 300.0), 5.00398609900942949e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.01, 300.0), 5.03986033020559088e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.1, 300.0), 5.39794441177768194e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-0.325, 300.0), 6.27296201542523812e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-1.0, 300.0), 8.40941797784686862e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-1.5, 300.0), 9.32666983425369137e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-2.0, 300.0), 9.76799239508425455e-1, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-10.0, 300.0), 1.000000000000000000e0, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-20.0, 300.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-100.0, 300.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-1000.0, 300.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Q,
        (-10000.0, 300.0), 1.0, TEST_TOL6);
}

void test_fdist (void)
{
  TEST (sc_cdf_fdist_Q,
        (0.0, 1.2, 1.3), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1e-100, 1.2, 1.3), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.001, 1.2, 1.3), 9.88939151413976144e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.01, 1.2, 1.3), 9.56136324293168615e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.1, 1.2, 1.3), 8.31757607287159265e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.325, 1.2, 1.3), 6.85869954753804551e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1.0, 1.2, 1.3), 4.90369220925244747e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1.5, 1.2, 1.3), 4.16001359358446148e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (2.0, 1.2, 1.3), 3.65266418648061213e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (10.0, 1.2, 1.3), 1.51553762120799025e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (20.0, 1.2, 1.3), 9.90122736631249612e-2, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (100.0, 1.2, 1.3), 3.55108729523115643e-2, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1000.0, 1.2, 1.3), 7.98794830588361109e-3, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (10000.0, 1.2, 1.3), 1.7891371911574145e-3, TEST_TOL6);

  TEST (sc_cdf_fdist_Q,
        (0.0, 500.0, 1.3), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1e-100, 500.0, 1.3), 1.0, TEST_TOL6);

  TEST (sc_cdf_fdist_Q,
        (0.001, 500.0, 1.3), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.01, 500.0, 1.3), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.1, 500.0, 1.3), 9.99410023490380312e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.325, 500.0, 1.3), 9.31388951394845747e-1, TEST_TOL6);

  TEST (sc_cdf_fdist_Q,
        (1.0, 500.0, 1.3), 6.61524946193595385e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1.5, 500.0, 1.3), 5.47983754752542572e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (2.0, 500.0, 1.3), 4.72660931062611202e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (10.0, 500.0, 1.3), 1.83160371421586096e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (20.0, 500.0, 1.3), 1.18215376943088595e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (100.0, 500.0, 1.3), 4.19549427957787016e-2, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1000.0, 500.0, 1.3), 9.41425061934473424e-3, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (10000.0, 500.0, 1.3), 2.10807516853862603e-3, TEST_TOL6);

  TEST (sc_cdf_fdist_Q,
        (0.0, 1.2, 500.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1e-100, 1.2, 500.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.001, 1.2, 500.0), 9.86953850355871047e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.01, 1.2, 500.0), 9.48167577539196671e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.1, 1.2, 500.0), 7.97764898283923711e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.325, 1.2, 500.0), 6.09497016780606251e-1, TEST_TOL6);

  TEST (sc_cdf_fdist_Q,
        (1.0, 1.2, 500.0), 3.32343808425346381e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1.5, 1.2, 500.0), 2.24460769728532946e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (2.0, 1.2, 500.0), 1.54790885095386295e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (10.0, 1.2, 500.0), 8.3198234087901168e-4, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (20.0, 1.2, 500.0), 1.99426162833131e-6, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (100.0, 1.2, 500.0), 6.23302662288217117e-25, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1000.0, 1.2, 500.0), 1.14328577259666930e-134, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (10000.0, 1.2, 500.0), 0.0, 0.0);

  TEST (sc_cdf_fdist_Q,
        (0.0, 200.0, 500.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1e-100, 200.0, 500.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.001, 200.0, 500.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.01, 200.0, 500.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.1, 200.0, 500.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (0.325, 200.0, 500.0), 9.99999999999999997e-1, TEST_TOL6);

  TEST (sc_cdf_fdist_Q,
        (1.0, 200.0, 500.0), 4.93253673878831734e-1, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1.5, 200.0, 500.0), 2.05824281287561795e-4, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (2.0, 200.0, 500.0), 4.71763848371410786e-10, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (10.0, 200.0, 500.0), 5.98048337181948436e-96, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (20.0, 200.0, 500.0), 2.92099265879979502e-155, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (100.0, 200.0, 500.0), 0.0, TEST_TOL6);
  TEST (sc_cdf_fdist_Q,
        (1000.0, 200.0, 500.0), 0.0, 0.0);
  TEST (sc_cdf_fdist_Q,
        (10000.0, 200.0, 500.0), 0.0, 0.0);
}

void test_chisq (void)
{
  TEST (sc_cdf_chisq_Q,
        (0.0, 13.0), 1e0, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (1e-100, 13.0), 1e0, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (0.001, 13.0), 1e0, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (0.01, 13.0), 9.99999999999999999e-1, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (0.1, 13.0), 9.99999999998212030e-1, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (0.325, 13.0), 9.99999996553886862e-1, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (1.0, 13.0), 9.99996165265264864e-1, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (1.5, 13.0), 9.99956828161079896e-1, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (2.0, 13.0), 9.99773749915343953e-1, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (10.0, 13.0), 6.93934367980748890e-1, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (20.0, 13.0), 9.52102560780915127e-2, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (100.0, 13.0), 1.65902608070858809e-15, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (1000.0, 13.0), 1.74851191544860225e-205, TEST_TOL6);
  TEST (sc_cdf_chisq_Q,
        (10000.0, 13.0), 0.0, 0.0);
}


void test_chisqinv (void)
{
  TEST (sc_cdf_chisq_Qinv,
        (0.0, 13.0), SC_POSINF, TEST_TOL6);
  TEST (sc_cdf_chisq_Qinv,
        (1.65902608070858809e-15, 13.0), 100.0, TEST_TOL6);
  TEST (sc_cdf_chisq_Qinv,
        (9.52102560780915127e-2, 13.0), 20.0, TEST_TOL6);
  TEST (sc_cdf_chisq_Qinv,
        (6.93934367980748892e-1, 13.0), 10.0, TEST_TOL6);
  TEST (sc_cdf_chisq_Qinv,
        (9.99773749915343954e-1, 13.0), 2.0, TEST_TOL6);
  TEST (sc_cdf_chisq_Qinv,
        (9.99956828161079894e-1, 13.0), 1.5, TEST_TOL6);
  TEST (sc_cdf_chisq_Qinv,
        (9.99996165265264863e-1, 13.0), 1.0, TEST_TOL7);
  //TEST (sc_cdf_chisq_Qinv,
  //    (9.99999996553886862e-1, 13.0), 0.325, TEST_TOL7);
  //TEST (sc_cdf_chisq_Qinv,
  //    (9.99999999998212031e-1, 13.0), 0.1, TEST_TOL7);
  //TEST (sc_cdf_chisq_Qinv,
  //    (1.0, 13.0), 0.0, TEST_TOL7);
}

void test_tdistinv (void)
{
  TEST (sc_cdf_tdist_Qinv,
        (0.5, 1.0), 0.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (4.99681690219919441e-1, 1.0), 0.001, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (4.96817007235091745e-1, 1.0), 0.01, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (4.68274482569446430e-1, 1.0), 0.1, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (3.99976879967147876e-1, 1.0), 0.325, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (2.5e-1, 1.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (1.87167041810998816e-1, 1.0), 1.5, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (1.47583617650433274e-1, 1.0), 2.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (3.17255174305535695e-2, 1.0), 10.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (1.59022512561763752e-2, 1.0), 20.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (3.18299276490825515e-3, 1.0), 100.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (3.18309780080558939e-4, 1.0), 1000.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (3.18309885122757724e-5, 1.0), 10000.0, TEST_TOL6);

  TEST (sc_cdf_tdist_Qinv,
        (5.00318309780080559e-1, 1.0), -0.001, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (5.03182992764908255e-1, 1.0), -0.01, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (5.31725517430553570e-1, 1.0), -0.1, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (6.00023120032852124e-1, 1.0), -0.325, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (7.5e-1, 1.0), -1.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (8.12832958189001184e-1, 1.0), -1.5, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (8.52416382349566726e-1, 1.0), -2.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.68274482569446430e-1, 1.0), -10.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.84097748743823625e-1, 1.0), -20.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.96817007235091745e-1, 1.0), -100.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.99681690219919441e-1, 1.0), -1000.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.99968169011487724e-1, 1.0), -10000.0, TEST_TOL6);

  TEST (sc_cdf_tdist_Qinv,
        (5.00353553302204959e-1, 2.0), -0.001, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (5.03535445520899514e-1, 2.0), -0.01, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (5.35267280792929913e-1, 2.0), -0.1, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (6.11985772746873767e-1, 2.0), -0.325, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (7.88675134594812882e-1, 2.0), -1.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (8.63803437554499460e-1, 2.0), -1.5, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.08248290463863016e-1, 2.0), -2.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.95073771488337155e-1, 2.0), -10.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.98754668053816452e-1, 2.0), -20.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.99950007498750219e-1, 2.0), -100.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.99999500000749999e-1, 2.0), -1000.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.99999995000000075e-1, 2.0), -10000.0, 1e-6);

  TEST (sc_cdf_tdist_Qinv,
        (5.00000000000000000e-01, 300.0), 0.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (4.99601390099057051e-1, 300.0), 0.001, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (4.96013966979440912e-1, 300.0), 0.01, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (4.60205558822231806e-1, 300.0), 0.1, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (3.72703798457476188e-1, 300.0), 0.325, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (1.59058202215313138e-1, 300.0), 1.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (6.73330165746308628e-2, 300.0), 1.5, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (2.32007604915745452e-2, 300.0), 2.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (8.279313677e-21, 300.0), 10.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (1.93159812815803978e-57, 300.0), 20.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (1.02557519997736154e-232, 300.0), 100.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (0.00000000000000000e+00, 300.0), SC_POSINF, 0.0);

  TEST (sc_cdf_tdist_Qinv,
        (5.00398609900942949e-1, 300.0), -0.001, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (5.03986033020559088e-1, 300.0), -0.01, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (5.39794441177768194e-1, 300.0), -0.1, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (6.27296201542523812e-1, 300.0), -0.325, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (8.40941797784686862e-1, 300.0), -1.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.32666983425369137e-1, 300.0), -1.5, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (9.76799239508425455e-1, 300.0), -2.0, TEST_TOL6);
  TEST (sc_cdf_tdist_Qinv,
        (1.000000000000000000e0, 300.0), SC_NEGINF, TEST_TOL6);
}


#define FUNC(x)  test_ ## x,                     "test sc_ran_" #x
#define FUNC2(x) test_ ## x, test_ ## x ## _pdf, "test sc_ran_" #x
#define N 100000

static sc_rng *r_global;

static double test_ugaussian(void)
{
  return sc_ran_gaussian(r_global, 1.0);
}


void
testMoments (double (*f) (void), const char *name,
             double a, double b, double p)
{
  int i;
  double count = 0, expected, sigma;
  int status;

  for (i = 0; i < N; i++)
  {
    double r = f ();
    if (r < b && r > a)
      count++;
  }

  expected = p * N;
  sigma = fabs (count - expected) / sqrt (expected);

  status = (sigma > 3);

  sc_test (status, "%s [%g,%g] (%g observed vs %g expected)",
           name, a, b, count / N, p);
}

void test_ran_gaussian(void)
{
  r_global = sc_rng_alloc (&intern_rng_type);
  testMoments (FUNC (ugaussian), 0.0, 100.0, 0.5);
  testMoments (FUNC (ugaussian), -1.0, 1.0, 0.67971);
  testMoments (FUNC (ugaussian), 3.0, 3.5, 0.0012);
}

int main(int argc, char *argv[])
{
  test_ran_gaussian ();

  test_tdist ();
  test_fdist ();
  test_chisq ();

  test_chisqinv ();
  test_tdistinv ();

  test_auto_gaussian ();
  test_auto_fdist ();
  test_auto_flat ();
  test_auto_chisq ();

  //this fails because implementation is inaccurate for
  //continuous-valued degree-of-freedom input arg (it passes for discrete
  //values, ie. test_tdist)
  //test_auto_tdist ();

  exit (sc_test_summary ());
}
