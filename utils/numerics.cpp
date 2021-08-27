/**
 * @brief Wrappers for core math routines from open-sources: VNL and CEPHES.
 */
/*
 * Original Author: Dennis Jen and Silvester Czanner
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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#define export // obsolete feature "export template" used in these header files

#include <vnl/vnl_det.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>

#include <vnl/algo/vnl_cholesky.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_nonlinear_minimizer.h>

#undef export

#include "fs_vnl/fs_cost_function.h"
#include "fs_vnl/fs_lbfgs.h"
#include "fs_vnl/fs_powell.h"

#include "cephes.h"
#include "diag.h"
#include "error.h"
#include "numerics.h"

// used for calculating the incomplete beta function
static double d_huge();
static double d_epsilon();
static double gamma_log(double x);
static double beta(double x, double y);

// used for calculating the incomplete gamma function
static double normal_01_cdf(double x);
static double d_min(double x, double y);

// this is needed by mrimorph to essentially subclass
void (*user_call_func)(float[]) = NULL;

static void ConvertFromFloatToVNLDouble(const float *iaSource,
                                        vnl_vector< double > &oDestination,
                                        const int icSourceAndDestination)
{
  // indexing p at 1 because of NR legacy
  for (int nDest = 0; nDest < icSourceAndDestination; nDest++) {
    oDestination(nDest) = static_cast< double >(iaSource[nDest + 1]);
  }
}

static void ConvertFromVNLDoubleToFloat(const vnl_vector< double > iSource,
                                        float *oaDestination,
                                        const int icSourceAndDestination)
{
  // indexing p at 1 because of NR legacy
  for (int i = 0; i < icSourceAndDestination; i++) {
    oaDestination[i + 1] = static_cast< float >(iSource(i));
  }
}

static void ConvertFromFloatToVNLDouble(float **iaSource,
                                        vnl_matrix< double > &oDestination,
                                        const int icSourceAndDestination)
{
  // indexing p at 1 because of NR legacy
  for (int row = 0; row < icSourceAndDestination; row++) {
    for (int col = 0; col < icSourceAndDestination; col++) {
      oDestination(row, col) = static_cast< double >(iaSource[row + 1][col + 1]);
    }
  }
}

static void ConvertFromVNLDoubleToFloat(const vnl_matrix< double > iSource,
                                        float **oaDestination,
                                        const int icSourceAndDestination)
{
  // indexing p at 1 because of NR legacy
  for (int row = 0; row < icSourceAndDestination; row++) {
    for (int col = 0; col < icSourceAndDestination; col++) {
      oaDestination[row + 1][col + 1] = static_cast< float >(iSource(row, col));
    }
  }
}

/**
 * D_HUGE returns a "huge" real value, usually the largest legal real.
 * Discussion:
 *
 * HUGE_VAL is the largest representable legal real number, and is usually
 * defined in math.h, or sometimes in stdlib.h.
 *
 * Original Author: John Burkardt
 * http://www.scs.fsu.edu/~burkardt/cpp_src/prob/prob.html
 *
 * @return a "huge" real value.
 */
static double d_huge() { return HUGE_VAL; }

/**
 * D_EPSILON returns the round off unit for floating arithmetic.
 *
 * Discussion:
 *
 *    D_EPSILON is a number R which is a power of 2 with the property that,
 *    to the precision of the computer's arithmetic,
 *      1 < 1 + R
 *    but
 *      1 = ( 1 + R / 2 )
 *
 * Original Author: John Burkardt
 * http://www.scs.fsu.edu/~burkardt/cpp_src/prob/prob.html
 *
 * @return floating point round-off unit.
 */
static double d_epsilon()
{
  double r;

  r = 1.0;

  while (1.0 < (double)(1.0 + r)) {
    r = r / 2.0;
  }

  return (2.0 * r);
}

/**
 * GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
 *
 * Discussion:
 *   The program uses rational functions that theoretically approximate
 *   LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
 *   approximation for 12 < X is from Hart et al, while approximations
 *   for X < 12.0 are similar to those in Cody and Hillstrom, but are
 *
 *   The accuracy achieved depends on the arithmetic system, the compiler,
 *   intrinsic functions, and proper selection of the machine-dependent
 *   constants.
 *
 * Original Author:
 *   William Cody and L. Stoltz
 *   Argonne National Laboratory
 *   C++ translation by John Burkardt.
 *   http://www.scs.fsu.edu/~burkardt/cpp_src/prob/prob.html
 *
 * Reference:
 *   William Cody and Kenneth Hillstrom,
 *   Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
 *   Mathematics of Computation,
 *   Volume 21, 1967, pages 198-203.
 *
 *   Kenneth Hillstrom,
 *   ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
 *   May 1969.
 *
 *   Hart, Ward Cheney, Charles Lawson, Maehly,
 *   Charles Mesztenyi, John Rice, Thacher, Witzgall,
 *   Computer Approximations,
 *   Wiley and sons, New York, 1968.
 *
 * Machine-dependent constants:
 * BETA   - radix for the floating-point representation.
 * MAXEXP - the smallest positive power of BETA that overflows.
 * XBIG   - largest argument for which LN(GAMMA(X)) is representable
 *          in the machine, i.e., the solution to the equation
 *            LN(GAMMA(XBIG)) = BETA**MAXEXP.
 * FRTBIG - Rough estimate of the fourth root of XBIG
 *
 * Approximate values for some important machines are:
 *                           BETA      MAXEXP         XBIG
 * CRAY-1        (S.P.)        2        8191       9.62E+2461
 * Cyber 180/855
 *   under NOS   (S.P.)        2        1070       1.72E+319
 * IEEE (IBM/XT,
 *   SUN, etc.)  (S.P.)        2         128       4.08E+36
 * IEEE (IBM/XT,
 *   SUN, etc.)  (D.P.)        2        1024       2.55D+305
 * IBM 3033      (D.P.)       16          63       4.29D+73
 * VAX D-Format  (D.P.)        2         127       2.05D+36
 * VAX G-Format  (D.P.)        2        1023       1.28D+305
 *
 *                          FRTBIG
 *
 * CRAY-1        (S.P.)   3.13E+615
 * Cyber 180/855
 *   under NOS   (S.P.)   6.44E+79
 * IEEE (IBM/XT,
 *   SUN, etc.)  (S.P.)   1.42E+9
 * IEEE (IBM/XT,
 *   SUN, etc.)  (D.P.)   2.25D+76
 * IBM 3033      (D.P.)   2.56D+18
 * VAX D-Format  (D.P.)   1.20D+9
 * VAX G-Format  (D.P.)   1.89D+76
 *
 * @param X The argument of the Gamma function.  X must be positive.
 *
 * @return The logarithm of the Gamma function of X.
 *   If X <= 0.0, or if overflow would occur, the program returns the
 *   value XINF, the largest representable floating point number.
 */
static double gamma_log(double x)
{
  double c[7] = {-1.910444077728E-03,
                 8.4171387781295E-04,
                 -5.952379913043012E-04,
                 7.93650793500350248E-04,
                 -2.777777777777681622553E-03,
                 8.333333333333333331554247E-02,
                 5.7083835261E-03};
  double corr;
  double d1 = -5.772156649015328605195174E-01;
  double d2 = 4.227843350984671393993777E-01;
  double d4 = 1.791759469228055000094023;
  double frtbig = 1.42E+09;
  int i;
  double p1[8] = {4.945235359296727046734888,
                  2.018112620856775083915565E+02,
                  2.290838373831346393026739E+03,
                  1.131967205903380828685045E+04,
                  2.855724635671635335736389E+04,
                  3.848496228443793359990269E+04,
                  2.637748787624195437963534E+04,
                  7.225813979700288197698961E+03};
  double p2[8] = {4.974607845568932035012064,
                  5.424138599891070494101986E+02,
                  1.550693864978364947665077E+04,
                  1.847932904445632425417223E+05,
                  1.088204769468828767498470E+06,
                  3.338152967987029735917223E+06,
                  5.106661678927352456275255E+06,
                  3.074109054850539556250927E+06};
  double p4[8] = {1.474502166059939948905062E+04,
                  2.426813369486704502836312E+06,
                  1.214755574045093227939592E+08,
                  2.663432449630976949898078E+09,
                  2.940378956634553899906876E+010,
                  1.702665737765398868392998E+011,
                  4.926125793377430887588120E+011,
                  5.606251856223951465078242E+011};
  double pnt68 = 0.6796875;
  double q1[8] = {6.748212550303777196073036E+01,
                  1.113332393857199323513008E+03,
                  7.738757056935398733233834E+03,
                  2.763987074403340708898585E+04,
                  5.499310206226157329794414E+04,
                  6.161122180066002127833352E+04,
                  3.635127591501940507276287E+04,
                  8.785536302431013170870835E+03};
  double q2[8] = {1.830328399370592604055942E+02,
                  7.765049321445005871323047E+03,
                  1.331903827966074194402448E+05,
                  1.136705821321969608938755E+06,
                  5.267964117437946917577538E+06,
                  1.346701454311101692290052E+07,
                  1.782736530353274213975932E+07,
                  9.533095591844353613395747E+06};
  double q4[8] = {2.690530175870899333379843E+03,
                  6.393885654300092398984238E+05,
                  4.135599930241388052042842E+07,
                  1.120872109616147941376570E+09,
                  1.488613728678813811542398E+010,
                  1.016803586272438228077304E+011,
                  3.417476345507377132798597E+011,
                  4.463158187419713286462081E+011};
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double xbig = 4.08E+36;
  double xden;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double xsq;
  //
  //  Return immediately if the argument is out of range.
  //
  if (x <= 0.0 || xbig < x) {
    return d_huge();
  }

  if (x <= d_epsilon()) {
    res = -log(x);
  }
  else if (x <= 1.5) {
    if (x < pnt68) {
      corr = -log(x);
      xm1 = x;
    }
    else {
      corr = 0.0;
      xm1 = (x - 0.5) - 0.5;
    }

    if (x <= 0.5 || pnt68 <= x) {
      xden = 1.0;
      xnum = 0.0;

      for (i = 0; i < 8; i++) {
        xnum = xnum * xm1 + p1[i];
        xden = xden * xm1 + q1[i];
      }

      res = corr + (xm1 * (d1 + xm1 * (xnum / xden)));
    }
    else {
      xm2 = (x - 0.5) - 0.5;
      xden = 1.0;
      xnum = 0.0;
      for (i = 0; i < 8; i++) {
        xnum = xnum * xm2 + p2[i];
        xden = xden * xm2 + q2[i];
      }

      res = corr + xm2 * (d2 + xm2 * (xnum / xden));
    }
  }
  else if (x <= 4.0) {
    xm2 = x - 2.0;
    xden = 1.0;
    xnum = 0.0;
    for (i = 0; i < 8; i++) {
      xnum = xnum * xm2 + p2[i];
      xden = xden * xm2 + q2[i];
    }

    res = xm2 * (d2 + xm2 * (xnum / xden));
  }
  else if (x <= 12.0) {
    xm4 = x - 4.0;
    xden = -1.0;
    xnum = 0.0;
    for (i = 0; i < 8; i++) {
      xnum = xnum * xm4 + p4[i];
      xden = xden * xm4 + q4[i];
    }

    res = d4 + xm4 * (xnum / xden);
  }
  else {
    res = 0.0;

    if (x <= frtbig) {
      res = c[6];
      xsq = x * x;

      for (i = 0; i < 6; i++) {
        res = res / xsq + c[i];
      }
    }

    res = res / x;
    corr = log(x);
    res = res + sqrtpi - 0.5 * corr;
    res = res + x * (corr - 1.0);
  }

  return res;
}

/**
 * BETA returns the value of the Beta function.
 *
 * Formula:
 *
 *   BETA(X,Y) = ( GAMMA(X) * GAMMA(Y) ) / GAMMA(X+Y)
 *
 * Properties:
 *
 *    BETA(X,Y) = BETA(Y,X).
 *    BETA(X,Y) = Integral ( 0 <= T <= 1 ) T**(X-1) (1-T)**(Y-1) dT.
 *    BETA(X,Y) = GAMMA(X) * GAMMA(Y) / GAMMA(X+Y)
 *
 * Restrictions:
 *
 *   Both X and Y must be greater than 0.
 *
 * Original Author: John Burkardt
 * http: *www.scs.fsu.edu/~burkardt/cpp_src/prob/prob.html
 *
 * @param x Parameter that define the Beta function.  Must be greater than 0.
 * @param y Parameter that define the Beta function.  Must be greater than 0.
 * @return Value of the Beta function.
 */
static double beta(double x, double y)
{
  if (x <= 0.0 || y <= 0.0) {
    std::cout << "\n";
    std::cout << "BETA - Fatal error!\n";
    std::cout << "  Both X and Y must be greater than 0.\n";
    exit(1);
  }

  return (exp(gamma_log(x) + gamma_log(y) - gamma_log(x + y)));
}

/**
 * BETA_INC returns the value of the incomplete Beta function.
 *
 *  Discussion:
 *
 *    This calculation requires an iteration.  In some cases, the iteration
 *    may not converge rapidly, or may become inaccurate.
 *
 *  Formula:
 *
 *    BETA_INC(A,B,X)
 *
 *      =   Integral ( 0 <= T <= X ) T**(A-1) (1-T)**(B-1) dT
 *        / Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT
 *
 *      =   Integral ( 0 <= T <= X ) T**(A-1) (1-T)**(B-1) dT
 *        / BETA(A,B)
 *
 *  Original Author: Majumder and Bhattacharjee.
 *
 *    C++ translation by John Burkardt.
 *    http://www.scs.fsu.edu/~burkardt/cpp_src/prob/prob.html
 *
 *  Reference:
 *
 *    Majumder and Bhattacharjee,
 *    Algorithm AS63,
 *    Applied Statistics,
 *    1973, volume 22, number 3.
 *
 * @param a Parameter of the function, 0.0 < a.
 * @param b Parameter of the function, 0.0 < b.
 * @param x Argument of the function.  Normally, 0.0 <= x <= 1.0.
 * @return Value of the function.
 */
float OpenBetaIncomplete(float a, float b, float x)
{
  double cx;
  int i;
  int it;
  int it_max = 5000;
  bool indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double tol = 1.0E-07;
  double value;
  double xx;

  if (a <= 0.0) {
    std::cout << "\n";
    std::cout << "BETA_INC - Fatal error!\n";
    std::cout << "  A <= 0.\n";
    exit(1);
  }

  if (b <= 0.0) {
    std::cout << "\n";
    std::cout << "BETA_INC - Fatal error!\n";
    std::cout << "  B <= 0.\n";
    exit(1);
  }

  if (x <= 0.0) {
    value = 0.0;
    return value;
  }
  else if (1.0 <= x) {
    value = 1.0;
    return value;
  }
  //
  //  Change tail if necessary and determine S.
  //
  psq = a + b;

  if (a < (a + b) * x) {
    xx = 1.0 - x;
    cx = x;
    pp = b;
    qq = a;
    indx = true;
  }
  else {
    xx = x;
    cx = 1.0 - x;
    pp = a;
    qq = b;
    indx = false;
  }

  term = 1.0;
  i = 1;
  value = 1.0;

  ns = (int)(qq + cx * (a + b));
  //
  //  Use Soper's reduction formulas.
  //
  rx = xx / cx;

  temp = qq - (double)i;
  if (ns == 0) {
    rx = xx;
  }

  it = 0;

  for (;;) {
    it = it + 1;

    if (it_max < it) {
      std::cout << "\n";
      std::cout << "BETA_INC - Fatal error!\n";
      std::cout << "  Maximum number of iterations exceeded!\n";
      std::cout << "  IT_MAX = " << it_max << "\n";
#if 0
      exit ( 1 );
#else
      return value;
#endif
    }

    term = term * temp * rx / (pp + (double)(i));
    value = value + term;
    temp = fabs(term);

    if (temp <= tol && temp <= tol * value) {
      break;
    }

    i = i + 1;
    ns = ns - 1;

    if (0 <= ns) {
      temp = qq - (double)i;
      if (ns == 0) {
        rx = xx;
      }
    }
    else {
      temp = psq;
      psq = psq + 1.0;
    }
  }
  //
  //  Finish calculation.
  //
  value = value * exp(pp * log(xx) + (qq - 1.0) * log(cx)) / (beta(a, b) * pp);

  if (indx) {
    value = 1.0 - value;
  }

  return value;
}

/**
 * D_MIN returns the minimum of two double precision values.
 *
 * Original Author: John Burkardt
 * http://www.scs.fsu.edu/~burkardt/cpp_src/prob/prob.html
 *
 * @param x Quantity to compare.
 * @param y Quantity to compare.
 * @return Minimum of x and y.
 */
static double d_min(double x, double y)
{
  if (y < x) {
    return y;
  }
  else {
    return x;
  }
}

/**
 * NORMAL_01_CDF evaluates the Normal 01 CDF.
 *
 * Original Author: John Burkardt
 * http://www.scs.fsu.edu/~burkardt/cpp_src/prob/prob.html
 *
 * Reference:
 *   A G Adams,
 *   Areas Under the Normal Curve,
 *   Algorithm 39,
 *   Computer j.,
 *   Volume 12, pages 197-198, 1969.
 *
 * @param X The argument of the CDF.
 * @return The value of the CDF.
 */
double normal_01_cdf(double x)
{
  double a1 = 0.398942280444;
  double a2 = 0.399903438504;
  double a3 = 5.75885480458;
  double a4 = 29.8213557808;
  double a5 = 2.62433121679;
  double a6 = 48.6959930692;
  double a7 = 5.92885724438;
  double b0 = 0.398942280385;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364;
  double b5 = 0.151679116635;
  double b6 = 5.29330324926;
  double b7 = 4.8385912808;
  double b8 = 15.1508972451;
  double b9 = 0.742380924027;
  double b10 = 30.789933034;
  double b11 = 3.99019417011;
  double cdf;
  double q;
  double y;
  //
  //  |X| <= 1.28.
  //
  if (fabs(x) <= 1.28) {
    y = 0.5 * x * x;

    q = 0.5 - fabs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))));
    //
    //  1.28 < |X| <= 12.7
    //
  }
  else if (fabs(x) <= 12.7) {
    y = 0.5 * x * x;

    q = exp(-y) * b0 /
        (fabs(x) - b1 +
         b2 / (fabs(x) + b3 + b4 / (fabs(x) - b5 + b6 / (fabs(x) + b7 - b8 / (fabs(x) + b9 + b10 / (fabs(x) + b11))))));
    //
    //  12.7 < |X|
    //
  }
  else {
    q = 0.0;
  }
  //
  //  Take account of negative X.
  //
  if (x < 0.0) {
    cdf = q;
  }
  else {
    cdf = 1.0 - q;
  }

  return cdf;
}

/**
 * GAMMA_INC computes the incomplete Gamma function.
 *
 * Definition:
 *
 *    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T**(P-1) EXP(-T) DT / GAMMA(P).
 *
 * Discussion:
 *
 *    GAMMA_INC(P,       0) = 0,
 *    GAMMA_INC(P,Infinity) = 1.
 *
 * Original Author: B L Shea
 * C++ translation by John Burkardt.
 * http://www.scs.fsu.edu/~burkardt/cpp_src/prob/prob.html
 *
 * Reference:
 *
 *    B L Shea,
 *    Chi-squared and Incomplete Gamma Integral,
 *    Algorithm AS239,
 *    Applied Statistics,
 *    Volume 37, Number 3, 1988, pages 466-473.
 *
 * @param p The exponent parameter, 0.0 < p.
 * @param x Integral limit parameter.  If x is less than or equal to 0,
 * GAMMA_INC is returned as 0.
 * @return The value of the function.
 */
float OpenGammaIncomplete(float p, float x)
{
  double a;
  double arg;
  double b;
  double c;
  double exp_arg_min = -88.0;
  double overflow = 1.0E+37;
  double plimit = 1000.0;
  double pn1;
  double pn2;
  double pn3;
  double pn4;
  double pn5;
  double pn6;
  double rn;
  double value;
  double tol = 1.0E-07;
  double xbig = 1.0E+08;

  value = 0.0;

  if (p <= 0.0) {
    std::cout << " \n";
    std::cout << "GAMMA_INC - Fatal error!\n";
    std::cout << "  Parameter P <= 0.\n";
    exit(1);
  }

  if (x <= 0.0) {
    value = 0.0;
    return value;
  }
  //
  //  Use a normal approximation if PLIMIT < P.
  //
  if (plimit < p) {
    pn1 = 3.0 * sqrt(p) * (pow(x / p, 1.0 / 3.0) + 1.0 / (9.0 * p) - 1.0);
    value = normal_01_cdf(pn1);
    return value;
  }
  //
  //  Is X extremely large compared to P?
  //
  if (xbig < x) {
    value = 1.0;
    return value;
  }
  //
  //  Use Pearson's series expansion.
  //  (P is not large enough to force overflow in the log of Gamma.
  //
  if (x <= 1.0 || x < p) {
    arg = p * log(x) - x - gamma_log(p + 1.0);
    c = 1.0;
    value = 1.0;
    a = p;

    for (;;) {
      a = a + 1.0;
      c = c * x / a;
      value = value + c;

      if (c <= tol) {
        break;
      }
    }

    arg = arg + log(value);

    if (exp_arg_min <= arg) {
      value = exp(arg);
    }
    else {
      value = 0.0;
    }
  }
  //
  //  Use a continued fraction expansion.
  //
  else {
    arg = p * log(x) - x - gamma_log(p);
    a = 1.0 - p;
    b = a + x + 1.0;
    c = 0.0;
    pn1 = 1.0;
    pn2 = x;
    pn3 = x + 1.0;
    pn4 = x * b;
    value = pn3 / pn4;

    for (;;) {
      a = a + 1.0;
      b = b + 2.0;
      c = c + 1.0;
      pn5 = b * pn3 - a * c * pn1;
      pn6 = b * pn4 - a * c * pn2;

      if (0.0 < fabs(pn6)) {
        rn = pn5 / pn6;

        if (fabs(value - rn) <= d_min(tol, tol * rn)) {
          arg = arg + log(value);

          if (exp_arg_min <= arg) {
            value = 1.0 - exp(arg);
          }
          else {
            value = 1.0;
          }

          return value;
        }
        value = rn;
      }
      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;
      //
      //  Rescale terms in continued fraction if terms are large.
      //
      if (overflow <= fabs(pn5)) {
        pn1 = pn1 / overflow;
        pn2 = pn2 / overflow;
        pn3 = pn3 / overflow;
        pn4 = pn4 / overflow;
      }
    }
  }

  return value;
}

int OpenDFPMin(float p[],
                          int n,
                          float iTolerance,
                          int *oIterations,
                          float *oFinalFunctionReturn,
                          float (*iFunction)(float[]),
                          void (*iDerivativeFunction)(float[], float[]),
                          void (*iStepFunction)(int, float, void *, float *),
                          void *iStepFunctionParams,
                          void (*iUserCallBackFunction)(float[]))
{
  int returnCode;
  fs_cost_function costFunction(iFunction, iDerivativeFunction, n);

  vnl_vector< double > finalParameters(n);
  ConvertFromFloatToVNLDouble(p, finalParameters, n);

  fs_lbfgs minimizer(costFunction);

  fs_lbfgs_observer observer;
  observer.setStepFunction(iStepFunction, iStepFunctionParams);
  observer.setUserCallbackFunction(iUserCallBackFunction);

  minimizer.setObserver(&observer);

  minimizer.memory = 7;

  // increased accuracy
  minimizer.line_search_accuracy = 0.1;

  const double functionTolerance = 0.001;
  minimizer.set_x_tolerance(functionTolerance);
  minimizer.set_f_tolerance(0.1);

  // I had to set the gradient convergence tolerance or else
  //       it refused to converge on the second test -- this  impacts the
  //       convergence!
  minimizer.set_g_tolerance(static_cast< double >(iTolerance));

  // keep iterating if the error isn't significantly better
  bool isSuccess = true;
  float previousValue = *oFinalFunctionReturn;
  bool shouldContinue = true;
  while (shouldContinue && isSuccess) {
    isSuccess = minimizer.minimize(finalParameters);
    float currentValue = minimizer.get_end_error();

    if ((previousValue - currentValue) / previousValue <= iTolerance) {
      shouldContinue = false;
    }
  }

  if (isSuccess) {
    // success
    *oIterations = observer.getNumberOfOptimalUpdates();
    *oFinalFunctionReturn = minimizer.get_end_error();
    ConvertFromVNLDoubleToFloat(finalParameters, p, n);
    returnCode = 0;
  }
  else {
    returnCode = minimizer.get_failure_code();

    if (returnCode == vnl_nonlinear_minimizer::ERROR_FAILURE) {
      *oIterations = observer.getNumberOfOptimalUpdates();
      *oFinalFunctionReturn = minimizer.get_end_error();

      ConvertFromVNLDoubleToFloat(finalParameters, p, n);
    }
    else if (returnCode == vnl_nonlinear_minimizer::ERROR_DODGY_INPUT) {
      ErrorExit(ERROR_BADPARM, "quasi-Newton minimization (lbfsg) dodgy input");
    }
    else if (returnCode == vnl_nonlinear_minimizer::FAILED_TOO_MANY_ITERATIONS ||
             returnCode == vnl_nonlinear_minimizer::FAILED_FTOL_TOO_SMALL ||
             returnCode == vnl_nonlinear_minimizer::FAILED_XTOL_TOO_SMALL ||
             returnCode == vnl_nonlinear_minimizer::FAILED_GTOL_TOO_SMALL) {
      ErrorExit(ERROR_BADPARM, "quasi-Newton minimization (lbfsg) tolerances too small");
    }
  }
  return (returnCode);
}

/**
 * Provides the eigen values and vectors for symmetric matrices.
 */
int OpenEigenSystem(float *iaData, int icData, float *oEigenValues, float *oEigenVectors)
{
  vnl_matrix< float > vnlMatrix(iaData, icData, icData);

  vnl_symmetric_eigensystem< float > eigenSystem(vnlMatrix);

  eigenSystem.V.copy_out(oEigenVectors);
  for (int i = 0; i < icData; i++) {
    oEigenValues[i] = eigenSystem.get_eigenvalue(i);
  }

  // check for errors?
  return NO_ERROR;
}

/**
 * Provides the eigen values and vectors for symmetric matrices.
 */
int OpenNonSymmetricEigenSystem(float *iaData, int icData, float *oEigenValues, float *oEigenVectors)
{
  // convert the data into a double
  double data[icData * icData];
  for (int row = 0; row < icData; row++) {
    for (int column = 0; column < icData; column++) {
      int index = column + row * icData;
      data[index] = iaData[index];
    }
  }

  vnl_matrix< double > vnlMatrix(data, icData, icData);

  // this does the actual eigensystem calculation
  vnl_real_eigensystem eigenSystem(vnlMatrix);

  // copy the eigen vectors (doubles) out into our float array
  for (int row = 0; row < icData; row++) {
    for (int column = 0; column < icData; column++) {
      int index = column + row * icData;
      oEigenVectors[index] = eigenSystem.Vreal[row][column];
    }
  }

  for (int row = 0; row < icData; row++) {
    oEigenValues[row] = eigenSystem.D[row].real();
  }

  // check for errors?
  return NO_ERROR;
}

void OpenPowell(float iaParams[],
                           float **ioInitialDirection,
                           int icParams,
                           float iTolerance,
                           int *oIterations,
                           float *oFinalFunctionReturn,
                           float (*iFunction)(float[]))
{
  fs_cost_function costFunction(iFunction);
  fs_powell minimizer(&costFunction);

  minimizer.set_x_tolerance(static_cast< double >(iTolerance));

  vnl_vector< double > finalParameters(icParams);
  ConvertFromFloatToVNLDouble(iaParams, finalParameters, icParams);

  vnl_matrix< double > initialDirection(icParams, icParams, vnl_matrix_identity);

  if (ioInitialDirection != NULL) {
    ConvertFromFloatToVNLDouble(ioInitialDirection, initialDirection, icParams);
  }

  int returnCode = minimizer.minimize(finalParameters, &initialDirection);

  // exit if failure
  if (returnCode == vnl_nonlinear_minimizer::FAILED_TOO_MANY_ITERATIONS) {
    ErrorExit(ERROR_BADPARM, "powell exceeding maximum iterations.");
  }
  else if (returnCode == fs_powell::ERROR_FAILURE || returnCode == fs_powell::ERROR_DODGY_INPUT ||
           returnCode == fs_powell::FAILED_FTOL_TOO_SMALL || returnCode == fs_powell::FAILED_XTOL_TOO_SMALL ||
           returnCode == fs_powell::FAILED_GTOL_TOO_SMALL) {
    ErrorExit(ERROR_BADPARM, "powell error.");
  }
  else {
    // success
    *oIterations = minimizer.get_num_iterations();
    *oFinalFunctionReturn = minimizer.get_end_error();
    ConvertFromVNLDoubleToFloat(finalParameters, iaParams, icParams);
    ConvertFromVNLDoubleToFloat(initialDirection, ioInitialDirection, icParams);
  }
}

/*------------------------------------------------------------------------
  OpenPowell2() - this is a mod of OpenPowell() above. There are several
  changes:
    1. It sets the ftol tolerance instead of the xtol. ftol
       is the fraction of the cost that the difference in successive costs
      must drop below to stop.
    2. You can set the max number iterations.
    3. If Gdiag > 0, then prints out info as it optimizes (setenv DIAG_NO 1)
    4. Returns 0 if no error, or 1 if error. The other version just exited.
    5. Sets all params even if an error is returned.
    6. Allows user to set the tolerance on the 1D Min
  I made a new function because BF uses OpenPowell() in a lot of places.
  It would be better to have more options on this function.
  Note: each "iteration" is a loop thru a 1D min for each parameter.
  ------------------------------------------------------------------------*/
int OpenPowell2(float iaParams[],
                           float **ioInitialDirection,
                           int icParams,
                           float iTolerance,
                           float iLinMinTol,
                           int MaxIterations,
                           int *oIterations,
                           float *oFinalFunctionReturn,
                           float (*iFunction)(float[]))
{
  fs_cost_function costFunction(iFunction);
  fs_powell minimizer(&costFunction);

  vnl_vector< double > finalParameters(icParams);
  ConvertFromFloatToVNLDouble(iaParams, finalParameters, icParams);

  vnl_matrix< double > initialDirection(icParams, icParams, vnl_matrix_identity);

  if (ioInitialDirection != NULL) {
    ConvertFromFloatToVNLDouble(ioInitialDirection, initialDirection, icParams);
  }

  if (Gdiag_no > 0) {
    minimizer.set_trace(1);
    minimizer.set_verbose(1);
  }
  minimizer.set_linmin_xtol(iLinMinTol);
  // minimizer.set_x_tolerance(iLinMinTol);

  minimizer.set_f_tolerance(iTolerance);
  minimizer.set_max_function_evals(MaxIterations);

  int returnCode = minimizer.minimize(finalParameters, &initialDirection);

  *oIterations = minimizer.get_num_iterations();
  *oFinalFunctionReturn = minimizer.get_end_error();
  ConvertFromVNLDoubleToFloat(finalParameters, iaParams, icParams);
  ConvertFromVNLDoubleToFloat(initialDirection, ioInitialDirection, icParams);

  if (returnCode == vnl_nonlinear_minimizer::FAILED_TOO_MANY_ITERATIONS) {
    printf("powell exceeded maximum iterations\n");
    return (1);
  }
  else if (returnCode == fs_powell::ERROR_FAILURE || returnCode == fs_powell::ERROR_DODGY_INPUT ||
           returnCode == fs_powell::FAILED_FTOL_TOO_SMALL || returnCode == fs_powell::FAILED_XTOL_TOO_SMALL ||
           returnCode == fs_powell::FAILED_GTOL_TOO_SMALL) {
    printf("powell error %d\n", returnCode);
    return (1);
  }

  // success
  return (0);
}

int OpenLUMatrixInverse(MATRIX *iMatrix, MATRIX *oInverse)
{
  // NO_ERROR from error.h
  int errorCode = NO_ERROR;

  vnl_matrix< float > vnlMatrix(iMatrix->data, iMatrix->rows, iMatrix->cols);

  unsigned int r = vnlMatrix.rows();
  if (r <= 4 && r == vnlMatrix.cols()) {
    if (r == 1) {
      if (vnlMatrix(0, 0) == 0.0)
        errorCode = ERROR_BADPARM;
      else
        oInverse->data[0] = 1.0 / vnlMatrix(0, 0);
    }
    else if (r == 2) {
      vnl_matrix_fixed< float, 2, 2 > m(vnlMatrix);
      if (vnl_det(m) == 0.0)
        errorCode = ERROR_BADPARM;
      else
        vnl_inverse(m).copy_out(oInverse->data);
    }
    else if (r == 3) {
      vnl_matrix_fixed< float, 3, 3 > m(vnlMatrix);
      if (vnl_det(m) == 0.0)
        errorCode = ERROR_BADPARM;
      else
        vnl_inverse(m).copy_out(oInverse->data);
    }
    else {
      vnl_matrix_fixed< float, 4, 4 > m(vnlMatrix);
      if (vnl_det(m) == 0.0)
        errorCode = ERROR_BADPARM;
      else
        vnl_inverse(m).copy_out(oInverse->data);
    }
  }

  else  // > 4x4 matrices
  {
    // the svd matrix inversion failed a test case, whereas qr passes, so we're
    // going to use the qr generated inverse
    vnl_qr< float > vnlMatrixInverter(vnlMatrix);

    // determinant of 0 means that it's singular
    if (vnlMatrixInverter.determinant() != 0.0) {
      vnl_matrix< float > inverse = vnlMatrixInverter.inverse();
      inverse.copy_out(oInverse->data);
    }
    else {
      errorCode = ERROR_BADPARM;
    }
  }

  return errorCode;
}

/*!
  \fn int OpenQRdecomposition(const MATRIX *iMatrix, MATRIX *oQ, MATRIX *oR)
  \brief Performs QR decomposition on the input matrix. The output
  matricies must be allocated already (Q is nrowsxnows where nrows is
  the number of rows in the input matrix, and R is same size as the
  input matrix). 
 */
int OpenQRdecomposition(const MATRIX *iMatrix, MATRIX *oQ, MATRIX *oR)
{
  if(iMatrix==NULL){
    printf("QRdecomposition(): input matrix is NULL\n");
    return(1);
  }
  if(oQ==NULL){
    printf("QRdecomposition(): output Q matrix is NULL\n");
    return(1);
  }
  if(oR==NULL){
    printf("QRdecomposition(): output R matrix is NULL\n");
    return(1);
  }
  if(oQ->rows != iMatrix->rows){
    printf("QRdecomposition(): Q rows != input Matrix rows\n");
    return(1);
  }
  if(oQ->cols != oQ->rows){
    printf("QRdecomposition(): Q cols != Q rows\n");
    return(1);
  }
  if(oR->rows != iMatrix->rows){
    printf("QRdecomposition(): R rows != input Matrix rows\n");
    return(1);
  }
  if(oR->cols != iMatrix->cols){
    printf("QRdecomposition(): R cols != input Matrix cols\n");
    return(1);
  }

  vnl_matrix< float >  vnlMatrix(iMatrix->data, iMatrix->rows, iMatrix->cols);
  vnl_qr<float> qr( vnlMatrix );
  qr.Q().copy_out(oQ->data);
  qr.R().copy_out(oR->data);
  return(0);
}

/**
 * Returns the determinant of matrix.
 * @param iMatrix
 * @return Returns 0 if the matrix is non-square.
 */
float OpenMatrixDeterminant(MATRIX *iMatrix)
{
  float determinant = 0.0;

  if (iMatrix->rows == iMatrix->cols) {
    vnl_matrix< float > vnlMatrix(iMatrix->data, iMatrix->rows, iMatrix->cols);
    if (vnlMatrix.rows() == 1)
      determinant = vnlMatrix(0, 0);
    else if (vnlMatrix.rows() == 2)
      determinant = vnl_det(vnl_matrix_fixed< float, 2, 2 >(vnlMatrix));
    else if (vnlMatrix.rows() == 3)
      determinant = vnl_det(vnl_matrix_fixed< float, 3, 3 >(vnlMatrix));
    else if (vnlMatrix.rows() == 4)
      determinant = vnl_det(vnl_matrix_fixed< float, 4, 4 >(vnlMatrix));
    else
      determinant = vnl_determinant< float >(vnlMatrix);
  }

  return determinant;
}

/**
 * Perform an SV Decomposition.
 * @param ioA Input matrix and output value of u.
 * @param oW Ouput diagonal vector.
 * @param oV Ouput V matrix.
 */
int OpenSvdcmp(MATRIX *ioA, VECTOR *oW, MATRIX *oV)
{
  int errorCode = NO_ERROR;

  if (ioA->rows >= ioA->cols) {
    vnl_matrix< float > vnlX(ioA->data, ioA->rows, ioA->cols);
    vnl_svd< float > svdMatrix(vnlX);

    if (svdMatrix.valid()) {
      svdMatrix.U().copy_out(ioA->data);
      svdMatrix.W().diagonal().copy_out(oW->data);
      svdMatrix.V().copy_out(oV->data);
    }
    else {
      errorCode = ERROR_BADPARM;
    }
  }
  else {
    errorCode = ERROR_BADPARM;
  }

  return errorCode;
}

/**
 * Generates a random number between 0 and 1.
 * The sequence of numbers generated
 * for a given seed will be the same.  If you want a different sequence, use a
 * different seed.  The behaviour of this function is meant to mimic that of
 * the ran1 algorithm in numerical recipes.
 */
float OpenRan1(long *iSeed)
{

#ifdef HAVE_OPENMP
  if (omp_get_thread_num() != 0) {
    fprintf(stderr, "%s:%d OpenRan1 called from non-0 thread but this is not conducive to reproducible behavior, nor is this code thread-safe\n", __FILE__, __LINE__);
    exit(1);
  }
#endif
  
  static const double MIN = 0.0;
  static const double MAX = 1.0;

  static long mSeed;
  static vnl_random mVnlRandom(mSeed);

  if (mSeed != *iSeed) {
    mSeed = *iSeed;
    mVnlRandom.reseed(mSeed);
  }

  float randomNumber = mVnlRandom.drand64(MIN, MAX);

  return randomNumber;
}

/**
 * Generates the second derivatives needed by splint.
 * @param iYStartDerivative The derivative at the beginning of the function.
 * @param iYEndDerivative The derivative at the last point of the function.
 */
void OpenSpline(
    float iaX[], float iaY[], int icXY, float iYStartDerivative, float iYEndDerivative, float oaYSecondDerivatives[])
{
  float *secondDerivatives = SplineCubicSet(
      icXY, iaX, iaY, SPLINE_USE_FIRST_DERIVATIVE, iYStartDerivative, SPLINE_USE_FIRST_DERIVATIVE, iYEndDerivative);

  for (int i = 0; i < icXY; i++) {
    oaYSecondDerivatives[i] = secondDerivatives[i];
  }
}

/**
 * Cubic spline interpolation.  This is to be used in conjunction with spline,
 * which generates the second derivatives needed by this function.  Be aware
 * that areas outside of where the control points are defined will be
 * interpolated as the closest edge point, rather than the usual extrapolated
 * point.
 */
void OpenSplint(
    float iaX[], float iaY[], float iaYSecondDerivatives[], int icYSecondDerivatives, float iX, float *oY)
{
  float firstDerivative = 0.0;
  float secondDerivative = 0.0;

  *oY = SplineCubicValue(icYSecondDerivatives, iaX, iX, iaY, iaYSecondDerivatives, &firstDerivative, &secondDerivative);
}

/**
 * Computes the second derivatives of a piecewise cubic spline.
 *
 * For data interpolation, the user must call SPLINE_SET to determine
 * the second derivative data, passing in the data to be interpolated,
 * and the desired boundary conditions.
 *
 * The data to be interpolated, plus the SPLINE_SET output, defines
 * the spline.  The user may then call SPLINE_VAL to evaluate the
 * spline at any point.

 * The cubic spline is a piecewise cubic polynomial.  The intervals
 * are determined by the "knots" or abscissas of the data to be
 * interpolated.  The cubic spline has continous first and second
 * derivatives over the entire interval of interpolation.
 *
 * For any point T in the interval T(IVAL), T(IVAL+1), the form of
 * the spline is
 *
 *   SPL(T) = A(IVAL)
 *          + B(IVAL) * ( T - T(IVAL) )
 *          + C(IVAL) * ( T - T(IVAL) )**2
 *          + D(IVAL) * ( T - T(IVAL) )**3
 *
 * If we assume that we know the values Y(*) and YPP(*), which represent
 * the values and second derivatives of the spline at each knot, then
 * the coefficients can be computed as:
 *
 *   A(IVAL) = Y(IVAL)
 *   B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
 *     - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
 *   C(IVAL) = YPP(IVAL) / 2
 *   D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )\
 *
 * Since the first derivative of the spline is
 *
 *   SPL'(T) =     B(IVAL)
 *           + 2 * C(IVAL) * ( T - T(IVAL) )
 *           + 3 * D(IVAL) * ( T - T(IVAL) )**2,
 *
 * the requirement that the first derivative be continuous at interior
 * knot I results in a total of N-2 equations, of the form:
 *
 *   B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
 *   + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
 *
 * or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
 *
 *   ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
 *   - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
 *   + YPP(IVAL-1) * H(IVAL-1)
 *   + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
 *   =
 *   ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
 *   - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
 *
 * or
 *
 *   YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
 *   + YPP(IVAL) * H(IVAL)
 *   =
 *   6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
 *   - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
 *
 * Boundary conditions must be applied at the first and last knots.
 * The resulting tridiagonal system can be solved for the YPP values.
 *
 * Original Author: John Burkardt
 * http://www.csit.fsu.edu/~burkardt/cpp_src/spline/spline.html
 *
 * @param n The number of data points.  N must be at least 2.
 * In the special case where N = 2 and IBCBEG = IBCEND = 0, the
 * spline will actually be linear.
 *
 * @param t The knot values, size n, that is, the points were data is
 * specified.  The knot values should be distinct, and increasing.
 *
 * @param y The data values, size n, to be interpolated.
 *
 * @param ibcbeg Left boundary condition flag:
 *   SPLINE_USE_QUADRATIC=0: the cubic spline should be a quadratic over the
 *    first interval;
 *   SPLINE_USE_FIRST_DERIVATIVE=1: the first derivative at the left endpoint
 *    should be YBCBEG;
 *   SPLINE_USE_SECOND_DERIVATIVE=2: the second derivative at the
 *    left endpoint should be YBCBEG.
 *
 * @param ybcbeg The values to be used in the boundary
 * conditions if IBCBEG is equal to 1 or 2.
 *
 * @param ibcend Right boundary condition flag:
 *   SPLINE_USE_QUADRATIC=0: the cubic spline should be a quadratic over the
 *    last interval;
 *   SPLINE_USE_FIRST_DERIVATIVE=1: the first derivative at the right endpoint
 *    should be YBCEND;
 *   SPLINE_USE_SECOND_DERIVATIVE=2: the second derivative at the right
 *    endpoint should be YBCEND.
 *
 * @param ybcend The values to be used in the boundary conditions if IBCEND is
 * equal to 1 or 2.
 *
 * @return The second derivatives of the cubic spline, size n.
 */
float *SplineCubicSet(int n, float t[], float y[], int ibcbeg, float ybcbeg, int ibcend, float ybcend)
{
  float *a;
  float *b;
  int i;
  float *ypp;
  //
  //  Check.
  //
  if (n <= 1) {
    std::cout << "\n";
    std::cout << "SPLINE_CUBIC_SET - Fatal error!\n";
    std::cout << "  The number of data points N must be at least 2.\n";
    std::cout << "  The input value is " << n << ".\n";
    exit(1);
  }

  for (i = 0; i < n - 1; i++) {
    if (t[i + 1] <= t[i]) {
      std::cout << "\n";
      std::cout << "SPLINE_CUBIC_SET - Fatal error!\n";
      std::cout << "  The knots must be strictly increasing, but\n";
      std::cout << "  T(" << i << ") = " << t[i] << "\n";
      std::cout << "  T(" << i + 1 << ") = " << t[i + 1] << "\n";
      exit(1);
    }
  }

  a = new float[3 * n];
  b = new float[n];
  //
  //  Set up the first equation.
  //
  if (ibcbeg == SPLINE_USE_QUADRATIC) {
    b[0] = 0.0;
    a[1 + 0 * 3] = 1.0;
    a[0 + 1 * 3] = -1.0;
  }
  else if (ibcbeg == SPLINE_USE_FIRST_DERIVATIVE) {
    b[0] = (y[1] - y[0]) / (t[1] - t[0]) - ybcbeg;
    a[1 + 0 * 3] = (t[1] - t[0]) / 3.0;
    a[0 + 1 * 3] = (t[1] - t[0]) / 6.0;
  }
  else if (ibcbeg == SPLINE_USE_SECOND_DERIVATIVE) {
    b[0] = ybcbeg;
    a[1 + 0 * 3] = 1.0;
    a[0 + 1 * 3] = 0.0;
  }
  else {
    std::cout << "\n";
    std::cout << "SPLINE_CUBIC_SET - Fatal error!\n";
    std::cout << "  IBCBEG must be 0, 1 or 2.\n";
    std::cout << "  The input value is " << ibcbeg << ".\n";
    delete[] a;
    delete[] b;
    exit(1);
  }
  //
  //  Set up the intermediate equations.
  //
  for (i = 1; i < n - 1; i++) {
    b[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]) - (y[i] - y[i - 1]) / (t[i] - t[i - 1]);
    a[2 + (i - 1) * 3] = (t[i] - t[i - 1]) / 6.0;
    a[1 + i * 3] = (t[i + 1] - t[i - 1]) / 3.0;
    a[0 + (i + 1) * 3] = (t[i + 1] - t[i]) / 6.0;
  }
  //
  //  Set up the last equation.
  //
  if (ibcend == SPLINE_USE_QUADRATIC) {
    b[n - 1] = 0.0;
    a[2 + (n - 2) * 3] = -1.0;
    a[1 + (n - 1) * 3] = 1.0;
  }
  else if (ibcend == SPLINE_USE_FIRST_DERIVATIVE) {
    b[n - 1] = ybcend - (y[n - 1] - y[n - 2]) / (t[n - 1] - t[n - 2]);
    a[2 + (n - 2) * 3] = (t[n - 1] - t[n - 2]) / 6.0;
    a[1 + (n - 1) * 3] = (t[n - 1] - t[n - 2]) / 3.0;
  }
  else if (ibcend == SPLINE_USE_SECOND_DERIVATIVE) {
    b[n - 1] = ybcend;
    a[2 + (n - 2) * 3] = 0.0;
    a[1 + (n - 1) * 3] = 1.0;
  }
  else {
    std::cout << "\n";
    std::cout << "SPLINE_CUBIC_SET - Fatal error!\n";
    std::cout << "  IBCEND must be 0, 1 or 2.\n";
    std::cout << "  The input value is " << ibcend << ".\n";
    delete[] a;
    delete[] b;
    exit(1);
  }
  //
  //  Solve the linear system.
  //
  if (n == 2 && ibcbeg == 0 && ibcend == 0) {
    ypp = new float[2];

    ypp[0] = 0.0;
    ypp[1] = 0.0;
  }
  else {
    ypp = d3_np_fs(n, a, b);

    if (!ypp) {
      std::cout << "\n";
      std::cout << "SPLINE_CUBIC_SET - Fatal error!\n";
      std::cout << "  The linear system could not be solved.\n";
      delete[] a;
      delete[] b;
      exit(1);
    }
  }

  delete[] a;
  delete[] b;
  return ypp;
}

/**
 * This function evaluates a piecewise cubic spline at a point.
 *
 * SPLINE_CUBIC_SET must have already been called to define the values of YPP.
 *
 * For any point T in the interval T(IVAL), T(IVAL+1), the form of
 * the spline is
 *
 *   SPL(T) = A
 *          + B * ( T - T(IVAL) )
 *          + C * ( T - T(IVAL) )**2
 *          + D * ( T - T(IVAL) )**3
 *
 * Here:
 *   A = Y(IVAL)
 *   B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
 *     - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
 *   C = YPP(IVAL) / 2
 *   D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
 *
 * Original Author: John Burkardt
 * http://www.csit.fsu.edu/~burkardt/cpp_src/spline/spline.html
 *
 * @param n The number of knots.
 * @param y The data values, size n, at the knots.
 * @param t The knot values, size n.
 * @param tval A point, typically between T[0] and T[N-1], at which the spline
 * is to be evalulated.  If TVAL lies outside this range,
 * extrapolation is used.
 * @param y The data values, size n, at the knots.
 * @param ypp The second derivatives, size n, of the spline at the knots.
 * @param ypval The derivative of the spline at TVAL.
 * @param yppval The second derivative of the spline at TVAL.
 * @return The value of the spline at TVAL.
 */
float SplineCubicValue(int n, float t[], float tval, float y[], float ypp[], float *ypval, float *yppval)
{
  float yval;
  int ival = n - 2;

  if (tval <= t[0]) {
    // enforce constant function outside of domain
    yval = y[0];
  }
  else if (tval >= t[n - 1]) {
    // enforce constant function outside of domain
    yval = y[n - 1];
    *ypval = 0.0;
    *ypval = 0.0;
  }
  else {
    // we're not outside of the range, so interpolate
    for (int i = 0; i < n - 1; i++) {
      if (tval < t[i + 1]) {
        ival = i;
        break;
      }
    }

    //  In the interval I, the polynomial is in terms of a normalized
    //  coordinate between 0 and 1.
    float dt = tval - t[ival];
    float h = t[ival + 1] - t[ival];

    yval = y[ival] + dt * ((y[ival + 1] - y[ival]) / h - (ypp[ival + 1] / 6.0 + ypp[ival] / 3.0) * h +
                           dt * (0.5 * ypp[ival] + dt * ((ypp[ival + 1] - ypp[ival]) / (6.0 * h))));

    *ypval = (y[ival + 1] - y[ival]) / h - (ypp[ival + 1] / 6.0 + ypp[ival] / 3.0) * h +
             dt * (ypp[ival] + dt * (0.5 * (ypp[ival + 1] - ypp[ival]) / h));

    *yppval = ypp[ival] + dt * (ypp[ival + 1] - ypp[ival]) / h;
  }

  return yval;
}

/**
 * D3_NP_FS factors and solves a D3 system.
 *
 * The D3 storage format is used for a tridiagonal matrix.
 * The superdiagonal is stored in entries (1,2:N), the diagonal in
 * entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
 * original matrix is "collapsed" vertically into the array.
 *
 * This algorithm requires that each diagonal entry be nonzero.
 * It does not use pivoting, and so can fail on systems that
 * are actually nonsingular.
 *
 * Example:
 *
 * Here is how a D3 matrix of order 5 would be stored:
 *
 *    *  A12 A23 A34 A45
 *   A11 A22 A33 A44 A55
 *   A21 A32 A43 A54  *
 *
 * Original Author: John Burkardt
 * http://www.csit.fsu.edu/~burkardt/cpp_src/spline/spline.html
 *
 * @param n The order of the linear system.
 *
 * @param a, 3*n size.  On input, the nonzero diagonals of the linear system.
 * On output, the data in these vectors has been overwritten by factorization
 * information.
 *
 * @param b The right hand side. size n.
 *
 * @return The solution of the linear system.
 * This is NULL if there was an error
 * because one of the diagonal entries was zero.  Size n.
 */
float *d3_np_fs(int n, float a[], float b[])
{
  int i;
  float *x;
  float xmult;
  //
  //  Check.
  //
  for (i = 0; i < n; i++) {
    if (a[1 + i * 3] == 0.0) {
      return NULL;
    }
  }

  x = new float[n];

  for (i = 0; i < n; i++) {
    x[i] = b[i];
  }

  for (i = 1; i < n; i++) {
    xmult = a[2 + (i - 1) * 3] / a[1 + (i - 1) * 3];
    a[1 + i * 3] = a[1 + i * 3] - xmult * a[0 + i * 3];
    x[i] = x[i] - xmult * x[i - 1];
  }

  x[n - 1] = x[n - 1] / a[1 + (n - 1) * 3];
  for (i = n - 2; 0 <= i; i--) {
    x[i] = (x[i] - a[0 + (i + 1) * 3] * x[i + 1]) / a[1 + i * 3];
  }

  return x;
}

#define NR_END 1
#define FREE_ARG char *

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v = (float *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
  if (!v) ErrorExit(ERROR_NOMEMORY, "could not allocate vector(%ld, %ld)", nl, nh);
  return v - nl + NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float *)));
  if (!m) ErrorExit(ERROR_NOMEMORY, "could not allocate matrix(%ld, %ld)", nrow, ncol);
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (float *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
  if (!m[nrl]) ErrorExit(ERROR_NOMEMORY, "could not allocate matrix(%ld, %ld) array", nrow, ncol);
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */ { free((FREE_ARG)(v + nl - NR_END)); }

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

/*  AUTHOR : Martin Reuter
    CREATED: 12/28/2009
    Computes Denman and Beavers square root iteration up to an epsilon in imax steps

    Input: 4x4 affine transformation matrix M
    Algorithm computes sqrt of 3x3 part and then adjusts the translation part
    Output: 4x4 sqrt of M
*/
vnl_matrix_fixed< double, 4, 4 > MatrixSqrt(const vnl_matrix_fixed< double, 4, 4 > &m)
{
  // assert(m.rows() == 4 && m.cols() == 4);

  vnl_matrix_fixed< double, 3, 3 > R;  // = m.extract(3,3,0,0);
  for (int rr = 0; rr < 3; rr++)
    for (int cc = 0; cc < 3; cc++) {
      R[rr][cc] = m[rr][cc];
    }

  // Denman and Beavers square root iteration

  int imax = 100;
  double eps = 0.00001;  // important to be small to guarantee symmetry,
                         // but even adding two zeros did not show differences in tests
  double err = 1000;
  // cout << "using square root iteartion (" << imax << ")"<< endl;
  vnl_matrix_fixed< double, 3, 3 > Yn(R);
  vnl_matrix_fixed< double, 3, 3 > Zn;
  Zn.set_identity();
  vnl_matrix_fixed< double, 3, 3 > Yni, Zni;
  vnl_matrix_fixed< double, 3, 3 > Ysq;

  int count = 0;
  while (count < imax && err > eps) {
    count++;

    // store invrse here (we change Yn below)
    Yni = vnl_inverse(Yn);
    Zni = vnl_inverse(Zn);

    // add inverse:
    Yn += Zni;
    Zn += Yni;

    Yn *= 0.5;
    Zn *= 0.5;

    Ysq = Yn * Yn;
    Ysq -= R;
    err = Ysq.absolute_value_max();
    // cout << " iteration " << count << "  err: "<< err << endl;
  }

  if (count > imax) {
    std::cerr << "Matrix Sqrt did not converge in " << imax << " steps!" << std::endl;
    std::cerr << "   ERROR: " << err << std::endl;
    // std::assert(err <= eps);
    exit(1);
  }

  // compute new T
  // Rh1 = R + I
  vnl_matrix_fixed< double, 3, 3 > Rh1(Yn);
  Rh1[0][0] += 1;
  Rh1[1][1] += 1;
  Rh1[2][2] += 1;

  vnl_vector_fixed< double, 3 > T;
  T[0] = m[0][3];
  T[1] = m[1][3];
  T[2] = m[2][3];

  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  vnl_vector_fixed< double, 3 > Th = vnl_inverse(Rh1) * T;  // vnl_svd<double>(Rh1).solve(T);

  // put everything together:
  vnl_matrix_fixed< double, 4, 4 > msqrt;
  msqrt[0][3] = Th[0];
  msqrt[1][3] = Th[1];
  msqrt[2][3] = Th[2];
  msqrt[3][0] = 0.0;
  msqrt[3][1] = 0.0;
  msqrt[3][2] = 0.0;
  msqrt[3][3] = 1.0;
  for (int c = 0; c < 3; c++)
    for (int r = 0; r < 3; r++) msqrt[r][c] = Yn[r][c];

  //    bool test = true;
  //    if (test)
  //    {
  //       vnl_matrix < double > ms2 = msqrt * msqrt;
  //       ms2 -= m;
  //       double sum = ms2.absolute_value_max();
  //       if (sum > eps)
  //       {
  //          cerr << " Error : " << sum << endl;
  // 				 cerr << " sqrt(M): " << endl << msqrt << endl;
  //          cerr << endl;
  //          assert(1==2);
  //       }
  //    }

  return msqrt;
}

MATRIX *MatrixSqrt(MATRIX *m, MATRIX *sqrtm)
{
  int i, j;
  vnl_matrix_fixed< double, 4, 4 > vnl_m;

  if (m->rows != 4 || m->cols != 4) {
    std::cerr << " Numerics MatrixSqrt m must be 4x4 " << std::endl;
    exit(1);
  }

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) vnl_m[i][j] = m->rptr[i + 1][j + 1];

  vnl_matrix_fixed< double, 4, 4 > vnl_msqrt = MatrixSqrt(vnl_m);

  if (sqrtm == NULL) sqrtm = MatrixAlloc(4, 4, MATRIX_REAL);

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) sqrtm->rptr[i + 1][j + 1] = vnl_msqrt[i][j];

  return sqrtm;
}

/*
   --------------------------------------------------------------------------

   Author: Silvester Czanner

   Created 7/26/2006

   Description: GSL library replacement for Freesurfer

   --------------------------------------------------------------------------
*/

//##################### FUNC PROTOS ######################################

static const unsigned long int SC_MASK_LO = 0x00ffffffUL; /*2^24 - 1 */
static const unsigned long int SC_MASK_HI = ~0x00ffffffUL;
static const unsigned long int two_over_24 = 16777216; /*2^24 */

unsigned long int sc_inc_status(sc_status_t *status);

double sc_uni_pos(const sc_rng *r);
double sc_uni(const sc_rng *r);

double sc_gamma_frac(const sc_rng *r, const double a);
double sc_gamma_large(const sc_rng *r, const double a);
double sc_gamma_int(const sc_rng *r, const unsigned int a);
double sc_ugauss(const sc_rng *r);

//##################### FUNCTIONS #######################################

//########### VNL ###################################################
/*!
  \fn int sc_linalg_cholesky_decomp(MATRIX *U)
  \params MATRIX *U - input and output matrix
  Computes cholesky decomposition C of U in-place.
  U = C'*C. C is upper triangular (same as matlab)
  Note: when using as a temporal filter, want lower tri.
  Returns 0 on success and nonzero otherwise (eg,
  when U is non-positive def).
 */
int sc_linalg_cholesky_decomp(MATRIX *U)
{
  int i, j, err;

  vnl_matrix< double > P(U->rows, U->cols);
  vnl_matrix< double > C(U->rows, U->cols);

  for (i = 0; i < U->rows; i++)
    for (j = 0; j < U->cols; j++) P(i, j) = *MATRIX_RELT(U, i + 1, j + 1);

  vnl_cholesky chol(P);
  err = chol.rank_deficiency();
  if (err) return (err);

  C = chol.upper_triangle();

  for (i = 0; i < U->rows; i++)
    for (j = 0; j < U->cols; j++) *MATRIX_RELT(U, i + 1, j + 1) = C(i, j);

  return (0);
}

//######### UTILS #######################################################

static void sc_err_msg(const char *msg)
{
  printf("ERROR: %s", msg);
  exit(1);
}

//########################## RNG ###############################

unsigned long int sc_inc_status(sc_status_t *status)
{
  int i, j;
  long int step;

  i = status->i;
  j = status->j;
  step = status->u[j] - status->u[i] - status->carry;

  if (step & SC_MASK_HI) {
    status->carry = 1;
    step &= SC_MASK_LO;
  }
  else {
    status->carry = 0;
  }
  status->u[i] = step;

  if (i == 0) {
    i = 23;
  }
  else {
    i--;
  }
  status->i = i;

  if (j == 0) {
    j = 23;
  }
  else {
    j--;
  }
  status->j = j;

  return (step);
}

unsigned long int sc_rng_get(const sc_rng *r)
{
  sc_status_t *status;
  int i, skip;
  unsigned long int ret;

  status = (sc_status_t *)r->status;
  skip = status->skip;
  ret = sc_inc_status(status);

  status->n++;
  if (status->n == 24) {
    status->n = 0;
    for (i = 0; i < skip; i++) sc_inc_status(status);
  }

  return (ret);
}

void sc_rng_set(sc_rng *r, unsigned long int seed_in)
{
  int const_lux;
  sc_status_t *status;
  int i;
  unsigned long int tmp;
  long int seed;

  const_lux = 389;
  status = (sc_status_t *)r->status;

  if (seed_in == 0) {
    seed = 314159265;
  }
  else {
    seed = seed_in;
  }

  for (i = 0; i < 24; i++) {
    tmp = seed / 53668;
    seed = 40014 * (seed - tmp * 53668) - tmp * 12211;
    if (seed < 0) {
      seed = seed + 2147483563;
    }
    status->u[i] = seed % two_over_24;
  }

  status->i = 23;
  status->j = 9;
  status->n = 0;
  status->skip = const_lux - 24;

  if (status->u[23] & SC_MASK_HI) {
    status->carry = 1;
  }
  else {
    status->carry = 0;
  }

  return;
}

sc_rng *sc_rng_alloc(const sc_rng_type *T)
{
  sc_rng *r;

  r = (sc_rng *)malloc(sizeof(sc_rng));
  if (r == NULL) {
    sc_err_msg("*sc_rng_alloc(): problem with allocation");
  }

  r->status = malloc(T->size);
  r->type = T;

  sc_rng_set(r, 0);

  return (r);
}

void sc_rng_free(sc_rng *r)
{
  if (r->status != NULL) free(r->status);

  if (r != NULL) free(r);

  return;
}

double sc_uni(const sc_rng *r)
{
  double u;

  u = (double)(sc_rng_get(r) / 16777216.0);

  return (u);
}

double sc_uni_pos(const sc_rng *r)
{
  double x;
  do {
    x = (double)(sc_rng_get(r) / 16777216.0);
  } while (x == 0);

  return (x);
}

//##################### RAN ###

double sc_ran_flat(const sc_rng *r, const double a, const double b)
{
  double u, ret;

  u = sc_uni(r);
  ret = a * (1 - u) + b * u;
  return (ret);
}

double sc_ran_gaussian(const sc_rng *r, const double sigma)
{
  double x, y, r2, ret;

  do {
    x = -1 + 2 * sc_uni(r);
    y = -1 + 2 * sc_uni(r);
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0);

  ret = sigma * y * sqrt(-2.0 * log(r2) / r2);
  return (ret);
}

double sc_ugauss(const sc_rng *r)
{
  double ret;

  ret = sc_ran_gaussian(r, 1.0);

  return (ret);
}

double sc_gamma_int(const sc_rng *r, const unsigned int a)
{
  unsigned int i;
  double prod, ret;

  if (a < 12) {
    prod = 1;
    for (i = 0; i < a; i++) {
      prod *= sc_uni_pos(r);
    }
    ret = -log(prod);
    return (ret);
  }
  else {
    ret = sc_gamma_large(r, (double)a);
    return (ret);
  }
}

double sc_gamma_large(const sc_rng *r, const double a)
{
  double sqa, y, v, ret;

  sqa = sqrt(2 * a - 1);
  do {
    do {
      y = tan(M_PI * sc_uni(r));
      ret = sqa * y + a - 1;
    } while (ret <= 0);
    v = sc_uni(r);
  } while (v > (1 + y * y) * exp((a - 1) * log(ret / (a - 1)) - sqa * y));

  return (ret);
}

double sc_gamma_frac(const sc_rng *r, const double a)
{
  double p, q, u, v, ret;
  p = M_E / (a + M_E);
  do {
    u = sc_uni(r);
    v = sc_uni_pos(r);

    if (u < p) {
      ret = exp((1 / a) * log(v));
      q = exp(-ret);
    }
    else {
      ret = 1 - log(v);
      q = exp((a - 1) * log(ret));
    }
  } while (sc_uni(r) >= q);

  return (ret);
}

double sc_ran_gamma(const sc_rng *r, const double a, const double b)
{
  unsigned int na;
  double ret;

  na = (unsigned int)floor(a);
  if (a == na) {
    ret = b * sc_gamma_int(r, na);
    return (ret);
  }
  else if (na == 0) {
    ret = b * sc_gamma_frac(r, a);
    return (ret);
  }
  else {
    ret = b * (sc_gamma_int(r, na) + sc_gamma_frac(r, a - na));
    return (ret);
  }
}

double sc_ran_fdist(const sc_rng *r, const double nu1, const double nu2)
{
  double Y1, Y2, ret;

  Y1 = sc_ran_gamma(r, nu1 / 2, 2.0);
  Y2 = sc_ran_gamma(r, nu2 / 2, 2.0);

  ret = (Y1 * nu2) / (Y2 * nu1);
  return (ret);
}

double sc_ran_chisq(const sc_rng *r, const double nu)
{
  double ret;

  ret = 2 * sc_ran_gamma(r, nu / 2, 1.0);
  return (ret);
}

double sc_ran_tdist(const sc_rng *r, const double nu)
{
  double Y1, Y2, Z, ret;

  if (nu <= 2) {
    Y1 = sc_ugauss(r);
    Y2 = sc_ran_chisq(r, nu);
    ret = Y1 / sqrt(Y2 / nu);
    return (ret);
  }
  else {
    do {
      Y1 = sc_ugauss(r);
      Y2 = sc_ran_exponential(r, 1 / (nu / 2 - 1));
      Z = Y1 * Y1 / (nu - 2);
    } while (1 - Z < 0 || exp(-Y2 - Z) > (1 - Z));

    ret = Y1 / sqrt((1 - 2 / nu) * (1 - Z));
    return (ret);
  }
}

double sc_ran_exponential(const sc_rng *r, const double mu)
{
  double u, ret;

  u = sc_uni_pos(r);
  ret = -mu * log(u);

  return (ret);
}

double sc_ran_binomial_pdf(unsigned int k, double p, unsigned int n)
{
  double ret;

  if (k == 0) {
    ret = bdtr(k, n, p);
  }
  else {
    ret = bdtr(k, n, p) - bdtr(k - 1, n, p);
  }

  return (ret);
}

//##################### CDF ###

double sc_cdf_flat_Q(double x, double a, double b)
{
  double ret = 0.0;

  if (x <= a) ret = 0.0;
  if (a < x && x < b) ret = (x - a) / (b - a);
  if (b <= x) ret = 1.0;

  return (1.0 - ret);
}

double sc_cdf_flat_Qinv(double Q, double a, double b)
{
  double ret = 0.0;

  if (Q == 0.0) ret = b;

  if (Q == 1.0) ret = a;

  if (Q > 0.0 && Q < 1.0) ret = Q * a + (1.0 - Q) * b;

  return (ret);
}

double sc_cdf_fdist_Q(double x, double nu1, double nu2)
{
  double ret = 0.0;

  ret = fdtrc(nu1, nu2, x);

  return (ret);
}

double sc_cdf_fdist_Qinv(double Q, double nu1, double nu2)
{
  double ret = 0.0;

  ret = fdtri(nu1, nu2, Q);

  return (ret);
}

double sc_cdf_tdist_Q(double x, double nu)
{
  double ret = 0.0;

  ret = 1.0 - stdtr(nu, x);

  return (ret);
}

double sc_cdf_tdist_Qinv(double Q, double nu)
{
  double ret = 0.0;

  ret = (-1) * stdtri(nu, Q);

  return (ret);
}

double sc_cdf_gaussian_Q(double x, double nu)
{
  double ret = 0.0;

  if (nu != 1.0) {
    printf("ERROR sc_cdf_gaussian_Q: sigma=%f != 1.0\n\n", nu);
    exit(1);
  }

  ret = 1 - ndtr(x);

  return (ret);
}

double sc_cdf_gaussian_Qinv(double Q, double nu)
{
  double ret = 0.0;

  if (nu != 1.0) {
    printf("ERROR sc_cdf_gaussian_Qinv: sigma=%f != 1.0\n\n", nu);
    exit(1);
  }

  ret = (-1) * ndtri(Q);

  return (ret);
}

double sc_cdf_chisq_Q(double x, double nu)
{
  double ret = 0.0;

  ret = chdtrc(nu, x);

  return (ret);
}

double sc_cdf_chisq_Qinv(double Q, double nu)
{
  double ret = 0.0;

  ret = chdtri(nu, Q);

  return (ret);
}
