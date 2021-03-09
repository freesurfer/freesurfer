/**
 * @brief A class to solve overconstrained system A X = B
 *
 *   it uses either least squares (standard regression)
 *   or a robust estimator (Tukey's Biweight with iterative reweighted least
 *   squares)
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

#include "Regression.h"
#include <iostream>
#include <cassert>
#include <math.h>
#include <limits>
#include <vector>
#include <fstream>
#include "RobustGaussian.h"

#define export // obsolete feature 'export template' used in these headers 
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_qr.h>
#undef export

#include <stdio.h>
#include <stdlib.h>
#include "error.h"

using namespace std;

template<class T>
vnl_vector<T> Regression<T>::getRobustEst(double sat, double sig)
{
  vnl_vector<T> w;
  return getRobustEstW(w, sat, sig);
}

// void print_mem_usage ()
// {
//  char buf[30];
//         snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
//         FILE* pf = fopen(buf, "r");
//         if (pf) {
//             unsigned size; //       total program size
//             //unsigned resident;//   resident set size
//             //unsigned share;//      shared pages
//             //unsigned text;//       text (code)
//             //unsigned lib;//        library
//             //unsigned data;//       data/stack
//             //unsigned dt;//         dirty pages (unused in Linux 2.6)
//             //fscanf(pf, "%u" /* %u %u %u %u %u"*/, &size/*, &resident, &share, &text, &lib, &data*/);
//             fscanf(pf, "%u" , &size);
//             //DOMSGCAT(MSTATS, std::setprecision(4) << size / (1024.0) << "MB mem used");
//      cout <<  size / (1024.0) << " MB mem used" << endl;
//         }
//         fclose(pf);
// //while (1)
// //{
// //   if ('n' == getchar())
// //       break;
// //}
// }

template<class T>
vnl_vector<T> Regression<T>::getRobustEstW(vnl_vector<T>& w, double sat,
    double sig)
{
  if (A)
    return getRobustEstWAB(w, sat, sig);
  else
    return vnl_vector<T>(1, getRobustEstWB(w, sat, sig));
}

template<class T>
double Regression<T>::getRobustEstWB(vnl_vector<T>& w, double sat, double sig)
{
  // constants
  int MAXIT = 20;
  //double EPS = 2e-12;
  double EPS = 2e-6;

  T muold = 0;
  T mu = 0;

  assert(A == NULL);

  // initialize with median:
//   {
//     double bb[B->rows];
//     for (int i = 1 ; i<=b->rows; i++)
//       bb[i-1] = B->rptr[i][1];
//     mu = RobustGaussian::median(bb,B->rows);
//   } // freeing bb

  mu = RobustGaussian<T>::median(b->data_block(), b->size());

  muold = mu + 100; // set out of reach of eps
  int count = 0;
  vnl_vector<T> r;
  w.fill(1.0);
  vnl_vector<T> rdsigma;
//   MATRIX * w = MatrixConstVal(1.0,B->rows,1,NULL);
//   MATRIX * r = MatrixConstVal(1.0,B->rows,1,NULL);
//   MATRIX * rdsigma = MatrixConstVal(1.0,B->rows,1,NULL);
  T sigma;
  double d1, d2, wi;

  //cout << endl<<endl<<" Values: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",B);

  while (fabs(muold - mu) > EPS && count < MAXIT)
  {
    count++;
    muold = mu;
    //cout << endl << "count " << count << " Myold " << muold << endl;
//     for (int i = 1 ; i<=B->rows; i++)
//       r->rptr[i][1] = muold - B->rptr[i][1];
    r = muold - *b;

    //cout << " residuals: " << endl;
    //   MatrixPrintFmt(stdout,"% 2.8f",r);

    sigma = getSigmaMAD(r);
    //cout << " sigma: " << sigma << endl;
    if (sigma < EPS) // if all r are the same (usually zero)
    {
      mu = muold;
//       for (int i=1 ; i<=w->rows; i++)
//         w->rptr[i][1] = 1;
      w.fill(1);
      break;
    }

//    for (int i = 1 ; i<=B->rows; i++)
//      rdsigma->rptr[i][1] = r->rptr[i][1] / sigma;
    rdsigma = (T) (1.0 / sigma) * r;
    //cout << " r/sigma: " << endl;
    //   MatrixPrintFmt(stdout,"% 2.8f",rdsigma);

    getSqrtTukeyDiaWeights(rdsigma, w, sat); // here we get sqrt of weights into w
    //cout << " weights: " << endl;
    //   MatrixPrintFmt(stdout,"% 2.8f",w);

    // compute new parameter mu (using weights)
    d1 = 0;
    d2 = 0;
    // d3 = 0;
    for (unsigned int i = 0; i < b->size(); i++)
    {
      wi = w[i] * w[i];
      // d1 += wi * r[i];  not sure if this is right, think it should be b (see next line)
      d1 += wi * b->operator[](i);
      d2 += wi;
      // d3 += wi * r[i] * r[i];
    }

    if (d2 < EPS) // weights zero
    {

      assert(d2 >= EPS);
    }
    mu = muold - (T) d1 / d2;

  }
//cout << "!!! final mu :  " << mu << endl;
//  MatrixFree(&r);
//  MatrixFree(&rdsigma);

  return mu;

}

/** Solves overconstrained system A p = b using
 M estimators (Beaton and Tukey's biweigt function).
 Returns vectors p and w (w by reference as it is large)
 where p is the robust solution and w the used weights.
 Method: iterative reweighted least squares
 */
template<class T>
vnl_vector<T> Regression<T>::getRobustEstWAB(vnl_vector<T>& wfinal, double sat,
    double sig)
{
  if (verbose > 1)
  {
    cout << "  Regression<T>::getRobustEstWAB( "<<sat<<" , "<<sig<<" ) " ;
    if (floatsvd) cout << "  FLOAT version " ;
    else cout << "  DOUBLE version " ;
    cout << endl;
  }
  
  // constants
  int MAXIT = 20;
  //double EPS = 2e-16;
  double EPS = 2e-12;

  // variables
  std::vector<T> err(MAXIT + 1);
  err[0] = numeric_limits<T>::infinity();
  err[1] = 1e20;
  double sigma;

  int arows = A->rows(); // large (voxels)
  int acols = A->cols(); // small (parameters)

  //pre-alocate vectors
  // init residuals (based on zero p, so r := b )
  vnl_vector<T> * r = new vnl_vector<T>(*b);
  vnl_vector<T> * p = new vnl_vector<T>(acols);
  vnl_vector<T> * w = new vnl_vector<T>(arows);

//  if (! w.valid() || ! r.valid() || !p.valild())
//     ErrorExit(ERROR_NO_MEMORY,"Regression<T>::getRobustEstWAB could not allocate memory for w,r,p") ;

  vnl_vector<T> *lastp = new vnl_vector<T>(acols);
  vnl_vector<T> *lastw = new vnl_vector<T>(arows);
//  if (! lastw.valid() || !lastr.valid() || !lastp.valid())
//     ErrorExit(ERROR_NO_MEMORY,"Regression<T>::getRobustEstWAB could not allocate memory for lastw,lastr,lastp") ;
  vnl_vector<T> *vtmp = NULL;

  int count = 0;
  int incr = 0;
  // iteration until we increase the error, we reach maxit or we have no error
  do
  {
    count++; //first = 1

    if (count > 1)
    {
      // store lastp (not necessary in first run)
      // instead of copy, we swap pointers:
      vtmp = lastp;
      lastp = p;
      p = vtmp;

      // store last weights
      // instead of copy, we swap pointers:
      vtmp = lastw;
      lastw = w;
      w = vtmp;
    }
    
    // normalize r and compute weights (or rather w = sqrt of weights)
    sigma = getSigmaMAD(*r);
    if (sigma < EPS) // e.g. if images are identical
    {
      cout << "  Sigma too small: " << sigma << " (identical images?)" << endl;
      w->fill(1.0);
    }
    else
    {
      //cout << "Sigma: " << sigma << endl;
      *r *= (1.0 / sigma);
      // here we get sqrt of weights into w
      getSqrtTukeyDiaWeights(*r, *w, sat);
    }
    // free residuals (to reduce max memory load) seems to have no effect on speed
    r->clear();

    // compute weighted least squares
    if (floatsvd)
      *p = getWeightedLSEstFloat(*w);
    else
      *p = getWeightedLSEst(*w);

    // compute new residuals
    *r = *b - (*A * *p);

    // and total errors (using new r)
    // err = sum (w r^2) / sum (w)
    T swr = 0;
    T sw = 0;
    for (unsigned int rr = 0; rr < r->size(); rr++)
    {
      T t1 = w->operator[](rr);
      T t2 = r->operator[](rr);
      t1 *= t1; // remember w is the sqrt of the weights
      t2 *= t2;
      sw += t1;
      swr += t1 * t2;
    }
    err[count] = swr / sw;
    //cout << "err [ " << count << " ] = " << err[count] << endl;
    if (err[count - 1] <= err[count])
      incr++;
  } while (incr < 1 && count < MAXIT && err[count] > EPS);

  delete (r); // won't be needed below

  vnl_vector<T> pfinal;
  if (err[count] > err[count - 1])
  {
    // take previous values (since actual values made the error to increase)
    // cout << " last step was no improvement, taking values : "<<  count-1 << endl;
    pfinal = *lastp;
    wfinal = *lastw;
    if (verbose > 1)
      cout << "     Step: " << count - 2 << " ERR: " << err[count - 1] << endl;
    lasterror = err[count - 1];
  }
  else
  {
    pfinal = *p;
    wfinal = *w;
    if (verbose > 1)
      cout << "     Step: " << count - 1 << " ERR: " << err[count] << endl;
    lasterror = err[count];
  }
  delete (p);
  delete (w);
  delete (lastw);
  delete (lastp);
      
    
  // compute statistics on weights:
  double d = 0.0;
  double dd = 0.0;
  double ddcount = 0;
  int zcount = 0;
  T val;
  for (unsigned int i = 0; i < wfinal.size(); i++)
  {
    val = wfinal[i];
    d += val;
    if (fabs(b->operator[](i)) > 0.00001)
    {
      dd += val;
      ddcount++;
      if (val < 0.1)
        zcount++;
    }
  }
  d /= wfinal.size();
  dd /= ddcount;
  if (verbose > 1)
    cout << "          weights average: " << dd << "  zero: "
        << (double) zcount / ddcount << flush;
  //"  on significant b vals ( " << ddcount << " ): " << dd <<endl;
  lastweight = dd;
  lastzero = (double) zcount / ddcount;

  return pfinal;
}


/** Solving \f$ p = [A^T W A]^{-1} A^T W b\f$     (with \f$ W = diag(w_i^2) \f$ )
 done by computing \f$ M := \sqrt{W} A\f$ and  \f$ v := \sqrt{W} b\f$
 then we have \f$ p = [ M^T M ]^{-1} M^T v  \f$
 or \f$ M p = v \f$, this we solve with QR decomposition (faster than svd).
 \param w vector representing a diagnoal matrix with the sqrt of the weights as elements
 */
template<class T>
vnl_vector<T> Regression<T>::getWeightedLSEst(const vnl_vector<T> & w)
{
  unsigned int rr, cc;

  assert(w.size() == A->rows());

  // compute wA  where w = diag(sqrt(W));
  vnl_matrix<T> wA(A->rows(), A->cols());
  for (rr = 0; rr < A->rows(); rr++)
    for (cc = 0; cc < A->cols(); cc++)
      wA(rr, cc) = A->operator()(rr, cc) * w(rr);

  vnl_qr<T>* QR = new vnl_qr<T>(wA);
  // I could maybe delete wA here?  

  // compute wb:
  vnl_vector<T> wb(b->size());
  //if (!wb.valid()) ErrorExit(ERROR_NO_MEMORY,"Regression<T>::getWeightedLSEst could not allocate memory for wb") ;
  for (rr = 0; rr < b->size(); rr++)
    wb[rr] = b->operator[](rr) * w[rr];

  // solve wA p = wb
  vnl_vector<T> p = QR->solve(wb);

  delete (QR);

  return p;
}

/** Uses FLOAT internaly.
 Solving \f$ p = [A^T W A]^{-1} A^T W b\f$     (with \f$ W = diag(w_i^2) \f$ )
 done by computing \f$ M := \sqrt{W} A\f$ and  \f$ v := \sqrt{W} b\f$
 then we have \f$ p = [ M^T M ]^{-1} M^T v  \f$
 or \f$ M p = v \f$, this we solve with QR decomposition (faster than svd).
 \param w vector representing a diagnoal matrix with the sqrt of the weights as elements
 */
template<class T>
vnl_vector<T> Regression<T>::getWeightedLSEstFloat(const vnl_vector<T> & w)
{
  unsigned int rr, cc;

  assert(w.size() == A->rows());

  // compute wA  where w = diag(sqrt(W));
  vnl_matrix<float> wA(A->rows(), A->cols());
  for (rr = 0; rr < A->rows(); rr++)
    for (cc = 0; cc < A->cols(); cc++)
      wA(rr, cc) = (float) (A->operator()(rr, cc) * w(rr));

  vnl_qr<float>* QR = new vnl_qr<float>(wA);
  // I could maybe delete wA here?  

  // compute wb:
  vnl_vector<float> wb(b->size());
  //if (!wb.valid()) ErrorExit(ERROR_NO_MEMORY,"Regression<T>::getWeightedLSEst could not allocate memory for wb") ;
  for (rr = 0; rr < b->size(); rr++)
    wb[rr] = b->operator[](rr) * w[rr];

  // solve
  vnl_vector<float> p = QR->solve(wb);
  delete (QR);

  // copy p back:
  vnl_vector<T> pd(p.size());
  for (rr = 0; rr < p.size(); rr++)
    pd[rr] = p[rr];

  return pd;
}

// template <class T>
// vnl_vector< T >  Regression<T>::getWeightedLSEst(const vnl_vector< T > & w)
// // w is a vector representing a diagnoal matrix with the sqrt of the weights as elements
// // solving p = [A^T W A]^{-1} A^T W b     (with W = diag(w_i^2) )
// // done by computing M := Sqrt(W) A and and  wb := sqrt(W) b
// // then we have p = [ M^T M ]^{-1} M^T wb    which is solved by computing the pseudo inverse
// // of psdi(M) and then p = psdi(M) * wb
// {
//   unsigned int rr, cc;
//   
//   assert(w.size() == A->rows());
// 
//   // compute wA  where w = diag(sqrt(W));
//   vnl_matrix < T > wA(A->rows(), A->cols());  
//   for (rr = 0; rr < A->rows();rr++)
//     for (cc = 0; cc < A->cols();cc++)
//       wA(rr,cc) = A->operator()(rr,cc) * w(rr);
//   
//   // compute pseudoInverse of wA:
//   // uses a LOT of memory!!! for a 520MB matrix A it will be close to 2Gig
//   vnl_svd< T >* svdMatrix = new vnl_svd< T >( wA );
//   if (! svdMatrix->valid() )
//   {
//     cerr << "    Regression<T>::getWeightedLSEst    could not compute pseudo inverse!" << endl;
//     exit(1);
//   }
//   wA = svdMatrix->pinverse(); 
//   delete (svdMatrix);
// 
// 
//   // compute wb:
//   vnl_vector < T > wb(b->size());
//   //if (!wb.valid()) ErrorExit(ERROR_NO_MEMORY,"Regression<T>::getWeightedLSEst could not allocate memory for wb") ;
//   for (rr = 0;rr<b->size();rr++)
//     wb[rr] = b->operator[](rr) * w[rr];
//   
//   // compute wAi * wb
//   vnl_vector < T > p = wA * wb;
// 
//   return p;
// }

// template <class T>
// vnl_vector< T >  Regression<T>::getWeightedLSEstFloat(const vnl_vector< T > & w)
// // uses FLOAT internaly
// // w is a vector representing a diagnoal matrix with the sqrt of the weights as elements
// // solving p = [A^T W A]^{-1} A^T W b     (with W = diag(w_i^2) )
// // done by computing M := Sqrt(W) A and and  wb := sqrt(W) b
// // then we have p = [ M^T M ]^{-1} M^T wb    which is solved by computing the pseudo inverse
// // of psdi(M) and then p = psdi(M) * wb
// {
//   unsigned int rr, cc;
//   
//   assert(w.size() == A->rows());
// 
//   // compute wA  where w = diag(sqrt(W));
//   vnl_matrix < float > wA(A->rows(), A->cols());
//   for (rr = 0; rr < A->rows();rr++)
//     for (cc = 0; cc < A->cols();cc++)
//       wA(rr,cc) = (float)(A->operator()(rr,cc) * w(rr));
//   
//   // compute pseudoInverse of wA:
//   // uses a LOT of memory!!! for a 520MB double matrix A it will be close to 2Gig 
//   // therefore using float should reduce the max memory usage by 1 gig
//   vnl_svd< float >* svdMatrix = new vnl_svd< float >( wA );
//   if (! svdMatrix->valid() )
//   {
//     cerr << "    Regression<T>::getWeightedLSEst    could not compute pseudo inverse!" << endl;
//     exit(1);
//   }
//   wA = svdMatrix->pinverse(); 
//   delete (svdMatrix);
// 
// 
//   // compute wb:
//   vnl_vector <float > wb(b->size());
//   //if (!wb.valid()) ErrorExit(ERROR_NO_MEMORY,"Regression<T>::getWeightedLSEst could not allocate memory for wb") ;
//   for (rr = 0;rr<b->size();rr++)
//     wb[rr] = b->operator[](rr) * w[rr];
//   
//   // compute wAi * wb
//   vnl_vector < float > p = wA * wb;
// 
//   // copy p back:
//    vnl_vector < T > pd(p.size());
//    for (rr = 0;rr<p.size();rr++)
//      pd[rr] = p[rr];
// 
//   return pd;
// }
// 

template<class T>
vnl_vector<T> Regression<T>::getLSEst()
{
  //cout << " Regression<T>::getLSEst " << endl;
  lastweight = -1;
  lastzero = -1;
  if (A == NULL) // LS solution is just the mean of B
  {
    assert(b!=NULL);
    T d = 0;
    for (unsigned int r = 0; r < b->size(); r++)
      d += b->operator[](r);
    vnl_vector<T> p(1);
    p(1) = d / b->size();

    return p;
  }

//    vnl_matrix< float > vnlX( ioA->data, ioA->rows, ioA->cols );
  vnl_svd<T> svdMatrix(*A);
  if (!svdMatrix.valid())
  {
    cerr << "    Regression<T>::getLSEst   could not compute pseudo inverse!"
        << endl;
    exit(1);
  }
  vnl_matrix<T>* Ai = new vnl_matrix<T>(svdMatrix.pinverse());
  vnl_vector<T> p = *Ai * *b;
  delete (Ai);

  // compute error:
  vnl_vector<T> R = (*A * p) - *b;
  double serror = 0;
  unsigned int rr;
  unsigned int n = R.size();
  for (rr = 0; rr < n; rr++)
  {
    serror += R[rr] * R[rr];
  }
//  cout << "     squared error: " << serror << endl;
  lasterror = serror;

  return p;
}

template<class T>
void getTukeyBiweight(const vnl_vector<T>& r, vnl_vector<T>& w, double sat)
{
  unsigned int n = r.size();
  assert(n==w.size());

  double aa, bb;
  for (unsigned int i = 0; i < n; i++)
  {
    if (fabs(r(i)) > sat)
      w(i) = (T) (sat * sat / 2.0);
    else
    {

      aa = r[i] / sat;
      bb = 1.0 - aa * aa;
      w[i] = (T) ((sat * sat / 2.0) * (1.0 - bb * bb * bb));
    }
  }
  //return w;
}

/** Computes sum of partial derivatives d_rho/d_sat(r_i)
 */
template<class T>
double Regression<T>::getTukeyPartialSat(const vnl_vector<T>& r, double sat)
{

  double sum = 0;
  double rt1;
  double rt2;
  double rt4;
  double sat3 = sat * sat * sat;
  double sat5 = sat3 * sat * sat;
  unsigned int rr;
  for (rr = 0; rr < r.size(); rr++)
  {
    rt1 = r[rr];
    if (fabs(rt1) > sat)
    {
      sum += sat;
    }
    else
    {
      rt2 = rt1 * rt1;
      rt4 = rt2 * rt2;
      sum += 3.0 * rt4 / sat3 - 2.0 * rt4 * rt2 / sat5;
    }
  }

  return sum;
}

/** Computes sqrt(weights) for a reweighted least squares approach.
 Returns elements of diag matrix W (as a column vector).
 Input: r (residuals) and sat
 Output: w sqrt of weights
 */
template<class T>
void Regression<T>::getSqrtTukeyDiaWeights(const vnl_vector<T>& r,
    vnl_vector<T>& w, double sat)
{
  //cout << " getTukeyDiaWeights  r size: " << r->rows << " , " << r->cols << endl;

  unsigned int n = r.size();
  assert(n == w.size());

  double t1;
  //double t2;
  unsigned int rr;
  //int ocount = 0;
  for (rr = 0; rr < n; rr++)
  {
    // cout << " fabs: " << fabs(r->rptr[rr][cc]) << " sat: " << sat << endl;
    if (fabs(r[rr]) >= sat)
    {
      w(rr) = 0.0;
      //ocount++;
    }
    else
    {
      t1 = r[rr] / sat;
      //t2 = 1.0 - t1 * t1;
      w(rr) = (T) (1.0 - t1 * t1); // returning sqrt
    }
  }
  //cout << " over threshold: " << ocount << " times ! " << endl;
  //return w;
}

template<class T>
T Regression<T>::VectorMedian(const vnl_vector<T>& v)
{
  unsigned int n = v.size();
  T* t = (T *) calloc(n, sizeof(T));
  if (t == NULL)
    ErrorExit(ERROR_NO_MEMORY,
        "Regression<T>::VectorMedian could not allocate memory for t");

  unsigned int r;
  // copy array to t
  for (r = 0; r < n; r++)
  {
    t[r] = v(r);
  }
  //for (int i = 0;i<n;i++) cout << " " << t[i];  cout << endl;

  T qs = RobustGaussian<T>::median(t, n);
  free(t);
  return qs;
}

/** Robust estimate for sigma (using median absolute deviation).
 \f$ 1.4826 med_i(|r_i - med_j(r_j)|) \f$
 */
template<class T>
T Regression<T>::getSigmaMAD(const vnl_vector<T>& v, T d)
{
  unsigned int n = v.size();
  T* t = (T *) calloc(n, sizeof(T));
  if (t == NULL)
    ErrorExit(ERROR_NO_MEMORY,
        "Regression<T>::getSigmaMAD could not allocate memory for t");

  unsigned int r;
  // copy array to t
  for (r = 0; r < n; r++)
  {
    t[r] = v[r];
  }
  //for (int i = 0;i<n;i++) cout << " " << t[i];  cout << endl;

  T qs = RobustGaussian<T>::mad(t, n, d);
  free(t);
  return qs;
}

/**
 \param fname is the filename without ending
 */
template<class T>
void Regression<T>::plotPartialSat(const std::string& fname)
{

  // residual is B (as b-Ap = b since p = 0 initially)

  // plot diffs
  string nbase = fname;
  int rf = nbase.rfind("/");
  if (rf != -1)
  {
    nbase = nbase.substr(rf + 1, nbase.length());
  }
  string fn = fname + ".plot";
  ofstream ofile(fn.c_str(), ios::out);
  bool png = false;
  if (png)
    ofile << "set terminal png medium size 800,600" << endl;
  else
    ofile << "set terminal postscript eps color" << endl;
  if (png)
    ofile << "set output \"" << nbase << ".png\"" << endl;
  else
    ofile << "set output \"" << nbase << ".eps\"" << endl;
  ofile << "plot ";
  ofile << " \"-\" notitle with lines 1" << endl;
  for (int j = 1; j <= 30; j++)
  {
    ofile << j << " " << getTukeyPartialSat(*b, j) << endl;
  }
  ofile << "e" << endl;
  ofile.close();

}
