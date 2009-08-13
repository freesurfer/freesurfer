/**
 * @file Regression.cpp
 * @brief A class to solve overconstrained system A X = B
 *
 *   it uses either least squares (standard regression)
 *   or a robust estimator (Tukey's Biweight with iterative reweighted least
 *   squares)
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/08/13 02:51:19 $
 *    $Revision: 1.15 $
 *
 * Copyright (C) 2008-2009
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include "Regression.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>
#include <fstream>
#include "RobustGaussian.h"

#ifdef __cplusplus
extern "C"
{
#endif
  #include <stdio.h>
  #include <stdlib.h>
  #include "error.h"
#ifdef __cplusplus
}
#endif



using namespace std;

MATRIX * Regression::getRobustEst(double sat, double sig)
{
  pair < MATRIX* , MATRIX* > pw = getRobustEstW(sat,sig);
  MatrixFree(&pw.second);
  return pw.first;
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

pair < MATRIX *, MATRIX *> Regression::getRobustEstW(double sat, double sig)
{
  if (A) return getRobustEstWAB(sat,sig);
  return getRobustEstWB(sat,sig);
}

pair < MATRIX *, MATRIX *> Regression::getRobustEstWB(double sat, double sig)
{
  // constants
  int MAXIT = 20;
  //double EPS = 2e-12;
  double EPS = 2e-6;

  double muold =0;
  double mu =0;

  assert(B->cols == 1);
  assert(A == NULL);

  // initialize with median:
  double bb[B->rows];
  for (int i = 1 ; i<=B->rows; i++)
    bb[i-1] = B->rptr[i][1];
  mu = RobustGaussian::median(bb,B->rows);

  muold = mu + 100; // set out of reach of eps
  int count = 0;
  MATRIX * w = MatrixConstVal(1.0,B->rows,1,NULL);
  MATRIX * r = MatrixConstVal(1.0,B->rows,1,NULL);
  MATRIX * rdsigma = MatrixConstVal(1.0,B->rows,1,NULL);
  double sigma,d1,d2,d3,wi;

  //cout << endl<<endl<<" Values: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",B);

  while (fabs(muold - mu) > EPS && count < MAXIT)
  {
    count++;
    muold = mu;
    //cout << endl << "count " << count << " Myold " << muold << endl;
    for (int i = 1 ; i<=B->rows; i++)
      r->rptr[i][1] = muold - B->rptr[i][1];

    //cout << " residuals: " << endl;
    //   MatrixPrintFmt(stdout,"% 2.8f",r);

    sigma = getSigmaMAD(r);
    //cout << " sigma: " << sigma << endl;
    if (sigma < EPS) // if all r are the same (usually zero)
    {
      mu = muold;
      for (int i=1 ; i<=w->rows; i++)
        w->rptr[i][1] = 1;
      break;
    }

    for (int i = 1 ; i<=B->rows; i++)
      rdsigma->rptr[i][1] = r->rptr[i][1] / sigma;
    //cout << " r/sigma: " << endl;
    //   MatrixPrintFmt(stdout,"% 2.8f",rdsigma);

    getSqrtTukeyDiaWeights(rdsigma, sat, w); // here we get sqrt of weights into w
    //cout << " weights: " << endl;
    //   MatrixPrintFmt(stdout,"% 2.8f",w);

    d1 = 0;
    d2 = 0;
    d3 = 0;
    for (int i = 1 ; i<=B->rows; i++)
    {
      wi = w->rptr[i][1] *w->rptr[i][1];
      d1 += wi * r->rptr[i][1];
      d2 += wi;
      d3 += wi * r->rptr[i][1] * r->rptr[i][1];
    }

    if (d2 < EPS) // weights zero
    {

      assert(d2 >= EPS);
    }
    mu = muold - d1/d2;

  }
//cout << "!!! final mu :  " << mu << endl;
  MatrixFree(&r);
  MatrixFree(&rdsigma);

  MATRIX * p = MatrixConstVal(mu,1,1,NULL);

  return pair < MATRIX*, MATRIX*> (p,w);

}

pair < MATRIX *, MATRIX *> Regression::getRobustEstWAB(double sat, double sig)
// solves overconstrained system A p = b using
// M estimators (Beaton and Tukey's biweigt function)
// returns pair < p, w  > of vectors
// where p is the robust solution and w the used weights
// Method: iterative reweighted least squares
{
//   cout << "  Regression::getRobustEstW( "<<sat<<" , "<<sig<<" ) " << endl;

  // constants
  int MAXIT = 20;
  //double EPS = 2e-16;
  double EPS = 2e-12;

  // variables
  vector < double > err(MAXIT+1);
  err[0] = numeric_limits<double>::infinity();
  err[1] = 1e20;
  int count = 1, rr, cc;
  double sigma;

  MATRIX * w = MatrixConstVal(1.0,A->rows,1,NULL);
  MATRIX * r = MatrixConstVal(1.0,A->rows,1,NULL);
  MATRIX * p = MatrixZero(A->cols, 1, NULL);
    
  if (w == NULL || r == NULL || p == NULL)
     ErrorExit(ERROR_NO_MEMORY,"Regression::getRobustEstWAB could not allocate memory for w,r,p") ;

  MATRIX * lastp = MatrixCopy(p,NULL);
  MATRIX * lastw = MatrixCopy(w,NULL);
  MATRIX * lastr = MatrixCopy(B,NULL);  // error lastr = b-A*p = b  , since p = 0 initially
  if (lastw == NULL || lastr == NULL || lastp == NULL)
     ErrorExit(ERROR_NO_MEMORY,"Regression::getRobustEstWAB could not allocate memory for lastw,lastr,lastp") ;

  MATRIX * wAi = NULL, *v = NULL;
  MATRIX * wA = MatrixAlloc(A->rows, A->cols, MATRIX_REAL);
  MATRIX * wr = MatrixAlloc(A->rows, 1, MATRIX_REAL);
  if (wA == NULL || wr == NULL )
     ErrorExit(ERROR_NO_MEMORY,"Regression::getRobustEstWAB could not allocate memory for wA,wr") ;

  // compute error (lastr = b-A*p)
  //v     = MatrixMultiply(A,p,NULL);
  //lastr = MatrixSubtract(B,v,lastr);
  //MatrixFree(&v); v = NULL;

  // iteration until we increase the error, we reach maxit or we have no error
  while ( err[count-1] > err[count] && count < MAXIT && err[count] > EPS )
  {
    // save previous values   (not necessary for first step)
    MatrixCopy(p, lastp);
    MatrixCopy(w, lastw);
    MatrixCopy(r, lastr);

    count++;

    //cout << " count: "<< count << endl;
    //cout << " memusage: " << endl;
    //print_mem_usage();
    
    
    // recompute weights
    if (count > 2)
    {
      sigma = getSigmaMAD(lastr);
      //cout << " sigma: " << sigma << endl;
      MatrixScalarMul(lastr,(float)1.0/sigma,r); // use r temporarily here
      getSqrtTukeyDiaWeights(r, sat, w); // here we get sqrt of weights into w
    }


    // solve WLS with the new weights
    // compute wA  where w = sqrt(W);
    for (rr = 1;rr<=wA->rows;rr++)
      for (cc = 1;cc<=wA->cols;cc++)
        wA->rptr[rr][cc] = A->rptr[rr][cc] * w->rptr[rr][1];

    // compute wr  = ((b-A*p).*w);
    for (rr = 1;rr<=wr->rows;rr++)
      wr->rptr[rr][1] = lastr->rptr[rr][1] * w->rptr[rr][1];

    // compute new p = lastp + pinv(wA)*wr
    //wAi = MatrixPseudoInverse(wA,wAi);
    if (wAi != NULL) MatrixFree(&wAi);
    wAi = MatrixSVDPseudoInverse(wA,NULL);

    if (wAi == NULL)
    {
      cerr << "    Regression::getRobustEstW  could not compute pseudo inverse!" << endl;
      //cerr << " wa: " << endl ;
      //   MatrixPrintFmt(stdout,"% 2.8f",wAi);
      //cerr << endl;
      assert(wAi != NULL);
    }
    //cout << " got inverse" << endl;
    // MatrixPrintFmt(stdout,"% 2.8f",wAi);
    // cout << endl;
    // MatrixPrintFmt(stdout,"% 2.8f",wr);
    // cout << endl;

    v   = MatrixMultiply(wAi,wr,v);
    // cout << " asdf" << endl;
    MatrixAdd(lastp,v,p);


    // compute new residual and total errors (using new p)
    // r = b - (A*p)
    MatrixMultiply(A,p,r);
    MatrixSubtract(B,r,r);
    // err = sum (w r^2) / sum (w)
    float swr = 0;
    float sw  = 0;
    int rr;
    for (rr = 1;rr<=r->rows;rr++)
    {
      float t1 = w->rptr[rr][1];
      float t2 = r->rptr[rr][1];
      t1 *= t1; // remember w is the sqrt of the weights
      t2 *= t2;
      sw  += t1;
      swr += t1*t2;
    }
    err[count] = swr/sw;

//    cout << " Step: " << count-1 << " ERR: "<< err[count]<< endl;
//    cout << " p  : "<< endl;
//    MatrixPrintFmt(stdout,"% 2.8f",p);
//    cout << endl;

    // do not save values here, to be able to use lastp when exiting earlier due to increasing error
    // MatrixCopy(p, lastp);
    // MatrixCopy(w, lastw);
    // MatrixCopy(r, lastr);

  }

  if (err[count]>err[count-1])
  {
    // take previous values (since actual values made the error to increase)
    // cout << " last step was no improvement, taking values : "<<  count-1 << endl;
    MatrixFree(&p);
    MatrixFree(&w);
    p = lastp;
    w = lastw;
    cout << "     Step: " << count-2 << " ERR: "<< err[count-1]<< endl;
    lasterror = err[count-1];
  }
  else
  {
    MatrixFree(&lastp);
    MatrixFree(&lastw);
    cout << "     Step: " << count-1 << " ERR: "<< err[count]<< endl;
    lasterror = err[count];
  }

  MatrixFree(&lastr);
  MatrixFree(&r);
  MatrixFree(&wA);
  MatrixFree(&wr);
  MatrixFree(&wAi);
  MatrixFree(&v);

//    cout << " p-final  : "<< endl;
//    MatrixPrintFmt(stdout,"% 2.8f",p);

  double d=0.0;
  double dd=0.0;
  double ddcount = 0;
  int zcount =  0;
  double val;
  for (int i = 1;i<=w->rows;i++)
  {
    val = (double) *MATRIX_RELT(w, i, 1);
    d +=val ;
    if (fabs(*MATRIX_RELT(B, i, 1)) > 0.00001)
    {
      dd+= val ;
      ddcount++;
      if (val < 0.1) zcount++;
    }
  }
  d /= w->rows;
  dd /= ddcount;
  cout << "          weights average: " << dd << "  zero: " << (double)zcount/ddcount << flush;
  //"  on significant b vals ( " << ddcount << " ): " << dd <<endl;
  lastweight = dd;
  lastzero   = (double)zcount/ddcount;
  return pair < MATRIX*, MATRIX*> (p,w);

}

MATRIX* Regression::getLSEst(MATRIX* outX)
{
  //cout << " Regression::getLSEst " << endl;
  lastweight = -1;
  lastzero   = -1;

  if (A==NULL) // LS solution is just the mean of B
  {
    if (!outX) outX = MatrixAlloc(1,1,MATRIX_REAL);
    assert(outX->rows == 1);
    assert(outX->cols == 1);
    assert(B->cols ==1);
    outX->rptr[1][1] = 0;
    for (int r=1;r<=B->rows;r++)
      outX->rptr[1][1] += B->rptr[r][1];
    outX->rptr[1][1] = outX->rptr[1][1] / B->rows;
    return outX;

  }

  MATRIX* Ai = MatrixSVDPseudoInverse(A,NULL);
  if (Ai == NULL)
  {
    cerr << "    Regression::getLSEst   could not compute pseudo inverse!" << endl;
    //cerr << "    A: " << endl ;
    //MatrixPrintFmt(stdout,"% 2.8f",A);
    //cerr << endl;
    assert(Ai != NULL);
  }
  outX   = MatrixMultiply(Ai,B,outX);
  MatrixFree(&Ai);

  // compute error:
  MATRIX* R = MatrixMultiply(A,outX,NULL);
  R = MatrixSubtract(R,B,R);
  double serror = 0;
  int cc, rr;
  for (cc = 1;cc<=R->cols;cc++)
    for (rr = 1;rr<=R->rows;rr++)
    {
      serror += R->rptr[rr][cc] * R->rptr[rr][cc];
    }
  MatrixFree(&R);
//  cout << "     squared error: " << serror << endl;
  lasterror = serror;

  return outX;
}

MATRIX* getTukeyBiweight(MATRIX* r, double sat, MATRIX* w)
{
  int n = r->rows;
  if (w == NULL) w = MatrixAlloc(n,1,  MATRIX_REAL);
  else (assert (n == w->rows));
  assert(r->cols ==1 );
  if (w == NULL) 
     ErrorExit(ERROR_NO_MEMORY,"Regression::getTukeyBiweight could not allocate memory for w") ;

  double a , b;
  for (int i = 0;i<n;i++)
  {
    if (fabs(r->rptr[i][1]) > sat) w->rptr[i][1] = sat*sat/2.0;
    else
    {

      a = r->rptr[i][1]/sat;
      b = 1.0 - a * a;
      w->rptr[i][1] = (sat*sat/2.0)*(1.0-b*b*b);
    }
  }
  return w;
}

double Regression::getTukeyPartialSat(MATRIX* r, double sat)
// computes sum of partial derivatives d_rho/d_sat(r_i)
{

  assert(r->cols ==1 );

  double sum = 0;
  double rt2;
  double rt4;
  double sat3 = sat *sat*sat;
  double sat5 = sat3 *sat*sat;
  int rr;
  for (rr = 1;rr<=r->rows;rr++)
  {
    if (fabs(r->rptr[rr][1]) > sat)
    {
      sum += sat;
    }
    else
    {
      rt2 = r->rptr[rr][1] * r->rptr[rr][1];
      rt4 = rt2 *rt2;
      sum += 3.0 * rt4/sat3 - 2.0 * rt4 *rt2 / sat5;
    }
  }

  return sum;
}

MATRIX * Regression::getSqrtTukeyDiaWeights(MATRIX * r, double sat, MATRIX * w)
// computes sqrt(weights) for a reweighted least squares approach
// returns elements of diag matrix W (as a column vector)
{
  //cout << " getTukeyDiaWeights  r size: " << r->rows << " , " << r->cols << endl;

  assert(r->cols == 1);
  int n = r->rows;
  if (w == NULL) w = MatrixAlloc(n,1,  MATRIX_REAL);
  else (assert (r->rows == w->rows));
  assert(r->cols ==1 );
  if (w == NULL) 
     ErrorExit(ERROR_NO_MEMORY,"Regression::getTukeyDiaWeights could not allocate memory for w") ;

  double t1, t2;
  int rr,cc;
  int count = 1;
  int ocount = 0;
  for (rr = 1;rr<=r->rows;rr++)
    for (cc = 1;cc<=r->cols;cc++)
    {
      // cout << " fabs: " << fabs(r->rptr[rr][cc]) << " sat: " << sat << endl;
      if (fabs(r->rptr[rr][cc]) > sat)
      {
        w->rptr[count][1] = 0;
        ocount++;
      }
      else
      {
        t1 = r->rptr[rr][cc]/sat;
        t2 = 1.0 - t1 *t1;
        //w->rptr[count][1] = t2 * t2;
        w->rptr[count][1] = t2 ; // returning sqrt
      }
      count++;
    }
  //cout << " over threshold: " << ocount << " times ! " << endl;
  return w;
}


double Regression::VectorMedian(MATRIX *v)
{
  int n = v->rows * v->cols;
  double* t = (double *)calloc(n, sizeof(double));
  if (t == NULL) 
     ErrorExit(ERROR_NO_MEMORY,"Regression::VectorMedian could not allocate memory for t") ;

  int r,c,cc;
  // copy array to t
  cc = 0;
  for (r=1; r <= v->rows; r++)
    for (c=1; c <= v->cols; c++)
    {
      t[cc] = v->rptr[r][c];
      cc++;
    }
  //for (int i = 0;i<n;i++) cout << " " << t[i];  cout << endl;

  double qs = RobustGaussian::median(t,n) ;
  free(t);
  return qs;

}

double Regression::getSigmaMAD(MATRIX *v, double d)
// robust estimate for sigma (using median absolute deviation)
// 1.4826 med_i(|r_i - med_j(r_j)|)
{
  int n = v->rows * v->cols;
  double* t = (double *)calloc(n, sizeof(double));
  if (t == NULL) 
     ErrorExit(ERROR_NO_MEMORY,"Regression::getSigmaMAD could not allocate memory for t") ;

  int r,c,cc;
  // copy array to t
  cc = 0;
  for (r=1; r <= v->rows; r++)
    for (c=1; c <= v->cols; c++)
    {
      t[cc] = v->rptr[r][c];
      cc++;
    }
  //for (int i = 0;i<n;i++) cout << " " << t[i];  cout << endl;

  double qs = RobustGaussian::mad(t,n,d) ;
  free(t);
  return qs;
}


void Regression::plotPartialSat(const std::string& fname)
// fname is the filename without ending
{

  // residual is B (as b-Ap = b since p = 0 initially)

  // plot diffs
  string nbase = fname;
  int rf = nbase.rfind("/");
  if (rf != -1)
  {
    nbase = nbase.substr(rf+1,nbase.length());
  }
  string fn = fname+".plot";
  ofstream ofile(fn.c_str(),ios::out);
  bool png = false;
  if (png) ofile << "set terminal png medium size 800,600" << endl;
  else ofile << "set terminal postscript eps color" << endl;
  if (png) ofile << "set output \""<< nbase <<".png\"" << endl;
  else ofile << "set output \""<< nbase <<".eps\"" << endl;
  ofile << "plot ";
  ofile << " \"-\" notitle with lines 1" << endl;
  for (int j = 1;j<=30; j++)
  {
    ofile << j << " " << getTukeyPartialSat(B, j) << endl;
  }
  ofile << "e" << endl;
  ofile.close();

}
