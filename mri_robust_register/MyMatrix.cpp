/**
 * @file MyMatrix.cpp
 * @brief A static class with Matrix operations
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/08/13 02:51:19 $
 *    $Revision: 1.2 $
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
#include <cassert>

#include "MyMatrix.h"
#include "Quaternion.h"

using namespace std;


double MyMatrix::AffineTransDistSq(MATRIX * a, MATRIX * b, double r)
// computes squared distance between a and b (4x4 affine transformations)
// D^2 = 1/5 r^2 Tr(A^tA) +  ||T_d||^2
// where T_d is the translation and A the Affine
// of the matrix d = a - b
// r is the radius specifying the volume of interest
// (this distance is used in Jenkinson 1999 RMS deviation - tech report
//    www.fmrib.ox.ac.uk/analysis/techrep )
// the center of the brain should be at the origin
{

  MATRIX* drigid = MatrixCopy(a,NULL);
  if (b) drigid = MatrixSubtract(drigid,b,drigid);
  else
  {
    MATRIX *id = MatrixIdentity(4,NULL);
    drigid = MatrixSubtract(drigid,id,drigid);
    MatrixFree(&id);
  }

  double EPS = 0.000001;
  assert(drigid->rows ==4 && drigid->cols == 4);
  assert(fabs(drigid->rptr[4][1]) < EPS);
  assert(fabs(drigid->rptr[4][2]) < EPS);
  assert(fabs(drigid->rptr[4][3]) < EPS);

  //cout << " drigid: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",drigid);

  // translation norm quadrat:
  double tdq = 0;
  for (int i=1; i <= 3; i++)
  {
    tdq += drigid->rptr[i][4] * drigid->rptr[i][4];
    drigid->rptr[i][4] = 0.0;
    drigid->rptr[4][i] = 0.0;
  }
  drigid->rptr[4][4] = 0.0;

  //cout << " trans dist2: " << tdq << endl;
  MATRIX* dt = MatrixTranspose(drigid, NULL);
  drigid = MatrixMultiply(dt,drigid,drigid);
  MatrixFree(&dt);

  // Trace of A^t A
  double tr = 0.0;
  for (int i=1; i <= 3; i++)
  {
    tr += drigid->rptr[i][i];
  }

  MatrixFree(&drigid);

  return (1.0/5.0) * r*r* tr + tdq;

}

double MyMatrix::RigidTransDistSq(MATRIX * a, MATRIX * b)
// computes squared distance between a and b (4x4 rigid transformations)
// D^2 = ||T_d||^2 + |log R_d|^2
// where T_d is the translation and R_d the rotation
// of the matrix d that rigidly transforms a to b
{

  MATRIX* drigid;
  if (!b) drigid = MatrixCopy(a,NULL);
  else
  {
    drigid = MatrixInverse(a,NULL);
    drigid = MatrixMultiply(b,drigid,drigid);
  }

  double EPS = 0.000001;
  assert(drigid->rows ==4 && drigid->cols == 4);
  assert(fabs(drigid->rptr[4][1]) < EPS);
  assert(fabs(drigid->rptr[4][2]) < EPS);
  assert(fabs(drigid->rptr[4][3]) < EPS);
  assert(fabs(drigid->rptr[4][4]-1) < EPS);

  //cout << " drigid: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",drigid);

  // translation norm quadrat:
  double tdq = 0;
  for (int r=1; r <= 3; r++)
  {
    tdq += drigid->rptr[r][4] * drigid->rptr[r][4];
  }

  //cout << " trans dist2: " << tdq << endl;

  // rotation norm:
  double rd = RotMatrixLogNorm(drigid);
  //cout << " rd: " << rd << endl;

  MatrixFree(&drigid);

  return rd*rd + tdq;

}

double MyMatrix::getFrobeniusDiff(MATRIX *m1, MATRIX *m2)
{

  double s,ss = 0.0;
  assert(m1->rows == m2->rows);
  assert(m1->cols == m2->cols);

  for (int r = 1;r<=m1->rows;r++)
    for (int c = 1;c<=m2->cols;c++)
    {
      s = *MATRIX_RELT(m1, r, c) - *MATRIX_RELT(m2, r, c);
      ss += s * s;
    }
  ss = sqrt(ss);
  return ss;
}


MATRIX * MyMatrix::MatrixSqrt (MATRIX * m, MATRIX *msqrt)
{
  assert(m->rows == 4 && m->cols == 4);
  msqrt = MatrixIdentity(4,msqrt);
  MATRIX* R =  MatrixAlloc(3,3,MATRIX_REAL);
  for (int rr = 1; rr<=3; rr++)
    for (int cc = 1; cc<=3; cc++)
    {
      *MATRIX_RELT(R, rr, cc) = *MATRIX_RELT(m, rr, cc);
    }


  //Denman and Beavers square root iteration

  int imax = 100;
  double eps = 0.0001;
  double err = 1000;
  //cout << "using square root iteartion (" << imax << ")"<< endl;
  MATRIX * Yn  = MatrixCopy(R,NULL);
  MATRIX * Zn  = MatrixIdentity(3,NULL);
  MATRIX * Zni = NULL;
  MATRIX * Yni = NULL;
  MATRIX * Ysq = NULL;
  int count = 0;
  while (count<imax && err > eps)
  {
    count++;
    Yni = MatrixInverse(Yn,Yni);
    Zni = MatrixInverse(Zn,Zni);
    assert(Yni && Zni);

    Yn = MatrixAdd(Yn,Zni,Yn);
    Zn = MatrixAdd(Zn,Yni,Zn);

    Yn = MatrixScalarMul(Yn,0.5,Yn);
    Zn = MatrixScalarMul(Zn,0.5,Zn);
    //cout << " matrix " << i << endl;
    //MatrixPrintFmt(stdout,"% 2.8f",Yn);
    //cout << endl;

    Ysq = MatrixMultiply(Yn,Yn,Ysq);
    Ysq = MatrixSubtract(Ysq,R,Ysq);
    err = 0;
    for (int c=1; c<4; c++)
      for (int r=1; r<4; r++)
        err += fabs(*MATRIX_RELT(Ysq, r, c)) ;

  }

  if (count > imax)
  {
    cerr << "Matrix Sqrt did not converge in " << imax << " steps!" << endl;
    cerr << "   ERROR: " << err << endl;
    assert(err <= eps);
  }

  MATRIX * Rh = Yn;
  //cout << "rh : " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",Rh);
  //cout << endl;
  MatrixFree(&Zni);
  MatrixFree(&Yni);
  MatrixFree(&Zn);
  MatrixFree(&Ysq);


  // compute new T
  MATRIX* Rh1 = MatrixCopy(Rh,NULL);
  *MATRIX_RELT(Rh1, 1, 1) =*MATRIX_RELT(Rh1, 1, 1) +1;
  *MATRIX_RELT(Rh1, 2, 2) =*MATRIX_RELT(Rh1, 2, 2) +1;
  *MATRIX_RELT(Rh1, 3, 3) =*MATRIX_RELT(Rh1, 3, 3) +1;

  VECTOR * T = MatrixAlloc(3,1,MATRIX_REAL);
  *MATRIX_RELT(T, 1, 1) =*MATRIX_RELT(m, 1, 4);
  *MATRIX_RELT(T ,2, 1) =*MATRIX_RELT(m, 2, 4);
  *MATRIX_RELT(T, 3, 1) =*MATRIX_RELT(m, 3, 4);

  MATRIX* Rh1i = MatrixInverse(Rh1,NULL);
  assert(Rh1i);

  VECTOR * Th = MatrixMultiply(Rh1i,T,NULL);

  *MATRIX_RELT(msqrt, 1, 4) =*MATRIX_RELT(Th, 1, 1) ;
  *MATRIX_RELT(msqrt, 2, 4) =*MATRIX_RELT(Th, 2, 1) ;
  *MATRIX_RELT(msqrt, 3, 4) =*MATRIX_RELT(Th, 3, 1) ;
  *MATRIX_RELT(msqrt, 4, 1) =0 ;
  *MATRIX_RELT(msqrt, 4, 2) =0 ;
  *MATRIX_RELT(msqrt, 4, 3) =0 ;
  *MATRIX_RELT(msqrt, 4, 4) =1 ;
  for (int c=1; c<4; c++)
    for (int r=1; r<4; r++)
      *MATRIX_RELT(msqrt, r, c) = *MATRIX_RELT(Rh, r, c) ;

  MatrixFree(&Th);
  MatrixFree(&Rh1i);
  MatrixFree(&T);
  MatrixFree(&Rh1);
  MatrixFree(&Rh);
  MatrixFree(&R);

//    bool test = true;
//    if (test)
//    {
//       MATRIX* ms2 = MatrixMultiply(msqrt,msqrt,NULL);
//       ms2 = MatrixSubtract(ms2,m,ms2);
//       double sum = 0;
//       for (int c=1; c<=4; c++)
//       for (int r=1; r<=4; r++)
//         sum += fabs(*MATRIX_RELT(ms2, r, c)) ;
//       if (sum > 0.0001)
//       {
//          cerr << " Error : " << sum << endl;
//          //MatrixPrintFmt(stdout,"% 2.8f",ms2);
//          cerr << endl;
//   assert(1==2);
//       }
//       MatrixFree(&ms2);
//    }

  return msqrt;
}


MATRIX* MyMatrix::getMatrix(std::vector < double > d, int r, int c, MATRIX* m)
// convert double array to matrix
{
  if (c==-1) c=r; // quadratic

  assert(r*c == (int)d.size());
  if (!m) m = MatrixAlloc(r,c,MATRIX_REAL);
  assert(m->rows ==r);
  assert(m->cols ==c);

  int rr,cc, count=0;
  for (rr = 1 ; rr <= r ; rr++)
    for (cc = 1 ; cc <= c  ; cc++)
    {
      *MATRIX_RELT(m, rr,cc) = d[count];
      count++;
    }

  return m;
}

MATRIX * MyMatrix::aff2mat(MATRIX * aff, MATRIX *outM)
// converts affine vector (12x1)
// into an affine matrix (homogeneous coord) 4x4
{
  if (outM == NULL) outM = MatrixAlloc(4, 4, MATRIX_REAL);
  MatrixIdentity(4,outM);

  int count = 1;
  for (int rr = 1;rr<=3;rr++)
    for (int cc = 1;cc<=4;cc++)
    {
      *MATRIX_RELT(outM, rr, cc) = *MATRIX_RELT(outM, rr, cc) +  *MATRIX_RELT(aff, count, 1);
      count++;
    }

  return outM;
}

MATRIX * MyMatrix::getHalfRT (MATRIX * m, MATRIX *mhalf)
// m must be rotation and translation only (rigid!)
{
  if (mhalf) MatrixFree(&mhalf);

  float d = MatrixDeterminant(m);
  assert(fabs(d-1) < 0.000001);

  Quaternion q;
  q.importMatrix(*MATRIX_RELT(m, 1, 1),*MATRIX_RELT(m, 1, 2),*MATRIX_RELT(m, 1, 3),
                 *MATRIX_RELT(m, 2, 1),*MATRIX_RELT(m, 2, 2),*MATRIX_RELT(m, 2, 3),
                 *MATRIX_RELT(m, 3, 1),*MATRIX_RELT(m, 3, 2),*MATRIX_RELT(m, 3, 3));
  //cout << "q: "<< q << endl;
  Quaternion qh = q.getHalfRotation();
  //cout << "qh: " << qh << endl;
  mhalf  = getMatrix(qh.getRotMatrix3dh(),4,4);
  MATRIX* Rh1 = MatrixIdentity(3,NULL);
  for (int rr = 1; rr<4;rr++)
    for (int cc = 1; cc<4;cc++)
      *MATRIX_RELT(Rh1, rr, cc) = *MATRIX_RELT(Rh1, rr, cc) + *MATRIX_RELT(mhalf, rr, cc);

  VECTOR * T = MatrixAlloc(3,1,MATRIX_REAL);
  *MATRIX_RELT(T, 1, 1) = *MATRIX_RELT(m, 1, 4);
  *MATRIX_RELT(T ,2, 1) = *MATRIX_RELT(m, 2, 4);
  *MATRIX_RELT(T, 3, 1) = *MATRIX_RELT(m, 3, 4);

  //cout << " rh1" << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",Rh1);

  MATRIX* Rh1i = MatrixInverse(Rh1,NULL);
  assert(Rh1i);

  VECTOR * Th = MatrixMultiply(Rh1i,T,NULL);

  *MATRIX_RELT(mhalf, 1, 4) =*MATRIX_RELT(Th, 1, 1) ;
  *MATRIX_RELT(mhalf, 2, 4) =*MATRIX_RELT(Th, 2, 1) ;
  *MATRIX_RELT(mhalf, 3, 4) =*MATRIX_RELT(Th, 3, 1) ;
  *MATRIX_RELT(mhalf, 4, 1) =0 ;
  *MATRIX_RELT(mhalf, 4, 2) =0 ;
  *MATRIX_RELT(mhalf, 4, 3) =0 ;
  *MATRIX_RELT(mhalf, 4, 4) =1 ;

  MatrixFree(&Th);
  MatrixFree(&Rh1i);
  MatrixFree(&T);
  MatrixFree(&Rh1);
  return mhalf;
}

double  MyMatrix::RotMatrixLogNorm(MATRIX * m)
// computes Frobenius norm of log of rot matrix
// this is equivalent to geodesic distance on rot matrices
// will look only at first three rows and colums and
// expects a rotation matrix there
{
  // assert we have no stretching only rot (and trans)
  float det = MatrixDeterminant(m);
  //cout << " det: " << det << endl;
  if (fabs(det-1.0) > 0.001)
  {
    cerr << "There is streching! det: " << det << endl;
    assert (fabs(det-1.0) < 0.001);
  }

  double trace = 0.0;
  for (int n=1; n <= 3; n++) trace += m->rptr[n][n];
  //cout << " trace : " << trace << endl;
  trace = 0.5*(trace-1.0);
  if (trace > 1.0) trace = 1.0;
  if (trace < -1.0) trace = -1.0;
  //cout << "  0.5*(trace-1): " << trace << endl;
  double theta = acos(trace); // gives [0..pi]

  return sqrt(2.0) * theta;

}

double MyMatrix::RotMatrixGeoDist(MATRIX * a, MATRIX *b)
{

  if (!b) return RotMatrixLogNorm(a);

  // if not 3x3, fetch first 3x3
  // and construct a^T b
  MATRIX *at =  MatrixAlloc(3,3,MATRIX_REAL);
  MATRIX *blocal =  MatrixAlloc(3,3,MATRIX_REAL);
  assert (a->rows >= 3 && a->cols >= 3);
  assert (b->rows >= 3 && b->cols >= 3);
  for (int r=1; r <= 3; r++)
    for (int c=1; c <= 3; c++)
    {
      at->rptr[r][c] = a->rptr[c][r];
      blocal->rptr[r][c] = b->rptr[r][c];
    }

  blocal = MatrixMultiply(at,blocal,blocal);

  double dist = RotMatrixLogNorm(blocal);

  MatrixFree(&at);
  MatrixFree(&blocal);

  return dist;

}

