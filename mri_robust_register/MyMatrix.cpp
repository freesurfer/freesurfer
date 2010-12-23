/**
 * @file MyMatrix.cpp
 * @brief A static class with Matrix operations
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/12/23 00:08:40 $
 *    $Revision: 1.8 $
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
#include <vnl/vnl_inverse.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/algo/vnl_complex_eigensystem.h>
#include <vnl/vnl_real.h>
#include <vnl/vnl_imag.h>
#include <vnl/algo/vnl_qr.h>

using namespace std;


////// conversion stuff  (note vnl stuff is double, MATRIX is float !!!)

MATRIX* MyMatrix::convertVNL2MATRIX(const vnl_vector <double > & v, MATRIX* outM)
{
   if (outM == NULL) outM = MatrixAlloc(v.size(),1,MATRIX_REAL);
	 assert ((int) v.size() == outM->rows);
	 assert (1 == outM->cols);
	 
	 for (unsigned int r = 0; r< v.size(); r++)
	    outM->rptr[r+1][1] = v[r];
	 return outM;
}

MATRIX* MyMatrix::convertVNL2MATRIX(const vnl_matrix <double > & m, MATRIX* outM)
{
   if (outM == NULL) outM = MatrixAlloc(m.rows(),m.cols(),MATRIX_REAL);
	 assert ((int) m.rows() == outM->rows);
	 assert ((int) m.cols() == outM->cols);
	 
	 for (unsigned int r = 0; r< m.rows(); r++)
	 for (unsigned int c = 0; c< m.cols(); c++)
	    outM->rptr[r+1][c+1] = m[r][c];
	 return outM;
}

vnl_vector <double > MyMatrix::convertVECTOR2VNL(VECTOR* m)
{
   assert(m->cols == 1);
   vnl_vector <double>  ret(m->rows);
	 
	 for (int r = 0; r< m->rows; r++)
	    ret[r] = m->rptr[r+1][1];
	 return ret;
}

vnl_matrix <double > MyMatrix::convertMATRIX2VNL(MATRIX* m)
{
   vnl_matrix <double> ret(m->rows, m->cols);
	 
	 for (int r = 0; r< m->rows; r++)
	 for (int c = 0; c< m->cols; c++)
	    ret[r][c] = m->rptr[r+1][c+1];
	 return ret;
}


////// VNL stuff


std::pair < vnl_matrix < double >, vnl_matrix < double > > 
   MyMatrix::MatrixSqrtAndInv(const vnl_matrix < double >& m)
{
  assert(m.rows() == 4 && m.cols() == 4);
	
	// extract R and T from M:
	vnl_matrix_fixed < double,3,3 > R;// = m.extract(3,3,0,0);
  for (int rr = 0; rr<3; rr++)
    for (int cc = 0; cc<3; cc++)
    {
      R[rr][cc] = m[rr][cc];
    }
  vnl_vector_fixed < double,3 > T;
	T[0] = m[0][3]; T[1] = m[1][3]; T[2] = m[2][3];
	// compute rotation and translation of M^{-1}:
	vnl_matrix_fixed < double,3,3 > Ri=vnl_inverse(R);
  vnl_vector_fixed < double,3 > Ti=- Ri * T;
	// T, Ri and Ti are needed later 

  //Compute sqrt(R) with
  //Denman and Beavers square root iteration

  int imax = 100;
  double eps = 0.00001;  // important to be small to guarantee symmetry,
	                       // but even adding two zeros did not show any 
												 // differences in tests
	double err = 1000;
  //cout << "using square root iteartion (" << imax << ")"<< endl;
	vnl_matrix_fixed < double,3,3 > Yn(R);
	vnl_matrix_fixed < double,3,3 > Zn; Zn.set_identity();
	vnl_matrix_fixed < double,3,3 > Yni,Zni;
	vnl_matrix_fixed < double,3,3 > Ysq;	

  int count = 0;
  while (count<imax && err > eps)
  {
    count++;

    //store invrse here (we change Yn below)
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
		//cout << " iteration " << count << "  err: "<< err << endl;
  }
  // now Yn is sqrt(R) AND Zn is sqrt(R)^-1
	
  if (count > imax)
  {
    cerr << "Matrix Sqrt did not converge in " << imax << " steps!" << endl;
    cerr << "   ERROR: " << err << endl;
    assert(err <= eps);
  }

  // construct sqrt(M) = sqrt(R) x + Th (affine)
  // compute new Th
	// Rh1 = R + I
	vnl_matrix_fixed < double,3,3 > Rh1(Yn);
	Rh1[0][0] += 1; Rh1[1][1] +=1; Rh1[2][2] += 1;
  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  vnl_vector_fixed < double,3 > Th = vnl_inverse(Rh1) * T; //vnl_svd<double>(Rh1).solve(T);
  // put everything together:
  vnl_matrix < double > msqrt(4,4);
  msqrt[0][3] = Th[0];
	msqrt[1][3] = Th[1];
	msqrt[2][3] = Th[2];
	msqrt[3][0] = 0.0; msqrt[3][1] = 0.0; msqrt[3][2] = 0.0; msqrt[3][3] = 1.0;
  for (int c=0; c<3; c++)
    for (int r=0; r<3; r++)
		  msqrt[r][c] = Yn[r][c];

  // construct sqrt(M)^-1 = sqrt(R)-1 x + Thm (affine)
  // compute new Thm
	// Rh1 = R + I
	vnl_matrix_fixed < double,3,3 > Rhi1(Zn);
	Rhi1[0][0] += 1; Rhi1[1][1] +=1; Rhi1[2][2] += 1;
  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  vnl_vector_fixed < double,3 > Thi = vnl_inverse(Rhi1) * Ti; //vnl_svd<double>(Rh1).solve(T);
  // put everything together:
  vnl_matrix < double > msqrti(4,4);
  msqrti[0][3] = Thi[0];
	msqrti[1][3] = Thi[1];
	msqrti[2][3] = Thi[2];
	msqrti[3][0] = 0.0; msqrti[3][1] = 0.0; msqrti[3][2] = 0.0; msqrti[3][3] = 1.0;
  for (int c=0; c<3; c++)
    for (int r=0; r<3; r++)
		  msqrti[r][c] = Zn[r][c];

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

  return std::pair < vnl_matrix < double >, vnl_matrix < double > > (msqrt,msqrti);
}

vnl_matrix < double > MyMatrix::MatrixSqrt(const vnl_matrix < double >& m)
// for now separate the translation, else we get defective m
// (where we cannot do the eigendecomposition trick)
// in the future use schur decomposition (as in LAPACK) and
// then replicate the sqrtm from matlab
{
  assert(m.rows() == 4 && m.cols() == 4);

	vnl_matrix < double > R(3,3);// = m.extract(3,3,0,0);
  for (int rr = 0; rr<3; rr++)
    for (int cc = 0; cc<3; cc++)
    {
      R[rr][cc] = m[rr][cc];
    }

//cout << endl << endl;
//cout << " M: " << endl << m << endl;
//cout << " R: " << endl << R << endl;

  vnl_matrix < double > rsqrt (3,3,0.0);
  vnl_complex_eigensystem esys(R, rsqrt); //complex part is zero
	
//cout << " V' " << endl << esys.R << endl;
//cout << " D " << endl << esys.W << endl;

	
	vnl_diag_matrix < vcl_complex < double > > Wsqrt(3);
	for (unsigned int i=0;i<3;i++)
	{
	  Wsqrt[i] = sqrt(esys.W[i]);
	}
//cout << " Wsqrt " << endl << Wsqrt << endl;

  vnl_matrix < vcl_complex < double > > Rcomp (3,3);

//  esys.R.inplace_transpose(); // store evec in columns
//	vnl_matrix_fixed < vcl_complex < double >, 4, 4 > Rt = esys.R;
//  Rcomp = (esys.R * Wsqrt) * vnl_inverse(Rt);
  vnl_qr < vcl_complex < double > > QR(esys.R);
	Rcomp =  QR.solve(Wsqrt*esys.R);
	Rcomp.inplace_transpose();
	
//cout << " Rcomp " << endl << Rcomp << endl;
 
	
	rsqrt = vnl_real(Rcomp);
	
	
  // compute new T
	// Rh1 = R + I
	vnl_matrix_fixed < double,3,3 > Rh1(rsqrt);
	Rh1[0][0] += 1; Rh1[1][1] +=1; Rh1[2][2] += 1;

  vnl_vector_fixed < double,3 > T;
	T[0] = m[0][3]; T[1] = m[1][3]; T[2] = m[2][3];

  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  vnl_vector_fixed < double,3 > Th = vnl_inverse(Rh1) * T; //vnl_svd<double>(Rh1).solve(T);

  // put everything together:
  vnl_matrix < double > msqrt(4,4);
  msqrt[0][3] = Th[0];
	msqrt[1][3] = Th[1];
	msqrt[2][3] = Th[2];
	msqrt[3][0] = 0.0; msqrt[3][1] = 0.0; msqrt[3][2] = 0.0; msqrt[3][3] = 1.0;
  for (int c=0; c<3; c++)
    for (int r=0; r<3; r++)
		  msqrt[r][c] = rsqrt[r][c];

//cout << " msqrt " << endl << msqrt << endl;
	

  bool test = true;
  double eps = 0.00000000000001;  
  if (test)
  {
	  double fnorm = vnl_imag(Rcomp).frobenius_norm();
		//cout << " fnorm (img) = " << fnorm << endl;
	  if (fnorm > eps)
		{
		   cerr << " Error complex result?: " << fnorm << endl << vnl_imag(Rcomp) << endl;
			 assert(1==2);
		}
	
    vnl_matrix < double > ms2 = msqrt * msqrt;
    ms2 -= m;
    double sum = ms2.absolute_value_max();
		//cout << " max = " << sum << endl;
    if (sum > eps)
    {
      cerr << " Error : " << sum << endl;
      cerr << " sqrt(M): " << endl << msqrt << endl;
      cerr << endl;
      assert(1==2);
    }
  }
	
  return msqrt;
}

vnl_matrix < double > MyMatrix::MatrixSqrtIter(const vnl_matrix < double >& m)
{
  assert(m.rows() == 4 && m.cols() == 4);
	
	vnl_matrix_fixed < double,3,3 > R;// = m.extract(3,3,0,0);
  for (int rr = 0; rr<3; rr++)
    for (int cc = 0; cc<3; cc++)
    {
      R[rr][cc] = m[rr][cc];
    }


  //Denman and Beavers square root iteration

  int imax = 100;
  double eps = 0.00001;  // important to be small to guarantee symmetry,
	                       // but even adding two zeros did not show any 
												 // differences in tests
	double err = 1000;
  //cout << "using square root iteartion (" << imax << ")"<< endl;
	vnl_matrix_fixed < double,3,3 > Yn(R);
	vnl_matrix_fixed < double,3,3 > Zn; Zn.set_identity();
	vnl_matrix_fixed < double,3,3 > Yni,Zni;
	vnl_matrix_fixed < double,3,3 > Ysq;	

  int count = 0;
  while (count<imax && err > eps)
  {
    count++;

    //store invrse here (we change Yn below)
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
		//cout << " iteration " << count << "  err: "<< err << endl;
  }

  if (count > imax)
  {
    cerr << "Matrix Sqrt did not converge in " << imax << " steps!" << endl;
    cerr << "   ERROR: " << err << endl;
    assert(err <= eps);
  }

  // compute new T
	// Rh1 = R + I
	vnl_matrix_fixed < double,3,3 > Rh1(Yn);
	Rh1[0][0] += 1; Rh1[1][1] +=1; Rh1[2][2] += 1;

  vnl_vector_fixed < double,3 > T;
	T[0] = m[0][3]; T[1] = m[1][3]; T[2] = m[2][3];

  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  vnl_vector_fixed < double,3 > Th = vnl_inverse(Rh1) * T; //vnl_svd<double>(Rh1).solve(T);

  // put everything together:
  vnl_matrix < double > msqrt(4,4);
  msqrt[0][3] = Th[0];
	msqrt[1][3] = Th[1];
	msqrt[2][3] = Th[2];
	msqrt[3][0] = 0.0; msqrt[3][1] = 0.0; msqrt[3][2] = 0.0; msqrt[3][3] = 1.0;
  for (int c=0; c<3; c++)
    for (int r=0; r<3; r++)
		  msqrt[r][c] = Yn[r][c];


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

double MyMatrix::AffineTransDistSq(const vnl_matrix < double >&a, const vnl_matrix < double >&b, double r)
// computes squared distance between a and b (4x4 affine transformations)
// D^2 = 1/5 r^2 Tr(A^tA) +  ||T_d||^2
// where T_d is the translation and A the Affine
// of the matrix d = a - b
// r is the radius specifying the volume of interest
// (this distance is used in Jenkinson 1999 RMS deviation - tech report
//    www.fmrib.ox.ac.uk/analysis/techrep )
// the center of the brain should be at the origin
{

  vnl_matrix <double> drigid(a);
  if (!b.empty()) drigid -= b;
	else
	{
	   vnl_matrix<double> id(drigid.rows(),drigid.cols());
		 id.set_identity();
		 drigid -= id;
	}

  double EPS = 0.000001;
  assert(drigid.rows() ==4 && drigid.cols() == 4);
  assert(fabs(drigid[3][0]) < EPS);
  assert(fabs(drigid[3][1]) < EPS);
  assert(fabs(drigid[3][2]) < EPS);

  //cout << " drigid: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",drigid);

  // translation norm quadrat:
  double tdq = 0;
  for (int i=0; i < 3; i++)
  {
    tdq += drigid[i][3] * drigid[i][3];
    drigid[i][3] = 0.0;
    drigid[3][i] = 0.0;
  }
  drigid[3][3] = 0.0;

  //cout << " trans dist2: " << tdq << endl;
	

  drigid = drigid.transpose() * drigid;

  // Trace of A^t A
  double tr = 0.0;
  for (int i=0; i < 3; i++)
  {
    tr += drigid[i][i];
  }


  return (1.0/5.0) * r*r* tr + tdq;

}

double MyMatrix::RigidTransDistSq(const vnl_matrix < double >&a, const vnl_matrix < double >&b)
// computes squared distance between a and b (4x4 rigid transformations)
// D^2 = ||T_d||^2 + |log R_d|^2
// where T_d is the translation and R_d the rotation
// of the matrix d that rigidly transforms a to b
{

  vnl_matrix <double> drigid;
  if (b.empty()) drigid = a;
  else
  {
    drigid = vnl_matrix_inverse<double>(a);
    drigid = b * drigid;
  }

  double EPS = 0.000001;
  assert(drigid.rows() ==4 && drigid.cols() == 4);
  assert(fabs(drigid[3][0]) < EPS);
  assert(fabs(drigid[3][1]) < EPS);
  assert(fabs(drigid[3][2]) < EPS);
  assert(fabs(drigid[3][3]-1) < EPS);

  //cout << " drigid: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",drigid);

  // translation norm quadrat:
  double tdq = 0;
  for (int r=0; r < 3; r++)
  {
    tdq += drigid[r][3] * drigid[r][3];
  }

  //cout << " trans dist2: " << tdq << endl;

  // rotation norm:
  double rd = RotMatrixLogNorm(drigid);
  //cout << " rd: " << rd << endl;


  return rd*rd + tdq;

}

double MyMatrix::getFrobeniusDiff(const vnl_matrix < double >&m1, const vnl_matrix < double >&m2)
{

  assert(m1.rows() == m2.rows());
  assert(m1.cols() == m2.cols());
	
	return (m1-m2).fro_norm();
}

double  MyMatrix::RotMatrixLogNorm(const vnl_matrix_fixed < double, 4, 4 >& m)
// computes Frobenius norm of log of rot matrix
// this is equivalent to geodesic distance on rot matrices
// will look only at first three rows and colums and
// expects a rotation matrix there
{
  // assert we have no stretching only rot (and trans)
  double det = vnl_determinant(m);;
  //cout << " det: " << det << endl;
  if (fabs(det-1.0) > 0.001)
  {
    cerr << "There is streching! det: " << det << endl;
    assert (fabs(det-1.0) < 0.001);
  }

  double trace = 0.0;
  for (int n=0; n < 3; n++) trace += m[n][n];
  //cout << " trace : " << trace << endl;
  trace = 0.5*(trace-1.0);
  if (trace > 1.0) trace = 1.0;
  if (trace < -1.0) trace = -1.0;
  //cout << "  0.5*(trace-1): " << trace << endl;
  double theta = acos(trace); // gives [0..pi]

  return sqrt(2.0) * theta;

}

vnl_matrix < double > MyMatrix::getVNLMatrix(std::vector < double > d, int r)
// convert double array to matrix
{
  int c = (int)d.size() / r;
  assert(r*c == (int)d.size());
	
  vnl_matrix < double > m(r,c);
	
  int rr,cc, count=0;
  for (rr = 0 ; rr < r ; rr++)
    for (cc = 0 ; cc < c  ; cc++)
    {
      m[rr][cc] = d[count];
      count++;
    }

  return m;
}

void MyMatrix::getRTfromM(const vnl_matrix_fixed < double , 4 , 4 > &m,
	                              vnl_matrix_fixed < double , 3 , 3 > &r, 
														    vnl_vector_fixed < double, 3 >      &t)
{
	r = m.extract(3,3);
	t = m.extract(3,1,0,3).get_column(0);
}

vnl_matrix_fixed < double , 4 , 4 >  MyMatrix::getMfromRT(
	                             const vnl_matrix_fixed < double , 3 , 3 > &r, 
															 const vnl_vector_fixed < double, 3 >      &t)
{
  vnl_matrix_fixed < double, 4, 4 > m;
	m.set_identity();
	for (unsigned int rr = 0;rr<3;rr++)
	{
	  for (unsigned int cc = 0;cc<3;cc++)
		   m[rr][cc] = r[rr][cc];
		m[rr][3] = t[rr];
  }
	return m;
}

LTA* MyMatrix::VOXmatrix2LTA(const vnl_matrix_fixed < double, 4 , 4 >& m, MRI* src, MRI* dst)
{
  LTA* ret =  LTAalloc(1,src);
  ret->xforms[0].m_L = convertVNL2MATRIX(m,ret->xforms[0].m_L);
  ret->xforms[0].m_L = MRIvoxelXformToRasXform (src,dst,ret->xforms[0].m_L,ret->xforms[0].m_L) ;
  ret->type = LINEAR_RAS_TO_RAS ;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

LTA* MyMatrix::RASmatrix2LTA(const vnl_matrix_fixed < double, 4 , 4 >& m, MRI* src, MRI* dst)
{
  LTA* ret =  LTAalloc(1,src);
  ret->xforms[0].m_L = convertVNL2MATRIX(m,ret->xforms[0].m_L);
  ret->type = LINEAR_RAS_TO_RAS ;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

////// MATRIX stuff

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

pair < MATRIX *, VECTOR * > MyMatrix::getRTfromM(MATRIX * M, MATRIX * R, VECTOR * T)
{
   // check dimenstions:
   assert (M->rows == 4);
   assert (M->cols == 4);
   if (R == NULL) R = MatrixAlloc(3,3,M->type);
   if (T == NULL) T = VectorAlloc(3,M->type);
   assert(R->rows ==3);
   assert(R->cols ==3);
   assert(T->rows ==3);
   assert(T->cols ==1);
   
   // check M
   double eps = 0.000001;
   assert (fabs(M->rptr[4][1]) < eps);
   assert (fabs(M->rptr[4][2]) < eps);
   assert (fabs(M->rptr[4][3]) < eps);
   assert (fabs(M->rptr[4][4] - 1.0) < eps);
   
    for (int c=1; c<4;c++)
    {
      for (int r=1; r<4;r++)
      {
        R->rptr[r][c] = M->rptr[r][c];
      }
      T->rptr[c][1] = M->rptr[c][4];
   }
   
   return pair <MATRIX *, VECTOR *> (R,T);
}

MATRIX * MyMatrix::getMfromRT(MATRIX * R, VECTOR * T, MATRIX * M)
{
   int type;
   if (R != NULL) type = R->type;
   else if ( T != NULL) type = T->type;
   else assert(R != NULL || T != NULL);
   
   if (M == NULL ) M = MatrixAlloc(4,4,type);
   if (R == NULL) R = MatrixIdentity(3,NULL);
   if (T == NULL) T = MatrixZero(3,1,NULL);
   
   // check dimensions:
   assert (M->rows == 4);
   assert (M->cols == 4);
   assert(R->rows ==3);
   assert(R->cols ==3);
   assert(T->rows ==3);
   assert(T->cols ==1);
   
   
    for (int c=1; c<4;c++)
    {
      for (int r=1; r<4;r++)
      {
        M->rptr[r][c] = R->rptr[r][c];
      }
      M->rptr[c][4] = T->rptr[c][1];
      M->rptr[4][c] = 0;
    }
    M->rptr[4][4] = 1;
   
   return M;
}
LTA* MyMatrix::VOXmatrix2LTA(MATRIX * m, MRI* src, MRI* dst)
{
  LTA* ret =  LTAalloc(1,src);
  ret->xforms[0].m_L = MRIvoxelXformToRasXform (src,dst,m,NULL) ;
  ret->type = LINEAR_RAS_TO_RAS ;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

LTA* MyMatrix::RASmatrix2LTA(MATRIX * m, MRI* src, MRI* dst)
{
  LTA* ret =  LTAalloc(1,src);
  ret->xforms[0].m_L = MatrixCopy(m,ret->xforms[0].m_L) ;
  ret->type = LINEAR_RAS_TO_RAS ;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

