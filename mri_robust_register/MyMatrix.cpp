/**
 * @brief A static class with Matrix operations
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
#include <cassert>
#include <algorithm>
#include "MyMatrix.h"
#include "Quaternion.h"
#include "utils.h" //nint

#define export // obsolete feature 'export template' used in these headers 
#include <vnl/vnl_inverse.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/algo/vnl_complex_eigensystem.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_real.h>
#include <vnl/vnl_imag.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/vnl_complexify.h>
#include <vcl_iostream.h>
#include <vnl/vnl_matlab_print.h>
#undef export

#include "matrix.h"

#define sign(x) (( x > 0 ) - ( x < 0 ))

using namespace std;

extern "C" {

// complex Schur decomposition
extern void zgees_(char *jobvs, char *sort, long (*select)(), long *n,
    complex<double> *a, long *lda, long *sdim, complex<double> *w,
    complex<double> *vs, long *ldvs, complex<double> *work, long *lwork,
    double *rwork, long *bwork, long *info);

// complex Schur reordering
extern void ztrsen_(char *job, char *compq, int *select, long *n,
    complex<double> *t, long *ldt, complex<double> *q, long *ldq,
    complex<double> *w, long *m, double *s, double *sep, complex<double> *work,
    long *lwork, long *info);

}

// void PrintComplex (const vnl_matrix < vcl_complex < double > > &M)
// {
// 
//   char signum[2]={ '+','-'};
// 
//   cout << " [ ";
//   
//   for (unsigned int r=0;r<M.rows();r++)
//   {
//     for (unsigned int c=0;c<M.cols();c++)
//       cout << " ( " << M[r][c].real() << " " << signum[(M[r][c].imag()<0)]<< " " << M[r][c].imag() << "i ) " ;
//     if (r < M.rows() -1 ) cout << " ; " << endl;
//   }
//   cout << " ] " << endl;
// 
// }

// extern "C" {
// /*: Computes Schur Decomposistion of nxn complex matrix */
// extern int v3p_netlib_zgees_(
//   char v3p_netlib_const *jobvs,
//   char v3p_netlib_const *sort,
//   v3p_netlib_integer (*select)(),
//   v3p_netlib_integer v3p_netlib_const *n,
//   v3p_netlib_doublecomplex *a,
//   v3p_netlib_integer v3p_netlib_const *lda,
//   v3p_netlib_integer *sdim,
//   v3p_netlib_doublecomplex *w,
//   v3p_netlib_doublecomplex *vs,
//   v3p_netlib_integer v3p_netlib_const *ldvs,
//   v3p_netlib_doublecomplex *work,
//   v3p_netlib_integer *lwork,
//   v3p_netlib_doublereal *rwork,
//   v3p_netlib_integer *bwork,
//   v3p_netlib_integer *info
//   );
// }

////// conversion stuff  (note vnl stuff is double, MATRIX is float !!!)

MATRIX* MyMatrix::convertVNL2MATRIX(const vnl_vector<double> & v, MATRIX* outM)
{
  if (outM == NULL)
    outM = MatrixAlloc(v.size(), 1, MATRIX_REAL);
  assert((int) v.size() == outM->rows);
  assert(1 == outM->cols);

  for (unsigned int r = 0; r < v.size(); r++)
    outM->rptr[r + 1][1] = v[r];
  return outM;
}

MATRIX* MyMatrix::convertVNL2MATRIX(const vnl_matrix<double> & m, MATRIX* outM)
{
  if (outM == NULL)
    outM = MatrixAlloc(m.rows(), m.cols(), MATRIX_REAL);
  assert((int) m.rows() == outM->rows);
  assert((int) m.cols() == outM->cols);

  for (unsigned int r = 0; r < m.rows(); r++)
    for (unsigned int c = 0; c < m.cols(); c++)
      outM->rptr[r + 1][c + 1] = m[r][c];
  return outM;
}

vnl_vector<double> MyMatrix::convertVECTOR2VNL(VECTOR* m)
{
  assert(m->cols == 1);
  vnl_vector<double> ret(m->rows);

  for (int r = 0; r < m->rows; r++)
    ret[r] = m->rptr[r + 1][1];
  return ret;
}

vnl_matrix<double> MyMatrix::convertMATRIX2VNL(MATRIX* m)
{
  vnl_matrix<double> ret(m->rows, m->cols);

  for (int r = 0; r < m->rows; r++)
    for (int c = 0; c < m->cols; c++)
      ret[r][c] = m->rptr[r + 1][c + 1];
  return ret;
}

////// VNL stuff

std::pair<vnl_matrix<double>, vnl_matrix<double> > MyMatrix::MatrixSqrtAndInvIter(
    const vnl_matrix<double>& m)
//Compute sqrt(R) with
//Denman and Beavers square root iteration
{
  assert(m.rows() == 4 && m.cols() == 4);

  // extract R and T from M:
  vnl_matrix_fixed<double, 3, 3> R; // = m.extract(3,3,0,0);
  for (int rr = 0; rr < 3; rr++)
    for (int cc = 0; cc < 3; cc++)
    {
      R[rr][cc] = m[rr][cc];
    }
  vnl_vector_fixed<double, 3> T;
  T[0] = m[0][3];
  T[1] = m[1][3];
  T[2] = m[2][3];
  // compute rotation and translation of M^{-1}:
  vnl_matrix_fixed<double, 3, 3> Ri = vnl_inverse(R);
  vnl_vector_fixed<double, 3> Ti = -Ri * T;
  // T, Ri and Ti are needed later 

  //Compute sqrt(R) with
  //Denman and Beavers square root iteration

  int imax = 100;
  double eps = 0.00001;  // important to be small to guarantee symmetry,
                         // but even adding two zeros did not show any 
                         // differences in tests
  double err = 1000;
  //cout << "using square root iteartion (" << imax << ")"<< endl;
  vnl_matrix_fixed<double, 3, 3> Yn(R);
  vnl_matrix_fixed<double, 3, 3> Zn;
  Zn.set_identity();
  vnl_matrix_fixed<double, 3, 3> Yni, Zni;
  vnl_matrix_fixed<double, 3, 3> Ysq;

  int count = 0;
  while (count < imax && err > eps)
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
  vnl_matrix_fixed<double, 3, 3> Rh1(Yn);
  Rh1[0][0] += 1;
  Rh1[1][1] += 1;
  Rh1[2][2] += 1;
  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  vnl_vector_fixed<double, 3> Th = vnl_inverse(Rh1) * T; //vnl_svd<double>(Rh1).solve(T);
  // put everything together:
  vnl_matrix<double> msqrt(4, 4);
  msqrt[0][3] = Th[0];
  msqrt[1][3] = Th[1];
  msqrt[2][3] = Th[2];
  msqrt[3][0] = 0.0;
  msqrt[3][1] = 0.0;
  msqrt[3][2] = 0.0;
  msqrt[3][3] = 1.0;
  for (int c = 0; c < 3; c++)
    for (int r = 0; r < 3; r++)
      msqrt[r][c] = Yn[r][c];

  // construct sqrt(M)^-1 = sqrt(R)-1 x + Thm (affine)
  // compute new Thm
  // Rh1 = R + I
  vnl_matrix_fixed<double, 3, 3> Rhi1(Zn);
  Rhi1[0][0] += 1;
  Rhi1[1][1] += 1;
  Rhi1[2][2] += 1;
  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  vnl_vector_fixed<double, 3> Thi = vnl_inverse(Rhi1) * Ti; //vnl_svd<double>(Rh1).solve(T);
  // put everything together:
  vnl_matrix<double> msqrti(4, 4);
  msqrti[0][3] = Thi[0];
  msqrti[1][3] = Thi[1];
  msqrti[2][3] = Thi[2];
  msqrti[3][0] = 0.0;
  msqrti[3][1] = 0.0;
  msqrti[3][2] = 0.0;
  msqrti[3][3] = 1.0;
  for (int c = 0; c < 3; c++)
    for (int r = 0; r < 3; r++)
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
//          cerr << " sqrt(M): " << endl << msqrt << endl;
//          cerr << endl;
//          assert(1==2);
//       }
//    }

  return std::pair<vnl_matrix<double>, vnl_matrix<double> >(msqrt, msqrti);
}

vnl_matrix<double> MyMatrix::PadeApproximantOfDegree(
    const vnl_matrix<double> & A, unsigned int m)
//PADEAPPROXIMANTOFDEGREE  Pade approximant to exponential.
//   F = PADEAPPROXIMANTOFDEGREE(A,m) is the degree m diagonal
//   Pade approximant to EXP(A), where m = 3, 5, 7, 9 or 13.
//   Series are evaluated in decreasing order of powers, which is
//   in approx. increasing order of maximum norms of the terms.
{
  assert(A.cols() == A.rows());
  int n = A.cols();
  std::vector<double> c = getPadeCoefficients(m);

  vnl_matrix<double> F(n, n);

  // Evaluate Pade approximant.
  switch (m)
  {

  case 3:
  case 5:
  case 7:
  case 9:
  {
    int s = (int) ceil((m + 1.0) / 2.0);
    std::vector<vnl_matrix<double> > Apowers(s, vnl_matrix<double>(n, n));
    Apowers[0].set_identity();
    Apowers[1] = A * A;
    for (int j = 2; j < s; j++)
      Apowers[j] = Apowers[j - 1] * Apowers[1];

    vnl_matrix<double> U(n, n, 0.0);
    vnl_matrix<double> V(n, n, 0.0);

    for (int j = m; j >= 1; j -= 2)
      U = U + c[j] * Apowers[j / 2];

    U = A * U;
    for (int j = m - 1; j >= 0; j -= 2)
      V = V + c[j] * Apowers[(j + 1) / 2];

//        F = (-U+V)\(U+V);!! 
    F = vnl_svd<double>(-U + V).solve(U + V);

  }
    break;
  case 13:
  {
    // For optimal evaluation need different formula for m >= 12.
    vnl_matrix<double> A2(A * A);
    vnl_matrix<double> A4(A2 * A2);
    vnl_matrix<double> A6(A2 * A4);
    vnl_matrix<double> id(n, n);
    id.set_identity();
    vnl_matrix<double> U(
        A
            * (A6 * (c[13] * A6 + c[11] * A4 + c[9] * A2) + c[7] * A6
                + c[5] * A4 + c[3] * A2 + c[1] * id));

    vnl_matrix<double> V(
        A6 * (c[12] * A6 + c[10] * A4 + c[8] * A2) + c[6] * A6 + c[4] * A4
            + c[2] * A2 + c[0] * id);
//      F = (-U+V)\(U+V);!!
    F = vnl_svd<double>(-U + V).solve(U + V);
  }
    break;
  }

  return F;
}

std::vector<double> MyMatrix::getPadeCoefficients(unsigned int m)
// GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
//    C = GETPADECOEFFICIENTS returns coefficients of numerator
//    of [m/m] Pade approximant, where m = 3,5,7,9,13.
{
  double * c = NULL;
  switch (m)
  {
  case 3:
  {
    double d[] =
    { 120.0, 60.0, 12.0, 1.0 };
    c = d;
    break;
  }
  case 5:
  {
    double d[] =
    { 30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0 };
    c = d;
    break;
  }
  case 7:
  {
    double d[] =
    { 17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0 };
    c = d;
    break;
  }
  case 9:
  {
    double d[] =
    { 17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0,
        2162160.0, 110880.0, 3960.0, 90.0, 1.0 };
    c = d;
    break;
  }
  case 13:
  {
    double d[] =
    { 64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
        1187353796428800.0, 129060195264000.0, 10559470521600.0, 670442572800.0,
        33522128640.0, 1323241920.0, 40840800.0, 960960.0, 16380.0, 182.0, 1.0 };
    c = d;
    break;
  }
  }
  std::vector<double> cvec(c, c + m + 1);
  //delete(c); ??
  return cvec;
}

void MyMatrix::expmchk(const vnl_matrix<double>& A,
    std::vector<unsigned int> &m_vals, std::vector<double> &theta)
//EXPMCHK Check the class of input A and
//    initialize m_vals and theta accordingly.
{
//        classA = class(A);
//        switch classA
//            case 'double'
  unsigned int i[] =
  { 3, 5, 7, 9, 13 };
  m_vals.assign(i, i + 5);

  // theta_m for m=1:13.
  double d[] =
  {
  //3.650024139523051e-008
  //5.317232856892575e-004
      1.495585217958292e-002,// m_vals = 3
      //8.536352760102745e-002
      2.539398330063230e-001,// m_vals = 5
      //5.414660951208968e-001
      9.504178996162932e-001,// m_vals = 7
      //1.473163964234804e+000
      2.097847961257068e+000,// m_vals = 9
      //2.811644121620263e+000
      //3.602330066265032e+000
      //4.458935413036850e+000
      5.371920351148152e+000 };           // m_vals = 13
  theta.assign(d, d + 5);

//             case 'single'
//                 m_vals = [3 5 7];
//                 % theta_m for m=1:7.
//                 theta = [%8.457278879935396e-004
//                          %8.093024012430565e-002
//                           4.258730016922831e-001  % m_vals = 3
//                          %1.049003250386875e+000
//                           1.880152677804762e+000  % m_vals = 5
//                          %2.854332750593825e+000
//                           3.925724783138660e+000];% m_vals = 7
//             otherwise
//                 error('MATLAB:expm:inputType','Input must be single or double.')
//         end

}

vnl_matrix<double> MyMatrix::MatrixExp(const vnl_matrix<double>& m)
//EXPM   Matrix exponential.
//  EXPM(X) is the matrix exponential of X.  EXPM is computed using
//   a scaling and squaring algorithm with a Pade approximation.
//
//   Although it is not computed this way, if X has a full set
//   of eigenvectors V with corresponding eigenvalues D, then
//   [V,D] = EIG(X) and EXPM(X) = V*diag(exp(diag(D)))/V.
//
//   EXP(X) computes the exponential of X element-by-element.
//
//   Reference:
//   N. J. Higham, The scaling and squaring method for the matrix
//   exponential revisited. SIAM J. Matrix Anal. Appl.,
//   26(4) (2005), pp. 1179-1193.
{

  assert(m.cols() == m.rows());
  int n = m.cols();

  std::vector<unsigned int> m_vals;
  std::vector<double> theta;
  expmchk(m, m_vals, theta); // Initialization
  double normA = m.operator_one_norm();

  vnl_matrix<double> F(n, n);

  if (normA <= theta.back())
  {
    // no scaling and squaring is required.
    for (unsigned int i = 0; i < m_vals.size(); i++)
      if (normA <= theta[i])
      {
        F = PadeApproximantOfDegree(m, m_vals[i]);
        break;
      }
  }
  else
  {
    int s;
    double t = frexp(normA / theta.back(), &s); //[t s] = log2(normA/theta(end));
    s = s - (t == 0.5); // adjust s if normA/theta(end) is a power of 2.
    double ss = 1.0 / pow(2.0, s);
    vnl_matrix<double> A(ss * m); //scaling
    F = PadeApproximantOfDegree(A, m_vals.back());
    for (int i = 0; i < s; i++)
      F = F * F;  // Squaring
  }

  return F;

}

// =========================================== MATRIX LOG =================================================

vnl_matrix<double> MyMatrix::MatrixLog(const vnl_matrix<double>& A,
    int maxlogiter)
// LOGM   Matrix logarithm.
//   L = LOGM(A) is the principal matrix logarithm of A, the inverse of EXPM(A).
//   L is the unique logarithm for which every eigenvalue has imaginary part
//   lying strictly between -PI and PI.  If A is singular or has any eigenvalues
//   on the negative real axis, then the principal logarithm is undefined,
//   a non-principal logarithm is computed, and a warning message is printed.
{

//cout << endl << "======================================================="<< endl << endl;

  //cout<< " MyMatrix::MatrixLog " << endl;

  assert(A.cols() == A.rows());
  int n = A.cols();

  if (n == 1)
    return vnl_matrix<double>(1, 1, log(A[0][0]));

  // First form complex Schur form (if m not already upper triangular, maybe check?)
  vnl_matrix<vcl_complex< double > > U(n,n);
  vnl_matrix<vcl_complex< double > > T(n,n);
//  if isequal(A,triu(A))
//    T = A; U = eye(n);
//  else
  SchurComplex(A, U, T);
  //cout << " First Schur U: " << U << endl;
  //cout << " First Schur T: " << T << endl;

  // Check eigenvalues
  double eps = 0.000000001;
  for (int i = 0; i < n; i++)
    if (T[i][i].real() <= 0.0 && abs(T[i][i].imag()) < eps)
    {
      cerr << "MatrixLog Error: " << endl;
      cerr << " Principal matrix logarithm is not defined for A with " << endl;
      cerr << "          nonpositive real eigenvalues. " << endl;
      assert(1==2);
    }

  // Handle special case of diagonal T.
// if isequal(T,tril(T))
//    F = U*diag(feval(fun,diag(T),0,varargin{:}))*U';
//    exitflag = 0; output = struct('terms',ones(n,1),'ind',{1:n});
//    return
// end
// 

  // Determine reordering of Schur form into block form.
  double delta = 0.1;
  std::vector<int> ord1(blocking(T, delta));
  //cout << " ord1: "; Print(ord1);
  std::vector<int> ord;
  std::vector<std::vector<int> > ind;
  swapping(ord1, ord, ind);     // Gives the blocking in ord and ind.
  //cout << " ord: "; Print(ord);
  //cout << " ind: "<< endl;
  //for (unsigned int iii = 0;iii<ind.size(); iii++)
  //{
  //  cout << " " << iii << " " ; Print (ind[iii]);
  //}  
  // Since ORDSCHUR puts highest index top left.
  int mord = *std::max_element(ord.begin(), ord.end());
  //cout << " mord : " << mord << endl;       
  for (unsigned int i = 0; i < ord.size(); i++)
    ord[i] = mord - ord[i];
  //cout << " ord new: " ; Print(ord);
  OrdSchurComplex(U, T, ord, U, T);

  // Calculate F(T)
  vnl_matrix<vcl_complex< double > > F(n,n,0.0);
  int m = ind.size();
  for (int col = 0; col < m; col++)
  {
    //cout << "==================col: " << col << endl;
    std::vector<int> jj = ind[col];
    //cout << " jj: "; Print(jj);
//    if prnt == 2 && max(j) > min(j)
//       fprintf('Evaluating function of block (%g:%g)\n', min(j), max(j))
//    end
//    if isequal(fun,@fun_log)

    //cout << " submatrix t(jj,jj): "<< endl << getSubMatrix(T,jj,jj) << endl;
    //cout << " mlog_isst: " << endl << MatrixLog_isst(getSubMatrix(T,jj,jj),maxlogiter) << endl;
    setSubMatrix(F, jj, jj,
        MatrixLog_isst(getSubMatrix(T, jj, jj), maxlogiter));
    //cout << " F " << F << endl;
//    else
//       [F(j,j), terms(col)] = funm_atom(T(j,j),fun,tol,maxterms,prnt>1,...
//                                        varargin{:});
//    end
// 
    for (int row = col - 1; row >= 0; row--)
    {
      //cout << "=======row : " << row << endl;
      std::vector<int> ii = ind[row];
      //cout << " ii: "; Print(ii);
      if (ii.size() == 1 && jj.size() == 1)
      {
        // Scalar case.
        int i = ii[0];
        int j = jj[0];
        //cout << " scalar case ( " << i << " , " << j << " )" << endl;
        //k = i+1:j-1;
//          if (i+1 > j-1)  // empty k below
//          {
//            //cout << " !!! empty k  !! "  << endl;
//            F[i][j] = (T[i][j]*(F[i][i] - F[j][j]))/(T(i,i)-T(j,j));
//          }
//          else
//          {
//            for (int k = i+1; k<= j-1; k++)
//            {
//              vcl_complex < double > temp = T[i][j]*(F[i][i] - F[j][j]) + F[i][k]*T[k][j] - T[i][k]*F[k][j];
//              F[i][j] = temp/(T(i,i)-T(j,j));
//            }
//          }

        vcl_complex< double > temp1(0);
        vcl_complex< double > temp2(0);
        for (int k = i + 1; k <= j - 1; k++)
        {
          temp1 += F[i][k] * T[k][j];
          temp2 += T[i][k] * F[k][j];
        }
        F[i][j] = (T[i][j] * (F[i][i] - F[j][j]) + temp1 - temp2)
            / (T(i, i) - T(j, j));

      }
      else
      {
        //cout << " matrix case  " << endl;

        // concatenate indices: (matlab) k = cat(2,ind{row+1:col-1});
        std::vector<int> kk;
        for (int k = row + 1; k < col; k++)
          kk.insert(kk.end(), ind[k].begin(), ind[k].end());
        //cout << " kk : " ; Print(kk);

        //rhs = F(i,i)*T(i,j) - T(i,j)*F(j,j) + F(i,k)*T(k,j) - T(i,k)*F(k,j);
        vnl_matrix<vcl_complex< double > > rhs(getSubMatrix(F,ii,ii) *getSubMatrix(T,ii,jj));
        rhs -= getSubMatrix(T, ii, jj) * getSubMatrix(F, jj, jj);
        //cout << " rhs1 : " << rhs << endl;
        if (kk.size() > 0)
        {
          rhs += getSubMatrix(F, ii, kk) * getSubMatrix(T, kk, jj);
          rhs -= getSubMatrix(T, ii, kk) * getSubMatrix(F, kk, jj);
        }
        //cout << " rhs2 : " << rhs << endl;

        //F(i,j) = sylv_tri(T(i,i),-T(j,j),rhs);
        setSubMatrix(F, ii, jj,
            sylv_tri(getSubMatrix(T, ii, ii), -getSubMatrix(T, jj, jj), rhs));
        //cout << " F " << F << endl;

      }
    } // rows in upper triangle
  }

//      cout << endl;
//      cout.precision(16);
//      cout << " T   = ";PrintComplex(T);
//      cout << " F   = ";PrintComplex(F);
//      cout.precision(6);

  F = U * F * U.conjugate_transpose();

  //cout << " log(A) ? " << F << endl;

  double eps_matlab = pow(2.0, -52);
  double imagone = vnl_imag(F).operator_one_norm();
  if (imagone > 10 * n * eps_matlab * F.operator_one_norm())
  {
    cerr << " MatrixLog Error: " << endl;
    cerr << "  Result too imaginary to ignore! ( " << imagone << " )" << endl;
    cerr.precision(16);
    cerr << " A = " << A << endl;
    cerr << " T = " << T << endl;
    cerr << " U = " << U << endl;
    cerr << " fro( real(U * T * U^*) - A) = "
        << (vnl_real(U * T * U.conjugate_transpose()) - A).frobenius_norm()
        << endl << endl;
    cerr << " log(A) = " << F << endl;
    cerr << " fro( exp(real(log(A))) - A)  = "
        << (MatrixExp(vnl_real(F)) - A).frobenius_norm() << endl << endl;
    exit(1);
  }

  return vnl_real(F);

}

void MyMatrix::OrdSchurComplexLogical(const vnl_matrix<vcl_complex< double > > & U ,
const vnl_matrix < vcl_complex < double > > & T,
const std::vector < int >& select,
vnl_matrix < vcl_complex < double > > & US,
vnl_matrix < vcl_complex < double > > & TS)
// wrapping ZTRSEN from LAPACk:
//  reorders the Schur factorization of a complex matrix
//  A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in
//  the leading positions on the diagonal of the upper triangular matrix
//  T, and the leading columns of Q form an orthonormal basis of the
//  corresponding right invariant subspace.
{

  //cout << " MyMatrix::OrdSchurComplex ( logical select )  NOT TESTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

  assert(T.rows() == T.cols());
  assert(U.rows() == U.cols());
  assert(U.rows() == T.rows());

  //cout << " T before : " << T << endl;

  int * selecti = new int[select.size()];
  for (unsigned int i = 0;i<select.size();i++)
  selecti[i] = select[i];

  long n = T.rows();
  TS = T.transpose();// transpose for fortran
  US = U.transpose();// transpose for fortran
  vnl_vector<vcl_complex<double > > W(n);
  long m;
  long lwork = n;
  vnl_vector < vcl_complex < double > > work(lwork);
  long info = 100;
  ztrsen_ (
  (char*)"N",     // JOB (no condition numbers)
  (char*)"V",// COMPQ (update the schur vectors also)
  selecti,// SELECT (logical array)
  &n,// N
  TS.data_block(),// T (input/output)
  &n,// LDT
  US.data_block(),// Q (input/output)
  &n,// LDQ
  W.data_block(),// W (output) reordered eigenvalues
  &m,// m (output) dimension of specified invariant subspace
  0,// S
  0,// SEP
  work.data_block(),// WORK
  &lwork,// LWORK 
  &info);
  delete[] selecti;
  TS.inplace_transpose();// switch back ..
  US.inplace_transpose();// from fortran ordering
  //cout << " info: " << info << endl;
  //cout << " T after " << TS << endl;

  if (info != 0)
  {
    // These return codes are taken from ztrsen.f:
    //*          = 0:  successful exit
    //*          < 0:  if INFO = -i, the i-th argument had an illegal value.

    cerr << __FILE__ ": info = " << info << ", something went wrong:\n";
    if (info < 0)
    {
      cerr << __FILE__ ": (internal error) the " << (-info) << "th argument had an illegal value\n";
    }
    else
    {
      cerr << __FILE__ ": unknown error\n";
    }

    assert(info==0);
  }
}

void MyMatrix::OrdSchurComplex(const vnl_matrix<vcl_complex< double > > & U,
const vnl_matrix < vcl_complex < double > > & T,
const std::vector < int > & clusters,
vnl_matrix < vcl_complex < double > > & US,
vnl_matrix < vcl_complex < double > > & TS)
// Reorders the Schur factorization 
//    X = U*T*U' of a matrix X so that a selected cluster of eigenvalues 
//    appears in the leading (upper left) diagonal of the 
//    triangular Schur matrix T, and the corresponding invariant 
//    subspace is spanned by the leading columns of U. 
// Can reorder multiple clusters at once.  Given a vector 
//    CLUSTERS of cluster indices, commensurate with E = EIG(T), and such
//    that all eigenvalues with the same CLUSTERS value form one cluster,
//    [US,TS] = ORDSCHUR(U,T,CLUSTERS) will sort the specified clusters in
//    descending order along the diagonal of TS, the cluster with highest 
//    index appearing in the upper left corner.
{

  assert(U.cols() == U.rows());
  assert(T.cols() == T.rows());
  assert(U.cols() == T.cols());
  assert(clusters.size() == T.cols());

  //cout << " U " << endl << U << endl;
  //cout << " T " << endl << T << endl;

  int n = T.cols();

  US = U;
  TS = T;

  int maxi = *std::max_element( clusters.begin(), clusters.end());
  std::vector < int > myclusters(clusters);
  //cout << " mycluster before "; Print(myclusters);
  //bool reorder = false;
  for (int i = 0; i<=maxi; i++)
  {

    // create logical 'select' vector of all elelemts where 
    //  mycluster == i 
    int foundi = 0;
    std::vector < int > select(n);
    for (int j = 0;j<n;j++)
    {
      if (myclusters[j] == i) // found one
      {
        select[j] = 1;
        foundi++;
      }
      else select[j] = 0; // don't reorder this one      
    }

    //cout << "select " ; Print(select);
    //cout << "foundi " << foundi << endl;

    if (foundi==0) continue;// if there was a gap in clusters idx, continue with next

    int startvals = 0;
    while ( startvals < n && select[startvals]) startvals++;
    //cout << "startvals " << startvals << endl;
    if (foundi==startvals) continue;// all selected are allready the first values

    //cout << " doing schur reorder... " << endl;
    // here we have found at least one value to reorder
    //reorder = true;
    OrdSchurComplexLogical(US,TS,select,US,TS);

    // also reorder myclusters accordingly:
    std::vector < int > my2(myclusters);
    int count = 0;
    for (int j = 0; j<n; j++ )
    {
      // setting i cluster to the front:
      if (j<foundi) myclusters[j] = i;
      else
      {
        // skipping i clusters
        while (count < n && my2[count] == i) count++;
        assert (my2[count] != i );
        // copy rest
        myclusters[j] = my2[count];
        count++;
      }
    }

    //cout << "mycluster after: " ; Print(myclusters);

  }

  // returning parameters US and TS
  //cout << " US " << endl << US << endl;
  //cout << " TS " << endl << TS << endl;
  //if (reorder)   exit(1);
}

vnl_matrix<vcl_complex< double > > MyMatrix::getSubMatrix(const vnl_matrix < vcl_complex < double > > & A,
    const std::vector < int > & rows,
    const std::vector < int > & cols)
// returns the submatrix of A consising of rows x cols (where these are indices)
{

  vnl_matrix < vcl_complex < double > > M (rows.size(),cols.size());

  for (unsigned int r = 0; r<rows.size(); r++)
  for (unsigned int c = 0; c<cols.size(); c++)
  M[r][c] = A[rows[r]][cols[c]];

  return M;

}

void MyMatrix::setSubMatrix(vnl_matrix<vcl_complex< double > > & A,
const std::vector < int > & rows,
const std::vector < int > & cols,
const vnl_matrix < vcl_complex < double > > & B)
// sets the submatrix of A specified by rows and cols to B
{
  assert(B.rows() == rows.size());
  assert(B.cols() == cols.size());

  for (unsigned int r = 0; r<rows.size(); r++)
  for (unsigned int c = 0; c<cols.size(); c++)
  {
    // also a should be large enough
    assert((int)A.rows()>rows[r]);
    assert((int)A.cols()>cols[c]);
    A[rows[r]][cols[c]] = B[r][c];
  }

}

vnl_matrix<vcl_complex<double > > MyMatrix::MatrixLog_isst(const vnl_matrix < vcl_complex < double > >& A, int maxlogiter)
// Log of triangular matrix by inverse scaling and squaring method.
//   X = LOGM_ISST(m, MAXLOGITER) computes the (principal) logarithm of an
//   upper triangular matrix m, for a matrix with no nonpositive real
//   eigenvalues, using the inverse scaling and squaring method with Pade
//   approximation.  At most MAXLOGITER square roots are computed.
//   [X, ITER] = LOGM_ISST(m, MAXLOGITER, PRNT) returns the number
//   ITER of square roots computed and prints this information if
//   PRNT is nonzero.
//
//   References:
//   S. H. Cheng, N. J. Higham, C. S. Kenney, and A. J. Laub, Approximating the
//      logarithm of a matrix to specified accuracy, SIAM J. Matrix Anal. Appl.,
//      22(4):1112-1125, 2001.
//   N. J. Higham, Evaluating Pade approximants of the matrix logarithm,
//      SIAM J. Matrix Anal. Appl., 22(4):1126-1135, 2001.
{
  int n = A.cols();

  if (n == 1) return vnl_matrix < vcl_complex < double > > (1,1,log(A[0][0]));

  vnl_matrix < vcl_complex < double > > T(A),R(A);

  double phi;
  bool prnt = false; // true;
  int iter, j,i,k;
//  vcl_complex < double > s;
  vnl_diag_matrix < vcl_complex < double > > I(n,1);//  I.set_identity();

  for (iter = 0;iter <maxlogiter; iter++)
  {

    phi = (T-I).frobenius_norm();

    if (phi <= 0.25)
    {
      if (prnt) cout <<"MatrixLog_isst computed "<< iter << " square roots. \n";
      break;
    }
    if ( iter == maxlogiter-1 )
    {
      cerr << "MatrixLog_isst:tooManyIterations Too many square roots." << endl;
      assert(iter < maxlogiter-1);
      exit (1);
    }

    // Compute upper triangular square root R of T, a column at a time.
    for (j= 0; j< n; j++)
    {
      R[j][j] = sqrt(T[j][j]);
      for (i=j-1; i >=0; i--)
      {
//           if (i+1>=j)  // empty k below => zero s
//             R[i][j] = (T[i][j])/(R[i][i] + R[j][j]);
//           else
//           for (k=i+1; k<j; k++)
//           {
//             //cout << i << " " << j << " " << k << endl;
//             s = R[i][k]*R[k][j];
//             R[i][j] = (T[i][j] - s)/(R[i][i] + R[j][j]);
//           }

        vcl_complex < double > s(0);
        for (k=i+1; k<j; k++)
        {
          //cout << i << " " << j << " " << k << endl;
          s += R[i][k]*R[k][j];
        }
        R[i][j] = (T[i][j] - s)/(R[i][i] + R[j][j]);

      }
    }

    T = R;
  }
  //cout << " T in logisst: " << endl << T << endl;
  //cout << " iter : " << iter << endl;
  //cout << " log pf : " << MatrixLog_pf(T-I,8) << endl;
  return vcl_complex < double > (pow(2.0,iter),0)*MatrixLog_pf(T-I,8);

}

vnl_matrix<vcl_complex<double > > MyMatrix::MatrixLog_pf(const vnl_matrix < vcl_complex <double > >& A, unsigned int m)
//LOGM_PF   Pade approximation to matrix log by partial fraction expansion.
//          Y = LOGM_PF(A,m) approximates LOG(EYE(SIZE(A))+A) by an [m/m]
//          Pade approximation.
{
  int n = A.cols();
  vnl_vector < double > nodes(n);
  vnl_vector < double > wts(n);
  gauss_legendre(m,nodes,wts);
  //cout << " nodes   " << nodes << endl;
  //cout << " weights "  << wts << endl;

  // Convert from [-1,1] to [0,1].
  nodes = (nodes + 1)/2.0;
  wts = wts/2.0;

  vnl_matrix < vcl_complex < double > > S(n,n,0.0);
  vnl_matrix < vcl_complex < double > > id(n,n); id.set_identity();

  vnl_matrix < vcl_complex < double > > At(A.transpose());
  for (unsigned int j=0;j<m; j++ )
  {
    S = S + vcl_complex<double>(wts[j],0.0) * (vnl_svd< vcl_complex < double > >(id+vcl_complex<double>(nodes[j],0.0)*At).solve(At));
    //S = S + wts[j]*(A/(id + nodes[j]*A)); //!!! rightsolv requires transpose
  }

  return S.transpose();
}

void MyMatrix::gauss_legendre(int n, vnl_vector<double>& x,
    vnl_vector<double> & w)
//GAUSS_LEGENDRE  Nodes and weights for Gauss-Legendre quadrature.
//                [X,W] = GAUSS_LEGENDRE(N) computes the nodes X and weights W
//                for n-point Gauss-Legendre quadrature.
//
// Reference:
// G. H. Golub and J. H. Welsch, Calculation of Gauss quadrature rules
// Math. Comp., 23(106):221-230, 1969.
{
  vnl_matrix<double> M(n, n, 0.0);
  double d;
  for (int i = 1; i < n; i++)
  {
    d = i / sqrt(4.0 * i * i - 1);
    M[i - 1][i] = d;
    M[i][i - 1] = d;
  }
  //cout << " offidag : " << M << endl;
  vnl_symmetric_eigensystem<double> eigs(M);

  //cout << " V : " << eigs.V << endl;
  //cout << " D : " << eigs.D << endl;

  x = eigs.D.diagonal();
  w = eigs.V.get_row(0);
  for (int i = 0; i < n; i++)
    w[i] = 2.0 * w[i] * w[i];

//  i = 1:n-1;
//  v = i./sqrt((2*i).^2-1);
//  [V,D] = eig( diag(v,-1)+diag(v,1) );
//  x = diag(D);
//  w = 2*(V(1,:)'.^2);
}

std::vector<int> MyMatrix::blocking(const vnl_matrix<vcl_complex< double > > &A,double delta)
//BLOCKING  Produce blocking pattern for block Parlett recurrence in FUNM.
//   M = BLOCKING(A, DELTA) accepts an upper triangular matrix
//   A and produces a blocking pattern, specified by the vector M,
//   for the block Parlett recurrence.
//   M(i) is the index of the block into which A(i,i) should be placed,
//   for i=1:LENGTH(A).
//   DELTA is a gap parameter (default 0.1) used to determine the blocking.
//
//   For A coming from a real matrix it should be posible to take
//   advantage of the symmetry about the real axis.  This code does not.
{
  int n = A.cols();

  vnl_vector < vcl_complex < double > > a(n);
  for (int i = 0; i< n; i++) a[i] = A[i][i];

  //cout << " A diag: " << a << endl;

  std::vector < int > m(n,-1);
  int maxM = 0;

  for (int i = 0; i<n;i++)
  {

    if (m[i] == -1)
    {
      m[i] = maxM;    // If a(i) hasn`t been assigned to a set
      maxM = maxM +1;// then make a new set and assign a(i) to it.
    }

    for (int j = i+1; j<n; j++) // run through remaining values
    {
      if (m[i] != m[j] )    // If a(i) and a(j) are not in same set.
      {
        if ( abs(a[i]-a[j]) <= delta )
        {
          // create union of these two sets
          if ( m[j] == -1 )
          m[j] = m[i];// If a(j) hasn`t been assigned to a
                      // set, assign it to the same set as a(i).
          else
          {
            int p = max(m[i],m[j]);
            int q = min(m[i],m[j]);
            for (int k = 0; k<n; k++)
            {
              if (m[k] == p)
              m[k]= q; // If a(j) has been assigned to a set
                       // place all the elements in the set
                       // containing a(j) into the set
                       // containing a(i) (or vice versa).
              else if (m[k] > p)
              m[k]= m[k]-1;
            }
            maxM = maxM - 1;
            // Tidying up. As we have deleted set
            // p we reduce the index of the sets
            // > p by 1.
          }
        }
      }
    }

  }

  return m;
}

void MyMatrix::swapping(const std::vector<int> &m, std::vector<int> &mm,
    std::vector<std::vector<int> > &ind)
//SWAPPING  Choose confluent permutation ordered by average index.
//   [MM,IND] = SWAPPING(M) takes a vector M containing the integers
//   1:K (some repeated if K < LENGTH(M)), where M(J) is the index of
//   the block into which the element T(J,J) of a Schur form T
//   should be placed.
//   It constructs a vector MM (a permutation of M) such that T(J,J)
//   will be located in the MM(J)'th block counting from the (1,1) position.
//   The algorithm used is to order the blocks by ascending
//   average index in M, which is a heuristic for minimizing the number
//   of swaps required to achieve this confluent permutation.
//   The cell array vector IND defines the resulting block form:
//   IND{i} contains the indices of the i'th block in the permuted form.
{
  int mmax = *std::max_element(m.begin(), m.end()) + 1;
  mm.clear();
  mm.resize(m.size());
  std::vector<std::pair<double, int> > g(mmax);
  std::vector<int> h(mmax);

  int lengthp, sump;
  for (int i = 0; i < mmax; i++)
  {
    lengthp = 0;
    sump = 0;
    for (unsigned int j = 0; j < m.size(); j++)
    {
      if (m[j] == i)
      {
        lengthp++;
        sump += j; // i+1 as in matlab will not change the order
      }
    }

    h[i] = lengthp;
    g[i].first = ((double) sump) / h[i];
    g[i].second = i;

    //p = find(m==i);
    //h(i) = length(p);
    //g(i) = sum(p)/h(i);
  }

//   cout << " g(pre): " ; Print(g);
// [x,y] = sort(g);
  std::sort(g.begin(), g.end());
//   cout << " g(post): " ; Print(g);

// h = [0 cumsum(h(y))];
  h = cumsum0(h, g);

  //cout << " h " ; Print(h);

  ind.clear();
  ind.resize(mmax);
  for (int i = 0; i < mmax; i++)
  {
    for (unsigned int j = 0; j < m.size(); j++)
      if (m[j] == g[i].second)
        mm[j] = i;

    ind[i].resize(h[i + 1] - h[i]);
    int count = 0;
    for (int j = h[i]; j < h[i + 1]; j++)
    {
      ind[i][count] = j;
      count++;
    }
  }
  // returning paramters mm and ind
}

std::vector<int> MyMatrix::cumsum0(const std::vector<int> &v,
    const std::vector<std::pair<double, int> > & w)
// helper for swapping
// computes the cumulative sum of v with the following modifications:
// 1. a zero is prepended
// 2. the second parameter of vector w is used for reordering the array first
// so if y is the array of all w.second this is equiv to matlab:
// cumsum0 = [0 cumsum(v(y))];
// note that w.first values are ignored, the pair is passed to make calling easier
{
  int n = (int) v.size();
  std::vector<int> cs(n + 1);
  cs[0] = 0;

  int sum = 0;
  for (int i = 0; i < n; i++)
  {
    sum = sum + v[w[i].second];
    cs[i + 1] = sum;
  }
  return cs;
}

vnl_matrix<vcl_complex< double > >
MyMatrix::sylv_tri(const vnl_matrix < vcl_complex < double > > & T,
    const vnl_matrix < vcl_complex < double > > & U,
    const vnl_matrix < vcl_complex < double > > & B)
//SYLV_TRI    Solve triangular Sylvester equation.
//   X = SYLV_TRI(T,U,B) solves the Sylvester equation
//   T*X + X*U = B, where T and U are square upper triangular matrices.
    {
      assert (T.cols() == T.rows());
      assert (U.cols() == U.rows());
      int m = T.cols();
      int n = U.cols();

      //cout << " sylv_tri " << endl;
      //cout << " T " << T << endl;
      //cout << " U " << U << endl;
      //cout << " B " << B << endl;

      vnl_matrix < vcl_complex < double > > M(m,n);

      // Forward substitution
      for (int i = 0;i<n;i++)
      {
        //cout << " col " << i << endl;

        vnl_diag_matrix < vcl_complex < double > > Uii(m,U[i][i]);
        vnl_matrix < vcl_complex < double > > A(T + Uii);
        //cout << " A " << A << endl;
        vnl_vector < vcl_complex < double > > y(B.get_column(i));
        if (i>0)
        y -= M.get_n_columns(0,i-1) * U.get_n_rows(0,i-1).get_column(i);
        //cout << " y " << Y << endl;
        M.set_column(i,vnl_svd<vcl_complex < double > >(A).solve(y));
        // matlab:
        // X(:,i) = (T + U[i][i]*Id) \ (B(:,i) - X(:,1:i-1)*U(1:i-1,i));
      }

      return M;
    }

// ==================================== MATRIX POWER ============================================================

/**
 Uses matrix log and exp to compute the matrix power
 @param[in]     A The input matrix
 @param[in]     d the exponent
 @return M to the power of d.
 */
vnl_matrix<double> MyMatrix::MatrixPow(const vnl_matrix<double>& A, double d)
// use log and exp to compute power
{
  return MatrixExp(d * MatrixLog(A));
}

// ==================================== MATRIX SQUAREROOT ============================================================

/**
 Computes the squareroot of a matrix \f$M = \sqrt{M} * \sqrt{M}\f$ 
 using complex SCHUR decomposition (same as matlab,
 although they seem to use slightly different LAPACK routines for the SCHUR decomposition)
 Note: there exists also a completely REAL algorithm (not implemented here), see
 Computing Real Square Roots of a Real Matrix
 by Nicholas J. Higham
 
 @param[in]     A The input matrix needs to be a square matrix.
 @return Square root of A.
 */
vnl_matrix<double> MyMatrix::MatrixSqrt(const vnl_matrix<double>& A)
{

//

//  vnl_matrix < double > M(4,4); M.set_identity();
////   M[0][0] = -149; M[0][1] = -50 ; M[0][2] = -154;
////   M[1][0] =  537; M[1][1] = 180 ; M[1][2] =  546;
////   M[2][0] =  -27; M[2][1] =  -9 ; M[2][2] =  -25;
//  M[0][3] = 0.09; M[1][3] = 0.19; M[2][3] = 0.47;

//   //cout << " matrix exponent: " << endl << MatrixExp(M) << endl;
//   cout << " matrix Log: " << endl << MatrixLog(M) << endl;
//   cout << " matrix Log: " << endl << MatrixPow(M,0.5) << endl;
// 
// exit(0);

//   cout.precision(18);
// //  cout << "A MATRIX : " <<  A << endl;
//   cout << " A = [ " << A.get_row(0) << " ; " << endl;
//   cout << "       " << A.get_row(1) << " ; " << endl;
//   cout << "       " << A.get_row(2) << " ; " << endl;
//   cout << "       " << A.get_row(3) << " ] " << endl;
//   
//   
//   cout.precision(5);

  //return MatrixPow(A,0.5); // for testing log and exp

  assert(A.cols() == A.rows());
  int n = A.cols();

  vnl_matrix<vcl_complex< double > > U(n,n);
  vnl_matrix<vcl_complex< double > > T(n,n);
  SchurComplex(A, U, T);

  vnl_matrix<vcl_complex< double > > Asqrt(n,n);

  if (isDiag(vnl_real(T)) && isDiag(vnl_imag(T)))
  {
    //cout << " isDiagonal ! " << endl;
    vnl_diag_matrix<vcl_complex< double > > R(n);
    for (int i = 0;i<n;i++)
    R[i] = sqrt(T[i][i]);  // Square root always exists. 

        Asqrt = U * R * U.conjugate_transpose();
      }
      else
      {

        // Compute upper triangular square root R of T, a column at a time.
        //cout << " upper triangular ! " << endl;

        vnl_matrix < vcl_complex < double > > R(n,n,0.0);
        // vcl_complex < double > s;
        for (int j= 0; j< n; j++)// column
        {
          R[j][j] = sqrt(T[j][j]);
          for (int i=j-1; i >=0; i--)
          {
            vcl_complex < double > s(0);
            for (int k=i+1; k<j; k++)
            {
              s += R[i][k]*R[k][j];
            }
            R[i][j] = (T[i][j] - s)/(R[i][i] + R[j][j]);
          }
        }

//     cout << endl;
//     cout.precision(16);
//     cout << " T = ";PrintComplex(T);
//     cout << " R = ";PrintComplex(R);
//     cout << " R*R = ";PrintComplex(R*R);
//     cout << " R*R - T = ";PrintComplex(R*R-T);
//     cout << " fro(R*R-T) = " << (R*R-T).frobenius_norm() << endl << endl ;
//     cout.precision(6);

        Asqrt = U * R * U.conjugate_transpose();

        //cout << " sqrt " << endl << Asqrt << endl;
      }

  bool test = true;
  if (test)
  {
    // test if imaginary values too large:
    double eps_matlab = pow(2.0, -52);
    double imagone = vnl_imag(Asqrt).operator_one_norm();
    if (imagone > 10 * n * eps_matlab * Asqrt.operator_one_norm())
    //if ( imagone > 0.000000001 )
    {
      cerr << "MyMatrix::MatrixSqrt Internal Error: " << endl;
      cerr << " Result too imaginary to ignore! ( " << imagone << " )" << endl;
      cerr << " Debug Info: " << endl;
      vnl_matlab_print(vcl_cerr,A,"A",vnl_matlab_print_format_long);
      cerr << endl;
      vnl_matlab_print(vcl_cerr,T,"T",vnl_matlab_print_format_long);
      cerr << endl;
      vnl_matlab_print(vcl_cerr,U,"U",vnl_matlab_print_format_long);
      cerr << endl;
      cerr << " fro( real(U * T * U^*) - A) = "
          << (vnl_real(U * T * U.conjugate_transpose()) - A).frobenius_norm()
          << endl << endl;
      vnl_matlab_print(vcl_cerr,Asqrt,"Asqrt",vnl_matlab_print_format_long);
      cerr << endl;
      cerr << " fro( real(sqrt(A)) ^2 - A)  = "
          << (vnl_real(Asqrt) * vnl_real(Asqrt) - A).frobenius_norm() << endl
          << endl;
      exit(1);
    }
    // double fnorm = vnl_imag(T).frobenius_norm();
    ////cout << " fnorm (img) = " << fnorm << endl;
    //if (fnorm > eps)
    //{
    //   cout << " Warning complex result?: " << fnorm << endl << vnl_imag(T) << endl;
    //     assert(1==2);
    //  }

    // test if sqrt^2==A
    double eps = 1e-10;
    vnl_matrix<double> msqrt = vnl_real(Asqrt);
    vnl_matrix<double> ms2 = msqrt * msqrt;
    ms2 -= A;
    //double sum = ms2.absolute_value_max();
    double fnorm = ms2.frobenius_norm();
    if (fnorm > eps)
    {
      cerr << "MyMatrix::MatrixSqrt Internal Error: " << endl;
      cerr << "  Result not close enough to square root!" << endl;
      cerr << "  Frobenius norm of difference: " << fnorm << endl;
      vnl_matlab_print(vcl_cerr,A,"A",vnl_matlab_print_format_long);
      cerr << endl;
      vnl_matlab_print(vcl_cerr,msqrt,"Asqrt",vnl_matlab_print_format_long);
      cerr << endl;
      vnl_matlab_print(vcl_cerr,ms2,"Diff",vnl_matlab_print_format_long);
      cerr << endl;
      exit(1);
    }
  }

  return vnl_real(Asqrt);

}

/**
 Checks if a matrix is diagonal. It returns true if
 the sum of squared differences (SSD) of the off diagonal elements
 is below a given epsilon.
 
 @param[in]     A The input matrix needs to be square.
 @param[in]     eps the threshold for the SSD.
 @return True if SSD in off diagonal elements is below eps.
 */
bool MyMatrix::isDiag(const vnl_matrix<double>& A, double eps)
{

  int i, j;
  int n = A.rows();

  double sse = 0.0;

  for (i = 1; i < n; i++)
    for (j = i - 1; j >= 0; j--)
    {
      sse += A[i][j] * A[i][j];
      sse += A[j][i] * A[j][i];
    }

  return sse < eps;

}

/**
 Computes the squareroot of a matrix \f$M = \sqrt{M} * \sqrt{M}\f$ 
 using its eigen-decomposition (here via svd).
 This works only if the matrix is diagonalizable.
 For now separate the translation, else we easily get defective m (not diagonalizable)
 (where we cannot do the eigendecomposition trick).
 
 @param[in]     m The input matrix needs to be a 4x4 affine matrix.
 @return Square root of m.
 */
vnl_matrix<double> MyMatrix::MatrixSqrtEigs(const vnl_matrix<double>& m)
{
  assert(m.rows() == 4 && m.cols() == 4);

  vnl_matrix<double> R(3, 3);  // = m.extract(3,3,0,0);
  for (int rr = 0; rr < 3; rr++)
    for (int cc = 0; cc < 3; cc++)
    {
      R[rr][cc] = m[rr][cc];
    }

//cout << endl << endl;
//cout << " M: " << endl << m << endl;
//cout << " R: " << endl << R << endl;

  vnl_matrix<double> rsqrt(3, 3, 0.0);
  vnl_complex_eigensystem esys(R, rsqrt); //complex part is zero

//cout << " V' " << endl << esys.R << endl;
//cout << " D " << endl << esys.W << endl;

  vnl_diag_matrix<vcl_complex< double > > Wsqrt(3);
  for (unsigned int i = 0; i < 3; i++)
  {
    Wsqrt[i] = sqrt(esys.W[i]);
  }
//cout << " Wsqrt " << endl << Wsqrt << endl;

  vnl_matrix<vcl_complex< double > > Rcomp (3,3);

//  esys.R.inplace_transpose(); // store evec in columns
//  vnl_matrix_fixed < vcl_complex < double >, 4, 4 > Rt = esys.R;
//  Rcomp = (esys.R * Wsqrt) * vnl_inverse(Rt);
  vnl_qr<vcl_complex< double > > QR(esys.R);
  Rcomp = QR.solve(Wsqrt * esys.R);
  Rcomp.inplace_transpose();

//cout << " Rcomp " << endl << Rcomp << endl;

  rsqrt = vnl_real(Rcomp);

  // compute new T
  // Rh1 = R + I
  vnl_matrix_fixed<double, 3, 3> Rh1(rsqrt);
  Rh1[0][0] += 1;
  Rh1[1][1] += 1;
  Rh1[2][2] += 1;

  vnl_vector_fixed<double, 3> T;
  T[0] = m[0][3];
  T[1] = m[1][3];
  T[2] = m[2][3];

  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  vnl_vector_fixed<double, 3> Th = vnl_inverse(Rh1) * T; //vnl_svd<double>(Rh1).solve(T);

  // put everything together:
  vnl_matrix<double> msqrt(4, 4);
  msqrt[0][3] = Th[0];
  msqrt[1][3] = Th[1];
  msqrt[2][3] = Th[2];
  msqrt[3][0] = 0.0;
  msqrt[3][1] = 0.0;
  msqrt[3][2] = 0.0;
  msqrt[3][3] = 1.0;
  for (int c = 0; c < 3; c++)
    for (int r = 0; r < 3; r++)
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
      cout << " Warning complex result?: " << fnorm << endl << vnl_imag(Rcomp)
          << endl;
      return MatrixSqrtIter(m);
      //assert(1==2);
    }

    vnl_matrix<double> ms2 = msqrt * msqrt;
    ms2 -= m;
    double sum = ms2.absolute_value_max();
    //cout << " max = " << sum << endl;
    if (sum > eps)
    {
      cout << " Warning : " << sum << endl;
      cout << " sqrt(M): " << endl << msqrt << endl;
      return MatrixSqrtIter(m);
    }
  }

  return msqrt;
}

/**
 Computes the sqrt of a matrix \f$M = \sqrt{M} * \sqrt{M}\f$
 using Schur decomposition on the 3x3 linear map after splitting off the translation.
 Note, [smith:01] say that sqrt(affine) is not necessarily affine
 
 @param[in]     m The input matrix needs to be a 4x4 affine matrix.
 @return Square root of m.
 */
vnl_matrix<double> MyMatrix::MatrixSqrtAffine(const vnl_matrix<double>& m)
{
  assert(m.rows() == 4 && m.cols() == 4);

  vnl_matrix_fixed<double, 3, 3> R; // = m.extract(3,3,0,0);
  for (int rr = 0; rr < 3; rr++)
    for (int cc = 0; cc < 3; cc++)
    {
      R[rr][cc] = m[rr][cc];
    }

  vnl_matrix_fixed<double, 3, 3> sqrtR(MatrixSqrt(R));

  vnl_vector_fixed<double, 3> T;
  T[0] = m[0][3];
  T[1] = m[1][3];
  T[2] = m[2][3];

  // Rh1 = sqrt(R) + I
  vnl_matrix_fixed<double, 3, 3> Rh1(sqrtR);
  Rh1[0][0] += 1;
  Rh1[1][1] += 1;
  Rh1[2][2] += 1;
  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  // can Rh1 always be inverted ??? not sure...
  vnl_vector_fixed<double, 3> Th = vnl_inverse(Rh1) * T; //vnl_svd<double>(Rh1).solve(T);

  // put everything together:
  vnl_matrix<double> msqrt(4, 4);
  msqrt[0][3] = Th[0];
  msqrt[1][3] = Th[1];
  msqrt[2][3] = Th[2];
  msqrt[3][0] = 0.0;
  msqrt[3][1] = 0.0;
  msqrt[3][2] = 0.0;
  msqrt[3][3] = 1.0;
  for (int c = 0; c < 3; c++)
    for (int r = 0; r < 3; r++)
      msqrt[r][c] = sqrtR[r][c];

  return msqrt;
}

/**
 Computes the matrix square root using an iterative algorithm by Denman and Beavers.
 The square root here is defined so that \f$M = \sqrt{M} * \sqrt{M}\f$.
 Note that symmetry cannot be guaranteed (termination may be before reaching optimum at
 machine precision).
 
 @param[in]     m The input matrix needs to be 4x4.
 @return Approximated square root of m.
 */
vnl_matrix<double> MyMatrix::MatrixSqrtIter(const vnl_matrix<double>& m)
//Denman and Beavers square root iteration
{
  assert(m.rows() == 4 && m.cols() == 4);

  vnl_matrix_fixed<double, 3, 3> R; // = m.extract(3,3,0,0);
  for (int rr = 0; rr < 3; rr++)
    for (int cc = 0; cc < 3; cc++)
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
  vnl_matrix_fixed<double, 3, 3> Yn(R);
  vnl_matrix_fixed<double, 3, 3> Zn;
  Zn.set_identity();
  vnl_matrix_fixed<double, 3, 3> Yni, Zni;
  vnl_matrix_fixed<double, 3, 3> Ysq;

  int count = 0;
  while (count < imax && err > eps)
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
  vnl_matrix_fixed<double, 3, 3> Rh1(Yn);
  Rh1[0][0] += 1;
  Rh1[1][1] += 1;
  Rh1[2][2] += 1;

  vnl_vector_fixed<double, 3> T;
  T[0] = m[0][3];
  T[1] = m[1][3];
  T[2] = m[2][3];

  // solve T = Rh1 * Th   <=>   Th = Rh1^-1 * T
  vnl_vector_fixed<double, 3> Th = vnl_inverse(Rh1) * T; //vnl_svd<double>(Rh1).solve(T);

  // put everything together:
  vnl_matrix<double> msqrt(4, 4);
  msqrt[0][3] = Th[0];
  msqrt[1][3] = Th[1];
  msqrt[2][3] = Th[2];
  msqrt[3][0] = 0.0;
  msqrt[3][1] = 0.0;
  msqrt[3][2] = 0.0;
  msqrt[3][3] = 1.0;
  for (int c = 0; c < 3; c++)
    for (int r = 0; r < 3; r++)
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
//          cerr << " sqrt(M): " << endl << msqrt << endl;
//          cerr << endl;
//          assert(1==2);
//       }
//    }

  return msqrt;
}

void MyMatrix::SchurComplex(const vnl_matrix<double> & M,
    vnl_matrix<vcl_complex< double > > & U,
    vnl_matrix < vcl_complex < double > > & T)
    {
      assert (M.rows() == M.cols());
      long n = M.rows();
      long sdim = 0;
      long lwork = 10*n;
      vnl_vector < vcl_complex < double > > work(lwork);
      vnl_vector < double > rwork(n);
      long info = 0;

      //vnl_matrix<vcl_complex<double> > A(vnl_complexify(M));
        T = vnl_complexify(M);
        vnl_vector < vcl_complex < double > > W(n);
        U.set_size(n,n);

        T.inplace_transpose();// for fortran call switch col and rows
        zgees_ (
            (char*)"V",// JOBVS
            (char*)"N",// SORT
            0,// SELECT
            &n,// N
            T.data_block(), &n,// T, LDA
            &sdim,// SDIM
            W.data_block(),// W
            U.data_block(),&n ,// U, LDVS
            work.data_block(), &lwork,// WORK, LWORK
            rwork.data_block(),// RWORK 
            0,// BWORK
            &info);// info

        T.inplace_transpose();// switch back ..
        U.inplace_transpose();// from fortran ordering

        if (info != 0)
        {
          // These return codes are taken from zgees.f:
          //*          = 0:  successful exit
          //*          < 0:  if INFO = -i, the i-th argument had an illegal value.
          //*          = 1,...,N:
          //*                The QR iteration failed. 
          //*          > N:  =N+1: EV could not be reoredered (ill conditioned)
          //*                =N+2: after reordering, roundoff changed values of
          //*                      some complex eigenvalues so that leading
          //*                      eigenvalues in the Schur form no
          //*                      longer satisfy SELECT=.TRUE.  This could also
          //*                      be caused due to scaling.
          cerr << __FILE__ ": info = " << info << ", something went wrong:\n";
          if (info < 0)
          {
            cerr << __FILE__ ": (internal error) the " << (-info) << "th argument had an illegal value\n";
          }
          else if (1 <= info && info <= n)
          {
            cerr << __FILE__ ": the QR iteration failed, but the last " << (n - info) << " eigenvalues may be correct\n";
          }
          else if (info == n+1)
          {
            cerr << __FILE__ ": EV could not be reordered (ill conditioned)\n";
          }
          else if (info == n+2)
          {
            cerr << __FILE__ ": roundoff error -- maybe due to poor scaling\n";
          }
          else
          {
            cerr << __FILE__ ": unknown error\n";
          }

          assert(info==0);
        }

      }

/**
  Computes mean of rotation matrices 
  */
vnl_matrix_fixed < double, 3, 3> MyMatrix::RotationMean(const std::vector < vnl_matrix_fixed < double, 3, 3 > > &vm)
{
  assert(vm.size() > 0);

  // check if inputs are indeed rotation matrices:
  double maxnorm =0.0;
  for (unsigned int i = 0; i<vm.size(); i++)
  {
    vnl_matrix_fixed < double , 3, 3> Id = vm[i] * vm[i].transpose();
    Id[0][0] -= 1.0; Id[1][1] -= 1.0; Id[2][2] -= 1.0;
    double fn = Id.frobenius_norm() ;
    //std::cout << " matrix " << i << "  frobenius norm: " <<  fn<<std::endl;
    if ( fn > maxnorm) maxnorm = fn;
  }
  std::cout << "Input rotation's max deviation from rotation is: " << maxnorm << std::endl;
  if (maxnorm > 0.0001) 
  {
    std::cerr << " MyMatrix::RotationMean input(s) not rotation" << std::endl;
    exit(1);
  }


  vnl_matrix_fixed<double, 3, 3> meanr(0.0);

 // http://brml.technion.ac.il/publications_files/1307386738.pdf
 // this methods averages the rotation vectors (axis where length = angle)
  std::vector < double > meanRV (3,0.0);
  for (unsigned int i = 0; i<vm.size(); i++)
  {
    Quaternion Q;
    Q.importMatrix(vm[i][0][0],vm[i][0][1],vm[i][0][2],
                   vm[i][1][0],vm[i][1][1],vm[i][1][2],
                   vm[i][2][0],vm[i][2][1],vm[i][2][2] );
    std::vector < double > v = Q.getRotVec();
    //std::cout << "RotVec " << i << " = " << v[0] <<" "<< v[1] << " " << v[2] << std::endl;
    meanRV[0] += v[0];
    meanRV[1] += v[1];
    meanRV[2] += v[2];
  }
  meanRV[0] /= vm.size();
  meanRV[1] /= vm.size();
  meanRV[2] /= vm.size();
  //std::cout << "Mean Rot Vec " << meanRV[0] <<" " << meanRV[1] << " "<< meanRV[2] << std::endl;
  Quaternion Q;
  Q.importRotVec(meanRV[0], meanRV[1], meanRV[2]);
  std::vector < double > rmat = Q.getRotMatrix3d();
  meanr[0][0] = rmat[0]; meanr[0][1] = rmat[1] ; meanr[0][2] = rmat[2];
  meanr[1][0] = rmat[3]; meanr[1][1] = rmat[4] ; meanr[1][2] = rmat[5];
  meanr[2][0] = rmat[6]; meanr[2][1] = rmat[7] ; meanr[2][2] = rmat[8];
  //vnl_matlab_print(vcl_cout,meanr,"meanr",vnl_matlab_print_format_long);std::cout << std::endl;
  return meanr;  

  // This older version only averages the rot matrices (mean is minimizer of Frobenius norm
  // in embedding space) and then maps the result back to a rotation matrix
  // using polar decomposition. This is similar to averaging points on a sphere 
  // in Euclidean space and then projecting the average back.
  for (unsigned int i = 0;i<vm.size(); i++)
    meanr += vm[i];
  // Project meanr back to SO(3) via polar decomposition:  
  meanr = (1.0 / vm.size()) * meanr;
  vnl_matlab_print(vcl_cout,meanr,"meanr",vnl_matlab_print_format_long);std::cout << std::endl;
  vnl_matrix<double> PolR(3, 3), PolS(3, 3);
  MyMatrix::PolarDecomposition(meanr, PolR, PolS);
  meanr = PolR;
    
  return meanr;

}


/**
 Computes geometric mean of matrices in vector vm.
  If n=-1 set n number of elements in vm (default)
  else use n for the division (usefull, if some elements are the identity matrix and not
  passed in the array).
  geo = exp ( sum_i log(mv[i]) / n )
 */
vnl_matrix < double > MyMatrix::GeometricMean(const std::vector < vnl_matrix < double > > & vm, int n)
{
  assert(vm.size() > 0);
  if (n==-1) n=(int)vm.size();
  vnl_matrix < double > geo(MatrixLog(vm[0]));
  for (unsigned int i = 1;i<vm.size();i++)
    geo = geo + MatrixLog(vm[i]);
  geo = MatrixExp((1.0/n) * geo);
  return geo;
}

/**
 Computes geometric mean of matrices in vector vm.
  If n=-1 set n number of elements in vm (default)
  else use n for the division (usefull, if some elements are the identity matrix and not
  passed in the array).
  geo = exp ( sum_i log(mv[i]) / n )
 */
vnl_matrix < double > MyMatrix::GeometricMean(const std::vector < vnl_matrix_fixed < double, 3, 3 > > & vm, int n)
{
  assert(vm.size() > 0);
  if (n==-1) n=(int)vm.size();
  vnl_matrix < double > geo(MatrixLog(vm[0]));
  for (unsigned int i = 1;i<vm.size();i++)
    geo = geo + MatrixLog(vm[i]);
  geo = MatrixExp((1.0/n) * geo);
  return geo;
}


// vnl_matrix < double > MyMatrix::GeometricMean(const std::vector < vnl_matrix < double > > & vm, int n)
// // Computes ( mv[0] * mv[1] ... * mv[m] )^{1/n}
// // by first computing the n-th square root of each factor and then multiplying,
// // to ensure the solution is the one that lies between the inputs
// // (imagine n=7 with six pos 45 degree rotations and one identity,
// // the result X should be 0 < X < 45 degree rotation, not negative, 
// // which happens by multiplying first)
// // If n=-1 set n number of elements in vm (default)
// // else use n for the squareroot (usefull, if some elements are the identity matrix and not
// // passed in the array).
// {
//   assert(vm.size() > 0);
// 
//   if (n==-1) n=(int)vm.size();
//   
//   double d = 1.0/n;
//   vnl_matrix < double > geo(MatrixPow(vm[0],d));
//   
//   // compute n-th sqrt or the rest
//   for (unsigned int i=1;i<vm.size();i++)
//    geo *= MatrixPow(vm[i],d);
// 
//   return geo;
// }

/**
 Computes polar decomposition of sqyare matrix A (using singular value decomposition):
 A = R S into a rotation matrix R and a symmetric matrix S
 @param[in]  A  Matrix to be decomposed
 @param[out] R  Rotation Matrix
 @param[out] S  Symmetric Matrix
 */
void MyMatrix::PolarDecomposition(const vnl_matrix<double> &A,
    vnl_matrix<double> &R, vnl_matrix<double> &S)
{
  assert(A.cols() == A.rows());
  vnl_svd<double> svd(A);
  if (svd.valid())
  {
    R = svd.U() * svd.V().transpose();
    S = svd.V() * svd.W() * svd.V().transpose();
    //vnl_matlab_print(vcl_cout,R,"R",vnl_matlab_print_format_long);cout << endl;
    //vnl_matlab_print(vcl_cout,S,"S",vnl_matlab_print_format_long);cout << endl;
  }
  else
  {
    cerr << "MyMatrix PolarDecomposition ERROR: SVD not possible?" << endl;
    exit(1);
  }

}

/**
 Computes polar decomposition of sqyare matrix A (using singular value decomposition):
 A = R S D into a rotation matrix R and a shear matrix S and diagonal matrix D
 @param[in]  A  Matrix to be decomposed
 @param[out] R  Rotation Matrix
 @param[out] S  Shear Matrix
 @param[out] D  Diagonal Matrix
 */
void MyMatrix::Polar2Decomposition(const vnl_matrix<double> &A,
    vnl_matrix<double> &R, vnl_matrix<double> &S, vnl_diag_matrix<double> &D)
{
  assert(A.cols() == A.rows());
  PolarDecomposition(A, R, S);

  // further decompose S into shear * diag(scales)
  D.set_size(A.cols());
  for (unsigned int c = 0; c < A.cols(); c++)
  {
    D[c] = S[c][c];
    S.set_column(c, S.get_column(c) / D[c]);
  }
  //vnl_matlab_print(vcl_cout,R,"Rot",vnl_matlab_print_format_long);cout << endl;
  //vnl_matlab_print(vcl_cout,S,"Shear",vnl_matlab_print_format_long);cout << endl;
  //vnl_matlab_print(vcl_cout,D,"Scale",vnl_matlab_print_format_long);cout << endl;  
}

/**
 Computes squared distance between a and b (4x4 affine transformations)
 \f$ D^2 = 1/5 r^2 Tr(A^tA) +  ||T_d||^2 \f$
 where \f$T_d\f$ is the translation and \f$A\f$ the Affine
 of the matrix \f$d = a - b\f$.
 r is the radius specifying the volume of interest.
 This distance is used in Jenkinson 1999 RMS deviation - tech report
 www.fmrib.ox.ac.uk/analysis/techrep .
 The center of the brain should be at the origin.
 */
double MyMatrix::AffineTransDistSq(const vnl_matrix<double>&a,
    const vnl_matrix<double>&b, double r)
{

  vnl_matrix<double> drigid(a);
  if (!b.empty())
    drigid -= b;
  else
  {
    vnl_matrix<double> id(drigid.rows(), drigid.cols());
    id.set_identity();
    drigid -= id;
  }

  assert(drigid.rows() ==4 && drigid.cols() == 4);
  assert(fabs(drigid[3][0]) < 0.000001);
  assert(fabs(drigid[3][1]) < 0.000001);
  assert(fabs(drigid[3][2]) < 0.000001);

  //cout << " drigid: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",drigid);

  // translation norm quadrat:
  double tdq = 0;
  for (int i = 0; i < 3; i++)
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
  for (int i = 0; i < 3; i++)
  {
    tr += drigid[i][i];
  }

  return (1.0 / 5.0) * r * r * tr + tdq;

}

/**
 Computes squared distance between a and b (4x4 rigid transformations)
 \f$ D^2 = ||T_d||^2 + |log R_d|^2 \f$
 where \f$ T_df$ is the translation and f$R_df$ the rotation
 of the matrix d that rigidly transforms a to b
 */
double MyMatrix::RigidTransDistSq(const vnl_matrix<double>&a,
    const vnl_matrix<double>&b)
{

  vnl_matrix<double> drigid;
  if (b.empty())
    drigid = a;
  else
  {
    drigid = vnl_matrix_inverse<double>(a);
    drigid = b * drigid;
  }

  assert(drigid.rows() ==4 && drigid.cols() == 4);
  assert(fabs(drigid[3][0]) < 0.000001);
  assert(fabs(drigid[3][1]) < 0.000001);
  assert(fabs(drigid[3][2]) < 0.000001);
  assert(fabs(drigid[3][3]-1) < 0.000001);

  //cout << " drigid: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",drigid);

  // translation norm quadrat:
  double tdq = 0;
  for (int r = 0; r < 3; r++)
  {
    tdq += drigid[r][3] * drigid[r][3];
  }

  //cout << " trans dist2: " << tdq << endl;

  // rotation norm:
  double rd = RotMatrixLogNorm(drigid);
  //cout << " rd: " << rd << endl;

  return rd * rd + tdq;

}

double MyMatrix::getFrobeniusDiff(const vnl_matrix<double>&m1,
    const vnl_matrix<double>&m2)
{

  assert(m1.rows() == m2.rows());
  assert(m1.cols() == m2.cols());

  return (m1 - m2).fro_norm();
}

/**
 The Frobenius norm of the log of a rotation matrix is the geodesic
 distance on rotation matrices.
 The procedure will look only at first three rows and colums and
 expects a rotation matrix there.
 */
double MyMatrix::RotMatrixLogNorm(const vnl_matrix_fixed<double, 4, 4>& m)
{

  // for rotation matrices:
  // tr R = 1 + 2 cos (theta)  , |theta| < pi
  // log R  = 0 for theta =0
  //        = (theta / (2 sin(theta)) ) * ( R - R^T)
  // see also Moakher , Means and averaging in the group of rotations
  

  // assert we have no stretching only rot (and trans)
  double det = vnl_determinant(m);
  //cout << " det: " << det << endl;
  if (fabs(det - 1.0) > 0.001)
  {
    cerr << "There is streching! det: " << det << endl;
    assert(fabs(det-1.0) < 0.001);
  }

  double trace = 0.0;
  for (int n = 0; n < 3; n++)
    trace += m[n][n];
  //cout << " trace : " << trace << endl;
  trace = 0.5 * (trace - 1.0);
  if (trace > 1.0)
    trace = 1.0;
  if (trace < -1.0)
    trace = -1.0;
  //cout << "  0.5*(trace-1): " << trace << endl;
  double theta = acos(trace); // gives [0..pi]

  return sqrt(2.0) * theta;

}

/**
 Assuming tri-linear interpolation, how large
 is the smoothing (averaging) introduced by this transform?
 Source and target geometries are taken from the LTA.
 */
double MyMatrix::getResampSmoothing(const LTA * lta)
{

  MATRIX *V2V;

  if (lta->type == LINEAR_VOX_TO_VOX)
    V2V = MatrixCopy(lta->xforms[0].m_L, NULL);
  else if (lta->type == LINEAR_RAS_TO_RAS)
  {
    V2V = MatrixCopy(lta->xforms[0].m_L, NULL);
    //           sI2R
    //     src -------> RAS
    //      |?           | mod (input)
    //      V            V
    //     dst <------- RAS
    //           dR2I
    MATRIX *sI2R = vg_i_to_r(&lta->xforms[0].src);
    MATRIX *dR2I = vg_r_to_i(&lta->xforms[0].dst);
    MATRIX *tmp = 0;
    if (sI2R == 0 || dR2I == 0)
    {
      cerr
          << "MyMatrix::getResampSmoothing: passed volumes did not have the info on i_to_r or r_to_i."
          << endl;
      exit(1);
    }
    tmp = MatrixMultiply(V2V, sI2R, NULL);
    MatrixMultiply(dR2I, tmp, V2V);
    MatrixFree(&tmp);
    MatrixFree(&sI2R);
    MatrixFree(&dR2I);
  }
  else
  {
    cerr << "MyMatrix::getResampSmoothing  LTA type " << lta->type
        << " not supported!" << endl;
    exit(1);
  }

  // convert to double:
  vnl_matrix<double> V = convertMATRIX2VNL(V2V);
  //vnl_matlab_print(vcl_cout,V,"V",vnl_matlab_print_format_long);cout << endl;
  vnl_matrix<double> Vinv = vnl_inverse(V);
  //vnl_matlab_print(vcl_cout,Vinv,"Vinv",vnl_matlab_print_format_long);cout << endl;

  assert(Vinv.rows() ==4 && Vinv.cols()==4);

  MatrixFree(&V2V);

  int width = lta->xforms[0].dst.width;
  int height = lta->xforms[0].dst.height;
  int depth = lta->xforms[0].dst.depth;
  vnl_vector<double> x(4), y(4);
  x[3] = 1.0;
  y[3] = 1.0;
  double row = 0.0, slice = 0.0, full = 0.0;
  for (int d = 0; d < depth; d++)
  {
    y[2] = d;
    slice = 0.0;
    for (int h = 0; h < height; h++)
    {
      y[1] = h;
      row = 0.0;
      for (int w = 0; w < width; w++)
      {
        y[0] = w;
        x = Vinv * y;

        row += fabs(nint(x[0]) - x[0]) + fabs(nint(x[1]) - x[1])
            + fabs(nint(x[2]) - x[2]);
        //cout << " Y " << y << " X " << x << endl;
        //cout << " err " << row << endl;
        //exit(1);
      }
      slice += row / width;
    }
    full += slice / height;
  }

  full /= depth;
  //full /= 1.5; // the maximum

  return full;

}

vnl_matrix<double> MyMatrix::getVNLMatrix(std::vector<double> d, int r)
// convert double array to matrix
{
  int c = (int) d.size() / r;
  assert(r*c == (int)d.size());

  vnl_matrix<double> m(r, c);

  int rr, cc, count = 0;
  for (rr = 0; rr < r; rr++)
    for (cc = 0; cc < c; cc++)
    {
      m[rr][cc] = d[count];
      count++;
    }

  return m;
}

void MyMatrix::getRTfromM(const vnl_matrix_fixed<double, 4, 4> &m,
    vnl_matrix_fixed<double, 3, 3> &r, vnl_vector_fixed<double, 3> &t)
{
  r = m.extract(3, 3);
  t = m.extract(3, 1, 0, 3).get_column(0);
}

vnl_matrix_fixed<double, 4, 4> MyMatrix::getMfromRT(
    const vnl_matrix_fixed<double, 3, 3> &r,
    const vnl_vector_fixed<double, 3> &t)
{
  vnl_matrix_fixed<double, 4, 4> m;
  m.set_identity();
  for (unsigned int rr = 0; rr < 3; rr++)
  {
    for (unsigned int cc = 0; cc < 3; cc++)
      m[rr][cc] = r[rr][cc];
    m[rr][3] = t[rr];
  }
  return m;
}

/**
 Returns rotation matrix with 1 (or -1 based on sign) in the 
 row specified by i,j,k for each column respectively.
 i,j,k permutations of 1,2,3 (possibly signed, so that det = 1).
 */
vnl_matrix_fixed<double, 4, 4> MyMatrix::getRot(int i, int j, int k)
{

  vnl_matrix_fixed<double, 4, 4> M(0.0);
  M[3][3] = 1.0;

  if (abs(i) + abs(j) + abs(k) != 6)
  {
    cout << "MyMatrix::getRot( " << i << " , " << j << " , " << k
        << " ) ERROR: pass permutation of 1 2 3" << endl;
    exit(1);
  }
  if (i * j * k < 0)
  {
    cout << "MyMatrix::getRot( " << i << " , " << j << " , " << k
        << " ) ERROR: determinant negative" << endl;
    exit(1);
  }

  M[i - 1][0] = sign(i);
  M[j - 1][1] = sign(j);
  M[k - 1][2] = sign(k);

  return M;

}

/** 
 Parametrizes all 24 possible rotations (multiples of 90 degrees in all directions)
 so i = 0..23
 */
vnl_matrix_fixed<double, 4, 4> MyMatrix::getRot(int i)
{

  int j = i / 4; // 0..5
  int k = i % 4; // 0..3

  int a = 1, b = 2, c = 3;

  switch (j)
  {
  // case 0: a=1; b=2; c=3; break;
  case 1:
    a = 1;
    b = 3;
    c = 2;
    break;
  case 2:
    a = 2;
    b = 1;
    c = 3;
    break;
  case 3:
    a = 2;
    b = 3;
    c = 1;
    break;
  case 4:
    a = 3;
    b = 1;
    c = 2;
    break;
  case 5:
    a = 3;
    b = 2;
    c = 1;
    break;
  }

  switch (k)
  {
  //case 0: a= a; b= b; c= c; break;
  case 1:
    b = -b;
    c = -c;
    break;
  case 2:
    a = -a;
    c = -c;
    break;
  case 3:
    a = -a;
    b = -b;
    break;
  }

  return getRot(a, b, c);

}

LTA* MyMatrix::VOXmatrix2LTA(const vnl_matrix_fixed<double, 4, 4>& m, MRI* src,
    MRI* dst)
{
  LTA* ret = LTAalloc(1, src);
  ret->xforms[0].m_L = convertVNL2MATRIX(m, ret->xforms[0].m_L);
//  ret->xforms[0].m_L = MRIvoxelXformToRasXform (src,dst,ret->xforms[0].m_L,ret->xforms[0].m_L) ;
//  ret->type = LINEAR_RAS_TO_RAS ;
  ret->type = LINEAR_VOX_TO_VOX;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

LTA* MyMatrix::RASmatrix2LTA(const vnl_matrix_fixed<double, 4, 4>& m, MRI* src,
    MRI* dst)
{
  LTA* ret = LTAalloc(1, src);
  ret->xforms[0].m_L = convertVNL2MATRIX(m, ret->xforms[0].m_L);
  ret->type = LINEAR_RAS_TO_RAS;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

vnl_matrix<double> MyMatrix::LTA2VOXmatrix(LTA * lta)
{

  vnl_matrix<double> M;
  if (lta->type == LINEAR_RAS_TO_RAS)
  {
    LINEAR_TRANSFORM *lt = &lta->xforms[0];
    MATRIX *m_L = MatrixCopy(lt->m_L, NULL);
    //           sI2R
    //     src -------> RAS
    //      |?           | mod (input)
    //      V            V
    //     dst <------- RAS
    //           dR2I
    MATRIX *sI2R = vg_i_to_r(&lt->src);
    MATRIX *dR2I = vg_r_to_i(&lt->dst);
    MATRIX *tmp = 0;
    if (sI2R == 0 || dR2I == 0)
    {
      cout
          << "MyMatrix::LTA2VOXmatrix ERROR : passed volumes did not have the info on i_to_r or r_to_i."
          << endl;
      exit(1);
    }
    tmp = MatrixMultiply(m_L, sI2R, NULL);
    MatrixMultiply(dR2I, tmp, m_L);
    MatrixFree(&tmp);
    MatrixFree(&sI2R);
    MatrixFree(&dR2I);
    M = convertMATRIX2VNL(m_L);
    MatrixFree(&m_L);
  }
  else if (lta->type == LINEAR_VOX_TO_VOX)
  {
    M = convertMATRIX2VNL(lta->xforms[0].m_L);
  }
  else
  {
    cout << " MyMatrix::LTA2VOXmatrix ERROR lta type not supported !" << endl;
    exit(1);
  }

  return M;
}

vnl_matrix<double> MyMatrix::LTA2RASmatrix(LTA * lta)
{

  vnl_matrix<double> M;
  if (lta->type == LINEAR_VOX_TO_VOX)
  {
    LINEAR_TRANSFORM *lt = &lta->xforms[0];
    MATRIX *m_L = MatrixCopy(lt->m_L, NULL);
    //           sR2I
    //     src <------- RAS
    //      | mod        | ?
    //      V            V
    //     dst -------> RAS
    //           dI2R
    MATRIX *sR2I = vg_r_to_i(&lt->src);
    MATRIX *dI2R = vg_i_to_r(&lt->dst);
    MATRIX *tmp = 0;
    if (sR2I == 0 || dI2R == 0)
    {
      cout
          << "MyMatrix::LTA2RASmatrix ERROR: passed volumes did not have the info on i_to_r or r_to_i."
          << endl;
      exit(1);
    }
    tmp = MatrixMultiply(m_L, sR2I, NULL);
    MatrixMultiply(dI2R, tmp, m_L);
    MatrixFree(&tmp);
    MatrixFree(&sR2I);
    MatrixFree(&dI2R);
    M = convertMATRIX2VNL(m_L);
    MatrixFree(&m_L);
  }
  else if (lta->type == LINEAR_RAS_TO_RAS)
  {
    M = convertMATRIX2VNL(lta->xforms[0].m_L);
  }
  else
  {
    cout << " MyMatrix::LTA2VOXmatrix ERROR lta type not supported !" << endl;
    exit(1);
  }

  return M;
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

  MATRIX* drigid = MatrixCopy(a, NULL);
  if (b)
    drigid = MatrixSubtract(drigid, b, drigid);
  else
  {
    MATRIX *id = MatrixIdentity(4, NULL);
    drigid = MatrixSubtract(drigid, id, drigid);
    MatrixFree(&id);
  }

  assert(drigid->rows ==4 && drigid->cols == 4);
  assert(fabs(drigid->rptr[4][1]) < 0.000001);
  assert(fabs(drigid->rptr[4][2]) < 0.000001);
  assert(fabs(drigid->rptr[4][3]) < 0.000001);

  //cout << " drigid: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",drigid);

  // translation norm quadrat:
  double tdq = 0;
  for (int i = 1; i <= 3; i++)
  {
    tdq += drigid->rptr[i][4] * drigid->rptr[i][4];
    drigid->rptr[i][4] = 0.0;
    drigid->rptr[4][i] = 0.0;
  }
  drigid->rptr[4][4] = 0.0;

  //cout << " trans dist2: " << tdq << endl;
  MATRIX* dt = MatrixTranspose(drigid, NULL);
  drigid = MatrixMultiply(dt, drigid, drigid);
  MatrixFree(&dt);

  // Trace of A^t A
  double tr = 0.0;
  for (int i = 1; i <= 3; i++)
  {
    tr += drigid->rptr[i][i];
  }

  MatrixFree(&drigid);

  return (1.0 / 5.0) * r * r * tr + tdq;

}

double MyMatrix::RigidTransDistSq(MATRIX * a, MATRIX * b)
// computes squared distance between a and b (4x4 rigid transformations)
// D^2 = ||T_d||^2 + |log R_d|^2
// where T_d is the translation and R_d the rotation
// of the matrix d that rigidly transforms a to b
{

  MATRIX* drigid;
  if (!b)
    drigid = MatrixCopy(a, NULL);
  else
  {
    drigid = MatrixInverse(a, NULL);
    drigid = MatrixMultiply(b, drigid, drigid);
  }

  assert(drigid->rows ==4 && drigid->cols == 4);
  assert(fabs(drigid->rptr[4][1]) < 0.000001);
  assert(fabs(drigid->rptr[4][2]) < 0.000001);
  assert(fabs(drigid->rptr[4][3]) < 0.000001);
  assert(fabs(drigid->rptr[4][4]-1) < 0.000001);

  //cout << " drigid: " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",drigid);

  // translation norm quadrat:
  double tdq = 0;
  for (int r = 1; r <= 3; r++)
  {
    tdq += drigid->rptr[r][4] * drigid->rptr[r][4];
  }

  //cout << " trans dist2: " << tdq << endl;

  // rotation norm:
  double rd = RotMatrixLogNorm(drigid);
  //cout << " rd: " << rd << endl;

  MatrixFree(&drigid);

  return rd * rd + tdq;

}

double MyMatrix::getFrobeniusDiff(MATRIX *m1, MATRIX *m2)
{

  double s, ss = 0.0;
  assert(m1->rows == m2->rows);
  assert(m1->cols == m2->cols);

  for (int r = 1; r <= m1->rows; r++)
    for (int c = 1; c <= m2->cols; c++)
    {
      s = *MATRIX_RELT(m1, r, c) - *MATRIX_RELT(m2, r, c);
      ss += s * s;
    }
  ss = sqrt(ss);
  return ss;
}

MATRIX * MyMatrix::MatrixSqrtIter(MATRIX * m, MATRIX *msqrt)
{
  assert(m->rows == 4 && m->cols == 4);
  msqrt = MatrixIdentity(4, msqrt);
  MATRIX* R = MatrixAlloc(3, 3, MATRIX_REAL);
  for (int rr = 1; rr <= 3; rr++)
    for (int cc = 1; cc <= 3; cc++)
    {
      *MATRIX_RELT(R, rr, cc) = *MATRIX_RELT(m, rr, cc);
    }

  //Denman and Beavers square root iteration

  int imax = 100;
  double eps = 0.0001;
  double err = 1000;
  //cout << "using square root iteartion (" << imax << ")"<< endl;
  MATRIX * Yn = MatrixCopy(R, NULL);
  MATRIX * Zn = MatrixIdentity(3, NULL);
  MATRIX * Zni = NULL;
  MATRIX * Yni = NULL;
  MATRIX * Ysq = NULL;
  int count = 0;
  while (count < imax && err > eps)
  {
    count++;
    Yni = MatrixInverse(Yn, Yni);
    Zni = MatrixInverse(Zn, Zni);
    assert(Yni && Zni);

    Yn = MatrixAdd(Yn, Zni, Yn);
    Zn = MatrixAdd(Zn, Yni, Zn);

    Yn = MatrixScalarMul(Yn, 0.5, Yn);
    Zn = MatrixScalarMul(Zn, 0.5, Zn);
    //cout << " matrix " << i << endl;
    //MatrixPrintFmt(stdout,"% 2.8f",Yn);
    //cout << endl;

    Ysq = MatrixMultiply(Yn, Yn, Ysq);
    Ysq = MatrixSubtract(Ysq, R, Ysq);
    err = 0;
    for (int c = 1; c < 4; c++)
      for (int r = 1; r < 4; r++)
        err += fabs(*MATRIX_RELT(Ysq, r, c));

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
  MATRIX* Rh1 = MatrixCopy(Rh, NULL);
  *MATRIX_RELT(Rh1, 1, 1) = *MATRIX_RELT(Rh1, 1, 1) + 1;
  *MATRIX_RELT(Rh1, 2, 2) = *MATRIX_RELT(Rh1, 2, 2) + 1;
  *MATRIX_RELT(Rh1, 3, 3) = *MATRIX_RELT(Rh1, 3, 3) + 1;

  VECTOR * T = MatrixAlloc(3, 1, MATRIX_REAL);
  *MATRIX_RELT(T, 1, 1) = *MATRIX_RELT(m, 1, 4);
  *MATRIX_RELT(T ,2, 1) = *MATRIX_RELT(m, 2, 4);
  *MATRIX_RELT(T, 3, 1) = *MATRIX_RELT(m, 3, 4);

  MATRIX* Rh1i = MatrixInverse(Rh1, NULL);
  assert(Rh1i);

  VECTOR * Th = MatrixMultiply(Rh1i, T, NULL);

  *MATRIX_RELT(msqrt, 1, 4) = *MATRIX_RELT(Th, 1, 1);
  *MATRIX_RELT(msqrt, 2, 4) = *MATRIX_RELT(Th, 2, 1);
  *MATRIX_RELT(msqrt, 3, 4) = *MATRIX_RELT(Th, 3, 1);
  *MATRIX_RELT(msqrt, 4, 1) = 0;
  *MATRIX_RELT(msqrt, 4, 2) = 0;
  *MATRIX_RELT(msqrt, 4, 3) = 0;
  *MATRIX_RELT(msqrt, 4, 4) = 1;
  for (int c = 1; c < 4; c++)
    for (int r = 1; r < 4; r++)
      *MATRIX_RELT(msqrt, r, c) = *MATRIX_RELT(Rh, r, c);

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

/**
 Converts std double vector (r x c) into a matrix (FS MATRIX)
 */
MATRIX* MyMatrix::getMatrix(std::vector<double> d, int r, int c, MATRIX* m)
{
  if (c == -1)
    c = r; // quadratic

  assert(r*c == (int)d.size());
  if (!m)
    m = MatrixAlloc(r, c, MATRIX_REAL);
  assert(m->rows ==r);
  assert(m->cols ==c);

  int rr, cc, count = 0;
  for (rr = 1; rr <= r; rr++)
    for (cc = 1; cc <= c; cc++)
    {
      *MATRIX_RELT(m, rr,cc) = d[count];
      count++;
    }

  return m;
}

/**
 Converts affine vector (12x1) into an affine matrix (homogeneous coord) 4x4 (FS MATRIX)
 */
MATRIX * MyMatrix::aff2mat(MATRIX * aff, MATRIX *outM)
{
  if (outM == NULL)
    outM = MatrixAlloc(4, 4, MATRIX_REAL);
  MatrixIdentity(4, outM);

  int count = 1;
  for (int rr = 1; rr <= 3; rr++)
    for (int cc = 1; cc <= 4; cc++)
    {
      *MATRIX_RELT(outM, rr, cc) = *MATRIX_RELT(outM, rr, cc)
          + *MATRIX_RELT(aff, count, 1);
      count++;
    }

  return outM;
}

/**
 Computes half the rotation and translation (FS MATRIX).
 m must be rotation and translation only (rigid!)
 */
MATRIX * MyMatrix::getHalfRT(MATRIX * m, MATRIX *mhalf)
{
  if (mhalf)
    MatrixFree(&mhalf);

  float d = MatrixDeterminant(m);
  d = fabs(d-1);
  assert(d < 0.000001);

  Quaternion q;
  q.importMatrix(*MATRIX_RELT(m, 1, 1), *MATRIX_RELT(m, 1, 2),
      *MATRIX_RELT(m, 1, 3), *MATRIX_RELT(m, 2, 1), *MATRIX_RELT(m, 2, 2),
      *MATRIX_RELT(m, 2, 3), *MATRIX_RELT(m, 3, 1), *MATRIX_RELT(m, 3, 2),
      *MATRIX_RELT(m, 3, 3));
  //cout << "q: "<< q << endl;
  Quaternion qh = q.getHalfRotation();
  //cout << "qh: " << qh << endl;
  mhalf = getMatrix(qh.getRotMatrix3dh(), 4, 4);
  MATRIX* Rh1 = MatrixIdentity(3, NULL);
  for (int rr = 1; rr < 4; rr++)
    for (int cc = 1; cc < 4; cc++)
      *MATRIX_RELT(Rh1, rr, cc) = *MATRIX_RELT(Rh1, rr, cc)
          + *MATRIX_RELT(mhalf, rr, cc);

  VECTOR * T = MatrixAlloc(3, 1, MATRIX_REAL);
  *MATRIX_RELT(T, 1, 1) = *MATRIX_RELT(m, 1, 4);
  *MATRIX_RELT(T ,2, 1) = *MATRIX_RELT(m, 2, 4);
  *MATRIX_RELT(T, 3, 1) = *MATRIX_RELT(m, 3, 4);

  //cout << " rh1" << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",Rh1);

  MATRIX* Rh1i = MatrixInverse(Rh1, NULL);
  assert(Rh1i);

  VECTOR * Th = MatrixMultiply(Rh1i, T, NULL);

  *MATRIX_RELT(mhalf, 1, 4) = *MATRIX_RELT(Th, 1, 1);
  *MATRIX_RELT(mhalf, 2, 4) = *MATRIX_RELT(Th, 2, 1);
  *MATRIX_RELT(mhalf, 3, 4) = *MATRIX_RELT(Th, 3, 1);
  *MATRIX_RELT(mhalf, 4, 1) = 0;
  *MATRIX_RELT(mhalf, 4, 2) = 0;
  *MATRIX_RELT(mhalf, 4, 3) = 0;
  *MATRIX_RELT(mhalf, 4, 4) = 1;

  MatrixFree(&Th);
  MatrixFree(&Rh1i);
  MatrixFree(&T);
  MatrixFree(&Rh1);
  return mhalf;
}

/**
 Computes Frobenius norm of log of rot matrix (FS MATRIX).
 This is equivalent to geodesic distance on rot matrices.
 The procedure will look only at first three rows and colums and
 expects a rotation matrix there.
 */
double MyMatrix::RotMatrixLogNorm(MATRIX * m)
{
  // assert we have no stretching only rot (and trans)
  float det = MatrixDeterminant(m);
  //cout << " det: " << det << endl;
  if (fabs(det - 1.0) > 0.001)
  {
    cerr << "There is streching! det: " << det << endl;
    assert(fabs(det-1.0) < 0.001);
  }

  double trace = 0.0;
  for (int n = 1; n <= 3; n++)
    trace += m->rptr[n][n];
  //cout << " trace : " << trace << endl;
  trace = 0.5 * (trace - 1.0);
  if (trace > 1.0)
    trace = 1.0;
  if (trace < -1.0)
    trace = -1.0;
  //cout << "  0.5*(trace-1): " << trace << endl;
  double theta = acos(trace); // gives [0..pi]

  return sqrt(2.0) * theta;

}

double MyMatrix::RotMatrixGeoDist(MATRIX * a, MATRIX *b)
{

  if (!b)
    return RotMatrixLogNorm(a);

  // if not 3x3, fetch first 3x3
  // and construct a^T b
  MATRIX *at = MatrixAlloc(3, 3, MATRIX_REAL);
  MATRIX *blocal = MatrixAlloc(3, 3, MATRIX_REAL);
  assert(a->rows >= 3 && a->cols >= 3);
  assert(b->rows >= 3 && b->cols >= 3);
  for (int r = 1; r <= 3; r++)
    for (int c = 1; c <= 3; c++)
    {
      at->rptr[r][c] = a->rptr[c][r];
      blocal->rptr[r][c] = b->rptr[r][c];
    }

  blocal = MatrixMultiply(at, blocal, blocal);

  double dist = RotMatrixLogNorm(blocal);

  MatrixFree(&at);
  MatrixFree(&blocal);

  return dist;

}

pair<MATRIX *, VECTOR *> MyMatrix::getRTfromM(MATRIX * M, MATRIX * R,
    VECTOR * T)
{
  // check dimenstions:
  assert(M->rows == 4);
  assert(M->cols == 4);
  if (R == NULL)
    R = MatrixAlloc(3, 3, M->type);
  if (T == NULL)
    T = VectorAlloc(3,M->type);
  assert(R->rows ==3);
  assert(R->cols ==3);
  assert(T->rows ==3);
  assert(T->cols ==1);

  // check M
  assert(fabs(M->rptr[4][1]) < 0.000001);
  assert(fabs(M->rptr[4][2]) < 0.000001);
  assert(fabs(M->rptr[4][3]) < 0.000001);
  assert(fabs(M->rptr[4][4] - 1.0) < 0.000001);

  for (int c = 1; c < 4; c++)
  {
    for (int r = 1; r < 4; r++)
    {
      R->rptr[r][c] = M->rptr[r][c];
    }
    T->rptr[c][1] = M->rptr[c][4];
  }

  return pair<MATRIX *, VECTOR *>(R, T);
}

MATRIX * MyMatrix::getMfromRT(MATRIX * R, VECTOR * T, MATRIX * M)
{
  int type = 0;
  if (R != NULL)
    type = R->type;
  else if (T != NULL)
    type = T->type;
  else
    assert(R != NULL || T != NULL);

  if (M == NULL)
    M = MatrixAlloc(4, 4, type);
  if (R == NULL)
    R = MatrixIdentity(3, NULL);
  if (T == NULL)
    T = MatrixZero(3, 1, NULL);

  // check dimensions:
  assert(M->rows == 4);
  assert(M->cols == 4);
  assert(R->rows ==3);
  assert(R->cols ==3);
  assert(T->rows ==3);
  assert(T->cols ==1);

  for (int c = 1; c < 4; c++)
  {
    for (int r = 1; r < 4; r++)
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
  LTA* ret = LTAalloc(1, src);
  ret->xforms[0].m_L = MRIvoxelXformToRasXform(src, dst, m, NULL);
  ret->type = LINEAR_RAS_TO_RAS;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

LTA* MyMatrix::RASmatrix2LTA(MATRIX * m, MRI* src, MRI* dst)
{
  LTA* ret = LTAalloc(1, src);
  ret->xforms[0].m_L = MatrixCopy(m, ret->xforms[0].m_L);
  ret->type = LINEAR_RAS_TO_RAS;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

