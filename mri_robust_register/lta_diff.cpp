/**
 * @file  lta_diff.cpp
 * @brief A programm to compute differences of lta files (transforms)
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2012/12/19 01:50:52 $
 *    $Revision: 1.26 $
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
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>
#include <limits>
#include <vcl_iostream.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_matlab_print.h>

#include "Registration.h"
#include "MyMatrix.h"

// all other software are all in "C"
#ifdef __cplusplus
extern "C"
{
#endif
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#include "transform.h"

#ifdef __cplusplus
}
#endif

using namespace std;

//static char vcid[] = "$Id: lta_diff.cpp,v 1.26 2012/12/19 01:50:52 mreuter Exp $";
char *Progname = NULL;
void writeVox2Vox(LTA * lta)
{
  cout << " convet to vox 2 vox" << endl;
  LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  cout << " writing" << endl;
  LTAwrite(lta, "test-vox2vox.lta");
}

double cornerdiff(LTA* lta1, LTA* lta2)
{
  // get vox2vox using lta geometry info
  LTAchangeType(lta1, LINEAR_VOX_TO_VOX);
  LTAchangeType(lta2, LINEAR_VOX_TO_VOX);

  VECTOR * v_X = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  VECTOR * v_Y1 = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  VECTOR * v_Y2 = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  VECTOR_ELT(v_X,4)= 1;
  VECTOR_ELT(v_Y1,4)= 1;
  VECTOR_ELT(v_Y2,4)= 1;

  assert(lta1->xforms[0].src.depth == lta2->xforms[0].src.depth);
  assert(lta1->xforms[0].src.height == lta2->xforms[0].src.height);
  assert(lta1->xforms[0].src.width == lta2->xforms[0].src.width);

  int y3, y2, y1;
  double d = 0;
  for (y3 = 0; y3 < 2; y3++)
  {
    V3_Z(v_X)= y3 * (lta1->xforms[0].src.depth-1);
    for (y2 = 0; y2 < 2; y2++)
    {
      V3_Y(v_X) = y2 * (lta1->xforms[0].src.height-1);
      for (y1 = 0; y1 < 2; y1++)
      {
        V3_X(v_X) = y1* (lta1->xforms[0].src.width-1);
        MatrixMultiply(lta1->xforms[0].m_L, v_X, v_Y1);
        MatrixMultiply(lta2->xforms[0].m_L, v_X, v_Y2);
        double d1 = V3_X(v_Y1) - V3_X(v_Y2);
        double d2 = V3_Y(v_Y1) - V3_Y(v_Y2);
        double d3 = V3_Z(v_Y1) - V3_Z(v_Y2);
        d += sqrt(d1*d1 + d2*d2 + d3*d3);
        //cout << " corner : " << V3_X(v_X) << " , " <<  V3_Y(v_X) << " , " <<  V3_Z(v_X) << endl;
        //cout << "   mapped to "<< V3_X(v_Y1) << " , " <<  V3_Y(v_Y1) << " , " <<  V3_Z(v_Y1) << endl;
        //cout << "   mapped to "<< V3_X(v_Y2) << " , " <<  V3_Y(v_Y2) << " , " <<  V3_Z(v_Y2) << endl;
      }
    }
  }

  return d / 8.0;
}

double cornerdiff(LTA* lta1)
{
  // get vox2vox using lta geometry info
  LTAchangeType(lta1, LINEAR_VOX_TO_VOX);

  VECTOR * v_X = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  VECTOR * v_Y1 = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  VECTOR_ELT(v_X,4)= 1;
  VECTOR_ELT(v_Y1,4)= 1;

  int y3, y2, y1;
  double d = 0;
  double dmax = 0;
  for (y3 = 0; y3 < 2; y3++)
  {
    V3_Z(v_X)= y3 * (lta1->xforms[0].src.depth-1);
    for (y2 = 0; y2 < 2; y2++)
    {
      V3_Y(v_X) = y2 * (lta1->xforms[0].src.height-1);
      for (y1 = 0; y1 < 2; y1++)
      {
        V3_X(v_X) = y1* (lta1->xforms[0].src.width-1);
        MatrixMultiply(lta1->xforms[0].m_L, v_X, v_Y1);
        double d1 = V3_X(v_Y1) - V3_X(v_X);
        double d2 = V3_Y(v_Y1) - V3_Y(v_X);
        double d3 = V3_Z(v_Y1) - V3_Z(v_X);
        double dd = sqrt(d1*d1 + d2*d2 + d3*d3);
        //cout << " dd: " << dd << endl;
        if ( dd > dmax) dmax = dd;
        d += dd;
        //cout << " corner : " << V3_X(v_X) << " , " <<  V3_Y(v_X) << " , " <<  V3_Z(v_X) << endl;
        //cout << "   mapped to "<< V3_X(v_Y1) << " , " <<  V3_Y(v_Y1) << " , " <<  V3_Z(v_Y1) << endl;
        //cout << "   mapped to "<< V3_X(v_Y2) << " , " <<  V3_Y(v_Y2) << " , " <<  V3_Z(v_Y2) << endl;
      }
    }
  }
  cout << " dmax: " << dmax << endl;
  return d / 8.0;
}

double determinant(MATRIX * M1, MATRIX* M2)
{

//   MATRIX* M = MatrixAlloc(4,4,MATRIX_REAL);
  MATRIX* M = MatrixCopy(M1, NULL);
  if (M2 != NULL)
  {
    //cout << " inverting" << endl;
    //M = MatrixInverse(M1,M);
    M = MatrixMultiply(M2, M, M);
  }

  double d = MatrixDeterminant(M);
  MatrixFree(&M);
  return d;
}

void testQuaternion(const vnl_matrix<double>& Rot)
{

  Quaternion Q;
  Q.importMatrix(Rot[0][0],Rot[0][1],Rot[0][2],Rot[1][0],Rot[1][1],Rot[1][2],Rot[2][0],Rot[2][1],Rot[2][2]);
  cout << "Quaternion: " << Q << endl;
  
  vnl_vector < double > v(3);
  v[0] = 1.2; v[1]=-0.5; v[2] = .7;
  cout << " v    :  " << v << endl;
  cout << " v rot:  " << Rot * v << endl;
  
  std::vector < double > v1 = Q.rotate(v[0],v[1],v[2]);
  cout << " v rotQ: " << v1[0] <<" " << v1[1] << " " << v1[2] << endl;
  
  std::vector < double > m = Q.getRotMatrix3d();
  cout << " M = [ " << m[0] << " " << m[1] << " " << m[2] << endl << m[3] << " " << m[4] << " " << m[5] << endl << m[6] << " " << m[7] << " " << m[8] << " ]" << endl; 

exit(1);
}

void decompose(MATRIX * M1, MATRIX* M2)
{
  vnl_matrix<double> m = MyMatrix::convertMATRIX2VNL(M1);

  if (M2 != NULL)
  {
    //cout << " inverting" << endl;
    //M = MatrixInverse(M1,M);
    vnl_matrix<double> m2 = MyMatrix::convertMATRIX2VNL(M2);
    m = m * m2;
  }

  cout << " Decompose M1*M2 into Rot * Shear * Scale + Trans: " << endl
      << endl;
      
  vnl_matrix<double> Rot, Shear;
  vnl_diag_matrix<double> Scale;
  MyMatrix::Polar2Decomposition(m.extract(3, 3), Rot, Shear, Scale);
  vnl_matlab_print(vcl_cout,Rot,"Rot",vnl_matlab_print_format_long);
  cout << endl;
  
  Quaternion Q;
  Q.importMatrix(Rot[0][0],Rot[0][1],Rot[0][2],Rot[1][0],Rot[1][1],Rot[1][2],Rot[2][0],Rot[2][1],Rot[2][2]);
  std::vector < double > v = Q.getRotVec();
  cout << "RotVec = [ " << v[0] << " " << v[1] << " " << v[2] << " ] " << endl;
  cout << endl;
  
  cout << "RotAngle = " << Q.getRotAngle() << endl;
  cout << endl;
  
  vnl_matlab_print(vcl_cout,Shear,"Shear",vnl_matlab_print_format_long);
  cout << endl;
  
  vnl_matlab_print(vcl_cout,Scale,"Scale",vnl_matlab_print_format_long);
  cout << endl;
  
  vnl_vector<double> t = m.extract(3, 1, 0, 3).get_column(0);
  vnl_matlab_print(vcl_cout,t,"Trans",vnl_matlab_print_format_long);
  cout << endl;
  
  cout << "AbsTrans = " << t.two_norm() << endl;
  cout << endl;
  
  cout << "Determinant = " << vnl_determinant(m) << endl << endl;

}

double sphereDiff(MATRIX * M1, MATRIX* M2, double r)
{

  MATRIX* M = MatrixAlloc(4, 4, MATRIX_REAL);
  if (M2 == NULL)
    M = MatrixCopy(M1, M);
  else
  {
    M = MatrixInverse(M1, M);
    M = MatrixMultiply(M2, M, M);
  }

  // MatrixPrintFmt(stdout,"% 2.8f",M);

  double dmax = 0;
  double dmin = std::numeric_limits<double>::infinity();
  double davg = 0;
  int counter = 0;

  VECTOR * v_A = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  VECTOR * v_B = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  VECTOR_ELT(v_A,4)= 1;
  VECTOR_ELT(v_B,4)= 1;

  // poles
  V3_X(v_A)= 0;
  V3_Y(v_A)= 0;
  V3_Z(v_A)= r;
  MatrixMultiply(M, v_A, v_B);
  // MatrixPrintFmt(stdout,"% 2.8f",v_A);
  // MatrixPrintFmt(stdout,"% 2.8f",v_B);

  double d1 = V3_X(v_B)- V3_X(v_A);
  double d2 = V3_Y(v_B)- V3_Y(v_A);
  double d3 = V3_Z(v_B)- V3_Z(v_A);
  double dd = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
  dmax = dd;
  dmin = dd;
  davg = dd;
  counter++;
  V3_X(v_A)= 0;
  V3_Y(v_A)= 0;
  V3_Z(v_A)= -r;
  MatrixMultiply(M, v_A, v_B);
  d1 = V3_X(v_B)- V3_X(v_A);
  d2 = V3_Y(v_B)- V3_Y(v_A);
  d3 = V3_Z(v_B)- V3_Z(v_A);
  dd = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
  if (dd > dmax)
    dmax = dd;
  if (dd < dmin)
    dmin = dd;
  davg += dd;
  counter++;

  // sample sphere with radius 1
  int max = 10;
  for (int i = -max + 1; i < max; i++)
  {
    double angle1 = (i * M_PI * 0.5) / max; // from -pi/2 to +pi/2
    // radius:
    double r1 = cos(angle1);
    double h = sin(angle1);
    // circumference is 2pi *r1, we want 4*max samples at aequator (where r=1 and cc 2pi)
    int max2 = int(4.0 * max * r1);
    for (int j = 0; j < max2; j++)
    {
      double angle2 = (2.0 * M_PI * j) / max2; // from 0 to 2pi
      V3_X(v_A)= r*r1*cos(angle2);
      V3_Y(v_A)= r*r1*sin(angle2);
      V3_Z(v_A)= r*h;
      MatrixMultiply(M, v_A, v_B);
      d1 = V3_X(v_B)- V3_X(v_A);
      d2 = V3_Y(v_B)- V3_Y(v_A);
      d3 = V3_Z(v_B)- V3_Z(v_A);
      dd = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
      if (dd > dmax)
        dmax = dd;
      if (dd < dmin)
        dmin = dd;
      davg += dd;
      counter++;
    }
  }
  davg = davg / counter;
//   cout << " max: " << dmax << " min: " << dmin << " avg: " << davg << endl;
  MatrixFree(&M);
  return dmax;
}

void testSphereDiff()
{

  MATRIX* T = MatrixIdentity(4, NULL);
  T->rptr[1][4] = 10.0;
  cout << sphereDiff(T, NULL, 100) << endl;

  // rotation 90degree around z axis
  MATRIX* R = MatrixIdentity(4, NULL);
  R->rptr[1][1] = 0;
  R->rptr[1][2] = -1;
  R->rptr[2][1] = 1;
  R->rptr[2][2] = 0;
  cout << sphereDiff(R, NULL, 100) << endl;

  exit(1);
}

double interpolationError2D(double angle)
{
  int side = 256;

  //MRI* mri_error = MRIalloc(side, side,1,MRI_FLOAT);
  MATRIX *a = MatrixAllocRotation(3, angle, Z_ROTATION);

  VECTOR * v_X = VectorAlloc(3, MATRIX_REAL);  // input (src) coordinates
  VECTOR * v_Y = VectorAlloc(3, MATRIX_REAL);  // transformed (dst) coordinates
  double errorsum = 0, x, y;
  int xm, xp, ym, yp;
  double val, xmd, ymd, xpd, ypd;  // d's are distances
  V3_Z(v_Y)= 0;
  for (int y1 = 0; y1 < side; y1++)
    for (int y2 = 0; y2 < side; y2++)
    {
      V3_X(v_Y)= y1 - 128;
      V3_Y(v_Y) = y2 - 128;
      MatrixMultiply(a, v_Y, v_X);

      x = V3_X(v_X) + 128;
      y = V3_Y(v_X) + 128;

//         xm = MAX((int)x, 0) ;
//         xp = MIN(side-1, xm+1) ;
//         ym = MAX((int)y, 0) ;
//         yp = MIN(side-1, ym+1) ;
      xm = (int)floor(x);
      xp = xm+1;
      ym = (int)floor(y);
      yp = ym+1;

      xmd = x - (float)xm;
      ymd = y - (float)ym;
      xpd = (1.0f - xmd);
      ypd = (1.0f - ymd);

      //cout << "x: " << x << " xm: " << xm << " xp: " << xp << " xmd: " << xmd << " xpd: " << xpd << endl;
      //cout << "y: " << y <<" ym: " << ym << " yp: " << yp << " ymd: " << ymd << " ypd: " << ypd << endl;
      //assert(x>=0);
      assert (xmd >= 0 && xpd >= 0);
      assert (ymd >= 0 && ypd >= 0);

      val = 0;// sum of distance to each coordinate (use smallest)
      if (xmd < xpd) val += xmd;
      else val += xpd;
      if (ymd < ypd) val += ymd;
      else val += ypd;
      //MRIFvox(mri_error,y1,y2,0) = (float)(val) ;
      errorsum += val;

    }
    //MRIwrite(mri_error,"mri_error.mgz");
    //MRIfree(&mri_error);
  return errorsum;
}

double interpolationError(LTA* lta)
{
  // get vox2vox using lta geometry info
  LTAchangeType(lta, LINEAR_VOX_TO_VOX);

  // sample from dst back to src
  MATRIX *mAinv = MatrixInverse(lta->xforms[0].m_L, NULL);
  if (!mAinv)
    ErrorExit(ERROR_BADPARM, "interpolationError: xform is singular");
  int width = lta->xforms[0].dst.width;
  int height = lta->xforms[0].dst.height;
  int depth = lta->xforms[0].dst.depth;
  VECTOR * v_X = VectorAlloc(4, MATRIX_REAL);  // input (src) coordinates
  VECTOR * v_Y = VectorAlloc(4, MATRIX_REAL);  // transformed (dst) coordinates
  int y3, y2, y1;
  double x, y, z;
  int xm, xp, ym, yp, zm, zp;
  double val, xmd, ymd, zmd, xpd, ypd, zpd;  // d's are distances

  MRI* mri_error = MRIalloc(width, height, depth, MRI_FLOAT);

  double errorsum = 0;
  v_Y->rptr[4][1] = 1.0f;
  for (y3 = 0; y3 < depth; y3++)
  {
    V3_Z(v_Y)= y3;
    for (y2 = 0; y2 < height; y2++)
    {
      V3_Y(v_Y) = y2;
      for (y1 = 0; y1 < width; y1++)
      {
        V3_X(v_Y) = y1;
        MatrixMultiply(mAinv, v_Y, v_X);

        x = V3_X(v_X);
        y = V3_Y(v_X);
        z = V3_Z(v_X);

        xm = MAX((int)x, 0);
        xp = MIN(width-1, xm+1);
        ym = MAX((int)y, 0);
        yp = MIN(height-1, ym+1);
        zm = MAX((int)z, 0);
        zp = MIN(depth-1, zm+1);

        xmd = x - (float)xm;
        ymd = y - (float)ym;
        zmd = z - (float)zm;
        xpd = (1.0f - xmd);
        ypd = (1.0f - ymd);
        zpd = (1.0f - zmd);

        val = 0; // sum of distance to each coordinate (use smallest)
        if (xmd < xpd) val += xmd;
        else val += xpd;
        if (ymd < ypd) val += ymd;
        else val += ypd;
        if (zmd < zpd) val += zmd;
        else val += zpd;
        MRIFvox(mri_error,y1,y2,y3) = (float)(val);
        errorsum += val;
      }
    }
  }
  MRIwrite(mri_error, "mri_error.mgz");
  MRIfree(&mri_error);
  return errorsum;

}

int main(int argc, char *argv[])
{

//testSphereDiff();

//   // Default initialization
//   int nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
//   if (nargs && argc - nargs == 1) exit (0);
//   argc -= nargs;
  Progname = argv[0];
//   argc --;
//   argv++;
//  if (vcid)
//  {};

  if (argc < 3)
  {
    cout << endl;
    cout << argv[0] << " file1.lta file2.lta [dist-type] [norm-div] [invert]"
        << endl;
    cout << endl;
    cout << "    dist-type " << endl;
    cout << "       1            Rigid Transform Distance (||log(R)|| + ||T||)"
        << endl;
    cout << "       2  (default) Affine Transform Distance (RMS) " << endl;
    cout << "       3            8-corners mean distance after transform "
        << endl;
    cout << "       4            Max Displacement on Sphere " << endl;
    cout << "       5            Determinant (scaling) M1*M2" << endl;
    cout
        << "       6            Interpolation Smoothing (only for first transform)"
        << endl;
    cout << "       7            Decomposition RAS M1*M2 = Rot*Shear*Scaling"
        << endl;
    cout << "       8            Decomposition VOX M1*M2 = Rot*Shear*Scaling"
        << endl;
    cout << "                       pass 'identity.nofile' for second lta "
        << endl;
    cout
        << "    norm-div  (=1)  divide final distance by this (e.g. step adjustment)"
        << endl;
    cout
        << "    invert          invert LTA: 0 false (default), 1 invert first, 2 second"
        << endl;
    cout << endl;
    cout
        << " For the RMS and Sphere we use a ball with radius 10cm at the RAS center."
        << endl;
    cout << " Instead of 'file2.lta' you can specify 'identity.nofile'."
        << endl;
    cout << endl;
    exit(1);
  }
  string lta1f = argv[1];
  string lta2f = argv[2];

  double d = 1.0;
  int disttype = 2;
  if (argc > 3)
  {
    disttype = atoi(argv[3]);
    assert(double(disttype) == atof(argv[3]));
  }
  if (argc > 4)
    d = atof(argv[4]);

  int invert = 0;
  if (argc > 5)
    invert = atoi(argv[5]);

  if (disttype == 100)
  {
    int steps = 200;
    double div = 16.0;
    vector<double> theta(steps);
    vector<double> err(steps);
    for (int i = 0; i < steps; i++)
    {
      // 0.. PI/div in 20 steps
      // -PI/div ..0 is symmetric
      theta[i] = M_PI * (i + 1) / ((steps) * div);
      err[i] = interpolationError2D(theta[i]);
    }
    ostringstream ss;
    ss << "interror-rot16";
    string fn = ss.str() + ".plot";
    ofstream f(fn.c_str(), ios::out);

    f << "set terminal postscript eps color" << endl;
    f << "set title \"Interpolation error when rotating \"" << endl;
    f << "set output \"" << ss.str() << ".eps\"" << endl;
    f << "plot  \"-\" notitle with lines 1" << endl;
    for (int i = 0; i < steps; i++)
    {
      cout << theta[i] << " " << err[i] << endl;
      f << theta[i] << " " << err[i] << endl;
    }
    f << "e" << endl;
    exit(0);
  }

  LTA* lta1 = LTAreadEx(lta1f.c_str());
  LTA* lta2 = LTAreadEx(lta2f.c_str());

  if (!lta1 || !lta2)
  {
    cerr << "Could not open one of the LTA input files" << endl;
    exit(1);
  }

  if (invert == 1)
  {
    VOL_GEOM vgtmp;
    LT *lt;
    MATRIX *m_tmp = lta1->xforms[0].m_L;
    lta1->xforms[0].m_L = MatrixInverse(lta1->xforms[0].m_L, NULL);
    MatrixFree(&m_tmp);
    lt = &lta1->xforms[0];
    if (lt->dst.valid == 0 || lt->src.valid == 0)
    {
      fprintf(stderr,
          "WARNING:********************************************************\n");
      fprintf(stderr,
          "WARNING:dst or src volume is invalid.  Inverse likely wrong.\n");
      fprintf(stderr,
          "WARNING:********************************************************\n");
    }
    copyVolGeom(&lt->dst, &vgtmp);
    copyVolGeom(&lt->src, &lt->dst);
    copyVolGeom(&vgtmp, &lt->src);
  }
  if (invert == 2)
  {
    VOL_GEOM vgtmp;
    LT *lt;
    MATRIX *m_tmp = lta2->xforms[0].m_L;
    lta2->xforms[0].m_L = MatrixInverse(lta2->xforms[0].m_L, NULL);
    MatrixFree(&m_tmp);
    lt = &lta2->xforms[0];
    if (lt->dst.valid == 0 || lt->src.valid == 0)
    {
      fprintf(stderr,
          "WARNING:********************************************************\n");
      fprintf(stderr,
          "WARNING:dst or src volume is invalid.  Inverse likely wrong.\n");
      fprintf(stderr,
          "WARNING:********************************************************\n");
    }
    copyVolGeom(&lt->dst, &vgtmp);
    copyVolGeom(&lt->src, &lt->dst);
    copyVolGeom(&vgtmp, &lt->src);
  }

  double dist = -1;
  //Registration R;
  LTAchangeType(lta1,LINEAR_VOX_TO_VOX);
  MATRIX* VOX1 = MatrixCopy(lta1->xforms[0].m_L,NULL);
  LTAchangeType(lta1, LINEAR_RAS_TO_RAS);
  MATRIX* RAS1 = MatrixCopy(lta1->xforms[0].m_L,NULL);
  MATRIX* RAS2 = NULL;
  MATRIX* VOX2 = NULL;
  if (lta2f != "identity.nofile")
  {
    LTAchangeType(lta2,LINEAR_VOX_TO_VOX);
    VOX2 = MatrixCopy(lta2->xforms[0].m_L,NULL);
    LTAchangeType(lta2, LINEAR_RAS_TO_RAS);
    RAS2 = MatrixCopy(lta2->xforms[0].m_L,NULL);
  }

  switch (disttype)
  {
  case 1:
    dist = sqrt(MyMatrix::RigidTransDistSq(RAS1, RAS2)) / d;
    break;
  case 2:
    dist = sqrt(MyMatrix::AffineTransDistSq(RAS1, RAS2)) / d;
    break;
  case 3:
    if (lta2f == "identity.nofile")
      dist = cornerdiff(lta1) / d;
    else
      dist = cornerdiff(lta1, lta2) / d;
    break;
  case 4:
    dist = sphereDiff(RAS1, RAS2, 100) / d;
    break;
  case 5:
    dist = determinant(RAS1, RAS2) / d;
    break;
  case 6:
    dist = MyMatrix::getResampSmoothing(lta1) / d;
    break;
  case 7:
    decompose(RAS1, RAS2);
    exit(0);
    break;
  case 8:
    decompose(VOX1, VOX2);
    exit(0);
    break;
  default:
    cerr << "ERROR: dist-type " << disttype << " unknown!" << endl;
    exit(1);
    break;
  }
  if (disttype < 7)
    cout << dist << endl;
    
  MatrixFree(&VOX1);
  MatrixFree(&RAS1);
  if (VOX2) MatrixFree(&VOX2);
  if (RAS2) MatrixFree(&RAS2);
}
