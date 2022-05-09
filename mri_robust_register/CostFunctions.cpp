/**
 * @brief A class that makes available many different cost functions for images
 *   and to combine multiple volumes by mean or median
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

#include "CostFunctions.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

#include "RobustGaussian.h"

#include "error.h"
#include "macros.h"
#include "mrimorph.h"
#include "matrix.h"

#define export // obsolete feature
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_det.h>
#undef export

using namespace std;

std::pair<float, float> CostFunctions::minmax(MRI *i)
{
  MRIiterator it1(i);
  it1.begin();
  float min = (*it1);
  float max = min;
  for (it1++; !it1.isEnd(); it1++)
  {
    if ((*it1) < min)
      min = (*it1);
    if ((*it1) > max)
      max = (*it1);
  }
  return std::pair<float, float>(min, max);
}

double CostFunctions::mean(MRI *mri, int frame)
{
/*  int count = 0;
  double d = 0.0;
  MRIiterator it1(mri);
  for (it1.begin(); !it1.isEnd(); it1++)
  {
    d += (*it1);
    count++;
  }
  return (float) (d / count);*/
  

// currently only for one frame at a time
// future, return vector < double > for all frames and do 
// mapping only once, then each process needs it's own vector < double > finally 
// sum across processess (basically reduction on a vector)

  double d = 0.0;
  unsigned int ocount = 0;
  unsigned int count = 0;

  int z;
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:d)  
#endif
  for (z = 0; z < mri->depth; z++)
  {
    int y,x;
    double v;
    for (y = 0; y < mri->height; y++)
    {
      for (x = 0; x < mri->width; x++)
      {
         //MRIsampleVolumeFrame(mri, x, y, z, frame, &v);
          v = MRIgetVoxVal(mri,x,y,z,frame);
          if (v == -1)
          {
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
            ocount++;
            continue;
          }
          d += v;  
          
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
         count++; // can be removed, but needs to agree with N below!!!
      }
    }
  }
  
  unsigned int n =  mri->width * mri->height * mri->depth - ocount;
  assert(n==count);
  return d/n;


}


double CostFunctions::mean(MRI * mri,
    const vnl_matrix_fixed<double, 4, 4>& Mi,
    int frame,
    int d1, int d2, int d3)
{
  return mean(mri,Mi,frame,d1,d2,d3,0,mri->width,0,mri->height,0,mri->depth);
}

double CostFunctions::mean(MRI * mri,
    const vnl_matrix_fixed<double, 4, 4>& Mi,
    int frame,
    int d1, int d2, int d3,
    int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
// currently only for one frame at a time
// future, return vector < double > for all frames and do 
// mapping only once, then each process needs it's own vector < double > finally 
// sum across processess (basically reduction on a vector)

  mri->outside_val = -1;
  if (frame <0 || frame >= mri->nframes )
  {
    std::cerr << " ERROR: CostFunctions::mean frame " << frame << " outside range !" << std::endl;
    exit(1);
  }
  int dt[3] = { xmax, ymax, zmax };
  dt[0] = dt[0] -d1 +1;
  dt[1] = dt[1] -d2 +1;
  dt[2] = dt[2] -d3 +1;

  int z;
  double d = 0.0;
  unsigned int ocount = 0;
  unsigned int count = 0;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:d)  
#endif
  for (z = zmin; z < dt[2]; z += d3)
  {
    int x, y;
    double xx, yx, zx;
    double v;
    double xz,yz,zz;
    double xy,yy,zy;
  
    xz = Mi[0][2] * z + Mi[0][3];
    yz = Mi[1][2] * z + Mi[1][3];
    zz = Mi[2][2] * z + Mi[2][3];
    for (y = ymin; y < dt[1]; y += d2)
    {
      xy = Mi[0][1] * y + xz;
      yy = Mi[1][1] * y + yz;
      zy = Mi[2][1] * y + zz;
      for (x = xmin; x < dt[0]; x += d1)
      {

        xx = Mi[0][0] * x + xy;
        yx = Mi[1][0] * x + yy;
        zx = Mi[2][0] * x + zy;

        //for (f = 0; f < dt[3]; f++)
        //{
          MRIsampleVolumeFrame(mri, xx, yx, zx, frame, &v);
          if (v == -1)
          {
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
            ocount++;
            continue;
          }
          d += v;  
          
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
           count++; // can be removed, but needs to agree with N below!!!
        //}
      }
    }
  }
  
  unsigned int n = (((xmax-xmin-1)/d1)+1) * (((ymax-ymin-1)/d2)+1)  * (((zmax-zmin-1)/d3)+1) -ocount;
  assert(n==count);
  return d/n;

}



double CostFunctions::var(MRI *i, int frame)
{
  double m = mean(i,frame);
  double d = 0.0;
  double dd;
  int count = 0;
  MRIiterator it1(i);
  for (it1.begin(); !it1.isEnd(); it1++)
  {
    dd = (*it1) - m;
    d += dd * dd;
    count++;
  }
  return  (d / count);
}


double CostFunctions::norm(MRI *mri, int frame)
{
  double d = 0.0;
  int z;
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:d)  
#endif
  for (z = 0; z < mri->depth; z++)
  {
    int y,x;
    double v;
    for (y = 0; y < mri->height; y++)
    {
      for (x = 0; x < mri->width; x++)
      {
          v = MRIgetVoxVal(mri,x,y,z,frame);
          d += v * v;  
      }
    }
  }
  return sqrt(d);
}

float CostFunctions::median(MRI *i)
{
  int n = i->width * i->height * i->depth;
  float* t = (float *) calloc(n, sizeof(float));

  int cc = 0;
  // copy array to t
  MRIiterator it1(i);
  for (it1.begin(); !it1.isEnd(); it1++)
  {
    t[cc] = (*it1);
    cc++;
  }

  float qs = RobustGaussian<float>::median(t, n);
  free(t);
  return qs;
}

float CostFunctions::mad(MRI *i, float d)
{
  int n = i->width * i->height * i->depth;
  float* t = (float *) calloc(n, sizeof(float));

  int cc = 0;
  // copy array to t
  MRIiterator it1(i);
  for (it1.begin(); !it1.isEnd(); it1++)
  {
    t[cc] = (*it1);
    cc++;
  }

  float qs = RobustGaussian<float>::mad(t, n, d);
  free(t);
  return qs;
}

/*double CostFunctions::leastSquares(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3)
{
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  double dt[4] = { mriT->width, mriT->height, mriT->depth , mriT->nframes };
  int x, y, z, f;
  double dd, d=0;
  double xs, ys, zs;
  double xt, yt, zt;
  double vs, vt;
  for (z = 0; z < dt[2] - d3 + 1; z += d3)
  {
    for (y = 0; y < dt[1] - d2 + 1; y += d2)
    {
      for (x = 0; x < dt[0] - d1 + 1; x += d1)
      {

        xt = Mti[0][0] * x + Mti[0][1] * y + Mti[0][2] * z + Mti[0][3];
        yt = Mti[1][0] * x + Mti[1][1] * y + Mti[1][2] * z + Mti[1][3];
        zt = Mti[2][0] * x + Mti[2][1] * y + Mti[2][2] * z + Mti[2][3];

        xs = Msi[0][0] * x + Msi[0][1] * y + Msi[0][2] * z + Msi[0][3];
        ys = Msi[1][0] * x + Msi[1][1] * y + Msi[1][2] * z + Msi[1][3];
        zs = Msi[2][0] * x + Msi[2][1] * y + Msi[2][2] * z + Msi[2][3];

        for (f = 0; f < dt[3]; f++)
        {
          MRIsampleVolumeFrame(mriS, xs, ys, zs, f, &vs);
          if (vs == -1) continue;
          MRIsampleVolumeFrame(mriT, xt, yt, zt, f, &vt);
          if (vt == -1) continue;
        
          dd = vs-vt;
          d += dd*dd;  
        }
      }
    }
  }

  return d;

}*/

double CostFunctions::leastSquares(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3)
{
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  int dt[4] = { mriT->width, mriT->height, mriT->depth , mriT->nframes };
  dt[0] = dt[0] -d1 +1;
  dt[1] = dt[1] -d2 +1;
  dt[2] = dt[2] -d3 +1;

  int z;
  double d = 0.0;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:d)  
#endif
  for (z = 0; z < dt[2]; z += d3)
  {
    int x, y, f;
    double dd;
    double xs, ys, zs;
    double xt, yt, zt;
    double vs, vt;
    double xtz,ytz,ztz,xsz,ysz,zsz;
    double xty,yty,zty,xsy,ysy,zsy;
  
    xtz = Mti[0][2] * z + Mti[0][3];
    ytz = Mti[1][2] * z + Mti[1][3];
    ztz = Mti[2][2] * z + Mti[2][3];
    xsz = Msi[0][2] * z + Msi[0][3];
    ysz = Msi[1][2] * z + Msi[1][3];
    zsz = Msi[2][2] * z + Msi[2][3];
    for (y = 0; y < dt[1]; y += d2)
    {
      xty = Mti[0][1] * y + xtz;
      yty = Mti[1][1] * y + ytz;
      zty = Mti[2][1] * y + ztz;
      xsy = Msi[0][1] * y + xsz;
      ysy = Msi[1][1] * y + ysz;
      zsy = Msi[2][1] * y + zsz;
      for (x = 0; x < dt[0]; x += d1)
      {

        xt = Mti[0][0] * x + xty;
        yt = Mti[1][0] * x + yty;
        zt = Mti[2][0] * x + zty;
        xs = Msi[0][0] * x + xsy;
        ys = Msi[1][0] * x + ysy;
        zs = Msi[2][0] * x + zsy;

        for (f = 0; f < dt[3]; f++)
        {
          MRIsampleVolumeFrame(mriS, xs, ys, zs, f, &vs);
          if (vs == -1) continue;
          MRIsampleVolumeFrame(mriT, xt, yt, zt, f, &vt);
          if (vt == -1) continue;
        
          dd = vs-vt;
          d += dd*dd;  
        }
      }
    }
  }

  return d;

}

double CostFunctions::leastSquares(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3, 
    const double &s1, const double &s2)
{
  // white background needs to be tested, 
  // in my test case there were black stripes, 
  // at boundary, leading to gray stripes after
  // downsampling, leading to big problems during alingment
  double sout = mriS->outside_val;
  double tout = mriT->outside_val;
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  int dt[4] = { mriT->width, mriT->height, mriT->depth , mriT->nframes };
  dt[0] = dt[0] -d1 +1;
  dt[1] = dt[1] -d2 +1;
  dt[2] = dt[2] -d3 +1;

  int z;
  double d = 0.0;
  //double oepss = sout / 255.0;
  //double oepst = tout / 255.0;
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:d)  
#endif
  for (z = 0; z < dt[2]; z += d3)
  {
    int x, y, f;
    double dd;
    double xs, ys, zs;
    double xt, yt, zt;
    double vs, vt;
    double xtz,ytz,ztz,xsz,ysz,zsz;
    double xty,yty,zty,xsy,ysy,zsy;
  
    xtz = Mti[0][2] * z + Mti[0][3];
    ytz = Mti[1][2] * z + Mti[1][3];
    ztz = Mti[2][2] * z + Mti[2][3];
    xsz = Msi[0][2] * z + Msi[0][3];
    ysz = Msi[1][2] * z + Msi[1][3];
    zsz = Msi[2][2] * z + Msi[2][3];
    for (y = 0; y < dt[1]; y += d2)
    {
      xty = Mti[0][1] * y + xtz;
      yty = Mti[1][1] * y + ytz;
      zty = Mti[2][1] * y + ztz;
      xsy = Msi[0][1] * y + xsz;
      ysy = Msi[1][1] * y + ysz;
      zsy = Msi[2][1] * y + zsz;
      for (x = 0; x < dt[0]; x += d1)
      {

        xt = Mti[0][0] * x + xty;
        yt = Mti[1][0] * x + yty;
        zt = Mti[2][0] * x + zty;
        xs = Msi[0][0] * x + xsy;
        ys = Msi[1][0] * x + ysy;
        zs = Msi[2][0] * x + zsy;

        for (f = 0; f < dt[3]; f++)
        {
          MRIsampleVolumeFrame(mriS, xs, ys, zs, f, &vs);
          if (vs == -1) continue;
          //if (fabs (vs -sout) < oepss) continue;
          MRIsampleVolumeFrame(mriT, xt, yt, zt, f, &vt);
          if (vt == -1) continue;
          //if (fabs (vs -sout) < oepss) vs = 0;
          //if (fabs(vt -tout) < oepst) vt= 0;
        
          dd = (s1*vs)-(s2*vt);
          d += dd*dd;  
        }
      }
    }
  }
  mriS->outside_val = sout;
  mriT->outside_val = tout;
  return d;

}

double CostFunctions::absDiff(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3, 
    const double &s1, const double &s2)
{
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  int dt[4] = { mriT->width, mriT->height, mriT->depth , mriT->nframes };
  dt[0] = dt[0] -d1 +1;
  dt[1] = dt[1] -d2 +1;
  dt[2] = dt[2] -d3 +1;

  int z;
  double d = 0.0;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:d)  
#endif
  for (z = 0; z < dt[2]; z += d3)
  {
    int x, y, f;
    double dd;
    double xs, ys, zs;
    double xt, yt, zt;
    double vs, vt;
    double xtz,ytz,ztz,xsz,ysz,zsz;
    double xty,yty,zty,xsy,ysy,zsy;
  
    xtz = Mti[0][2] * z + Mti[0][3];
    ytz = Mti[1][2] * z + Mti[1][3];
    ztz = Mti[2][2] * z + Mti[2][3];
    xsz = Msi[0][2] * z + Msi[0][3];
    ysz = Msi[1][2] * z + Msi[1][3];
    zsz = Msi[2][2] * z + Msi[2][3];
    for (y = 0; y < dt[1]; y += d2)
    {
      xty = Mti[0][1] * y + xtz;
      yty = Mti[1][1] * y + ytz;
      zty = Mti[2][1] * y + ztz;
      xsy = Msi[0][1] * y + xsz;
      ysy = Msi[1][1] * y + ysz;
      zsy = Msi[2][1] * y + zsz;
      for (x = 0; x < dt[0]; x += d1)
      {

        xt = Mti[0][0] * x + xty;
        yt = Mti[1][0] * x + yty;
        zt = Mti[2][0] * x + zty;
        xs = Msi[0][0] * x + xsy;
        ys = Msi[1][0] * x + ysy;
        zs = Msi[2][0] * x + zsy;

        for (f = 0; f < dt[3]; f++)
        {
          MRIsampleVolumeFrame(mriS, xs, ys, zs, f, &vs);
          if (vs == -1) continue;
          MRIsampleVolumeFrame(mriT, xt, yt, zt, f, &vt);
          if (vt == -1) continue;
        
          dd = (s1*vs)-(s2*vt);
          d += fabs(dd);  
        }
      }
    }
  }

  return d;

}


// needs to be updated (see parallel versions above)
double CostFunctions::leastSquares(MRI_BSPLINE *mriS, MRI_BSPLINE* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3)
{
  mriS->coeff->outside_val = -1;
  mriT->coeff->outside_val = -1;
  double dt[4] = { (double)mriT->coeff->width, (double)mriT->coeff->height, (double)mriT->coeff->depth , (double)mriT->coeff->nframes };
  int x, y, z, f;
  double dd, d=0;
  double xs, ys, zs;
  double xt, yt, zt;
  double vs, vt;
  for (z = 0; z < dt[2] - d3 + 1; z += d3)
  {
    for (y = 0; y < dt[1] - d2 + 1; y += d2)
    {
      for (x = 0; x < dt[0] - d1 + 1; x += d1)
      {

        xt = Mti[0][0] * x + Mti[0][1] * y + Mti[0][2] * z + Mti[0][3];
        yt = Mti[1][0] * x + Mti[1][1] * y + Mti[1][2] * z + Mti[1][3];
        zt = Mti[2][0] * x + Mti[2][1] * y + Mti[2][2] * z + Mti[2][3];

        xs = Msi[0][0] * x + Msi[0][1] * y + Msi[0][2] * z + Msi[0][3];
        ys = Msi[1][0] * x + Msi[1][1] * y + Msi[1][2] * z + Msi[1][3];
        zs = Msi[2][0] * x + Msi[2][1] * y + Msi[2][2] * z + Msi[2][3];

        for (f = 0; f < dt[3]; f++)
        {
          MRIsampleBSpline(mriS, xs, ys, zs, f, &vs);
          if (vs == -1) continue;
          MRIsampleBSpline(mriT, xt, yt, zt, f, &vt);
          if (vt == -1) continue;
        
          dd = vs-vt;
          d += dd*dd;  
        }
      }
    }
  }

  return d;

}

/*float CostFunctions::leastSquares(MRI * i1, MRI * i2)
{
  assert(i1 != NULL);

  if (i2)
  {
    assert(i1->width == i2->width);
    assert(i1->height == i2->height);
    assert(i1->depth == i2->depth);
  }

  double d = 0;
  //cout << "CostFunctions::leastSquares chunk data" <<endl;
  double dd;
  if (i2 == NULL)
  {
    MRIiterator it1(i1);
    for (it1.begin(); !it1.isEnd(); it1++)
      d += (*it1) * (*it1);
  }
  else
  {
    MRIiterator it1(i1);
    MRIiterator it2(i2);
    it2.begin();
    //assert(i1->type == i2->type);
    for (it1.begin(); !it1.isEnd(); it1++)
    {
      //cout << "it1: " << *it1 << " it2: " << *it2 << endl;
      dd = (double) (*it1) - (double) (*it2);
      d += dd * dd;
      it2++;
    }

  }
  //cout << " d: " << d << endl;
  return (float) d;
}*/


double CostFunctions::localNCC(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3 )
{
  int blockradius = 2;
  int bw = 2*blockradius+1;
  //int bw2 = bw*bw;
  //int bw3 = bw*bw2;
  
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  
  // because we may need to iterate several times, map once 
  // to avoid additional resampling
  MRI * nmriS = mapMRI(mriS,Msi,1,1,1);
  MRI * nmriT = mapMRI(mriT,Mti,1,1,1);

  int dt[4] = { mriT->width, mriT->height, mriT->depth , mriT->nframes };
  dt[0] = dt[0] -d1 - blockradius;
  dt[1] = dt[1] -d2 - blockradius;
  dt[2] = dt[2] -d3 - blockradius;

  int z;
  double d = 0.0;
  unsigned int dcount = 0;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:d)  
#endif
  for (z = blockradius; z < dt[2]; z += d3)
  {
    int x, y, f,i,j,k;
    double dd;
    //double xs, ys, zs;
    //double xt, yt, zt;
    double vs, vt;
    
    // precompute xyz coords for all z in block (include translation)
    int nz = z - blockradius;
    int nzm = nz + bw;
    int ny, nx, nym, nxm;
    
    double smean, tmean, svar, tvar,tdd,sdd;
    int counter;
    
    for (y = blockradius; y < dt[1]; y += d2)
    {
      ny = y - blockradius;
      nym = ny + bw;
      for (x = blockradius; x < dt[0]; x += d1)
      {
        nx = x - blockradius;
        nxm = nx+bw;
        
        for (f = 0; f < dt[3]; f++)
        {
          // loop over block to compute means
          smean = 0.0;
          tmean = 0.0;
          counter = 0;
          for (i=nz;i<nzm;i++)
          for (j=ny;j<nym;j++)
          for (k=nx;k<nxm;k++)
          {
            //MRIsampleVolumeFrame(mriS, xsx[k]+xsy[j]+xsz[i], ysx[k]+ysy[j]+ysz[i], zsx[k]+zsy[j]+zsz[i], f, &vs);
            vs = MRIgetVoxVal(nmriS,k,j,i,f);
            if (vs == -1) goto nextcoord;
            //MRIsampleVolumeFrame(mriT, xtx[k]+xty[j]+xtz[i], ytx[k]+yty[j]+ytz[i], ztx[k]+zty[j]+ztz[i], f, &vt);
            vt = MRIgetVoxVal(nmriT,k,j,i,f);
            if (vt == -1) goto nextcoord;
            tmean += vt;
            smean += vs;
            //tdata[counter]= vt;
            //sdata[counter]= vs;
            counter++;
          }
          tmean /= counter;
          smean /= counter;
          
          // loop over block again
          dd = 0.0;
          svar = 0.0;
          tvar = 0.0;
          for (i=nz;i<nzm;i++)
          for (j=ny;j<nym;j++)
          for (k=nx;k<nxm;k++)
          {    
            vs = MRIgetVoxVal(nmriS,k,j,i,f);
            if (vs == -1) goto nextcoord;
            //MRIsampleVolumeFrame(mriT, xtx[k]+xty[j]+xtz[i], ytx[k]+yty[j]+ytz[i], ztx[k]+zty[j]+ztz[i], f, &vt);
            vt = MRIgetVoxVal(nmriT,k,j,i,f);
            if (vt == -1) goto nextcoord;
            tdd = vt - tmean;
            sdd = vs - smean;
            dd += tdd * sdd;
            svar += tdd * tdd;
            tvar += sdd * sdd;
          }
          //svar /= bw3;
          //tvar /= bw3;
            
          //d += dd / (bw3 * svar * tvar); 
          if (svar ==0 || tvar == 0) goto nextcoord;
          svar = sqrt( svar);
          tvar = sqrt( tvar);
          d += fabs(dd / (svar * tvar)); 
          
          if (std::isnan(d))
          {
            cout << "d is nan " << dd << " " << svar << " " << tvar << endl;
            exit (1);
          }
          
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
          dcount++;
        }
// I hate goto's but this is the easiest to break out of all the nested loops above
nextcoord:
      i=0; // needs something here 
      }
    }
    
  }

  MRIfree(&nmriS);
  MRIfree(&nmriT);

  return d/dcount;

}

MRI * CostFunctions::mapMRI(MRI *mri, const vnl_matrix_fixed<double, 4, 4>& Mi,
    int d1, int d2, int d3)
{
  int width  = ((mri->width  -1) / d1 ) +1;
  int height = ((mri->height -1) / d2 ) +1;
  int depth  = ((mri->depth  -1) / d3 ) +1;
  MRI * nmri = MRIallocSequence(width,height,depth,MRI_FLOAT,mri->nframes);
  nmri->outside_val = mri->outside_val;
  
  int iz;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) 
#endif
  for (iz = 0; iz < depth; iz ++)
  {
    int z = iz * d3; 
    int ix,x, iy,y, f;
    double xt, yt, zt;
    double v;
    double xtz,ytz,ztz;
    double xty,yty,zty;
  
    xtz = Mi[0][2] * z + Mi[0][3];
    ytz = Mi[1][2] * z + Mi[1][3];
    ztz = Mi[2][2] * z + Mi[2][3];
    
    for (iy = 0; iy < height; iy++)
    {
      y = iy * d2;
      xty = Mi[0][1] * y + xtz;
      yty = Mi[1][1] * y + ytz;
      zty = Mi[2][1] * y + ztz;
  
      for (ix = 0; ix < width; ix++)
      {
        x = ix * d1;
        xt = Mi[0][0] * x + xty;
        yt = Mi[1][0] * x + yty;
        zt = Mi[2][0] * x + zty;

        for (f = 0; f < mri->nframes; f++)
        {
          MRIsampleVolumeFrame(mri, xt, yt, zt, f, &v);
          MRIsetVoxVal(nmri,ix,iy,iz,f,v);
        }
      }
    }
  }
  return nmri;
}


double CostFunctions::NCC(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3)
{
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  
  // because we need to iterate several times, map once 
  // to avoid additional resampling
  MRI * nmriS = mapMRI(mriS,Msi,d1,d2,d3);
  MRI * nmriT = mapMRI(mriT,Mti,d1,d2,d3);
  
  // compute means for each frame
  vector < double > meanS (mriS->nframes);
  vector < double > meanT (mriT->nframes);
  assert(mriS->nframes == mriT->nframes);
  for (int f = 0; f<mriT->nframes; f++)
  {
    meanS[f] = mean(nmriS,f);
    meanT[f] = mean(nmriT,f);
  }

  assert ( mriT->nframes ==1 ); // for now  

  int z;
  double d = 0.0;
  double sigS = 0.0;
  double sigT = 0.0;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:d)  
#endif
  for (z = 0; z < nmriT->depth; z++)
  {
    int x, y, f;
    double ds, dt, vs, vt;
    
    for (y = 0; y < nmriT->height; y++)
    {
      for (x = 0; x < nmriT->width; x++)
      {

        for (f = 0; f < nmriT->nframes; f++)
        {
          //MRIsampleVolumeFrame(nmriS, x, y, z, f, &vs);
          vs = MRIgetVoxVal(nmriS,x,y,z,f);
          if (vs == -1) continue;
          //MRIsampleVolumeFrame(nmriT, x, y, z, f, &vt);
          vt = MRIgetVoxVal(nmriT,x,y,z,f);
          if (vt == -1) continue;
        
           
          ds = (vs-meanS[f]);
          dt = (vt-meanT[f]);
          d += ds*dt;
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
          sigS += ds*ds;
          
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
          sigT += dt*dt;
        }
      }
    }
  }

  sigS = sqrt(sigS);
  sigT = sqrt(sigT);

  MRIfree(&nmriS);
  MRIfree(&nmriT);

  return fabs(d / (sigS * sigT));

}


double CostFunctions::tukeyBiweight(MRI * i1, MRI * i2, double sat)
{
  assert(i1 != NULL);

  if (i2)
  {
    assert(i1->width == i2->width);
    assert(i1->height == i2->height);
    assert(i1->depth == i2->depth);
  }

  int n = i1->width * i1->height * i1->depth;
  float* diff = (float *) calloc(n, sizeof(float));

  int cc = 0;
  if (i2 == NULL)
  {
    MRIiterator it1(i1);
    for (it1.begin(); !it1.isEnd(); it1++)
      diff[cc] = (*it1);
    cc++;
  }
  else
  {
    MRIiterator it1(i1);
    MRIiterator it2(i2);
    it2.begin();
    //assert(i1->type == i2->type);
    for (it1.begin(); !it1.isEnd(); it1++)
    {
      //cout << "it1: " << *it1 << " it2: " << *it2 << endl;
      diff[cc] = (*it1) - (*it2);
      //if (isnan(diff[cc])) cout << "it1: " << *it1 << " it2: " << *it2 << endl;
      ////if (diff[cc] != 0.0) cout << "it1: " << *it1 << " it2: " << *it2 << endl;
      cc++;
      it2++;
    }

  }

  // the problem with looking at the sigma of the difference is
  // a) that it can be zero (since most voxels are in the background and
  // using uchar makes them and also most differences zero)
  // b) for uchar images it can switch from from 0 to 1 to 2
  // producing jumps in the cost function
  // more thought needs to go into this in the future
  exit(1);

  float sigma = RobustGaussian<float>::mad(diff, n);
  //if (sigma == 0.0)
  sigma = 1.4826;
  cout << "sigma: " << sigma << endl;
  double d = 0;
  for (int i = 0; i < n; i++)
    d += rhoTukeyBiweight(diff[i] / sigma, sat);

  free(diff);
  cout << " d: " << d << endl;
  return d;
}



/*double CostFunctions::tukeyBiweight(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3, double sat)
{
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  int dt[4] = { mriT->width, mriT->height, mriT->depth , mriT->nframes };
  dt[0] = dt[0] -d1 +1;
  dt[1] = dt[1] -d2 +1;
  dt[2] = dt[2] -d3 +1;

  int z;
  double d = 0.0;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)  reduction(+:d) 
#endif
  for (z = 0; z < dt[2]; z += d3)
  {
    int x, y, f;
    double dd;
    double xs, ys, zs;
    double xt, yt, zt;
    double vs, vt;
    double xtz,ytz,ztz,xsz,ysz,zsz;
    double xty,yty,zty,xsy,ysy,zsy;
  
    xtz = Mti[0][2] * z + Mti[0][3];
    ytz = Mti[1][2] * z + Mti[1][3];
    ztz = Mti[2][2] * z + Mti[2][3];
    xsz = Msi[0][2] * z + Msi[0][3];
    ysz = Msi[1][2] * z + Msi[1][3];
    zsz = Msi[2][2] * z + Msi[2][3];
    for (y = 0; y < dt[1]; y += d2)
    {
      xty = Mti[0][1] * y + xtz;
      yty = Mti[1][1] * y + ytz;
      zty = Mti[2][1] * y + ztz;
      xsy = Msi[0][1] * y + xsz;
      ysy = Msi[1][1] * y + ysz;
      zsy = Msi[2][1] * y + zsz;
      for (x = 0; x < dt[0]; x += d1)
      {

        xt = Mti[0][0] * x + xty;
        yt = Mti[1][0] * x + yty;
        zt = Mti[2][0] * x + zty;
        xs = Msi[0][0] * x + xsy;
        ys = Msi[1][0] * x + ysy;
        zs = Msi[2][0] * x + zsy;

        for (f = 0; f < dt[3]; f++)
        {
          MRIsampleVolumeFrame(mriS, xs, ys, zs, f, &vs);
          if (vs == -1) continue;
          MRIsampleVolumeFrame(mriT, xt, yt, zt, f, &vt);
          if (vt == -1) continue;
        
          dd = vs-vt;
          d += rhoTukeyBiweight(dd,sat);  
        }
      }
    }
  }

  return d;

}*/

double CostFunctions::tukeyBiweight(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3, double sat)
{
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  int dt[4] = { mriT->width, mriT->height, mriT->depth , mriT->nframes };
  dt[0] = dt[0] -d1 +1;
  dt[1] = dt[1] -d2 +1;
  dt[2] = dt[2] -d3 +1;

  int z;

  unsigned int n = ((mriT->width+1)/d1) * ((mriT->height+1)/d2) * ((mriT->depth+1)/d3) * mriT->nframes;
  float* diff = (float *) calloc(n, sizeof(float));
  unsigned int pn = ((mriT->width+1)/d1) * ((mriT->height+1)/d2) *mriT->nframes;


#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)  
#endif
  for (z = 0; z < dt[2]; z += d3)
  {
    int x, y, f;
    double dd;
    double xs, ys, zs;
    double xt, yt, zt;
    double vs, vt;
    double xtz,ytz,ztz,xsz,ysz,zsz;
    double xty,yty,zty,xsy,ysy,zsy;
    
    unsigned int count = 0;
    unsigned int zpos = (z/d3) * pn;
    
    xtz = Mti[0][2] * z + Mti[0][3];
    ytz = Mti[1][2] * z + Mti[1][3];
    ztz = Mti[2][2] * z + Mti[2][3];
    xsz = Msi[0][2] * z + Msi[0][3];
    ysz = Msi[1][2] * z + Msi[1][3];
    zsz = Msi[2][2] * z + Msi[2][3];
    for (y = 0; y < dt[1]; y += d2)
    {
      xty = Mti[0][1] * y + xtz;
      yty = Mti[1][1] * y + ytz;
      zty = Mti[2][1] * y + ztz;
      xsy = Msi[0][1] * y + xsz;
      ysy = Msi[1][1] * y + ysz;
      zsy = Msi[2][1] * y + zsz;
      for (x = 0; x < dt[0]; x += d1)
      {

        xt = Mti[0][0] * x + xty;
        yt = Mti[1][0] * x + yty;
        zt = Mti[2][0] * x + zty;
        xs = Msi[0][0] * x + xsy;
        ys = Msi[1][0] * x + ysy;
        zs = Msi[2][0] * x + zsy;

        for (f = 0; f < dt[3]; f++)
        {
          MRIsampleVolumeFrame(mriS, xs, ys, zs, f, &vs);
          if (vs == -1) continue;
          MRIsampleVolumeFrame(mriT, xt, yt, zt, f, &vt);
          if (vt == -1) continue;
        
          dd = vs-vt;
          //d += rhoTukeyBiweight(dd,sat);  
          diff[zpos+count] = (float) dd;
          count++;
        }
      }
    }
  }
  
 // float sigma = RobustGaussian<float>::mad(diff, n);
  //if (sigma == 0.0)
  //sigma = 1.4826;
  //cout << "sigma: " << sigma << endl;
  //if (sigma < 1.0) sigma = 1;
   
  double d = 0.0;
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:d)
#endif
  for (unsigned int i = 0; i < n; i++)
  {
   // d += rhoTukeyBiweight(diff[i] / sigma, sat);
    d += rhoTukeyBiweight(diff[i] , sat);
  }
  free(diff);
  return d;

}


double CostFunctions::normalizedCorrelation(MRI * i1, MRI * i2)
{

  assert(i1->width == i2->width);
  assert(i1->height == i2->height);
  assert(i1->depth == i2->depth);

  double d = 0;
  float d1 = 0;
  float d2 = 0;
  double dd1 = 0;
  double dd2 = 0;

  if (i1->ischunked && i2->ischunked)
  {
    MRIiterator it1(i1);
    MRIiterator it2(i2);
    for (it1.begin(); !it1.isEnd(); it1++)
    {
      d1 = (*it1);
      d2 = (*it2);
      d += d1 * d2;
      dd1 += d1 * d1;
      dd2 += d2 * d2;
      it2++;
    }
  }
  else
  {
    for (int z = 0; z < i1->depth; z++)
      for (int y = 0; y < i1->height; y++)
        for (int x = 0; x < i1->width; x++)
        {
          d1 = MRIgetVoxVal(i1, x, y, z, 0);
          d2 = MRIgetVoxVal(i2, x, y, z, 0);
          d += d1 * d2;
          dd1 += d1 * d1;
          dd2 += d2 * d2;
        }
  }
  return (d / (sqrt(dd1) * sqrt(dd2)));
}


double CostFunctions::SBCost(MRI * mriS, MRI * mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti,
    int d1, int d2, int d3)
{

  // matlab:
  //N = length(vJ);
  //vI = vI - mean(vI); vI = vI / norm(vI);
  //vJ = vJ - mean(vJ); vJ = vJ / norm(vJ);
  //[~, ind] = sort(vI + sign(vI'*vJ)*vJ, 'descend');
  //Ic = cumsum(vI(ind));
  //Jc = cumsum(vJ(ind));
  //cf = -max( (Ic(1:end-1).^2 + Jc(1:end-1).^2) ./ ((1:(N-1)).*((N-1):-1:1))');

  // assert same sizes
  if (mriS->width != mriT->width || mriS->height != mriT->height || mriS->depth != mriT->depth)
  {
    cerr << "CostFunctions::SBCost image dimensions differ!" << endl;
    exit(1);
  }
  
  // currently works only on single frame
  if (mriS->nframes >1 || mriT->nframes>1)
  {
    cerr << "CostFunctions::SBCost only works on single frame inputs!" << endl;
    exit(1);
  }
  
  // map and subsample
  MRI * nmriS = mapMRI(mriS,Msi,d1,d2,d3);
  MRI * nmriT = mapMRI(mriT,Mti,d1,d2,d3);

  // get total size
  unsigned int n  =  nmriS->width * nmriS->height * nmriS->depth;

  // get means
  double m1 = mean(nmriS,0);
  double m2 = mean(nmriT,0);
  
  // convert to stl vector, de-mean and compute norm
  unsigned int counter = 0;
  std::vector < double > iv1(n), iv2(n);
  double norm1 = 0.0;
  double norm2 = 0.0;
  double f1,f2 = 0;
  for (int z = 0; z < nmriS->depth; z++)
    for (int y = 0; y < nmriS->height; y++)
      for (int x = 0; x < nmriS->width; x++)
      {
        f1 = MRIgetVoxVal(nmriS, x, y, z, 0)-m1;
        f2 = MRIgetVoxVal(nmriT, x, y, z, 0)-m2;
        iv1[counter] = f1;
        iv2[counter] = f2;
        norm1 += f1*f1;
        norm2 += f2*f2;
        counter++;
      }
  norm1 = sqrt(norm1);
  norm2 = sqrt(norm2);
      
  // free mem (mapped images not needed as we have std::vector now)
  MRIfree(&nmriS);
  MRIfree(&nmriT);

  // normalize and compute sign
  double prod = 0.0;
  for (unsigned int i = 0; i<n; i++)
  {
    iv1[i] /= norm1;
    iv2[i] /= norm2;
    prod += iv1[i] * iv2[i];    
  }
  int sign = std::copysign(1,prod);

  // compute sorting vector (into iv1 to save space)
  for (unsigned int i = 0; i<n; i++)
  {
    iv1[i] += sign * iv2[i];
  }
  
  // sorting
  // initialize original index locations
  std::vector<size_t> idx(iv1.size());
  iota(idx.begin(), idx.end(), 0);
  // sort indexes based on comparing values in iv1
  // using std::stable_sort to avoid index re-orderings
  // when v contains elements of equal values 
  // stable_sort(idx.begin(), idx.end(),
  sort(idx.begin(), idx.end(),
       [&iv1](size_t i1, size_t i2) {return iv1[i1] > iv1[i2];});

  // cumsums
  std::vector < double > c1(n), c2(n);
  f1=0.0; 
  f2=0.0;
  for (unsigned int i = 0; i<n; i++)
  {
    f1 += iv1[idx[i]];
    f2 += iv2[idx[i]];
    c1[i] = f1;
    c2[i] = f2;
  }
  
  // max
  double mm = 0.0; // max will be positive
  for (unsigned int i = 0; i<n-1; i++)
  {
    f1 = (c1[i]*c1[i] + c2[i]*c2[i]) / ((i+1.0) * (n-1.0-i));
    if (f1 > mm) mm = f1;
  }
  return -mm*n;
}


double CostFunctions::moment(MRI *i, int x, int y, int z)
{
  double dd = 0.0;
  for (int d = 0; d < i->depth; d++)
    for (int h = 0; h < i->height; h++)
      for (int w = 0; w < i->width; w++)
      {
        dd += pow(d + 1.0, z) * pow(h + 1.0, y) * pow(w + 1.0, x)
            * MRIgetVoxVal(i, (int) w, (int) h, (int) d, 0);
      }
  return dd;
}

std::vector<double> CostFunctions::centroid(MRI *i)
// M_100/M_000 , M_010/M_000 , M_001 / M_000
// now ignore outside_vals in centroid computation (for white backgrounds)
{
  //cout << "CostFunctions::centroid" << endl;
  std::vector<double> dd(3, 0.0);
  double n = 0;
  double val;
  double eps = i->outside_val/255.0;
  for (int d = 0; d < i->depth; d++)
    for (int h = 0; h < i->height; h++)
      for (int w = 0; w < i->width; w++)
      {
        val = MRIgetVoxVal(i, w, h, d, 0);
        if (fabs(val-i->outside_val) < eps) val = 0.0; // invisible
        n += val;
        dd[0] += (w + 1) * val;
        dd[1] += (h + 1) * val;
        dd[2] += (d + 1) * val;
      }
  dd[0] = dd[0] / n;
  dd[1] = dd[1] / n;
  dd[2] = dd[2] / n;
  if (std::isnan(dd[0] + dd[1] + dd[2]))
  {
    cerr << "CostFunctions::centroid is NAN (empty image? n = " << n << " )"
        << endl;
    exit(1);
  }

  return dd;
}

vnl_matrix_fixed<double, 3, 3> CostFunctions::orientation(MRI *i)
// M_ijk below is the moment(i,j,k)
{
  // compute mean
  std::vector<double> dd(3, 0.0);
  double n = 0;
  float val;
  int wp1, hp1, dp1;

//  MATRIX* cov = MatrixAlloc(3,3,MATRIX_REAL);
//  cov = MatrixZero(3,3,cov);
  vnl_matrix_fixed<double, 3, 3> cov(0.0);

  for (int d = 0; d < i->depth; d++)
  {
    dp1 = d + 1;
    for (int h = 0; h < i->height; h++)
    {
      hp1 = h + 1;
      for (int w = 0; w < i->width; w++)
      {
        wp1 = w + 1;
        val = MRIgetVoxVal(i, w, h, d, 0);

        n += val; // M_000
        dd[0] += (wp1) * val; // M_100
        dd[1] += (hp1) * val; // M_010
        dd[2] += (dp1) * val; // M_001

        cov[0][0] += wp1 * wp1 * val; // M_200
        cov[1][1] += hp1 * hp1 * val; // M_020
        cov[2][2] += dp1 * dp1 * val; // M_002
        cov[0][1] += wp1 * hp1 * val; // M_110
        cov[1][2] += hp1 * dp1 * val; // M_011
        cov[0][2] += wp1 * dp1 * val; // M_101

      }
    }
  }
  // compute means:
  dd[0] = dd[0] / n;
  dd[1] = dd[1] / n;
  dd[2] = dd[2] / n;

  // compute cov matrix
  cov[0][0] = (cov[0][0] / n) - dd[0] * dd[0]; // M_200/M_000 - mean_x^2
  cov[1][1] = (cov[1][1] / n) - dd[1] * dd[1];
  cov[2][2] = (cov[2][2] / n) - dd[2] * dd[2];
  cov[0][1] = (cov[0][1] / n) - dd[0] * dd[1]; // M_110/M_000 - mean_x * mean_y
  cov[1][2] = (cov[1][2] / n) - dd[1] * dd[2];
  cov[0][2] = (cov[0][2] / n) - dd[0] * dd[2];

  // mirror off diagonal elements: cov is symmetric
  cov[1][0] = cov[0][1];
  cov[2][1] = cov[1][2];
  cov[2][0] = cov[0][2];

  // compute Eigenvectors
  vnl_symmetric_eigensystem<double> SymEig(cov);
  // sort:
  unsigned int smallest = 0;
  unsigned int largest = 0;
  for (uint i = 1; i < 3; i++)
  {
    if (SymEig.D[largest] < SymEig.D[i])
      largest = i;
    if (SymEig.D[smallest] > SymEig.D[i])
      smallest = i;
  }
  unsigned int largetosmall[3];
  largetosmall[0] = largest;
  largetosmall[1] = 3 - smallest - largest;
  largetosmall[2] = smallest;
  // should be sorted when using this lookup:
  assert(SymEig.D[largetosmall[0]] >= SymEig.D[largetosmall[1]]);
  assert(SymEig.D[largetosmall[1]] >= SymEig.D[largetosmall[2]]);
  vnl_matrix_fixed<double, 3, 3> evec;
  for (uint i = 0; i < 3; i++)
    evec.set_column(i, SymEig.V.get_column(largetosmall[i]));

//   float eval[3];
//   MATRIX* evec = MatrixEigenSystem(cov, eval ,NULL) ;
//   // should be sorted:
//   assert(eval[0] >= eval[1]);
//   assert(eval[1] >= eval[2]);

  // make det positive:
  //double d = MatrixDeterminant(evec);
  double d = vnl_det(evec);
  //vnl_matlab_print(vcl_cerr,evec,"evec",vnl_matlab_print_format_long);
  //cout << " det = " << d << endl;
  if (d < 0)
  {
    //cout << "Orientation: neg. determinant ..  fixing" << endl;
    for (int r = 0; r < 3; r++)
      evec[r][0] = -evec[r][0];
    //vnl_matlab_print(vcl_cerr,evec,"evec2",vnl_matlab_print_format_long);
  }

  //cout << " evals: " << eval[0] << " " << eval[1] << " " << eval[2] << endl;
  //cout << " evecs: " << endl;
  // MatrixPrintFmt(stdout,"% 2.8f",evec);

  //MatrixFree(&cov);
  return evec;

}
