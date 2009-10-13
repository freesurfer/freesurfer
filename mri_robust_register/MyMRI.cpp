/**
 * @file  MyMRI.cpp
 * @brief A class for MRI utils
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/10/13 20:08:00 $
 *    $Revision: 1.3 $
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
#include <iostream>
#include "MyMRI.h"

#ifdef __cplusplus
extern "C"
{
#endif
#include "limits.h"
#include "error.h"
#include "macros.h"
#include "mrimorph.h"

#ifdef __cplusplus
}
#endif

using namespace std;

// mainly extraced from mri_convert
MRI* MyMRI::makeConform(MRI *mri, MRI *out, bool fixvoxel, bool fixtype)
{

  out = MRIcopy(mri,out);

  if (mri->type == MRI_UCHAR && mri->xsize == 1 && mri->ysize == 1 && mri->zsize ==1
      && mri->thick == 1 && mri->ps == 1) return out;


  double conform_size = 1;
  int conform_width = findRightSize(mri, conform_size);

  MRI * temp = MRIallocHeader(mri->width, mri->height, mri->depth, mri->type);
  MRIcopyHeader(mri, temp);
  temp->width = temp->height = temp->depth = conform_width;
  temp->imnr0 = 1;
  temp->imnr1 = conform_width;
  temp->type = MRI_UCHAR;
  temp->thick = conform_size;
  temp->ps = conform_size;
  temp->xsize = temp->ysize = temp->zsize = conform_size;

  if (fixvoxel)
  {
    cout << "Making input conform to 1mm voxels" << endl;
    printf("Original Data has (%g, %g, %g) mm size and (%d, %d, %d) voxels.\n",
           mri->xsize, mri->ysize, mri->zsize,
           mri->width, mri->height, mri->depth);
    printf("Data is conformed to %g mm size and %d voxels for all directions\n",
           conform_size, conform_width);
  }
  temp->xstart = temp->ystart = temp->zstart = - conform_width/2;
  temp->xend = temp->yend = temp->zend = conform_width/2;
  temp->x_r = -1.0;
  temp->x_a =  0.0;
  temp->x_s =  0.0;
  temp->y_r =  0.0;
  temp->y_a =  0.0;
  temp->y_s = -1.0;
  temp->z_r =  0.0;
  temp->z_a =  1.0;
  temp->z_s =  0.0;

  /* ----- change type if necessary ----- */
  int no_scale_flag = FALSE;
  if (mri->type != temp->type && fixtype)
  {
    printf("changing data type from %d to %d (noscale = %d)...\n",
           mri->type,temp->type,no_scale_flag);
    MRI * mri2  = MRISeqchangeType(out, temp->type, 0.0, 0.999, no_scale_flag);
    if (mri2 == NULL)
    {
      printf("ERROR: MRISeqchangeType\n");
      exit(1);
    }
    MRIfree(&out);
    out = mri2;
  }

  /* ----- reslice if necessary ----- */
  if ((mri->xsize != temp->xsize ||
       mri->ysize != temp->ysize ||
       mri->zsize != temp->zsize ||
       mri->width != temp->width ||
       mri->height != temp->height ||
       mri->depth != temp->depth ||
       mri->x_r != temp->x_r ||
       mri->x_a != temp->x_a ||
       mri->x_s != temp->x_s ||
       mri->y_r != temp->y_r ||
       mri->y_a != temp->y_a ||
       mri->y_s != temp->y_s ||
       mri->z_r != temp->z_r ||
       mri->z_a != temp->z_a ||
       mri->z_s != temp->z_s ||
       mri->c_r != temp->c_r ||
       mri->c_a != temp->c_a ||
       mri->c_s != temp->c_s) && fixvoxel)
  {
    printf("Reslicing using ");

    int resample_type_val = SAMPLE_TRILINEAR;

    switch (resample_type_val)
    {
    case SAMPLE_TRILINEAR:
      printf("trilinear interpolation \n");
      break;
    case SAMPLE_NEAREST:
      printf("nearest \n");
      break;
    case SAMPLE_SINC:
      printf("sinc \n");
      break;
    case SAMPLE_CUBIC:
      printf("cubic \n");
      break;
    case SAMPLE_WEIGHTED:
      printf("weighted \n");
      break;
    }
    MRI * mri2 = MRIresample(out, temp, resample_type_val);
    if (mri2 == NULL)
    {
      cerr << "makeConform: MRIresample did not return MRI" << endl;
      exit(1);
    }
    MRIfree(&out);
    out = mri2;
  }

  MRIfree(&temp);
  return out;


}


MRI * MyMRI::MRIvalscale(MRI *mri_src, MRI *mri_dst, double s)
// recently found also: MRIscalarMul in mri.h (but has no clipping for short?)
{

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  int      width, height, depth, x, y, z ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  short    *ps_src, *ps_dst ;
  BUFTYPE  *pb_src, *pb_dst ;
  float    *pf_src, *pf_dst, val ;

  switch (mri_src->type)
  {
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pf_src = &MRIFvox(mri_src, 0, y, z) ;
        pf_dst = &MRIFvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = *pf_src++ ;
          val *= s;
          *pf_dst++ = val ;
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        ps_src = &MRISvox(mri_src, 0, y, z) ;
        ps_dst = &MRISvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = (float)(*ps_src++) ;
          val *= s;
          if (val < SHRT_MIN) val = SHRT_MIN;
          if (val > SHRT_MAX) val = SHRT_MAX;
          *ps_dst++ = (short)nint(val) ;
        }
      }
    }
    break ;
  case MRI_UCHAR:
    assert(s > 0);
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pb_src = &MRIvox(mri_src, 0, y, z) ;
        pb_dst = &MRIvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = (float)*pb_src++ ;
          val *= s;
          if (val > 255) val = 255;
          *pb_dst++ = (BUFTYPE)nint(val) ;
        }
      }
    }
    break ;
  default:
    ErrorReturn(mri_dst,
                (ERROR_UNSUPPORTED, "MRIvalScale: unsupported type %d",
                 mri_src->type)) ;
  }

  return(mri_dst) ;

}
MRI * MyMRI::convolute(MRI * mri, MRI * filter, int dir)
// dir 1,2,3  : x,y,z
// filter should be dimension (x,1,1) with odd length x
{
  assert(filter->height ==1 && filter->depth ==1);

  int d[3];
  d[0] = mri->width;
  d[1] = mri->height;
  d[2] = mri->depth;
  //cout << " sizeorig: " << d[0] << " " << d[1] << " " << d[2] << endl;
  int dm1 = dir-1;
  d[dm1] = d[dm1] - filter->width + 1;
  //cout << " sizetarget: " << d[0] << " " << d[1] << " " << d[2] << endl;
  MRI * result = MRIalloc(d[0], d[1], d[2], MRI_FLOAT);
  if (result == NULL) 
        ErrorExit(ERROR_NO_MEMORY,"MyMRI::convolute could not allocate memory for result") ;
  
  MRIclear(result);
  int dd,ff,a,b;
  int ip, frev;
  for (dd = 0; dd < d[dm1]; dd++)
    for (ff = 0; ff < filter->width; ff++)
    {
      ip = dd + ff ;
      frev = filter->width - 1 - ff;
      for ( a = 0; a<d[dir%3]; a++)
        for ( b = 0; b<d[(dir+1)%3]; b++)
        {
          if (mri->type == MRI_FLOAT)
          {
            //cout << " working on MRI_float " << endl;
            if (dir == 1)
            {
              //cout << " reslt( " << dd << " " << a << " " << b << " )   frev= " << frev << " ip: " << ip << endl;
              MRIFvox(result, dd, a, b) += MRIFvox(filter, frev, 0, 0) * MRIFvox(mri, ip, a, b);
            }
            else if (dir ==2)
              MRIFvox(result, b, dd, a) += MRIFvox(filter, frev, 0, 0) * MRIFvox(mri, b, ip, a);
            else if (dir == 3)
              MRIFvox(result, a, b, dd) += MRIFvox(filter, frev, 0, 0) * MRIFvox(mri, a, b, ip);
            else assert(dir > 0 && dir < 4);
          }
          else if (mri->type == MRI_UCHAR || mri->type == MRI_INT || mri->type == MRI_SHORT || mri->type == MRI_LONG)
          {
//          if (dir == 1)
//   {
//      //cout << " reslt( " << dd << " " << a << " " << b << " )   frev= " << frev << " ip: " << ip << endl;
//      MRIFvox(result, dd, a, b) += MRIFvox(filter, frev, 0, 0) * (int)MRIvox(mri, ip, a, b);
//   }
//   else if (dir ==2)
//      MRIFvox(result, b, dd, a) += MRIFvox(filter, frev, 0, 0) * (int)MRIvox(mri, b, ip, a);
//   else if (dir == 3)
//      MRIFvox(result, a, b, dd) += MRIFvox(filter, frev, 0, 0) * (int)MRIvox(mri, a, b, ip);
//   else assert(dir > 0 && dir < 4);

            //cout << " reslt( " << dd << " " << a << " " << b << " )   frev= " << frev << " ip: " << ip << endl;
            if (dir == 1)
            {
              //cout << " mri " <<  MRIgetVoxVal(mri, ip, a, b,0) << endl;
              MRIFvox(result, dd, a, b) += MRIFvox(filter, frev, 0, 0) * MRIgetVoxVal(mri, ip, a, b,0);
            }
            else if (dir ==2)
              MRIFvox(result, b, dd, a) += MRIFvox(filter, frev, 0, 0) * MRIgetVoxVal(mri, b, ip, a,0);
            else if (dir == 3)
              MRIFvox(result, a, b, dd) += MRIFvox(filter, frev, 0, 0) * MRIgetVoxVal(mri, a, b, ip,0);
            else assert(dir > 0 && dir < 4);

          }
          else // cannot deal with type
            assert(1==2);
        }
    }


  return result;

}

MRI * MyMRI::getPrefilter()
{
  MRI *mri_prefilter ;
  mri_prefilter = MRIalloc(5,1,1, MRI_FLOAT);
  MRIFvox(mri_prefilter, 0, 0, 0) =  0.03504 ;
  MRIFvox(mri_prefilter, 1, 0, 0) =  0.24878 ;
  MRIFvox(mri_prefilter, 2, 0, 0) =  0.43234 ;
  MRIFvox(mri_prefilter, 3, 0, 0) =  0.24878 ;
  MRIFvox(mri_prefilter, 4, 0, 0) =  0.03504 ;

  return mri_prefilter;
}

MRI * MyMRI::getDerfilter()
{
  MRI *mri_derfilter ;
  mri_derfilter = MRIalloc(5,1,1, MRI_FLOAT);
  MRIFvox(mri_derfilter, 0, 0, 0) =  0.10689 ;
  MRIFvox(mri_derfilter, 1, 0, 0) =  0.28461 ;
  MRIFvox(mri_derfilter, 2, 0, 0) =  0.0 ;
  MRIFvox(mri_derfilter, 3, 0, 0) =  -0.28461 ;
  MRIFvox(mri_derfilter, 4, 0, 0) =  -0.10689 ;
  return mri_derfilter;
}

MRI * MyMRI::subSample(MRI * mri)
{
  int w = (mri->width +1) / 2;
  int h = (mri->height+1) / 2;
  int d = (mri->depth +1) / 2;
//   int d = mri->depth;
  MRI* mri_sub = MRIalloc(w,h,d,mri->type);
  int x,y,z;
  for (z = 0;z<d;z++)
    for (y = 0;y<h;y++)
      for (x = 0;x<w;x++)
        MRIsetVoxVal(mri_sub,x,y,z,0,MRIgetVoxVal(mri,2*x,2*y,2*z,0));

  return mri_sub;

}

MRI * MyMRI::getBlur(MRI* mriS)
{
  MRI *mri_prefilter = MyMRI::getPrefilter();
  MRI *tmp1 = MyMRI::convolute(mriS,mri_prefilter,1);
  MRI *tmp2 = MyMRI::convolute(tmp1,mri_prefilter,2);
  MRI *tmp3 = MyMRI::convolute(tmp2,mri_prefilter,3);
  MRIfree(&tmp1);
  MRIfree(&tmp2);
  MRIfree(&mri_prefilter);
  return tmp3;
}

MRI * MyMRI::getPartial(MRI* mriS, int dir)
// dir 1,2,3  = x,y,z
{
  assert(dir > 0);
  assert(dir < 4);

  // construct convolution masks:
  MRI *mri_prefilter = getPrefilter();
  MRI *mri_derfilter = getDerfilter();

  // convolute with derivative filter in dir axis
  // and with prefilter along the other two axis

  //int whd[3] = {MRI_WIDTH ,MRI_HEIGHT, MRI_DEPTH};
  MRI* mtmp = MRIcopyFrame(mriS, NULL, 0, 0) ;
  MRI* mtmp2;
  //int klen = mri_prefilter->width ;
  for (int i =1;i<=3; i++)
  {
    if (i==dir)
    {
      //MRIconvolve1d(mtmp, mtmp, &MRIFvox(mri_derfilter, 0, 0, 0), klen, whd[i-1], 0, 0) ;
      mtmp2 = convolute(mtmp,mri_derfilter,i);
      MRIfree(&mtmp);
      mtmp = mtmp2;
    }
    else
    {
      //MRIconvolve1d(mtmp, mtmp, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[i-1], 0, 0) ;
      mtmp2 = convolute(mtmp,mri_prefilter,i);
      MRIfree(&mtmp);
      mtmp = mtmp2;
    }
  }
//  MRI *mri_dst = MRIclone(mriS, NULL) ;
//  MRIcopyFrame(mtmp, mri_dst, 0, 0) ;    /* convert it back to UCHAR */
//  MRIfree(&mtmp);
//  return mri_dst;
  MRIfree(&mri_prefilter);
  MRIfree(&mri_derfilter);
  return mtmp;
}

bool MyMRI::getPartials(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz, MRI* &outblur)
{

  assert(outfx == NULL && outfy == NULL && outfz == NULL && outblur == NULL );

  // construct convolution masks:
  MRI *mri_prefilter = getPrefilter();
  MRI *mri_derfilter = getDerfilter();

  MRI* mdz   = convolute(mri,mri_derfilter,3);
  MRI* mbz   = convolute(mri,mri_prefilter,3);

  MRI* mdzby = convolute(mdz,mri_prefilter,2);
  MRI* mbzby = convolute(mbz,mri_prefilter,2);
  MRI* mbzdy = convolute(mbz,mri_derfilter,2);
  MRIfree(&mdz);
  MRIfree(&mbz);

  outfx = convolute(mbzby,mri_derfilter,1);
  outfy = convolute(mbzdy,mri_prefilter,1);
  outfz = convolute(mdzby,mri_prefilter,1);
  outblur = convolute(mbzby,mri_prefilter,1);

  //cout << " size fx: " << outfx->width << " " << outfx->height << " " << outfx->depth << endl;

  MRIfree(&mdzby);
  MRIfree(&mbzby);
  MRIfree(&mbzdy);

  MRIfree(&mri_prefilter);
  MRIfree(&mri_derfilter);

  return true;
}

MRI * MyMRI::getBlur2(MRI* mri)
{
  MRI* outblur = MRIgaussianSmooth(mri, 1, 1,NULL);
//  mri_kernel = MRIgaussian1d(1, -1) ;
// MRI* outblur = MRIconvolveGaussian(mri, NULL, mri_kernel);
//  MRIfree(&mri_kernel);
  return outblur;
}

bool MyMRI::getPartials2(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz, MRI* &outblur)
{

  assert(outfx == NULL && outfy == NULL && outfz == NULL && outblur == NULL );

  outblur = getBlur2(mri);
  outfx   = MRIxDerivative(outblur,NULL);
  outfy   = MRIyDerivative(outblur,NULL);
  outfz   = MRIzDerivative(outblur,NULL);

  return true;
}



// this function is called when conform is done
// copied from mri_convert
int MyMRI::findRightSize(MRI *mri, float conform_size)
{
  // user gave the conform_size
  double xsize, ysize, zsize;
  double fwidth, fheight, fdepth, fmax;
  int conform_width;

  xsize = mri->xsize;
  ysize = mri->ysize;
  zsize = mri->zsize;

  // now decide the conformed_width
  // calculate the size in mm for all three directions
  fwidth = mri->xsize*mri->width;
  fheight = mri->ysize*mri->height;
  fdepth = mri->zsize*mri->depth;
  // pick the largest
  if (fwidth> fheight)
    fmax = (fwidth > fdepth) ? fwidth : fdepth;
  else
    fmax = (fdepth > fheight) ? fdepth : fheight;
  // get the width with conform_size
  conform_width = (int) ceil(fmax/conform_size);

  // just to make sure that if smaller than 256, use 256 anyway
  if (conform_width < 256)
    conform_width = 256;
  // conform_width >= 256.   allow 10% leeway
  else if ((conform_width -256.)/256. < 0.1)
    conform_width = 256;

  // if more than 256, warn users
  if (conform_width > 256)
  {
    fprintf(stderr, "WARNING =================="
            "++++++++++++++++++++++++"
            "=======================================\n");
    fprintf(stderr, "The physical sizes are "
            "(%.2f mm, %.2f mm, %.2f mm), "
            "which cannot fit in 256^3 mm^3 volume.\n",
            fwidth, fheight, fdepth);
    fprintf(stderr, "The resulting volume will have %d slices.\n",
            conform_width);
    fprintf(stderr, "If you find problems, please let us know "
            "(freesurfer@nmr.mgh.harvard.edu).\n");
    fprintf(stderr, "=================================================="
            "++++++++++++++++++++++++"
            "===============\n\n");
  }
  return conform_width;
}

MATRIX* MyMRI::MRIgetZslice(MRI * mri, int slice)
// extract a z-slice from mri
// return as matrix
{
  assert(slice >=0 && slice < mri->depth);
  MATRIX* m = MatrixAlloc(mri->height,mri->width,MATRIX_REAL);
  int x,y;
  for (y = 0 ; y < mri->height ; y++)
    for (x = 0 ; x < mri->width  ; x++)
    {
      if (mri->type == MRI_FLOAT)
        *MATRIX_RELT(m, y+1, x+1) = MRIFvox(mri, x, y, slice);
      else if (mri->type == MRI_UCHAR || mri->type == MRI_INT)
        *MATRIX_RELT(m, y+1, x+1) = (int)MRIvox(mri,x,y,slice);
      else (assert (1==2));
    }

  return m;
}

