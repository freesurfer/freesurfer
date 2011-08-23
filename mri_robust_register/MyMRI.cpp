/**
 * @file  MyMRI.cpp
 * @brief A class for MRI utils
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2011/08/23 18:53:40 $
 *    $Revision: 1.8.2.1 $
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
#include <cassert>
#include <iostream>
#include "MyMRI.h"
#include "MyMatrix.h"

#ifdef __cplusplus
extern "C"
{
#endif
#include "limits.h"
#include "error.h"
#include "macros.h"
#include "mrimorph.h"
#include "histo.h"


#ifdef __cplusplus
}
#endif

using namespace std;


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
//  MRIFvox(mri_derfilter, 0, 0, 0) =  0.10689 ;
//  MRIFvox(mri_derfilter, 1, 0, 0) =  0.28461 ;
//  MRIFvox(mri_derfilter, 2, 0, 0) =  0.0 ;
//  MRIFvox(mri_derfilter, 3, 0, 0) =  -0.28461 ;
//  MRIFvox(mri_derfilter, 4, 0, 0) =  -0.10689 ;
  MRIFvox(mri_derfilter, 0, 0, 0) =  -0.10689 ;
  MRIFvox(mri_derfilter, 1, 0, 0) =  -0.28461 ;
  MRIFvox(mri_derfilter, 2, 0, 0) =   0.0 ;
  MRIFvox(mri_derfilter, 3, 0, 0) =   0.28461 ;
  MRIFvox(mri_derfilter, 4, 0, 0) =   0.10689 ;
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

MRI * MyMRI::getBlur(MRI* mriS, MRI* mriT)
{
  MRI *mri_prefilter = MyMRI::getPrefilter();
  mriT = MRIconvolveGaussian(mriS,mriT,mri_prefilter);
  MRIfree(&mri_prefilter);
  return mriT;
//   MRI *mri_prefilter = MyMRI::getPrefilter();
//   MRI *tmp1 = MyMRI::convolute(mriS,mri_prefilter,1);
//   MRI *tmp2 = MyMRI::convolute(tmp1,mri_prefilter,2);
//   MRI *tmp3 = MyMRI::convolute(tmp2,mri_prefilter,3);
//   MRIfree(&tmp1);
//   MRIfree(&tmp2);
//   MRIfree(&mri_prefilter);
//   return tmp3;
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

  //double size = mri->width *mri->height * mri->depth *sizeof(float) / (1024.0 * 1024.0);
  //cout.precision(4);
  //cout << " MyMRI::getPartials allocation  max: " << 5 * size << " Mb  return = " << 4*size << " Mb" << endl;

  // construct convolution masks:
  MRI *mri_prefilter = getPrefilter();
  MRI *mri_derfilter = getDerfilter();
  int klen = mri_prefilter->width ;
  int whd[3] = {MRI_WIDTH ,MRI_HEIGHT, MRI_DEPTH};

  MRI* mdz =  MRIconvolve1d(mri, NULL, &MRIFvox(mri_derfilter, 0, 0, 0), klen, whd[3-1], 0, 0) ;
  //MRI* mdz   = convolute(mri,mri_derfilter,3);

  MRI* mbz =  MRIconvolve1d(mri, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[3-1], 0, 0) ;
  //MRI* mbz   = convolute(mri,mri_prefilter,3);

  MRI* mdzby = MRIconvolve1d(mdz, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[2-1], 0, 0) ;
  //MRI* mdzby = convolute(mdz,mri_prefilter,2);
  
  MRIfree(&mdz);
  outfz = MRIconvolve1d(mdzby, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[1-1], 0, 0) ;
  //outfz = convolute(mdzby,mri_prefilter,1);
  MRIfree(&mdzby);
  
  MRI* mbzby = MRIconvolve1d(mbz, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[2-1], 0, 0) ;
  //MRI* mbzby = convolute(mbz,mri_prefilter,2);
  MRI* mbzdy = MRIconvolve1d(mbz, NULL, &MRIFvox(mri_derfilter, 0, 0, 0), klen, whd[2-1], 0, 0) ;
  //MRI* mbzdy = convolute(mbz,mri_derfilter,2);
  MRIfree(&mbz);
  outfy = MRIconvolve1d(mbzdy, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[1-1], 0, 0) ;
  //outfy = convolute(mbzdy,mri_prefilter,1);
  MRIfree(&mbzdy);
  
  outfx = MRIconvolve1d(mbzby, NULL, &MRIFvox(mri_derfilter, 0, 0, 0), klen, whd[1-1], 0, 0) ;
  //outfx = convolute(mbzby,mri_derfilter,1);
  outblur = MRIconvolve1d(mbzby, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[1-1], 0, 0) ;
  //outblur = convolute(mbzby,mri_prefilter,1);
  MRIfree(&mbzby);

  //cout << " size fx: " << outfx->width << " " << outfx->height << " " << outfx->depth << endl;

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



std::vector < int >  MyMRI::findRightSize(MRI *mri, float conform_size, bool conform)
// determines width, height, depth to cover mri with conform_size isotropic voxels
//  bool conform makes cube image (w=h=d) and adjust to min of 256
{

  double xsize, ysize, zsize;
  double fwidth, fheight, fdepth, fmax;
  int conform_width;

  xsize = mri->xsize;
  ysize = mri->ysize;
  zsize = mri->zsize;

  // now decide the conformed_width
  // calculate the size in mm for all three directions
  fwidth  = mri->xsize*mri->width;
  fheight = mri->ysize*mri->height;
  fdepth  = mri->zsize*mri->depth;
  
  vector < int > ret(3);
  double eps = 0.0001; // to prevent ceil(2.0*64 / 2.0) = ceil(64.000000000001) = 65
  ret[0] = (int) ceil((fwidth/conform_size)-eps);
  ret[1] = (int) ceil((fheight/conform_size)-eps);
  ret[2] = (int) ceil((fdepth/conform_size)-eps);
  
  //cout << " zsize: " << zsize << " depth: " << mri->depth << " fdepth: " << fdepth << " new depth: " << ret[2] << endl;
  
  if (! conform) return ret;
  
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

//   // if more than 256, warn users
//   if (conform_width > 256)
//   {
// //     fprintf(stderr, "WARNING =================="
// //             "++++++++++++++++++++++++"
// //             "=======================================\n");
// //     fprintf(stderr, "The physical sizes are "
// //             "(%.2f mm, %.2f mm, %.2f mm), "
// //             "which cannot fit in 256^3 mm^3 volume.\n",
// //             fwidth, fheight, fdepth);
// //     fprintf(stderr, "The resulting volume will have %d slices.\n",
// //             conform_width);
// //     fprintf(stderr, "If you find problems, please let us know "
// //             "(freesurfer@nmr.mgh.harvard.edu).\n");
// //     fprintf(stderr, "=================================================="
// //             "++++++++++++++++++++++++"
// //             "===============\n\n");
//       cout << " ... Will not fit into 256^3 volume, will have " << conform_width << " slices" << endl;
//   }

  
  ret[0] = conform_width;
  ret[1] = conform_width;
  ret[2] = conform_width;

  return ret;
}

bool MyMRI::isConform(MRI * mri)
{

  return ( mri->xsize == 1   && mri->ysize == 1   && mri->zsize == 1 &&
           mri->width == 256 && mri->height == 256 && mri->depth == 256 );

}

bool MyMRI::isIsotropic(MRI * mri)
{

  return ( mri->xsize == mri->ysize  && mri->xsize ==  mri->zsize );

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

MRI * MyMRI::MRIlinearTransform(MRI* mriS, MRI* mriT, const vnl_matrix_fixed < double,4,4 >& m)
{
   MATRIX * mm = MyMatrix::convertVNL2MATRIX(m,NULL);
   mriT = ::MRIlinearTransform(mriS,mriT,mm);
   MatrixFree(&mm);
   return mriT;
}

vnl_matrix_fixed < double, 4, 4 > MyMRI::MRIvoxelXformToRasXform(MRI * mri_src, MRI * mri_dst, const vnl_matrix_fixed < double,4,4> &m_vox)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras;

  if (!mri_dst)
    mri_dst = mri_src ;  /* assume they will be in the same space */

  m_ras_to_voxel = MRIgetRasToVoxelXform(mri_src) ;
  m_voxel_to_ras = MRIgetVoxelToRasXform(mri_dst) ;

  //m_tmp = MatrixMultiply(m_voxel_xform, m_ras_to_voxel, NULL) ;
  //m_ras_xform = MatrixMultiply(m_voxel_to_ras, m_tmp, m_ras_xform) ;
  
  vnl_matrix_fixed < double, 4, 4 > m_ras = MyMatrix::convertMATRIX2VNL(m_voxel_to_ras) *
       m_vox * MyMatrix::convertMATRIX2VNL(m_ras_to_voxel);


  MatrixFree(&m_voxel_to_ras);
  MatrixFree(&m_ras_to_voxel);

  return(m_ras) ;


}

MRI * MyMRI::entropyImage(MRI* mri, int radius, int sigma )
{
  int width  = mri->width;
  int height = mri->height;
  int depth  = mri->depth;

  MRI * mriIn = mri;
  if (mri->type != MRI_UCHAR)
  {
    int no_scale_flag = FALSE;
    printf("changing data type from %d to %d (noscale = %d)...\n",
           mri->type,MRI_UCHAR,no_scale_flag);
    mriIn  = MRISeqchangeType(mri, MRI_UCHAR, 0.0, 0.999, no_scale_flag);
    if (mriIn == NULL)
    {
      printf("ERROR: MRISeqchangeType\n");
      exit(1);
    }
   // MRIwrite(mriIn,"mriIn.mgz");
  }

  MRI * entI =  MRIalloc(width,height,depth, MRI_FLOAT);
  entI = MRIcopyHeader(mri,entI);

	//float patchGaussFactor = gaussFactor( patchSigma );
	float minEntropy = 1e8; 
  float maxEntropy = -1.0;
  
  //HISTOGRAM *hpdf = HISTOgaussianPDF(NULL,0,sigma,3*radius*radius);
  int twosig2 = 2*sigma*sigma;
  float gfactor = 1.0/sqrt(M_PI*twosig2);
  
  int x,y,z;
  for (z=0;z<depth;z++)
  {
    if ((z+1)%10 == 0) cout << "#" << flush;
    else  cout << "." << flush;
    for (y=0;y<height;y++)
    {
      for (x=0;x<width;x++)
      {
        // compute neighborhood bounds
        int xMin = x - radius; 
        int xMax = x + radius;
        int yMin = y - radius;
        int yMax = y + radius;
        int zMin = z - radius;
        int zMax = z + radius;
        if (xMin < 0) xMin = 0;
        if (xMax > width - 1) xMax = width - 1;
        if (yMin < 0) yMin = 0;
        if (yMax > height - 1) yMax = height - 1;
        if (zMin < 0) zMin = 0;
        if (zMax > depth - 1) zMax = depth - 1;

        // compute histogram of pixel values in neighborhood
        HISTOGRAM *h = HISTOalloc(64);
        HISTOclear(h,h);
        HISTOinit( h,64,0,252 );
        int count = 0; //, outsideCount = 0;
        for (int zp = zMin; zp <= zMax; zp++)
        {
          for (int yp = yMin; yp <= yMax; yp++)
          {
            for (int xp = xMin; xp <= xMax; xp++)
            {
//            if (mask.data( xp, yp )) {
              int index = MRIvox(mriIn, xp, yp, zp) / 4;
              int dx = xp - x;
              int dy = yp - y;
              int dz = zp - z;
              int distSqd =  (dx * dx + dy * dy + dz * dz);
              //float weight = gauss( distSqd, patchGaussFactor ); // fix(faster): store in table?
             // cout << " min: " << hpdf->min << " max: " << hpdf->max << "  distsqd: " << distSqd << endl;
             // assert(distSqd <= hpdf->max);
             // assert(distSqd >= hpdf->min);
             
              //float weight = hpdf->counts[distSqd];
              float weight = gfactor * exp(- (float)distSqd/twosig2);
              h->counts[ index ] += weight;
              if (index - 1 >= 0)
                h->counts[ index - 1 ] += weight * 0.5f;
              if (index + 1 < 64)
                h->counts[ index + 1 ] += weight * 0.5f;
              count++;
//            } else {
//              outsideCount++;
//            }
            }
          }
        }

        // normalize intensity?

        // if have enough values
        if (count)
        {
  
//          // normalize the histogram
//          float factor = 1.0f / histogram.sum();
//          multiply( histogram, factor, histogram );

          // compute entropy from histogram
          float e = (float) HISTOgetEntropy(h);
          if (e<0 || isnan(e) )
          {
            cout << " count: " << count << "  e: " << e << endl;
            cout << " h bin size: " << h->bin_size << endl;
            cout << " hmin: " << h->min << endl;
            cout << " hmax: " << h->max << endl;
            cout << " total: " << HISTOtotal(h) << endl;
            
            for (int b = 0 ; b < h->nbins ; b++)
              cout << " ["<< b << "] : " << h->bins[b] << "  " << h->counts[b] << endl;
            exit(1);
          }
          MRIsetVoxVal(entI,x,y,z,0,e);
          if (e < minEntropy)
            minEntropy = e;
          if (e > maxEntropy)
            maxEntropy = e;
        }
        else
        {
          assert (count > 0 );
          cout << " count : " << count << endl;
          MRIsetVoxVal(entI,x,y,z,0,0.0);
        
        }
      }
    }
  }

//  // clear unititialzed values to the minimum value
//  for (int y = 0; y < height; y++) {
//    for (int x = 0; x < width; x++) {
//      if (entImage.data( x, y ) == -1)
//        entImage.data( x, y ) = minEntropy;
//    }
//  }
//
//  // display stats
//  float min = 0, mean = 0, max = 0;
//  imageStats( entImage, min, mean, max );
//  disp( 1, "entropy: %f / %f / %f", min, mean, max );
  cout << " Min Entropy: " << minEntropy << "  Max Entropy: " << maxEntropy << endl;

  // transform to [0, 255] range
//  aptr<ImageGrayU> entImageGray = toUChar( entImage );
  
 
  int no_scale_flag = FALSE;
//  printf("changing data type from %d to %d (noscale = %d)...\n",
 //          mri->type,MRI_UCHAR,no_scale_flag);
  MRI* mriEnt = MRISeqchangeType(entI, MRI_UCHAR, 0.0, 0.999, no_scale_flag);
  if (mriEnt == NULL)
  {
    printf("ERROR: MRISeqchangeType\n");
    exit(1);
  }
  //MRIwrite(entI,"mriEntFloat.mgz");
  MRIfree(&entI);
  
  if (mri != mriIn)
    MRIfree(&mriIn);
  
  return mriEnt;

}
