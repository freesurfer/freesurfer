/**
 * @brief A class for MRI utils
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
#include <iostream>
#include <algorithm>
#include <vector>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "MyMRI.h"
#include "MyMatrix.h"
#include "CostFunctions.h"
#include "RobustGaussian.h"

#include "limits.h"
#include "error.h"
#include "macros.h"
#include "mrimorph.h"
#include "histo.h"

using namespace std;

MRI * MyMRI::MRInorm255(MRI *mri_src, MRI *mri_dst)
// normalizes so that min =0 and max=255
{
  std::pair<float, float> mm = CostFunctions::minmax(mri_src);
  return MyMRI::MRIvalscale(mri_src, mri_dst, 255 / (mm.second - mm.first), mm.first);
}

MRI * MyMRI::MRIvalscale(MRI *mri_src, MRI *mri_dst, double s, double b)
// recently found also: MRIscalarMul in mri.h (but has no clipping for short?)
{

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL);

  int width, height, depth, nf, x, y, z, f;
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  nf = mri_src->nframes;
  short *ps_src, *ps_dst;
  BUFTYPE *pb_src, *pb_dst;
  float *pf_src, *pf_dst, val;

  switch (mri_src->type)
  {
  case MRI_FLOAT:
    for (f = 0; f < nf ; f++ )
    for (z = 0; z < depth; z++)
    {
      for (y = 0; y < height; y++)
      {
        pf_src = &MRIFseq_vox(mri_src, 0, y, z, f) ;
        pf_dst = &MRIFseq_vox(mri_dst, 0, y, z, f);
        for (x = 0; x < width; x++)
        {
          val = *pf_src++;
          val -= b;
          val *= s;
          *pf_dst++ = val;
        }
      }
    }
    break;
    case MRI_SHORT:
    for (f = 0; f < nf ; f++ )
    for (z = 0; z < depth; z++)
    {
      for (y = 0; y < height; y++)
      {
        ps_src = &MRISseq_vox(mri_src, 0, y, z, f);
        ps_dst = &MRISseq_vox(mri_dst, 0, y, z, f);
        for (x = 0; x < width; x++)
        {
          val = (float)(*ps_src++);
          val -= b;
          val *= s;
          if (val < SHRT_MIN) val = SHRT_MIN;
          if (val > SHRT_MAX) val = SHRT_MAX;
          *ps_dst++ = (short)nint(val);
        }
      }
    }
    break;
    case MRI_USHRT:
    for (f = 0; f < nf ; f++ )
    for (z = 0; z < depth; z++)
    {
      for (y = 0; y < height; y++)
      {
        unsigned short *ps_src2 = &MRIUSseq_vox(mri_src, 0, y, z, f);
        unsigned short *ps_dst2 = &MRIUSseq_vox(mri_dst, 0, y, z, f);
        for (x = 0; x < width; x++)
        {
          val = (float)(*ps_src2++);
          val -= b;
          val *= s;
          if (val < 0) val = 0;
          if (val > USHRT_MAX) val = USHRT_MAX;
          *ps_dst2++ = (unsigned short)nint(val);
        }
      }
    }
    break;
    case MRI_UCHAR:
    assert(s > 0);
    for (f = 0; f < nf ; f++ )
    for (z = 0; z < depth; z++)
    {
      for (y = 0; y < height; y++)
      {
        pb_src = &MRIseq_vox(mri_src, 0, y, z, f);
        pb_dst = &MRIseq_vox(mri_dst, 0, y, z, f);
        for (x = 0; x < width; x++)
        {
          val = (float)*pb_src++;
          val -= b;
          val *= s;
          if (val > 255) val = 255;
          *pb_dst++ = (BUFTYPE)nint(val);
        }
      }
    }
    break;
    default:
    ErrorReturn(mri_dst,(ERROR_UNSUPPORTED, "MRIvalScale: unsupported type %d",mri_src->type)) ;
    break;
  }

  return (mri_dst);

}

MRI * MyMRI::getPrefilter()
{
  MRI *mri_prefilter;
  mri_prefilter = MRIalloc(5, 1, 1, MRI_FLOAT);
  MRIFvox(mri_prefilter, 0, 0, 0) = 0.03504;
  MRIFvox(mri_prefilter, 1, 0, 0) = 0.24878;
  MRIFvox(mri_prefilter, 2, 0, 0) = 0.43234;
  MRIFvox(mri_prefilter, 3, 0, 0) = 0.24878;
  MRIFvox(mri_prefilter, 4, 0, 0) = 0.03504;

  return mri_prefilter;
}

MRI * MyMRI::getDerfilter()
{
  MRI *mri_derfilter;
  mri_derfilter = MRIalloc(5, 1, 1, MRI_FLOAT);
  MRIFvox(mri_derfilter, 0, 0, 0) = -0.10689;
  MRIFvox(mri_derfilter, 1, 0, 0) = -0.28461;
  MRIFvox(mri_derfilter, 2, 0, 0) = 0.0;
  MRIFvox(mri_derfilter, 3, 0, 0) = 0.28461;
  MRIFvox(mri_derfilter, 4, 0, 0) = 0.10689;
  return mri_derfilter;
}

MRI * MyMRI::subSample(MRI * mri_src, MRI * mri_dst, bool fixheader,
    int randpos)
// if randpos == -1, do not use random locations, else start at randpos;
// this only makes sense when images are smoothed before
// else you need to average neighboring voxel as in MRIdownsample2 (mri.h)
{

//  int w = (mri_src->width +1) / 2;
//  int h = (mri_src->height+1) / 2;
//  int d = (mri_src->depth +1) / 2;
//  if (fixheader)
//  {
//    w = mri_src->width  / 2;
//    h = mri_src->height / 2;
//    d = mri_src->depth  / 2;
//  }

  // for random position within 2x2x2 box, we need each box to exist in image:
  int w = mri_src->width / 2;
  int h = mri_src->height / 2;
  int d = mri_src->depth / 2;

  if (mri_src->depth == 1)
    d = 1; // 2d image

  if (!mri_dst)
  {
    mri_dst = MRIallocSequence(w, h, d, mri_src->type, mri_src->nframes);
    MRIcopyHeader(mri_src, mri_dst);
  }

  int x, y, z, f;
  int dx, dy, dz = 0;
  for (z = 0; z < d; z++)
    for (y = 0; y < h; y++)
      for (x = 0; x < w; x++)
      {
        if (randpos >= 0)
        {
          // random offset 0 or 1:
          dx = (int) (2.0 * MyMRI::getRand(randpos));
          randpos++;
          dy = (int) (2.0 * MyMRI::getRand(randpos));
          randpos++;
          if (d > 1)
          {
            dz = (int) (2.0 * MyMRI::getRand(randpos));
            randpos++;
          }
          for (f = 0; f<mri_src->nframes; f++)
            MRIsetVoxVal(mri_dst, x, y, z, f,
                MRIgetVoxVal(mri_src, 2 * x + dx, 2 * y + dy, 2 * z + dz, f));
        }
        else
          for (f = 0; f<mri_src->nframes; f++)
            MRIsetVoxVal(mri_dst, x, y, z, f,
                MRIgetVoxVal(mri_src, 2 * x, 2 * y, 2 * z, f));
      }

  if (fixheader) // adjusts header of dest so that RAS coordinates agree
  {
    mri_dst->imnr0 = mri_src->imnr0;
    mri_dst->imnr1 = mri_src->imnr0 + mri_dst->depth - 1;
    mri_dst->xsize = mri_src->xsize * 2;
    mri_dst->ysize = mri_src->ysize * 2;
    mri_dst->zsize = mri_src->zsize * 2;
    mri_dst->thick = mri_src->thick * 2;
    mri_dst->ps = mri_src->ps * 2;

    // adjust cras
    //printf("COMPUTING new CRAS\n") ;
    VECTOR* C = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(C,1)= mri_src->width/2+0.5;
    VECTOR_ELT(C,2)= mri_src->height/2+0.5;
    VECTOR_ELT(C,3)= mri_src->depth/2+0.5;
    VECTOR_ELT(C,4)= 1.0;
    MATRIX* V2R = extract_i_to_r(mri_src);
    MATRIX* P = MatrixMultiply(V2R, C, NULL);
    mri_dst->c_r = P->rptr[1][1];
    mri_dst->c_a = P->rptr[2][1];
    mri_dst->c_s = P->rptr[3][1];
    MatrixFree(&P);
    MatrixFree(&V2R);
    VectorFree(&C);

    MRIreInitCache(mri_dst);
  }

  return mri_dst;

}

MRI * MyMRI::getBlur(MRI* mriS, MRI* mriT)
{
  MRI *mri_prefilter = MyMRI::getPrefilter();
  mriT = MRIconvolveGaussian(mriS, mriT, mri_prefilter);
  MRIfree(&mri_prefilter);
  return mriT;
}

bool MyMRI::getPartials(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz,
    MRI* &outblur)
{

  assert(outfx == NULL && outfy == NULL && outfz == NULL && outblur == NULL);

  //double size = mri->width *mri->height * mri->depth *sizeof(float) / (1024.0 * 1024.0);
  //cout.precision(4);
  //cout << " MyMRI::getPartials allocation  max: " << 5 * size << " Mb  return = " << 4*size << " Mb" << endl;

  // construct convolution masks:
  MRI *mri_prefilter = getPrefilter();
  MRI *mri_derfilter = getDerfilter();
  int klen = mri_prefilter->width;
  int whd[3] =
  { MRI_WIDTH, MRI_HEIGHT, MRI_DEPTH };
  bool is2d = (mri->depth == 1);

  outfx = MRIclone(mri,NULL);
  outfy = MRIclone(mri,NULL);
  if (!is2d)
    outfz = MRIclone(mri,NULL);
  else 
    outfz = NULL;
  outblur = MRIclone(mri,NULL);
  for (int f = 0; f < mri->nframes; f++)
  {
    // compute outfx
    MRI* mdx = MRIconvolve1d(mri, NULL, &MRIFvox(mri_derfilter, 0, 0, 0), klen, whd[1-1], f, 0) ;
    MRI* mdxby = MRIconvolve1d(mdx, NULL, &MRIFvox(mri_prefilter, 0, 0, 0) , klen, whd[2-1], 0, 0);
    MRIfree(&mdx);
    if (is2d)
      MRIcopyFrame(mdxby,outfx,0,f);
    else
      MRIconvolve1d(mdxby, outfx, &MRIFvox(mri_prefilter, 0, 0, 0) , klen, whd[3-1], 0, f);
    MRIfree(&mdxby);

    // compute outfy 
    MRI* mbx = MRIconvolve1d(mri, NULL, &MRIFvox(mri_prefilter, 0, 0, 0) , klen, whd[1-1], f, 0);
    MRI* mbxdy = MRIconvolve1d(mbx, NULL, &MRIFvox(mri_derfilter, 0, 0, 0) , klen, whd[2-1], 0, 0);
    MRI* mbxby = MRIconvolve1d(mbx, NULL, &MRIFvox(mri_prefilter, 0, 0, 0) , klen, whd[2-1], 0, 0);
    MRIfree(&mbx);
    if (is2d)
      MRIcopyFrame(mbxdy,outfy,0,f);
    else
      MRIconvolve1d(mbxdy, outfy, &MRIFvox(mri_prefilter, 0, 0, 0) , klen, whd[3-1], 0, f);
    MRIfree(&mbxdy);

    // compute outfz (using mbxby from above)
    if (!is2d)
      outfz = MRIconvolve1d(mbxby, outfz, &MRIFvox(mri_derfilter, 0, 0, 0) , klen, whd[3-1], 0, f);

    // compute outblur
    if (is2d)
      MRIcopyFrame(mbxby,outblur,0,f);
    else
      MRIconvolve1d(mbxby, outblur, &MRIFvox(mri_prefilter, 0, 0, 0) , klen, whd[3-1], 0, f);
    MRIfree(&mbxby);

  }

//   MRI* mdz;
//   MRI* mbz;
//   if (mri->depth > 1)
//   {
//     mdz =  MRIconvolve1d(mri, NULL, &MRIFvox(mri_derfilter, 0, 0, 0), klen, whd[3-1], 0, 0) ;
//     //MRI* mdz   = convolute(mri,mri_derfilter,3);
// 
//     mbz =  MRIconvolve1d(mri, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[3-1], 0, 0) ;
//     //MRI* mbz   = convolute(mri,mri_prefilter,3);
//   }
//   else
//   {
//     mdz = MRIcopy(mri,NULL);
//     mbz = MRIcopy(mri,NULL);
//   }
// 
//   MRI* mdzby = MRIconvolve1d(mdz, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[2-1], 0, 0) ;
//   //MRI* mdzby = convolute(mdz,mri_prefilter,2);
//   
//   MRIfree(&mdz);
//   outfz = MRIconvolve1d(mdzby, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[1-1], 0, 0) ;
//   //outfz = convolute(mdzby,mri_prefilter,1);
//   MRIfree(&mdzby);
//   
//   MRI* mbzby = MRIconvolve1d(mbz, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[2-1], 0, 0) ;
//   //MRI* mbzby = convolute(mbz,mri_prefilter,2);
//   MRI* mbzdy = MRIconvolve1d(mbz, NULL, &MRIFvox(mri_derfilter, 0, 0, 0), klen, whd[2-1], 0, 0) ;
//   //MRI* mbzdy = convolute(mbz,mri_derfilter,2);
//   MRIfree(&mbz);
//   outfy = MRIconvolve1d(mbzdy, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[1-1], 0, 0) ;
//   //outfy = convolute(mbzdy,mri_prefilter,1);
//   MRIfree(&mbzdy);
//   
//   outfx = MRIconvolve1d(mbzby, NULL, &MRIFvox(mri_derfilter, 0, 0, 0), klen, whd[1-1], 0, 0) ;
//   //outfx = convolute(mbzby,mri_derfilter,1);
//   outblur = MRIconvolve1d(mbzby, NULL, &MRIFvox(mri_prefilter, 0, 0, 0), klen, whd[1-1], 0, 0) ;
//   //outblur = convolute(mbzby,mri_prefilter,1);
//   MRIfree(&mbzby);

  //cout << " size fx: " << outfx->width << " " << outfx->height << " " << outfx->depth << endl;

  MRIfree(&mri_prefilter);
  MRIfree(&mri_derfilter);

  return true;
}


std::vector<int> MyMRI::findRightSize(MRI *mri, float conform_size,
    bool conform)
// determines width, height, depth to cover mri with conform_size isotropic voxels
//  bool conform makes cube image (w=h=d) and adjust to min of 256
{
  //string s = ""; if (conform) s = " , 256-conform" ;
  //cout << " MyMRI::findRightSize(mri, target_voxel_size=" << conform_size << s << " )" << endl;
  double xsize, ysize, zsize;
  double fwidth, fheight, fdepth, fmax;
  int conform_width;

  xsize = fabs(mri->xsize);
  ysize = fabs(mri->ysize);
  zsize = fabs(mri->zsize);

  // now decide the conformed_width
  // calculate the size in mm for all three directions
  fwidth  = xsize * mri->width;
  fheight = ysize * mri->height;
  fdepth  = zsize * mri->depth;
  
  //cout << "in: " <<mri->width << " " << mri->height << " " << mri->depth << endl;
  //cout << "out: " << fwidth << " " << fheight << " " << fdepth << endl;

  vector<int> ret(3);
  double eps = 0.0001; // to prevent ceil(2.0*64 / 2.0) = ceil(64.000000000001) = 65
  ret[0] = (int) ceil((fwidth / conform_size) - eps);
  ret[1] = (int) ceil((fheight / conform_size) - eps);
  ret[2] = (int) ceil((fdepth / conform_size) - eps);
  if (mri->depth == 1)
    ret[2] = 1; // 2d image

  //cout << " zsize: " << zsize << " depth: " << mri->depth << " fdepth: " << fdepth << " new depth: " << ret[2] << endl;

  if (!conform)
    return ret;

  // pick the largest image side length
  fmax = fwidth;
  if (fheight > fmax)
    fmax = fheight;
  if (mri->depth > 1 && fdepth > fmax)
    fmax = fdepth;

  // get the width with conform_size
  conform_width = (int) ceil(fmax / conform_size);

  // just to make sure that if smaller than 256, use 256 anyway
  if (conform_width < 256)
    conform_width = 256;
  // conform_width >= 256.   allow 10% leeway
  else if ((conform_width - 256.) / 256. < 0.1)
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
  if (mri->depth == 1)
    ret[2] = 1; // 2d image

  return ret;
}

bool MyMRI::isConform(MRI * mri)
{

  return (mri->xsize == 1 && mri->ysize == 1 && mri->zsize == 1
      && mri->width == 256 && mri->height == 256 && mri->depth == 256);

}

bool MyMRI::isIsotropic(MRI * mri)
{
  if (mri->depth == 1) return (mri->xsize == mri->ysize);
  
  return (mri->xsize == mri->ysize && mri->xsize == mri->zsize);

}

MATRIX* MyMRI::MRIgetZslice(MRI * mri, int slice, int frame)
// extract a z-slice from mri
// return as matrix
{
  assert(slice >=0 && slice < mri->depth);
  MATRIX* m = MatrixAlloc(mri->height, mri->width, MATRIX_REAL);
  int x, y;
  for (y = 0; y < mri->height; y++)
    for (x = 0; x < mri->width; x++)
    {
      if (mri->type == MRI_FLOAT)
        *MATRIX_RELT(m, y+1, x+1) = MRIFseq_vox(mri, x, y, slice, frame) ;
        else if (mri->type == MRI_UCHAR || mri->type == MRI_INT)
        *MATRIX_RELT(m, y+1, x+1) = (int)MRIseq_vox(mri,x,y,slice,frame);
        else (assert (1==2));
      }

  return m;
}

MRI * MyMRI::MRIlinearTransform(MRI* mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& m)
{
  MATRIX * mm = MyMatrix::convertVNL2MATRIX(m.as_matrix(), NULL);
  mriT = ::MRIlinearTransform(mriS, mriT, mm);
  MatrixFree(&mm);
  return mriT;
}

vnl_matrix_fixed<double, 4, 4> MyMRI::MRIvoxelXformToRasXform(MRI * mri_src,
    MRI * mri_dst, const vnl_matrix_fixed<double, 4, 4> &m_vox)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras;

  if (!mri_dst)
    mri_dst = mri_src; /* assume they will be in the same space */

  m_ras_to_voxel = MRIgetRasToVoxelXform(mri_src);
  m_voxel_to_ras = MRIgetVoxelToRasXform(mri_dst);

  //m_tmp = MatrixMultiply(m_voxel_xform, m_ras_to_voxel, NULL) ;
  //m_ras_xform = MatrixMultiply(m_voxel_to_ras, m_tmp, m_ras_xform) ;

  vnl_matrix_fixed<double, 4, 4> m_ras = MyMatrix::convertMATRIX2VNL(
      m_voxel_to_ras) * m_vox * MyMatrix::convertMATRIX2VNL(m_ras_to_voxel);

  MatrixFree(&m_voxel_to_ras);
  MatrixFree(&m_ras_to_voxel);

  return (m_ras);

}

MRI * MyMRI::gaussianCube(int size)
{
  double sigma = size / (4 * sqrt(2 * log(2)));
  cout << "MyMRI::gaussianCube( size: " << size << " )  sigma: " << sigma
      << endl;

  MRI * g = MRIalloc(size, size, size, MRI_FLOAT);
  int x, y, z;
  double xs, ys, zs;
  int hsize = (size - 1) / 2;
  sigma = 2.0 * sigma * sigma;
  double sum = 0.0;
  double dtmp;
  for (z = 0; z < size; z++)
  {
    zs = z - hsize;
    zs = zs * zs / sigma;
    for (y = 0; y < size; y++)
    {
      ys = y - hsize;
      ys = ys * ys / sigma;
      for (x = 0; x < size; x++)
      {
        xs = x - hsize;
        xs = xs * xs / sigma;
        dtmp = exp(-xs - ys - zs);
        MRIsetVoxVal(g, x, y, z, 0, dtmp);
        sum += dtmp;
      }
    }
  }
  // normalize
  for (z = 0; z < size; z++)
    for (y = 0; y < size; y++)
      for (x = 0; x < size; x++)
        MRIsetVoxVal(g, x, y, z, 0, MRIgetVoxVal(g, x, y, z, 0) / (float) sum);

  return g;
}

MRI * MyMRI::setTypeUCHAR(MRI * mri)
{
  int no_scale_flag = FALSE;
  printf("       -- changing data type from %d to %d (noscale = %d)...\n",
      mri->type, MRI_UCHAR , no_scale_flag);
    // not sure what happens if mri is already UCHAR?
  MRI * mri2 = MRISeqchangeType(mri, MRI_UCHAR, 0.0, 0.999, no_scale_flag);
  if (mri2 == NULL)
  {
    printf("ERROR: MRISeqchangeType in MyMRI::setTypeUCHAR\n");
    exit(1);
  }
  return mri2;
}

void MyMRI::setMaxOutsideVal(MRI * mri)
{
  float maxval = MRIgetVoxVal (mri,0,0,0,0);
  int w,h,d;
  for (d = 0; d<mri->depth; d++)
  {
    for (h = 0; h<mri->height; h++)
    {
      for (w = 0; w<mri->width; w++)
      {
         const float & val = MRIgetVoxVal (mri, w, h, d,0);
         if (val > maxval) maxval = val;
      }
    }
  }
  mri->outside_val = maxval;
}


float MyMRI::getBackground(MRI * mri)
// simple approach : count min and max values
// this will fail if background is gray
{
  float minval = MRIgetVoxVal (mri,0,0,0,0);
  float maxval = minval;
  int w,h,d;
  for (d = 0; d<mri->depth; d++)
  {
    for (h = 0; h<mri->height; h++)
    {
      for (w = 0; w<mri->width; w++)
      {
         const float & val = MRIgetVoxVal (mri, w, h, d,0);
         if (val > maxval) maxval = val;
         if (val < minval) minval = val;
      }
    }
  }

  int mins = 0;
  int maxs = 0;
  float eps = (maxval-minval)/255.0;
  for (d = 0; d<mri->depth; d++)
  {
    for (h = 0; h<mri->height; h++)
    {
      for (w = 0; w<mri->width; w++)
      {
         const float & val = MRIgetVoxVal (mri, w, h, d,0);
         if (maxval-val<eps) maxs++;
         if (val-minval<eps) mins++;
      }
    }
  }
  //std::cout << "min: " << minval << " (# " << mins<<" ),   max: " << maxval << " (# " <<maxs<<" ) " << std::endl; 
  if (mins > maxs) return minval;
  else return maxval;
  
}

double MyMRI::entropyPatch(MRI * mri, int x, int y, int z, int radius,
    int nbins, MRI* kernel, bool ball)
// expects uchar image between 0..(nbins-1) intensity values
{
  //int width  = mri->width;
  //int height = mri->height;
  //int depth  = mri->depth;

  // compute neighborhood bounds
  int xMin = x - radius;
  int xMax = x + radius;
  int yMin = y - radius;
  int yMax = y + radius;
  int zMin = z - radius;
  int zMax = z + radius;
  //int xr = xMin;
  //int yr = yMin;
  //int zr = zMin;
  //if (xMin < 0) xMin = 0;
  // if (xMax > width - 1) xMax = width - 1;
  // if (yMin < 0) yMin = 0;
  // if (yMax > height - 1) yMax = height - 1;
//  if (zMin < 0) zMin = 0;
//  if (zMax > depth - 1) zMax = depth - 1;

  // compute histogram of pixel values in neighborhood
  double histo[nbins];
  for (int o = 0; o < nbins; o++)
    histo[o] = 0.0f;

  int minI = MRIvox(mri, x, y, z) ;
  int maxI = minI;
  int index, dx, dy, dz, zp, yp, xp;
  double r2 = radius * radius;
  //cout << " minI: " << minI << " maxI: " << maxI << " float: " << MRIgetVoxVal(mri,x,y,z,0) << endl;
  for (zp = zMin; zp <= zMax; zp++)
  {
    for (yp = yMin; yp <= yMax; yp++)
    {
      for (xp = xMin; xp <= xMax; xp++)
      {
        //cout <<  xp << " " << yp << " " << zp << endl;
        dx = xp - xMin - radius;
        dy = yp - yMin - radius;
        dz = zp - zMin - radius;
        if (ball)
        {
          if (dx * dx + dy * dy + dz * dz > r2)
            continue;
        }
        index = MRIvox(mri, xp, yp, zp) ;
        if (index < minI)
          minI = index;
        if (index > maxI)
          maxI = index;
        //assert(index>=0);
        //assert(index<nbins);
        histo[index] += MRIFvox(kernel,dx,dy,dz) ;
      }
    }
  }

        // empty or all same value:
  if (minI == maxI)
  {
    //cout << " minI == maxI == " << minI << endl;
    // HISTOfree(&h) ;
    return 0.0;
  }

//  // smooth histo
//  HISTOGRAM *hs = HISTOalloc(nbins) ;
//  hs->bin_size = h->bin_size ;
//  hs->min = h->min ;
//  hs->max = h->max ;
//  float smooth [5] = { 0.05, 0.25, 0.4, 0.25, 0.05 };
//  float total;
//  float alltotal = 0.0;
//  int pos,b,b1,kx;
//  for (b = 0 ; b < nbins ; b++)
//  {
//    for (total = 0.0f, pos = 0 ; pos < 5 ; pos++)
//    {
//      kx = pos - 2 ;
//      b1 = b + kx ;
//      if (b1 >= nbins || b1 < 0)
//        continue ;
//
//      total += smooth[pos] * (float)h->counts[b1] ;
//    }
//    hs->counts[b] = total ;
//    hs->bins[b] = h->bins[b] ;
//    alltotal += total;
//  }

  // entropy (on normalized histo)
  double alltotal = 0.0f;
  for (int o = 0; o < nbins; o++)
    alltotal += histo[o];

  double entropy = 0.0, temp;
  for (int b = 0; b < nbins; b++)
  {
    //if (h->counts[b] > 0)
    if (histo[b] > 0)
    {
      //temp = (double)hs->counts[b]/alltotal ;
      temp = (double) histo[b] / alltotal;
      // shannon
      entropy -= temp * log(temp);
      // burg
      //entropy += log(temp);
      // Renyi, alpha = 2 : ent = - log( sum(histo(idx).^2) );
      //entropy += temp*temp;
    }
  }
  // Renyi
  //entropy = -log(entropy);

  return entropy;
}

// MRI * MyMRI::entropyImage(MRI* mri, int radius )
// {
// 
//   int nbins = 64;
// 
//   int width  = mri->width;
//   int height = mri->height;
//   int depth  = mri->depth;
// 
//   // convert to uchar (0..255)
//   MRI * mriIn;
//   if (mri->type != MRI_UCHAR)
//   {
//     int no_scale_flag = FALSE;
//     printf("changing data type from %d to %d (noscale = %d)...\n",
//            mri->type,MRI_UCHAR,no_scale_flag);
//     mriIn  = MRISeqchangeType(mri, MRI_UCHAR, 0.0, 0.999, no_scale_flag);
//     if (mriIn == NULL)
//     {
//       printf("ERROR: MRISeqchangeType\n");
//       exit(1);
//     }
//    // MRIwrite(mriIn,"mriIn.mgz");
//   }
//   else
//   {
//      mriIn = MRIcopy(mri,NULL);
//   }
// 
//   // scale factor for image to fit bin size (0...(nbins-1))
//   float factor = nbins/256.0;
//   //cout << "factor :  " << factor << endl;
//   int x,y,z;
//   for (z=0;z<depth;z++)
//     for (y=0;y<height;y++)
//       for (x=0;x<width;x++)
//           MRIsetVoxVal(mriIn,x,y,z,0,(int)(MRIgetVoxVal(mriIn,x,y,z,0)*factor));
//   
//   //MRIwrite(mriIn,"mriin.mgz");
//   
//   
//   // get gaussian kernel
//   MRI * kernel = gaussianCube(2*radius+1);
//   //MRIwrite(mriIn,"mrikernel.mgz");
// 
//   
//   // get Entropy image
//   MRI * entI =  MRIalloc(width,height,depth, MRI_FLOAT);
//   entI = MRIcopyHeader(mri,entI);
//   entI->type = MRI_FLOAT;
// 
//   float e;
//   float minEntropy = 1000000;
//   float maxEntropy = 0;
//   int zinterval = depth/10;
//   int zinterval2 = depth/50;
//   int xmax = width-radius-1;
//   int ymax = height-radius-1;
//   int zmax = depth-radius-1;
//   for (z=0;z<depth;z++)
//   {
//     if ((z+1)%zinterval == 0) cout << (z+1)/zinterval << flush;
//     else  if ((z+1)%zinterval2 ==0) cout << "." << flush;
//         
//     for (y=0;y<height;y++)
//     {
//       for (x=0;x<width;x++)
//       {
//           //zero padding
//           if (x < radius || y< radius || z<radius || x>xmax || y>ymax || z>zmax)
//           {
//             MRIsetVoxVal(entI,x,y,z,0,0.0);
//             continue;
//           }
//             
//           // inside
//           e = (float) entropyPatch(mriIn,x,y,z,radius,nbins,kernel);
//           if (e<0 || isnan(e) )
//           {
//              cerr<< " ERROR in MyMRI::entropyImage: entropy is negative or nan ? " << endl;
//           }
//           MRIsetVoxVal(entI,x,y,z,0,e);
//           if (e < minEntropy)
//             minEntropy = e;
//           if (e > maxEntropy)
//             maxEntropy = e;
//       }
//     }
//   }
// 
//   cout << endl << " Min Entropy: " << minEntropy << "  Max Entropy: " << maxEntropy << endl;
// 
//   // scale entropy image to 0..255 and make uchar
// //  int no_scale_flag = FALSE;
// //  MRI* mriEnt = MRISeqchangeType(entI, MRI_UCHAR, 0.0, 0.999, no_scale_flag);
// //  if (mriEnt == NULL)
// //  {
// //    printf("ERROR: MRISeqchangeType\n");
// //    exit(1);
// //  }
// //
// //  //MRIwrite(entI,"mriEntFloat.mgz");
// 
//   // cleanup
// //  MRIfree(&entI);
//   MRIfree(&mriIn);
//   MRIfree(&kernel);
//   
// //  return mriEnt;
//   return entI;
// }

double MyMRI::noiseVar(MRI * mri)
{
  cout << "MyMRI::noiseVar : " << flush;
  int width = mri->width - 1;
  int height = mri->height - 1;
  int depth = mri->depth - 1;
  int insize = (width - 1) * (height - 1) * (depth - 1);
  int x = 0, y = 0, z = 0;

  double epsi = 0.0;

//   MRI * mrifloat = NULL;
//   if (mri->type != MRI_FLOAT)
//   {
//     mrifloat = MRISeqchangeType(mri, MRI_FLOAT, 0.0, 1, TRUE);
//     if (mrifloat == NULL)
//     {
//       printf("ERROR: MRISeqchangeType\n");
//       exit(1);
//     }
//   }
//   else
//     mrifloat = mri;
// 
//   MRI *mri_6filter = MRIalloc(3,1,1, MRI_FLOAT);
//   MRIFvox(mri_6filter, 0, 0, 0) =  1.0/2.0 ;
//   MRIFvox(mri_6filter, 1, 0, 0) =  -1.0 ;
//   MRIFvox(mri_6filter, 2, 0, 0) =  1.0/2.0 ;
//   MRI * mri6 = MRIconvolveGaussian(mrifloat,NULL,mri_6filter);
//   MRIfree(&mri_6filter);
//   
//   MRIwrite(mri6,"mri6.mgz");

  double sum = 0.0;
  for (z = 1; z < depth; z++)
  {
    //if ((z+1)%zinterval == 0) cout << (z+1)/zinterval << flush;
    //else  if ((z+1)%zinterval2 ==0) cout << "." << flush;
    for (y = 1; y < height; y++)
    {
      for (x = 1; x < width; x++)
      {
        epsi = MRIgetVoxVal(mri, x - 1, y, z, 0)
            + MRIgetVoxVal(mri, x + 1, y, z, 0);
        epsi += MRIgetVoxVal(mri, x, y - 1, z, 0)
            + MRIgetVoxVal(mri, x, y + 1, z, 0);
        epsi += MRIgetVoxVal(mri, x, y, z - 1, 0)
            + MRIgetVoxVal(mri, x, y, z + 1, 0);
        epsi = MRIgetVoxVal(mri, x, y, z, 0) - (epsi / 6.0);
        sum += epsi * epsi;
      }
    }
  }

//  if (mri != mrifloat)
//    MRIfree(&mrifloat);

  sum = (sum / insize) * (6.0 / 7.0);
  cout.precision(12);
  cout << sum << endl;
  return sum;
}

MRI * MyMRI::nlsdImage(MRI * mri, int prad, int nrad)
{
  if (mri->nframes > 1)
  {
    cerr << "MyMRI::entropyImage multiple frames not supported!" << endl;
    exit(1);
  }

  int width = mri->width;
  int height = mri->height;
  int depth = mri->depth;
  int x = 0, y = 0, z = 0, xx = 0, yy = 0, zz = 0, xxx = 0, yyy = 0, zzz = 0;
  //int insize = width*height*depth;

  // allocate non-local shape descriptor image
  MRI * nlsdI = MRIalloc(width, height, depth, MRI_FLOAT);
  nlsdI = MRIcopyHeader(mri, nlsdI);
  nlsdI->type = MRI_FLOAT;
  MRIclear(nlsdI);

  int radius = prad + nrad;
  int xmax = width - radius - 1;
  int ymax = height - radius - 1;
  int zmax = depth - radius - 1;
  double wij = 0.0;
  double d = 0.0;
  double val = 0.0;
  int nsize1 = (nrad * 2) + 1;
  int nsize = nsize1 * nsize1 * nsize1;
  int count = 0;

  double nvar = noiseVar(mri);
  double sigma2 = sqrt(2) * nvar; // needs to be estimated from full image
  //double s2 = sqrt(2);

#ifdef HAVE_OPENMP
#pragma omp parallel for firstprivate(y,x,val,zz,yy,xx,zzz,yyy,xxx,d,wij,count) shared(depth,height,width,radius,nrad,prad,mri,nlsdI,nsize,sigma2) schedule(static,1)
#endif
  for (z = 0; z < depth; z++)
  {

    if (z % 10 == 0)
    {
      cout << "." << flush;
    }
    //if ((z+1)%zinterval == 0) cout << (z+1)/zinterval << flush;
    //else  if ((z+1)%zinterval2 ==0) cout << "." << flush;
    for (y = 0; y < height; y++)
    {
      for (x = 0; x < width; x++)
      {
        //zero padding at boarder:
        if (x < radius || y < radius || z < radius || x > xmax || y > ymax
            || z > zmax)
        {
          MRIsetVoxVal(nlsdI, x, y, z, 0, 0.0);
          continue;
        }

        val = 0.0;
        double * weights = new double[nsize];
        count = 0;
        // run over (non-local) neighborhood
        for (zz = -nrad; zz <= nrad; zz++)
          for (yy = -nrad; yy <= nrad; yy++)
            for (xx = -nrad; xx <= nrad; xx++)
            {

              // run over local neighborhood
              wij = 0.0;
              for (zzz = -prad; zzz <= prad; zzz++)
                for (yyy = -prad; yyy <= prad; yyy++)
                  for (xxx = -prad; xxx <= prad; xxx++)
                  {
                    d = MRIgetVoxVal(mri, x + xxx, y + yyy, z + zzz, 0)
                        - MRIgetVoxVal(mri, x + xx + xxx, y + yy + yyy,
                            z + zz + zzz, 0);
                    wij += d * d;
                  }
              wij = exp(-wij / sigma2);
              weights[count] = wij;
              count++;
              val += wij;

            }

        val /= nsize;
        d = 0.0;
        for (zz = 0; zz < nsize; zz++)
        {
          d += weights[zz] - val;
        }

        delete[] weights;

        MRIsetVoxVal(nlsdI, x, y, z, 0, d);

      }
    }
  }
  cout << endl;
  return nlsdI;
}

// MRI * MyMRI::nlsdImage(MRI * mri, int prad, int nrad)
// {
// 
//   int width  = mri->width;
//   int height = mri->height;
//   int depth  = mri->depth;
//   int x=0,y=0,z=0,xx=0,yy=0,zz=0;
//   double val = 0.0;
//   int nthreads=1, tid=0;
// #ifdef HAVE_OPENMP
// #pragma omp parallel 
//   {
//     nthreads = omp_get_num_threads();
//   }
// #endif
// 
//   // allocate non-local shape descriptor image
//   MRI * nlsdI =  MRIalloc(width,height,depth, MRI_FLOAT);
//   nlsdI = MRIcopyHeader(mri,nlsdI);
//   nlsdI->type = MRI_FLOAT;
//   MRIclear(nlsdI);
// 
//   // allocate tmp mri images (one for each thread)
//   MRI * mritmps[_MAX_FS_THREADS];
//   for (tid=0;tid<nthreads;tid++)
//   {
//     mritmps[tid] =  MRIalloc(width,height,depth, MRI_FLOAT);
//     if (mritmps[tid] == NULL)
//     {
//       printf("ERROR MyMRI::nlsdImage: mritmp allocation failed (memory overflow)!\n");
//       exit(1);
//     }
//   }
//   tid = 0;
//   
//   // lookup table for neighorhood offset coords:
//   int nsize = (nrad*2)+1;
//   nsize = nsize * nsize * nsize;
//   int nxoffset[nsize];
//   int nyoffset[nsize];
//   int nzoffset[nsize];
//   int count = 0;
//   for (zz = -nrad; zz <= nrad; zz++)
//   for (yy = -nrad; yy <= nrad; yy++)
//   for (xx = -nrad; xx <= nrad; xx++)
//   {
//     nxoffset[count] = xx;
//     nyoffset[count] = yy;
//     nzoffset[count] = zz;
//     count++;
//   }
//   
//   // pre compute sqrt(2) sigma^2
//   double nvar = noiseVar(mri);
//   double sigma2 = sqrt(2) * nvar; // needs to be estimated from full image
//   
//   // pre compute average filter for patch size:
//   int psize1 = (prad*2)+1;
//   MRI *avg_filter = MRIalloc(psize1,1,1, MRI_FLOAT);
//   for (int i = 0;i<psize1; i++) MRIFvox(avg_filter, i, 0, 0) = -1;
//   
//   
//   // run over neighborhood offsets (single loop for better paralellization)
//   tid = 0;
// #ifdef HAVE_OPENMP
// #pragma omp parallel for firstprivate(tid,z,y,x,zz,yy,xx,val) shared(depth,height,width,nxoffset,nyoffset,nzoffset,nlsdI) schedule(static,1)
// #endif
//   for (count = 0; count < nsize; count++)
//   {
// #ifdef HAVE_OPENMP
//     tid = omp_get_thread_num();
// #endif
// 
// #ifdef HAVE_OPENMP
//   #pragma omp critical
// #endif  
//     cout << " " << count << " ( " << tid << " ) " << endl;
//     // 1. subtract shifted image from original at current neighbor hood offset
//     for (z=0; z < depth;  z++)
//     {
//       zz = nzoffset[count] +z;
//       for (y=0; y < height; y++)
//       {
//         yy = nyoffset[count] +y;
//         for (x=0; x < width;  x++)
//         {
//           xx = nxoffset[count] + x;
//           val = MRIgetVoxVal(mri,x,y,z,0);
//           if ( xx >= 0 && yy >= 0 && zz >= 0 && xx < width && yy < height && zz < depth)
//             val -= MRIgetVoxVal(mri,xx,yy,zz,0);
//           // 2. compute pointwise square
//           MRIsetVoxVal(mritmps[tid],x,y,z,0,val*val);
//         }
//       }
//     }
//     
//     // 3. filter with average filter in patch box
//     MRI * mrit2 = MRIconvolveGaussian(mritmps[tid],NULL,avg_filter);
// 
//     // 4. compute (fixed j offset)  wij = exp ( - sum(squared differences) / (sqrt(2) sigma^2)
//     for (z=0; z < depth;  z++)
//     for (y=0; y < height; y++)
//     for (x=0; x < width;  x++)
//     {
//       val = exp( MRIgetVoxVal(mrit2,x,y,z,0) / sigma2);
//       //MRIsetVoxVal(mritemp[tid],x,y,z,0,val);
// #ifdef HAVE_OPENMP
//   #pragma omp critical
// #endif  
//       MRIFvox(nlsdI, x, y, z) += val;
//     }
//   
//     
//   
//     MRIfree(&mrit2);
//   }
//   
//   nlsdI = MRIscalarMul(nlsdI,nlsdI,1.0/nsize);
//   return nlsdI;
// }

MRI * MyMRI::entropyImage(MRI* mri, int radius, bool ball, bool correction, MRI * mask)
{

  if (mri->nframes > 1)
  {
    cerr << "MyMRI::entropyImage multiple frames not supported!" << endl;
    exit(1);
  }
  
  if (correction && radius > 1)
  {
    cout << "WARNING: using correction mode with radius > 1 may take a very long time !" << endl;
    cout << "         Correction usually allows small radis (e.g. 1 in 3D)." << endl;
  }

  //int nbins = 64;
  int nbins = 256;

  int width = mri->width;
  int height = mri->height;
  int depth = mri->depth;
  int x, y, z;
  bool is2d = (mri->depth == 1);

  int insize = width * height * depth;
  unsigned char * mriIn = new unsigned char[insize];
  if (!mriIn)
  {
    cout << "ERROR: Memory Overflow for mriIn in myMRI::entropyImage" << endl;
    exit(1);
  }
  unsigned char * mriMask = NULL;
  MRI* localMask = NULL;
  if (mask)
  {
    mriMask = new unsigned char[insize];
    if (!mriMask)
    {
      cout << "ERROR: Memory Overflow for mriMask in myMRI::entropyImage" << endl;
      exit(1);
    }
    MATRIX* m = MRIgetVoxelToVoxelXform(mask,mri);
    localMask = MRIcopy(mri,NULL);
    localMask = MRIlinearTransformInterp(mask,localMask,m,SAMPLE_NEAREST);
    MatrixFree(&m);
  }
  

  // -------------------------------------------------------------------------------------
  // scale image to fit bins (probably faster w/o using histogram)
//   double minIn = MRIgetVoxVal(mri,0,0,0,0);
//   double maxIn = minIn;
//   float tempIn;
//   for (z=0; z<depth; z++)
//   for (y=0; y<height;y++)
//   for (x=0; x<width; x++)
//   {
//     tempIn = MRIgetVoxVal(mri,x,y,z,0);
//     if (tempIn > maxIn) maxIn=tempIn;
//     if (tempIn < minIn) minIn=tempIn;
//   }
//   double factor = ((double)nbins) / (maxIn - minIn + 0.0001);
//   int count = 0;
//   double temp;
//   for (z=0; z<depth; z++)
//   for (y=0; y<height;y++)
//   for (x=0; x<width; x++)
//   {
//     temp = floor(factor * (MRIgetVoxVal(mri,x,y,z,0) - minIn));
//     assert(temp >= 0 && temp < nbins);
//     mriIn[count] = (unsigned char) temp;
//     count++;
//   }

  // -------------------------------------------------------------------------------------
  // convert to uchar (0..255)
  MRI * mriuchar = NULL;
  if (mri->type != MRI_UCHAR)
  {
    int no_scale_flag = FALSE;
    printf("  - changing data type from %d to %d (noscale = %d)...\n", mri->type,
        MRI_UCHAR, no_scale_flag);
    mriuchar = MRISeqchangeType(mri, MRI_UCHAR, 0.0, 0.999, no_scale_flag);
    if (mriuchar == NULL)
    {
      printf("ERROR: MRISeqchangeType\n");
      exit(1);
    }
    // MRIwrite(mriIn,"mriIn.mgz");
  }
  else
  {
    mriuchar = mri;
  }

  double factor = nbins / 256.0;
  double temp;
  int count = 0;
  //cout << "factor :  " << factor << endl;
  for (z = 0; z < depth; z++)
    for (y = 0; y < height; y++)
      for (x = 0; x < width; x++)
      {
        temp = (int) (MRIgetVoxVal(mriuchar, x, y, z, 0) * factor);
        assert(temp >= 0 && temp < nbins);
        mriIn[count] = (unsigned char) temp;
        
        // if mask:
        if (localMask)
        {
          temp = (int) MRIgetVoxVal(localMask, x, y, z, 0);
          if (temp > 0) temp = 1;
          if (temp < 0) temp = 0;
          mriMask[count] = (unsigned char) temp;
        }
        count++;
      }
  if (mriuchar != mri)
    MRIfree(&mriuchar);
  if (localMask)
    MRIfree(&localMask);

  // -------------------------------------------------------------------------------------
  // get gaussian kernel
  int ssize = 2 * radius + 1;
  double sigma = ssize / (4 * sqrt(2 * log(2)));
  cout << "  - compute Gaussian cube( size: " << ssize << " )  sigma: " << sigma
      << endl;
  int zssize = ssize;
  if (is2d) zssize = 1;
  int gsize = ssize * ssize * zssize;
  double* g = new double[gsize];
  if (!g)
  {
    cout << "ERROR: Memory Overflow for g in myMRI::entropyImage" << endl;
    exit(1);
  }
  double sum = 0.0;
  count = 0;
  double dtmp, zs, ys, xs;
  sigma = 2.0 * sigma * sigma;
  zs = 0;
  for (z = 0; z < zssize; z++)
  {
    if (!is2d)
    { 
      zs = z - radius;
      zs = zs * zs / sigma;
    }
    for (y = 0; y < ssize; y++)
    {
      ys = y - radius;
      ys = ys * ys / sigma;
      for (x = 0; x < ssize; x++)
      {
        xs = x - radius;
        xs = xs * xs / sigma;
        dtmp = exp(-xs - ys - zs);
        //assert(count < gsize);
        g[count] = dtmp;
        count++;
        sum += dtmp;
      }
    }
  }
  // normalize
  for (z = 0; z < gsize; z++)
    g[z] /= sum;

  // get Entropy image
  MRI * entI = MRIalloc(width, height, depth, MRI_FLOAT);
  entI = MRIcopyHeader(mri, entI);
  entI->type = MRI_FLOAT;

  //double minEntropy = 1000000;
  //double maxEntropy = 0;
  //int zinterval = depth/10;
  //int zinterval2 = depth/50;
  int xmax = width - radius - 1;
  int ymax = height - radius - 1;
  int zmax = depth - radius - 1;
  if (correction)
  {
    xmax--;
    ymax--;
    zmax--;
  }
  if (is2d) zmax = 0;
  double histo[nbins];
  double histosum = 0, entropy = 0, etmp = 0;
  int o = 0, zz = 0, yy = 0, xx = 0;
  count = 0;
  int pos = 0;
//  unsigned char index=0;
  x = 0;
  y = 0;
  int r2 = radius * radius;
  int wtimesh = width * height;
  //cout << depth << " " << height << " " << width << endl;
  int Nthreads = 1;
#ifdef HAVE_OPENMP
  Nthreads = omp_get_max_threads();
#endif  
  if (correction)
    cout << "  - compute entropy image now with correction (threads="<<Nthreads<<")... " << endl;
  else
    cout << "  - compute entropy image now (threads="<<Nthreads<<")... " << endl;


#ifdef HAVE_OPENMP
#pragma omp parallel for firstprivate(y,x,o,histo,count,zz,yy,xx,pos,histosum,entropy,etmp) shared(depth,height,width,radius,entI,nbins,ssize,xmax,ymax,zmax,mriIn,g) schedule(static,1)
#endif
  for (z = 0; z < depth; z++)
  {
    //if ((z+1)%zinterval == 0) cout << (z+1)/zinterval << flush;
    //else  if ((z+1)%zinterval2 ==0) cout << "." << flush;

    for (y = 0; y < height; y++)
    {
      for (x = 0; x < width; x++)
      {
        //zero padding at boarder:
        if (depth > 1)
        {
          if (x < radius || y < radius || z < radius || x > xmax || y > ymax
              || z > zmax)
          {
            MRIsetVoxVal(entI, x, y, z, 0, 0.0);
            continue;
          }
        }
        else // 2D
        {
          if (x < radius || y < radius || x > xmax || y > ymax)
          {
            MRIsetVoxVal(entI, x, y, z, 0, 0.0);
            continue;
          }
        }

        // outside mask, continue:
        if (mriMask)
        {
          int mpos = (z * height + y ) * width + x;
          if (mriMask[mpos] == 0 )
          {
            MRIsetVoxVal(entI, x, y, z, 0, 0.0);
            continue;
          }       
        }
        
        
        // inside compute entropy in box or ball:

        // init histo
        for (o = 0; o < nbins; o++)
          histo[o] = 0.0f;

        // compute histo
        count = 0;
        int zstart = -radius;
        int zend = radius + 1;
        int z2 = 0, y2 = 0, x2 = 0;
        int yp = 0, zp = 0;
        if (depth == 1) //is2d
        {
          zstart = 0;
          zend = 1;
        }
        for (zz = zstart; zz < zend; zz++)
        {
          if (ball)
            z2 = zz * zz - r2;
          zp = (z + zz) * height;
          for (yy = -radius; yy <= radius; yy++)
          {
            if (ball)
              y2 = z2 + yy * yy;
            yp = ((y + yy) + zp) * width;
            for (xx = -radius; xx <= radius; xx++)
            {
              if (ball)
              {
                x2 = y2 + xx * xx;
                if (x2 > 0)
                {
                  //cout << " continue : " << zz << " " << yy << " " << xx << " " << radius << " " << x2 << endl;
                  count++;
                  continue;
                }
              }
              // update the histogram
              pos = (x + xx) + yp;
              if (correction)
              {
                pos = (x + xx) + yp;                
                if (is2d)
                {
                  get3Dcorrection(histo, mriIn[pos], mriIn[pos + 1], mriIn[pos + width], mriIn[pos + width + 1], nbins);
                }
                else
                {
                  unsigned int voxel1, voxel2, voxel3, voxel4, voxel5, voxel6, voxel7, voxel8;
                  voxel1 = mriIn[pos];
                  voxel2 = mriIn[pos + 1];
                  voxel3 = mriIn[pos + width];
                  voxel4 = mriIn[pos + wtimesh];
                  voxel5 = mriIn[pos + width + 1];
                  voxel6 = mriIn[pos + wtimesh + 1];
                  voxel7 = mriIn[pos + wtimesh + width];
                  voxel8 = mriIn[pos + wtimesh + width + 1];
                  // histo update happens in get3Dcorrection
                  // TET1        
                  get3Dcorrection(histo, voxel1, voxel2, voxel3, voxel4, nbins);
                  // TET2
                  get3Dcorrection(histo, voxel5, voxel2, voxel3, voxel8, nbins);
                  // TET3
                  get3Dcorrection(histo, voxel7, voxel8, voxel3, voxel4, nbins);
                  // TET4
                  get3Dcorrection(histo, voxel6, voxel2, voxel8, voxel4, nbins);
                  // TET5
                  get3Dcorrection(histo, voxel8, voxel2, voxel3, voxel4, nbins);
                }
              }
              else
              {
                pos = (x + xx) + yp;
                histo[mriIn[pos]] += g[count];
                count++;
              }
            }
          }
        }

        histosum = 0.0;
        for (o = 0; o < nbins; o++)
          histosum += histo[o];

        entropy = 0.0f;
        for (o = 0; o < nbins; o++)
          if (histo[o] > 0)
          {
            etmp = histo[o] / histosum;
            entropy -= etmp * log(etmp);
          }

        if (entropy < 0 || std::isnan(entropy))
        {
          cerr << " ERROR in MyMRI::entropyImage: entropy is negative or nan ? "
              << endl;
        }
        MRIsetVoxVal(entI, x, y, z, 0, (float) entropy);
        //if (entropy < minEntropy)
        //  minEntropy = entropy;
        //if (entropy > maxEntropy)
        //  maxEntropy = entropy;

      }
    }
  }

  //cout << endl << " Min Entropy: " << minEntropy << "  Max Entropy: " << maxEntropy << endl;

  // scale entropy image to 0..255 and make uchar
//  int no_scale_flag = FALSE;
//  MRI* mriEnt = MRISeqchangeType(entI, MRI_UCHAR, 0.0, 0.999, no_scale_flag);
//  if (mriEnt == NULL)
//  {
//    printf("ERROR: MRISeqchangeType\n");
//    exit(1);
//  }
//
//  //MRIwrite(entI,"mriEntFloat.mgz");

  // cleanup
  delete[] g;
  delete[] mriIn;
  if (mriMask) delete[] mriMask;

//  MRIwrite(entI,"enttest.mgz");
//  return mriEnt;

  if (correction)
  {
    // correct for half voxel shift (half voxel or half mm?)
    MATRIX * T1 = MatrixIdentity(4, NULL);
    *MATRIX_RELT(T1, 1, 4) = 0.5;
    *MATRIX_RELT(T1, 2, 4) = 0.5;
    *MATRIX_RELT(T1, 3, 4) = 0.5;
    entI = MRIlinearTransformInterp(entI, NULL, T1, SAMPLE_CUBIC_BSPLINE);
    MatrixFree(&T1);
  }

  return entI;
}

void MyMRI::get3Dcorrection(double* histo, unsigned int v1, unsigned int v2,
    unsigned int v3, unsigned int v4, unsigned int intRange)
{
  unsigned int vals[] =
  { v1, v2, v3, v4 };
  vector<unsigned int> vect(vals, vals + 4);
  sort(vect.begin(), vect.end());
  unsigned int I3 = vect[0];
  unsigned int I2 = vect[1];
  unsigned int I1 = vect[2];
  unsigned int I0 = vect[3];
  unsigned int d = I3;
  unsigned int c = I2 - d;
  unsigned int b = I1 - d;
  unsigned int a = I0 - d;
  vector<double> tempPDF(intRange, 0.0);
  unsigned int i;
  unsigned int vecLen;

  if (a == 0 && b == 0 && c == 0)
  {
    tempPDF[I3] = 1;
  }
  else
  {

    if (I3 != I2)
    {
      vecLen = I2 - I3 + 1;
      vector<double> tempVec(vecLen, 0.0);
      for (i = 0; i < vecLen; i++)
        tempVec[i] = I3 + i;
      tempVec[0] = tempVec[0] + 0.25;
      tempVec[vecLen - 1] = tempVec[vecLen - 1] - 0.25;

      double PDFvec[vecLen];
      for (i = 0; i < vecLen; i++)
        PDFvec[i] = 3.0f / a * (tempVec[i] - I3) / b * (tempVec[i] - I3) / c;
      PDFvec[0] = 0.5 * PDFvec[0];
      PDFvec[vecLen - 1] = 0.5 * PDFvec[vecLen - 1];

      for (i = I3; i <= I2; i++)
        tempPDF[i] = PDFvec[i - I3];
    }

    if (I2 != I1)
    {
      vecLen = I1 - I2 + 1;
      vector<double> tempVec(vecLen, 0.0);
      for (i = 0; i < vecLen; i++)
        tempVec[i] = I2 + i;
      tempVec[0] = tempVec[0] + 0.25;
      tempVec[vecLen - 1] = tempVec[vecLen - 1] - 0.25;

      double PDFvec[vecLen];
      for (i = 0; i < vecLen; i++)
        PDFvec[i] =
            (3.0f / a)
                * ((((I1 - tempVec[i]) / (b - c)) * ((tempVec[i] - I3) / b))
                    + (((I0 - tempVec[i]) / (a - c))
                        * ((tempVec[i] - I2) / (b - c))));
      PDFvec[0] = 0.5 * PDFvec[0];
      PDFvec[vecLen - 1] = 0.5 * PDFvec[vecLen - 1];

      for (i = I2; i <= I1; i++)
        tempPDF[i] = PDFvec[i - I2];
    }

    if (I1 != I0)
    {
      vecLen = I0 - I1 + 1;
      vector<double> tempVec(vecLen, 0.0);
      for (i = 0; i < vecLen; i++)
        tempVec[i] = I1 + i;
      tempVec[0] = tempVec[0] + 0.25;
      tempVec[vecLen - 1] = tempVec[vecLen - 1] - 0.25;

      double PDFvec[vecLen];
      for (i = 0; i < vecLen; i++)
        PDFvec[i] = (3.0f / a)
            * (((I0 - tempVec[i]) / (a - c)) * ((I0 - tempVec[i]) / (a - b)));
      PDFvec[0] = 0.5 * PDFvec[0];
      PDFvec[vecLen - 1] = 0.5 * PDFvec[vecLen - 1];

      for (i = I1; i <= I0; i++)
        tempPDF[i] = PDFvec[i - I1];
    }

    double sum = 0;
    for (i = 0; i < intRange; i++)
      sum += tempPDF[i];
    for (i = 0; i < intRange; i++)
      tempPDF[i] = tempPDF[i] / sum;
  }

  for (i = 0; i < intRange; i++)
    histo[i] = histo[i] + tempPDF[i];

}


MRI * MyMRI::getNormalizedImage(MRI *mri)
{
  MRI * mriN = MRIallocSequence(mri->width,mri->height, mri->depth, MRI_FLOAT, mri->nframes);
  MRIcopyHeader(mri,mriN);
  mriN->type = MRI_FLOAT;
  unsigned int N = mri->width * mri->height * mri->depth;
  
  // do each frame independently
  for (int f = 0; f< mri->nframes; f++)
  {
  
    // compute mean:
    double mean = 0.0;
    int y;
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:mean)  
#endif
    for (y = 0; y < mri->height; y++)
    {
      int z,x;
      double v;
      for (z = 0; z < mri->depth; z++)
      {
        for (x = 0; x < mri->width; x++)
        {
          v = MRIgetVoxVal(mri,x,y,z,f);
          mean += v;  
        }
      }
    }
    mean = mean / N;
    
    // compute std:
    double std = 0.0;
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+:std)  
#endif
    for (y = 0; y < mri->height; y++)
    {
      int z,x;
      double v;
      for (z = 0; z < mri->depth; z++)
      {
        for (x = 0; x < mri->width; x++)
        {
          v = MRIgetVoxVal(mri,x,y,z,f) - mean;
          std += v*v;  
        }
      }
    }
    std = sqrt(std / (N-1.5));
    
    // fill output image
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) 
#endif
    for (y = 0; y < mri->height; y++)
    {
      int z,x;
      double v;
      for (z = 0; z < mri->depth; z++)
       {
        for (x = 0; x < mri->width; x++)
        {
          v = (MRIgetVoxVal(mri,x,y,z,f) - mean) / std;
          MRIsetVoxVal(mriN,x,y,z,f,v);
        }
      }
    }
    

  } // for each frame
  
  return mriN;
}
  
MRI * MyMRI::getNormalizedImage(MRI *mri, int blockradius)
{
  double stdeps = 0.0001;
  
  MRI * mriN = MRIallocSequence(mri->width,mri->height, mri->depth, MRI_FLOAT, mri->nframes);
  MRIcopyHeader(mri,mriN);
  mriN->type = MRI_FLOAT;
  
  int bw = 2*blockradius+1;
  int bw3 = bw * bw * bw;
  int dt[4] = { mri->width, mri->height, mri->depth , mri->nframes };
  dt[0] = dt[0] - 1 - blockradius;
  dt[1] = dt[1] - 1 - blockradius;
  dt[2] = dt[2] - 1 - blockradius;
  int blockradiusz = blockradius;
  int bwz = bw;
  if (mri->depth == 1) // 2D case
  {
    blockradiusz = 0;
    bwz = 1;
    bw3 = bw * bw;
    dt[2] = 1;
  }

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) 
#endif
  for (int y = blockradius; y < dt[1]; y++)
  {
    int x,z,f,i,j,k;
    
    double v, mean, std;
    
    // precompute xyz coords for all z in block (include translation)
    int ny, nx, nz, nym, nxm, nzm;
    ny = y - blockradius;
    nym = ny + bw;
    
    for (x = blockradius; x < dt[0]; x++)
    {
      nx = x - blockradius;
      nxm = nx+bw;
      
      for (z = blockradiusz; z < dt[2]; z++)
      {
        nz = z - blockradiusz;
        nzm = nz + bwz;
        
        for (f = 0; f < dt[3]; f++)
        {
          // loop over block to compute mean
          mean = 0.0;
          for (i=nz;i<nzm;i++)
          for (j=ny;j<nym;j++)
          for (k=nx;k<nxm;k++)
          {
            v = MRIgetVoxVal(mri,k,j,i,f);
            mean += v;
          }
          mean /= bw3;
          
          // compute std for block
          std = 0.0;
          for (i=nz;i<nzm;i++)
          for (j=ny;j<nym;j++)
          for (k=nx;k<nxm;k++)
          {
            v = MRIgetVoxVal(mri,k,j,i,f) - mean;
            std += v*v;
          }
          std = sqrt(std/(bw3 - 1.5));  // is this also true in 2D?
          
          if (std < stdeps) std = 1.0;
          
          v = (MRIgetVoxVal(mri,x,y,z,f) - mean) / std;
          MRIsetVoxVal(mriN,x,y,z,f,v);
        }
      }
    }
  }
  return mriN;
}

MRI * MyMRI::meanFilter(MRI *mri, int blockradius)
{
  MRI * mriN = MRIallocSequence(mri->width,mri->height, mri->depth, MRI_FLOAT, mri->nframes);
  MRIcopyHeader(mri,mriN);
  mriN->type = MRI_FLOAT;
  
  int bw = 2*blockradius+1;
  int bw3 = bw * bw * bw;
  int dt[4] = { mri->width, mri->height, mri->depth , mri->nframes };
  dt[0] = dt[0] - 1 - blockradius;
  dt[1] = dt[1] - 1 - blockradius;
  dt[2] = dt[2] - 1 - blockradius;
  int blockradiusz = blockradius;
  int bwz = bw;
  if (mri->depth == 1) // 2D case
  {
    blockradiusz = 0;
    bwz = 1;
    bw3 = bw * bw;
    dt[2] = 1;
  }

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) 
#endif
  for (int y = blockradius; y < dt[1]; y++)
  {
    int x,z,f,i,j,k;
    
    double mean;
    
    // precompute xyz coords for all z in block (include translation)
    int ny, nx, nz, nym, nxm, nzm;
    
    ny = y - blockradius;
    nym = ny + bw;
    for (x = blockradius; x < dt[0]; x++)
    {
      nx = x - blockradius;
      nxm = nx+bw;
      
      for (z = blockradiusz; z < dt[2]; z++)
      {
        nz = z - blockradiusz;
        nzm = nz + bwz;
        
        for (f = 0; f < dt[3]; f++)
        {
          // loop over block to compute mean
          mean = 0.0;
          for (i=nz;i<nzm;i++)
          for (j=ny;j<nym;j++)
          for (k=nx;k<nxm;k++)
          {
            mean += MRIgetVoxVal(mri,k,j,i,f);
          }
          mean /= bw3;
                    
          MRIsetVoxVal(mriN,x,y,z,f,mean);
        }
      }
    }
  }
  return mriN;
}

MRI * MyMRI::medianFilter(MRI *mri, int blockradius)
{
  MRI * mriN = MRIallocSequence(mri->width,mri->height, mri->depth, mri->type, mri->nframes);
  MRIcopyHeader(mri,mriN);
  //mriN->type = MRI_FLOAT;
  
  int bw = 2*blockradius+1;
  int bw3 = bw * bw * bw;
  int dt[4] = { mri->width, mri->height, mri->depth , mri->nframes };
  dt[0] = dt[0] - 1 - blockradius;
  dt[1] = dt[1] - 1 - blockradius;
  dt[2] = dt[2] - 1 - blockradius;
  int blockradiusz = blockradius;
  int bwz = bw;
  if (mri->depth == 1) // 2D case
  {
    blockradiusz = 0;
    bwz = 1;
    bw3 = bw * bw;
    dt[2] = 1;
  }

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) 
#endif
  for (int y = blockradius; y < dt[1]; y++)
  {
    int x,z,f,i,j,k,count;
    
    double median;
    float* t = (float *) calloc(bw3, sizeof(float));
    
    // precompute xyz coords for all z in block (include translation)
    int ny, nx,nz, nym, nxm, nzm;
    ny = y - blockradius;
    nym = ny + bw;
    
    for (x = blockradius; x < dt[0]; x++)
    {
      nx = x - blockradius;
      nxm = nx+bw;
      
      for (z = blockradiusz; z < dt[2]; z++)
      {
        nz = z - blockradiusz;
        nzm = nz + bwz;
        
        for (f = 0; f < dt[3]; f++)
        {
          // loop over block to compute mean
          count = 0;
          for (i=nz;i<nzm;i++)
          for (j=ny;j<nym;j++)
          for (k=nx;k<nxm;k++)
          {
            
            t[count] = MRIgetVoxVal(mri,k,j,i,f);
            count++;
          }
          median = RobustGaussian<float>::median(t, bw3);
                    
          MRIsetVoxVal(mriN,x,y,z,f,median);
        }
      }
    }
    free(t);
  }
  return mriN;
}
