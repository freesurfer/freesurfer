#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include "mri.h"
#include "error.h"
#include "histo.h"
#include "mri_conform.h"

extern int errno;

MRI *conform_type(MRI *mri);
MRI *conform_voxels(MRI *mri);
MRI *conform_size(MRI *mri);
MRI *conform_direction(MRI *mri);

MATRIX *MRIgetConformMatrix(MRI *mri)
{

  MRI *templ;
  MATRIX *m_resample ;

  if(mri->ras_good_flag == 0)
  {
    mri->x_r = -1.0;  mri->x_a =  0.0;  mri->x_s =  0.0;
    mri->y_r =  0.0;  mri->y_a =  0.0;  mri->y_s = -1.0;
    mri->z_r =  0.0;  mri->z_a =  1.0;  mri->z_s =  0.0;
 } 

  templ = MRIallocHeader(256, 256, 256, MRI_UCHAR);

  templ->imnr0 = 1;
  templ->imnr1 = 256;
  templ->thick = 1.0;
  templ->ps = 1.0;
  templ->xsize = templ->ysize = templ->zsize = 1.0;
  templ->xstart = templ->ystart = templ->zstart = -128.0;
  templ->xend = templ->yend = templ->zend = 128.0;
  templ->x_r = -1.0;      templ->x_a =  0.0;      templ->x_s =  0.0;
  templ->y_r =  0.0;      templ->y_a =  0.0;      templ->y_s = -1.0;
  templ->z_r =  0.0;      templ->z_a =  1.0;      templ->z_s =  0.0;
  templ->c_r =  0.0;      templ->c_a =  0.0;      templ->c_s =  0.0;
  // definition of the conform one is that c_(r,a,s) = 0 at (128,128,128).
  templ->slice_direction = MRI_CORONAL;
  templ->tr = mri->tr ; templ->te = mri->te ; 
  templ->flip_angle = mri->flip_angle ; templ->ti = mri->ti ; 


  m_resample = MRIgetResampleMatrix(mri, templ);

  MRIfree(&templ) ;

  return(m_resample);
}

MRI *MRIconform(MRI *mri)
{

  MRI *templ, *mri2, *res;

  res = MRIcopy(mri, NULL); /* don't mess with the input */

  if(res->ras_good_flag == 0)
  {
    res->x_r = -1.0;  res->x_a =  0.0;  res->x_s =  0.0;
    res->y_r =  0.0;  res->y_a =  0.0;  res->y_s = -1.0;
    res->z_r =  0.0;  res->z_a =  1.0;  res->z_s =  0.0;
 } 

  templ = MRIallocHeader(256, 256, 256, MRI_UCHAR);

  templ->imnr0 = 1;
  templ->imnr1 = 256;
  templ->thick = 1.0;
  templ->ps = 1.0;
  templ->xsize = templ->ysize = templ->zsize = 1.0;
  templ->xstart = templ->ystart = templ->zstart = -128.0;
  templ->xend = templ->yend = templ->zend = 128.0;
  templ->x_r = -1.0;      templ->x_a =  0.0;      templ->x_s =  0.0;
  templ->y_r =  0.0;      templ->y_a =  0.0;      templ->y_s = -1.0;
  templ->z_r =  0.0;      templ->z_a =  1.0;      templ->z_s =  0.0;
  templ->c_r = res->c_r;  templ->c_a = res->c_a;  templ->c_s = res->c_s;
  templ->slice_direction = MRI_CORONAL;
  templ->tr = mri->tr ; templ->te = mri->te ; 
  templ->flip_angle = mri->flip_angle ; templ->ti = mri->ti ; 

  /* ----- change type if necessary ----- */
  if(res->type != templ->type)
  {
    mri2 = MRIchangeType(res, templ->type, 0.0, 0.999, FALSE);
    MRIfree(&res);
    if(mri2 == NULL)
      return(NULL);
    res = mri2;
  }


  /* ----- reslice if necessary ----- */
  if(res->xsize != templ->xsize || res->ysize != templ->ysize   || res->zsize != templ->zsize ||
     res->width != templ->width || res->height != templ->height || res->depth != templ->depth ||
     res->x_r != templ->x_r || res->x_a != templ->x_a || res->x_s != templ->x_s ||
     res->y_r != templ->y_r || res->y_a != templ->y_a || res->y_s != templ->y_s ||
     res->z_r != templ->z_r || res->z_a != templ->z_a || res->z_s != templ->z_s)
  {
    mri2 = MRIresample(res, templ, RESAMPLE_INTERPOLATE);
    MRIfree(&res);
    if(mri2 == NULL)
      return(NULL);
    res = mri2;
  }

  return(res);

#if 0

  MRI *temp;
  int copied_flag = 0;

  /* conform data type if needed */
  if(mri->type != MRI_UCHAR)
  {
    temp = conform_type(mri);
    MRIfree(&mri);
    mri = temp;
    copied_flag = 1;
  }


  /* conform voxel sizes if needed */
  if(!(mri->xsize == 1 && mri->ysize == 1 && mri->zsize == 1))
  {
    temp = conform_voxels(mri);
    MRIfree(&mri);
    mri = temp;
    copied_flag = 1;
  }

  /* conform volume size if needed (256x256x256) */
  if(mri->width != 256 || mri->height != 256 || mri->depth != 256)
  {
    temp = conform_size(mri);
    MRIfree(&mri);
    mri = temp;
    copied_flag = 1;
  }

  /* conform slice direction if needed */
  if(mri->slice_direction != MRI_CORONAL)
  {
    temp = conform_direction(mri);
    MRIfree(&mri);
    mri = temp;
    copied_flag = 1;
  }

  /* don't return the same structure that was passed */
  if(!copied_flag)
  {
    temp = mri;
    mri = MRIcopy(temp, NULL);
  }

  return(mri);

#endif

}  /*  end MRIconform()  */

#define N_BINS 256
#define FRACTION_SCALE 0.0001

MRI *conform_type(MRI *mri)
{

  MRI *mri2;
  int x, y, z;
  float min = 0.0, max = 0.0;
  float scale;
  float this = 0.0;
  int counts[N_BINS];
  float bin_size;
  int i;
  int bin;
  int nv, nv_needed;
  float e1, e2;
  float e1d, e2d;
  float that;

  for(i = 0;i < N_BINS;i++)
    counts[i] = 0;

  if(mri->slices == NULL)
  {
    mri2 = MRIcopy(mri, NULL);
    mri2->type = MRI_UCHAR;
    return(mri2);
  }

  /* pixel value scaling goes here */

  if(mri->type == MRI_UCHAR)
    min = max = (float)MRIvox(mri, 0, 0, 0);
  if(mri->type == MRI_INT)
    min = max = (float)MRIIvox(mri, 0, 0, 0);
  if(mri->type == MRI_LONG)
    min = max = (float)MRILvox(mri, 0, 0, 0);
  if(mri->type == MRI_FLOAT)
    min = max = (float)MRIFvox(mri, 0, 0, 0);
  if(mri->type == MRI_SHORT)
    min = max = (float)MRISvox(mri, 0, 0, 0);

  for(x = 0;x < mri->width;x++)
    for(y = 0;y < mri->height;y++)
      for(z = 0;z < mri->depth;z++)
      {
      if(mri->type == MRI_UCHAR)
        this = (float)MRIvox(mri, x, y, z);
      if(mri->type == MRI_INT)
        this = (float)MRIIvox(mri, x, y, z);
      if(mri->type == MRI_LONG)
        this = (float)MRILvox(mri, x, y, z);
      if(mri->type == MRI_FLOAT)
        this = (float)MRIFvox(mri, x, y, z);
      if(mri->type == MRI_SHORT)
        this = (float)MRISvox(mri, x, y, z);

      if(this > max)
        max = this;
      if(this < min)
        min = this;
        
      }

  bin_size = (float)(max - min) / (float)N_BINS;

  for(x = 0;x < mri->width;x++)
    for(y = 0;y < mri->height;y++)
      for(z = 0;z < mri->depth;z++)
      {
      if(mri->type == MRI_UCHAR)
        this = (float)MRIvox(mri, x, y, z);
      if(mri->type == MRI_INT)
        this = (float)MRIIvox(mri, x, y, z);
      if(mri->type == MRI_LONG)
        this = (float)MRILvox(mri, x, y, z);
      if(mri->type == MRI_FLOAT)
        this = (float)MRIFvox(mri, x, y, z);
      if(mri->type == MRI_SHORT)
        this = (float)MRISvox(mri, x, y, z);

      if(this == max)
        bin = N_BINS - 1;
      else
        bin = (int)((this - (float)min) / bin_size);
      if(bin < 0)
      {
        printf("bin < 0\n");
        bin = 0;
      }
      if(bin > N_BINS - 1)
      {
        printf("bin > N_BINS - 1\n");
        bin = N_BINS - 1;
      }

      counts[bin]++;

      }

  nv_needed = (int)(mri->height * mri->width * mri->depth * FRACTION_SCALE);

  nv = 0;
  for(i = 0;i < N_BINS && nv < nv_needed;i++)
    nv += counts[i];

  if(i == -1)
  {
    errno = 0;
    ErrorExit(ERROR_BADPARM, "MRIconform (value scaling): histogram is too thin for\na clipping fraction of %g",
              FRACTION_SCALE);
  }

  e1 = (i-1) * bin_size + min;

  nv = 0;
  for(i = N_BINS-1;i >= 0 && nv < nv_needed;i--)
    nv += counts[i];

  if(i == N_BINS)
  {
    errno = 0;
    ErrorExit(ERROR_BADPARM, "MRIconform (value scaling): histogram is too thin for\na clipping fraction of %g",
              FRACTION_SCALE);
  }

  e2 = (i+1) * bin_size + min;

  e1d = FRACTION_SCALE * 255;
  e2d = 255 - e1d;

  scale = (e2d - e1d) / (e2 - e1);

  mri2 = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_UCHAR, 1);
  MRIcopyHeader(mri, mri2);

  for(x = 0;x < mri->width;x++)
    for(y = 0;y < mri->height;y++)
      for(z = 0;z < mri->depth;z++)
      {
      if(mri->type == MRI_UCHAR)
        this = (float)MRIvox(mri, x, y, z);
      if(mri->type == MRI_INT)
        this = (float)MRIIvox(mri, x, y, z);
      if(mri->type == MRI_LONG)
        this = (float)MRILvox(mri, x, y, z);
      if(mri->type == MRI_FLOAT)
        this = (float)MRIFvox(mri, x, y, z);
      if(mri->type == MRI_SHORT)
        this = (float)MRISvox(mri, x, y, z);

      that = scale * (this - e1) + e1d;

      if(that < 0)
        MRIvox(mri2, x, y, z) = 0;
      else if (that > 255)
        MRIvox(mri2, x, y, z) = 255;
      else
        MRIvox(mri2, x, y, z) = (unsigned char)that;

      }

  return(mri2);

} /* end conform_type() */

MRI *conform_voxels(MRI *mri)
{

  MRI *mri2;
  float fovx, fovy, fovz;
  float max_fov;
  float min_size;
  int x, y, z;
  int x0, y0, z0;
  float xs, ys, zs;
  int xi, yi, zi;
  float xf, yf, zf, xf2, yf2, zf2;
  float ras_scale;

  fovx = mri->width * mri->xsize;
  fovy = mri->height * mri->ysize;
  fovz = mri->depth * mri->zsize;

  min_size = 1;
  max_fov = 256.;

  if(mri->slices == NULL)
  {
    mri2 = MRIcopy(mri, NULL);
    mri2->fov = max_fov;
    mri2->width = mri2->height = mri2->depth = 256;
    mri2->xsize = mri2->ysize = mri2->zsize = min_size;
    return(mri2);
  }

  mri2 = MRIallocSequence(256, 256, 256, MRI_UCHAR, 1);

  mri2->xdir = mri->xdir;
  mri2->ydir = mri->ydir;
  mri2->zdir = mri->zdir;
  mri2->slice_direction = mri->slice_direction;

  mri2->fov = max_fov;
  mri2->xstart = mri2->ystart = mri2->zstart = -128 * min_size;
  mri2->xend = mri2->yend = mri2->zend = 128 * min_size;

  x0 = (int)(mri->width * (-mri->xstart / fovx));
  y0 = (int)(mri->height * (-mri->ystart / fovy));
  z0 = (int)(mri->depth * (-mri->zstart / fovz));

  for(x = 0;x < 256;x++)
    for(y = 0;y < 256;y++)
      for(z = 0;z < 256;z++)
      {

        xs = ((float)(x - 128) * mri2->xsize / mri->xsize + (float)x0);
        ys = ((float)(y - 128) * mri2->ysize / mri->ysize + (float)y0);
        zs = ((float)(z - 128) * mri2->zsize / mri->zsize + (float)z0);

        xi = (int)xs;
        yi = (int)ys;
        zi = (int)zs;

        if(xi < 1 || yi < 1 || zi < 1 || xi > mri->width - 2 || yi > mri->height - 2 || zi > mri->depth - 2)
          MRIvox(mri2, x, y, z) = 0;
        else
        {
          xf = xs - xi;
          yf = ys - yi;
          zf = zs - zi;
          xf2 = 1.-xf;
          yf2 = 1.-yf;
          zf2 = 1.-zf;

          MRIvox(mri2, x, y, z) = (xf)  * (yf)  * (zf)  *  (float)MRIvox(mri, xi + 1, yi + 1, zi + 1) + 
                                  (xf)  * (yf)  * (zf2) *  (float)MRIvox(mri, xi + 1, yi + 1, zi    ) + 
                                  (xf)  * (yf2) * (zf)  *  (float)MRIvox(mri, xi + 1, yi,     zi + 1) + 
                                  (xf)  * (yf2) * (zf2) *  (float)MRIvox(mri, xi + 1, yi,     zi    ) + 
                                  (xf2) * (yf)  * (zf)  *  (float)MRIvox(mri, xi,     yi + 1, zi + 1) + 
                                  (xf2) * (yf)  * (zf2) *  (float)MRIvox(mri, xi,     yi + 1, zi    ) + 
                                  (xf2) * (yf2) * (zf)  *  (float)MRIvox(mri, xi,     yi,     zi + 1) + 
                                  (xf2) * (yf2) * (zf2) *  (float)MRIvox(mri, xi,     yi,     zi    );
        }

      }

  if(mri->ras_good_flag)
  {
    ras_scale = sqrt(mri->x_r * mri->x_r + mri->x_a * mri->x_a + mri->x_s * mri->x_s);
    mri2->x_r = mri->x_r / ras_scale;
    mri2->x_a = mri->x_a / ras_scale;
    mri2->x_s = mri->x_s / ras_scale;
    ras_scale = sqrt(mri->y_r * mri->y_r + mri->y_a * mri->y_a + mri->y_s * mri->y_s);
    mri2->y_r = mri->y_r / ras_scale;
    mri2->y_a = mri->y_a / ras_scale;
    mri2->y_s = mri->y_s / ras_scale;
    ras_scale = sqrt(mri->z_r * mri->z_r + mri->z_a * mri->z_a + mri->z_s * mri->z_s);
    mri2->z_r = mri->z_r / ras_scale;
    mri2->z_a = mri->z_a / ras_scale;
    mri2->z_s = mri->z_s / ras_scale;

    mri2->c_r = mri->c_r;
    mri2->c_a = mri->c_a;
    mri2->c_s = mri->c_s;

    mri2->ras_good_flag = 1;

  }

  return(mri2);

} /* end conform_voxels() */

MRI *conform_size(MRI *mri)
{

  MRI *mri2;
  int pre_x, pre_y, pre_z;
  int i, j;

  if(mri->slices == NULL)
  {
    mri2 = MRIcopy(mri, NULL);
    mri2->width = mri2->height = mri2->depth = 256;
    return(mri2);
  }
  else
    mri2 = MRIalloc(256, 256, 256, mri->type);

  MRIcopyHeader(mri, mri2);

  pre_x = (int)((256 - mri->width) / 2);
  pre_y = (int)((256 - mri->height) / 2);
  pre_z = (int)((256 - mri->depth) / 2);

  for(i = 0;i < 256;i++)
    for(j = 0;j < 256;j++)
      memset(mri2->slices[i][j], 0x00, 256);

 for(i = 0;i < mri->depth;i++)
    for(j = 0;j < mri->height;j++)
      memcpy(&(mri2->slices[i+pre_z][j+pre_y][pre_x]), &(mri->slices[i][j][0]), mri->width);

  return(mri2);

} /* end conform_size() */

MRI *conform_direction(MRI *mri)
{

  MRI *mri2;

  if(mri->slices == NULL)
    mri2 = MRIcopy(mri, NULL);
  else
  {
    mri2 = MRIreorder(mri, NULL, mri->xdir, mri->ydir, mri->zdir);
  }

  mri2->xdir = XDIM;
  mri2->ydir = YDIM;
  mri2->zdir = ZDIM;
  mri2->slice_direction = MRI_CORONAL;

  return(mri2);

} /* end conform_direction() */

/*  EOF  */
