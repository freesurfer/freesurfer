#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "mri.h"
#include "error.h"
#include "histo.h"
#include "mri_conform.h"

MRI *conform_type(MRI *mri);
MRI *conform_voxels(MRI *mri);
MRI *conform_direction(MRI *mri);

MRI *MRIconform(MRI *mri)
{

  MRI *mri2, *mri3, *mri4;

  mri3 = conform_type(mri);

  mri4 = conform_voxels(mri3);
  MRIfree(&mri3);
  mri2 = conform_direction(mri4);
  MRIfree(&mri4);

  return(mri2);

}  /*  end MRIconform()  */

#define N_BINS 256
#define FRACTION_SCALE 0.01

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

  if(mri->type == MRI_UCHAR)
  {
    mri2 = MRIcopy(mri, NULL);
    return(mri2);
  }

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
        printf("bin < 0\n");
      if(bin > N_BINS - 1)
        printf("bin > N_BINS - 1\n");

      counts[bin]++;

      }

  nv_needed = (int)(mri->height * mri->width * mri->depth * FRACTION_SCALE);

  nv = 0;
  for(i = 0;i < N_BINS && nv < nv_needed;i++)
    nv += counts[i];

  if(i == -1)
    ErrorExit(ERROR_BADPARM, "MRIconform (value scaling): histogram is too thin for\na clipping fraction of %g",
              FRACTION_SCALE);

  e1 = (i-1) * bin_size + min;

  nv = 0;
  for(i = N_BINS-1;i >= 0 && nv < nv_needed;i--)
    nv += counts[i];

  if(i == N_BINS)
    ErrorExit(ERROR_BADPARM, "MRIconform (value scaling): histogram is too thin for\na clipping fraction of %g",
              FRACTION_SCALE);

  e2 = (i+1) * bin_size + min;

  scale = 255.0 / (e2 - e1);

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

      if(this <= e1)
        MRIvox(mri2, x, y, z) = 0;
      else if(this >= e2)
        MRIvox(mri2, x, y, z) = 255;
      else
        MRIvox(mri2, x, y, z) = (unsigned char)(scale * (this - e1));

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

  fovx = mri->width * mri->xsize;
  fovy = mri->height * mri->ysize;
  fovz = mri->depth * mri->zsize;

  min_size = 1;
  max_fov = 256.;
/*
  max_fov = (fovx > fovy ? fovx : fovy);
  max_fov = (fovz > max_fov ? fovx : max_fov);
  min_size = max_fov / 256.0;
*/
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

  return(mri2);

} /* end conform_voxels() */

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
