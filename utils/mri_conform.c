#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "mri.h"
#include "error.h"
#include "mri_conform.h"

/*  if the pixel data is sorted by value, these values tell what
  fraction of the pixels are to be clipped to the maximum and
  minimum new pixel values -- the fraction of pixels between
  the two are scaled linearly
*/
#define  BOTTOM_CLIP  0.001
#define TOP_CLIP  0.999


static MRI *interpolate_and_pad(MRI *mri);
static short s_partial_quicksort(short list[], int n, float cutoff_fraction);
static unsigned char *s_scale(short list[], int list_length);
static float f_partial_quicksort(float list[], int n, float cutoff_fraction);
static unsigned char *f_scale(float list[], int list_length);

MRI *MRIconform(MRI *mri, void *p_data, int xdim, int ydim, int zdim)
{

  int n_pixels = mri->width * mri->height * mri->depth;
  float bottom_cutoff, top_cutoff;
  unsigned char *uchar_list;
  int i, j, k;
  MRI *mri2, *mri3, *mri4;

  if(p_data == NULL)
    {
    mri2 = MRIallocHeader(mri->width, mri->height, 256, MRI_UCHAR);
    MRIcopyHeader(mri, mri2);
    mri2->depth = 256;
    mri2->type = MRI_UCHAR;

    if(mri2 == NULL)
      ErrorReturn(NULL, (ERROR_NO_MEMORY, "MRIconform: can't allocate sequence"));

    mri2->zsize = mri2->xsize;
    mri2->zstart = mri2->xstart;
    mri2->zend = mri2->xend;

    mri2->imnr1 = 256;

    mri2->slice_direction = MRI_CORONAL;
    return(mri2);

    }
  else
    {

    if(mri->type == MRI_FLOAT)
      uchar_list = f_scale((float *)p_data, n_pixels);
    else if(mri->type == MRI_SHORT)
      uchar_list = s_scale((short *)p_data, n_pixels);
    else
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIconform: unsupported data type %d\n", mri->type));

    mri2 = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_UCHAR, 1);
    if(mri2 == NULL)
      ErrorReturn(NULL, (ERROR_NO_MEMORY, "MRIconform: can't allocate sequence"));
    MRIcopyHeader(mri, mri2);
    mri2->depth = mri->depth;

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
          MRIvox(mri2, i, j, k) = uchar_list[i + mri->width * (j + k*mri->height)];

    mri3 = interpolate_and_pad(mri2);
    MRIfree(&mri2);

    mri4 = MRIreorder(mri3, NULL, xdim, ydim, zdim);
    mri4->slice_direction = MRI_CORONAL;
    MRIfree(&mri3);

    return(mri4);
    }

}  /*  end MRIconform()  */

static MRI *interpolate_and_pad(MRI *mri)
{

  MRI *mri2;
  int i, j, k;
  float findex, fifract;
  int fiint;

  mri2 = MRIallocSequence(mri->width, mri->height, 256, MRI_UCHAR, 1);
  MRIcopyHeader(mri, mri2);
  mri2->zsize = mri2->xsize;
  mri2->zstart = mri2->xstart;
  mri2->zend = mri2->xend;
  mri2->imnr1 = 256;

  for(i = 0;i < 256;i++)
    {
    findex = (i-128) * mri2->zsize / mri->zsize + mri->depth / 2;
    fiint = (int)(floor(findex));
    fifract = findex - floor(findex);
    if(fiint < 0 || fiint > mri->depth - 2)
      {
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->width;k++)
          MRIvox(mri2, k, j, i) = 0;
      }
    else
      {
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->width;k++)
          MRIvox(mri2, k, j, i) = (1. - fifract) * (float)MRIvox(mri, k, j, fiint) + fifract * (float)MRIvox(mri, k, j, fiint+1);
      }
    }

  return(mri2);

}  /*  end interpolate_and_pad()  */

/*  quicksort algorithm, sorting only the partition needed to find the
  (n*cutoff_fraction)th value
  see Sedgewick, Algorithms in C, pp. 126-130
*/
static short s_partial_quicksort(short list[], int n, float cutoff_fraction)
{

  int cutoff_index = (int)(cutoff_fraction * n);

  int i, j, l, r;
  short v, t;

  l = 1;
  r = n-1;

  while(r > l)
    {
    v = list[r]; i = l-1; j = r;
    for(;;)
      {
      while(list[++i] < v);
      while(list[--j] > v);
      if(i >= j)
        break;
      t = list[i]; list[i] = list[j]; list[j] = t;
      }
    t = list[i]; list[i] = list[r]; list[r] = t;
    if(i >= cutoff_index)
      r = i-1;
    if(i <= cutoff_index)
      l = i+1;
    }

  return(list[l-1]);


}  /*  end  s_partial_quicksort()  */

static unsigned char *s_scale(short list[], int list_length)
{

  unsigned char *new_list;
  int i;
  short scale_bottom, scale_top;
  float scale_factor;
  short *list_copy;

  list_copy = (short *)malloc(list_length * 2);
  memcpy(list_copy, list, list_length * 2);
  scale_bottom = s_partial_quicksort(list_copy, list_length, BOTTOM_CLIP);
  scale_top = s_partial_quicksort(list_copy, list_length, TOP_CLIP);
  free(list_copy);
  scale_factor = 255 / (float)(scale_top - scale_bottom);

  new_list = (unsigned char *)malloc(list_length);

  for(i = 0;i < list_length;i++)
    {
    if(list[i] > scale_top)
      list[i] = scale_top;
    if(list[i] < scale_bottom)
      list[i] = scale_bottom;
    new_list[i] = (unsigned char)((list[i] - scale_bottom) * scale_factor);
    }  
  return(new_list);

}  /*  end s_scale()  */

/*  quicksort algorithm, sorting only the partition needed to find the
  (n*cutoff_fraction)th value
  see Sedgewick, Algorithms in C, pp. 126-130
*/
static float f_partial_quicksort(float list[], int n, float cutoff_fraction)
{

  int cutoff_index = (int)(cutoff_fraction * n);

  int i, j, l, r;
  float v, t;

  l = 1;
  r = n-1;

  while(r > l)
    {
    v = list[r]; i = l-1; j = r;
    for(;;)
      {
      while(list[++i] < v);
      while(list[--j] > v);
      if(i >= j)
        break;
      t = list[i]; list[i] = list[j]; list[j] = t;
      }
    t = list[i]; list[i] = list[r]; list[r] = t;
    if(i >= cutoff_index)
      r = i-1;
    if(i <= cutoff_index)
      l = i+1;
    }

  return(list[l-1]);


}  /*  end  f_partial_quicksort()  */

static unsigned char *f_scale(float list[], int list_length)
{

  unsigned char *new_list;
  int i;
  float scale_bottom, scale_top, scale_factor;
  float *list_copy;

  list_copy = (float *)malloc(list_length * 4);
  memcpy(list_copy, list, list_length * 4);
  scale_bottom = f_partial_quicksort(list_copy, list_length, BOTTOM_CLIP);
  scale_top = f_partial_quicksort(list_copy, list_length, TOP_CLIP);
  free(list_copy);
  scale_factor = 255 / (scale_top - scale_bottom);

  new_list = (unsigned char *)malloc(list_length);

  for(i = 0;i < list_length;i++)
    {
    if(list[i] > scale_top)
      list[i] = scale_top;
    if(list[i] < scale_bottom)
      list[i] = scale_bottom;
    new_list[i] = (unsigned char)((list[i] - scale_bottom) * scale_factor);
    }
  
  return(new_list);

}  /*  end f_scale()  */

/*  EOF  */
