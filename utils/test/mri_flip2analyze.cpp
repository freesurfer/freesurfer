/*
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


/**
 * @author Yasunari Tosa
 * @date   Wed Dec 22 13:32:49 2004
 *
 * @brief  flip axis direction to analyze orient 0 direction
 *
 *
 */

#include <iostream>
#include <iomanip>

extern "C"
{
#include "mri.h"

  char *Progname = "mri_flip2analyze";
}

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    cout << "Usage: mri_flip2analyze <srcvol> <dstvol>" << endl;
    cout << "   if you want analyze format, specify the destination with .img extension" << endl;
    return 0;
  }

  MRI *src = MRIread(argv[1]);
  ///////////////////////////////////////////////////////////////////////////////////////
  // decide which one of the six ways direction cosines look like
  //     0. a 0 0   1. a 0 0  2. 0 0 c  3. 0 b 0  4. 0 b 0  5. 0 0 c
  //        0 b 0      0 0 c     a 0 0     a 0 0     0 0 c     0 b 0
  //        0 0 c      0 b 0     0 b 0     0 0 c     a 0 0     a 0 0  where a, b, c = +/- 1
  //
  // analyze orient = 0 is  LAS
  //       -1 0 0
  //        0 1 0
  //        0 0 1
  /////////////////////////////////////////////////////////////////////////////////////

  int which_xras = ( fabs(src->x_r) > fabs(src->x_a) ) ?
                   ( ( fabs(src->x_r) > fabs(src->x_s) ) ?  0 : 2 ) :  ( (fabs(src->x_a ) > fabs(src->x_s)) ? 1 : 2);
  int which_yras = ( fabs(src->y_r) > fabs(src->y_a) ) ?
                   ( ( fabs(src->y_r) > fabs(src->y_s) ) ?  0 : 2 ) :  ( (fabs(src->y_a ) > fabs(src->y_s)) ? 1 : 2) ;

  float dstXsize = 0.;
  float dstYsize = 0.;
  float dstZsize = 0.;

  // fix axis size in the destination
  int index = 0;
  int dstWidth = 0;
  int dstHeight = 0;
  int dstDepth = 0;
  switch (which_xras)
  {
  case 0:
      dstXsize = src->xsize;
    if (which_yras == 1)
    {
      index = 0;
      dstYsize = src->ysize;
      dstZsize = src->zsize;
      dstWidth = src->width;
      dstHeight = src->height;
      dstDepth = src->depth;
    }
    else
    {
      index = 1;
      dstYsize = src->zsize;
      dstZsize = src->ysize;
      dstWidth = src->width;
      dstHeight = src->depth;
      dstDepth = src->height;
    }
    break;
  case 1:
    dstYsize = src->xsize;
    if (which_yras == 2)
    {
      index = 2;
      dstXsize = src->zsize;
      dstZsize = src->ysize;
      dstWidth = src->depth;
      dstHeight = src->width;
      dstDepth = src->height;
    }
    else
    {
      index = 3;
      dstXsize = src->ysize;
      dstZsize = src->zsize;
      dstWidth = src->height;
      dstHeight = src->width;
      dstDepth = src->depth;
    }
    break;
  case 2:
    dstZsize = src->xsize;
    if (which_yras == 0)
    {
      index = 4;
      dstXsize = src->ysize;
      dstYsize = src->zsize;
      dstWidth = src->height;
      dstHeight = src->depth;
      dstDepth = src->width;
    }
    else
    {
      index = 5;
      dstXsize = src->zsize;
      dstYsize = src->ysize;
      dstWidth = src->depth;
      dstHeight = src->height;
      dstDepth = src->width;
    }
    break;
  }
  MRI *dst = MRIallocSequence(dstWidth, dstHeight, dstDepth, src->type, src->nframes);
  dst->xsize = dstXsize;
  dst->ysize = dstYsize;
  dst->zsize = dstZsize;
  //set to Analyze orient = 0 orientation
  dst->x_r = -1;
  dst->x_a = 0;
  dst->x_s = 0;
  dst->y_r = 0;
  dst->y_a = 1;
  dst->y_s = 0;
  dst->z_r = 0;
  dst->z_a = 0;
  dst->z_s = 1;
  dst->c_r = 0;
  dst->c_a = 0;
  dst->c_s = 0;
  // now set the voxel values
  int dstx;
  int dsty;
  int dstz;
  for (int frame=0; frame < src->nframes; frame++)
    for (int z = 0; z < src->depth; z++)
      for (int y = 0; y < src->height; y++)
        for (int x = 0; x < src->width; x++)
        {
          switch (index)
          {
            ///////////////////////////////////////////
          case 0:
            // a 0 0
            // 0 b 0
            // 0 0 c
            dstx = (src->x_r > 0.) ? (src->width-x-1): x;   // X
            dsty = (src->y_a > 0.) ? y : (src->height-y-1); // Y
            dstz = (src->z_s > 0.) ? z : (src->depth-z-1);  // Z
            break;
          case 1:
            // a 0 0
            // 0 0 c
            // 0 b 0
            dstx = (src->x_r > 0.) ? (src->width-x-1): x;   // X
            dsty = (src->z_a > 0.) ? z : (src->depth-z-1);  // Z'
            dstz = (src->y_s > 0.) ? y : (src->height-y-1); // Y'
            break;
            ////////////////////////////////////////////
          case 2:
            // 0 0 c
            // a 0 0
            // 0 b 0
            dstx = (src->z_r > 0.) ? (src->depth-z-1): z;   // Z"
            dsty = (src->x_a > 0.) ? x : (src->width-x-1);  // X'
            dstz = (src->y_s > 0.) ? y : (src->height-y-1); // Y'
            break;
          case 3:
            // 0 b 0
            // a 0 0
            // 0 0 c
            dstx = (src->y_r > 0.) ? (src->height-y-1): y; // Y"
            dsty = (src->x_a > 0.) ? x : (src->width-x-1); // X'
            dstz = (src->z_s > 0.) ? z : (src->depth-z-1); // Z
            break;
            ////////////////////////////////////////////
          case 4:
            // 0 b 0
            // 0 0 c
            // a 0 0
            dstx = (src->y_r > 0.) ? (src->height-y-1) : y; // Y"
            dsty = (src->z_a > 0.) ? z : (src->depth-z-1);  // Z'
            dstz = (src->x_s > 0.) ? x : (src->width-x-1);  // X"
            break;
          case 5:
            // 0 0 c
            // 0 b 0
            // a 0 0
            dstx = (src->z_r > 0.) ? (src->depth-z-1) : z;  // Z"
            dsty = (src->y_a > 0.) ? y : (src->height-y-1); // Y
            dstz = (src->x_s > 0.) ? x : (src->width-x -1); // X"
            break;
          default:
            cerr << "NO such case is allowed. index is out of range" << endl;
            return -1;
          }
          MRIsetVoxVal(dst, dstx, dsty, dstz, frame, MRIgetVoxVal(src, x, y, z, frame));
        }

  MRIwrite(dst, argv[2]);
  MRIfree(&src);
  MRIfree(&dst);
}
