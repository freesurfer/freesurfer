/**
 * @brief Creates a distance transform.
 *
 */
/*
 * Original Author: Florent Segonne
 *    $Author$
 *    $Date$
 *    $Revision$
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

#include "fastmarching.h"

MRI *MRIextractDistanceMap(MRI *mri_src, MRI *mri_dst, int label, float max_distance, int mode, MRI *mri_mask)
{
  MRI *mri_distance = NULL;

  // int  free_mri = 0 ;

  // if (mri_src->type != MRI_FLOAT)
  // {
  //   mri_src = MRIchangeType(mri_src, MRI_FLOAT, 0, 1, 1) ;
  //   free_mri = 1 ;
  // }

  // make sure that the max distance is greater than 0
  if (max_distance <= 0) {
    max_distance = 2 * MAX(MAX(mri_src->width, mri_src->height), mri_src->depth);
  }

  if (mri_dst == NULL) {
    mri_dst = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
  }

  mri_distance = mri_dst;

  // make sure that the member variables are all the same
  if (mri_dst->width != mri_src->width || mri_dst->height != mri_src->height || mri_dst->depth != mri_src->depth ||
      mri_dst->type != MRI_FLOAT) {
    if (mri_dst->width != mri_src->width)
      fprintf(stderr,
              "ERROR : incompatible structure with mri_src:\n"
              "mri_dst->width=%d != mri_src->width=%d\n",
              mri_dst->width,
              mri_src->width);

    if (mri_dst->height != mri_src->height)
      fprintf(stderr,
              "ERROR : incompatible structure with mri_src:\n"
              "mri_dst->height=%d != mri_src->height=%d\n",
              mri_dst->height,
              mri_src->height);

    if (mri_dst->depth != mri_src->depth)
      fprintf(stderr,
              "ERROR : incompatible structure with mri_src:\n"
              "mri_dst->depth=%d != mri_src->depth=%d\n",
              mri_dst->depth,
              mri_src->depth);

    if (mri_dst->type != MRI_FLOAT)
      fprintf(stderr,
              "ERROR : incompatible structure with mri_dst:\n"
              "mri_dst->type=%d != MRI_FLOAT\n",
              mri_dst->type);
  }
  else {
    // set values to zero
    for (int z = 0; z < mri_dst->depth; z++)
      for (int y = 0; y < mri_dst->height; y++)
        for (int x = 0; x < mri_dst->width; x++) MRIFvox(mri_dst, x, y, z) = 0.0f;

    const int outside = 1;
    const int inside = 2;
    const int both = 3;

    // positive inside and positive outside
    const int bothUnsigned = 4;

    if (mode == outside || mode == both || mode == bothUnsigned) {
      FastMarching< +1 > fastmarching_out(mri_distance, mri_mask);
      fastmarching_out.SetLimit(max_distance);
      fastmarching_out.InitFromMRI(mri_src, label);
      fastmarching_out.Run(max_distance);
    }

    if (mode == inside || mode == both || mode == bothUnsigned) {
      FastMarching< -1 > fastmarching_in(mri_distance, mri_mask);
      fastmarching_in.SetLimit(-max_distance);
      fastmarching_in.InitFromMRI(mri_src, label);
      fastmarching_in.Run(-max_distance);
    }

    // if unsigned, then we want than the inside and out to be positive
    if (mode == bothUnsigned) {
      for (int z = 0; z < mri_distance->depth; z++)
        for (int y = 0; y < mri_distance->height; y++)
          for (int x = 0; x < mri_distance->width; x++)
            MRIFvox(mri_distance, x, y, z) = fabs(MRIFvox(mri_distance, x, y, z));
    }
  }

  return mri_distance;
}
