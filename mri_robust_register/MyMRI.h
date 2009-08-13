/**
 * @file  MyMRI.h
 * @brief A class for MRI utils as used for registration
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/08/13 02:51:19 $
 *    $Revision: 1.2 $
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
//
// written by Martin Reuter
// Aug. 12th ,2009
//

#ifndef MyMRI_H
#define MyMRI_H

#ifdef __cplusplus
extern "C"
{
#endif
#include "matrix.h"
#include "mri.h"
#ifdef __cplusplus
}
#endif

#include <utility>
#include <string>
#include <vector>

class MyMRI
{
public:
  static MRI* makeConform(MRI *mri, MRI *out, bool fixvoxel = true, bool fixtype = true);

  static MRI *  MRIvalscale(MRI *mri_src, MRI *mri_dst, double s);
  static MRI * convolute(MRI * mri, MRI * filter, int dir);
  static MRI * getPrefilter();
  static MRI * getDerfilter();
  static MRI * subSample(MRI * mri);
  static MRI * getBlur(MRI* mriS);
  static MRI * getPartial(MRI* mriS, int dir);
  static bool  getPartials(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz, MRI* &outblur);
  static MRI * getBlur2(MRI* mri);
  static bool  getPartials2(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz, MRI* &outblur);

  static int findRightSize(MRI *mri, float conform_size);

  static MATRIX* MRIgetZslice(MRI * mri, int slice);
};


#endif
