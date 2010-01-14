/**
 * @file  MyMRI.h
 * @brief A class for MRI utils as used for registration
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/01/14 19:41:04 $
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
#include <vnl/vnl_matrix_fixed.h>

class MyMRI
{
public:
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

  static std::vector < int > findRightSize(MRI *mri, float conform_size, bool conform);
	static bool isConform(MRI *mri);
	static bool isIsotropic(MRI *mri);

  static MATRIX* MRIgetZslice(MRI * mri, int slice);
	
  static vnl_matrix_fixed < double, 4, 4 > MRIvoxelXformToRasXform(MRI * src, MRI * trg, const vnl_matrix_fixed < double, 4, 4> &vox);
	static MRI* MRIlinearTransform(MRI* mriS, MRI* mriT, const vnl_matrix_fixed < double, 4, 4 >& m);
};


#endif
