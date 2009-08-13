//
// MyMRI is a class for several MRI operations
//    as used for registration
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
