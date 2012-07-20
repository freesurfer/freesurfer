/**
 * @file  MyMRI.h
 * @brief A class for MRI utils as used for registration
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2012/07/20 23:31:17 $
 *    $Revision: 1.12 $
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
  static MRI * MRIvalscale(MRI *mri_src, MRI *mri_dst, double s, double b = 0.0);
  static MRI * MRInorm255(MRI *mri_src, MRI *mri_dst);
  static bool  getPartials(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz, MRI* &outblur);
  static MRI * getBlur(MRI* mriS, MRI* mriT);
  static MRI * subSample(MRI * mri_src, MRI * mri_dst = NULL, bool fixheader = false);
  static MRI * entropyImage( MRI * mri, int radius, int sigma ); 
  static MRI * entropyImage( MRI * mri, int radius, bool ball = false ); 
  static MRI * nlsdImage(MRI * mri, int prad, int nrad);
  static double noiseVar(MRI * mri);

  static vnl_matrix_fixed < double, 4, 4 > MRIvoxelXformToRasXform(MRI * src, MRI * trg, const vnl_matrix_fixed < double, 4, 4> &vox);
	static MRI* MRIlinearTransform(MRI* mriS, MRI* mriT, const vnl_matrix_fixed < double, 4, 4 >& m);

	static bool isIsotropic(MRI *mri);
  static std::vector < int > findRightSize(MRI *mri, float conform_size, bool conform);
  
  static MRI * gaussianCube(int size);
  
private:
  static MRI * getPrefilter();
  static MRI * getDerfilter();
	
	// these are not really used any longer (maybe delete)?
  static MRI * getPartial(MRI* mriS, int dir);
//  static MRI * getBlur2(MRI* mri);
//  static bool  getPartials2(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz, MRI* &outblur);
  static MRI * convolute(MRI * mri, MRI * filter, int dir);

	static bool isConform(MRI *mri);

  static MATRIX* MRIgetZslice(MRI * mri, int slice);
  static double entropyPatch(MRI * mri, int x, int y, int z, int radius, int nbins, MRI* kernel,bool ball = false);
	
};


#endif
