/**
 * @brief A class for MRI utils as used for registration
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
//
// written by Martin Reuter
// Aug. 12th ,2009
//
#ifndef MyMRI_H
#define MyMRI_H

#include "matrix.h"
#include "mri.h"

#include <utility>
#include <string>
#include <vector>

#define export // obsolete feature 'export template' used in these headers 
#include <vnl/vnl_matrix_fixed.h>
#undef export

/** \class MyMRI
 * \brief Static class for some MRI operations
 */
class MyMRI
{
public:
  //! Scale intensities by s, offset b
  static MRI * MRIvalscale(MRI *mri_src, MRI *mri_dst, double s,
      double b = 0.0);
  //! Normalize intensities to 0..255
  static MRI * MRInorm255(MRI *mri_src, MRI *mri_dst);
  //! Compute partial derivatives and blurr images
  static bool getPartials(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz,
      MRI* &outblur);
  //! Get blurred image
  static MRI * getBlur(MRI* mriS, MRI* mriT);
  //! Subsample the image
  static MRI * subSample(MRI * mri_src, MRI * mri_dst = NULL, bool fixheader =
      true, int randpos = 0);
  //! Get Entropy image (old)
  static MRI * entropyImage(MRI * mri, int radius, int sigma);
  //! Get Entropy image using box or ball
  static MRI * entropyImage(MRI * mri, int radius, bool ball = false,
      bool correction = false, MRI * mask=NULL);
  //! NLSD not tested
  static MRI * nlsdImage(MRI * mri, int prad, int nrad);

  //! subtracts mean and divides by standard deviation (global)
  static MRI * getNormalizedImage(MRI *mri);
  
  //! subtracts local mean and std
  static MRI * getNormalizedImage(MRI *mri, int boxradius);

  //! mean Filter
  static MRI * meanFilter(MRI *mri, int boxradius);

  //! median Filter
  static MRI * medianFilter(MRI *mri, int boxradius);

  static double noiseVar(MRI * mri);

  static vnl_matrix_fixed<double, 4, 4> MRIvoxelXformToRasXform(MRI * src,
      MRI * trg, const vnl_matrix_fixed<double, 4, 4> &vox);
  static MRI* MRIlinearTransform(MRI* mriS, MRI* mriT,
      const vnl_matrix_fixed<double, 4, 4>& m);

  //! Check if image is isotropic
  static bool isIsotropic(MRI *mri);

  static std::vector<int> findRightSize(MRI *mri, float conform_size,
      bool conform);

  //! Gaussian kernel in cube of given size
  static MRI * gaussianCube(int size);

  //! Get (next) random number (uniform 0..1)
  static double getRand(int & randpos);
  
  //! Get background intensity (vox val with largest level set)
  static float getBackground(MRI * mri);
  
  //! Set outside value in header to the max value of image
  static void setMaxOutsideVal(MRI * mri);
  
  //! Change type to UCHAR
  static MRI* setTypeUCHAR(MRI * mri);

private:
  static MRI * getPrefilter();
  static MRI * getDerfilter();

  static bool isConform(MRI *mri);

  static MATRIX* MRIgetZslice(MRI * mri, int slice, int frame = 0);
  static double entropyPatch(MRI * mri, int x, int y, int z, int radius,
      int nbins, MRI* kernel, bool ball = false);
  static void get3Dcorrection(double* histo, unsigned int v1, unsigned int v2,
      unsigned int v3, unsigned int v4, unsigned int intRange);

};

inline double MyMRI::getRand(int & randpos)
// adjusts randpos if > 100
{
  static double uniform[101] =
      { 0.825960894329359, 0.391799656122844, 0.606822280998442,
          0.923192169644637, 0.772662867376799, 0.394913419972908,
          0.562648884384293, 0.821213453107628, 0.077777173551299,
          0.536306917657496, 0.514107239245956, 0.761227511728186,
          0.978425890910315, 0.157030594951506, 0.147496203406248,
          0.783403101626463, 0.637363106003275, 0.257349279177570,
          0.109403726553149, 0.136954287985675, 0.654171911591805,
          0.204476823923606, 0.164700938803686, 0.566376934840801,
          0.854095286894540, 0.430581644078462, 0.397193805113614,
          0.045536084142315, 0.140620297119782, 0.724179303888238,
          0.210184975302471, 0.346103835004574, 0.226683913788046,
          0.308131967172401, 0.451070768412025, 0.831843821825446,
          0.086760750364257, 0.854135129404520, 0.742590231477529,
          0.053858310273945, 0.122791324135168, 0.526967625411493,
          0.320307444448230, 0.520062337463421, 0.878893684505441,
          0.443226585521493, 0.320014638636649, 0.868171615160687,
          0.764138797526944, 0.864611801159897, 0.083244000805187,
          0.340050247951733, 0.581607039757426, 0.291315250124114,
          0.467637373824935, 0.596640965216293, 0.897387480194169,
          0.121159463814373, 0.018178277226624, 0.747666956017887,
          0.670108792666926, 0.092391423049263, 0.917528192569663,
          0.347581829497044, 0.857229173664249, 0.632564501155682,
          0.363593396666436, 0.385827512659221, 0.003355759148503,
          0.969283884272977, 0.459201707491767, 0.800715562489883,
          0.296945204007638, 0.619917383791628, 0.748749776544515,
          0.394902341568457, 0.354479046607957, 0.867388084434840,
          0.627567204750354, 0.984294463747630, 0.824472893444168,
          0.295272747204669, 0.927251152590456, 0.119266587069099,
          0.241265513320426, 0.527775796652046, 0.060203502196487,
          0.835363100840541, 0.148398316485924, 0.121273963075904,
          0.683207184266160, 0.500002874704566, 0.561265939626174,
          0.234256657698203, 0.486854859925734, 0.071141206079346,
          0.442630693187859, 0.327200604299592, 0.226827609433358,
          0.971944076189183, 0.302612670611030 };
  if (randpos > 100)
    randpos = randpos % 101;
  return uniform[randpos];
}

#endif
