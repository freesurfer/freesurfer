/**
 * @file MultiRegistration.h
 * @brief A class to handle registration of multiple files
 *
 * MultiRegistration is a class to compute a robust registration
 *  of several images. It makes use routines from Registration 
 *
 * written by Martin Reuter
 *  Aug. 12th ,2009
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2012/09/21 23:05:15 $
 *    $Revision: 1.20 $
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

#ifndef MultiRegistration_H
#define MultiRegistration_H

#include <utility>
#include <string>
#include <vector>
//#include <iostream>

#include "RegRobust.h"

#ifdef __cplusplus
extern "C"
{
#endif
#include "matrix.h"
#include "mri.h"
#include "mriBSpline.h"
#include "transform.h"
#ifdef __cplusplus
}
#endif

/** \class MultiRegistration
 * \brief Class for co-registering several images (same modality)
 */
class MultiRegistration
{
public:
  MultiRegistration() :
      outdir("./"), transonly(false), rigid(true), robust(true), sat(4.685), satit(
          false), debug(0), iscale(false), subsamplesize(-1), highit(-1), fixvoxel(
          false), keeptype(false), average(1), doubleprec(false), backupweights(
          false), sampletype(SAMPLE_CUBIC_BSPLINE), mri_mean(NULL)
  {
  }
  ;

  MultiRegistration(const std::vector<std::string> mov) :
      outdir("./"), transonly(false), rigid(true), robust(true), sat(4.685), satit(
          false), debug(0), iscale(false), subsamplesize(-1), highit(-1), fixvoxel(
          false), keeptype(false), average(1), doubleprec(false), backupweights(
          false), sampletype(SAMPLE_CUBIC_BSPLINE), mri_mean(NULL)
  {
    loadMovables(mov);
  }
  ;

  ~MultiRegistration()
  {
    clear();
  }
  ;

  //! Initialize co-registration based on tpi
  bool initialXforms(int tpi, bool fixtp, int regmaxres, int regitmax,
      double regeps);
  //! Iteratively estimate template in mid-space
  bool computeTemplate(int avitmax, double aveps, int regitmax, double regeps);
  //! Half way template (special case for two inputs)
  bool halfWayTemplate(int regmaxres, int regitmax, double regeps,
      bool vox2vox);

  //! Write the tempalate image
  bool writeMean(const std::string& mean);
  //! Write conform tempalate image
  bool writeConformMean(const std::string& cmean);
  //! Write all LTAs
  bool writeLTAs(const std::vector<std::string> & nltas, bool vox2vox,
      const std::string & mean);
  //! Write all mapped movables
  bool writeWarps(const std::vector<std::string>& nwarps);
  //! Write all intensity scales
  bool writeIntensities(const std::vector<std::string>& nintens);
  //! Write all weights
  bool writeWeights(const std::vector<std::string>& nweights,
      bool oneminusweights);

  //! Load all inputs
  int loadMovables(const std::vector<std::string> pmov);
  //! Load initial transforms
  int loadLTAs(const std::vector<std::string> nltas);
  //! Load initial intensity scales
  int loadIntensities(const std::vector<std::string> nintens);
  //! Clear everything
  void clear();

  //! Get a random seed from input images
  unsigned int getSeed();

  //! Set output directory
  void setOutdir(const std::string & s)
  {
    outdir = s;
  }
  ;
  //! Restrict to translation only
  void setTransonly(bool r)
  {
    transonly = r;
  }
  ;
  //! Run rigid registration
  void setRigid(bool r)
  {
    rigid = r;
  }
  ;
  //! Toggle robustness
  void setRobust(bool r)
  {
    robust = r;
  }
  ;
  //! Specify saturation for outlier sensitivity
  void setSaturation(double d)
  {
    sat = d;
  }
  ;
  //! Switch on automatic saturation estimation
  void setSatit(bool b)
  {
    satit = b;
  }
  ;
  //! Set debug level
  void setDebug(int d)
  {
    debug = d;
  }
  ;
  //! Toggle global intensity scaling
  void setIscale(bool i)
  {
    iscale = i;
  }
  ;
  //! Toggle fixed voxel???
  void setFixVoxel(bool i)
  {
    fixvoxel = i;
  }
  ;
  //! Specify if we keep input type 
  void setKeepType(bool i)
  {
    keeptype = i;
  }
  ;
  //! Specify method to create average (1 mean, 2 median..)
  void setAverage(int i)
  {
    average = i;
  }
  ;
  //! Specify size to start subsampling
  void setSubsamplesize(int sss)
  {
    subsamplesize = sss;
  }
  ;
  //! Specify iteration number on highest resolution
  void setHighit(int hit)
  {
    highit = hit;
  }
  ;
  //! Specify precision for registration
  void setDoublePrec(bool b)
  {
    doubleprec = b;
  }
  ;
  //! Specify if weights are keept
  void setBackupWeights(bool b)
  {
    backupweights = b;
  }
  ;

  //! Sample type when creating averages
  void setSampleType(int st)
  {
    switch (st)
    {
    case SAMPLE_TRILINEAR:
    case SAMPLE_CUBIC_BSPLINE:
    case SAMPLE_NEAREST:
      break;
    default:
      std::cout << "ERROR MultiRegistration:setSampleType: " << st
          << " not supported type!" << std::endl;
      exit(1);
    }
    sampletype = st;
  }
  //! Get the sample type
  int getSampleType()
  {
    return sampletype;
  }
  ;

  //! Maps mov based on ltas (also iscale) and then averages them
  bool mapAndAverageMov(int itdebug);

  //! (not tested)
  MRI * averageConformSet(int itdebug = 0);

  //! Creates ltas based on centroids and maps and averages       
  bool initialAverageSet();

  //! Averages a set of images (assumed to be aligned)
  static MRI* averageSet(const std::vector<MRI *>& set, MRI* mean, int method,
      double sat);

private:

  void normalizeIntensities(void);

  void initRegistration(RegRobust & R);

  // copy of input filenames
  std::vector<std::string> mov;
  std::vector<std::string> iltas;
  std::vector<std::string> iintens;

  // Parameter:
  std::string outdir;
  bool transonly;
  bool rigid;
  bool robust;
  double sat;
  bool satit;
  int debug;
  bool iscale;
  int subsamplesize;
  int highit;

  bool fixvoxel;
  bool keeptype;
  int average;
  bool doubleprec;
  bool backupweights;
  int sampletype;

  // DATA
  std::vector<MRI*> mri_mov;
  std::vector<MRI_BSPLINE*> mri_bsplines;
  std::vector<LTA*> ltas;
  std::vector<MRI*> mri_warps;
  std::vector<MRI*> mri_weights;
  std::vector<double> intensities;
  MRI * mri_mean;

};

#endif
