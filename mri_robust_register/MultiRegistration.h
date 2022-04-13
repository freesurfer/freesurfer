/**
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

#ifndef MultiRegistration_H
#define MultiRegistration_H

#include <utility>
#include <string>
#include <vector>
//#include <iostream>

#include "RegRobust.h"

#include "matrix.h"
#include "mri.h"
#include "mriBSpline.h"
#include "transform.h"

/** \class MultiRegistration
 * \brief Class for co-registering several images (same modality)
 */
class MultiRegistration
{
public:
  MultiRegistration() :
      outdir("./"), transonly(false), rigid(true), robust(true), sat(4.685),
          satit(false), debug(0), iscale(false), iscaleonly(false),
          nomulti(false), subsamplesize(-1), highit(-1), fixvoxel(false),
          keeptype(false), average(1), doubleprec(false), backupweights(false),
	sampletype(SAMPLE_CUBIC_BSPLINE), crascenter(false), resthresh(0.01), mri_mean(NULL)
  {
  }

  MultiRegistration(const std::vector<std::string> mov) :
      outdir("./"), transonly(false), rigid(true), robust(true), sat(4.685),
          satit(false), debug(0), iscale(false), iscaleonly(false),
          nomulti(false), subsamplesize(-1), highit(-1), fixvoxel(false),
          keeptype(false), average(1), doubleprec(false), backupweights(false),
          sampletype(SAMPLE_CUBIC_BSPLINE), crascenter(false), resthresh(0.01), mri_mean(NULL)
  {
    loadMovables(mov);
  }

  ~MultiRegistration()
  {
    clear();
  }
  
  //! Print parameter settings
  void printParams()
  {
    std::cout << std::boolalpha << std::endl;
    std::cout << " MultiRegistration Parameters " << std::endl << std::endl;
    std::cout << " Outdir:        " << outdir << std::endl;
    std::cout << " TransOnly:     " << transonly << std::endl;
    std::cout << " Rigid:         " << rigid << std::endl;
    std::cout << " Robust:        " << robust << std::endl;
    if (satit)
      std::cout << " Satit:         " << satit << std::endl;
    else
      std::cout << " Sat:           " << sat << std::endl;
    std::cout << " Iscale:        " << iscale<< std::endl;
    std::cout << " IscaleOnly:    " << iscaleonly << std::endl;
    std::cout << " NoMulti:       " << nomulti << std::endl;
    std::cout << " SubsampleSize: " << subsamplesize << std::endl;
    std::cout << " HighIt:        " << highit << std::endl;
    std::cout << " FixVoxel:      " << fixvoxel << std::endl;
    std::cout << " KeepType:      " << keeptype << std::endl;
    std::cout << " Average:       " << average << std::endl;
    std::cout << " DoublePrec:    " << doubleprec << std::endl;
    std::cout << " BackupWeights: " << backupweights << std::endl;
    std::cout << " SampleType:    " << sampletype<< std::endl;
    std::cout << " CRASCenter:    " << crascenter<< std::endl;
    std::cout << " Resthresh:     " << resthresh << std::endl;
    std::cout << " Debug:         " << debug << std::endl;
    std::cout <<  std::noboolalpha << std::endl;
  
  }
  

  //! Initialize co-registration based on tpi
  bool initialXforms(int tpi, bool fixtp, int regmaxres, int regitmax, double regeps);
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
  //! Write header-adjusted movables.
  bool writeMapMovHdr(const std::vector<std::string>& mapmovhdr);
  //! Write all intensity scales
  bool writeIntensities(const std::vector<std::string>& nintens);
  //! Write all weights
  bool writeWeights(const std::vector<std::string>& nweights,
      bool oneminusweights);

  //! Load all inputs
  int loadMovables(const std::vector<std::string> &pmov, const std::vector<std::string>& masks = std::vector<std::string>());
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

  //! Restrict to translation only
  void setTransonly(bool r)
  {
    transonly = r;
  }

  //! Run rigid registration
  void setRigid(bool r)
  {
    rigid = r;
  }

  //! Toggle robustness
  void setRobust(bool r)
  {
    robust = r;
  }

  //! Specify saturation for outlier sensitivity
  void setSaturation(double d)
  {
    sat = d;
  }

  //! Switch on automatic saturation estimation
  void setSatit(bool b)
  {
    satit = b;
  }

  //! Set debug level
  void setDebug(int d)
  {
    debug = d;
  }

  //! Toggle global intensity scaling
  void setIscale(bool i)
  {
    iscale = i;
  }

  //! Toggle only intensity scaling (no spacial transform)
  void setIscaleOnly(bool i)
  {
    iscaleonly = i;
  }

  //! Toggle multi-res or only high res
  void setNoMulti(bool i)
  {
    nomulti = i;
  }

  //! Toggle fixed voxel???
  void setFixVoxel(bool i)
  {
    fixvoxel = i;
  }

  //! Specify if we keep input type 
  void setKeepType(bool i)
  {
    keeptype = i;
  }

  //! Specify method to create average (1 mean, 2 median..)
  void setAverage(int i)
  {
    average = i;
  }

  //! Specify size to start subsampling
  void setSubsamplesize(int sss)
  {
    subsamplesize = sss;
  }

  //! Specify iteration number on highest resolution
  void setHighit(int hit)
  {
    highit = hit;
  }

  //! Specify precision for registration
  void setDoublePrec(bool b)
  {
    doubleprec = b;
  }

  //! Specify if weights are keept
  void setBackupWeights(bool b)
  {
    backupweights = b;
  }

  //! Specify voxel threshold, default is 0.001
  void setResthresh(float thresh)
  {
    resthresh = thresh;
  }

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

  //! If true: center at average CRAS, else at avg. barycenter
  void useCRAS(bool b)
  {
    crascenter=b;
  }

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

  vnl_matrix_fixed<double, 3, 3> getAverageCosines();
  MRI * createTemplateGeo();

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
  bool iscaleonly;
  bool nomulti;
  int subsamplesize;
  int highit;

  bool fixvoxel;
  bool keeptype;
  int average;
  bool doubleprec;
  bool backupweights;
  int sampletype;
  bool crascenter;
  float resthresh;

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
