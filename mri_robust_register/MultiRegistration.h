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
 *    $Date: 2011/10/07 22:28:51 $
 *    $Revision: 1.15 $
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

#include "Registration.h"

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

class MultiRegistration
{
public:
  MultiRegistration():outdir("./"),transonly(false),rigid(true),robust(true),sat(4.685),satit(false),
	                     debug(0),iscale(false),subsamplesize(-1),highit(-1),fixvoxel(false),
											 keeptype(false),average(1),doubleprec(false),sampletype(SAMPLE_CUBIC_BSPLINE),
                       mri_mean(NULL)
		{};
  MultiRegistration(const std::vector < std::string > mov):outdir("./"),transonly(false),
	                     rigid(true),robust(true),sat(4.685),satit(false),debug(0),iscale(false),
											 subsamplesize(-1),highit(-1),fixvoxel(false),keeptype(false),average(1),doubleprec(false),
											 sampletype(SAMPLE_CUBIC_BSPLINE),mri_mean(NULL)
  { loadMovables(mov);};
		
  ~MultiRegistration()
  {clear();};
	 
  bool initialXforms(int tpi, bool fixtp, int regmaxres, int regitmax, double regeps);
  bool computeTemplate(int avitmax, double aveps, int regitmax, double regeps);
  bool halfWayTemplate(int regmaxres, int regitmax, double regeps, bool vox2vox);

  bool writeMean(const std::string& mean);
  bool writeConformMean(const std::string& cmean);
  bool writeLTAs(const std::vector < std::string > & nltas, bool vox2vox, const std::string & mean);
  bool writeWarps(const std::vector <  std::string >& nwarps);
  bool writeIntensities(const std::vector < std::string >& nintens);
  bool writeWeights(const std::vector < std::string >& nweights, bool oneminusweights);


  int loadMovables(const std::vector < std::string > mov);
  int loadLTAs(const std::vector < std::string > nltas);
  int loadIntensities(const std::vector < std::string > nintens);
  void clear();
  
  unsigned int getSeed();
	 
  // Set parameters:
  void setOutdir(const std::string & s)
  {
    outdir = s;
  }; 
  void setTransonly(bool r)
  {
    transonly = r;
  };
  void setRigid(bool r)
  {
    rigid = r;
  };
  void setRobust(bool r)
  {
    robust = r;
  };
  void setSaturation(double d)
  {
    sat = d;
  };
  void setSatit(bool b)
  {
    satit = b;
  };
  void setDebug(int d)
  {
    debug = d;
  };
  void setIscale(bool i)
  {
    iscale = i;
  };
  void setFixVoxel(bool i)
  {
    fixvoxel = i;
  };
  void setKeepType(bool i)
  {
    keeptype = i;
  };
  void setAverage(int i)
  {
    average = i;
  };
  void setSubsamplesize (int sss)
  {
    subsamplesize = sss;
  };
  void setHighit (int hit)
  {
    highit = hit;
  };
  void setDoublePrec(bool b)
  {
    doubleprec = b;
  }
  // sample type when creating averages
  void setSampleType(int st)
  {
    switch (st)
    {
    case SAMPLE_TRILINEAR:
    case SAMPLE_CUBIC_BSPLINE:
      break;
    default:
      std::cout << "ERROR Registration:setSampleType: " << st << " not supported type!" << std::endl;
      exit(1);
    }
    sampletype = st;
  }
  int getSampleType() {return sampletype;};
  
  // maps mov based on ltas (also iscale) and averages:
  bool mapAndAverageMov(int itdebug);
  
  MRI * averageConformSet(int itdebug = 0);
	
	// creates ltas based on centroids and maps and averages			 
  bool initialAverageSet();

  // averages a set of images (assumed to be aligned)
  static MRI* averageSet(const std::vector < MRI * >& set,
                       MRI* mean, int method, double sat);


private:

  void normalizeIntensities(void);

  void initRegistration(Registration & R);
  
  // copy of input filenames
  std::vector <std::string > mov;
  std::vector <std::string > iltas;
  std::vector <std::string > iintens;
	
  // copy of output filenames
  //std::string mean;
  // std::vector <std::string> nltas;
  //std::vector <std::string> nweights;
  // std::vector <std::string> nwarps;
	
  // Parameter:
  std::string outdir;
  bool   transonly;
  bool   rigid;
  bool   robust;
  double sat;
  bool   satit;
  int    debug;
  bool   iscale;
  int    subsamplesize;
  int    highit;
	
  bool   fixvoxel;
  bool   keeptype;
  int    average;
  bool   doubleprec;
  int    sampletype;
	
  // DATA
  std::vector < MRI* > mri_mov;
  std::vector < MRI_BSPLINE* > mri_bsplines;
  std::vector < LTA* > ltas;
  std::vector < MRI* > mri_warps;
  std::vector < MRI* > mri_weights;
  std::vector < double > intensities;
  MRI * mri_mean;

};


#endif
