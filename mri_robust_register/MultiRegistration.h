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
 *    $Date: 2011/02/08 22:31:43 $
 *    $Revision: 1.13 $
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
#include "transform.h"
#ifdef __cplusplus
}
#endif

class MultiRegistration
{
public:
  MultiRegistration():outdir("./"),transonly(false),rigid(true),robust(true),sat(4.685),satit(false),
	                     debug(0),iscale(false),subsamplesize(-1),highit(-1),fixvoxel(false),
											 keeptype(false),average(1),doubleprec(false),mri_mean(NULL)
		{};
  MultiRegistration(const std::vector < std::string > mov):outdir("./"),transonly(false),
	                     rigid(true),robust(true),sat(4.685),satit(false),debug(0),iscale(false),
											 subsamplesize(-1),highit(-1),fixvoxel(false),keeptype(false),average(1),doubleprec(false),
											 mri_mean(NULL)
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
	
  bool averageSet(int itdebug = 0, int interp = SAMPLE_TRILINEAR);
  MRI * averageConformSet(int itdebug = 0);
	
  static MRI* averageSet(const std::vector < MRI * >& set,
                       MRI* mean, int method, double sat);
											 
  static MRI* initialAverageSet(const std::vector < MRI * >& set,
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
	
  // DATA
  std::vector < MRI* > mri_mov;
  std::vector < LTA* > ltas;
  std::vector < MRI* > mri_warps;
  std::vector < MRI* > mri_weights;
  std::vector < double > intensities;
  MRI * mri_mean;

};


#endif
