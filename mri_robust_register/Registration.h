/**
 * @file Registration.h
 * @brief A class to compute a robust symmetric registration
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/08/13 02:51:19 $
 *    $Revision: 1.18 $
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

// written by Martin Reuter
// Nov. 4th ,2008
//

#ifndef Registration_H
#define Registration_H

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
//#include <iostream>

class Registration
{
public:
  Registration(): sat(-1),iscale(false),transonly(false),rigid(true),
      robust(true),rtype(1),subsamplesize(-1),debug(0),initorient(false),
      mri_source(NULL),mri_target(NULL), Minit(NULL),Mfinal(NULL),
      mri_weights(NULL), mov2weights(NULL),dst2weights(NULL),
      lastp(NULL), mri_indexing(NULL)
  {};
  Registration(MRI * s, MRI *t): sat(-1),iscale(false),transonly(false),rigid(true),
      robust(true),rtype(1),subsamplesize(-1),debug(0),initorient(false),
      mri_source(MRIcopy(s,NULL)),mri_target(MRIcopy(t,NULL)),
      Minit(NULL),Mfinal(NULL),mri_weights(NULL),
      mov2weights(NULL),dst2weights(NULL),lastp(NULL),
      mri_indexing(NULL)
  {};

  virtual ~Registration();

  void clear(); // initialize registration (keep source and target and gauss pyramid)
  void freeGPS() {freeGaussianPyramid(gpS);};
  void freeGPT() {freeGaussianPyramid(gpT);};

  // Set parameters:
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
  void setDebug(int d)
  {
    debug = d;
  };
  void setIscale(bool i)
  {
    iscale = i;
  };
  void setRtype(int r)
  {
    rtype = r;
  };
  void setMinit(MATRIX* m)
  {
    Minit = MatrixCopy(m,Minit);
  };
  void setSource (MRI * s, bool fixvoxel = false, bool fixtype = false);
  void setTarget (MRI * t, bool fixvoxel = false, bool fixtype = false);
  void setSubsamplesize (int sss)
  {
    subsamplesize = sss;
  };
  void setName(const std::string &n);
  void setInitOrient(bool io)
  {
    initorient = io;
  };

  bool isIscale()
  {
    return iscale;
  };
  std::string  getName()
  {
    return name;
  };
  MRI * getWeights()
  {
    return mri_weights;
  };
  std::pair <MATRIX*, MATRIX*> getHalfWayMaps()
  {
    std::pair <MATRIX*, MATRIX*> md2w(mov2weights,dst2weights);
    return md2w;
  };

  // compute registration
  virtual std::pair <MATRIX*, double> computeIterativeRegistration( int n,double epsit,MRI * mriS=NULL, MRI* mriT=NULL, MATRIX* Minit = NULL, double iscaleinit = 1.0);
  std::pair <MATRIX*, double> computeIterativeRegSat( int n,double epsit,MRI * mriS=NULL, MRI* mriT=NULL, MATRIX* Minit = NULL, double iscaleinit = 1.0);
  std::pair <MATRIX*, double> computeMultiresRegistration (int stopres, int n,double epsit, MRI * mriS= NULL, MRI* mriT= NULL, MATRIX* Minit = NULL, double iscaleinit = 1.0);

  bool warpSource(const std::string & fname, MATRIX* M = NULL, double is = -1);
  bool warpSource(MRI* orig, MRI* target, const std::string &fname, MATRIX* M = NULL, double is = -1);

  // testing
  void testRobust(const std::string & fname, int testno);

  double computeSatEstimate (int reslevel, int n,double epsit, MRI * mriS=NULL, MRI* mriT=NULL, MATRIX* mi=NULL, double scaleinit=1.0 );

protected:

  //   returns weights:
  std::pair < MATRIX*, MRI* > computeRegistrationStepW(MRI * mriS = NULL, MRI* mriT=NULL);
  //   returns param vector:
  MATRIX* computeRegistrationStepP(MRI * mriS, MRI* mriT);
  //   returns 4x4 matrix and iscale:
  std::pair <MATRIX*, double> computeRegistrationStep(MRI * mriS = NULL, MRI* mriT = NULL);

  //conversion
  std::pair < MATRIX*, double > convertP2Md(MATRIX* p);

  // initial registration using moments
  MATRIX * initializeTransform(MRI *mri_in, MRI *mri_ref);
  int init_scaling(MRI *mri_in, MRI *mri_ref, MATRIX *m_L); // NOT TESTED !!!!!!!

  double sat;
  bool iscale;
  bool transonly;
  bool rigid;
  bool robust;
  int rtype;
  int subsamplesize;
  std::string name;
  std::string nbase;
  int debug;
  bool initorient;
  //bool outweights;
  //std::string weightsname;

  MRI * mri_source;
  std::vector < MRI* > gpS;
  MRI * mri_target;
  std::vector < MRI* > gpT;
  MATRIX * Minit;
  MATRIX * Mfinal;
  double iscalefinal;



private:

  // construct Ab and R:
  MATRIX* constructR(MATRIX* p);
  std::pair < MATRIX*, VECTOR* > constructAb(MRI *mriS, MRI *mriT);
  std::pair < MATRIX*, VECTOR* > constructAb2(MRI *mriS, MRI *mriT);


  // conversions
  MATRIX * rt2mat(MATRIX * r, MATRIX * t, MATRIX *outM); // uses global rtype flag
  MATRIX * p2mat(MATRIX * p6, MATRIX *outM); // calls rt2mat (uses global rtype)

  // gaussian pyramid:
  std::vector < MRI* > buildGaussianPyramid (MRI * mri_in, int n);
  void freeGaussianPyramid(std::vector< MRI* >& p);

  MRI * mri_weights;
  MATRIX * mov2weights;
  MATRIX * dst2weights;

  // help vars
  MATRIX* lastp;
  double zeroweights;

  MRI * mri_indexing;
};


#endif
