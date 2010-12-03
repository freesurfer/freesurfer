/**
 * @file Registration.h
 * @brief A class to compute a robust symmetric registration
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/12/03 02:48:25 $
 *    $Revision: 1.30 $
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
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>

class Registration
{
 template <class T> friend class RegistrationStep;
public:
  Registration(): sat(-1),iscale(false),transonly(false),rigid(true),
      robust(true),rtype(1),subsamplesize(-1),debug(0),verbose(1),initorient(false),
      inittransform(true),highit(-1),mri_source(NULL),mri_target(NULL),
      iscalefinal(1.0),doubleprec(false),wlimit(0.175),symmetry(true),resample(false),
			mri_weights(NULL), mri_hweights(NULL),mri_indexing(NULL)
  {};
  Registration(MRI * s, MRI *t): sat(-1),iscale(false),transonly(false),rigid(true),
      robust(true),rtype(1),subsamplesize(-1),debug(0),verbose(1),initorient(false),
      inittransform(true),highit(-1),mri_source(MRIcopy(s,NULL)),mri_target(MRIcopy(t,NULL)),
      iscalefinal(1.0),doubleprec(false),wlimit(0.175),symmetry(true),resample(false),
			mri_weights(NULL),mri_hweights(NULL),mri_indexing(NULL)
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
  void setHighit(int i)
  {
    highit = i;
  };
  void setDebug(int d)
  {
    debug = d;
    if (d>0) verbose = 2;
  };
  void setVerbose(int i)
	// 0 very quiet, 1 default, 2 detail
  {
    if (i < 0) i=0;
    else if (i > 2) i=2;
    verbose = i;
  };
  void setIscale(bool i)
  {
    iscale = i;
  };
  void setRtype(int r)
  {
    rtype = r;
  };
  void setMinit(const vnl_matrix < double > & m)
  {
    Minit =m; // this is for the original volumes (not in resampled space!)
  };
  void setSource (MRI * s, bool conform = false, bool keeptype = false);
  void setTarget (MRI * t, bool conform = false, bool keeptype = false);
  void setSourceAndTarget(MRI * s, MRI * t, bool keeptype = false);
  void setSubsamplesize (int sss)
  {
    subsamplesize = sss;
  };
  void setName(const std::string &n);
  void setInitOrient(bool io)
  {
    initorient = io;
  };
  void setInitTransform(bool it)
  {
    inittransform = it;
  };
  void setDoublePrec(bool b)
  {
    doubleprec = b;
  };
  void setWLimit( double d)
  {
    wlimit = d;
  };
  void setSymmetry( bool b)
  {
    symmetry = b;
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
  MRI * getHWeights()
  {
    return mri_hweights;
  };
  std::pair <vnl_matrix_fixed< double, 4, 4 > , vnl_matrix_fixed < double, 4, 4 > > getHalfWayMaps();

  double findSaturation (MRI * mriS=NULL, MRI* mriT=NULL, const vnl_matrix < double >& Minit = vnl_matrix < double >(), double iscaleinit = 1.0 );
  std::pair <MATRIX*, double> computeIterativeRegSat( int n,double epsit,MRI * mriS=NULL, MRI* mriT=NULL, MATRIX* Minit = NULL, double iscaleinit = 1.0);

  // compute registration
  virtual void computeIterativeRegistration( int n,double epsit,MRI * mriS=NULL, MRI* mriT=NULL, const vnl_matrix < double > &Minit = vnl_matrix<double>(), double iscaleinit = 1.0);
  void computeMultiresRegistration (int stopres, int n,double epsit, MRI * mriS= NULL, MRI* mriT= NULL, const vnl_matrix < double > &Minit = vnl_matrix<double>(), double iscaleinit = 1.0);

  // get final transform (might be different from mfinal due to possible resampling)
  vnl_matrix_fixed < double , 4 , 4 >  getFinalVox2Vox ();
  double getFinalIscale() {return iscalefinal;};

  double estimateIScale(MRI *mriS, MRI *mriT);

//   bool warpSource(const std::string & fname, MATRIX* M = NULL, double is = -1);
//   bool warpSource(MRI* orig, MRI* target, const std::string &fname, MATRIX* M = NULL, double is = -1);

  // testing
  void testRobust(const std::string & fname, int testno);

  double computeSatEstimate (int reslevel, int n,double epsit, MRI * mriS=NULL, MRI* mriT=NULL, MATRIX* mi=NULL, double scaleinit=1.0 );

  // return centroid (of input, not of resampled internal version)
  std::vector < double > getCentroidS();
  std::vector < double > getCentroidT();
	

protected:

  //   returns weights:
 // std::pair < vnl_vector <double >, MRI* > computeRegistrationStepW(MRI * mriS = NULL, MRI* mriT=NULL);
  //   returns param vector:
//  vnl_vector <double > computeRegistrationStepP(MRI * mriS, MRI* mriT);
  //   returns 4x4 matrix and iscale:
//  std::pair < vnl_matrix_fixed <double,4,4 >, double> computeRegistrationStep(MRI * mriS = NULL, MRI* mriT = NULL);

  //conversion
  std::pair < vnl_matrix_fixed <double,4,4 >, double > convertP2Md(const vnl_vector <double >& p);
  std::pair < MATRIX*, double > convertP2MATRIXd(MATRIX* p);
	
	//MRI * applyParams(MRI * mri_in, const vnl_vector<double>& p, MRI * mri_dst=NULL, bool inverse=false);

  // transform minit into resampled space
	vnl_matrix < double > getMinitResampled();

  // initial registration using moments
  vnl_matrix_fixed < double,4,4>  initializeTransform(MRI *mri_in, MRI *mri_ref);
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
  int verbose;
  bool initorient;
  bool inittransform;
  int highit;

  MRI * mri_source;
  std::vector < MRI* > gpS;
  std::vector < double > centroidS;
  MRI * mri_target;
  std::vector < MRI* > gpT;
  std::vector < double > centroidT;
  vnl_matrix < double >  Minit;
  vnl_matrix < double >  Mfinal;
  double iscalefinal;
  bool doubleprec;
  double wlimit;
  bool symmetry;

  bool resample;
  vnl_matrix < double >  Rsrc;
  vnl_matrix < double >  Rtrg;

private:

  // construct Ab and R:
  //MATRIX* constructR(MATRIX* p);
  //std::pair < MATRIX*, VECTOR* > constructAb(MRI *mriS, MRI *mriT);
  //std::pair < vnl_matrix <double>, vnl_vector <double> > constructAb(MRI *mriS, MRI *mriT);
//  void constructAb(MRI *mriS, MRI *mriT, vnl_matrix < double > &A, vnl_vector < double > &b);
  std::pair < MATRIX*, VECTOR* > constructAb2(MRI *mriS, MRI *mriT);

  // conversions
  MATRIX * rt2mat(MATRIX * r, MATRIX * t, MATRIX *outM); // uses global rtype flag
  MATRIX * p2mat(MATRIX * p6, MATRIX *outM); // calls rt2mat (uses global rtype)

  bool needReslice(MRI *mri, double vsize = -1, int xdim =-1, int ydim=-1, int zdim=-1, bool fixtype = true);
  std::pair< MRI* , vnl_matrix_fixed < double, 4, 4> > makeIsotropic(MRI *mri, MRI *out, double vsize = -1, int xdim =-1, int ydim=-1, int zdim=-1, bool fixtype = true);

  void findSatMultiRes(const vnl_matrix < double > &mi, double scaleinit );

  // gaussian pyramid:
  std::vector < MRI* > buildGaussianPyramid (MRI * mri_in, int min = 16);
  void freeGaussianPyramid(std::vector< MRI* >& p);

  MRI * mri_weights;
  MRI * mri_hweights;
  vnl_matrix < double> mov2weights;
  vnl_matrix < double> dst2weights;
  double wcheck; // set from computeRegistrationStepW
  double wchecksqrt; // set from computeRegistrationStepW
//  double zeroweights;// set from computeRegistrationStepW

  // help vars
//	vnl_vector < double > lastp;

  MRI * mri_indexing;
};


#endif
