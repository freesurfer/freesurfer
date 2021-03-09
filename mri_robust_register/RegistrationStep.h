/**
 * @brief A class to compute a robust symmetric registration (single step)
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

// written by Martin Reuter
// Jul. 16 ,2010
//
#ifndef RegistrationStep_H
#define RegistrationStep_H

#include "mri.h"

#include <utility>
#include <vector>
#include <cassert>
#include <iostream>
#include <string>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include "MyMRI.h"
#include "MyMatrix.h"
#include "Regression.h"
#include "Quaternion.h"
#include "Transformation.h"
#include "RegRobust.h"

template<class T>
class RegistrationStep
{
public:

  //! Constructor sets member from registration
  RegistrationStep(const RegRobust & R) :
      sat(R.sat), iscale(R.iscale), transonly(R.transonly), rigid(R.rigid), isoscale(
          R.isoscale), trans(R.trans), costfun(R.costfun), rtype(1), subsamplesize(
          R.subsamplesize), debug(R.debug), verbose(R.verbose), floatsvd(false), iscalefinal(
          R.iscalefinal), mri_weights(NULL), mri_indexing(NULL)
  {
  }

  //! Destructor to cleanup index image and weights
  ~RegistrationStep()
  {
    if (mri_indexing)
      MRIfree(&mri_indexing);
    if (mri_weights)
      MRIfree(&mri_weights);
  }

  //! Compute a single registration step
  std::pair<vnl_matrix_fixed<double, 4, 4>, double> computeRegistrationStep(
      MRI * mriS, MRI* mriT);

  double getwcheck()
  {
    return wcheck;
  }

  double getwchecksqrt()
  {
    return wchecksqrt;
  }

  double getzeroweights()
  {
    return zeroweights;
  }

  //! Return pointer to the weights image
  MRI * getWeights()
  {
    return mri_weights;
  }
  //??? who takes care of freeing? Currently we do.

  //! Return 4x4 transformation matrix and intensity scale parameter
  std::pair<vnl_matrix_fixed<double, 4, 4>, double> getMd()
  {
    return Md;
  }

  void setFloatSVD(bool fsvd)
  {
    floatsvd = fsvd;
  }
  // only makes sense for T=double;

  // only public because of resampling testing in Registration.cpp
  // should be made protected at some point.
  void constructAb(MRI *mriS, MRI *mriT, vnl_matrix<T> &A, vnl_vector<T> &b);

  // called from computeRegistrationStepW
  // and externally from RegPowell (not anymore, now use transformation model)
  //static std::pair < vnl_matrix_fixed <double,4,4 >, double > convertP2Md(const vnl_vector < T >& p,bool iscale,int rtype);
  //static std::pair < vnl_matrix_fixed <double,4,4 >, double > convertP2Md2(const vnl_vector < T >& p,bool iscale,int rtype);

protected:

  vnl_matrix<T> constructR(const vnl_vector<T> & p);

private:
// in:

  double sat;
  bool iscale;
  bool transonly;
  bool rigid;
  bool isoscale;
  Transformation* trans;
  Registration::Cost costfun;
  int rtype;
  int subsamplesize;
  int debug;
  int verbose;
  bool floatsvd; // should be removed
  double iscalefinal; // from the last step, used in constructAB

// out:

  std::pair<vnl_matrix_fixed<double, 4, 4>, double> Md;
  MRI* mri_weights;

  double wcheck; // set from computeRegistrationStepW
  double wchecksqrt; // set from computeRegistrationStepW
  double zeroweights; // set from computeRegistrationStepW

//internal
  MRI * mri_indexing;
  vnl_vector<T> pvec;

};

/** Computes Registration Single Step
 The mri's have to be in same space.
 Returns transformation matrix (and created internal float MRI with the weights, if robust, else weights ==NULL).
 Member parameter rtype only for rigid (2: affine restriction to rigid, 1: use rigid from robust-paper)
 */
template<class T>
std::pair<vnl_matrix_fixed<double, 4, 4>, double> RegistrationStep<T>::computeRegistrationStep(
    MRI * mriS, MRI* mriT)
{

  if (trans == NULL)
  {
    cerr<< "ERROR in RegistrationStep: no transform specified ..." << endl;
    exit(1);
  }

  vnl_matrix<T> A;
  vnl_vector<T> b;

  if (rigid && rtype == 2)
  {
    if (verbose > 1)
      std::cout << "rigid and rtype 2 !" << std::endl;
    assert(rtype !=2);

    // compute non rigid A
    rigid = false;
    constructAb(mriS, mriT, A, b);
    rigid = true;
    // now restrict A  (= A R(lastp) )
    vnl_matrix<T> R;
    if (pvec.size() > 0)
      R = constructR(pvec); // construct from last param estimate
    else // construct from identity:
    {
      int l = 6;
      if (!rigid)
        l = 12;
      if (iscale)
        l++;
      vnl_vector<T> tempp(l, 0.0);
      R = constructR(tempp);
      //MatrixPrintFmt(stdout,"% 2.8f",R);exit(1);
    }
    A = A * R.transpose();
  }
  else
  {
    //std::cout << "Rtype  " << rtype << std::endl;

    constructAb(mriS, mriT, A, b);
  }

  if (verbose > 1)
    std::cout << "   - checking A and b for nan ..." << std::flush;
  if (!A.is_finite() || !b.is_finite())
  {
    std::cerr << " A or b constain NAN or infinity values!!" << std::endl;
    exit(1);
  }

  if (verbose > 1)
    std::cout << "  DONE" << std::endl;

  Regression<T> R(A, b);
  R.setVerbose(verbose);
  R.setFloatSvd(floatsvd);
  if (costfun == Registration::ROB)
  {
    vnl_vector<T> w;
    if (verbose > 1)
      std::cout << "   - compute robust estimate ( sat " << sat << " )..."
          << std::flush;
    if (sat < 0)
      pvec = R.getRobustEstW(w);
    else
      pvec = R.getRobustEstW(w, sat);

    A.clear();
    b.clear();

    if (verbose > 1)
      std::cout << "  DONE" << std::endl;

//    std::cout << " pvec  : "<< std::endl;
//    std::cout.precision(16);
//    for (unsigned int iii=0;iii<pvec.size(); iii++)
//      std::cout << pvec[iii] << std::endl;;
//    std::cout.precision(8);

    // transform weights vector back to 3d (mri real)
    if (mri_weights
        && (mri_weights->width != mriS->width
            || mri_weights->height != mriS->height
            || mri_weights->depth != mriS->depth))
      MRIfree(&mri_weights);

    if (!mri_weights)
    {
      mri_weights = MRIallocSequence(mriS->width, mriS->height, mriS->depth, MRI_FLOAT, mriS->nframes);
      MRIcopyHeader(mriS, mri_weights);
      mri_weights->type = MRI_FLOAT;
      MRIsetResolution(mri_weights, mriS->xsize, mriS->ysize, mriS->zsize);
      mri_weights->outside_val = 1.0;
    }

    int x, y, z, f;
    unsigned int count = 0;
    long int val;
    wcheck = 0.0;
    wchecksqrt = 0.0;
    // sigma = max(widht,height,depth) / 6;
    double sigma = mriS->width;
    if (mriS->height > sigma)
      sigma = mriS->height;
    if (mriS->depth > sigma)
      sigma = mriS->depth;
    sigma = sigma / 6.0;
    double sigma22 = 2.0 * sigma * sigma;
    double factor = 1.0 / sqrt(M_PI * sigma22);
    double dsum = 0.0;
    //MRI * gmri = MRIalloc(mriS->width,mriS->height,mriS->depth, MRI_FLOAT);
    for (f = 0; f < mriS->nframes; f++)
      for (z = 0; z < mriS->depth; z++)
        for (x = 0; x < mriS->width; x++)
          for (y = 0; y < mriS->height; y++)
          {
            val = MRILseq_vox(mri_indexing,x,y,z,f) ;
             //if (val < 0) std::cout << " val: " << val << endl;
              if (val < 0) MRIFseq_vox(mri_weights, x, y, z, f) = 1.0;  // anything special set to ignore
              if (val == -4) MRIFseq_vox(mri_weights, x, y, z, f) = 0.0;// background in only one image (0.0 = label outlier)

              //if (val == -10) MRIFseq_vox(mri_weights, x, y, z, f) = -0.5;    // init value (border)
              //else if (val == -1) MRIFseq_vox(mri_weights, x, y, z, f) = -0.6;// zero element (skipped)
              //else if (val == -2) MRIFseq_vox(mri_weights, x, y, z, f) = -0.7;// nan element (skipped)
              //else if (val == -4) MRIFseq_vox(mri_weights, x, y, z, f) = -0.9;// background in only one image
              //else if (val == -5) MRIFseq_vox(mri_weights, x, y, z, f) = -1;  // backgroun/outside in both images

              if (val >= 0.0)
              {
                //std::cout << "val: " << val << "  xyz: " << x << " " << y << " " << z << " " << std::flush;
                assert(val < (int)w.size());
                assert(val >=0);
                double wtemp = w[val] * w[val];
                MRIFseq_vox(mri_weights, x, y, z, f) = wtemp;
                // compute distance to center:
                double xx = x-0.5*mriS->width;
                double yy = y-0.5*mriS->height;
                double zz = z-0.5*mriS->depth;
                double distance2 = xx*xx+yy*yy+zz*zz;
                double gauss = factor * exp(- distance2 / sigma22 );
                dsum += gauss;
                wcheck += gauss * (1.0 - wtemp);
                wchecksqrt += gauss * (1.0 - w[val]);//!!!!! historical, better not use the square root (use wcheck)
                //std::cout << " w^2 : "<< wtemp << "  gauss: " << gauss << "  wcheck+= " << gauss * (1.0 - wtemp) << endl;
                count++;
              }
            }
          //cout << std::endl;
//    MRIwrite(gmri,"mri_gauss.mgz");
//    MRIwrite(mri_indexing, "mri_indexing.mgz");
//    MRIwrite(mri_weights, "mri_weights.mgz");
//    MRIfree(&gmri);
    assert(count == w.size());

    wcheck = wcheck / dsum;
    wchecksqrt = wchecksqrt / dsum;
    if (verbose > 1)
      std::cout << "   - Weight Check: " << wcheck << "  wsqrt: " << wchecksqrt
          << std::endl;
//    if (wcheck > 0.5) 
//    {
//       std::cerr << " Too many voxels in the center are removed! Try to set a larger SAT value! " << std::endl;
//       exit(1);
//    }

  }
  else
  {
    if (verbose > 1)
      std::cout << "   - compute least squares estimate ..." << std::flush;
    pvec = R.getLSEst();

    A.clear();
    b.clear();
    if (verbose > 1)
      std::cout << "  DONE" << std::endl;
    // no weights in this case
    if (mri_weights)
      MRIfree(&mri_weights);
  }

//  zeroweights = R.getLastZeroWeightPercent();
  zeroweights = R.getLastWeightPercent(); // does not need pointers A and B to be valid

//  R.plotPartialSat(name);

//   if (mriS->depth ==1 || mriT->depth ==1)
//   {
//     if (mriS->depth != mriT->depth) { cout << " ERROR both src and trg need to be 2d or 3d!" << endl; exit(1);}
//     Md = convertP2Md2(pvec,iscale,rtype);
//   }
//   else
//     Md = convertP2Md(pvec,iscale,rtype);

  Md.second = 0.0;
  if (iscale)
    Md.second = pvec[pvec.size() - 1];
  trans->setParameters(pvec);
  Md.first = trans->getMatrix();

  return Md;
}

/** Constructs matrix A and vector b for robust regression
   (see Reuter et. al, Neuroimage 2010)
 */
template<class T>
void RegistrationStep<T>::constructAb(MRI *mriS, MRI *mriT, vnl_matrix<T>& A,
    vnl_vector<T>&b)
{

  if (verbose > 1)
    std::cout << "   - constructAb: " << std::endl;

  if (mriS->nframes == 0) mriS->nframes = 1;
  if (mriT->nframes == 0) mriT->nframes = 1;

  assert(mriT != NULL);
  assert(mriS != NULL);
  assert(mriS->width == mriT->width);
  assert(mriS->height== mriT->height);
  assert(mriS->depth == mriT->depth);
  assert(mriS->nframes == mriT->nframes);
  assert(mriS->type == mriT->type);

  bool is2d = false;
  //cout << "Sd: " << mriS->depth << " Td: " << mriT->depth << endl;
  if (mriS->depth == 1 || mriT->depth == 1)
  {
    if (mriT->depth != mriS->depth)
    {
      cout << "ERROR: both source and target need to be 2D or 3D" << endl;
      exit(1);
    }
    is2d = true;
  }

  // Allocate and initialize indexing volume
  int z, y, x, f;
  long int ss = mriS->width * mriS->height * mriS->depth * mriS->nframes;
  if (mri_indexing)
    MRIfree(&mri_indexing);
  int itype = MRI_INT;
  if (ss > std::numeric_limits<int>::max())
  {
    if (verbose > 1)
    {
      std::cout << "     -- using LONG for indexing ... " << std::endl;
      double mu = ((double) ss) * sizeof(long int) / (1024.0 * 1024.0);
      std::cout << "     -- allocating " << mu << "Mb mem for indexing ... "
          << std::flush;
    }
    itype = MRI_LONG;
  }
  else
  {
    if (verbose > 1)
    {
      std::cout << "     -- using INT for indexing ... " << std::endl;
      double mu = ((double) ss) * sizeof(int) / (1024.0 * 1024.0);
      std::cout << "     -- allocating " << mu << "Mb mem for indexing ... "
          << std::flush;
    }
    itype = MRI_INT;
  }
  mri_indexing = MRIallocSequence(mriS->width, mriS->height, mriS->depth, itype, mriS->nframes);
  mri_indexing->outside_val = -10;
  if (mri_indexing == NULL)
    ErrorExit(ERROR_NO_MEMORY,
        "Registration::constructAB could not allocate memory for mri_indexing");
  if (verbose > 1)
    std::cout << " done!" << std::endl;
  // initialize with -10
  for (z = 0; z < mriS->depth; z++)
    for (x = 0; x < mriS->width; x++)
      for (y = 0; y < mriS->height; y++)
        for (f = 0; f < mriS->nframes; f++)
          if (itype == MRI_LONG)
            MRILseq_vox(mri_indexing,x,y,z,f) = -10;
          else
            MRIIseq_vox(mri_indexing,x,y,z,f) = -10;

  // determine if we will subsample below:
  bool dosubsample = false;
  if (subsamplesize > 0)
    dosubsample = (mriS->width > subsamplesize && mriS->height > subsamplesize
        && (mriS->depth > subsamplesize || mriS->depth == 1));

  // we will need the derivatives (fx1,fy1,fz1), smoothed image (ft1) and average (SpTh)
  if (verbose > 1)
    std::cout << "     -- compute derivatives ... " << std::flush;
  MRI *SpTh = MRIallocSequence(mriS->width, mriS->height, mriS->depth, MRI_FLOAT, mriS->nframes);
  SpTh = MRIadd(mriS, mriT, SpTh);
  SpTh = MRIscalarMul(SpTh, SpTh, 0.5);
  MRI *fx1 = NULL, *fy1 = NULL, *fz1 = NULL, *ft1 = NULL;
  MyMRI::getPartials(SpTh, fx1, fy1, fz1, ft1);
  MRIfree(&SpTh);
  MRI *SmT = MRIallocSequence(mriS->width, mriS->height, mriS->depth, MRI_FLOAT, mriS->nframes);
  SmT = MRIsubtract(mriS, mriT, SmT);
  SmT = MyMRI::getBlur(SmT, SmT);

  if (verbose > 1)
    std::cout << " done!" << std::endl;
  //MRIwrite(fx1,"fx.mgz");
  //MRIwrite(fy1,"fy.mgz");
  //MRIwrite(fz1,"fz.mgz");
  //MRIwrite(ft1,"ft.mgz");

  // subsample if wanted
  MRI *fx, *fy, *fz, *ft;
  if (dosubsample)
  {
    if (verbose > 1)
      std::cout << "     -- subsample ... " << std::flush;

    // by default the subsample routine uses random offsets, we compute the coordinates below again
    // when looping through the subsampled image (for index creation and for constructing the matrices)
    fx = MyMRI::subSample(fx1, NULL, false, 0);
    MRIfree(&fx1);
    fy = MyMRI::subSample(fy1, NULL, false, 0);
    MRIfree(&fy1);
    if (fz1)
    {
      fz = MyMRI::subSample(fz1, NULL, false, 0);
      MRIfree(&fz1);
    }
    else
      fz = NULL;
    ft = MyMRI::subSample(ft1, NULL, false, 0);
    MRIfree(&ft1);

    MRI *SmTt = SmT;
    SmT = MyMRI::subSample(SmTt, NULL, false, 0);
    MRIfree(&SmTt);

    if (verbose > 1)
      std::cout << " done! " << std::endl;
  }
  else //just rename
  {
    fx = fx1;
    fy = fy1;
    fz = fz1;
    ft = ft1;
  }

//cout << " size fx : " << fx->width << " , " << fx->height << " , " << fx->depth << std::endl;
//cout << " size src: " << mriS->width << " , " << mriS->height << " , " << mriS->depth << std::endl;

  // compute 'counti': the number of rows needed (zero elements need to be removed)
  int n = fx->width * fx->height * fx->depth * fx->nframes;
  if (verbose > 1)
    std::cout << "     -- size " << fx->width << " x " << fx->height << " x "
        << fx->depth << " x " << fx->nframes << " = " << n << std::flush;
  long int counti = 0;
  double eps = 0.00001;
  double oepss = eps+mriS->outside_val/255.0;
  double oepst = eps+mriT->outside_val/255.0;
  int fxd = fx->depth;
  int fxw = fx->width;
  int fxh = fx->height;
  int fxf = fx->nframes;
  int fxstart = 0;
  int xp1, yp1, zp1;
  int ocount = 0, ncount = 0, zcount = 0;
  float fzval = eps/2.0;
  int dx, dy, dz;
  int randpos = 0;
  for (z = fxstart; z < fxd; z++)
    for (x = fxstart; x < fxw; x++)
      for (y = fxstart; y < fxh; y++)
      {
        // check if position is outside either source or target:
        if (dosubsample)
        {
          // dx,dy and dz need to agree with the subsampling above
          dx = (int)(2.0*MyMRI::getRand(randpos));
          randpos++;
          dy = (int)(2.0*MyMRI::getRand(randpos));
          randpos++;
          xp1 = 2*x+dx;
          yp1 = 2*y+dy;
          if (is2d) zp1 = z;
          else
          {
            dz = (int)(2.0*MyMRI::getRand(randpos));
            randpos++;
            zp1 = 2*z+dz;
          }
        }
        else 
        {
          xp1 = x;
          yp1 = y;
          zp1 = z; 
        }
        assert(xp1 < mriS->width);
        assert(yp1 < mriS->height);
        assert(zp1 < mriS->depth);
//        if ( MRIgetVoxVal(mriS,xp1,yp1,zp1,0) == mriS->outside_val || MRIgetVoxVal(mriT,xp1,yp1,zp1,0) == mriT->outside_val )
        if ( fabs(MRIgetVoxVal(mriS,xp1,yp1,zp1,0)- mriS->outside_val) <= oepss || fabs(MRIgetVoxVal(mriT,xp1,yp1,zp1,0)- mriT->outside_val)<=oepst )
        {
          //std::cout << "voxel outside (" << xp1 << " " << yp1 << " " << zp1 << " )  mriS: " <<MRIFvox(mriS,xp1,yp1,zp1) << "  mriT: " << MRIFvox(mriT,xp1,yp1,zp1)  << "  ovalS: " << mriS->outside_val << "  ovalT: " << mriT->outside_val<< std::endl;
          ocount+=fxf; // will be outside in all frames then
          continue;
        }

        // nan and zero values will also be skipped 
        for (f=0;f<fxf;f++)
        {
          const float & ftval = MRIFseq_vox(ft, x, y, z, f);
          const float & fxval = MRIFseq_vox(fx, x, y, z, f);
          const float & fyval = MRIFseq_vox(fy, x, y, z, f);
          if (!is2d)
            fzval = MRIFseq_vox(fz, x, y, z, f) ;
                
          if (isnan(fxval) || isnan(fyval) || isnan(fzval) || isnan(ftval) )
          {
            //if (verbose > 0) std::cout << " found a nan value!!!" << std::endl;
            ncount++;
            continue;
          }
          if (fabs(fxval) < eps  && fabs(fyval) < eps && fabs(fzval) < eps )
          {
            //if (verbose > 0) std::cout << " found a zero element !!!" << std::endl;
            zcount++;
            continue;
          }
          counti++; // found another good voxel
         }
       }
      
  if (verbose > 1 && n > counti)
    std::cout << "  need only: " << counti << std::endl;

  if (counti == 0)
  {
    std::cerr << std::endl;
    std::cerr
        << " ERROR: All entries are zero! Images do not overlap (anymore?)."
        << std::endl;
    std::cerr
        << "    This can have several reasons (i.e. different modalities, different "
        << std::endl;
    std::cerr
        << "    intensity scales, large non-linearities, too diff. voxel sizes ...)"
        << std::endl;
    std::cerr
        << "    Try calling with --noinit (if the original images are well aligned)"
        << std::endl;
    std::cerr << "    Maybe use --ixform <init.lta> with an approx. alignment"
        << std::endl;
    std::cerr << "    obtained from tkregister or another registration program."
        << std::endl;
    std::cerr << "    Or do some prior intensity correction? " << std::endl;
    std::cerr
        << "    You can also try to switch off symmetric reg. via --nosym "
        << std::endl;
    std::cerr << "    or manually play around with the --sat parameter. "
        << std::endl;
    std::cerr << std::endl;
    exit(1);
  }

  if (verbose > 1)
    cout << "     -- nans: " << ncount << " zeros: " << zcount << " outside: "
        << ocount << endl;

  // allocate the space for A and B
  int pnum = trans->getDOF();
  if (iscale)
    pnum++;
  //cout << " pnum: " << pnum << "  counti: " << counti<<  endl;
  double amu = ((double) counti * (pnum + 1)) * sizeof(T) / (1024.0 * 1024.0); // +1 =  rowpointer vector
  double bmu = (double) counti * sizeof(T) / (1024.0 * 1024.0);
  if (verbose > 1)
    std::cout << "     -- allocating " << amu + bmu << "Mb mem for A and b ... "
        << std::flush;
  bool OK = A.set_size(counti, pnum);
  OK = OK && b.set_size(counti);
  if (!OK)
  {
    std::cout << std::endl;
    ErrorExit(ERROR_NO_MEMORY,
        "Registration::constructAB could not allocate memory for A and b");
  }
  if (verbose > 1)
    std::cout << " done! " << std::endl;
  double maxmu = 5 * amu + 7 * bmu;
  string fstr = "";
  if (floatsvd)
  {
    maxmu = amu + 3 * bmu + 2 * (amu + bmu);
    fstr = "-float";
  }
  if (verbose > 1)
    std::cout << "         (MAX usage in SVD" << fstr << " will be > " << maxmu
        << "Mb mem + 6 MRI) " << std::endl;
  if (maxmu > 3800)
  {
    std::cout << "     -- WARNING: mem usage large: " << maxmu
        << "Mb mem + 6 MRI" << std::endl;
    //string fsvd;
    //if (doubleprec) fsvd = "remove --doubleprec and/or ";
    std::cout << "          Maybe use --subsample <int> " << std::endl;
  }

//        char ch;
//        std::cout << "Press a key to continue iterations: ";
//        std::cin  >> ch;

  // Loop and construct A and b
  long int count = 0;
  ocount = 0;
  randpos = 0;
  fzval = eps/2.0;

  for (z = fxstart; z < fxd; z++)
    for (x = fxstart; x < fxw; x++)
      for (y = fxstart; y < fxh; y++)
      {
        // check if position is outside either source or target:
        if (dosubsample)
        {
          // dx,dy and dz need to agree with the subsampling above
          dx = (int)(2.0*MyMRI::getRand(randpos));
          randpos++;
          dy = (int)(2.0*MyMRI::getRand(randpos));
          randpos++;
          xp1 = 2*x+dx;
          yp1 = 2*y+dy;
          if (is2d) zp1 = z;
          else
          {
            dz = (int)(2.0*MyMRI::getRand(randpos));
            randpos++;
            zp1 = 2*z+dz;
          }
        }
        else // if not subsampled
        {
          xp1 = x;
          yp1 = y;
          zp1 = z;
        }
        assert(xp1 < mriS->width);
        assert(yp1 < mriS->height);
        assert(zp1 < mriS->depth);
        const float & mriSval = MRIgetVoxVal(mriS,xp1,yp1,zp1,0);
        const float & mriTval = MRIgetVoxVal(mriT,xp1,yp1,zp1,0);
        
//        if ( mriSval == mriS->outside_val || mriTval == mriT->outside_val )
        if ( fabs(mriSval- mriS->outside_val) <= oepss || fabs(mriTval- mriT->outside_val)<=oepst )
        {
          //std::cout << "voxel outside (" << xp1 << " " << yp1 << " " << zp1 << " )  mriS: " <<MRIFvox(mriS,xp1,yp1,zp1) << "  mriT: " << MRIFvox(mriT,xp1,yp1,zp1)  << "  ovalS: " << mriS->outside_val << "  ovalT: " << mriT->outside_val<< std::endl;
          int outval = -4;
          if (fabs(mriSval - mriS->outside_val)<=oepss && fabs(mriTval- mriT->outside_val) <= oepst )
            outval = -5;
          for (f=0;f<fxf;f++)
            MRILseq_vox(mri_indexing, xp1, yp1, zp1,f) = outval;
          ocount+=fxf;
          //cout << " " << ocount << flush;
          continue;
        }
        
        // loop through all frames
        for (f=0; f<fxf; f++)
        {
          const float & ftval = MRIFseq_vox(ft, x, y, z, f);
          const float & fxval = MRIFseq_vox(fx, x, y, z, f);
          const float & fyval = MRIFseq_vox(fy, x, y, z, f);
          if (!is2d) fzval = MRIFseq_vox(fz, x, y, z, f);
          
          // skip nans or zeros
          if (isnan(fxval) || isnan(fyval) || isnan(fzval) || isnan(ftval) )
          {
            //if (verbose > 0) std::cout << " found a nan value!!!" << std::endl;
            MRILseq_vox(mri_indexing, xp1, yp1, zp1, f) = -2;
            continue;
          }
          if (fabs(fxval) < eps && fabs(fyval) < eps && fabs(fzval) < eps )
          {
            //if (verbose > 0) std::cout << " found a zero element !!!" << std::endl;
            MRILseq_vox(mri_indexing, xp1, yp1, zp1, f) = -1;
            continue;
          }

          assert(counti > count);

          MRILseq_vox(mri_indexing, xp1, yp1, zp1, f) = count;

          //cout << "x: " << x << " y: " << y << " z: " << z << " count: "<< count << std::endl;
          //cout << " " << count << " mrifx: " << MRIFvox(mri_fx, x, y, z) << " mrifx int: " << (int)MRIvox(mri_fx,x,y,z) <<endl;

          // new: now use transformation model to get the gradient vector
          vnl_vector < double > grad = trans->getGradient(x,fxval,y,fyval,z,fzval);
          int dof = grad.size();
          for (int pno = 0; pno < dof; pno++)
          {
            A[count][pno] = grad[pno];
          }

//         if (transonly)
//         {
//           if (is2d)
//           {
//             vnl_vector < double > grad = Trans->getGradient(x,fxval,y,fyval,z,fzval);
//             dof = grad.size();
//             for (int pno = 0; pno < dof; pno++)
//             {
//               A[count][pno] =  grad[pno];
//             }
//            // A[count][0] = fxval;
//            // A[count][1] = fyval;
//            // dof = 2;
//           }
//           else
//           {
//             A[count][0] = fxval;
//             A[count][1] = fyval;
//             A[count][2] = fzval;
//             dof = 3;
//           }
//         }
//         else if (rigid)
//         {
//           if (is2d)
//           {
//             vnl_vector < double > grad = Trans->getGradient(x,fxval,y,fyval,z,fzval);
//             dof = grad.size();
//             for (int pno = 0; pno < dof; pno++)
//             {
//               A[count][pno] =  grad[pno];
//             }
//             //A[count][0] =  fxval;
//             //A[count][1] =  fyval;
//             //A[count][2] = (fyval*xp1 - fxval*yp1);
//             //dof = 3;     
//           }
//           else
//           {
//             A[count][0] =  fxval;
//             A[count][1] =  fyval;
//             A[count][2] =  fzval;
//             A[count][3] = (fzval*yp1 - fyval*zp1);
//             A[count][4] = (fxval*zp1 - fzval*xp1);
//             A[count][5] = (fyval*xp1 - fxval*yp1);
//             dof = 6;
//           }
//           
//         }
//         else if (isoscale)
//         {
//           if (is2d) // [ p -q ; q p ] + T
//           {
//             vnl_vector < double > grad = Trans->getGradient(x,fxval,y,fyval,z,fzval);
//             dof = grad.size();
//             for (int pno = 0; pno < dof; pno++)
//             {
//               A[count][pno] =  grad[pno];
//             }
//             //A[count][0] =  fxval;
//             //A[count][1] =  fyval;
//             //A[count][2] =  fxval*xp1 + fyval*yp1;
//             //A[count][3] = -fxval*yp1 + fyval*xp1;
//             //dof = 4;
//           }
//           else
//           {
//             A[count][0] =  fxval;
//             A[count][1] =  fyval;
//             A[count][2] =  fzval;
//             A[count][3] = (fzval*yp1 - fyval*zp1);
//             A[count][4] = (fxval*zp1 - fzval*xp1);
//             A[count][5] = (fyval*xp1 - fxval*yp1);
//             A[count][6] = (fxval*xp1 + fyval*yp1);
//             dof = 7;
//             cerr << " Isoscale in 3D not implemented yet, use ridig or affine" <<endl;
//             exit(1);
//           }        
//         }
//         else // affine
//         {
//           if (is2d)
//           {
//             //A[count][0]  = fxval*xp1;
//             //A[count][1]  = fxval*yp1;
//             //A[count][2]  = fxval;
//             //A[count][3]  = fyval*xp1;
//             //A[count][4]  = fyval*yp1;
//             //A[count][5]  = fyval;
//             //dof = 6;
//             
//             vnl_vector < double > grad = Trans->getGradient(x,fxval,y,fyval,z,fzval);
//             dof = grad.size();
//             for (int pno = 0; pno < dof; pno++)
//             {
//               A[count][pno] =  grad[pno];
//             }
//           }
//           else
//           {
//             A[count][0]  = fxval*xp1;
//             A[count][1]  = fxval*yp1;
//             A[count][2]  = fxval*zp1;
//             A[count][3]  = fxval;
//             A[count][4]  = fyval*xp1;
//             A[count][5]  = fyval*yp1;
//             A[count][6]  = fyval*zp1;
//             A[count][7]  = fyval;
//             A[count][8]  = fzval*xp1;
//             A[count][9]  = fzval*yp1;
//             A[count][10] = fzval*zp1;
//             A[count][11] = fzval;
//             dof = 12;
//           }
//         }

          // ISCALE
          // intensity model: R(s,IS,IT) = exp(-0.5 s) IT - exp(0.5 s) IS
          //                  R'  = -0.5 ( exp(-0.5 s) IT + exp(0.5 s) IS)
          //   ft = 0.5 ( exp(-0.5s) IT + exp(0.5s) IS)  (average of intensity adjusted images)
          if (iscale) A[count][dof] = ftval;

          // A p = b = IS - IT
          b[count] = MRIFseq_vox(SmT, x, y, z, f);

          count++;// start with 0 above

        }
      }
        //cout << " ocount : " << ocount << endl;    
        //cout << " counti: " << counti << " count : " << count<< endl;    
  assert(counti == count);

//   vnl_matlab_print(vcl_cerr,A,"A",vnl_matlab_print_format_long);std::cerr << std::endl;    
//   vnl_matlab_print(vcl_cerr,b,"b",vnl_matlab_print_format_long);std::cerr << std::endl;    

  // free remaining MRI    
  MRIfree(&fx);
  MRIfree(&fy);
  if (fz)
    MRIfree(&fz);
  MRIfree(&ft);
  MRIfree(&SmT);
//MRIwrite(mri_indexing,"mriindexing2.mgz");
//exit(1);
  return;
}

// template <class T>
// pair < vnl_matrix_fixed <double,4,4 >, double > RegistrationStep<T>::convertP2Md(const vnl_vector < T >& p, bool iscale, int rtype)
// // rtype : use restriction (if 2) or rigid from robust paper
// // returns registration as 4x4 matrix M, and iscale
// {
// //   std::cout << " RegistrationStep<T>::convertP2Md(MATRIX* p) (p->rows: " << p->rows << " )" << std::flush;
//   std::pair < vnl_matrix_fixed <double,4,4 >, double> ret; ret.second = 0.0;
// 
//   int psize = p.size();
//   
//   if (iscale)
//   {
//     //std::cout << " has intensity " << std::endl;
//     // last is intensity scale    
//     // ISCALECHANGE:
//     psize--;
//     ret.second =  (double) p[psize];
//   }
//   
//   if (rtype == 1)
//   {
//     Transformation * Trans = NULL;
//     if (psize == 12)     Trans = new Transform3dAffine(p);
//     else if (psize == 7) Trans = new Transform3dIsoscale(p);
//     else if (psize == 6) Trans = new Transform3dRigid(p);
//     else if (psize == 3) Trans = new Transform3dTranslate(p);
//     else
//     {
//       std::cerr << " ERROR: unknown 3d transformation (type 1) with " << psize << " DOF" << std::endl;
//       assert(1==2);
//     }
//     
//     ret.first = Trans->getMatrix();
//     free(Trans);
//     return ret;
//   }
//   else if (rtype == 2)
//   {
//     Transformation * Trans = NULL;
//     if (psize == 12)     Trans = new Transform3dAffine2(p);
//     else if (psize == 7) Trans = new Transform3dIsoscale2(p);
//     else if (psize == 6) Trans = new Transform3dRigid2(p);
//     else if (psize == 3) Trans = new Transform3dTranslate(p);
//     else 
//     {
//       std::cerr << " ERROR: unknown 3d transformation (type 2) with " << psize << " DOF" << std::endl;
//       assert(1==2);
//     }
//     
//     ret.first = Trans->getMatrix();
//     free(Trans);
//     return ret;  
//   }
//   else
//   {
//     std::cerr << " ERROR: unknown rtype (should never get here) ..." << std::endl;
//     assert(1==2);
//   }
//   
//   
//   
//   // now transformation parameters:
//   
//   if (psize == 12) // affine
//   {
//     if (rtype == 1)
//     {
//       // affine, just the 12 parameters as matrix add-ons
//       ret.first.set_identity();
//       int count = 0;
//       for (int rr = 0;rr<3;rr++)
//       for (int cc = 0;cc<4;cc++)
//       {
//         ret.first[rr][cc] +=  p[count];
//         count++;
//       }
//     }
//     else if (rtype == 2)
//     {
//       // M = T*shear*Scale*Rot
//       
//       //Rot
//       Quaternion q;
//       q.importZYXAngles(-p[5], p[4], -p[3]); // same as spm now
//       vnl_matrix < double > rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(),3);
//       //Scale
//       vnl_matrix < double > smat(3,3,0.0);
//       smat[0][0] = p[6]; smat[1][1] = p[7]; smat[2][2] = p[8];
//       //Shear
//       vnl_matrix <double > zmat(3,3); zmat.set_identity();
//       zmat[0][1] = p[9]; zmat[0][2] = p[10]; zmat[1][2] = p[11];
//       // product 3x3
//       vnl_matrix <double> M3 = zmat * smat * rmat;
//       // consturct 4x4 with translation also:
//       int rr, cc;
//       for (rr=0;rr<3;rr++)
//       {
//         for (cc=0;cc<3;cc++) // copy M3
//           ret.first[rr][cc] =M3[rr][cc];
// 
//         // copy translation into 4th column
//         ret.first[rr][3] = p[rr];
//         // set 4th row to zero
//         ret.first[3][rr] = 0.0;
//       }
//       //except 4,4
//       ret.first[3][3] = 1.0;
//       
//       
//     }
//     else assert(1==2);
//   } 
//   else if (psize == 7) // rigid and isotropic scale: tx,ty,tz,r1,r2,r3,s
//   {
//     
//     // M = T*(Scale*Rot)
// 
//     //Rotation
//     Quaternion q;
//     if (rtype == 1)
//     {
//       q.importRotVec(p[3],p[4],p[5]);
//     }
//     else if (rtype == 2)
//     {
//       // first convert rotation to quaternion (clockwise)
//       //q.importZYXAngles(-p[5], -p[4], -p[3]);
//       q.importZYXAngles(-p[5],p[4],-p[3]); // same as spm now
//     }           
//     else assert(1==2);    
//     vnl_matrix < double > rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(),3);
// 
//     // scale
//     rmat = ((double)p[6]) * rmat;
//       
//     // copy 
//     ret.first.set_identity();
//     int rr, cc;
//     for (rr=0;rr<3;rr++)
//       for (cc=0;cc<3;cc++) // copy 
//         ret.first[rr][cc] = rmat[rr][cc];
// 
//     // translation
//     ret.first[0][3] = p[0];
//     ret.first[1][3] = p[1];
//     ret.first[2][3] = p[2];
//   
//   }
//   else if (psize == 6) // rigid tx,ty,tz,r1,r2,r3
//   {
//     // converts rot vector (3x1) and translation vector (3x1)
//     // into an affine matrix (homogeneous coord) 4x4
//     // if global rtype ==1 r1,r2,r3 are as in robust paper (axis, and length is angle)
//     // if global rtype ==2 then r1,r2,r3 are angles around x,y,z axis (order 1zrot,2yrot,3xrot)
//     Quaternion q;
//     if (rtype == 2)
//     {
//       // first convert rotation to quaternion (clockwise)
//       //q.importZYXAngles(-p[5], -p[4], -p[3]);
//       q.importZYXAngles(-p[5], p[4], -p[3]); // same as spm now
//     }
//     else if (rtype == 1)
//     {
//       // first convert rotation to quaternion;
//       q.importRotVec(p[3],p[4],p[5]);
//     }
//     else assert (1==2);
//     // then to rotation matrix
//     vnl_matrix < double > rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(),3);
//     
//     int rr, cc;
//     for (rr=0;rr<3;rr++)
//     {
//       for (cc=0;cc<3;cc++) // copy rot-matrix
//         ret.first[rr][cc] = rmat[rr][cc];
// 
//       // copy translation into 4th column
//       ret.first[rr][3] = p[rr];
//       // set 4th row to zero
//       ret.first[3][rr] = 0.0;
//     }
//     //except 4,4
//     ret.first[3][3] = 1.0;
//   }
//   else if (psize == 3) // translation only
//   {
//     ret.first.set_identity();
//     ret.first[0][3] = p[0];
//     ret.first[1][3] = p[1];
//     ret.first[2][3] = p[2];
//   }
//   else
//   {
//     cerr << " transformation neither 3,6,7 nor 12 dof : " << psize <<" ??" << std::endl;
//     assert(1==2);
//   }
// 
// //   std::cout << " -- DONE " << std::endl;
//   return ret;
// }
// 
// template <class T>
// pair < vnl_matrix_fixed <double,4,4 >, double > RegistrationStep<T>::convertP2Md2(const vnl_vector < T >& p, bool iscale, int rtype)
// // rtype : use restriction (if 2) or rigid from robust paper
// // returns registration as 4x4 matrix M, and iscale
// {
// //  if (iscale) std::cout << " RegistrationStep<T>::convertP2Md2(p,iscale=true,"<<rtype<<") (p length: " << p.size() << " )" << std::flush;
// //  else  std::cout << " RegistrationStep<T>::convertP2Md2(p,iscale=false,"<<rtype<<") (p length: " << p.size() << " )" << std::flush;
// //  std::cout<< endl; vnl_matlab_print(vcl_cout,p,"p",vnl_matlab_print_format_long);std::cout << std::endl;
//   std::pair < vnl_matrix_fixed <double,4,4 >, double> ret; ret.second = 0.0;
// 
//   int psize = p.size();
//   
//   if (iscale) // iscale
//   {
//     //std::cout << " has intensity " << std::endl;
//     // last is intensity scale    
//     // ISCALECHANGE:
//     psize--;
//     ret.second =  (double) p[psize];
//   }
// 
//   if (rtype == 1)
//   {
//     Transformation * Trans = NULL;
//     if (psize == 6)      Trans = new Transform2dAffine(p);
//     else if (psize == 4) Trans = new Transform2dIsoscale(p);
//     else if (psize == 3) Trans = new Transform2dRigid(p);
//     else if (psize == 2) Trans = new Transform2dTranslate(p);
//     else
//     {
//       std::cerr << " ERROR: unknown 2d transformation (type 1) with " << psize << " DOF" << std::endl;
//       assert(1==2);
//     }
//     
//     ret.first = Trans->getMatrix();
//     free(Trans);
//     return ret;
//   }
//   else if (rtype == 2)
//   {
//     Transformation * Trans = NULL;
//     if (psize == 6)      Trans = new Transform2dAffine2(p);
//     else if (psize == 4) Trans = new Transform2dIsoscale2(p);
//     else if (psize == 3) Trans = new Transform2dRigid(p);
//     else if (psize == 2) Trans = new Transform2dTranslate(p);
//     else
//     {
//       std::cerr << " ERROR: unknown 2d transformation (type 2) with " << psize << " DOF" << std::endl;
//       assert(1==2);
//     }
//     
//     ret.first = Trans->getMatrix();
//     free(Trans);
//     return ret;
//   }
//   else
//   {
//     cerr << " ERROR: unknown rtype (should never get here) ..." << std::endl;
//     assert(1==2);
//   }
//   
// 
// 
//   
//   // now transformation parameters:
//   
//   if (psize == 6) //AFFINE 2D
//   {
//     if (rtype == 1)
//     {
//       // affine, just the 6 parameters as matrix add-ons
//       /*ret.first.set_identity();
//       ret.first[0][0] += p[0];
//       ret.first[0][1] += p[1];
//       ret.first[0][3] += p[2];
//       ret.first[1][0] += p[3];
//       ret.first[1][1] += p[4];
//       ret.first[1][3] += p[5];*/
//       Transform2dAffine A2d(p);
//       ret.first = A2d.getMatrix();
//     }
//     else if (rtype == 2)
//     {
//       // M = T*shear*Scale*Rot
//       // Translation
//       vnl_vector_fixed <double,3 > t(p[0],p[1],0);   
//       //Rot
//       Quaternion q;
//       q.importZYXAngles(-p[2],0,0);
//       vnl_matrix < double > rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(),3);
//       //Scale
//       vnl_matrix < double > smat(3,3,0.0);
//       smat[0][0] = p[3]; smat[1][1] = p[4]; smat[2][2] = 1;
//       //Shear
//       vnl_matrix <double > zmat(3,3); zmat.set_identity();
//       zmat[0][1] = p[5]; 
//       // product 3x3
//       vnl_matrix <double> M3 = zmat * smat * rmat;
//       // consturct 4x4 with translation also:
//       int rr, cc;
//       for (rr=0;rr<3;rr++)
//       {
//         for (cc=0;cc<3;cc++) // copy M3
//           ret.first[rr][cc] =M3[rr][cc];
// 
//         // copy translation into 4th column
//         ret.first[rr][3] = t[rr];
//         // set 4th row to zero
//         ret.first[3][rr] = 0.0;
//       }
//       //except 4,4
//       ret.first[3][3] = 1.0;
//       
//       
//     }
//     else assert(1==2);
//   } 
//   else if (psize == 4) // rigid and isotropic scaling
//   {
//     if (rtype == 1)
//     {
// //       // matrix add-ons
// //       ret.first.set_identity();
// //       ret.first[0][0] += p[2];
// //       ret.first[0][1] += -p[3];
// //       ret.first[0][3] += p[0];
// //       ret.first[1][0] += p[3];
// //       ret.first[1][1] += p[2];
// //       ret.first[1][3] += p[1];
//       Transform2dIsoscale A2d(p);
//       ret.first = A2d.getMatrix();
//     }
//     else if (rtype == 2)
//     {
//       // M = T*(Scale*Rot)
//       //Rot
//       Quaternion q;
//       q.importZYXAngles(-p[2],0,0);
//       vnl_matrix < double > rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(),3);
// 
//       // scale
//       ret.first.set_identity();
//       ret.first[0][0] = p[3] * rmat[0][0];
//       ret.first[0][1] = p[3] * rmat[0][1];
//       ret.first[1][0] = p[3] * rmat[1][0];
//       ret.first[1][1] = p[3] * rmat[1][1];
//      
//       // translation
//       ret.first[0][3] = p[0];
//       ret.first[1][3] = p[1];
//     }
//     else assert(1==2);
//   
//   
//   }
//   else if (psize == 3) // rigid (rot and trans xy)
//   {
//     if (rtype == 1)
//     {
// //       // matrix add-ons
// //       ret.first.set_identity();
// //       ret.first[0][0] += p[2];
// //       ret.first[0][1] += -p[3];
// //       ret.first[0][3] += p[0];
// //       ret.first[1][0] += p[3];
// //       ret.first[1][1] += p[2];
// //       ret.first[1][3] += p[1];
//       Transform2dRigid A2d(p);
//       ret.first = A2d.getMatrix();
//     }
//     else if (rtype == 2)
//     {
//      // rigid: first 2 translation, next rotation (as a vector)
//     // split translation and rotation:
//     vnl_vector_fixed <double,3 > t(p[0],p[1],0);
//     double r = p[2];
// 
//     // converts rot vector (3x1) and translation vector (3x1)
//     // into an affine matrix (homogeneous coord) 4x4
//     // if global rtype ==1 r1,r2,r3 are as in robust paper (axis, and length is angle)
//     // if global rtype ==2 then r1,r2,r3 are angles around x,y,z axis (order 1zrot,2yrot,3xrot)
//     vnl_matrix < double > rmat;
//     Quaternion q;
//     if (rtype == 2)
//     {
//       // first convert rotation to quaternion (clockwise)
//       //q.importZYXAngles(-r[2], -r[1], -r[0]);
//       q.importZYXAngles(-r, 0, 0); // same as spm now
//     }
//     else if (rtype == 1)
//     {
//       // first convert rotation to quaternion;
//       q.importRotVec(0,0,r);
//       //cout << " r: " << r << endl;
//       //cout << " q: " << q << endl;
//     }
//     else assert (1==2);
//     // then to rotation matrix
//     rmat = MyMatrix::getVNLMatrix(q.getRotMatrix3d(),3);
//     
//     int rr, cc;
//     for (rr=0;rr<3;rr++)
//     {
//       for (cc=0;cc<3;cc++) // copy rot-matrix
//         ret.first[rr][cc] = rmat[rr][cc];
// 
//       // copy translation into 4th column
//       ret.first[rr][3] = t[rr];
//       // set 4th row to zero
//       ret.first[3][rr] = 0.0;
//     }
//     //except 4,4
//     ret.first[3][3] = 1.0;
//     }
//   }
//   else if (psize == 2) // translation only
//   {
//       Transform2dTranslate A2d(p);
//       ret.first = A2d.getMatrix();
// //    ret.first.set_identity();
// //    ret.first[0][3] = p[0];
// //    ret.first[1][3] = p[1];
// //    ret.first[2][3] = 0;
//   }
//   else
//   {
//     cerr << " transformation neither 6,4,3 nor 2 dof : " << psize <<" ??" << std::endl;
//     assert(1==2);
//   }
// //  std::cout<< endl; vnl_matlab_print(vcl_cout,ret.first,"Mt",vnl_matlab_print_format_long);std::cout << std::endl;
// 
// //   std::cout << " -- DONE " << std::endl;
//   return ret;
// }

/** Construct restriction matrix (to restrict the affine problem to less parameters)
 ! Experimental !
 if p->rows == 6 use only rigid
 if p->rows == 7 use also intensity scale
 if p->rows == 3 use only trans
 if p->rows == 4 use only trans + intensity
 */
template<class T>
vnl_matrix<T> RegistrationStep<T>::constructR(const vnl_vector<T> & p)
{
  assert(p.size() > 0);
  assert(p.size() == 6 || p.size()==7);

  int adim = 12;
  if (iscale)
  {
    assert(p.size() == 7 || p.size() ==4);
    adim++;
  }
  vnl_matrix<T> R(p.size(), adim, 0.0);

  // translation p0,p1,p2 map to m3,m7,m11, (counting from zero)
  R[0][3] = 1.0;
  R[1][7] = 1.0;
  R[2][11] = 1.0;

  // iscale (p6 -> m12)
  if (p.size() == 7)
    R[6][12] = 1.0;
  if (p.size() == 4)
    R[3][12] = 1.0;

  if (p.size() <= 4)
    return R;

  // rotation derivatives (dm_i/dp_i)
  double s4 = sin(p[3]);
  double c4 = cos(p[3]);
  double s5 = sin(p[4]);
  double c5 = cos(p[4]);
  double s6 = sin(p[5]);
  double c6 = cos(p[5]);

  R[4][0] = (T) (-s5 * c6);
  R[5][0] = (T) (-c5 * s6);

  R[4][1] = (T) (-s5 * s6);
  R[5][1] = (T) (c5 * c6);

  R[4][2] = (T) (-c5);

  R[3][4] = (T) (c4 * s5 * c6 + s4 * s6);
  R[4][4] = (T) (s4 * c5 * c6);
  R[5][4] = (T) (-s4 * s5 * s6 - c4 * c6);

  R[3][5] = (T) (c4 * s5 * s6 - s4 * c6);
  R[4][5] = (T) (s4 * c5 * s6);
  R[5][5] = (T) (s4 * s5 * c6 - c4 * s6);

  R[3][6] = (T) (c4 * c5);
  R[4][6] = (T) (-s4 * s5);

  R[3][8] = (T) (-s4 * s5 * c6 + c4 * s6);
  R[4][8] = (T) (c4 * c5 * c6);
  R[5][8] = (T) (-c4 * s5 * s6 + s4 * c6);

  R[3][9] = (T) (-s4 * s5 * s6 - c4 * c6);
  R[4][9] = (T) (c4 * c5 * s6);
  R[5][9] = (T) (c4 * s5 * c6 + s4 * s6);

  R[3][10] = (T) (-s4 * c5);
  R[4][10] = (T) (-c4 * s5);

  return R;
}

#endif
