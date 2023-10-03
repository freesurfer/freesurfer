/**
 * @brief A class to compute a registration using robust regression
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
// Sep. 4th ,2012
//
#ifndef RegRobust_H
#define RegRobust_H

#include "Registration.h"
#if ITK_VERSION_MAJOR >= 5
#include <iostream>
#include <vcl_compiler.h>
#else
#include <vcl_iostream.h>
#endif

// forward declaration
template<class T> class RegistrationStep;

/** \class RegRobust
 * \brief Class for robust registration 
 */
class RegRobust: public Registration
{
  template<class T> friend class RegistrationStep;
public:
  RegRobust() :
      Registration(), sat(-1), wlimit(0.16), mri_weights(NULL), mri_hweights(
          NULL), mri_indexing(NULL)
  {
  }
  
  virtual ~RegRobust();
  virtual MRI * getHalfWayGeom()
  {
    return mri_weights;
  }
  
  virtual void clear();
  //! Return the weights in target space
  virtual MRI * getWeights()
  {
    return mri_weights;
  }
  

  double estimateIScale(MRI *mriS, MRI *mriT);
  //! Estimate saturation parameter (higher sensitivity to outliers if sat small)
  double findSaturation();
  //! Set saturation parameter
  void setSaturation(double d)
  {
    sat = d;
  }
  
  //! Set weight limit for saturation estimation
  void setWLimit(double d)
  {
    wlimit = d;
  }

  //! Get Name of Registration class
  virtual std::string getClassName() {return "RegRobust";}
  
protected:
  virtual void computeIterativeRegistrationFull(int n, double epsit, MRI * mriS,
      MRI* mriT, const vnl_matrix<double> &Minit, double iscaleinit);
  //! To call the actual registration step
  template<class T> void iterativeRegistrationHelper(int nmax, double epsit,
      MRI * mriS, MRI* mriT, const vnl_matrix<double>& m, double scaleinit);

private:
  void findSatMultiRes(const vnl_matrix<double> &mi, double scaleinit);
  double wcheck; // set from computeRegistrationStepW
  double wchecksqrt; // set from computeRegistrationStepW

  // PRIVATE DATA
  double sat;
  double wlimit;
  MRI * mri_weights;
  MRI * mri_hweights;
  MRI * mri_indexing;

};

/** This is a template function to avoid code duplication.
 Computes iterative registration (recomputing A and b in each step).
 Retruns 4x4 matrix Mfinal and iscalefinal (class member).
 The caller needs to retrieve any really final transform with getFinalVox2Vox.
 */
template<class T>
void RegRobust::iterativeRegistrationHelper(int nmax, double epsit, MRI * mriS,
    MRI* mriT, const vnl_matrix<double>& m, double scaleinit)
{
  if (!mriS)
    mriS = mri_source;
  if (!mriT)
    mriT = mri_target;

  assert(mriS && mriT);
  assert(nmax>0);

  std::pair<vnl_matrix_fixed<double, 4, 4>, double> cmd(
      vnl_matrix_fixed<double, 4, 4>(), 0.0);
  std::pair<vnl_matrix_fixed<double, 4, 4>, double> fmd(
      vnl_matrix_fixed<double, 4, 4>(), 0.0);

  // check if mi (inital transform) is passed
  if (!m.empty())
    fmd.first = m;
  else if (!Minit.empty())
    fmd.first = getMinitResampled();
  else
    fmd.first = initializeTransform(mriS, mriT);

  // ISCALECHANGE:
  // intensity model: R(s,IS,IT) = exp(-0.5 s) IT - exp(0.5 s) IS
  //                  R'  = -0.5 ( exp(-0.5 s) IT + exp(0.5 s) IS)
  //  thus iscalefinal= exp(0.5 s) * (1/exp(-0.5 s)) = exp(s)
  //     ==> s = log (iscalefinal)
  iscalefinal = 1.0;
  if (scaleinit != 1.0)
    iscalefinal = scaleinit;
  else
    iscalefinal = iscaleinit;
  fmd.second = log(iscalefinal);

  if (verbose > 1)
  {
    std::cout << "   - initial iscale: " << iscalefinal << std::endl;
    std::cout << "   - initial transform:\n";
    vnl_matlab_print(vcl_cout,fmd.first,"Minit",vnl_matlab_print_format_long);
    std::cout << std::endl;
  }

  //std::cout << "mris width " << mriS->width << std::endl;
  MRI* mri_Swarp = NULL;
  MRI* mri_Twarp = NULL;
  vnl_matrix_fixed<double, 4, 4> mh2;
  vnl_matrix_fixed<double, 4, 4> mh;
  vnl_matrix_fixed<double, 4, 4> mhi;
  mhi.set_identity();
//  vnl_matrix_fixed<double , 4, 4> mi;

  double diff = 100.0;
  double idiff = 0.0;
  double ieps = 0.001; // exp(ieps) = 1.0010005, so stop if below 0.1% change

  int i = 1;

  // here create RegistrationStep of the specific type: double or float:
  RegistrationStep<T> RStep(*this);
  RStep.setFloatSVD(true); // even with double, use float for SVD to save some memory
  T tt;

  converged = false;
  while (!converged && i <= nmax)
  {
    if (verbose > 0)
      std::cout << " Iteration(" << typeid(tt).name() << "): " << i
          << std::endl;
    if (verbose == 1 && subsamplesize > 0)
      if (mriS->width > subsamplesize && mriS->height > subsamplesize
          && (mriS->depth > subsamplesize || mriS->depth == 1))
        std::cout << " (subsample " << subsamplesize << ") " << std::flush;
    if (verbose > 1)
      std::cout << std::endl;

    // map source (and if symmetric also target) to new space for next iteration
    // get mri_?warp and half way maps (mh, mhi):
    mapToNewSpace(fmd.first, iscalefinal, mriS, mriT, mri_Swarp, mri_Twarp, mh,
        mhi);

    if (debug > 0)
    {
      // write hw images before next registration step: 
      MRIwrite(mri_Swarp, (name + "-mriS-hw.mgz").c_str());
      MRIwrite(mri_Twarp, (name + "-mriT-hw.mgz").c_str());
      char ch;
      std::cout << "HW inputs written as" << std::endl;
      std::cout << (name + "-mriS-hw.mgz") << std::endl;
      std::cout << (name + "-mriT-hw.mgz") << std::endl;
      std::cout << "Press a key to run registration: ";
      std::cin >> ch;
    }

    // ==========================================================================
    //
    // compute Registration
    //
    if (verbose > 1)
      std::cout << "   - compute new registration" << std::endl;
    cmd = RStep.computeRegistrationStep(mri_Swarp, mri_Twarp);
    wcheck = RStep.getwcheck();
    wchecksqrt = RStep.getwchecksqrt();
    // ==========================================================================

    // store M and d
    if (verbose > 1)
    {
      std::cout << "   - recieved matrix update " << std::endl;
      vnl_matlab_print(vcl_cout,cmd.first,"Mupdate",vnl_matlab_print_format_long);
      std::cout << std::endl;
      std::cout << "   - store old transform" << std::endl;
      vnl_matlab_print(vcl_cout,fmd.first,"Mold",vnl_matlab_print_format_long);
      std::cout <<  std::endl;
    }
    vnl_matrix_fixed<double, 4, 4> fmdtmp(fmd.first);
    if (symmetry)
    {
      mh2 = vnl_inverse(mhi); // M = mh2 * mh
      // new M = mh2 * cm * mh
      fmd.first = (mh2 * cmd.first) * mh;
    }
    else
      fmd.first = cmd.first * fmd.first; // was '* mh' which should be the same
    if (verbose > 1)
    {
      std::cout << "   - updated full transform" << std::endl;
      vnl_matlab_print(vcl_cout,fmd.first,"Mnew",vnl_matlab_print_format_long);
      std::cout <<  std::endl;
    }

    // ISCALECHANGE:
    if (iscale)
    {
      fmd.second -= cmd.second; // adjust log
      iscalefinal = exp(fmd.second); // compute full factor (source to target)
      idiff = fabs(cmd.second);
      std::ostringstream istar;
      if (idiff <= ieps)
        istar << " <= " << ieps << "  :-)";
      if (verbose > 0)
        std::cout << "     -- intensity log diff: abs(" << cmd.second << ") "
            << istar.str() << std::endl;
    }

    if (!rigid)
      diff = MyMatrix::getFrobeniusDiff(fmd.first.as_matrix(), fmdtmp.as_matrix());
    else
      diff = sqrt(MyMatrix::RigidTransDistSq(fmd.first.as_matrix(), fmdtmp.as_matrix()));
    if (verbose > 1)
      std::cout << "     -- old diff. to prev. transform: " << diff
          << std::endl;
    diff = sqrt(MyMatrix::AffineTransDistSq(fmd.first.as_matrix(), fmdtmp.as_matrix(), 100));
    std::ostringstream star;
    if (diff <= epsit)
      star << "  <= " << epsit << "   :-)";
    else if (i == nmax)
      star << " max it: " << nmax << " reached!";
    if (verbose > 0)
      std::cout << "     -- diff. to prev. transform: " << diff << star.str()
          << std::endl;
    //std::cout << " intens: " << fmd.second << std::endl;
    i++;
    converged = (diff <= epsit && idiff <= ieps);

    if (debug > 0)
    {
      // write weights
      if (costfun == ROB)
      {
        if (mri_hweights)
          MRIfree(&mri_hweights);
        assert(RStep.getWeights());
        mri_hweights = MRIcopy(RStep.getWeights(), NULL);
        MRIwrite(mri_hweights, (name + "-mriS-hw-weights.mgz").c_str());
        char ch;
        std::cout << "Weights written for this iteration" << std::endl;
        std::cout << "Press a key to continue iterations: ";
        std::cin >> ch;
      }
    } // end if debug

  } // end while loop

  // adjust converged (don't tell the outside about intensiy problems)
  converged = (diff <= epsit);

  if (mri_hweights)
    MRIfree(&mri_hweights);
  if (costfun == ROB)
  {
    assert(RStep.getWeights());
    mri_hweights = MRIcopy(RStep.getWeights(), NULL);
  }

  //   DEBUG OUTPUT
  if (debug > 0)
  {
    // write weights and warped images after last step:

    MRIwrite(mri_Swarp, (name + "-mriS-mapped.mgz").c_str());
    MRIwrite(mri_Twarp, (name + "-mriT-mapped.mgz").c_str());
    MRI* salign = MRIclone(mriS, NULL);
    salign = MyMRI::MRIlinearTransform(mri_Swarp, salign, cmd.first);
    MRIwrite(salign, (name + "-mriS-align.mgz").c_str());
    MRIfree(&salign);
  } // end if debug

  // store weights (mapped to target space):
  if (mri_weights)
    MRIfree(&mri_weights);
  if (mri_hweights)
  {
    if (debug > 0)
    {
      // in the half-way space:
      std::string n = name + std::string("-mriS-weights.mgz");
      MRIwrite(mri_hweights, n.c_str());
    }
    // remove negative weights (markers) set to 1
   /* int x, y, z;
    for (z = 0; z < mri_hweights->depth; z++)
      for (x = 0; x < mri_hweights->width; x++)
        for (y = 0; y < mri_hweights->height; y++)
        {
          if (MRIFvox(mri_hweights,x,y,z) < 0) MRIFvox(mri_hweights,x,y,z) = 1;
        }*/

    mri_weights = MRIalloc(mriT->width, mriT->height, mriT->depth, MRI_FLOAT);
    MRIcopyHeader(mriT, mri_weights);
    mri_weights->type = MRI_FLOAT;
    mri_weights->outside_val = 1; 
    if (symmetry && ! iscaleonly)
      mri_weights = MyMRI::MRIlinearTransform(mri_hweights, mri_weights, mh2);
    else
      mri_weights = MRIcopy(mri_hweights, mri_weights);
  }

  //vnl_matlab_print(vcl_cerr,fmd.first,"fmd",vnl_matlab_print_format_long);std::cerr << std::endl;
  //vnl_matlab_print(vcl_cerr,cmd.first,"cmd",vnl_matlab_print_format_long);std::cerr << std::endl;
  //vnl_matlab_print(vcl_cerr,mh,"mov1hw",vnl_matlab_print_format_long);std::cerr << std::endl;
  //vnl_matlab_print(vcl_cerr,mhi,"dst1hw",vnl_matlab_print_format_long);std::cerr << std::endl;

  // adjust half way maps to new midpoint based on final transform
  if (verbose > 1)
    std::cout << "     -- adjusting half-way maps " << std::endl;
  vnl_matrix_fixed<double, 4, 4> ch = MyMatrix::MatrixSqrt(fmd.first.as_matrix());
  // do not just assume c = ch*ch, rather c = ch2 * ch
  // for transforming target we need ch2^-1 = ch * c^-1
  vnl_matrix_fixed<double, 4, 4> ci = vnl_inverse(fmd.first);
  vnl_matrix_fixed<double, 4, 4> chi = ch * ci;
  mov2weights = ch.as_matrix();
  dst2weights = chi.as_matrix();

  //vnl_matlab_print(vcl_cerr,mov2weights,"mov2hw",vnl_matlab_print_format_long);std::cerr << std::endl;
  //vnl_matlab_print(vcl_cerr,dst2weights,"dst2hw",vnl_matlab_print_format_long);std::cerr << std::endl;

  MRIfree(&mri_Twarp);
  MRIfree(&mri_Swarp);

  Mfinal = fmd.first.as_matrix();

}

#endif

