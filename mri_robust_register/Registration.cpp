/**
 * @file Registration.cpp
 * @brief A class to compute a robust symmetric registration
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/08/17 22:55:21 $
 *    $Revision: 1.31 $
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

#include "Registration.h"
#include "Quaternion.h"
#include "MyMatrix.h"
#include "MyMRI.h"
#include "Regression.h"
#include "CostFunctions.h"

#include <limits>
#include <cassert>
#include <fstream>
#include <sstream>

#ifdef __cplusplus
extern "C"
{
#endif
#include "limits.h"
#include "error.h"
#include "macros.h"
#include "mrimorph.h"

#ifdef __cplusplus
}
#endif

using namespace std;

Registration::~Registration()
{ // we cleanup our private variables
  //std::cout << " Destroy Registration" << std::endl;
  if (mri_source) MRIfree(&mri_source);
  if (mri_target) MRIfree(&mri_target);
  if (Minit) MatrixFree(&Minit);
  if (Mfinal) MatrixFree(&Mfinal);
  if (lastp) MatrixFree(&lastp);
  if (mri_indexing) MRIfree(&mri_indexing);
  if (mri_weights) MRIfree(&mri_weights);
  if (gpS.size() > 0) freeGaussianPyramid(gpS);
  if (gpT.size() > 0) freeGaussianPyramid(gpT);
  if (mov2weights) MatrixFree(&mov2weights);
  if (dst2weights) MatrixFree(&dst2weights);
  //std::cout << " Done " << std::endl;
}


void Registration::clear() // initialize registration (keep source and target and gauss pyramid)
// initialize registration (keep source and target and gauss pyramid)
{
  transonly = false;
  rigid = true;
  robust = true;
  sat = -1;
  iscale = false;
  rtype = 1;
  subsamplesize = -1;
  initorient = false;
  debug = 0;
  if (Minit) MatrixFree(&Minit);
  if (Mfinal) MatrixFree(&Mfinal);
  if (lastp) MatrixFree(&lastp);
  if (mri_indexing) MRIfree(&mri_indexing);
  if (mri_weights) MRIfree(&mri_weights);
  mri_weights= NULL;
  if (mov2weights) MatrixFree(&mov2weights);
  mov2weights = NULL;
  if (dst2weights) MatrixFree(&dst2weights);
  dst2weights =NULL;
  if (mri_weights) MRIfree(&mri_weights);
}

pair < MATRIX*, MRI* > Registration::computeRegistrationStepW(MRI * mriS, MRI* mriT)
// computes Registration single step
// the mri's have to be in same space
// returns parameter vector and float MRI with the weights (if robust, else weights ==NULL)
//
// if rigid (only trans and rot) else use affine transform
// if robust  use M-estimator instead of ordinary least squares
// if iscale add parameter for intensity scaling
// rtype only for rigid (2: affine restriction to rigid, 1: use rigid from robust-paper)
{

//  MATRIX *A = constructA_rt(mriS, NULL);
//  MATRIX *b = constructb(mriS, mriT, NULL);


  //cout << " compute Registration!"<< endl;
  //cout << "   - construct AB: " << endl;
  pair < MATRIX*, VECTOR*> Ab;
//  // cannot compute restriction without last p
//  if (rtype == 2 && lastp == NULL) rtype = 1;
//  cout << " rtype: " << rtype << endl;
  if (rigid && rtype==2)
  {
    if (verbose > 1) cout << "rigid and rtype 2 !" << endl;
    // compute non rigid A
    rigid = false;
    Ab = constructAb(mriS,mriT);
    rigid = true;
    // now restrict A  (= A R(lastp) )
    MATRIX* R;
    if (lastp) R = constructR(lastp);
    else
    {
      int l = 6;
      if (!rigid) l = 12;
      if (iscale) l++;
      MATRIX* tm = MatrixAlloc(l,1,MATRIX_REAL);
      MatrixClear(tm);
      R = constructR(tm);
      MatrixFree(&tm);
      //MatrixPrintFmt(stdout,"% 2.8f",R);exit(1);

    }
    MATRIX* Rt = MatrixTranspose(R,NULL);
    MATRIX* nA = MatrixMultiply(Ab.first,Rt,NULL);
    //MatrixPrintFmt(stdout,"% 2.8f",nA);exit(1);
    MatrixFree(&Ab.first);
    Ab.first = nA;
    MatrixFree(&R);
    MatrixFree(&Rt);
  }
  else
  {
    //cout << "Rtype  " << rtype << endl;

    Ab = constructAb(mriS,mriT);
  }

  for (int rr = 1; rr<= Ab.first->rows; rr++)
  {
    for (int cc = 1; cc <= Ab.first->cols; cc++)
    {
      assert (!isnan(*MATRIX_RELT(Ab.first, rr, cc)));
    }
    assert (!isnan(*MATRIX_RELT(Ab.second, rr, 1)));
  }
  // DEBUG:
  //MatrixWriteTxt("A.txt",Ab.first);
  //MatrixWriteTxt("b.txt",Ab.second);

  pair < MATRIX*, MATRIX* > pwm(NULL,NULL);
  pair < MATRIX*, MRI* > pw(NULL,NULL);
  Regression R(Ab.first,Ab.second);
	R.setVerbose(verbose);
  if (robust)
  {
    if (verbose > 1) cout << "   - compute robust estimate ( sat "<<sat<<" )..." << flush;
    if (sat < 0) pwm = R.getRobustEstW();
    else pwm = R.getRobustEstW(sat);
    if (verbose > 1) cout << "  DONE" << endl;
    pw.first  = pwm.first;
    pw.second = NULL;

//    cout << " pw-final  : "<< endl;
//    MatrixPrintFmt(stdout,"% 2.8f",pw.first);

    // transform weights vector back to 3d (mri real)
    pw.second = MRIalloc(mriS->width, mriS->height, mriS->depth, MRI_FLOAT);
    MRIcopyHeader(mriS, pw.second) ;
    pw.second->type = MRI_FLOAT;
    MRIsetResolution(pw.second, mriS->xsize, mriS->ysize, mriS->zsize);
    int x,y,z;
    long int count = 1, val;
    for (z = 0 ; z < mriS->depth  ; z++)
      for (x = 0 ; x < mriS->width  ; x++)
        for (y = 0 ; y < mriS->height ; y++)
        {
          val = MRILvox(mri_indexing,x,y,z);
          if (val == 0) MRIFvox(pw.second, x, y, z) = -0.5;
          else if (val == -1) MRIFvox(pw.second, x, y, z) = -1;
          else
          {
            //cout << val << "  xyz: " << x << " " << y << " " << z << " " << flush;
            assert(val <= pwm.second->rows);
            MRIFvox(pw.second, x, y, z) = *MATRIX_RELT(pwm.second,val , 1);
            //cout << "d"<<*MATRIX_RELT(pwm.second,val , 1)<< " " << MRIFvox(pw.second, x, y, z)<< endl;
            count++;
          }
        }
    //cout << endl;
    //cout << " count-1: " << count-1 << " rows: " << pwm.second->rows << endl;
    assert(count-1 == pwm.second->rows);

    if (pwm.second != NULL) MatrixFree(&pwm.second);
  }
  else
  {
    if (verbose > 1) cout << "   - compute least squares estimate ..." << flush;
    pw.first = R.getLSEst();
    if (verbose > 1) cout << "  DONE" << endl;
    pw.second = NULL; // no weights in this case
  }

//  zeroweights = R.getLastZeroWeightPercent();
  zeroweights = R.getLastWeightPercent();

//  R.plotPartialSat(name);

  MatrixFree(&Ab.first);
  MatrixFree(&Ab.second);

//    cout << " pw-final  : "<< endl;
//    MatrixPrintFmt(stdout,"% 2.8f",pw.first);

  return pw;
}

MATRIX* Registration::computeRegistrationStepP(MRI * mriS, MRI* mriT)
// computes Registration single step
// retruns parameter Vector only
{
  pair < MATRIX*, MRI*> pw = computeRegistrationStepW(mriS,mriT);
  if (pw.second)
  {
    //cout << "mrisweight widht " << pw.second->width << endl;
    if (debug > 0)
    {
      string n = name+string("-mriS-weights.mgz");
      MRIwrite(pw.second,n.c_str());
    }
    MRIfree(&pw.second);
  }
  return pw.first;
}


pair <MATRIX*,double> Registration::computeRegistrationStep(MRI * mriS, MRI* mriT)
// computes Registration single step
// retruns 4x4 matrix and iscale value
{
//cout << "  Registration::computeRegistrationStep " << endl;
  MATRIX* p = computeRegistrationStepP(mriS,mriT);
//   cout << " prows: " << p->rows << endl;
  pair <MATRIX*, double> pd = convertP2Md(p);
  MatrixFree(&p);
  return pd;
}


// pair < MATRIX*, double > Registration::computeIterativeRegistration( int n, MRI * mriS, MRI* mriT, MATRIX* m, double scaleinit)
// // computes iterative registration (recomputing A and b in each step)
// // retruns 4x4 matrix and iscale value
// {
//  //  if (!mriS) mriS = mri_source;
//  //  if (!mriT) mriT = mri_target;
//
//    assert (mriS && mriT);
//
//    int MAX = n;
//    int i;
//    pair < MATRIX*, double > cmd(NULL,1.0);
//    pair < MATRIX*, double > fmd(NULL,scaleinit);
//    if (m) fmd.first = MatrixCopy(m,NULL);
//    else if (Minit) fmd.first = MatrixCopy(Minit,NULL);
//    else fmd.first = initializeTransform(mriS,mriT);
//
//   cout << "   - initial transform:\n" ;
//   MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
//
//    //cout << "mris widht " << mriS->width << endl;
//    MRI* mri_Swarp   = NULL;
//    MRI* mri_Twarp   = NULL;
//    MATRIX* p;
//
//    lastp = NULL;
//
//    for (i=1;i<=MAX;i++)
//    {
//       cout << " Iteration: " << i << endl;
//
//
//        // warp source to target
// //       cout << "   - warping source" << endl;
// //       mri_Swarp = MRIlinearTransform(mriS, mri_Swarp, fmd.first);
// //       mri_Twarp = MRIcopy(mriT,mri_Twarp);
//
//  // here maybe better to symmetrically warp both images SQRT(M)
//  // this keeps the problem symmetric
//          cout << "   - warping source and target (sqrt)" << endl;
//          MATRIX * mh  = MatrixSqrt(fmd.first);
//   // do not just assume m = mh*mh, instead compute B = mh^-1 m
//   //  (in fact we need B^-1= m^-1 mh for transforming target):
//   MATRIX * mi  = MatrixInverse(fmd.first,NULL);
//   MATRIX * mhi = MatrixMultiply(mi,mh,NULL);
//
//          if (mri_Swarp) MRIfree(&mri_Swarp);
//          mri_Swarp =  MRIlinearTransform(mriS,NULL, mh);
//          if (mri_Twarp) MRIfree(&mri_Twarp);
//          mri_Twarp =  MRIlinearTransform(mriT,NULL, mhi);
//   MatrixFree(&mh);
//   MatrixFree(&mi);
//   MatrixFree(&mhi);
//
//        // adjust intensity
//        if (iscale)
//        {
//           cout << "   - adjusting intensity ( "<< fmd.second << " ) " << endl;
//           //MRIvalscale(mri_Swarp,mri_Swarp,fmd.second);
//           MRIvalscale(mri_Swarp,mri_Swarp,(1.0+fmd.second)*0.5);
//           MRIvalscale(mri_Twarp,mri_Twarp,(1.0+ 1.0/fmd.second)*0.5);
//        }
//
//
//        // compute Registration
//        cout << "   - compute new registration" << endl;
//  //cout << "mriswarp widht " << mri_Swarp->width << endl;
//        p = computeRegistrationStepP(mri_Swarp,mri_Twarp);
//        cmd = convertP2Md(p);
//        if (lastp) MatrixFree(&lastp);
//        lastp = p;
//        p = NULL;
//
// //   DEBUG
//      MRIwrite(mri_Swarp,"mriS-warp.mgz");
//      MRIwrite(mri_Twarp,"mriT-warp.mgz");
//      MRI* salign = MRIlinearTransform(mri_Swarp, NULL,cmd.first);
//      MRIwrite(salign,"mriS-align.mgz");
//      MRIfree(&salign);
//
//        // store M and d
//        cout << "   - store transform" << endl;
//      //cout << endl << " current : Matrix: " << endl;
//      //MatrixPrintFmt(stdout,"% 2.8f",cmd.first);
//      //cout << " intens: " << cmd.second << endl;
//      MATRIX* fmdtmp = MatrixCopy(fmd.first,NULL);
//        fmd.first = MatrixMultiply(cmd.first,fmd.first,fmd.first);
//        fmd.second *= cmd.second;
//        if (cmd.first != NULL) MatrixFree(&cmd.first);
//      //cout << endl << " Matrix: " << endl;
//      //MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
//      cout << "       difference to prev. transform: " << getFrobeniusDiff(fmd.first, fmdtmp) << endl;
//      MatrixFree(&fmdtmp);
//      //cout << " intens: " << fmd.second << endl;
//
//    }
//
//   MRIfree(&mri_Twarp);
//   MRIfree(&mri_Swarp);
//
//   Mfinal = MatrixCopy(fmd.first, Mfinal);
//   iscalefinal = fmd.second;
//
//    return fmd;
// }

pair < MATRIX*, double > Registration::computeIterativeRegistration( int nmax,double epsit, MRI * mriS, MRI* mriT, MATRIX* m, double scaleinit)
// computes iterative registration (recomputing A and b in each step)
// retruns 4x4 matrix and iscale value
{
  if (!mriS) mriS = mri_source;
  if (!mriT) mriT = mri_target;

  assert (mriS && mriT);
  assert (nmax>0);

  pair < MATRIX*, double > cmd(NULL,1.0);
  pair < MATRIX*, double > fmd(NULL,scaleinit);

  // check if mi (inital transform) is passed
  if (m) fmd.first = MatrixCopy(m,NULL);
  else if (Minit) fmd.first = MatrixCopy(Minit,NULL);
  else fmd.first = initializeTransform(mriS,mriT);

  if (verbose > 1)
  {
    cout << "   - initial transform:\n" ;
    MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
  }

  //cout << "mris widht " << mriS->width << endl;
  MRI* mri_Swarp   = NULL;
  MRI* mri_Twarp   = NULL;
  pair < MATRIX*, MRI*> pw(NULL,NULL);
  MATRIX * mh2 = NULL;
  MATRIX * mh  = NULL;
  MATRIX * mhi = NULL;

  if (lastp) MatrixFree(&lastp);
  lastp = NULL;
  double diff = 100;
  int i = 1;
  MATRIX* mras = NULL;
  while (diff > epsit && i<=nmax)
  {
    if (verbose >0) cout << " Iteration: " << i << flush;
		if (verbose == 1 && subsamplesize > 0)
		  if (mriS->width > subsamplesize || mriS->height > subsamplesize || mriS->depth > subsamplesize)
			  cout << " (subsample " << subsamplesize << ") " << flush;
		if (verbose >1) cout << endl;

    // warp source to target
//       cout << "   - warping source" << endl;
//       mri_Swarp = MRIlinearTransform(mriS, mri_Swarp, fmd.first);
//       mri_Twarp = MRIcopy(mriT,mri_Twarp);

    // here maybe better to symmetrically warp both images SQRT(M)
    // this keeps the problem symmetric
    if (verbose >1) cout << "   - warping source and target (sqrt)" << endl;
    if (mh) MatrixFree(&mh);
    // Old: half way voxelxform
    mh  = MyMatrix::MatrixSqrt(fmd.first);
    // do not just assume m = mh*mh, rather m = mh2 * mh
    // for transforming target we need mh2^-1 = mh * m^-1
    MATRIX * mi  = MatrixInverse(fmd.first,NULL);
    //MATRIX * mhi = MatrixMultiply(mi,mh,NULL); //old
    mhi = MatrixMultiply(mh,mi,mhi);

//     // get half way rasxform: (not working : low resolution?)
//     mras = MRIvoxelXformToRasXform (mriS,mriT,fmd.first,mras);
//     mh   = MatrixSqrt(mras);
//     // do not just assume m = mh*mh, rather m = mh2 * mh
//     // for transforming target we need mh2^-1 = mh * m^-1
//     MATRIX * mi  = MatrixInverse(mras,NULL);
//     //MATRIX * mhi = MatrixMultiply(mi,mh,NULL); //old
//     mhi = MatrixMultiply(mh,mi,mhi);
//     // make vox2vox
//     mh  = MRIrasXformToVoxelXform (mriS,mriS,mh,mh);
//     mhi = MRIrasXformToVoxelXform (mriT,mriT,mhi,mhi);

    if (mri_Swarp) MRIfree(&mri_Swarp);
    mri_Swarp = MRIclone(mriS,NULL);
    mri_Swarp = MRIlinearTransform(mriS,mri_Swarp, mh);
    if (mri_Twarp) MRIfree(&mri_Twarp);
    mri_Twarp = MRIclone(mriS,NULL); // bring them to same space (just use src geometry)
    mri_Twarp = MRIlinearTransform(mriT,mri_Twarp, mhi);

    // adjust intensity
    if (iscale)
    {
      if (verbose >1) cout << "   - adjusting intensity ( "<< fmd.second << " ) " << endl;
      //MRIvalscale(mri_Swarp,mri_Swarp,fmd.second);
      MyMRI::MRIvalscale(mri_Swarp,mri_Swarp,(1.0+fmd.second)*0.5);
      MyMRI::MRIvalscale(mri_Twarp,mri_Twarp,(1.0+ 1.0/fmd.second)*0.5);
    }


    // compute Registration
    if (verbose >1) cout << "   - compute new registration" << endl;
    //p = computeRegistrationStepP(mri_Swarp,mri_Twarp);
    if (pw.second) MRIfree(&pw.second);
    pw = computeRegistrationStepW(mri_Swarp,mri_Twarp);

    if (cmd.first != NULL) MatrixFree(&cmd.first);
    cmd = convertP2Md(pw.first);
    if (lastp) MatrixFree(&lastp);
    lastp = pw.first;
    pw.first = NULL;


    // store M and d
    if (verbose >1) cout << "   - store transform" << endl;
    //cout << endl << " current : Matrix: " << endl;
    //MatrixPrintFmt(stdout,"% 2.8f",cmd.first);
    //cout << " intens: " << cmd.second << endl;
    MATRIX* fmdtmp = MatrixCopy(fmd.first,NULL);
    //fmd.first = MatrixMultiply(cmd.first,fmd.first,fmd.first); //old
    if (mh2) MatrixFree(&mh2);
    mh2 = MatrixInverse(mhi,NULL); // M = mh2 * mh
    // new M = mh2 * cm * mh
    fmd.first = MatrixMultiply(mh2,cmd.first,fmd.first);
    fmd.first = MatrixMultiply(fmd.first,mh,fmd.first);

    fmd.second *= cmd.second;
    //cout << endl << " Matrix: " << endl;
    //MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
    if (!rigid) diff = MyMatrix::getFrobeniusDiff(fmd.first, fmdtmp);
    else        diff = sqrt(MyMatrix::RigidTransDistSq(fmd.first, fmdtmp));
    if (verbose >1) cout << "     -- old diff. to prev. transform: " << diff << endl;
    diff = sqrt(MyMatrix::AffineTransDistSq(fmd.first, fmdtmp, 100));
		ostringstream star;
		if (diff < epsit) star << "  < " << epsit << "  :-)" ;
		else if (i == nmax) star<< " max it: " << nmax << " reached!";
    if (verbose >0) cout << "     -- diff. to prev. transform: " << diff << star.str() << endl;
    //cout << " intens: " << fmd.second << endl;

    MatrixFree(&fmdtmp);
    MatrixFree(&mi);
    //MatrixFree(&mh);
    //MatrixFree(&mhi);
    //MatrixFree(&mh2);
		i++;

  } // end while loop
  if (mras) MatrixFree(&mras);
  
  //   DEBUG OUTPUT
  if (debug > 0)
  {
    // write weights and warped images after last step:

    MRIwrite(mri_Swarp,(name+"-mriS-warp.mgz").c_str());
    MRIwrite(mri_Twarp,(name+"-mriT-warp.mgz").c_str());
    MRI* salign = MRIclone(mriS,NULL);
    salign = MRIlinearTransform(mri_Swarp, salign,cmd.first);
    MRIwrite(salign,(name+"-mriS-align.mgz").c_str());
    MRIfree(&salign);
    if (pw.second)
    {
      // in the half-way space:
      string n = name+string("-mriS-weights.mgz");
      MRIwrite(pw.second,n.c_str());
    }
  } // end if debug

  // store weights (mapped to target space):
  if (pw.second)
  {
    // remove negative weights (markers) set to 1
    int x,y,z;
    for (z = 0 ; z < pw.second->depth  ; z++)
      for (x = 0 ; x < pw.second->width  ; x++)
        for (y = 0 ; y < pw.second->height ; y++)
        {
          if (MRIFvox(pw.second,x,y,z) < 0) MRIFvox(pw.second,x,y,z) = 1;
        }
    MRI * mtmp = MRIalloc(mriT->width,mriT->height,mriT->depth,MRI_FLOAT);
    MRIcopyHeader(mriT,mtmp);
    mtmp->type = MRI_FLOAT;
    mtmp = MRIlinearTransform(pw.second,mtmp,mh2);
    //MRIwrite(mtmp,weightsname.c_str());
    MRIfree(&pw.second);
    pw.second = mtmp;
  }
  if (mri_weights) MRIfree(&mri_weights);
  mri_weights = pw.second;
  if (mov2weights) MatrixFree(&mov2weights);
  if (dst2weights) MatrixFree(&dst2weights);
  mov2weights = mh; // no freeing needed
  dst2weights = mhi;

  if (diff > epsit) // adjust mh and mhi to new midpoint
  {
    if (verbose >1) cout << "     -- adjusting half-way maps " << endl;
    MATRIX * ch = MyMatrix::MatrixSqrt(fmd.first);
    // do not just assume c = ch*ch, rather c = ch2 * ch
    // for transforming target we need ch2^-1 = ch * c^-1
    MATRIX * ci  = MatrixInverse(cmd.first,NULL);
    MATRIX * chi = MatrixMultiply(ch,ci,NULL);
    // append ch or chi to mh mhi
    mov2weights = MatrixMultiply(ch,mh,NULL);
    dst2weights = MatrixMultiply(chi,mhi,NULL);
    MatrixFree(&mh);
    MatrixFree(&mhi);
  }

  MRIfree(&mri_Twarp);
  MRIfree(&mri_Swarp);
  MatrixFree(&mh2);
  if (cmd.first != NULL) MatrixFree(&cmd.first);


  Mfinal = MatrixCopy(fmd.first, Mfinal);
  iscalefinal = fmd.second;

  return fmd;
}

pair < MATRIX*, double > Registration::computeIterativeRegSat( int n,double epsit, MRI * mriS, MRI* mriT, MATRIX* m, double scaleinit)
// tests trough many saturations:
{
  if (!mriS) mriS = mri_source;
  if (!mriT) mriT = mri_target;

  assert (mriS && mriT);

  pair < MATRIX*, double > cmd(NULL,1.0);
  pair < MATRIX*, double > fmd(NULL,scaleinit);

  // check if mi (inital transform) is passed
  if (m)          fmd.first = MatrixCopy(m,NULL);
  else if (Minit) fmd.first = MatrixCopy(Minit,NULL);
  else            fmd.first = initializeTransform(mriS,mriT) ;

  if (verbose >1) 
  {
    cout << "   - initial transform:\n" ;
    MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
  }

  double satval[20] = {20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
  vector < double > diffs (20);

  string nametmp = name;
  std::stringstream sout;

  for (int si = 0;si<(int)diffs.size();si++)
  {
    sat = satval[si];

    std::stringstream out;
    out << sat;
    name = nametmp+"-sat"+out.str();

    cmd = computeIterativeRegistration(n,epsit,mriS,mriT,fmd.first,fmd.second);

    if (!rigid) diffs[si] = MyMatrix::getFrobeniusDiff(fmd.first, cmd.first);
    else        diffs[si] = sqrt(MyMatrix::RigidTransDistSq(fmd.first, cmd.first));
      if (verbose >1) cout << "       difference on sat " << sat << " to prev. transform: " << diffs[si] << endl;

    fmd.second = cmd.second;
    MatrixFree(&fmd.first);
    fmd.first = cmd.first;

    // store transform
    LTA * lta = LTAalloc(1,mriS);
    lta->xforms[0].m_L = MRIvoxelXformToRasXform (mriS, mriT, fmd.first, lta->xforms[0].m_L) ;
    lta->type = LINEAR_RAS_TO_RAS ;
    getVolGeom(mriS, &lta->xforms[0].src);
    getVolGeom(mriT, &lta->xforms[0].dst);
    LTAwriteEx(lta, (name+".lta").c_str()) ;
    LTAfree(&lta);


  }

  name = nametmp;

  // plot diffs
  string fbase = name+"-sat";
  ofstream ofile((fbase+".plot").c_str(),ios::out);
  bool png = false;
  if (png) ofile << "set terminal png medium size 800,600" << endl;
  else ofile << "set terminal postscript eps color" << endl;
  if (png) ofile << "set output \""<< fbase <<".png\"" << endl;
  else ofile << "set output \""<< fbase <<".eps\"" << endl;
  ofile << "plot ";
  ofile << " \"-\" notitle with lines 1" << endl;
  for (int j = 0;j<(int)diffs.size(); j++)
  {
    ofile << -satval[j] << " " << diffs[j] << endl;
  }
  ofile << "e" << endl;
  ofile.close();


  return fmd;
}

pair < MATRIX*, double> Registration::computeMultiresRegistration (int stopres, int n,double epsit, MRI * mriS, MRI* mriT, MATRIX* mi, double scaleinit )
// stopres : stops on this resolution level (0 highest resoltuion ...)
// n: number of max iterations on each resolution
// epsit: epsilon to stop iterations
{
  if (verbose >0) cout << endl << endl << " Registration::computeMultiresRegistration " << endl;
  if (verbose >0) cout << "   - Gaussian Pyramid " << endl;
  if (!mriS) mriS = mri_source;
  else
  {
    if (gpS.size() > 0) freeGaussianPyramid(gpS);
    gpS = buildGaussianPyramid(mriS,100);
  }
  if (!mriT) mriT = mri_target;
  else
  {
    if (gpT.size() > 0) freeGaussianPyramid(gpT);
    gpT = buildGaussianPyramid(mriT,100);
  }

  if (gpS.size() ==0) gpS = buildGaussianPyramid(mriS,100);
  if (gpT.size() ==0) gpT = buildGaussianPyramid(mriT,100);
  assert(gpS.size() == gpT.size());

  int resolution = gpS.size();

// debug : save pyramid
//  for (uint i = 0;i<gpS.size();i++)
//  {
//   char fn[40];
//   sprintf(fn, "pyramid-%d.mgz", i+1);
//   MRIwrite(gpS[i],fn);
//  }

  MATRIX *m  = MatrixIdentity(4,NULL);

  // variables to store matrix m and scaling factor d:
  pair < MATRIX* , double > cmd;
  pair < MATRIX* , double > md(NULL,scaleinit);

  // check if mi (inital transform) is passed
  if (mi) md.first = MatrixCopy(mi,NULL);
  else if (Minit) md.first = MatrixCopy(Minit,NULL);
  else md.first = initializeTransform(mriS,mriT);

  if (verbose >0 ) 
  {
    cout << "   - initial transform:\n" ;
    MatrixPrintFmt(stdout,"% 2.8f",md.first);
		cout << "   - initial iscale: " << scaleinit <<endl;
  }

  // adjust minit to current (lowest) resolution:
  int rstart = 2;
  for (int r = 1; r<=resolution-rstart; r++)
    for (int rr = 1;rr<=3;rr++)
      md.first->rptr[rr][4]  = 0.5 *  md.first->rptr[rr][4];

  if (verbose >1 ) 
  {
    cout << "   - initial adjusted:\n" ;
    MatrixPrintFmt(stdout,"% 2.8f",md.first);
  }


//     MRI* falign = MRIlinearTransform(gpS[resolution-rstart], NULL,md.first);
//     MRIwrite(falign,"mriS-lowres-initaligned.mgz");
//     MRIfree(&falign);
//     MRIwrite(gpT[resolution-rstart],"mriT-lowres.mgz");


//  md.second = 1; //??
//  MRI* mri_Swarp   = NULL;
//  MRI* mri_Twarp   = NULL;

  bool iscaletmp = iscale;
  iscale = false; //disable intensity scaling on low resolutions

  for (int r = resolution-rstart;r>=stopres;r--)
  {
//    MRIwrite(gpS[r],"mriS-smooth.mgz");
//    MRIwrite(gpT[r],"mriT-smooth.mgz");

    if (verbose >0 ) cout << endl << "Resolution: " << r << endl;

    if (r<=1) iscale = iscaletmp; // set iscale if set by user

//        if (transonly)
//        {
//           MATRIX * mdh = MatrixCopy(md.first,NULL);
//    mdh->rptr[1][4] = 0.5 *mdh->rptr[1][4] ;
//    mdh->rptr[2][4] = 0.5 *mdh->rptr[2][4];
//    mdh->rptr[3][4] = 0.5 *mdh->rptr[3][4];
//
//    MATRIX * mdhi= MatrixCopy(md.first,NULL);
//    mdhi->rptr[1][4] = 0.5 *mdhi->rptr[1][4] ;
//    mdhi->rptr[2][4] = 0.5 *mdhi->rptr[2][4];
//    mdhi->rptr[3][4] = 0.5 *mdhi->rptr[3][4];
//
//           cout << "   - warping source and target (trans)" << endl;
//           if (mri_Swarp) MRIfree(&mri_Swarp);
//           mri_Swarp =  MRIlinearTransform(gpS[r],NULL, mdh);
//           if (mri_Twarp) MRIfree(&mri_Twarp);
//           mri_Twarp =  MRIlinearTransform(gpT[r],NULL, mdhi);
//    MatrixFree(&mdh);
//    MatrixFree(&mdhi);
//        }
//        else if (rigid)
//        {

//          //  cout << "   - warping source to target (rigid)" << endl;
//          // if (mri_Swarp) MRIfree(&mri_Swarp);
//          // mri_Swarp =  MRIlinearTransform(gpS[r],NULL, md.first);
//          // if (mri_Twarp) MRIfree(&mri_Twarp);
//          // mri_Twarp =  MRIcopy(gpT[r],NULL);

//
//           cout << "   - warping source and target (rigid)" << endl;
//           MATRIX * mh   = getHalfRT(md.first);
//           MATRIX * mi  = MatrixInverse(md.first,NULL);
//    //MATRIX * mhi  = MatrixInverse(mh,NULL);
//      MATRIX * mhi = MatrixMultiply(mi,mh,NULL);
//           if (mri_Swarp) MRIfree(&mri_Swarp);
//           mri_Swarp =  MRIlinearTransform(gpS[r],NULL, mh);
//           if (mri_Twarp) MRIfree(&mri_Twarp);
//           mri_Twarp =  MRIlinearTransform(gpT[r],NULL, mhi);
//    MatrixFree(&mh);
//    MatrixFree(&mi);
//    MatrixFree(&mhi);

//         }
//        else // affine
//        {
//          // warp source to target
//          // !!! here maybe better to symmetrically warp both images SQRT(M)!!!
//          MATRIX * mh = MatrixSqrt(md.first);
//    MATRIX * mi  = MatrixInverse(md.first,NULL);
//    MATRIX * mhi = MatrixMultiply(mi,mh,NULL);
//
//          cout << "   - warping source and target (sqrt)" << endl;
//          if (mri_Swarp) MRIfree(&mri_Swarp);
//          mri_Swarp =  MRIlinearTransform(gpS[r],NULL, mh);
//          if (mri_Twarp) MRIfree(&mri_Twarp);
//           mri_Twarp =  MRIlinearTransform(gpT[r],NULL, mhi);
//    MatrixFree(&mh);
//    MatrixFree(&mhi);
//    MatrixFree(&mi);
//      }

//        // adjust intensity
//        if (iscale)
//        {
//           cout << "   - adjusting intensity ( "<< md.second << " ) " << endl;
//           MRIvalscale(mri_Swarp,mri_Swarp,(1.0+md.second)*0.5);
//           MRIvalscale(mri_Twarp,mri_Twarp,(1.0+ 1.0/md.second)*0.5);
//        }

    // compute Registration
    if (verbose >1 ) cout << "- compute new iterative registration" << endl;
    if (cmd.first) MatrixFree(&cmd.first);
    cmd = computeIterativeRegistration(n,epsit,gpS[r],gpT[r],md.first,md.second);
//       cmd = computeIterativeRegSat(n,gpS[r],gpT[r],md.first,md.second);
    if (verbose > 1)
    {
      cout << endl << " current : Matrix: " << endl;
      MatrixPrintFmt(stdout,"% 2.8f",cmd.first);
      cout << " intens: " << cmd.second << endl;

      // adjust to highest level for output only:
      double tx = cmd.first->rptr[1][4];
      double ty = cmd.first->rptr[2][4];
      double tz = cmd.first->rptr[3][4];
      for (int ll = r; ll > 0; ll--)
      {
        tx *= 2;
        ty*=2;
        tz*=2;
      }
      cout << " equiv trans on highres: " << tx << " " << ty << " " << tz << endl;

      // save resampled version on this level:
//    MRI* salign = MRIlinearTransform(gpS[r], NULL,cmd.first);
//    MRIwrite(salign,"mriS-lowres-aligned.mgz");
//     MRIwrite(gpT[r],"mriT-lowres.mgz");
//     MRIfree(&salign);
    }

    if (r !=0) // adjust matrix to higher resolution level
    {
      for (int rr = 1; rr<=3; rr++)
      {
        // md.first->rptr[rr][4]  = 2.0 *  md.first->rptr[rr][4];
        cmd.first->rptr[rr][4] = 2.0 * cmd.first->rptr[rr][4];
      }
    }
//       m  = MatrixMultiply(cmd.first,md.first,m);
//       MatrixCopy(m,md.first);
    MatrixCopy(cmd.first,md.first);
    md.second = cmd.second;
    if (verbose > 1)
    {
      cout << endl << " Matrix: " << endl;
      MatrixPrintFmt(stdout,"% 2.8f",md.first);
      cout << " intens: " << md.second << endl;
    }
  } // resolution loop
	
    if (verbose == 1)
    {
      cout << endl << "   - final transform: " << endl;
      MatrixPrintFmt(stdout,"% 2.8f",md.first);
      cout << "   - final iscale: " << md.second << endl;
    }
	
	
  // cleanup
  if (cmd.first) MatrixFree(&cmd.first);
  MatrixFree(&m);
  //MRIfree(&mri_Swarp);
  //MRIfree(&mri_Twarp);

  Mfinal = MatrixCopy(md.first, Mfinal);
  iscalefinal = md.second;

  return md;



}

double  Registration::computeSatEstimate (int reslevel, int n,double epsit, MRI * mriS, MRI* mriT, MATRIX* mi, double scaleinit )
{
  cout << " Registration::computeSatEstimate " << endl;
  reslevel = 1;
  double PERCENT = 0.85;
  double EPS     = 0.01;

  pair < MATRIX*, double> fmd = computeMultiresRegistration(reslevel+1,n,epsit,mriS,mriT,mi,scaleinit);
  //     cout << endl << " Matrix: " << endl;
  //     MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
  //     cout << " intens: " << fmd.second << endl;

  // Md is allready adjusted to current reslevel

  cout <<  endl << "Compute Sat estimate on Resolution " << reslevel << endl;
  cout << "   - warping source and target (sqrt)" << endl;

  MATRIX * mh = MyMatrix::MatrixSqrt(fmd.first);
  // do not just assume m = mh*mh, rather m = mh2 * mh
  // for transforming target we need mh2^-1 = mh * m^-1
  MATRIX * mii  = MatrixInverse(fmd.first,NULL);
  MATRIX *mhi = MatrixMultiply(mh,mii,NULL);


  MRI* mri_Swarp = MRIclone(gpS[reslevel],NULL);
  mri_Swarp = MRIlinearTransform(gpS[reslevel],mri_Swarp, mh);
  MRI* mri_Twarp = MRIclone(gpS[reslevel],NULL); // bring them to same space (just use src geometry)
  mri_Twarp = MRIlinearTransform(gpT[reslevel],mri_Twarp, mhi);

  // adjust intensity
  if (iscale)
  {
    cout << "   - adjusting intensity ( "<< fmd.second << " ) " << endl;
    //MRIvalscale(mri_Swarp,mri_Swarp,fmd.second);
    MyMRI::MRIvalscale(mri_Swarp,mri_Swarp,(1.0+fmd.second)*0.5);
    MyMRI::MRIvalscale(mri_Twarp,mri_Twarp,(1.0+ 1.0/fmd.second)*0.5);
  }

  pair < MATRIX*, VECTOR*> Ab;
  Ab = constructAb(mri_Swarp,mri_Twarp);
  pair < MATRIX*, MATRIX* > pwm(NULL,NULL);
  Regression R(Ab.first,Ab.second);
	R.setVerbose(verbose);

  pair < double , double > interval(3,20);
  pair < double , double > wpercent;
  cout << "   - get weights left ( " << interval.first << " )" << endl;
  pwm = R.getRobustEstW(interval.first);
  cout << endl;
  MatrixFree(&pwm.first);
  MatrixFree(&pwm.second);
  wpercent.first = R.getLastWeightPercent();
  cout << "   - get weights right ( " << interval.second << " )" << endl;
  pwm = R.getRobustEstW(interval.second);
  cout << endl;
  MatrixFree(&pwm.first);
  MatrixFree(&pwm.second);
  wpercent.second = R.getLastWeightPercent();

  assert (wpercent.first < PERCENT); // should be at sat ==1, otherwise one could try to go smaller?
  assert (wpercent.second > PERCENT); // should be at sat ==20, otherwise one could try to go higher
  double m, mp;
  int count = 0;
  while (interval.second - interval.first > 0.1 && wpercent.second - wpercent.first > EPS)
  {
    cout << endl << " Interval : w[ " << interval.first << " , " << interval.second << " ] = [ " << wpercent.first << " , " << wpercent.second << " ] : " << wpercent.second - wpercent.first<< endl;
    count++;

    // m = (PERCENT - wpercent.first) * (interval.second - interval.first) /  (wpercent.second - wpercent.first) + interval.first;
    m = 0.5*(interval.second + interval.first);
    cout << "   new test: " << m  << endl;
    pwm = R.getRobustEstW(m);
    cout << endl;
    MatrixFree(&pwm.first);
    MatrixFree(&pwm.second);
    mp = R.getLastWeightPercent();
    cout << "    yields: " << mp << endl;

    if (mp < PERCENT)
    {
      interval.first = m;
      wpercent.first = mp;
    }
    else
    {
      interval.second = m;
      wpercent.second = mp;
    }
  }

  // cleanup
  if (mri_Twarp) MRIfree(&mri_Twarp);
  if (mri_Swarp) MRIfree(&mri_Swarp);
  MatrixFree(&mh);
  MatrixFree(&mii);
  MatrixFree(&mhi);

  cout << "Optimal sat ( at " << PERCENT << " ) : " << (interval.second + interval.first) *0.5 << endl;
  cout << " after " << count << " steps " << endl;
  return (interval.second + interval.first) *0.5;
}

// double  Registration::computeSatEstimate (int reslevel, int n,double epsit, MRI * mriS, MRI* mriT, MATRIX* mi, double scaleinit )
// {
//    cout << " Registration::computeSatEstimate " << endl;
//    if (!mriS) mriS = mri_source;
//    if (!mriT) mriT = mri_target;
//
//    cout << "   - building Gaussian Pyramid " << endl;
//    vector < MRI* > gpS = buildGaussianPyramid(mriS,100);
//    vector < MRI* > gpT = buildGaussianPyramid(mriT,100);
//    assert(gpS.size() == gpT.size());
//    int resolution = gpS.size();
//
//
//    int stoplevel = reslevel;
//
// // debug : save pyramid
// //  for (uint i = 0;i<gpS.size();i++)
// //  {
// //   char fn[40];
// //   sprintf(fn, "pyramid-%d.mgz", i+1);
// //   MRIwrite(gpS[i],fn);
// //  }
//
//
//    vector < pair < double , double >  >  satzero(20);
//    for (int i = 0 ; i<20;i++) satzero[i].first = double(i)+1.0;
//
//
//
// for (int s = 0 ; s < 20; s++)
// {
//
//    cout << endl << " SATUTRATION : " << satzero[s].first << endl;
//    sat = satzero[s].first;
//
//   MATRIX *m  = MatrixIdentity(4,NULL);
//
//   // variables to store matrix m and scaling factor d:
//   pair < MATRIX* , double > cmd;
//   pair < MATRIX* , double > md(NULL,scaleinit);
//
//   // check if mi (inital transform) is passed
//   if (mi)
//   {
//       md.first = MatrixCopy(mi,NULL);
//   }
//   else if (Minit)
//       md.first = MatrixCopy(Minit,NULL);
//   else
//   {
//     //  md.first = MatrixIdentity(4,NULL);
//     //   md.first = initialize_transform(mriS,mriT);
//     // use voxtovox as init:
//     md.first = MRIgetVoxelToVoxelXform(mriS,mriT) ;
//   }
//
//   if (debug > 0)
//   {
//      cout << "   - initial transform:\n" ;
//      MatrixPrintFmt(stdout,"% 2.8f",md.first);
//   }
//
//   // adjust minit to current (lowest) resolution:
//   int rstart = 2;
//   for (int r = 1; r<=resolution-rstart; r++)
//   for (int rr = 1;rr<=3;rr++)
//      md.first->rptr[rr][4]  = 0.5 *  md.first->rptr[rr][4];
//
//    if(debug >0)
//    {
//       cout << "   - initial adjusted:\n" ;
//       MatrixPrintFmt(stdout,"% 2.8f",md.first);
//    }
//
//
// //     MRI* falign = MRIlinearTransform(gpS[resolution-rstart], NULL,md.first);
// //     MRIwrite(falign,"mriS-lowres-initaligned.mgz");
// //     MRIfree(&falign);
// //     MRIwrite(gpT[resolution-rstart],"mriT-lowres.mgz");
//
//
// //  md.second = 1; //??
// //  MRI* mri_Swarp   = NULL;
// //  MRI* mri_Twarp   = NULL;
//    vector < MATRIX* > matvec (resolution-rstart+1,NULL);
//    vector < double >  intvec (resolution-rstart+1,1);
// //  for (int r = resolution-rstart;r>=0;r--)
//   for (int r = resolution-rstart;r>=stoplevel;r--)
//   {
//  //    MRIwrite(gpS[r],"mriS-smooth.mgz");
//  //    MRIwrite(gpT[r],"mriT-smooth.mgz");
//
//       cout << endl << "Resolution: " << r << endl;
//
//
//        // compute Registration
//        cout << "   - compute new registration" << endl;
//        if (cmd.first) MatrixFree(&cmd.first);
//        cmd = computeIterativeRegistration(n,epsit,gpS[r],gpT[r],md.first,md.second);
// //       cmd = computeIterativeRegSat(n,gpS[r],gpT[r],md.first,md.second);
//        matvec[r] = MatrixCopy(cmd.first,matvec[r]);
//        intvec[r] = cmd.second;
//
//     if(debug > 0)
//     {
//      cout << endl << " current : Matrix: " << endl;
//      MatrixPrintFmt(stdout,"% 2.8f",cmd.first);
//      cout << " intens: " << cmd.second << endl;
//
//      // adjust to highest level for output only:
//      double tx = cmd.first->rptr[1][4];
//      double ty = cmd.first->rptr[2][4];
//      double tz = cmd.first->rptr[3][4];
//      for (int ll = r; ll > 0; ll--)
//      { tx *= 2; ty*=2; tz*=2;}
//       cout << " equiv trans on highres: " << tx << " " << ty << " " << tz << endl;
//
//      // save resampled version on this level:
//  //    MRI* salign = MRIlinearTransform(gpS[r], NULL,cmd.first);
//  //    MRIwrite(salign,"mriS-lowres-aligned.mgz");
// //     MRIwrite(gpT[r],"mriT-lowres.mgz");
// //     MRIfree(&salign);
//      }
//
//       if (r !=0) // adjust matrix to higher resolution level
//       {
//          for (int rr = 1; rr<=3; rr++)
//   {
//            // md.first->rptr[rr][4]  = 2.0 *  md.first->rptr[rr][4];
//             cmd.first->rptr[rr][4] = 2.0 * cmd.first->rptr[rr][4];
//   }
//       }
// //       m  = MatrixMultiply(cmd.first,md.first,m);
// //       MatrixCopy(m,md.first);
//        MatrixCopy(cmd.first,md.first);
//        md.second = cmd.second;
//      if (debug > 0)
//      {
//         cout << endl << " Matrix: " << endl;
//         MatrixPrintFmt(stdout,"% 2.8f",md.first);
//         cout << " intens: " << md.second << endl;
//      }
//   }
//   // cleanup
//   if (cmd.first) MatrixFree(&cmd.first);
//   MatrixFree(&m);
//   if (md.first) MatrixFree(&md.first);
//
//   //Mfinal = MatrixCopy(md.first, Mfinal);
//   //iscalefinal = md.second;
//
//   satzero[s].second = zeroweights;
//   cout << " sat: " << satzero[s].first << "  zeroweights: " << satzero[s].second << endl;
//
// }
//    // plot diffs
//    string fbase = name;
//    int rf = fbase.rfind("/");
//  if (rf != -1)
//  {
//             fbase = fbase.substr(rf+1,fbase.length());
//  }
//
//    ofstream ofile((name+".plot").c_str(),ios::out);
//    bool png = false;
//    if (png) ofile << "set terminal png medium size 800,600" << endl;
//    else ofile << "set terminal postscript eps color" << endl;
//    if (png) ofile << "set output \""<< fbase <<".png\"" << endl;
//    else ofile << "set output \""<< fbase <<".eps\"" << endl;
//    ofile << "plot ";
//    ofile << " \"-\" notitle with lines 1" << endl;
//    for (int j = 0;j<(int)satzero.size(); j++)
//    {
//       ofile << satzero[j].first << " " << satzero[j].second << endl;
//    }
//    ofile << "e" << endl;
//    ofile.close();
//
//   freeGaussianPyramid(gpS);
//   freeGaussianPyramid(gpT);
//
//   return -1.0;
// }
//
void Registration::testRobust(const std::string& fname, int testno)
{


  cout << " testRobust " << fname << "  testno: " << testno << endl;

  char fn[fname.length()+1];
  for (uint i = 0;i<fname.length();i++)
    fn[i] = fname[i];
  fn[fname.length()] = '\0';

//  MRI * mri1 = MRIalloc(10, 10, 10, MRI_FLOAT);
//  MRI * mri2 = MRIalloc(10, 10, 10, MRI_FLOAT);
//   MRI * mri1 = MRIalloc(10, 10, 10, MRI_INT);
//   MRI * mri2 = MRIalloc(10, 10, 10, MRI_INT);
//      int x,y,z;
//      for (z = 0 ; z < mri1->depth  ; z++)
//      for (y = 0 ; y < mri1->height ; y++)
//      for (x = 0 ; x < mri1->width  ; x++)
//      {
//      //cout << " x: " << x << "  y: " << y << " z: " << z << endl;
//       //  MRIFvox(mri, x, y, z) = *MATRIX_RELT(ii, y+1, x+1);
//       //  MRIFvox(mri1, x, y, z) = (x+1)*(y+1);
//       //  MRIFvox(mri2, x, y, z) = (x+1)*(y+1) +1;
//         MRIvox(mri1, x, y, z) = (x+1)*(y+1);
//         MRIvox(mri2, x, y, z) = (x+1)*(y+1) +1;
//      }
//   cout << " compute Registration now:"<< endl;
//   pair <MATRIX *,MRI*> pwt = computeRegistration(mri1,mri2,true,true,false);
//   cout << " done" << endl;
//   MatrixPrintFmt(stdout,"% 2.8f",pwt.first);
//   exit(1);
//      MRIwrite(mri1,"test1.mgz");
//      MRIwrite(mri2,"test2.mgz");



  MRI * mri = MRIread(fn) ;
  cout << " Read MRI successfully!" << endl;

//    MATRIX*  testm  = MatrixAllocRotation(4,0.08,Z_ROTATION);
//     testm->rptr[1][4] = 2.2;
//     testm->rptr[2][4] = .3;
//     testm->rptr[3][4] = 1.4;
//
// cout << "Testmatrix: " << endl;
//    MatrixPrintFmt(stdout,"% 2.8f",testm);
// cout << endl;
// cout << "Half RT: " << endl;
//    MATRIX* th1 = getHalfRT(testm);
//    MatrixPrintFmt(stdout,"% 2.8f",th1);
// cout << endl;
// cout << " half*half: " << endl;
//    MATRIX* thm1 = MatrixMultiply(th1, th1,NULL);
//    MatrixPrintFmt(stdout,"% 2.8f",thm1);
// cout << endl;
//
//
// cout << "sqrt  : " << endl;
//    MATRIX* th2 = MatrixSqrt(testm);
//    MatrixPrintFmt(stdout,"% 2.8f",th2);
// cout << endl;
// cout << " half*half: " << endl;
//    MATRIX* thm2 = MatrixMultiply(th2, th2,NULL);
//    MatrixPrintFmt(stdout,"% 2.8f",thm2);
// cout << endl;
//
// exit(1);

  vector < MRI* > gpS = buildGaussianPyramid(mri,100);
  int level = gpS.size();
//  int level = 4;
//  MRIwrite(gpS[gpS.size()-level],"small.mgz");
  //cout << "sfasf" << endl;

  MATRIX* a = NULL, *ai = NULL;
  MRI* mriTs = NULL, *mriTt = NULL;

  double theta  ;
  double iscaleval = 1.0;
  switch (testno)
  {
  case 0 : // identity
    cout << "Test " << testno << " : Identity" << endl;
    a  = MatrixIdentity(4,a);
    ai = MatrixIdentity(4,ai);
    mriTs = MRIcopy(gpS[gpS.size()-level],NULL);
    mriTt = MRIcopy(gpS[gpS.size()-level],NULL);
    MRIwrite(mriTs,"idS.mgz");
    MRIwrite(mriTt,"idT.mgz");

    break;
  case 1 : // translation
    cout << "Test " << testno << " : Translation" << endl;
    a  = MatrixIdentity(4,a);
    a->rptr[1][4] = 2.2;
    a->rptr[2][4] = .3;
    a->rptr[3][4] = 1.4;
    ai = MatrixInverse(a,ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);
    MRIwrite(mriTs,"transS.mgz");
    MRIwrite(mriTt,"transT.mgz");

    break;
  case 2:   // rotation
    cout << "Test " << testno << " : Rotation" << endl;
    theta = 0.08;
    a  = MatrixAllocRotation(4,theta,Z_ROTATION);
    ai = MatrixInverse(a,ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);
    MRIwrite(mriTs,"rotS.mgz");
    MRIwrite(mriTt,"rotT.mgz");
    break;
  case 3:   // intensity
    cout << "Test " << testno << " : Intensity" << endl;
    a  = MatrixIdentity(4,a);
    ai = MatrixIdentity(4,ai);
    iscaleval = 0.8;
    mriTs = MRIcopy(gpS[gpS.size()-level], NULL);
    mriTt = MyMRI::MRIvalscale(gpS[gpS.size()-level], NULL, iscaleval);
    MRIwrite(mriTs,"iscaleS.mgz");
    MRIwrite(mriTt,"iscaleT.mgz");
    break;
  case 4:   // rotation and translation
    cout << "Test " << testno << " : Rotation and Translation" << endl;
    theta = 0.08;
    a  = MatrixAllocRotation(4,theta,Z_ROTATION);
    a->rptr[1][4] = 2.2;
    a->rptr[2][4] = .3;
    a->rptr[3][4] = 1.4;
    ai = MatrixInverse(a,ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);
    iscaleval = 0.8;
    mriTt = MyMRI::MRIvalscale(mriTt, NULL, iscaleval);
    MRIwrite(mriTs,"rottransS.mgz");
    MRIwrite(mriTt,"rottransT.mgz");
    break;
  case 5 : // translation and junk
    cout << "Test " << testno << " : Translation and Noise" << endl;
    a  = MatrixIdentity(4,a);
    a->rptr[1][4] = .2;
    a->rptr[2][4] = .3;
    a->rptr[3][4] = .4;
    ai = MatrixInverse(a,ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);

    for (int dd = 0;dd<mriTs->depth/3;dd++)
      for (int cc = 0;cc<mriTs->height/3;cc++)
        for (int rr = 0;rr<mriTs->width/3;rr++)
          MRIvox(mriTs, rr, cc, dd) = (rr)+(cc)+(dd)+1;

    MRIwrite(mriTs,"junktransS.mgz");
    MRIwrite(mriTt,"junktransT.mgz");
    break;
  case 6 : // rotation and junk
    cout << "Test " << testno << " : Rotation and Noise" << endl;
    theta = 0.02;
    a  = MatrixAllocRotation(4,theta,Z_ROTATION);
    ai = MatrixInverse(a,ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);

    //double nn = (mriTs->depth/3) +(mriTs->height/3) +(mriTs->width/3)+1;
    for (int dd = 0;dd<mriTs->depth/3;dd++)
      for (int cc = 0;cc<mriTs->height/3;cc++)
        for (int rr = 0;rr<mriTs->width/3;rr++)
          MRIvox(mriTs, rr, cc, dd) = ((rr)+(cc)+(dd)+1) ;

    MRIwrite(mriTs,"junkrotS.mgz");
    MRIwrite(mriTt,"junkrotT.mgz");

    break;
  case 7 : // skaling
    cout << "Test " << testno << " : Scaling" << endl;
    a  = MatrixIdentity(4,NULL);
    a->rptr[1][1] = 1.01;
    a->rptr[2][2] = 1.04;
    a->rptr[3][3] = 1.06;

    ai = MatrixInverse(a,ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);

    //double nn = (mriTs->depth/3) +(mriTs->height/3) +(mriTs->width/3)+1;
//    for (int dd = 0;dd<mriTs->depth/3;dd++)
//   for (int cc = 0;cc<mriTs->height/3;cc++)
//   for (int rr = 0;rr<mriTs->width/3;rr++)
//     MRIvox(mriTs, rr, cc, dd) = ((rr)+(cc)+(dd)+1) ;

    MRIwrite(mriTs,"scaleS.mgz");
    MRIwrite(mriTt,"scaleT.mgz");

    break;

  case 20: //error functions when rotating
  {
    int steps = 50;
    double div = 4.0;
    vector < double > theta(steps);
    vector < double > err(steps);
    vector < double > mls(steps);
    vector < double > mls2(steps);
    //level--;
    for (int i=0; i<steps; i++)
    {
      // 0.. PI/div in 20 steps
      // -PI/div ..0 is symmetric
      theta[i] = M_PI * i / ((steps-1)*div);

      a  = MatrixAllocRotation(4,0.5*theta[i],Z_ROTATION);
      a  = MatrixMultiply(a, extract_i_to_r(gpS[gpS.size()-level]), a);
      a  = MatrixMultiply(extract_r_to_i(gpS[gpS.size()-level]) , a, a);
      ai = MatrixInverse(a,ai);

      mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
      mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);
      //MRIwrite(mriTs,"test20-s.mgz");
      //MRIwrite(mriTt,"test20-t.mgz");
      MatrixFree(&a);
      MatrixFree(&ai);
      ai = NULL;
      pair < MATRIX*, VECTOR*> Ab;
      transonly = false;
      robust    = true;
      rigid     = true;
      iscale    = false;

      Ab = constructAb(mriTs, mriTt);
      pair < MATRIX*, MATRIX* > pwm(NULL,NULL);
      Regression R(Ab.first,Ab.second);
			R.setVerbose(verbose);
      sat = 5;

      cout << "   - compute robust estimate ( sat "<<sat<<" )..." << flush;
      pwm = R.getRobustEstW(sat);
      MatrixFree(&pwm.first);
      MatrixFree(&pwm.second);
      err[i] = R.getLastError();
      cout << "angle: " << theta[i] << "  error: " << err[i] << endl;
      MATRIX *m = R.getLSEst();
      MatrixFree(&m);
      mls[i] = R.getLastError();
      cout << "angle: " << theta[i] << "  mls: " << mls[i] << endl;
      MRI * mridiff = MRIalloc(mriTs->width, mriTs->height, mriTs->depth, MRI_FLOAT);
      mridiff = MRIsubtract(mriTs,mriTt,mridiff);
      double ddd = 0;
      for (int d = 0;d<mriTs->depth;d++)
        for (int h = 0;h<mriTs->height;h++)
          for (int w = 0;w<mriTs->width;w++)
            ddd += MRIgetVoxVal(mridiff,w,h,d,1) * MRIgetVoxVal(mridiff,w,h,d,1);
      mls2[i] = ddd;
      cout << "angle: " << theta[i] << "  mls: " << mls2[i] << endl;
      MatrixFree(&Ab.first);
      MatrixFree(&Ab.second);
    }

    ostringstream ss;
    ss << "r-error-rot4-l" << level;
    string fn = ss.str()+".plot";
    ofstream f(fn.c_str(),ios::out);

    f << "set terminal postscript eps color" << endl;
    f << "set title \"(Robust) error when rotating on level " << level <<"\"" << endl;
    f << "set output \""<< ss.str() << ".eps\"" << endl;
    f << "plot  \"-\" notitle with lines 1" << endl;
    for (int i=0; i<steps; i++)
    {
      cout << theta[i] << " " << err[i] << endl;
      f << theta[i] << " " << err[i] << endl;
    }
    f << "e" << endl;

    ostringstream ss2;
    ss2 << "ls-error-rot4-l" << level;
    string fn2 = ss2.str()+".plot";
    ofstream f2(fn2.c_str(),ios::out);

    f2 << "set terminal postscript eps color" << endl;
    f2 << "set title \"(LeastSquares) error when rotating on level " << level <<"\"" << endl;
    f2 << "set output \""<< ss2.str() << ".eps\"" << endl;
    f2 << "plot  \"-\" notitle with lines 1" << endl;
    for (int i=0; i<steps; i++)
    {
      cout << theta[i] << " " << mls[i] << endl;
      f2 << theta[i] << " " << mls[i] << endl;
    }
    f2 << "e" << endl;

    ostringstream ss3;
    ss3 << "ils-error-rot4-l" << level;
    string fn3 = ss3.str()+".plot";
    ofstream f3(fn3.c_str(),ios::out);

    f3 << "set terminal postscript eps color" << endl;
    f3 << "set title \"(IntensityLeastSquares) error when rotating on level " << level <<"\"" << endl;
    f3 << "set output \""<< ss3.str() << ".eps\"" << endl;
    f3 << "plot  \"-\" notitle with lines 1" << endl;
    for (int i=0; i<steps; i++)
    {
      cout << theta[i] << " " << mls2[i] << endl;
      f3 << theta[i] << " " << mls2[i] << endl;
    }
    f3 << "e" << endl;

    exit(0);
    break;
  }
  default:
    assert(1==2);
  }

  cout << " Transformed , now registering ..." << endl;

  int steps;
  steps = 3;
  rtype = 2;

//    transonly = true;
//    robust = false;
//    rigid = false;
//    iscale = false;
//    pair <MATRIX*, double> pwlst  = computeIterativeRegistration(steps,mriTs,mriTt);
//    robust = true;
//    pair <MATRIX*, double> pwt  = computeIterativeRegistration(steps,mriTs,mriTt);
//    iscale = true;
//    pair <MATRIX*, double> pwit   = computeIterativeRegistration(steps,mriTs,mriTt);

  transonly = false;
  rigid     = false;
  robust = true;
  iscale = true;
  sat = 5;

//  pair <MATRIX*, double> pwit    = computeIterativeRegistration(steps,mriTs,mriTt);
  pair <MATRIX*, double> pw    = computeMultiresRegistration(0,5,0.01,mriTs,mriTt);
  exit(0);
  robust = false;
//   pair <MATRIX*, double> pwls  = computeIterativeRegistration(steps,mriTs,mriTt);
//   pair <MATRIX*, double> pwls  = computeMultiresRegistration(mriTs,mriTt);
  iscale = true;
  robust = true;
//   pair <MATRIX*, double> pwi    = computeIterativeRegistration(steps,mriTs,mriTt);
//    pair <MATRIX*, double> pwi   = computeMultiresRegistration(mriTs,mriTt);
//
//    robust = false;
//    rigid = false;
//    iscale = false;
//    pair <MATRIX*, double> apwls = computeIterativeRegistration(steps,mriTs,mriTt);
//    robust = true;
//    pair <MATRIX*, double> apw   = computeIterativeRegistration(steps,mriTs,mriTt);
//    iscale = true;
//    pair <MATRIX*, double> apwi  = computeIterativeRegistration(steps,mriTs,mriTt);
//

  cout << endl << endl << " Actual Transformation: " << endl;
  MATRIX * aa = MatrixMultiply(ai,ai,NULL);
  MatrixPrintFmt(stdout,"% 2.8f",aa);
  cout << " iscale: " << iscaleval << endl;
  MatrixFree(&aa);

//   cout << endl << " Trans - Least Square: " << endl;
//   MatrixPrintFmt(stdout,"% 2.8f",pwlst.first);
//   MatrixFree(&pwlst.first);
//   cout << endl << " Trans - Robust M-Est: " << endl;
//   MatrixPrintFmt(stdout,"% 2.8f",pwt.first);
//   MatrixFree(&pwt.first);
//   cout << endl << " Trans - Robust M-Est - iscale: " << endl;
//   MatrixPrintFmt(stdout,"% 2.8f",pwit.first);
//   cout << " iscale: " << pwi.second << endl << endl;
//   MatrixFree(&pwit.first);

//  cout << endl << " Rigid - Least Square: " << endl;
//  MatrixPrintFmt(stdout,"% 2.8f",pwls.first);
//  MatrixFree(&pwls.first);

  cout << endl << " Rigid - Robust M-Est: " << endl;
  MatrixPrintFmt(stdout,"% 2.8f",pw.first);
  MatrixFree(&pw.first);

//  cout << endl << " Rigid - Robust M-Est (iterations only): " << endl;
//  MatrixPrintFmt(stdout,"% 2.8f",pwit.first);
//  MatrixFree(&pwit.first);

//   cout << endl << " Rigid - Robust M-Est - iscale: " << endl;
//   MatrixPrintFmt(stdout,"% 2.8f",pwi.first);
//   cout << " iscale: " << pwi.second << endl << endl;
//   MatrixFree(&pwi.first);
//
//   cout << endl << " Non Rigid - Least Square: " << endl;
//   MatrixPrintFmt(stdout,"% 2.8f",apwls.first);
//   MatrixFree(&apwls.first);
//
//   cout << endl << " Non Rigid - Robust M-Est: " << endl;
//   MatrixPrintFmt(stdout,"% 2.8f",apw.first);
//   MatrixFree(&apw.first);
//
//   cout << endl << " Non Rigid - Robust M-Est - iscale: " << endl;
//   MatrixPrintFmt(stdout,"% 2.8f",apwi.first);
//   cout << " iscale: " << apwi.second << endl << endl;
//   MatrixFree(&apwi.first);


  exit(0);
}

bool Registration::warpSource( const string &fname, MATRIX* M, double is)
{
  assert(mri_source);
  assert(mri_target);

  int nframes = mri_source->nframes;
  mri_source->nframes = 1 ;

  MRI *mri_aligned = MRIclone(mri_target,NULL);
  if (M)
  {
    // cout << "using M:" << endl ;
    // MatrixPrint(stdout,M) ;
    mri_aligned = MRIlinearTransform(mri_source, mri_aligned,M);
  }
  else if (Mfinal)  mri_aligned = MRIlinearTransform(mri_source, mri_aligned,Mfinal);
  else
  {
    cerr << "warpSource error: no matrix set!" << endl;
    MRIfree(&mri_aligned) ;
    return false;
  }

  // here also do scaling of intensity values
  if (is == -1) is = iscalefinal;
  if (is >0 && is != 1)
    mri_aligned = MyMRI::MRIvalscale(mri_aligned, mri_aligned, is);

  mri_source->nframes = nframes ;

//    sprintf(fname, "%s_after_final_alignment", parms.base_name) ;
//    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
//    sprintf(fname, "%s_target", parms.base_name) ;
//    MRIwriteImageViews(mri_dst, fname, IMAGE_SIZE) ;

  char cfname[STRLEN];
  strcpy(cfname, fname.c_str()) ;

  MRIwrite(mri_aligned, cfname) ;
  MRIfree(&mri_aligned) ;

  return true; // no error treatment so far
}

bool Registration::warpSource(MRI* orig, MRI* target, const string &fname, MATRIX* M, double is)
// warps the mri orig to target
{
  assert(orig);

  int nframes = orig->nframes;
  orig->nframes = 1 ;

  MATRIX* m_Lvox;
  if (M) m_Lvox = MatrixCopy(M,NULL);
  else m_Lvox   = MatrixCopy(Mfinal,NULL);

  /* convert it to RAS mm coordinates */
  MATRIX* m_L = MRIvoxelXformToRasXform(mri_source, mri_target, m_Lvox, NULL) ;

  MRI *mri_aligned = MRIclone(target,NULL); //cp header and alloc space
  mri_aligned = MRIapplyRASlinearTransform(orig, mri_aligned, m_L) ;

  // here also do scaling of intensity values
  if (is == -1) is = iscalefinal;
  if (is >0 && is != 1)
    mri_aligned = MyMRI::MRIvalscale(mri_aligned, mri_aligned, is);

  orig->nframes = nframes ;

//    sprintf(fname, "%s_after_final_alignment", parms.base_name) ;
//    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
//    sprintf(fname, "%s_target", parms.base_name) ;
//    MRIwriteImageViews(mri_dst, fname, IMAGE_SIZE) ;

  char cfname[STRLEN];
  strcpy(cfname, fname.c_str()) ;

  MRIwrite(mri_aligned, cfname) ;
  MRIfree(&mri_aligned) ;

  return true; // no error treatment so far
}

pair < MATRIX*, VECTOR* > Registration::constructAb(MRI *mriS, MRI *mriT)
// similar to robust paper
// (with symmetry and iscale)
{

  if (verbose > 1) cout << "   - constructAb: " << endl;

  assert(mriT != NULL);
  assert(mriS != NULL);
  assert(mriS->width == mriT->width);
  assert(mriS->height== mriT->height);
  assert(mriS->depth == mriT->depth);
  assert(mriS->type  == mriT->type);
  //assert(mriS->width == mask->width);
  //assert(mriS->height== mask->height);
  //assert(mriS->depth == mask->depth);
  //assert(mask->type == MRI_INT);
  //MRIclear(mask);

  int z,y,x;
  long int ss = mriS->width * mriS->height * mriS->depth;
  if (mri_indexing) MRIfree(&mri_indexing);
  if (ss > std::numeric_limits<int>::max())
  {
     if (verbose > 1) cout << "     -- using LONG for indexing ... " << flush;
     mri_indexing = MRIalloc(mriS->width, mriS->height, mriS->depth,MRI_LONG);
     if (mri_indexing == NULL) 
        ErrorExit(ERROR_NO_MEMORY,"Registration::constructAB could not allocate memory for mri_indexing") ;
     if (verbose > 1) cout << " done!" << endl;
  }
  else 
  {
     double mu = ((double)ss) * sizeof(int) / (1024.0 * 1024.0);
     if (verbose > 1) cout << "     -- allocating " << mu << "Mb mem for indexing ... " << flush;
     mri_indexing = MRIalloc(mriS->width, mriS->height, mriS->depth,MRI_INT);
     if (mri_indexing == NULL) 
        ErrorExit(ERROR_NO_MEMORY,"Registration::constructAB could not allocate memory for mri_indexing") ;
     if (verbose > 1) cout << " done!" << endl;
  }

  for (z = 0 ; z < mriS->depth ; z++)
    for (x = 0 ; x < mriS->width ; x++)
      for (y = 0 ; y < mriS->height ; y++)
        MRILvox(mri_indexing, x, y, z) = 0;


  bool dosubsample = false;
  if (subsamplesize > 0)
    dosubsample = (mriS->width > subsamplesize && mriS->height > subsamplesize && mriS->depth > subsamplesize);
//   bool dosubsample = true;

//    // remove 0 voxels:
//   for (z = 0 ; z < mriS->depth ; z++)
//   for (y = 0 ; y < mriS->height ; y++)
//   for (x = 0 ; x < mriS->width ; x++)
//   {
//      if (MRIgetVoxVal(mriS,x,y,z,0) <= 0)
//            MRIsetVoxVal(mriS,x,y,z,0,std::numeric_limits<float>::quiet_NaN());
//      if (MRIgetVoxVal(mriT,x,y,z,0) <= 0)
//            MRIsetVoxVal(mriT,x,y,z,0,std::numeric_limits<float>::quiet_NaN());
//
//   }

  // we will need the derivatives
  if (verbose > 1) cout << "     -- compute derivatives ... " << flush;
  MRI *Sfx=NULL,*Sfy=NULL,*Sfz=NULL,*Sbl=NULL;
  MRI *Tfx=NULL,*Tfy=NULL,*Tfz=NULL,*Tbl=NULL;
  MyMRI::getPartials(mriS,Sfx,Sfy,Sfz,Sbl);
  MyMRI::getPartials(mriT,Tfx,Tfy,Tfz,Tbl);

  MRI * fx1  = MRIadd(Sfx,Tfx,Tfx);
  MRIscalarMul(fx1,fx1,0.5);
  MRI * fy1  = MRIadd(Sfy,Tfy,Tfy);
  MRIscalarMul(fy1,fy1,0.5);
  MRI * fz1  = MRIadd(Sfz,Tfz,Tfz);
  MRIscalarMul(fz1,fz1,0.5);
 // MRI * ft1  = MRIsubtract(Tbl,Sbl,Tbl); //T-S = f1-f2 = - delta f from paper
  MRI * ft1  = MRIsubtract(Tbl,Sbl,NULL); //T-S = f1-f2 = - delta f from paper
  MRIfree(&Sfx);
  MRIfree(&Sfy);
  MRIfree(&Sfz);
  // Sbl will be needed below

  if (verbose > 1) cout << " done!" << endl;
  //MRIwrite(fx1,"fx.mgz");
  //MRIwrite(fy1,"fy.mgz");
  //MRIwrite(fz1,"fz.mgz");
  //MRIwrite(ft1,"ft.mgz");

  MRI * fx,* fy,* fz,* ft;
  if (dosubsample)
  {
    if (verbose > 1) cout << "     -- subsample ... "<< flush;

    fx = MyMRI::subSample(fx1);
    fy = MyMRI::subSample(fy1);
    fz = MyMRI::subSample(fz1);
    ft = MyMRI::subSample(ft1);
    MRIfree(&fx1);
    MRIfree(&fy1);
    MRIfree(&fz1);
    MRIfree(&ft1);

    if (verbose > 1) cout << " done! " << endl;
  }
  else
  {
    fx = fx1;
    fy = fy1;
    fz = fz1;
    ft = ft1;
  }

  // allocate the space
  int n = fx->width * fx->height * fx->depth;
  int pnum = 12;
  if (transonly)  pnum = 3;
  else if (rigid) pnum = 6;
  if (iscale) pnum++;
  pair <MATRIX*, VECTOR* > Ab;

  MATRIX* A = MatrixAlloc(n,pnum,MATRIX_REAL);
  VECTOR* b = MatrixAlloc(n,1,MATRIX_REAL);
     if (A == NULL || b == NULL) 
        ErrorExit(ERROR_NO_MEMORY,"Registration::constructAB could not allocate memory for A or b") ;

  if (verbose > 1) cout << "     -- size " << fx->width << " " << fx->height << " " << fx->depth << flush;

  long int count = 0;
  int xp1,yp1,zp1;
  double eps = 0.00001;
  for (z = 0 ; z < fx->depth ; z++)
    for (x = 0 ; x < fx->width ; x++)
      for (y = 0 ; y < fx->height ; y++)
      {
        if (isnan(MRIFvox(fx, x, y, z)) ||isnan(MRIFvox(fy, x, y, z)) || isnan(MRIFvox(fz, x, y, z)) || isnan(MRIFvox(ft, x, y, z)) )
        {
          if (verbose > 0) cout << " found a nan value!!!" << endl;
          continue;
        }

        if (dosubsample)
        {
          xp1 = 2*x+3;
          yp1 = 2*y+3;
          zp1 = 2*z+3;
        }
        else
        {
          xp1 = x+3;
          yp1 = y+3;
          zp1 = z+3; // if not subsampled
        }
        assert(xp1 < mriS->width);
        assert(yp1 < mriS->height);
        assert(zp1 < mriS->depth);


        if (fabs(MRIFvox(fx, x, y, z)) < eps  && fabs(MRIFvox(fy, x, y, z)) < eps &&  fabs(MRIFvox(fz, x, y, z)) < eps )
        {
          //cout << " found a zero row!!!" << endl;
          MRILvox(mri_indexing, xp1, yp1, zp1) = -1;
          continue;
        }

        count++; // start with 1

        if (xp1 >= mriS->width || yp1 >= mriS->height || zp1 >= mriS->depth)
        {

          cerr << " outside !!! " << xp1 << " " << yp1 << " " << zp1 << endl;
          assert(1==2);
        }

        MRILvox(mri_indexing, xp1, yp1, zp1) = count;

        //cout << "x: " << x << " y: " << y << " z: " << z << " coutn: "<< count << endl;
        //cout << " " << count << " mrifx: " << MRIFvox(mri_fx, x, y, z) << " mrifx int: " << (int)MRIvox(mri_fx,x,y,z) <<endl;
        if (transonly)
        {
          *MATRIX_RELT(A, count, 1) = MRIFvox(fx, x, y, z);
          *MATRIX_RELT(A, count, 2) = MRIFvox(fy, x, y, z);
          *MATRIX_RELT(A, count, 3) = MRIFvox(fz, x, y, z);
         // if (iscale) *MATRIX_RELT(A, count, 4) = MRIFvox(Sbl, x, y, z);
          if (iscale) *MATRIX_RELT(A, count, 4) = MRIFvox(Tbl, x, y, z);
        }
        else if (rigid)
        {
          *MATRIX_RELT(A, count, 1) = MRIFvox(fx, x, y, z);
          *MATRIX_RELT(A, count, 2) = MRIFvox(fy, x, y, z);
          *MATRIX_RELT(A, count, 3) = MRIFvox(fz, x, y, z);
          *MATRIX_RELT(A, count, 4) = MRIFvox(fz, x, y, z)*yp1 - MRIFvox(fy, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 5) = MRIFvox(fx, x, y, z)*zp1 - MRIFvox(fz, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 6) = MRIFvox(fy, x, y, z)*xp1 - MRIFvox(fx, x, y, z)*yp1;
//          if (iscale) *MATRIX_RELT(A, count, 7) = MRIFvox(Sbl, x, y, z);
          if (iscale) *MATRIX_RELT(A, count, 7) = MRIFvox(Tbl, x, y, z);
					
        }
        else // affine
        {
          *MATRIX_RELT(A, count, 1)  = MRIFvox(fx, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 2)  = MRIFvox(fx, x, y, z)*yp1;
          *MATRIX_RELT(A, count, 3)  = MRIFvox(fx, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 4)  = MRIFvox(fx, x, y, z);
          *MATRIX_RELT(A, count, 5)  = MRIFvox(fy, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 6)  = MRIFvox(fy, x, y, z)*yp1;
          *MATRIX_RELT(A, count, 7)  = MRIFvox(fy, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 8)  = MRIFvox(fy, x, y, z);
          *MATRIX_RELT(A, count, 9)  = MRIFvox(fz, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 10) = MRIFvox(fz, x, y, z)*yp1;
          *MATRIX_RELT(A, count, 11) = MRIFvox(fz, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 12) = MRIFvox(fz, x, y, z);
         // if (iscale) *MATRIX_RELT(A, count, 13) = MRIFvox(Sbl, x, y, z);
          if (iscale) *MATRIX_RELT(A, count, 13) = MRIFvox(Tbl, x, y, z);
        }

        *MATRIX_RELT(b, count, 1) = - MRIFvox(ft, x, y, z);

      }
      
  // free remaining MRI    
  MRIfree(&fx);
  MRIfree(&fy);
  MRIfree(&fz);
  MRIfree(&ft);
  MRIfree(&Sbl);
  if (Tbl) MRIfree(&Tbl);

  // adjust sizes
  Ab.first  = MatrixAlloc(count,pnum,MATRIX_REAL);
  Ab.second = MatrixAlloc(count,1,MATRIX_REAL);
  for (int rr = 1; rr<= count; rr++)
  {
    *MATRIX_RELT(Ab.second, rr, 1) = *MATRIX_RELT(b, rr, 1);
    for (int cc = 1; cc <= pnum; cc++)
    {
      *MATRIX_RELT(Ab.first, rr, cc) = *MATRIX_RELT(A, rr, cc);
      assert (!isnan(*MATRIX_RELT(Ab.first, rr, cc)));
    }
    assert (!isnan(*MATRIX_RELT(Ab.second, rr, 1)));
  }

  if (verbose > 1) cout << " ( " << count << " non-zero voxels )"<< endl;

  MatrixFree(&A);
  MatrixFree(&b);

  return Ab;
}

pair < MATRIX*, VECTOR* > Registration::constructAb2(MRI *mriS, MRI *mriT)
{
  if (verbose > 1) cout << " constructAb2 " << endl;
  assert(mriT != NULL);
  assert(mriS != NULL);
  assert(mriS->width == mriT->width);
  assert(mriS->height== mriT->height);
  assert(mriS->depth == mriT->depth);
  assert(mriS->type  == mriT->type);
  //assert(mriS->width == mask->width);
  //assert(mriS->height== mask->height);
  //assert(mriS->depth == mask->depth);
  //assert(mask->type == MRI_INT);
  //MRIclear(mask);

  // set <= 0 to NaN
  int z,y,x;
  for (z = 0 ; z < mriS->depth ; z++)
    for (y = 0 ; y < mriS->height ; y++)
      for (x = 0 ; x < mriS->width ; x++)
      {
        if (MRIgetVoxVal(mriS,x,y,z,0) <= 0)
          MRIsetVoxVal(mriS,x,y,z,0,std::numeric_limits<float>::quiet_NaN());
        if (MRIgetVoxVal(mriT,x,y,z,0) <= 0)
          MRIsetVoxVal(mriT,x,y,z,0,std::numeric_limits<float>::quiet_NaN());

      }

  // we will need the derivatives
  if (verbose > 1) cout << " compute derivatives ... " << flush;
  MRI *Sfx=NULL,*Sfy=NULL,*Sfz=NULL,*Sbl=NULL;
  MRI *Tfx=NULL,*Tfy=NULL,*Tfz=NULL,*Tbl=NULL;
  MyMRI::getPartials(mriS,Sfx,Sfy,Sfz,Sbl);
  MyMRI::getPartials(mriT,Tfx,Tfy,Tfz,Tbl);

  if (verbose > 1) cout << " done!" << endl;

  if (verbose > 1) cout << " Subsample ... "<< flush;

  MRI * fx  = MyMRI::subSample(Tfx);
  MRI * fy  = MyMRI::subSample(Tfy);
  MRI * fz  = MyMRI::subSample(Tfz);
  MRI * ssb = MyMRI::subSample(Sbl);
  MRI * stb = MyMRI::subSample(Tbl);

  MRIfree(&Sfx);
  MRIfree(&Sfy);
  MRIfree(&Sfz);
  MRIfree(&Sbl);
  MRIfree(&Tfx);
  MRIfree(&Tfy);
  MRIfree(&Tfz);
  MRIfree(&Tbl);

  if (verbose > 1) cout << " done! " << endl;


  // allocate the space
  int n = fx->width * fx->height * fx->depth;
  int pnum = 12;
  if (rigid) pnum =6;
  if (iscale) pnum++;

  pair <MATRIX*, VECTOR* > Ab;

  MATRIX* A = MatrixAlloc(n,pnum,MATRIX_REAL);
  VECTOR* b = MatrixAlloc(n,1,MATRIX_REAL);

  int count = 0;
  int xp1,yp1,zp1;
  for (z = 0 ; z < fx->depth ; z++)
    for (x = 0 ; x < fx->width ; x++)
      for (y = 0 ; y < fx->height ; y++)
      {
        if (isnan(MRIFvox(fx, x, y, z)) || isnan(MRIFvox(fy, x, y, z)) || isnan(MRIFvox(fz, x, y, z)) || isnan(MRIFvox(ssb, x, y, z))|| isnan(MRIFvox(stb, x, y, z)))
          continue;

        count++;
        xp1 = 2*x+3;
        yp1 = 2*y+3;
        zp1 = 2*z+3;
        //zp1 = z+3;
        //cout << "x: " << x << " y: " << y << " z: " << z << " coutn: "<< count << endl;
        //cout << " " << count << " mrifx: " << MRIFvox(mri_fx, x, y, z) << " mrifx int: " << (int)MRIvox(mri_fx,x,y,z) <<endl;
        if (rigid)
        {
          assert(!rigid);

          //if (iscale)
          //  *MATRIX_RELT(A, count, 7) = -MRIFvox(Tbl, x, y, z);
          if (iscale) *MATRIX_RELT(A, count, 7) = MRIFvox(stb, x, y, z);
        }
        else // affine
        {
          *MATRIX_RELT(A, count, 1)  = MRIFvox(fx, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 2)  = MRIFvox(fx, x, y, z)*yp1;
          *MATRIX_RELT(A, count, 3)  = MRIFvox(fx, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 4)  = MRIFvox(fx, x, y, z);
          *MATRIX_RELT(A, count, 5)  = MRIFvox(fy, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 6)  = MRIFvox(fy, x, y, z)*yp1;
          *MATRIX_RELT(A, count, 7)  = MRIFvox(fy, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 8)  = MRIFvox(fy, x, y, z);
          *MATRIX_RELT(A, count, 9)  = MRIFvox(fz, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 10) = MRIFvox(fz, x, y, z)*yp1;
          *MATRIX_RELT(A, count, 11) = MRIFvox(fz, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 12) = MRIFvox(fz, x, y, z);
          //if (iscale) *MATRIX_RELT(A, count, 13) = MRIFvox(ssb, x, y, z);
          if (iscale) *MATRIX_RELT(A, count, 13) = MRIFvox(stb, x, y, z);
        }

        *MATRIX_RELT(b, count, 1) =  (MRIFvox(ssb, x, y, z) - MRIFvox(stb, x, y, z));

      }




  // adjust sizes
  Ab.first  = MatrixAlloc(count,pnum,MATRIX_REAL);
  Ab.second = MatrixAlloc(count,1,MATRIX_REAL);
  for (int rr = 1; rr<= count; rr++)
  {
    *MATRIX_RELT(Ab.second, rr, 1) = *MATRIX_RELT(b, rr, 1);
    for (int cc = 1; cc <= pnum; cc++)
    {
      *MATRIX_RELT(Ab.first, rr, cc) = *MATRIX_RELT(A, rr, cc);
    }
  }

  if (verbose > 1) cout << " Considering: " << count << " non-zero voxels"<< endl;

  MatrixFree(&A);
  MatrixFree(&b);
  MRIfree(&fx);
  MRIfree(&fy);
  MRIfree(&fz);
  MRIfree(&ssb);
  MRIfree(&stb);
  return Ab;

}

MATRIX* Registration::constructR(MATRIX* p)
// Construct restriction matrix (to restrict the affine problem to less parameters)
// if p->rows == 6 use only rigid
// if p->rows == 7 use also intensity scale
// if p->rows == 3 use only trans
// if p->rows == 4 use only trans + intensity
{
  assert(p != NULL);
  assert((p->rows == 6 || p->rows==7) && p->cols ==1);

  int adim = 12;
  if (iscale)
  {
    assert(p->rows == 7 || p->rows ==4);
    adim++;
  }
  MATRIX* R = MatrixAlloc(p->rows,adim,MATRIX_REAL);
  MatrixClear(R);


  // translation p1,p2,p3 map to m4,m8,m12
  *MATRIX_RELT(R, 1,  4) = 1.0;
  *MATRIX_RELT(R, 2,  8) = 1.0;
  *MATRIX_RELT(R, 3, 12) = 1.0;

  // iscale (p7 -> m13)
  if (p->rows ==7) *MATRIX_RELT(R, 7, 13) = 1.0;
  if (p->rows ==4) *MATRIX_RELT(R, 4, 13) = 1.0;

  if (p->rows <=4) return R;

  // rotation derivatives (dm_i/dp_i)
  double s4 = sin(*MATRIX_RELT(p, 4, 1));
  double c4 = cos(*MATRIX_RELT(p, 4, 1));
  double s5 = sin(*MATRIX_RELT(p, 5, 1));
  double c5 = cos(*MATRIX_RELT(p, 5, 1));
  double s6 = sin(*MATRIX_RELT(p, 6, 1));
  double c6 = cos(*MATRIX_RELT(p, 6, 1));

  *MATRIX_RELT(R, 5,  1) = -s5*c6;
  *MATRIX_RELT(R, 6,  1) = -c5*s6;

  *MATRIX_RELT(R, 5,  2) = -s5*s6;
  *MATRIX_RELT(R, 6,  2) =  c5*c6;

  *MATRIX_RELT(R, 5,  3) = -c5;

  *MATRIX_RELT(R, 4,  5) =  c4*s5*c6+s4*s6;
  *MATRIX_RELT(R, 5,  5) =  s4*c5*c6;
  *MATRIX_RELT(R, 6,  5) = -s4*s5*s6-c4*c6;

  *MATRIX_RELT(R, 4,  6) =  c4*s5*s6-s4*c6;
  *MATRIX_RELT(R, 5,  6) =  s4*c5*s6;
  *MATRIX_RELT(R, 6,  6) =  s4*s5*c6-c4*s6;

  *MATRIX_RELT(R, 4,  7) =  c4*c5;
  *MATRIX_RELT(R, 5,  7) = -s4*s5;

  *MATRIX_RELT(R, 4,  9) = -s4*s5*c6+c4*s6;
  *MATRIX_RELT(R, 5,  9) =  c4*c5*c6;
  *MATRIX_RELT(R, 6,  9) = -c4*s5*s6+s4*c6;

  *MATRIX_RELT(R, 4, 10) = -s4*s5*s6-c4*c6;
  *MATRIX_RELT(R, 5, 10) =  c4*c5*s6;
  *MATRIX_RELT(R, 6, 10) =  c4*s5*c6+s4*s6;

  *MATRIX_RELT(R, 4, 11) = -s4*c5;
  *MATRIX_RELT(R, 5, 11) = -c4*s5;

  return R;
}


MATRIX * Registration::rt2mat(MATRIX * r, MATRIX * t, MATRIX *outM)
// converts rot vector (3x1) and translation vector (3x1)
// into an affine matrix (homogeneous coord) 4x4
// if global rtype ==1 r1,r2,r3 are as in robust paper (axis, and length is angle)
// if global rtype ==2 then r1,r2,r3 are angles around x,y,z axis (order 1zrot,2yrot,3xrot)
{
  if (outM == NULL)
    outM = MatrixAlloc(4, 4, MATRIX_REAL);

  assert(r->rows == 3 && r->cols == 1);
  assert(t->rows == 3 && t->cols == 1);

  MATRIX *rmat;

  if (rtype == 2)
  {
//      MATRIX* rx = MatrixAllocRotation(3,*MATRIX_RELT(r, 1, 1),X_ROTATION);
//      MATRIX* ry = MatrixAllocRotation(3,*MATRIX_RELT(r, 2, 1),Y_ROTATION);
//      MATRIX* rz = MatrixAllocRotation(3,*MATRIX_RELT(r, 3, 1),Z_ROTATION);
//      MATRIX* tmp= MatrixMultiply(rx,ry,NULL);
//      rmat = MatrixMultiply(tmp,rz,NULL);
//      MatrixFree(&rx); MatrixFree(&ry); MatrixFree(&rz); MatrixFree(&tmp);
//      MatrixPrintFmt(stdout,"% 2.8f",rmat);
//      cout << endl;


    // first convert rotation to quaternion (clockwise)
    Quaternion q;
    q.importZYXAngles(-*MATRIX_RELT(r, 3, 1), -*MATRIX_RELT(r, 2, 1), -*MATRIX_RELT(r, 1, 1));
    // then to rotation matrix
    rmat = MyMatrix::getMatrix(q.getRotMatrix3d(),3);
    //MatrixPrintFmt(stdout,"% 2.8f",rmat2);

  }
  else if (rtype ==1)
  {

    // first convert rotation to quaternion
    Quaternion q;
    q.importRotVec(*MATRIX_RELT(r, 1, 1),*MATRIX_RELT(r, 2, 1),*MATRIX_RELT(r, 3, 1));
    // then to rotation matrix
    rmat = MyMatrix::getMatrix(q.getRotMatrix3d(),3);

  }
  else assert (1==2);

  int rr, cc;
  for (rr=1;rr<=3;rr++)
  {
    for (cc=1;cc<=3;cc++) // copy rot-matrix
      *MATRIX_RELT(outM, rr, cc) = *MATRIX_RELT(rmat, rr, cc);

    // copy translation into 4th column
    *MATRIX_RELT(outM, rr, 4) = *MATRIX_RELT(t, rr, 1);
    // set 4th row to zero
    *MATRIX_RELT(outM, 4, rr) = 0.0;
  }
  //except 4,4
  *MATRIX_RELT(outM, 4, 4) = 1.0;

  MatrixFree(&rmat);

  return outM;
}


MATRIX * Registration::p2mat(MATRIX * p6, MATRIX *outM)
// converts trans-rot vector (6x1)
// into an affine matrix (homogeneous coord) 4x4
// if rtype ==2 , then p4,p5,p6 are angles around x,y,z axis
// else rtype ==1 they are as in robust paper
{
  assert(p6->rows == 6 && p6->cols == 1);
  MATRIX* t = MatrixAlloc(3, 1, MATRIX_REAL);
  MATRIX* r = MatrixAlloc(3, 1, MATRIX_REAL);
  for (int rr = 1;rr<=3;rr++)
  {
    *MATRIX_RELT(t, rr, 1) = *MATRIX_RELT(p6, rr, 1);
    *MATRIX_RELT(r, rr, 1) = *MATRIX_RELT(p6, rr+3, 1);
  }

  outM = rt2mat(r,t,outM);
  MatrixFree(&r);
  MatrixFree(&t);
  return outM;
}

pair < MATRIX*, double > Registration::convertP2Md(MATRIX* p)
// rtype : use restriction (if 2) or rigid from robust paper
// returns registration as 4x4 matrix M, and iscale
{
//   cout << " Registration::convertP2Md(MATRIX* p) (p->rows: " << p->rows << " )" << flush;
  pair < MATRIX*, double> ret(NULL,1.0);
  MATRIX* pt;

  if (p->rows == 4 ||p->rows == 7 || p->rows == 13) // iscale
  {
    //cout << " has intensity " << endl;
    // cut off intensity scale
		
    //ret.second = 1.0-*MATRIX_RELT(p, p->rows, 1);		
    ret.second = 1.0/(1.0+*MATRIX_RELT(p, p->rows, 1));
		
    pt=  MatrixAlloc(p->rows -1, 1,MATRIX_REAL);
    for (int rr = 1; rr<p->rows; rr++)
      *MATRIX_RELT(pt, rr, 1) = *MATRIX_RELT(p, rr, 1);
  }
  else pt = MatrixCopy(p,NULL);

  if (pt->rows == 12) ret.first = MyMatrix::aff2mat(pt,NULL);
  else if (pt->rows == 6)
  {
    ret.first = p2mat(pt,NULL);
  }
  else if (pt->rows ==3)
  {
    ret.first = MatrixIdentity(4,NULL);
    *MATRIX_RELT(ret.first, 1, 4) = *MATRIX_RELT(pt, 1, 1);
    *MATRIX_RELT(ret.first, 2, 4) = *MATRIX_RELT(pt, 2, 1);
    *MATRIX_RELT(ret.first, 3, 4) = *MATRIX_RELT(pt, 3, 1);
  }
  else
  {
    cerr << " parameter neither 3,6 nor 12 : " << pt->rows <<" ??" << endl;
    assert(1==2);
  }

  MatrixFree(&pt);
//   cout << " -- DONE " << endl;
  return ret;
}



vector < MRI* > Registration::buildGaussianPyramid (MRI * mri_in, int n)
{
  vector <MRI* > p (n);
  MRI * mri_tmp;
// if (mri_in->type == MRI_UCHAR) cout << " MRI_UCHAR" << endl;
// else cout << " type: " << mri_in->type << endl;



  MRI *mri_kernel ;
  mri_kernel = MRIgaussian1d(1.08, 5) ;
  //mri_kernel = MRIalloc(5,1,1, MRI_FLOAT);
  MRIFvox(mri_kernel, 0, 0, 0) =  0.0625 ;
  MRIFvox(mri_kernel, 1, 0, 0) =  0.25 ;
  MRIFvox(mri_kernel, 2, 0, 0) =  0.375 ;
  MRIFvox(mri_kernel, 3, 0, 0) =  0.25 ;
  MRIFvox(mri_kernel, 4, 0, 0) =  0.0625 ;

  mri_tmp = mri_in;
  int i;
  int min = 16; // stop when we are smaller than min
  p[0] = MRIconvolveGaussian(mri_in, NULL, mri_kernel);
  for (i = 1;i<n;i++)
  {
    mri_tmp = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel) ;
    p[i] = MRIdownsample2(mri_tmp,NULL);
    MRIfree(&mri_tmp);
    mri_tmp = p[i];
    if (p[i]->width < min || p[i]->height <min || p[i]->depth <min)
      break;
  }
  if (i<n) p.resize(i+1);

  MRIfree(&mri_kernel);

  return p;
}

void Registration::freeGaussianPyramid(vector< MRI* >& p)
{
  for (uint i = 0;i<p.size();i++)
    MRIfree(&p[i]);
  p.clear();
}


// ---------------------- Initial Transform using Moments -----------------------------

MATRIX * Registration::initializeTransform(MRI *mriS, MRI *mriT)
{
  cout << "   - computing initial Transform\n" ;

  MATRIX* Minit = MRIgetVoxelToVoxelXform(mriS,mriT) ;
  MRI * mri_tmp = MRIlinearTransform(mriS,NULL,Minit);

  // find centroids:
  vector < double > centroidS = CostFunctions::centroid(mri_tmp);
  vector < double > centroidT = CostFunctions::centroid(mriT);


  //bool initorient = false; // do not use orientation (can be off due to different cropping)
  //bool initorient = true;
  if (initorient)
  {
    // find orientation:
    MATRIX * evT = CostFunctions::orientation(mriT);
    MATRIX * evS = CostFunctions::orientation(mri_tmp);

    // adjust orientation (flip) to target
    int fcount = 0;
    for (int c=1;c<4;c++)
    {
      double d = 0;
      for (int r=1;r<4;r++)
        d += evT->rptr[r][c] * evS->rptr[r][c];
      if (fabs(d)<0.1)
      {
        cout << "        Over 84 degree difference ( " << 360*acos(fabs(d))/M_PI << " ) in column " << c << ", skipping pre-orientation" << endl;
        initorient =false;
        break;
      }
      if (d<0) // flip
      {
        fcount++;
        for (int r=1;r<4;r++)
          evS->rptr[r][c] = - evS->rptr[r][c];
      }

    }
    assert(fcount ==0 || fcount==2); // make sure det > 0 still

    if (debug)
    {
      cout << endl <<" evecs S: " << endl;
      MatrixPrintFmt(stdout,"% 2.8f",evS);

      cout << endl <<" evecs T: " << endl;
      MatrixPrintFmt(stdout,"% 2.8f",evT);
    }

    if (initorient) // use orientation info for intial alignment
    {
      cout << "     -- using orientation info\n" ;
      // make 4x4
      MATRIX * evS2 = MatrixIdentity(4,NULL);
      MATRIX * evT2 = MatrixIdentity(4,NULL);
      for (int c=1; c<4;c++)
        for (int r=1; r<4;r++)
        {
          evS2->rptr[r][c] = evS->rptr[r][c];
          evT2->rptr[r][c] = evT->rptr[r][c];
        }

      // first move Src centroid to origin
      Minit->rptr[1][4] -= centroidS[0];
      Minit->rptr[2][4] -= centroidS[1];
      Minit->rptr[3][4] -= centroidS[2];

      // Rotate evT * evS^-1
      MATRIX* evS2i = MatrixInverse(evS2,NULL);
      assert(evS2i);
      Minit = MatrixMultiply(evS2i,Minit,Minit);
      Minit = MatrixMultiply(evT2,Minit,Minit);

      // finally move to T centroid
      Minit->rptr[1][4] += centroidT[0];
      Minit->rptr[2][4] += centroidT[1];
      Minit->rptr[3][4] += centroidT[2];

      if (debug)
      {
        mri_tmp = MRIlinearTransform(mriS,mri_tmp,Minit);
        MRIwrite(mri_tmp,"init-align-rot.mgz");
      }
      MatrixFree(&evS2);
      MatrixFree(&evS2i);
      MatrixFree(&evT2);
    }
    if (evT) MatrixFree(&evT);
    if (evS) MatrixFree(&evS);

  }

  if (!initorient) // if orientation did not work
  {
    cout << "     -- using translation info\n" ;
    Minit->rptr[1][4] += centroidT[0]-centroidS[0];
    Minit->rptr[2][4] += centroidT[1]-centroidS[1];
    Minit->rptr[3][4] += centroidT[2]-centroidS[2];
    if (debug)
    {
      mri_tmp = MRIlinearTransform(mriS,mri_tmp,Minit);
      MRIwrite(mri_tmp,"init-align-trans.mgz");
    }
  }

  MRIfree(&mri_tmp);

  return Minit;
}



// NOT TESTED !!!!!!!
#define MAX_DX   1.2
#define MAX_DY   1.2
#define MAX_DZ   1.2
#define MIN_DX   (1.0/MAX_DX)
#define MIN_DY   (1.0/MAX_DY)
#define MIN_DZ   (1.0/MAX_DZ)
#define MAX_RATIO 1.2
int Registration::init_scaling(MRI *mri_in, MRI *mri_ref, MATRIX *m_L)
{
  MATRIX      *m_scaling ;
  float       sx, sy, sz, dx, dy, dz ;
  MRI_REGION  in_bbox, ref_bbox ;

  m_scaling = MatrixIdentity(4, NULL) ;

  MRIboundingBox(mri_in, 60, &in_bbox) ;
  MRIboundingBox(mri_ref, 60, &ref_bbox) ;
  sx = (float)ref_bbox.dx / (float)in_bbox.dx ;
  sy = (float)ref_bbox.dy / (float)in_bbox.dy ;
  sz = (float)ref_bbox.dz / (float)in_bbox.dz ;
  dx = (ref_bbox.x+ref_bbox.dx-1)/2 - (in_bbox.x+in_bbox.dx-1)/2 ;
  dy = (ref_bbox.y+ref_bbox.dy-1)/2 - (in_bbox.y+in_bbox.dy-1)/2 ;
  dz = (ref_bbox.z+ref_bbox.dz-1)/2 - (in_bbox.z+in_bbox.dz-1)/2 ;

  if (sx > MAX_DX)
    sx = MAX_DX ;
  if (sx < MIN_DX)
    sx = MIN_DX ;
  if (sy > MAX_DY)
    sy = MAX_DY ;
  if (sy < MIN_DY)
    sy = MIN_DY ;
  if (sz > MAX_DZ)
    sz = MAX_DZ ;
  if (sz < MIN_DZ)
    sz = MIN_DZ ;
//  if (Gdiag & DIAG_SHOW)
  fprintf(stderr, "initial scaling: (%2.2f, %2.2f, %2.2f) <-- "
          "(%d/%d,%d/%d,%d/%d)\n",
          sx,sy,sz, ref_bbox.dx, in_bbox.dx, ref_bbox.dy, in_bbox.dy,
          ref_bbox.dz, in_bbox.dz) ;
  *MATRIX_RELT(m_scaling, 1, 1) = sx ;
  *MATRIX_RELT(m_scaling, 2, 2) = sy ;
  *MATRIX_RELT(m_scaling, 3, 3) = sz ;

#if 0
  *MATRIX_RELT(m_L, 1, 4) = dx ;
  *MATRIX_RELT(m_L, 2, 4) = dy ;
  *MATRIX_RELT(m_L, 3, 4) = dz ;
#endif
  MatrixMultiply(m_scaling, m_L, m_L) ;
  return(NO_ERROR) ;
}

void Registration::setSource (MRI * s, bool fixvoxel, bool fixtype)
//local copy
{
  if (! (fixvoxel || fixtype)) mri_source = MRIcopy(s,mri_source);
  else mri_source = MyMRI::makeConform(s,NULL,fixvoxel,fixtype);
  if (gpS.size() > 0) freeGaussianPyramid(gpS);
  //cout << "mri_source" << mri_source << endl;
}

void Registration::setTarget (MRI * t, bool fixvoxel, bool fixtype)
//local copy
{
  if (! (fixvoxel || fixtype))   mri_target = MRIcopy(t,mri_target);
  else mri_target = MyMRI::makeConform(t,NULL,fixvoxel,fixtype);
  if (gpT.size() > 0) freeGaussianPyramid(gpT);
  //cout << "mri_target" << mri_target << endl;
}

void Registration::setName(const std::string &n)
// set name and nbase (base name without path)
{
  name = n;
  nbase = n;
  int rf = nbase.rfind("/");
  if (rf != -1)
  {
    nbase = nbase.substr(rf+1,nbase.length());
  }
}
