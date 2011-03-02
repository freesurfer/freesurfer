/**
 * @file RegPowell.cpp
 * @brief A class to compute a robust registration using Powell
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:24 $
 *    $Revision: 1.8 $
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
#include "RegPowell.h"
#include "CostFunctions.h"
#include "MyMatrix.h"
#include "MyMRI.h"

#include "numerics.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vnl/vnl_inverse.h>
#include "RegistrationStep.h"

using namespace std;

RegPowell* RegPowell::tocurrent = NULL;
MRI * RegPowell::scf = NULL;
MRI * RegPowell::tcf = NULL;
int RegPowell::pcount = 0;
vnl_matrix_fixed < double, 4, 4> RegPowell::mh1 = vnl_matrix_fixed < double, 4, 4>();
vnl_matrix_fixed < double, 4, 4> RegPowell::mh2 = vnl_matrix_fixed < double, 4, 4>();
int RegPowell::icount = 0;

float RegPowell::costFunction(float p[])// starting at 1
{
//cout << " costFunction (" << flush;


  // copy p into vector
  //MATRIX * v = MatrixAlloc(pcount,1,MATRIX_REAL);
	vnl_vector < double > v (pcount);
  for (int i = 0;i<pcount; i++)
  {
    v[i] = p[i+1];// p starting 1
//      cout << p[i] << " ";
  }
//    cout << ") = " << flush;

  // transform into matrix and iscale double
  //pair < vnl_matrix_fixed < double , 4, 4 > , double > Md = tocurrent->convertP2Md(v);
  pair < vnl_matrix_fixed < double , 4, 4 > , double > Md = RegistrationStep<double>::convertP2Md(v,tocurrent->rtype);

  // new full M = mh2 * cm * mh1
  //Md.first = MatrixMultiply(mh2,Md.first,Md.first);
  //Md.first = MatrixMultiply(Md.first,mh1,Md.first);
	Md.first = mh2 * Md.first * mh1;

  // compute new half way maps
  vnl_matrix_fixed < double , 4, 4 > mh  = MyMatrix::MatrixSqrt(Md.first);
  // do not just assume m = mh*mh, rather m = mh2 * mh
  // for transforming target we need mh2^-1 = mh * m^-1
  //MATRIX * mi  = MatrixInverse(Md.first,NULL);
  //MATRIX * mhi = MatrixMultiply(mh1,mi,NULL);
  vnl_matrix_fixed < double , 4, 4 > mhi = mh1 * vnl_inverse(Md.first);
	
  // map both to half way
  MRI* mri_Swarp = MRIclone(scf,NULL);
  mri_Swarp = MyMRI::MRIlinearTransform(scf,mri_Swarp, mh);
  MRI* mri_Twarp = MRIclone(scf,NULL); // bring them to same space (just use src geometry)
  mri_Twarp = MyMRI::MRIlinearTransform(tcf,mri_Twarp, mhi);

  // adjust intensity
  if (tocurrent->iscale)
  {
    //cout << "   - adjusting intensity ( "<< fmd.second << " ) " << endl;

    MyMRI::MRIvalscale(mri_Swarp,mri_Swarp,(1.0+Md.second)*0.5);
    MyMRI::MRIvalscale(mri_Twarp,mri_Twarp,(1.0+ 1.0/Md.second)*0.5);
  }

  // compute error
  float d;
  if (tocurrent->robust) d = (float) CostFunctions::tukeyBiweight(mri_Swarp,mri_Twarp,tocurrent->sat);
  else d = (float) CostFunctions::leastSquares(mri_Swarp,mri_Twarp);

  icount++;
  if (icount%100 == 0) cout << "*" << flush;
//    if (icount%20 == 0)
//    {
//       cout << endl << " iteration : " << icount << endl;
//       MatrixPrintFmt(stdout,"% 2.8f",v);
//       cout << " error: " << d  << endl;
//       int ii = icount / 20;
//       std::stringstream out;
//       out << ii;
//       MRIwrite(mri_Swarp,("shw-"+out.str()+".mgz").c_str());
//       MRIwrite(mri_Twarp,("thw-"+out.str()+".mgz").c_str());
//    }

  // clean up
  MRIfree(&mri_Swarp);
  MRIfree(&mri_Twarp);
  //MatrixFree(&v);
  //MatrixFree(&Md.first);
  //MatrixFree(&mhi);
  //MatrixFree(&mh);
  //MatrixFree(&mi);

  return d;
}

void RegPowell::computeIterativeRegistration( int nmax,double epsit,MRI * mriS, MRI* mriT, const vnl_matrix < double > & m , double iscaleinit)
// retruns 4x4 matrix and iscale value
{
  if (!mriS) mriS = mri_source;
  if (!mriT) mriT = mri_target;

  assert (mriS && mriT);

  tocurrent = this; // so that we can access this from static cost function

  pair < vnl_matrix_fixed < double, 4, 4> , double > fmd(vnl_matrix_fixed < double, 4 , 4> () ,iscaleinit);

  // check if mi (inital transform) is passed
  if (!m.empty()) fmd.first = m;
  else if (!Minit.empty()) fmd.first = Minit;
  else fmd.first = initializeTransform(mriS,mriT) ;

  if (debug > 0)
  {
    cout << "   - initial transform:\n" ;
    //MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
		cout << fmd.first << endl;
  }


  // here maybe better to symmetrically warp both images SQRT(M)
  // this keeps the problem symmetric
  cout << "   - warping source and target (sqrt)" << endl;
  //if (mh1) MatrixFree(&mh1);
  mh1 = MyMatrix::MatrixSqrt(fmd.first);
  // do not just assume m = mh*mh, rather m = mh2 * mh
  // for transforming target we need mh2^-1 = mh * m^-1
  //MATRIX * mi  = MatrixInverse(fmd.first,NULL);
  //MATRIX * mhi = MatrixMultiply(mh1,mi,NULL);
  vnl_matrix_fixed < double, 4, 4 > mhi = mh1 * vnl_inverse(fmd.first);
	
  //set static
  //if (mh2) MatrixFree(&mh2);
  //mh2 = MatrixInverse(mhi,NULL); // M = mh2 * mh1
  mh2 = vnl_inverse(mhi); // M = mh2 * mh1

  //if (mri_Swarp) MRIfree(&mri_Swarp);
  //mri_Swarp = MRIclone(mriS,NULL);
  //mri_Swarp = MRIlinearTransform(mriS,mri_Swarp, mh);
  //if (mri_Twarp) MRIfree(&mri_Twarp);
  //mri_Twarp = MRIclone(mriS,NULL); // bring them to same space (just use src geometry)
  //mri_Twarp = MRIlinearTransform(mriT,mri_Twarp, mhi);

 // //MatrixFree(&mh);
  //MatrixFree(&mhi);
 // MatrixFree(&mi);

  // adjust intensity later in first powell call

  scf = mriS;
  tcf = mriT;

  // create parameter vector:
  pcount = 3; // transolny
  if (rigid) pcount = 6;
  else pcount = 12;
  if (pcount==3) assert(transonly);
  if (iscale) pcount++;


  // compute Registration
  cout << "   - compute new registration ( " << pcount << " params )" << endl;
  float fret, fstart, min_sse;
  float* p   = ::vector(1, pcount+1) ;
  float** xi = ::matrix(1, pcount+1, 1, pcount+1) ;
  float tol = 1e-5; //-8
  int maxiter = 36;
  int iter;

  for (int i = 1;i<=pcount;i++)
  {
    p[i] = 0.0;

    for (int j = i;j<=pcount;j++)
    {
      xi[i][j] = 0.0;
      xi[j][i] = 0.0;
    }
    xi[i][i] = 1.0;
  }
  if (iscale) p[pcount] = iscaleinit;

  min_sse = costFunction(p) ;
  cout << " min_sse: " << min_sse << endl;
  icount = 0;
  //OpenPowell(p, xi, pcount, tol, &iter, &fret, costFunction);
  OpenPowell2(p, xi, pcount, tol,tol,maxiter, &iter, &fret, costFunction);
  cout << endl << "best alignment initial powell: " << fret << " (" << iter << " steps)" << endl;

  int count = 0;
  do
  {
    count++;
    // reinitialize powell directions
    for (int r = 1 ; r <= pcount ; r++)
    {
      for (int c = 1 ; c <= pcount ; c++)
      {
        xi[r][c] = r == c ? 1 : 0 ;
      }
    }

    fstart = fret ;
    icount = 0;
    OpenPowell(p, xi, pcount, tol, &iter, &fret, costFunction);
    cout << endl << " return from " << count << " powell!!!" << endl;


//    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    //   printf("best alignment after powell: %2.3f (%d steps)\n", fret, iter) ;
    cout << "best alignment after powell: " << fret << " (" << iter << " steps)" << endl;

    if ((fstart-fret)/fstart < tol)
      break ;
  }
  while (fret < fstart) ;

  vnl_vector < double > v(RegPowell::pcount);
  for (int i = 1;i<=RegPowell::pcount; i++)
    v[i-1] = p[i];
 // MatrixPrintFmt(stdout,"% 2.8f",v);
	cout << v << endl;

  //if (fmd.first) MatrixFree(&fmd.first);
  fmd = RegistrationStep<double>::convertP2Md(v,rtype);
  //MatrixFree(&v);
  //MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
	cout << fmd.first << endl;

  free_matrix(xi, 1, pcount+1, 1, pcount+1) ;
  free_vector(p, 1, pcount+1) ;

  // new full M = mh2 * cm * mh1
  //fmd.first = MatrixMultiply(mh2,fmd.first,fmd.first);
  //fmd.first = MatrixMultiply(fmd.first,mh1,fmd.first);
  fmd.first = mh2 * fmd.first * mh1;
  Mfinal = fmd.first;
  iscalefinal = fmd.second;

  cout << endl << " DONE " << endl;
  cout << endl << "Final Transform:" << endl;
  
  //MatrixPrintFmt(stdout,"% 2.8f",Mfinal);
	cout << Mfinal << endl;

  MRI* mri_Swarp = MRIclone(mriT,NULL);
  mri_Swarp = MyMRI::MRIlinearTransform(mriS,mri_Swarp,Mfinal);
  MRIwrite(mriS,"mriS.mgz");
  MRIwrite(mri_Swarp,"mriSwarp.mgz");
  MRIwrite(mriT,"mriT.mgz");

//  exit(1);

//  return fmd ;


//
//        //p = computeRegistrationStepP(mri_Swarp,mri_Twarp);
//        if (pw.second) MRIfree(&pw.second);
//        pw = computeRegistrationStepW(mri_Swarp,mri_Twarp);
//
//        if (cmd.first != NULL) MatrixFree(&cmd.first);
//        cmd = convertP2MATRIXd(pw.first);
//        if (lastp) MatrixFree(&lastp);
//        lastp = pw.first;
//        pw.first = NULL;
//
//
//        // store M and d
//        cout << "   - store transform" << endl;
//        //cout << endl << " current : Matrix: " << endl;
//        //MatrixPrintFmt(stdout,"% 2.8f",cmd.first);
//        //cout << " intens: " << cmd.second << endl;
//        MATRIX* fmdtmp = MatrixCopy(fmd.first,NULL);
//        //fmd.first = MatrixMultiply(cmd.first,fmd.first,fmd.first); //old
//        if(mh2) MatrixFree(&mh2);
//        mh2 = MatrixInverse(mhi,NULL); // M = mh2 * mh
//        // new M = mh2 * cm * mh
//        fmd.first = MatrixMultiply(mh2,cmd.first,fmd.first);
//        fmd.first = MatrixMultiply(fmd.first,mh,fmd.first);
//
//        fmd.second *= cmd.second;
//        //cout << endl << " Matrix: " << endl;
//        //MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
//        if (!rigid) diff = getFrobeniusDiff(fmd.first, fmdtmp);
//        else        diff = sqrt(RigidTransDistSq(fmd.first, fmdtmp));
//        cout << "     -- old difference to prev. transform: " << diff << endl;
//        diff = sqrt(AffineTransDistSq(fmd.first, fmdtmp, 100));
//        cout << "     -- difference to prev. transform: " << diff << endl;
//        //cout << " intens: " << fmd.second << endl;
//
//        MatrixFree(&fmdtmp);
//        MatrixFree(&mi);
//        //MatrixFree(&mh);
//        //MatrixFree(&mhi);
//        //MatrixFree(&mh2);
//
//    }
//
//    //   DEBUG OUTPUT
//    if (debug > 0)
//    {
//     // write weights and warped images after last step:
//
//        MRIwrite(mri_Swarp,(name+"-mriS-warp.mgz").c_str());
//        MRIwrite(mri_Twarp,(name+"-mriT-warp.mgz").c_str());
//        MRI* salign = MRIclone(mriS,NULL);
//        salign = MRIlinearTransform(mri_Swarp, salign,cmd.first);
//        MRIwrite(salign,(name+"-mriS-align.mgz").c_str());
//        MRIfree(&salign);
//        if (pw.second)
//        {
//      // in the half-way space:
//             string n = name+string("-mriS-weights.mgz");
//             MRIwrite(pw.second,n.c_str());
//        }
//     }
//
//    // store weights (mapped to target space):
//    if (pw.second)
//    {
//          // remove negative weights (markers) set to 1
//          int x,y,z;
//          for (z = 0 ; z < pw.second->depth  ; z++)
//          for (x = 0 ; x < pw.second->width  ; x++)
//          for (y = 0 ; y < pw.second->height ; y++)
//          {
//             if (MRIFvox(pw.second,x,y,z) < 0) MRIFvox(pw.second,x,y,z) = 1;
//          }
//    MRI * mtmp = MRIalloc(mriT->width,mriT->height,mriT->depth,MRI_FLOAT);
//   MRIcopyHeader(mriT,mtmp);
//   mtmp->type = MRI_FLOAT;
//    mtmp = MRIlinearTransform(pw.second,mtmp,mh2);
//          //MRIwrite(mtmp,weightsname.c_str());
//          MRIfree(&pw.second);
//   pw.second = mtmp;
//    }
//    if (mri_weights) MRIfree(&mri_weights);
//    mri_weights = pw.second;
//    if (mov2weights) MatrixFree(&mov2weights);
//    if (dst2weights) MatrixFree(&dst2weights);
//    mov2weights = mh; // no freeing needed
//    dst2weights = mhi;
//
//    if (diff > epsit) // adjust mh and mhi to new midpoint
//    {
//       cout << "     -- adjusting half-way maps " << endl;
//       MATRIX * ch = MatrixSqrt(fmd.first);
//       // do not just assume c = ch*ch, rather c = ch2 * ch
//       // for transforming target we need ch2^-1 = ch * c^-1
//       MATRIX * ci  = MatrixInverse(cmd.first,NULL);
//       MATRIX * chi = MatrixMultiply(ch,ci,NULL);
//       // append ch or chi to mh mhi
//       mov2weights = MatrixMultiply(ch,mh,NULL);
//       dst2weights = MatrixMultiply(chi,mhi,NULL);
//       MatrixFree(&mh);
//       MatrixFree(&mhi);
//    }
//
//    MRIfree(&mri_Twarp);
//    MRIfree(&mri_Swarp);
//    MatrixFree(&mh2);
//    if (cmd.first != NULL) MatrixFree(&cmd.first);
//
//
//    Mfinal = MatrixCopy(fmd.first, Mfinal);
//    iscalefinal = fmd.second;
//
//    return fmd;


}
