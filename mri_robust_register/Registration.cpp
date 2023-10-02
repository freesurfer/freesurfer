/**
 * @brief A class to compute a robust symmetric registration
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

#include "Registration.h"
#include "RegistrationStep.h"
#include "Quaternion.h"
#include "MyMatrix.h"
#include "MyMRI.h"
#include "Regression.h"
#include "CostFunctions.h"
#include "mriBSpline.h"

#include <limits>
#include <cassert>
#include <fstream>
#include <sstream>
#if ITK_VERSION_MAJOR >= 5
#include <iostream>
#include <vcl_compiler.h>
#else
#include <vcl_iostream.h>
#endif
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matlab_print.h>
#include <vnl/algo/vnl_determinant.h>

#include "limits.h"
#include "error.h"
#include "macros.h"
#include "mrimorph.h"

using namespace std;

/** Cleanup our private variables: mri_source, mri_target, Gaussian Pyramids and transform */
Registration::~Registration()
{ 
//std::cout << " Destroy Registration" << std::endl;
  if (mri_source)
    MRIfree(&mri_source);
  if (mri_target)
    MRIfree(&mri_target);
//  if (mri_indexing) MRIfree(&mri_indexing);
//  if (mri_weights)  MRIfree(&mri_weights);
//  if (mri_hweights) MRIfree(&mri_hweights);
  if (gpS.size() > 0)
    freeGaussianPyramid(gpS);
  if (gpT.size() > 0)
    freeGaussianPyramid(gpT);
  if (trans)
    delete trans;
  //std::cout << " Done " << std::endl;
}

/**
  Initialize registration (keep source and target and gauss pyramid)
 also keep Rsrc and Rtrg (resampling matrices, if exist).
 and keeps transformation type */
void Registration::clear()
{
//  if (trans) delete trans;
//  setRigid();
  costfun = ROB;
//  sat = -1;
  iscale = false;
//  rtype = 1;
  subsamplesize = -1;
  initorient = false;
  debug = 0;

//  if (mri_indexing) MRIfree(&mri_indexing);
//  if (mri_weights)  MRIfree(&mri_weights);
//  if (mri_hweights) MRIfree(&mri_hweights);
//  mri_weights= NULL;

  Minit.clear();
  Mfinal.clear();
  //lastp.clear();
  mov2weights.clear();
  dst2weights.clear();
}

/**
 Computes iterative registration as 4x4 matrix Mfinal and iscalefinal (class member).
 The caller needs to retrieve any really final transform with getFinalVox2Vox.
 */
void Registration::computeIterativeRegistration(int nmax, double epsit)
{

  // wraping private routine
  // parameters need to be set by caller vial set methods:
  computeIterativeRegistrationFull(nmax, epsit, mri_source, mri_target,
      getMinitResampled(), iscaleinit);

}

// void Registration::computeIterativeRegistration( int nmax,double epsit, MRI * mriS, MRI* mriT, const vnl_matrix < double >& m, double scaleinit)
// // private routine, as called from multiregistration (passing mriS and mriT...)
// // computes iterative registration (recomputing A and b in each step)
// // retruns 4x4 matrix Mfinal and iscalefinal (class member)
// // The caller needs to retrieve any really final transform with getFinalVox2Vox
// {
// 
//   // call helper to avoid code duplication:
//   
//   if (doubleprec)
//    iterativeRegistrationHelper<double>(nmax,epsit,mriS,mriT,m,scaleinit);
//   else
//    iterativeRegistrationHelper<float>(nmax,epsit,mriS,mriT,m,scaleinit);
//     
//   return;
//   
// }

// // update this to vnl:   or remove???
// pair < MATRIX*, double > Registration::computeIterativeRegSat( int n,double epsit, MRI * mriS, MRI* mriT, MATRIX* m, double scaleinit)
// // tests trough many saturations:
// {
// //   if (!mriS) mriS = mri_source;
// //   if (!mriT) mriT = mri_target;
// // 
// //   assert (mriS && mriT);
// // 
// //   pair < MATRIX*, double > cmd(NULL,1.0);
//    pair < MATRIX*, double > fmd(NULL,scaleinit);
// // 
// //   // check if mi (inital transform) is passed
// //   if (m)          fmd.first = MatrixCopy(m,NULL);
// //   else if (!Minit.empty()) fmd.first = MatrixCopy(Minit,NULL);
// //   else            fmd.first = initializeTransform(mriS,mriT) ;
// //
// //  if (scaleinit != 1.0) md.second = scaleinit;
// //  else md.second = iscaleinit;
// // 
// //   if (verbose >1) 
// //   {
// //     cout << "   - initial transform:\n" ;
// //     MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
// //   }
// // 
// //   double satval[20] = {20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
// //   vector < double > diffs (20);
// // 
// //   string nametmp = name;
// //   std::stringstream sout;
// // 
// //   for (int si = 0;si<(int)diffs.size();si++)
// //   {
// //     sat = satval[si];
// // 
// //     std::stringstream out;
// //     out << sat;
// //     name = nametmp+"-sat"+out.str();
// // 
// //     cmd = computeIterativeRegistration(n,epsit,mriS,mriT,fmd.first,fmd.second);
// // 
// //     if (!rigid) diffs[si] = MyMatrix::getFrobeniusDiff(fmd.first, cmd.first);
// //     else        diffs[si] = sqrt(MyMatrix::RigidTransDistSq(fmd.first, cmd.first));
// //       if (verbose >1) cout << "       difference on sat " << sat << " to prev. transform: " << diffs[si] << endl;
// // 
// //     fmd.second = cmd.second;
// //     MatrixFree(&fmd.first);
// //     fmd.first = cmd.first;
// // 
// //     // store transform
// //     LTA * lta = LTAalloc(1,mriS);
// //     lta->xforms[0].m_L = MRIvoxelXformToRasXform (mriS, mriT, fmd.first, lta->xforms[0].m_L) ;
// //     lta->type = LINEAR_RAS_TO_RAS ;
// //     getVolGeom(mriS, &lta->xforms[0].src);
// //     getVolGeom(mriT, &lta->xforms[0].dst);
// //     LTAwriteEx(lta, (name+".lta").c_str()) ;
// //     LTAfree(&lta);
// // 
// // 
// //   }
// // 
// //   name = nametmp;
// // 
// //   // plot diffs
// //   string fbase = name+"-sat";
// //   ofstream ofile((fbase+".plot").c_str(),ios::out);
// //   bool png = false;
// //   if (png) ofile << "set terminal png medium size 800,600" << endl;
// //   else ofile << "set terminal postscript eps color" << endl;
// //   if (png) ofile << "set output \""<< fbase <<".png\"" << endl;
// //   else ofile << "set output \""<< fbase <<".eps\"" << endl;
// //   ofile << "plot ";
// //   ofile << " \"-\" notitle with lines 1" << endl;
// //   for (int j = 0;j<(int)diffs.size(); j++)
// //   {
// //     ofile << -satval[j] << " " << diffs[j] << endl;
// //   }
// //   ofile << "e" << endl;
// //   ofile.close();
// // 
// // 
//    return fmd;
// }

// double Registration::findSaturation (MRI * mriS, MRI* mriT, const vnl_matrix < double > & mi , double scaleinit )
// {
//   if (verbose >0) cout << endl << endl << " Registration::findSaturation " << endl;
//   if (!mriS) mriS = mri_source;
//   if (!mriT) mriT = mri_target;
//   if (sat == -1) sat = 4.685; // set start value
//   
// //  if (! MyMRI::isConform(mriS) || ! MyMRI::isConform(mriT) )
// //  {
// //    cout << " WARNING: image(s) not conform!" << endl;
// //    cout << "     cannot estimate sat, will use default " << sat << endl;
// //    return sat;  
// //  }
// 
//   // if mriS and mriT have been passed, redo pyramid
//   if (mriS != mri_source)
//   {
//     if (gpS.size() > 0) freeGaussianPyramid(gpS);
//     gpS = buildGaussianPyramid(mriS,100); OLD!!!
//   }
//   if (mriT != mri_target)
//   {
//     if (gpT.size() > 0) freeGaussianPyramid(gpT);
//     gpT = buildGaussianPyramid(mriT,100); OLD!!!
//   }
// 
//   if (gpS.size() ==0) gpS = buildGaussianPyramid(mriS,100);
//   if (gpT.size() ==0) gpT = buildGaussianPyramid(mriT,100);
//   assert(gpS.size() == gpT.size());
// 
//   int resolution = gpS.size();
// 
//   vnl_matrix_fixed < double, 4, 4> m; m.set_identity();
// 
//   // variables to store matrix m and scaling factor d:
//   pair < vnl_matrix_fixed < double, 4, 4> , double > cmd;
//   pair < vnl_matrix_fixed < double, 4, 4> , double > md(vnl_matrix_fixed < double, 4, 4> (),scaleinit);
// 
//   // check if mi (inital transform) is passed
//   if (!mi.empty()) md.first = mi;
//   else if (!Minit.empty()) md.first = getMinitResampled();
//   else md.first = initializeTransform(mriS,mriT);
// 
//   if (verbose >1 ) 
//   {
//     cout << "   - initial transform:\n" ;
//     cout << md.first << endl;
//     cout << "   - initial iscale: " << scaleinit <<endl;
//   }
// 
//   // adjust md.first to current (lowest) resolution:
//   int rstart = 2; OLD!!!
//   for (int r = 1; r<=resolution-rstart; r++)
//     for (int rr = 0;rr<3;rr++)
//       md.first[rr][3]  = 0.5 *  md.first[rr][3];
//       
//   vnl_matrix_fixed < double, 4, 4> firstbackup = md.first;
// 
//   if (verbose >1 ) 
//   {
//     cout << "   - initial adjusted:\n" ;
//     cout << md.first << endl;
//   }
// 
//   bool iscaletmp = iscale;
// //  iscale = false; //disable intensity scaling on low resolutions
// 
//   int counter = 0;
//   for (int r = resolution-rstart;r>=2;r--)
//   {
// 
//     if (verbose >1 ) cout << endl << "Resolution: " << r << endl;
// 
// //    if (r==2) iscale = iscaletmp; // set iscale if set by user
// 
// 
//     // compute Registration
//     if (verbose >2 ) cout << "   - compute new iterative registration" << endl;
//     
//     int n = 3;
//     if (r==2) n = 1;
//     int vv = verbose;
//     if (verbose == 1) verbose = 0;
//     computeIterativeRegistration(n,0.05,gpS[r],gpT[r],md.first,md.second);
//     cmd.first = Mfinal;
//     cmd.second = iscalefinal;    
//     verbose = vv;
//     
//     if (verbose > 1)
//     {
//       cout << endl << " current : Matrix: " << endl;
//       cout << cmd.first << endl;
//       cout << " intens: " << cmd.second << endl;
// 
//       // adjust to highest level for output only:
//       double tx = cmd.first[0][3];
//       double ty = cmd.first[1][3];
//       double tz = cmd.first[2][3];
//       for (int ll = r; ll > 0; ll--)
//       {
//         tx *= 2;
//         ty*=2;
//         tz*=2;
//       }
//       cout << " equiv trans on highres: " << tx << " " << ty << " " << tz << endl;
//     }
//     
//     if (r == 2)
//     {
//        counter++;
//        if (verbose >1 ) cout << " Iteration: " << counter << endl;
//        if (debug)
//        {
//          // write out wcheck
//          string fn = getName() + "-wcheck-est.txt";
//          ofstream f(fn.c_str(),ios::out);
//          f << sat << " " << wcheck << endl;
//          f.close();  
//          string fn2 = getName() + "-wchecksqrt-est.txt";
//          ofstream f2(fn.c_str(),ios::out);
//          f2 << sat << " " << wchecksqrt << endl;
//          f2.close();  
//        }
//        
//        if (wcheck > wlimit)
//        {
//           sat = sat+0.5;
//           if (verbose > 1) cout << "   - Weight check " << wcheck << " > "<< wlimit  << " increasing sat: " << sat << endl;
//           md.first = firstbackup;
//           md.second = scaleinit;
//           r = resolution-rstart+1;
//           continue;
//         }
//     }
// 
//     if (r !=0) // adjust matrix to higher resolution level
//     {
//       for (int rr = 0; rr<3; rr++)
//       {
//         cmd.first[rr][3] = 2.0 * cmd.first[rr][3];
//       }
//     }
//     if (r == 2) // passed the test, adjust to highest
//     {
//       for (int rr = 0; rr<3; rr++)
//       {
//         cmd.first[rr][3] = 4.0 * cmd.first[rr][3];
//       }
//     }
//     md.first = cmd.first;
//     md.second = cmd.second;
//     if (verbose > 1)
//     {
//       cout << endl << " Matrix: " << endl;
//       cout << md.first << endl;
//       cout << " intens: " << md.second << endl;
//     }
//   } // resolution loop
//   
//   if (verbose > 1)
//   {
//     cout << endl << "   - current transform: " << endl;
//     cout << md.first << endl;
//     cout << "   - current iscale: " << md.second << endl;
//   }
//     
//   if (verbose > 0 )  cout << "   - final SAT: " << sat << " ( it: " << counter << " , weight check " << wcheck << " <= "<< wlimit << " )" <<  endl;
//     
//   iscale = iscaletmp;
// 
//   return sat;
// }

// void Registration::findSatMultiRes(const vnl_matrix < double > &mi, double scaleinit )
// // helper for findSaturation
// // basically the code from multiresoltuion
// // all kinds of stuff is initialized before (e.g. pyramid)
// {
//   // variables to store matrix m and scaling factor d:
//   pair < vnl_matrix_fixed < double, 4, 4> , double > cmd;
//   pair < vnl_matrix_fixed < double, 4, 4> , double > md(mi,scaleinit);
// 
//   // allow 2d case (depth == 1)
//   if ( gpS[0]->width < 16 || gpS[0]->height < 16 ||  (gpS[0]->depth < 16 && !gpS[0]->depth == 1))
//   {
//      ErrorExit(ERROR_BADFILE, "Input images must be larger than 16^3.\n") ;
//   }
//   int resolution = gpS.size();
//   assert(resolution >= 1); // otherwise we should have exited above
//   int rstart = 1;  // at least 16^3
//   
//   // stop if we get larger than 64^3 or if we reach highest resolution:
//   int stopres;
//   for (stopres = resolution-rstart; stopres>0; stopres--)
//   {
//      if (gpS[stopres]->width >= 64 || gpS[stopres]->height >= 64 || gpS[stopres]->depth >= 64) break;
//   }
//     
// //  bool iscaletmp = iscale;
// //  iscale = false; //disable intensity scaling on low resolutions
// 
//   for (int r = resolution-rstart;r>=stopres;r--)
//   {
//     if (verbose >1 ) 
//     { 
//       cout << endl << "Resolution: " << r << endl;
//       cout << " gpS ( " << gpS[r]->width <<" , " << gpS[r]->height << " , " << gpS[r]->depth << " )" << endl;
//       cout << " gpT ( " << gpT[r]->width <<" , " << gpT[r]->height << " , " << gpT[r]->depth << " )" << endl;
//     }
// 
// //    if (r==2) iscale = iscaletmp; // set iscale if set by user
// 
// 
//     // compute Registration
//     if (verbose >2 ) cout << "   - compute new iterative registration" << endl;
//     
//     int n = 3;
//     if (r==stopres) n = 1;
//     int vv = verbose;
//     if (verbose == 1) verbose = 0;
//     computeIterativeRegistration(n,0.05,gpS[r],gpT[r],md.first,md.second);
//     cmd.first = Mfinal;
//     cmd.second = iscalefinal;    
//     verbose = vv;
//     
//     if (verbose > 1)
//     {
//       cout << endl << " current Matrix: " << endl;
//       vnl_matlab_print(vcl_cout,cmd.first,"Tc",vnl_matlab_print_format_long);cout << endl;
//       cout << " intens: Ic = " << cmd.second << endl;
// 
//       // adjust to highest level for output only:
//       double tx = cmd.first[0][3];
//       double ty = cmd.first[1][3];
//       double tz = cmd.first[2][3];
//       for (int ll = r; ll > 0; ll--)
//       {
//         tx*=2;
//         ty*=2;
//         tz*=2;
//       }
//       cout << " equiv trans on highres: " << tx << " " << ty << " " << tz << endl;
//     }
//     
// //     if (r == stopres)
// //     {
// //        if (debug)
// //        {
// //          // write out wcheck
// //          string fn = getName() + "-wcheck-est.txt";
// //          ofstream f(fn.c_str(),ios::out);
// //          f << sat << " " << wcheck << endl;
// //          f.close();  
// //          string fn2 = getName() + "-wchecksqrt-est.txt";
// //          ofstream f2(fn.c_str(),ios::out);
// //          f2 << sat << " " << wchecksqrt << endl;
// //          f2.close();  
// //        }
// //        
// // //        if (wcheck > wlimit)
// // //        {
// // //           sat = sat+0.5;
// // //           if (verbose > 1) cout << "   - Weight check " << wcheck << " > "<< wlimit  << " increasing sat: " << sat << endl;
// // //           md.first = firstbackup;
// // //           md.second = scaleinit;
// // //           r = resolution-rstart+1;
// // //           continue;
// // //         }
// //     }
// 
//     if (r !=0) // adjust matrix to higher resolution level
//     {
//       for (int rr = 0; rr<3; rr++)
//       {
//         cmd.first[rr][3] = 2.0 * cmd.first[rr][3];
//       }
//     }
//     if (r == stopres) // passed the test, adjust to highest
//     {
//       for (; r>1; r--) // adjusted once already in the lines above
//       {
//         for (int rr = 0; rr<3; rr++)
//         {
//           cmd.first[rr][3] = 4.0 * cmd.first[rr][3];
//         }
//       }
//     }
//     md.first = cmd.first;
//     md.second = cmd.second;
//     if (verbose > 1)
//     {
//       cout << endl << " Matrix: " << endl;
//       vnl_matlab_print(vcl_cout,md.first,"T",vnl_matlab_print_format_long);cout << endl;
//       cout << " Intensity:  I = " << md.second << endl;
//     }
//   } // resolution loop
//   
// }
// 
// double Registration::findSaturation ( )
// {
//   if (verbose >0) cout << endl << endl << " Registration::findSaturation " << endl;
// //   if (!mriS) mriS = mri_source;
// //   if (!mriT) mriT = mri_target;
//   MRI * mriS = mri_source;
//   MRI * mriT = mri_target;
//   
//   vnl_matrix_fixed < double, 4, 4> m; m.set_identity();
// 
//   // variables to store matrix m and scaling factor d:
//   pair < vnl_matrix_fixed < double, 4, 4> , double > md(vnl_matrix_fixed < double, 4, 4> (),iscaleinit);
// 
//   if (!Minit.empty()) md.first = getMinitResampled();
//   else md.first = initializeTransform(mriS,mriT);
// 
// //   if (scaleinit != 1.0) md.second = scaleinit;
// //   else md.second = iscaleinit;
// 
//   int MINS = 16;
//   if (minsize > MINS) MINS = minsize; // use minsize, but at least 16
//   pair <int, int> limits = getGPLimits(mriS,mriT,MINS,maxsize);
//   if (gpS.size() ==0) gpS = buildGPLimits(mriS,limits);
//   if (gpT.size() ==0) gpT = buildGPLimits(mriT,limits);
//   assert(gpS.size() == gpT.size());
//   if ( gpS[0]->width < MINS || gpS[0]->height < MINS || (gpS[0]->depth < MINS && gpS[0]->depth != 1))
//   {
//      ErrorExit(ERROR_BADFILE, "Input images must be larger than 16^3.\n") ;
//   }
//   if ( gpT[0]->width < MINS || gpT[0]->height < MINS || (gpT[0]->depth < MINS && gpT[0]->depth != 1))
//   {
//      ErrorExit(ERROR_BADFILE, "Input images must be larger than 16^3.\n") ;
//   }  
//   
//   int resolution = gpS.size();
//   
//   assert(resolution >= 1); // otherwise we should have exited above
//   int rstart = 1;  // at least 16^3, last and most coarse image
//   
//   // stop if we get larger than 64^3
//   // should be as close as possible to 64
//   // the problem is that on other scales the wcheck limit gets meaningless
//   int stopres;
//   for (stopres = resolution-rstart; stopres>0; stopres--)
//   {
//      if (gpS[stopres]->width >= 64 || gpS[stopres]->height >= 64 || gpS[stopres]->depth >= 64) break;
//   }
//   
//   if ( gpS[stopres]->width < 32 || gpS[stopres]->height < 32 || gpS[stopres]->depth < 32)
//   {
//      cout << endl<< " ========================================================================" << endl;
//      cout << " WARNING: image might be too small (or ill shaped) for --satit to work." << endl;
//      cout << "          Try to manually specify --sat # if not satisfied with result! " << endl;
//      cout << " ========================================================================" << endl << endl;;
//   }
//   
//   cout << endl << "   - Max Resolution used: " << stopres << endl;
//   cout << "     -- gpS ( " << gpS[stopres]->width <<" , " << gpS[stopres]->height << " , " << gpS[stopres]->depth << " )" << endl;
//   cout << "     -- gpT ( " << gpT[stopres]->width <<" , " << gpT[stopres]->height << " , " << gpT[stopres]->depth << " )" << endl;
// 
// 
//   if (verbose >1 ) 
//   {
//     cout << "   - initial transform:\n" ;
//     vnl_matlab_print(vcl_cout,md.first,"Ti",vnl_matlab_print_format_long);cout << endl;
//     cout << "   - initial iscale:   Ii = " << md.second <<endl;
//   }
// 
//   // adjust md.first to current (lowest) resolution:
//   for (int r = 1; r<=resolution-rstart; r++)
//     for (int rr = 0;rr<3;rr++)
//       md.first[rr][3]  = 0.5 *  md.first[rr][3];
//       
//   vnl_matrix_fixed < double, 4, 4> firstbackup = md.first;
// 
//   if (verbose >1 ) 
//   {
//     cout << "   - initial adjusted:\n" ;
//     vnl_matlab_print(vcl_cout,md.first,"Tia",vnl_matlab_print_format_long);cout << endl;
//   }
// 
//   
//   // -------------------------------------------- RUN LOOP ----------------------------------
//   // 
//   cout << "   - running loop to estimate saturation parameter:\n" ;
//   double satdiff = 0.5; // stop if we get closer than this
//   double satmax  = 0;
//   double satmin  = 0;
//   double wmin    = -1; 
//   double wmax    = -1;
//   int counter    = 0;
//   while ( satmin == 0.0 || satmax == 0.0 || satmax-satmin > satdiff )
//   {
//      counter++;
//      if ( satmin == 0 && satmax == 0 ) sat = 16;
//      else if (satmin ==0) sat = 0.5 * satmax;
//      else if (satmax ==0) sat = 2 * satmin;
//      else sat = 0.5*(satmax+satmin);
//      if (verbose > 0)
//      {
//         if (counter > 1) cout << "         min sat: " << satmin << " ( " <<wmin <<" ), max sat: " << satmax << " ( " <<wmax <<" ), sat diff: " << satmax-satmin <<endl;
//         cout << "     -- Iteration: " << counter << "  trying sat: " << sat <<endl;
//      }
//      findSatMultiRes(md.first, md.second );
//      if (wcheck > wlimit)
//      {
//         satmin = sat;
//         wmin   = wcheck;
//      }
//      else
//      {
//        satmax = sat;
//        wmax = wcheck;
//      }
//      
//      // if sat low (sensitive) and still not many outliers
//      // quit, to prevent long (infinite) search
//      // e.g. if source and target are same image
//      if (sat < 6 && wcheck < 0.04)
//      {
//         satmax = sat;
//         satmin = sat;
//         wmax= wcheck;
//         wmin= wcheck;
//         break;
//      }
//   }
//   
//   // -------------------------------------------- SELECT FINAL ---------------------------------
//   // 
//   if (wmax <= wlimit )
//   {
//     sat = satmax;
//     wcheck = wmax;
//   }
//   else
//   {
//     assert (wmin <= wlimit);
//     sat = satmin;
//     wcheck = wmin;
//   }
//     
//   if (verbose > 0 )  cout << "   - final SAT: " << sat << " ( it: " << counter << " , weight check " << wcheck << " <= "<< wlimit << " )" <<  endl;
// 
//   if (debug)
//   {
//       // write out wcheck
//       string fn = getName() + "-wcheck-est.txt";
//       ofstream f(fn.c_str(),ios::out);
//       f << sat << " " << wcheck << endl;
//       f.close();  
//   }
// 
//   return sat;
// }

/**
 \param stopres stops on this resolution level (0 highest resolution ...)
 \param n number of max iterations on each resolution
 \param epsit epsilon to stop iterations
 */
void Registration::computeMultiresRegistration(int stopres, int n, double epsit)
{
  if (verbose > 0)
    cout << endl << endl << " Registration::computeMultiresRegistration "
        << endl;
//  if (verbose >0) cout << "   - Gaussian Pyramid " << endl;

//   if (!mriS) mriS = mri_source;
//   else
//   {
//     if (gpS.size() > 0) freeGaussianPyramid(gpS);
//     gpS = buildGaussianPyramid(mriS,16);
//   }
//   if (!mriT) mriT = mri_target;
//   else
//   {
//     if (gpT.size() > 0) freeGaussianPyramid(gpT);
//     gpT = buildGaussianPyramid(mriT,16);
//   }

  MRI * mriS = mri_source;
  MRI * mriT = mri_target;

  // variables to store matrix m and scaling factor d:
  pair<vnl_matrix_fixed<double, 4, 4>, double> cmd;
  pair<vnl_matrix_fixed<double, 4, 4>, double> md(
      vnl_matrix_fixed<double, 4, 4>(), iscaleinit);

  // check if mi (inital transform) is passed
  // if (!mi.empty()) md.first =mi;
  // else 
  if (!Minit.empty())
    md.first = getMinitResampled();
  else
    md.first = initializeTransform(mriS, mriT);

//  if (scaleinit != 1.0) md.second = scaleinit;
//  else md.second = iscaleinit;

  if (debug)
  {
    cout << " Debug: writing inital MOV resampled to DST ..." << endl;
    vnl_matlab_print(vcl_cout,md.first,"Minit",vnl_matlab_print_format_long);
    cout << endl;
    MRI * mri_tmp = MRIclone(mriT, NULL); // bring to same space as target (output after resampling)
    mri_tmp = MyMRI::MRIlinearTransform(mriS, mri_tmp, md.first);
    MyMRI::MRIvalscale(mri_tmp, mri_tmp, md.second);
    string fn = getName() + "-mapmov-init.mgz";
    cout << "    Filename: " << fn << endl;
    MRIwrite(mri_tmp, fn.c_str());
    MRIfree(&mri_tmp);
  }

  int MINS = 16;
  if (minsize > MINS)
    MINS = minsize; // use minsize, but at least 16

  ////freeGaussianPyramid(gpS);
  ////freeGaussianPyramid(gpT);
  //if (gpS.size() ==0) gpS = buildGaussianPyramid(mriS,MINS,maxsize);
  //if (gpT.size() ==0) gpT = buildGaussianPyramid(mriT,MINS,maxsize);
  pair<int, int> limits = getGPLimits(mriS, mriT, MINS, maxsize);
  if (gpS.size() == 0)
    gpS = buildGPLimits(mriS, limits);
  if (gpT.size() == 0)
    gpT = buildGPLimits(mriT, limits);
  assert(gpS.size() == gpT.size());
  if (gpT[0]->width < MINS || gpT[0]->height < MINS
      || (gpT[0]->depth < MINS && gpT[0]->depth != 1))
  {
    ErrorExit(ERROR_BADFILE, "Input images must be larger than 16^3.\n");
  }
  if (gpS[0]->width < 16 || gpS[0]->height < 16
      || (gpS[0]->depth < 16 && gpS[0]->depth != 1))
  {
    ErrorExit(ERROR_BADFILE, "Input images must be larger than 16^3.\n");
  }
  int resolution = gpS.size();
  //if ((int) gpT.size() < resolution) resolution = gpT.size();
  //gpS.resize(resolution);
  //gpT.resize(resolution);

  if (debug)
  {
    saveGaussianPyramid(gpS, "pyramidS");
    saveGaussianPyramid(gpT, "pyramidT");
  }

  if (verbose > 0)
  {
    cout << "   - initial transform:\n";
    vnl_matlab_print(vcl_cout,md.first,"Ti",vnl_matlab_print_format_long);
    cout << endl;
    cout << "   - initial iscale:  Ii =" << md.second << endl;
  }

  // adjust minit to current (lowest) resolution:
  //int rstart = 1;
  //for (int r = 1; r<=resolution-rstart; r++)
  for (int r = 1; r <= limits.second; r++)
    for (int rr = 0; rr < 3; rr++)
      md.first[rr][3] = 0.5 * md.first[rr][3];

  vnl_matrix_fixed<double, 4, 4> firstbackup = md.first;

  if (verbose > 1)
  {
    cout << "   - initial adjusted:\n";
    vnl_matlab_print(vcl_cout,md.first,"Tia",vnl_matlab_print_format_long);
    cout << endl;
  }

  if (debug)
  {
    cout << " Debug: writing inital MOV (lowres) resampled to DST (lowres)..."
        << endl;
    MRI * mri_tmp = MRIclone(gpT[resolution - 1], NULL); // bring to same space as target (output after resampling)
    mri_tmp = MyMRI::MRIlinearTransform(gpS[resolution - 1], mri_tmp, md.first);
    MyMRI::MRIvalscale(mri_tmp, mri_tmp, md.second);
    string fn = getName() + "-mapmov-low-init.mgz";
    MRIwrite(mri_tmp, fn.c_str());
    MRIfree(&mri_tmp);
    fn = getName() + "-mriT-low.mgz";
    MRIwrite(gpT[resolution - 1], fn.c_str());
  }

  bool iscaletmp = iscale;
//  iscale = false; //disable intensity scaling on low resolutions

  if (resolution - 1 < stopres)
  {
    ErrorExit(ERROR_BADFILE, "Input images have insufficient resoltuion.\n");
  }

  int m = n;
  converged = true;
  for (int r = resolution - 1; r >= stopres; r--)
  {
    //MRIwrite(gpS[r],"mriS-smooth.mgz");
    //MRIwrite(gpT[r],"mriT-smooth.mgz");

    if (verbose > 0)
    {
      cout << endl << "Resolution: " << r;
      cout << "  S( " << gpS[r]->width << " " << gpS[r]->height << " "
          << gpS[r]->depth << " )";
      cout << "  T( " << gpT[r]->width << " " << gpT[r]->height << " "
          << gpT[r]->depth << " )" << endl;
    }

    if (r <= 1)
      iscale = iscaletmp; // set iscale if set by user

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
    if (verbose > 1)
      cout << "- compute new iterative registration" << endl;
    //if (cmd.first) MatrixFree(&cmd.first);
    m = n;
    if (r == 0 && highit == 0)
    {
      cout << "  skipping (highit == 0) ..." << endl;
      break;
    }
    if (r == 0 && highit > 0)
      m = highit;
    computeIterativeRegistrationFull(m, epsit, gpS[r], gpT[r], md.first, md.second);
    cmd.first = Mfinal;
    cmd.second = iscalefinal;

    if (verbose > 1)
    {
      cout << endl << " current : Matrix: " << endl;
      //MatrixPrintFmt(stdout,"% 2.8f",cmd.first);
      vnl_matlab_print(vcl_cout,cmd.first,"Tc",vnl_matlab_print_format_long);
      cout << endl;
      cout << " intens:  Ic = " << cmd.second << endl;

      // adjust to highest level for output only:
      double tx = cmd.first[0][3];
      double ty = cmd.first[1][3];
      double tz = cmd.first[2][3];
      for (int ll = r; ll > 0; ll--)
      {
        tx *= 2;
        ty *= 2;
        tz *= 2;
      }
      cout << " equiv trans on highres: " << tx << " " << ty << " " << tz
          << endl;
    }

//     if (r == 2) // r==stopres (but stopres not determined here), better do this output in saturationestimation routine
//     {
//        // write out wcheck
//        if (debug)
//        {
//          string fn = getName() + "-wcheck.txt";
//          ofstream f(fn.c_str(),ios::out);
//          f << sat << " " << wcheck << endl;
//          f.close();  
//          string fn2 = getName() + "-wchecksqrt.txt";
//          ofstream f2(fn.c_str(),ios::out);
//          f2 << sat << " " << wchecksqrt << endl;
//          f2.close();  
//        }
//        if (verbose > 1)
//        {
//          cout << " wcheck : " << wcheck << "  wchecksqrt: " << wchecksqrt << endl;
//        }
//        
// //        if (wcheck > 0.3)
// //        {
// //           sat = sat+1;
// //           MatrixFree(&md.first);
// //           md.first = MatrixCopy(firstbackup,NULL);
// //           md.second = scaleinit;
// //           r = resolution-rstart+1;
// //           continue;
// //       }
//     }

    if (r != 0) // adjust matrix to higher resolution level
    {
      for (int rr = 0; rr < 3; rr++)
      {
        cmd.first[rr][3] = 2.0 * cmd.first[rr][3];
      }
    }
////       m  = MatrixMultiply(cmd.first,md.first,m);
////       MatrixCopy(m,md.first);
//    MatrixCopy(cmd.first,md.first);
    md.first = cmd.first;
    md.second = cmd.second;
    if (verbose > 1)
    {
      cout << endl << " Matrix: " << endl;
      vnl_matlab_print(vcl_cout,md.first,"Tcpost",vnl_matlab_print_format_long);
      cout << endl;
      cout << " intens: Icpost = " << md.second << endl;
    }
  } // resolution loop

  // adjust final transform to highest resolution (stopres + limits.first steps)
  if (stopres > 0)
    stopres -= 1; //because then we did one step inside of loop already
  // go through remaining stopres and also min pyramid steps
  for (int r = 0; r < stopres + limits.first; r++)
  {
    for (int rr = 0; rr < 3; rr++)
    {
      md.first[rr][3] = 2.0 * md.first[rr][3];
    }
  }

  if (verbose == 1)
  {
    cout << endl << "   - final transform: " << endl;
    vnl_matlab_print(vcl_cout,md.first,"Tf",vnl_matlab_print_format_long);
    cout << endl;
    cout << "   - final iscale:  If = " << md.second << endl;
  }

  if (!converged && verbose > 0)
  {
    cout << endl;
    cout << "**********************************************************"
        << endl;
    cout << "*" << endl;
    cout << "* WARNING: Registration did not converge in " << m << " steps! "
        << endl;
    cout << "*          Problem might be ill posed. " << endl;
    cout << "*          Please inspect output manually!" << endl;
    cout << "*" << endl;
    cout << "**********************************************************"
        << endl;
  }

  Mfinal = md.first;
  iscalefinal = md.second;

}

// // update this to vnl:
// double  Registration::computeSatEstimate (int reslevel, int n,double epsit, MRI * mriS, MRI* mriT, MATRIX* mi, double scaleinit )
// {
// //   cout << " Registration::computeSatEstimate " << endl;
// //   reslevel = 1;
// //   double PERCENT = 0.85;
// //   double EPS     = 0.01;
// // 
// //   pair < MATRIX*, double> fmd = computeMultiresRegistration(reslevel+1,n,epsit,mriS,mriT,mi,scaleinit);
// //   //     cout << endl << " Matrix: " << endl;
// //   //     MatrixPrintFmt(stdout,"% 2.8f",fmd.first);
// //   //     cout << " intens: " << fmd.second << endl;
// // 
// //   // Md is allready adjusted to current reslevel
// // 
// //   cout <<  endl << "Compute Sat estimate on Resolution " << reslevel << endl;
// //   cout << "   - warping source and target (sqrt)" << endl;
// // 
// //   MATRIX * mh = MyMatrix::MatrixSqrt(fmd.first);
// //   // do not just assume m = mh*mh, rather m = mh2 * mh
// //   // for transforming target we need mh2^-1 = mh * m^-1
// //   MATRIX * mii  = MatrixInverse(fmd.first,NULL);
// //   MATRIX *mhi = MatrixMultiply(mh,mii,NULL);
// // 
// // 
// //   MRI* mri_Swarp = MRIclone(gpS[reslevel],NULL);
// //   mri_Swarp = MRIlinearTransform(gpS[reslevel],mri_Swarp, mh);
// //   MRI* mri_Twarp = MRIclone(gpS[reslevel],NULL); // bring them to same space (just use src geometry)
// //   mri_Twarp = MRIlinearTransform(gpT[reslevel],mri_Twarp, mhi);
// // 
// //   // adjust intensity
// //   if (iscale)
// //   {
// //     cout << "   - adjusting intensity ( "<< fmd.second << " ) " << endl;
// //     //MRIvalscale(mri_Swarp,mri_Swarp,fmd.second);
// //     MyMRI::MRIvalscale(mri_Swarp,mri_Swarp,(1.0+fmd.second)*0.5);
// //     MyMRI::MRIvalscale(mri_Twarp,mri_Twarp,(1.0+ 1.0/fmd.second)*0.5);
// //   }
// //   vnl_matrix<double> A;
// //   vnl_vector<double> b;
// //   constructAb(mri_Swarp,mri_Twarp,A,b);
// //   pair < MATRIX*, MATRIX* > pwm(NULL,NULL);
// //   Regression<double> R(A,b);
// //   R.setVerbose(verbose);
// //   R.setFloatSvd(floatsvd);
// // 
// //   pair < double , double > interval(3,20);
// //   pair < double , double > wpercent;
// //   cout << "   - get weights left ( " << interval.first << " )" << endl;
// //   vnl_vector <double> p;
// //   p = R.getRobustEst(interval.first);
// //   cout << endl;
// //   wpercent.first = R.getLastWeightPercent();
// //   cout << "   - get weights right ( " << interval.second << " )" << endl;
// //   p = R.getRobustEst(interval.second);
// //   cout << endl;
// //   wpercent.second = R.getLastWeightPercent();
// // 
// //   assert (wpercent.first < PERCENT); // should be at sat ==1, otherwise one could try to go smaller?
// //   assert (wpercent.second > PERCENT); // should be at sat ==20, otherwise one could try to go higher
// //   double m, mp;
// //   int count = 0;
// //   while (interval.second - interval.first > 0.1 && wpercent.second - wpercent.first > EPS)
// //   {
// //     cout << endl << " Interval : w[ " << interval.first << " , " << interval.second << " ] = [ " << wpercent.first << " , " << wpercent.second << " ] : " << wpercent.second - wpercent.first<< endl;
// //     count++;
// // 
// //     // m = (PERCENT - wpercent.first) * (interval.second - interval.first) /  (wpercent.second - wpercent.first) + interval.first;
// //     m = 0.5*(interval.second + interval.first);
// //     cout << "   new test: " << m  << endl;
// //     p = R.getRobustEst(m);
// //     cout << endl;
// //     mp = R.getLastWeightPercent();
// //     cout << "    yields: " << mp << endl;
// // 
// //     if (mp < PERCENT)
// //     {
// //       interval.first = m;
// //       wpercent.first = mp;
// //     }
// //     else
// //     {
// //       interval.second = m;
// //       wpercent.second = mp;
// //     }
// //   }
// // 
// //   // cleanup
// //   if (mri_Twarp) MRIfree(&mri_Twarp);
// //   if (mri_Swarp) MRIfree(&mri_Swarp);
// //   MatrixFree(&mh);
// //   MatrixFree(&mii);
// //   MatrixFree(&mhi);
// // 
// //   cout << "Optimal sat ( at " << PERCENT << " ) : " << (interval.second + interval.first) *0.5 << endl;
// //   cout << " after " << count << " steps " << endl;
// //   return (interval.second + interval.first) *0.5;
// return -1;
// }
// 
// // double  Registration::computeSatEstimate (int reslevel, int n,double epsit, MRI * mriS, MRI* mriT, MATRIX* mi, double scaleinit )
// // {
// //    cout << " Registration::computeSatEstimate " << endl;
// //    if (!mriS) mriS = mri_source;
// //    if (!mriT) mriT = mri_target;
// //
// //    cout << "   - building Gaussian Pyramid " << endl;
// //    vector < MRI* > gpS = buildGaussianPyramid(mriS,100);
// //    vector < MRI* > gpT = buildGaussianPyramid(mriT,100);
// //    assert(gpS.size() == gpT.size());
// //    int resolution = gpS.size();
// //
// //
// //    int stoplevel = reslevel;
// //
// // // debug : save pyramid
// // //  for (uint i = 0;i<gpS.size();i++)
// // //  {
// // //   char fn[40];
// // //   sprintf(fn, "pyramid-%d.mgz", i+1);
// // //   MRIwrite(gpS[i],fn);
// // //  }
// //
// //
// //    vector < pair < double , double >  >  satzero(20);
// //    for (int i = 0 ; i<20;i++) satzero[i].first = double(i)+1.0;
// //
// //
// //
// // for (int s = 0 ; s < 20; s++)
// // {
// //
// //    cout << endl << " SATUTRATION : " << satzero[s].first << endl;
// //    sat = satzero[s].first;
// //
// //   MATRIX *m  = MatrixIdentity(4,NULL);
// //
// //   // variables to store matrix m and scaling factor d:
// //   pair < MATRIX* , double > cmd;
// //   pair < MATRIX* , double > md(NULL,scaleinit);
// //
// //   // check if mi (inital transform) is passed
// //   if (mi)
// //   {
// //       md.first = MatrixCopy(mi,NULL);
// //   }
// //   else if (Minit)
// //       md.first = MatrixCopy(Minit,NULL);
// //   else
// //   {
// //     //  md.first = MatrixIdentity(4,NULL);
// //     //   md.first = initialize_transform(mriS,mriT);
// //     // use voxtovox as init:
// //     md.first = MRIgetVoxelToVoxelXform(mriS,mriT) ;
// //   }
// //
// //   if (debug > 0)
// //   {
// //      cout << "   - initial transform:\n" ;
// //      MatrixPrintFmt(stdout,"% 2.8f",md.first);
// //   }
// //
// //   // adjust minit to current (lowest) resolution:
// //   int rstart = 2;
// //   for (int r = 1; r<=resolution-rstart; r++)
// //   for (int rr = 1;rr<=3;rr++)
// //      md.first->rptr[rr][4]  = 0.5 *  md.first->rptr[rr][4];
// //
// //    if(debug >0)
// //    {
// //       cout << "   - initial adjusted:\n" ;
// //       MatrixPrintFmt(stdout,"% 2.8f",md.first);
// //    }
// //
// //
// // //     MRI* falign = MRIlinearTransform(gpS[resolution-rstart], NULL,md.first);
// // //     MRIwrite(falign,"mriS-lowres-initaligned.mgz");
// // //     MRIfree(&falign);
// // //     MRIwrite(gpT[resolution-rstart],"mriT-lowres.mgz");
// //
// //
// // //  md.second = 1; //??
// // //  MRI* mri_Swarp   = NULL;
// // //  MRI* mri_Twarp   = NULL;
// //    vector < MATRIX* > matvec (resolution-rstart+1,NULL);
// //    vector < double >  intvec (resolution-rstart+1,1);
// // //  for (int r = resolution-rstart;r>=0;r--)
// //   for (int r = resolution-rstart;r>=stoplevel;r--)
// //   {
// //  //    MRIwrite(gpS[r],"mriS-smooth.mgz");
// //  //    MRIwrite(gpT[r],"mriT-smooth.mgz");
// //
// //       cout << endl << "Resolution: " << r << endl;
// //
// //
// //        // compute Registration
// //        cout << "   - compute new registration" << endl;
// //        if (cmd.first) MatrixFree(&cmd.first);
// //        cmd = computeIterativeRegistration(n,epsit,gpS[r],gpT[r],md.first,md.second);
// // //       cmd = computeIterativeRegSat(n,gpS[r],gpT[r],md.first,md.second);
// //        matvec[r] = MatrixCopy(cmd.first,matvec[r]);
// //        intvec[r] = cmd.second;
// //
// //     if(debug > 0)
// //     {
// //      cout << endl << " current : Matrix: " << endl;
// //      MatrixPrintFmt(stdout,"% 2.8f",cmd.first);
// //      cout << " intens: " << cmd.second << endl;
// //
// //      // adjust to highest level for output only:
// //      double tx = cmd.first->rptr[1][4];
// //      double ty = cmd.first->rptr[2][4];
// //      double tz = cmd.first->rptr[3][4];
// //      for (int ll = r; ll > 0; ll--)
// //      { tx *= 2; ty*=2; tz*=2;}
// //       cout << " equiv trans on highres: " << tx << " " << ty << " " << tz << endl;
// //
// //      // save resampled version on this level:
// //  //    MRI* salign = MRIlinearTransform(gpS[r], NULL,cmd.first);
// //  //    MRIwrite(salign,"mriS-lowres-aligned.mgz");
// // //     MRIwrite(gpT[r],"mriT-lowres.mgz");
// // //     MRIfree(&salign);
// //      }
// //
// //       if (r !=0) // adjust matrix to higher resolution level
// //       {
// //          for (int rr = 1; rr<=3; rr++)
// //   {
// //            // md.first->rptr[rr][4]  = 2.0 *  md.first->rptr[rr][4];
// //             cmd.first->rptr[rr][4] = 2.0 * cmd.first->rptr[rr][4];
// //   }
// //       }
// // //       m  = MatrixMultiply(cmd.first,md.first,m);
// // //       MatrixCopy(m,md.first);
// //        MatrixCopy(cmd.first,md.first);
// //        md.second = cmd.second;
// //      if (debug > 0)
// //      {
// //         cout << endl << " Matrix: " << endl;
// //         MatrixPrintFmt(stdout,"% 2.8f",md.first);
// //         cout << " intens: " << md.second << endl;
// //      }
// //   }
// //   // cleanup
// //   if (cmd.first) MatrixFree(&cmd.first);
// //   MatrixFree(&m);
// //   if (md.first) MatrixFree(&md.first);
// //
// //   //Mfinal = MatrixCopy(md.first, Mfinal);
// //   //iscalefinal = md.second;
// //
// //   satzero[s].second = zeroweights;
// //   cout << " sat: " << satzero[s].first << "  zeroweights: " << satzero[s].second << endl;
// //
// // }
// //    // plot diffs
// //    string fbase = name;
// //    int rf = fbase.rfind("/");
// //  if (rf != -1)
// //  {
// //             fbase = fbase.substr(rf+1,fbase.length());
// //  }
// //
// //    ofstream ofile((name+".plot").c_str(),ios::out);
// //    bool png = false;
// //    if (png) ofile << "set terminal png medium size 800,600" << endl;
// //    else ofile << "set terminal postscript eps color" << endl;
// //    if (png) ofile << "set output \""<< fbase <<".png\"" << endl;
// //    else ofile << "set output \""<< fbase <<".eps\"" << endl;
// //    ofile << "plot ";
// //    ofile << " \"-\" notitle with lines 1" << endl;
// //    for (int j = 0;j<(int)satzero.size(); j++)
// //    {
// //       ofile << satzero[j].first << " " << satzero[j].second << endl;
// //    }
// //    ofile << "e" << endl;
// //    ofile.close();
// //
// //   freeGaussianPyramid(gpS);
// //   freeGaussianPyramid(gpT);
// //
// //   return -1.0;
// // }
// //

void Registration::testRobust(const std::string& fname, int testno)
{

  cout << " testRobust " << fname << "  testno: " << testno << endl;

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

  MRI * mri = MRIread(fname.c_str());
  cout << " Read MRI " << fname.c_str() << " successfully!" << endl;

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

  vector<MRI*> gpS = buildGaussianPyramid(mri, 16, -1);
  //int level = gpS.size();
  int level = gpS.size() - 1;
//  MRIwrite(gpS[gpS.size()-level],"small.mgz");
  //cout << "sfasf" << endl;

  MATRIX* a = NULL, *ai = NULL;
  MRI* mriTs = NULL, *mriTt = NULL;

  double theta;
  double iscaleval = 1.0;
  switch (testno)
  {
  case 0: // identity
    cout << "Test " << testno << " : Identity" << endl;
    a = MatrixIdentity(4, a);
    ai = MatrixIdentity(4, ai);
    mriTs = MRIcopy(gpS[gpS.size() - level], NULL);
    mriTt = MRIcopy(gpS[gpS.size() - level], NULL);
    MRIwrite(mriTs, "idS.mgz");
    MRIwrite(mriTt, "idT.mgz");

    break;
  case 1: // translation
    cout << "Test " << testno << " : Translation" << endl;
    a = MatrixIdentity(4, a);
    a->rptr[1][4] = 2.2;
    a->rptr[2][4] = .3;
    a->rptr[3][4] = 1.4;
    ai = MatrixInverse(a, ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size() - level], mriTs, a,
        SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size() - level], mriTt, ai,
        SAMPLE_TRILINEAR);
    MRIwrite(mriTs, "transS.mgz");
    MRIwrite(mriTt, "transT.mgz");

    break;
  case 2:   // rotation
    cout << "Test " << testno << " : Rotation" << endl;
    theta = 0.08;
    a = MatrixAllocRotation(4, theta, Z_ROTATION);
    ai = MatrixInverse(a, ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size() - level], mriTs, a,
        SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size() - level], mriTt, ai,
        SAMPLE_TRILINEAR);
    MRIwrite(mriTs, "rotS.mgz");
    MRIwrite(mriTt, "rotT.mgz");
    break;
  case 3:   // intensity
  {
    cout << "Test " << testno << " : Intensity" << endl;
    a = MatrixIdentity(4, a);
    ai = MatrixIdentity(4, ai);
    iscaleval = 0.8;
    mriTs = MRIcopy(gpS[gpS.size() - level], NULL);
    mriTt = MyMRI::MRIvalscale(gpS[gpS.size() - level], NULL, iscaleval);
    MRI* tttt = mriTs;
    mriTs = mriTt;
    mriTt = tttt;
    MRIwrite(mriTs, "iscaleS.mgz");
    MRIwrite(mriTt, "iscaleT.mgz");
    break;
  }
  case 4:   // rotation and translation
    cout << "Test " << testno << " : Rotation and Translation" << endl;
    theta = 0.08;
    a = MatrixAllocRotation(4, theta, Z_ROTATION);
    a->rptr[1][4] = 2.2;
    a->rptr[2][4] = .3;
    a->rptr[3][4] = 1.4;
    ai = MatrixInverse(a, ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size() - level], mriTs, a,
        SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size() - level], mriTt, ai,
        SAMPLE_TRILINEAR);
    iscaleval = 0.8;
    mriTt = MyMRI::MRIvalscale(mriTt, NULL, iscaleval);
    MRIwrite(mriTs, "rottransS.mgz");
    MRIwrite(mriTt, "rottransT.mgz");
    break;
  case 5: // translation and junk
    cout << "Test " << testno << " : Translation and Noise" << endl;
    a = MatrixIdentity(4, a);
    a->rptr[1][4] = .2;
    a->rptr[2][4] = .3;
    a->rptr[3][4] = .4;
    ai = MatrixInverse(a, ai);

    mriTs = MRIlinearTransformInterp(gpS[gpS.size() - level], mriTs, a,
        SAMPLE_TRILINEAR);
    mriTt = MRIlinearTransformInterp(gpS[gpS.size() - level], mriTt, ai,
        SAMPLE_TRILINEAR);

    for (int dd = 0; dd < mriTs->depth / 3; dd++)
      for (int cc = 0; cc < mriTs->height / 3; cc++)
        for (int rr = 0; rr < mriTs->width / 3; rr++)
          MRIvox(mriTs, rr, cc, dd) = (rr)+(cc)+(dd)+1;

          MRIwrite(mriTs,"junktransS.mgz");
          MRIwrite(mriTt,"junktransT.mgz");
          break;
          case 6 : // rotation and junk
          cout << "Test " << testno << " : Rotation and Noise" << endl;
          theta = 0.02;
          a = MatrixAllocRotation(4,theta,Z_ROTATION);
          ai = MatrixInverse(a,ai);

          mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
          mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);

          //double nn = (mriTs->depth/3) +(mriTs->height/3) +(mriTs->width/3)+1;
          for (int dd = 0;dd<mriTs->depth/3;dd++)
          for (int cc = 0;cc<mriTs->height/3;cc++)
          for (int rr = 0;rr<mriTs->width/3;rr++)
          MRIvox(mriTs, rr, cc, dd) = ((rr)+(cc)+(dd)+1);

          MRIwrite(mriTs,"junkrotS.mgz");
          MRIwrite(mriTt,"junkrotT.mgz");

          break;
          case 7 :// scaling
          cout << "Test " << testno << " : Scaling" << endl;
          a = MatrixIdentity(4,NULL);
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
//   case 8 : //iscale only  (move to regRobust)
//   {
//     cout << "Test " << testno << " : Intensity" << endl;
//     a  = MatrixIdentity(4,a);
//     ai = MatrixIdentity(4,ai);
//     iscaleval = 0.8;
//     iscaleval = 0.2;
//     level = 3;
//     mriTs = MRIcopy(gpS[gpS.size()-level], NULL);
//     mriTt = MyMRI::MRIvalscale(gpS[gpS.size()-level], NULL, iscaleval);
//     MRIwrite(mriTs,"iscaleS.mgz");
//     MRIwrite(mriTt,"iscaleT.mgz");
//     
//     iscalefinal = 1;
//     MRI *snew = NULL;
//     MRI *tnew = NULL;
//     cout << " find: " << iscaleval << "  or : " << 1.0/iscaleval << endl;
//     for (int l = 0;l<10;l++)
//     { 
//       cout << " L " << l << " iscale: " << iscalefinal << endl;
//       snew = MyMRI::MRIvalscale(mriTs,snew, sqrt(iscalefinal));
//       tnew = MyMRI::MRIvalscale(mriTt,tnew, 1.0/sqrt(iscalefinal));      
//       estimateIScale(snew,tnew);
//       MRIfree(&snew);
//     }
//     cout << "final: " << iscalefinal << "  init: " << iscaleval << endl;
//     exit (0);
//     break;
//   }      
          case 9:// rotation and translation and iscale
          {
            cout << "Test " << testno << " : Rotation and Translation" << endl;
            theta = 0.08;
            a = MatrixAllocRotation(4,theta,Z_ROTATION);
            a->rptr[1][4] = 2.2;
            a->rptr[2][4] = .3;
            a->rptr[3][4] = 1.4;
            ai = MatrixInverse(a,ai);

            mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
            mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);
            iscaleval = 0.8;
            mriTt = MyMRI::MRIvalscale(mriTt, mriTt, iscaleval);
//     MRI* tttt = mriTs;
//     mriTs = mriTt;
//     mriTt = tttt;
            MRIwrite(mriTs,"rottransintS.mgz");
            MRIwrite(mriTt,"rottransintT.mgz");
            break;
          }
//   case 20: //error functions when rotating  (needs to be moved to RegRobust !!!)
//   {
//     int steps = 50;
//     double div = 4.0;
//     vector < double > theta(steps);
//     vector < double > err(steps);
//     vector < double > mls(steps);
//     vector < double > mls2(steps);
//     //level--;
//     RegistrationStep<float> RStep(*this);
//     for (int i=0; i<steps; i++)
//     {
//       // 0.. PI/div in 20 steps
//       // -PI/div ..0 is symmetric
//       theta[i] = M_PI * i / ((steps-1)*div);
// 
//       a  = MatrixAllocRotation(4,0.5*theta[i],Z_ROTATION);
//       a  = MatrixMultiply(a, extract_i_to_r(gpS[gpS.size()-level]), a);
//       a  = MatrixMultiply(extract_r_to_i(gpS[gpS.size()-level]) , a, a);
//       ai = MatrixInverse(a,ai);
// 
//       mriTs = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTs, a, SAMPLE_TRILINEAR);
//       mriTt = MRIlinearTransformInterp(gpS[gpS.size()-level],mriTt, ai, SAMPLE_TRILINEAR);
//       //MRIwrite(mriTs,"test20-s.mgz");
//       //MRIwrite(mriTt,"test20-t.mgz");
//       MatrixFree(&a);
//       MatrixFree(&ai);
//       ai = NULL;
// 
//       transonly = false;
//       costfun   = ROB;
//       rigid     = true;
//       iscale    = false;
// 
//       vnl_matrix< float > A;
//       vnl_vector< float > b;
//       RStep.constructAb(mriTs, mriTt,A,b);
//       pair < MATRIX*, MATRIX* > pwm(NULL,NULL);
//       Regression< float > R(A,b);
//       R.setVerbose(verbose);
// 
//       sat = 5;
// 
//       cout << "   - compute robust estimate ( sat "<<sat<<" )..." << flush;
//       R.getRobustEst(sat);
// 
//       err[i] = R.getLastError();
//       cout << "angle: " << theta[i] << "  error: " << err[i] << endl;
//       R.getLSEst();
//       mls[i] = R.getLastError();
//       cout << "angle: " << theta[i] << "  mls: " << mls[i] << endl;
//       MRI * mridiff = MRIalloc(mriTs->width, mriTs->height, mriTs->depth, MRI_FLOAT);
//       mridiff = MRIsubtract(mriTs,mriTt,mridiff);
//       double ddd = 0;
//       for (int d = 0;d<mriTs->depth;d++)
//         for (int h = 0;h<mriTs->height;h++)
//           for (int w = 0;w<mriTs->width;w++)
//             ddd += MRIgetVoxVal(mridiff,w,h,d,0) * MRIgetVoxVal(mridiff,w,h,d,0);
//             //ddd += MRIgetVoxVal(mridiff,w,h,d,1) * MRIgetVoxVal(mridiff,w,h,d,1);
//       mls2[i] = ddd;
//       cout << "angle: " << theta[i] << "  mls: " << mls2[i] << endl;
//       A.clear();
//       b.clear();
//     }
// 
//     ostringstream ss;
//     ss << "r-error-rot4-l" << level;
//     string fn = ss.str()+".plot";
//     ofstream f(fn.c_str(),ios::out);
// 
//     f << "set terminal postscript eps color" << endl;
//     f << "set title \"(Robust) error when rotating on level " << level <<"\"" << endl;
//     f << "set output \""<< ss.str() << ".eps\"" << endl;
//     f << "plot  \"-\" notitle with lines 1" << endl;
//     for (int i=0; i<steps; i++)
//     {
//       cout << theta[i] << " " << err[i] << endl;
//       f << theta[i] << " " << err[i] << endl;
//     }
//     f << "e" << endl;
// 
//     ostringstream ss2;
//     ss2 << "ls-error-rot4-l" << level;
//     string fn2 = ss2.str()+".plot";
//     ofstream f2(fn2.c_str(),ios::out);
// 
//     f2 << "set terminal postscript eps color" << endl;
//     f2 << "set title \"(LeastSquares) error when rotating on level " << level <<"\"" << endl;
//     f2 << "set output \""<< ss2.str() << ".eps\"" << endl;
//     f2 << "plot  \"-\" notitle with lines 1" << endl;
//     for (int i=0; i<steps; i++)
//     {
//       cout << theta[i] << " " << mls[i] << endl;
//       f2 << theta[i] << " " << mls[i] << endl;
//     }
//     f2 << "e" << endl;
// 
//     ostringstream ss3;
//     ss3 << "ils-error-rot4-l" << level;
//     string fn3 = ss3.str()+".plot";
//     ofstream f3(fn3.c_str(),ios::out);
// 
//     f3 << "set terminal postscript eps color" << endl;
//     f3 << "set title \"(IntensityLeastSquares) error when rotating on level " << level <<"\"" << endl;
//     f3 << "set output \""<< ss3.str() << ".eps\"" << endl;
//     f3 << "plot  \"-\" notitle with lines 1" << endl;
//     for (int i=0; i<steps; i++)
//     {
//       cout << theta[i] << " " << mls2[i] << endl;
//       f3 << theta[i] << " " << mls2[i] << endl;
//     }
//     f3 << "e" << endl;
// 
//     exit(0);
//     break;
//   }
          case 30:
          {
            cout << "Test " << testno << " : Save GP lowest" << endl;
            MRIwrite(gpS[gpS.size()-1],"gps-low.mgz");
            exit(0);
            break;
          }
          default:
          {
            cout << "Test " << testno << "unknown!" << endl;
            exit(1);
          }
        }

  cout << " Transformed , now registering ..." << endl;

//  rtype = 2;

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
  rigid = false;
  costfun = ROB;
  iscale = true;
  setVerbose(2);
//  sat = 5;

//  pair <MATRIX*, double> pwit    = computeIterativeRegistration(steps,mriTs,mriTt);
  pair<vnl_matrix_fixed<double, 4, 4>, double> pw;
  setSourceAndTarget(mriTs, mriTt);
  computeMultiresRegistration(0, 5, 0.01);
  pw.first = Mfinal;
  pw.second = iscalefinal;

  exit(0);
  costfun = LS;
//   pair <MATRIX*, double> pwls  = computeIterativeRegistration(steps,mriTs,mriTt);
//   pair <MATRIX*, double> pwls  = computeMultiresRegistration(mriTs,mriTt);
  iscale = true;
  costfun = ROB;
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
  MATRIX * aa = MatrixMultiply(ai, ai, NULL);
  MatrixPrintFmt(stdout, "% 2.8f", aa);
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
  //MatrixPrintFmt(stdout,"% 2.8f",pw.first);
  cout << pw.first << endl;
  //MatrixFree(&pw.first);

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

// can be adjusted to vnl_matrix, but maybe remove?
// bool Registration::warpSource( const string &fname, MATRIX* M, double is)
// {
//   assert(mri_source);
//   assert(mri_target);
// 
//   int nframes = mri_source->nframes;
//   mri_source->nframes = 1 ;
// 
//   MRI *mri_aligned = MRIclone(mri_target,NULL);
//   if (M)
//   {
// cout << "using M:" << endl ;
// MatrixPrint(stdout,M) ;
//     mri_aligned = MRIlinearTransform(mri_source, mri_aligned,M);
//   }
//   else if (Mfinal)  mri_aligned = MyMRI::MRIlinearTransform(mri_source, mri_aligned,Mfinal);
//   else
//   {
//     cerr << "warpSource error: no matrix set!" << endl;
//     MRIfree(&mri_aligned) ;
//     return false;
//   }
// 
// here also do scaling of intensity values
//   if (is == -1) is = iscalefinal;
//   if (is >0 && is != 1)
//     mri_aligned = MyMRI::MRIvalscale(mri_aligned, mri_aligned, is);
// 
//   mri_source->nframes = nframes ;
// 
//    sprintf(fname, "%s_after_final_alignment", parms.base_name) ;
//    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
//    sprintf(fname, "%s_target", parms.base_name) ;
//    MRIwriteImageViews(mri_dst, fname, IMAGE_SIZE) ;
// 
//   char cfname[STRLEN];
//   strcpy(cfname, fname.c_str()) ;
// 
//   MRIwrite(mri_aligned, cfname) ;
//   MRIfree(&mri_aligned) ;
// 
//   return true; // no error treatment so far
// }

// bool Registration::warpSource(MRI* orig, MRI* target, const string &fname, MATRIX* M, double is)
// // warps the mri orig to target
// {
//   assert(orig);
// 
//   int nframes = orig->nframes;
//   orig->nframes = 1 ;
// 
//   MATRIX* m_Lvox;
//   if (M) m_Lvox = MatrixCopy(M,NULL);
//   else m_Lvox   = MatrixCopy(Mfinal,NULL);
// 
//   /* convert it to RAS mm coordinates */
//   MATRIX* m_L = MRIvoxelXformToRasXform(mri_source, mri_target, m_Lvox, NULL) ;
// 
//   MRI *mri_aligned = MRIclone(target,NULL); //cp header and alloc space
//   mri_aligned = MRIapplyRASlinearTransform(orig, mri_aligned, m_L) ;
// 
//   // here also do scaling of intensity values
//   if (is == -1) is = iscalefinal;
//   if (is >0 && is != 1)
//     mri_aligned = MyMRI::MRIvalscale(mri_aligned, mri_aligned, is);
// 
//   orig->nframes = nframes ;
// 
// //    sprintf(fname, "%s_after_final_alignment", parms.base_name) ;
// //    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
// //    sprintf(fname, "%s_target", parms.base_name) ;
// //    MRIwriteImageViews(mri_dst, fname, IMAGE_SIZE) ;
// 
//   char cfname[STRLEN];
//   strcpy(cfname, fname.c_str()) ;
// 
//   MRIwrite(mri_aligned, cfname) ;
//   MRIfree(&mri_aligned) ;
// 
//   return true; // no error treatment so far
// }

// double Registration::estimateIScale(MRI *mriS, MRI *mriT)
// {
//  
//   if (verbose > 1) cout << "   - estimateIScale: " << endl;
// 
//   assert(mriT != NULL);
//   assert(mriS != NULL);
//   assert(mriS->width == mriT->width);
//   assert(mriS->height== mriT->height);
//   assert(mriS->depth == mriT->depth);
//   assert(mriS->type  == mriT->type);
//   //assert(mriS->width == mask->width);
//   //assert(mriS->height== mask->height);
//   //assert(mriS->depth == mask->depth);
//   //assert(mask->type == MRI_INT);
//   //MRIclear(mask);
// 
//   int z,y,x;
//   long int ss = mriS->width * mriS->height * mriS->depth;
//   if (mri_indexing) MRIfree(&mri_indexing);
//   if (ss > std::numeric_limits<int>::max())
//   {
//      if (verbose > 1) cout << "     -- using LONG for indexing ... " << flush;
//      mri_indexing = MRIalloc(mriS->width, mriS->height, mriS->depth,MRI_LONG);
//      if (mri_indexing == NULL) 
//         ErrorExit(ERROR_NO_MEMORY,"Registration::estimateIScale could not allocate memory for mri_indexing") ;
//      if (verbose > 1) cout << " done!" << endl;
//   }
//   else 
//   {
//      double mu = ((double)ss) * sizeof(int) / (1024.0 * 1024.0);
//      if (verbose > 1) cout << "     -- allocating " << mu << "Mb mem for indexing ... " << flush;
//      mri_indexing = MRIalloc(mriS->width, mriS->height, mriS->depth,MRI_INT);
//      if (mri_indexing == NULL) 
//         ErrorExit(ERROR_NO_MEMORY,"Registration::estimateIScale could not allocate memory for mri_indexing") ;
//      if (verbose > 1) cout << " done!" << endl;
//   }
// 
//   for (z = 0 ; z < mriS->depth ; z++)
//     for (x = 0 ; x < mriS->width ; x++)
//       for (y = 0 ; y < mriS->height ; y++)
//         MRILvox(mri_indexing, x, y, z) = 0;
// 
// 
//   bool dosubsample = false;
//   if (subsamplesize > 0)
//     dosubsample = (mriS->width > subsamplesize && mriS->height > subsamplesize && (mriS->depth > subsamplesize || mriS->depth == 1));
//   dosubsample = false; // needs to be fixed below!!! indeces are now randomized
// 
//   // we will need the blurred images (as float):
//   if (verbose > 1) cout << "     -- compute smoothie ... " << flush;
//   MRI *Sbl= MRIalloc(mriS->width, mriS->height, mriS->depth,MRI_FLOAT);
//   Sbl = MRIcopy(mriS,Sbl);
//   Sbl= MyMRI::getBlur(Sbl,Sbl);
//   MRI *Tbl= MRIalloc(mriT->width, mriT->height, mriT->depth,MRI_FLOAT);
//   Tbl = MRIcopy(mriT,Tbl);
//   Tbl= MyMRI::getBlur(Tbl,Tbl);
// 
//   if (verbose > 1) cout << " done!" << endl;
// 
//   if (dosubsample)
//   {
//     if (verbose > 1) cout << "     -- subsample ... "<< flush;
// 
//     MRI * Sblt = Sbl;
//     Sbl = MyMRI::subSample(Sblt);
//     MRIfree(&Sblt);
//     MRI * Tblt = Tbl;
//     Tbl = MyMRI::subSample(Tblt);
//     MRIfree(&Tblt);
// 
//     if (verbose > 1) cout << " done! " << endl;
//   }
// 
//   // compute 'counti': the number of rows needed (zero elements need to be removed)
//   int n = Sbl->width * Sbl->height * Sbl->depth;
//   if (verbose > 1) cout << "     -- size " << Sbl->width << " x " << Sbl->height << " x " << Sbl->depth << " = " << n << flush;
//   long int counti = 0;
//   double eps = 0.00001;
//   for (z = 0 ; z < Sbl->depth ; z++)
//     for (x = 0 ; x < Sbl->width ; x++)
//       for (y = 0 ; y < Sbl->height ; y++)
//       {
//         if (isnan(MRIFvox(Sbl, x, y, z)) ||isnan(MRIFvox(Tbl, x, y, z)))
//         {
//           if (verbose > 0) cout << " found a nan value!!!" << endl;
//           continue;
//         }
//         if (fabs(MRIFvox(Sbl, x, y, z)) < eps  && fabs(MRIFvox(Tbl, x, y, z)) < eps  )
//         {
//           //if (verbose > 0) cout << " found a zero element !!!" << endl;
//           continue;
//         }
//         counti++; // start with 1
//       }  
//   if (verbose > 1) cout << "  need only: " << counti << endl;
//   if (counti == 0)
//   {
//      cerr << endl;
//      cerr << " ERROR: All entries are zero! Images do not overlap (anymore?)." << endl;
//      cerr << "    This can have several reasons (i.e. different modalities, different "<< endl;
//      cerr << "    intensity scales, large non-linearities, too diff. voxel sizes ...)" << endl;
//      //cerr << "    Try calling with --noinit (if the original images are well aligned)" << endl;
//      cerr << "    Maybe use --transform <init.lta> with an approx. alignment" <<endl;
//      cerr << "    obtained from tkregister or another registration program." << endl;
//      cerr << "    Or do some prior intensity correction? " << endl;
//      cerr << endl;
//      exit(1);
//   }
//    
// 
//   // allocate the space for A and B
//   double abmu = ( (double)counti ) * sizeof(double) / (1024.0 * 1024.0);
//   if (verbose > 1) cout << "     -- allocating " << abmu << "Mb mem for A and b ... " << flush;
//   pair < vnl_matrix <double> , vnl_vector < double> > Ab( vnl_matrix <double>(counti,1), vnl_vector<double> (counti));
//   if (verbose > 1) cout << " done! " << endl;
// //      if (A == NULL || b == NULL) 
// //         ErrorExit(ERROR_NO_MEMORY,"Registration::estimateIScale could not allocate memory for A or b") ;
// 
//   if (verbose > 1) cout << "     -- size " << Sbl->width << " " << Sbl->height << " " << Sbl->depth << flush;
// 
//   long int count = 0;
//   int xp1,yp1,zp1;
//   for (z = 0 ; z < Sbl->depth ; z++)
//     for (x = 0 ; x < Sbl->width ; x++)
//       for (y = 0 ; y < Sbl->height ; y++)
//       {
//         if (isnan(MRIFvox(Sbl, x, y, z)) ||isnan(MRIFvox(Tbl, x, y, z))  )
//         {
//           if (verbose > 0) cout << " found a nan value!!!" << endl;
//           continue;
//         }
// 
//         if (dosubsample)
//         {
//           xp1 = 2*x;
//           yp1 = 2*y;
//           zp1 = 2*z;
//         }
//         else
//         {
//           xp1 = x;
//           yp1 = y;
//           zp1 = z; // if not subsampled
//         }
//         assert(xp1 < mriS->width);
//         assert(yp1 < mriS->height);
//         assert(zp1 < mriS->depth);
// 
// 
//         if (fabs(MRIFvox(Sbl, x, y, z)) < eps  && fabs(MRIFvox(Tbl, x, y, z)) < eps )
//         {
//           //cout << " found a zero row!!!" << endl;
//           MRILvox(mri_indexing, xp1, yp1, zp1) = -1;
//           continue;
//         }
// 
//         //count++; // start with 1
// 
//         if (xp1 >= mriS->width || yp1 >= mriS->height || zp1 >= mriS->depth)
//         {
// 
//           cerr << " outside !!! " << xp1 << " " << yp1 << " " << zp1 << endl;
//           assert(1==2);
//         }
// 
//         MRILvox(mri_indexing, xp1, yp1, zp1) = count;
// 
//         //Ab.first[count][0]  = 0.5 / iscalefinal *( MRIFvox(Tbl,x,y,z) + MRIFvox(Sbl, x, y, z)); 
//         //Ab.first[count][0]  = MRIFvox(Sbl, x, y, z);
// 
//         // intensity model: R(s,IS,IT) = exp(-0.5 s) IT - exp(0.5 s) IS
//         //                  R'  = -0.5 ( exp(-0.5 s) IT + exp(0.5 s) IS)
//         Ab.first[count][0]  = 0.5 * (MRIFvox(Tbl,x,y,z) + MRIFvox(Sbl, x, y, z));
// 
//         Ab.second[count] = -(MRIFvox(Tbl, x, y, z) - MRIFvox(Sbl, x, y, z));
//         
//         count++; // start with 0
// 
//       }
//       
//   // free remaining MRI    
//   if (Sbl) MRIfree(&Sbl);
//   if (Tbl) MRIfree(&Tbl);
// 
//   assert(counti == count);
// 
//   Regression<double> R(Ab.first,Ab.second);
//   R.setVerbose(verbose);
//   R.setFloatSvd(true); // even for double, the svd can be float, better switch to float all toghether
//   
//   vnl_vector<double> p( R.getRobustEst());
// 
//   double is = p[0];
//   double s = log(iscalefinal);
//   s = s - is;
//   iscalefinal = exp(s);
//   cout << " ISCALE: " << iscalefinal << " returned: " << is  << endl;
//   
//   return iscalefinal;
// }

//static int counter = 0;

// void Registration::constructAb(MRI *mriS, MRI *mriT,vnl_matrix < double >& A,vnl_vector<double>&b)
// // similar to robust paper
// // (with symmetry and iscale)
// {
// 
//   if (verbose > 1) cout << "   - constructAb: " << endl;
// 
//   assert(mriT != NULL);
//   assert(mriS != NULL);
//   assert(mriS->width == mriT->width);
//   assert(mriS->height== mriT->height);
//   assert(mriS->depth == mriT->depth);
//   assert(mriS->type  == mriT->type);
//   //assert(mriS->width == mask->width);
//   //assert(mriS->height== mask->height);
//   //assert(mriS->depth == mask->depth);
//   //assert(mask->type == MRI_INT);
//   //MRIclear(mask);
// 
//   int z,y,x;
//   long int ss = mriS->width * mriS->height * mriS->depth;
//   if (mri_indexing) MRIfree(&mri_indexing);
//   if (ss > std::numeric_limits<int>::max())
//   {
//      if (verbose > 1) cout << "     -- using LONG for indexing ... " << flush;
//      mri_indexing = MRIalloc(mriS->width, mriS->height, mriS->depth,MRI_LONG);
//      if (mri_indexing == NULL) 
//         ErrorExit(ERROR_NO_MEMORY,"Registration::constructAB could not allocate memory for mri_indexing") ;
//      if (verbose > 1) cout << " done!" << endl;
//   }
//   else 
//   {
//      double mu = ((double)ss) * sizeof(int) / (1024.0 * 1024.0);
//      if (verbose > 1) cout << "     -- allocating " << mu << "Mb mem for indexing ... " << flush;
//      mri_indexing = MRIalloc(mriS->width, mriS->height, mriS->depth,MRI_INT);
//      if (mri_indexing == NULL) 
//         ErrorExit(ERROR_NO_MEMORY,"Registration::constructAB could not allocate memory for mri_indexing") ;
//      if (verbose > 1) cout << " done!" << endl;
//   }
// 
//   for (z = 0 ; z < mriS->depth ; z++)
//     for (x = 0 ; x < mriS->width ; x++)
//       for (y = 0 ; y < mriS->height ; y++)
//         MRILvox(mri_indexing, x, y, z) = -10;
// 
// 
//   bool dosubsample = false;
//   if (subsamplesize > 0)
//     dosubsample = (mriS->width > subsamplesize && mriS->height > subsamplesize && mriS->depth > subsamplesize);
// //   bool dosubsample = true;
// 
// //    // remove 0 voxels:
// //   for (z = 0 ; z < mriS->depth ; z++)
// //   for (y = 0 ; y < mriS->height ; y++)
// //   for (x = 0 ; x < mriS->width ; x++)
// //   {
// //      if (MRIgetVoxVal(mriS,x,y,z,0) <= 0)
// //            MRIsetVoxVal(mriS,x,y,z,0,std::numeric_limits<float>::quiet_NaN());
// //      if (MRIgetVoxVal(mriT,x,y,z,0) <= 0)
// //            MRIsetVoxVal(mriT,x,y,z,0,std::numeric_limits<float>::quiet_NaN());
// //
// //   }
// 
//   // we will need the derivatives
//   if (verbose > 1) cout << "     -- compute derivatives ... " << flush;
// #if 0
//   MRI *Sfx=NULL,*Sfy=NULL,*Sfz=NULL,*Sbl=NULL;
//   MRI *Tfx=NULL,*Tfy=NULL,*Tfz=NULL,*Tbl=NULL;
//   MyMRI::getPartials(mriS,Sfx,Sfy,Sfz,Sbl);
//   MyMRI::getPartials(mriT,Tfx,Tfy,Tfz,Tbl);
//   MRI * SmT  = MRIsubtract(Sbl,Tbl,NULL); //S-T = f2-f1 =  delta f from paper
//   MRI * fx1  = MRIadd(Sfx,Tfx,Tfx); // store at Tfx location, renamed to fx1
//   MRIfree(&Sfx);
//   MRIscalarMul(fx1,fx1,0.5);
//   MRI * fy1  = MRIadd(Sfy,Tfy,Tfy); // store at Tfy location, renamed to fy1
//   MRIfree(&Sfy);
//   MRIscalarMul(fy1,fy1,0.5);
//   MRI * fz1  = MRIadd(Sfz,Tfz,Tfz); // store at Tfz location, renamed to fz1
//   MRIfree(&Sfz);
//   MRIscalarMul(fz1,fz1,0.5);
//   MRI * ft1  = MRIadd(Sbl,Tbl,Tbl); // store at Tbl location, renamed to ft1
//   MRIfree(&Sbl);
//   MRIscalarMul(ft1,ft1,0.5);
// #else
//   MRI *SpTh = MRIalloc(mriS->width,mriS->height,mriS->depth,MRI_FLOAT);
//   SpTh = MRIadd(mriS,mriT,SpTh);
//   SpTh = MRIscalarMul(SpTh,SpTh,0.5);
//   MRI *fx1 = NULL, *fy1 = NULL, *fz1 = NULL, *ft1 = NULL;
//   MyMRI::getPartials(SpTh,fx1,fy1,fz1,ft1);
//   MRI * SmT = MRIalloc(mriS->width,mriS->height,mriS->depth,MRI_FLOAT);
//   SmT = MRIsubtract(mriS,mriT,SmT);
//   SmT = MyMRI::getBlur(SmT,SmT);
// #endif
//   if (verbose > 1) cout << " done!" << endl;
//   //MRIwrite(fx1,"fx.mgz");
//   //MRIwrite(fy1,"fy.mgz");
//   //MRIwrite(fz1,"fz.mgz");
//   //MRIwrite(ft1,"ft.mgz");
// 
//   // subsample if wanted
//   MRI * fx,* fy,* fz,* ft;
//   if (dosubsample)
//   {
//     if (verbose > 1) cout << "     -- subsample ... "<< flush;
// 
//     fx = MyMRI::subSample(fx1);
//     MRIfree(&fx1);
//     fy = MyMRI::subSample(fy1);
//     MRIfree(&fy1);
//     fz = MyMRI::subSample(fz1);
//     MRIfree(&fz1);
//     ft = MyMRI::subSample(ft1);
//     MRIfree(&ft1);
//     
//     MRI * SmTt = SmT;
//     SmT = MyMRI::subSample(SmTt);
//     MRIfree(&SmTt);
//     
//     if (verbose > 1) cout << " done! " << endl;
//   }
//   else //just rename
//   {
//     fx = fx1;
//     fy = fy1;
//     fz = fz1;
//     ft = ft1;
//   }
// 
// //cout << " size fx : " << fx->width << " , " << fx->height << " , " << fx->depth << endl;
// //cout << " size src: " << mriS->width << " , " << mriS->height << " , " << mriS->depth << endl;
// 
//   // compute 'counti': the number of rows needed (zero elements need to be removed)
//   int n = fx->width * fx->height * fx->depth;
//   if (verbose > 1) cout << "     -- size " << fx->width << " x " << fx->height << " x " << fx->depth << " = " << n << flush;
//   long int counti = 0;
//   double eps = 0.00001;
//   int fxd = fx->depth ;
//   int fxw = fx->width ;
//   int fxh = fx->height ;
//   int fxstart = 0;
//   for (z = fxstart ; z < fxd ; z++)
//     for (x = fxstart ; x < fxw ; x++)
//       for (y = fxstart ; y < fxh ; y++)
//       {
//         if (isnan(MRIFvox(fx, x, y, z)) ||isnan(MRIFvox(fy, x, y, z)) || isnan(MRIFvox(fz, x, y, z)) || isnan(MRIFvox(ft, x, y, z)) )
//         {
//           if (verbose > 0) cout << " found a nan value!!!" << endl;
//           continue;
//         }
//         if (fabs(MRIFvox(fx, x, y, z)) < eps  && fabs(MRIFvox(fy, x, y, z)) < eps &&  fabs(MRIFvox(fz, x, y, z)) < eps )
//         {
//           //if (verbose > 0) cout << " found a zero element !!!" << endl;
//           continue;
//         }
//         counti++; // start with 1
//       }  
//   if (verbose > 1) cout << "  need only: " << counti << endl;
//   if (counti == 0)
//   {
//      cerr << endl;
//      cerr << " ERROR: All entries are zero! Images do not overlap." << endl;
//      cerr << "    Try calling with --noinit (if the original images are well aligned)" << endl;
//      cerr << "    Or use --transform <init.lta> with an approximate alignment" <<endl;
//      cerr << "    obtained from tkregister or another registration program." << endl << endl;
//      exit(1);
//   }
//    
// 
//   // allocate the space for A and B
//   int pnum = 12;
//   if (transonly)  pnum = 3;
//   else if (rigid) pnum = 6;
//   if (iscale) pnum++;
// 
//   double amu = ((double)counti*(pnum+1)) * sizeof(double) / (1024.0 * 1024.0); // +1 =  rowpointer vector
//   double bmu = (double)counti * sizeof(double) / (1024.0 * 1024.0);
//   if (verbose > 1) cout << "     -- allocating " << amu + bmu<< "Mb mem for A and b ... " << flush;
//   bool OK = A.set_size(counti,pnum);
//   OK = OK && b.set_size(counti);
//   if ( !OK )
//   {
//      cout << endl;
//      ErrorExit(ERROR_NO_MEMORY,"Registration::constructAB could not allocate memory for A and b") ;
//   
//   }
// //   MATRIX* A = MatrixAlloc(counti,pnum,MATRIX_REAL);
// //   VECTOR* b = MatrixAlloc(counti,1,MATRIX_REAL);
// //   if (A == NULL || b == NULL) 
// //   {
// //     cout << endl;
// //     ErrorExit(ERROR_NO_MEMORY,"Registration::constructAB could not allocate memory for A and b") ;
// //   }
//   if (verbose > 1) cout << " done! " << endl;
//   double maxmu = 5* amu + 7 * bmu;
//   if (floatsvd) maxmu = amu + 3*bmu + 2*(amu+bmu);
//   if (verbose > 1 ) cout << "         (MAX usage in SVD will be > " << maxmu << "Mb mem + 6 MRI) " << endl;
//   if (maxmu > 3800)
//   {
//     cout << "     -- WARNING: mem usage large: " << maxmu <<"Mb mem + 6 MRI" << endl;
//     string fsvd;
//     if (!floatsvd) fsvd = "--floatsvd and/or ";
//     cout << "          Maybe use "<<fsvd<< "--subsample <int> " << endl;
//   }
// 
//   // Loop and construct A and b
//   int xp1,yp1,zp1;
//   long int count = 0;
//   for (z = fxstart ; z < fxd ; z++)
//     for (x = fxstart ; x < fxw ; x++)
//       for (y = fxstart ; y < fxh ; y++)
//       {
//         if (isnan(MRIFvox(fx, x, y, z)) ||isnan(MRIFvox(fy, x, y, z)) || isnan(MRIFvox(fz, x, y, z)) || isnan(MRIFvox(ft, x, y, z)) )
//         {
//           //if (verbose > 0) cout << " found a nan value!!!" << endl;
//           continue;
//         }
// 
//         if (dosubsample)
//         {
//           //xp1 = 2*x+2;
//           //yp1 = 2*y+2;
//           //zp1 = 2*z+2;
//           xp1 = 2*x;
//           yp1 = 2*y;
//           zp1 = 2*z;
//         }
//         else // if not subsampled, shift only due to 5 tab derivative above
//         {
//           //xp1 = x+2;
//           //yp1 = y+2;
//           //zp1 = z+2; 
//           xp1 = x;
//           yp1 = y;
//           zp1 = z; 
//         }
//         assert(xp1 < mriS->width);
//         assert(yp1 < mriS->height);
//         assert(zp1 < mriS->depth);
// 
// 
//         if (fabs(MRIFvox(fx, x, y, z)) < eps  && fabs(MRIFvox(fy, x, y, z)) < eps &&  fabs(MRIFvox(fz, x, y, z)) < eps )
//         {
//           //if (verbose > 0) cout << " found a zero element !!!" << endl;
//           MRILvox(mri_indexing, xp1, yp1, zp1) = -1;
//           continue;
//         }
// 
//        // count++; // start with 1
// 
//         if (xp1 >= mriS->width || yp1 >= mriS->height || zp1 >= mriS->depth)
//         {
// 
//           cerr << " outside !!! " << xp1 << " " << yp1 << " " << zp1 << endl;
//           assert(1==2);
//         }
// 
//         MRILvox(mri_indexing, xp1, yp1, zp1) = count;
// 
//         //cout << "x: " << x << " y: " << y << " z: " << z << " coutn: "<< count << endl;
//         //cout << " " << count << " mrifx: " << MRIFvox(mri_fx, x, y, z) << " mrifx int: " << (int)MRIvox(mri_fx,x,y,z) <<endl;
//         int dof = 0;
//         if (transonly)
//         {
// //           *MATRIX_RELT(A, count, 1) = MRIFvox(fx, x, y, z);
// //           *MATRIX_RELT(A, count, 2) = MRIFvox(fy, x, y, z);
// //           *MATRIX_RELT(A, count, 3) = MRIFvox(fz, x, y, z);
//           A[count][0] = MRIFvox(fx, x, y, z);
//           A[count][1] = MRIFvox(fy, x, y, z);
//           A[count][2] = MRIFvox(fz, x, y, z);
//           dof = 3;
//         }
//         else if (rigid)
//         {
// //           *MATRIX_RELT(A, count, 1) = MRIFvox(fx, x, y, z);
// //           *MATRIX_RELT(A, count, 2) = MRIFvox(fy, x, y, z);
// //           *MATRIX_RELT(A, count, 3) = MRIFvox(fz, x, y, z);
// //           *MATRIX_RELT(A, count, 4) = MRIFvox(fz, x, y, z)*yp1 - MRIFvox(fy, x, y, z)*zp1;
// //           *MATRIX_RELT(A, count, 5) = MRIFvox(fx, x, y, z)*zp1 - MRIFvox(fz, x, y, z)*xp1;
// //           *MATRIX_RELT(A, count, 6) = MRIFvox(fy, x, y, z)*xp1 - MRIFvox(fx, x, y, z)*yp1;
//           A[count][0] =  MRIFvox(fx, x, y, z);
//           A[count][1] =  MRIFvox(fy, x, y, z);
//           A[count][2] =  MRIFvox(fz, x, y, z);
//           A[count][3] =  (MRIFvox(fz, x, y, z)*yp1 - MRIFvox(fy, x, y, z)*zp1);
//           A[count][4] =  (MRIFvox(fx, x, y, z)*zp1 - MRIFvox(fz, x, y, z)*xp1);
//           A[count][5] =  (MRIFvox(fy, x, y, z)*xp1 - MRIFvox(fx, x, y, z)*yp1);
//           dof = 6;
//           
//         }
//         else // affine
//         {
// //           *MATRIX_RELT(A, count, 1)  = MRIFvox(fx, x, y, z)*xp1;
// //           *MATRIX_RELT(A, count, 2)  = MRIFvox(fx, x, y, z)*yp1;
// //           *MATRIX_RELT(A, count, 3)  = MRIFvox(fx, x, y, z)*zp1;
// //           *MATRIX_RELT(A, count, 4)  = MRIFvox(fx, x, y, z);
// //           *MATRIX_RELT(A, count, 5)  = MRIFvox(fy, x, y, z)*xp1;
// //           *MATRIX_RELT(A, count, 6)  = MRIFvox(fy, x, y, z)*yp1;
// //           *MATRIX_RELT(A, count, 7)  = MRIFvox(fy, x, y, z)*zp1;
// //           *MATRIX_RELT(A, count, 8)  = MRIFvox(fy, x, y, z);
// //           *MATRIX_RELT(A, count, 9)  = MRIFvox(fz, x, y, z)*xp1;
// //           *MATRIX_RELT(A, count, 10) = MRIFvox(fz, x, y, z)*yp1;
// //           *MATRIX_RELT(A, count, 11) = MRIFvox(fz, x, y, z)*zp1;
// //           *MATRIX_RELT(A, count, 12) = MRIFvox(fz, x, y, z);
//           A[count][0]  = MRIFvox(fx, x, y, z)*xp1;
//           A[count][1]  = MRIFvox(fx, x, y, z)*yp1;
//           A[count][2]  = MRIFvox(fx, x, y, z)*zp1;
//           A[count][3]  = MRIFvox(fx, x, y, z);
//           A[count][4]  = MRIFvox(fy, x, y, z)*xp1;
//           A[count][5]  = MRIFvox(fy, x, y, z)*yp1;
//           A[count][6]  = MRIFvox(fy, x, y, z)*zp1;
//           A[count][7]  = MRIFvox(fy, x, y, z);
//           A[count][8]  = MRIFvox(fz, x, y, z)*xp1;
//           A[count][9]  = MRIFvox(fz, x, y, z)*yp1;
//           A[count][10] = MRIFvox(fz, x, y, z)*zp1;
//           A[count][11] = MRIFvox(fz, x, y, z);
//           dof = 12;
//         }
// 
// ////          if (iscale) *MATRIX_RELT(A, count, 7) = MRIFvox(Sbl, x, y, z);
//  //         if (iscale) *MATRIX_RELT(A, count, dof+1) = MRIFvox(Tbl, x, y, z);
//  
//  // !! ISCALECHANGE
//         // if (iscale) *MATRIX_RELT(A, count, dof+1) =  (0.5 / iscalefinal) * ( MRIFvox(Tbl, x, y, z) + MRIFvox(Sbl,x,y,z));
// 
//         // if (iscale) A[count][dof] =  (0.5 / iscalefinal) * ( MRIFvox(Tbl, x, y, z) + MRIFvox(Sbl,x,y,z));
//          if (iscale) A[count][dof] =  MRIFvox(ft, x, y, z) / iscalefinal;
// //         if (iscale) A[count][dof]  = MRIFvox(Sbl,x,y,z);
//          
// //         *MATRIX_RELT(b, count, 1) = - MRIFvox(ft, x, y, z); // ft = T-S => -ft = S-T
// //         b[count] = - MRIFvox(ft, x, y, z); // ft was = T-S => -ft = S-T
//          b[count] =  MRIFvox(SmT, x, y, z); // S-T
// 
//         count++; // start with 0 above
// 
//       }
//   assert(counti == count);
//       
//   // free remaining MRI    
//   MRIfree(&fx);
//   MRIfree(&fy);
//   MRIfree(&fz);
//   MRIfree(&ft);
//   MRIfree(&SmT);
//   //if (Sbl) MRIfree(&Sbl);
//   //if (Tbl) MRIfree(&Tbl);
// 
//   // Setup return pair
//  //  pair <MATRIX*, VECTOR* > Ab(A,b);
//   
// //   if (counter == 1) exit(1);
// //   counter++;
// //   MatrixWriteTxt((name+"A.txt").c_str(),A);
// //   MatrixWriteTxt((name+"b.txt").c_str(),b);
// //  if (counter == 1) exit(1);
//   
// //   // adjust sizes
// //   pair <MATRIX*, VECTOR* > Ab(NULL,NULL);
// //   double abmu2 = ((double)count*(pnum+1)) * sizeof(float) / (1024.0 * 1024.0);
// //   if (verbose > 1) cout << "     -- allocating another " << abmu2 << "Mb mem for A and b ... " << flush;
// //   Ab.first  = MatrixAlloc(count,pnum,MATRIX_REAL);
// //   Ab.second = MatrixAlloc(count,1,MATRIX_REAL);
// //   if (Ab.first == NULL || Ab.second == NULL) 
// //   {
// //     cout << endl;
// //     ErrorExit(ERROR_NO_MEMORY,"Registration::constructAB could not allocate memory for Ab.first Ab.second") ;
// //   }
// //   if (verbose > 1) cout << " done! " << endl;  
// //   for (int rr = 1; rr<= count; rr++)
// //   {
// //     *MATRIX_RELT(Ab.second, rr, 1) = *MATRIX_RELT(b, rr, 1);
// //     for (int cc = 1; cc <= pnum; cc++)
// //     {
// //       *MATRIX_RELT(Ab.first, rr, cc) = *MATRIX_RELT(A, rr, cc);
// //       assert (!isnan(*MATRIX_RELT(Ab.first, rr, cc)));
// //     }
// //     assert (!isnan(*MATRIX_RELT(Ab.second, rr, 1)));
// //   }
// //   MatrixFree(&A);
// //   MatrixFree(&b);
// 
//   return;
// }

pair<MATRIX*, VECTOR*> Registration::constructAb2(MRI *mriS, MRI *mriT)
{
  if (verbose > 1)
    cout << " constructAb2 " << endl;
  assert(mriT != NULL);
  assert(mriS != NULL);
  assert(mriS->width == mriT->width);
  assert(mriS->height== mriT->height);
  assert(mriS->depth == mriT->depth);
  assert(mriS->type == mriT->type);
  //assert(mriS->width == mask->width);
  //assert(mriS->height== mask->height);
  //assert(mriS->depth == mask->depth);
  //assert(mask->type == MRI_INT);
  //MRIclear(mask);

  // set <= 0 to NaN
  int z, y, x;
  for (z = 0; z < mriS->depth; z++)
    for (y = 0; y < mriS->height; y++)
      for (x = 0; x < mriS->width; x++)
      {
        if (MRIgetVoxVal(mriS, x, y, z, 0) <= 0)
          MRIsetVoxVal(mriS, x, y, z, 0,
              std::numeric_limits<float>::quiet_NaN());
        if (MRIgetVoxVal(mriT, x, y, z, 0) <= 0)
          MRIsetVoxVal(mriT, x, y, z, 0,
              std::numeric_limits<float>::quiet_NaN());

      }

  // we will need the derivatives
  if (verbose > 1)
    cout << " compute derivatives ... " << flush;
  MRI *Sfx = NULL, *Sfy = NULL, *Sfz = NULL, *Sbl = NULL;
  MRI *Tfx = NULL, *Tfy = NULL, *Tfz = NULL, *Tbl = NULL;
  MyMRI::getPartials(mriS, Sfx, Sfy, Sfz, Sbl);
  MyMRI::getPartials(mriT, Tfx, Tfy, Tfz, Tbl);

  if (verbose > 1)
    cout << " done!" << endl;

  if (verbose > 1)
    cout << " Subsample ... " << flush;

  MRI * fx = MyMRI::subSample(Tfx);
  MRI * fy = MyMRI::subSample(Tfy);
  MRI * fz = MyMRI::subSample(Tfz);
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

  if (verbose > 1)
    cout << " done! " << endl;

  // allocate the space
  int n = fx->width * fx->height * fx->depth;
  int pnum = 12;
  if (rigid)
    pnum = 6;
  if (iscale)
    pnum++;

  pair<MATRIX*, VECTOR*> Ab;

  MATRIX* A = MatrixAlloc(n, pnum, MATRIX_REAL);
  VECTOR* b = MatrixAlloc(n, 1, MATRIX_REAL);

  int count = 0;
  int xp1, yp1, zp1;
  for (z = 0; z < fx->depth; z++)
    for (x = 0; x < fx->width; x++)
      for (y = 0; y < fx->height; y++)
      {
        if (isnan(MRIFvox(fx, x, y, z) )|| isnan(MRIFvox(fy, x, y, z)) || isnan(MRIFvox(fz, x, y, z)) || isnan(MRIFvox(ssb, x, y, z))|| isnan(MRIFvox(stb, x, y, z)))
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
          *MATRIX_RELT(A, count, 1) = MRIFvox(fx, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 2) = MRIFvox(fx, x, y, z)*yp1;
          *MATRIX_RELT(A, count, 3) = MRIFvox(fx, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 4) = MRIFvox(fx, x, y, z);
          *MATRIX_RELT(A, count, 5) = MRIFvox(fy, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 6) = MRIFvox(fy, x, y, z)*yp1;
          *MATRIX_RELT(A, count, 7) = MRIFvox(fy, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 8) = MRIFvox(fy, x, y, z);
          *MATRIX_RELT(A, count, 9) = MRIFvox(fz, x, y, z)*xp1;
          *MATRIX_RELT(A, count, 10) = MRIFvox(fz, x, y, z)*yp1;
          *MATRIX_RELT(A, count, 11) = MRIFvox(fz, x, y, z)*zp1;
          *MATRIX_RELT(A, count, 12) = MRIFvox(fz, x, y, z);
          //if (iscale) *MATRIX_RELT(A, count, 13) = MRIFvox(ssb, x, y, z);
          if (iscale) *MATRIX_RELT(A, count, 13) = MRIFvox(stb, x, y, z);
        }

        *MATRIX_RELT(b, count, 1) = (MRIFvox(ssb, x, y, z) - MRIFvox(stb, x, y, z));

      }

      // adjust sizes
  Ab.first = MatrixAlloc(count, pnum, MATRIX_REAL);
  Ab.second = MatrixAlloc(count, 1, MATRIX_REAL);
  for (int rr = 1; rr <= count; rr++)
  {
    *MATRIX_RELT(Ab.second, rr, 1) = *MATRIX_RELT(b, rr, 1);
    for (int cc = 1; cc <= pnum; cc++)
    {
      *MATRIX_RELT(Ab.first, rr, cc) = *MATRIX_RELT(A, rr, cc);
    }
  }

  if (verbose > 1)
    cout << " Considering: " << count << " non-zero voxels" << endl;

  MatrixFree(&A);
  MatrixFree(&b);
  MRIfree(&fx);
  MRIfree(&fy);
  MRIfree(&fz);
  MRIfree(&ssb);
  MRIfree(&stb);
  return Ab;

}

// not needed anywhere
// MRI * Registration::applyParams(MRI * mri_in, const vnl_vector <double> & p, MRI * mri_dst, bool inverse)
// {
// 
//   if (!mri_dst)
//   {
//     mri_dst = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, mri_in->type,mri_in->nframes) ;
//     MRIcopyHeader(mri_in, mri_dst) ;
//   }
// 
//   assert (mri_in->width == mri_dst->width);
//   assert (mri_in->height== mri_dst->height);
//   assert (mri_in->depth == mri_dst->depth);
//   assert (mri_in->type  == mri_dst->type);
//   assert (mri_in->nframes == mri_dst->nframes);
//   assert (mri_in->nframes == 1);
// 
// 
//   // depending on number of params etc, apply:
//   bool doiscale = (p.size() == 4 || p.size()==7 || p.size() == 13);
//   assert (doiscale == iscale);
//   
//   
//   // here we use global rtype:
//   pair < vnl_matrix < double >, double> Md = convertP2Md(p);
//   
//   if (!inverse)
//   {
//     MATRIX * tmpmat = MyMatrix::convertVNL2MATRIX(Md.first,NULL);
//     mri_dst = MRIlinearTransform(mri_in,mri_dst, tmpmat);
//     Md.second = p[p.size()-1];
//     MatrixFree(&tmpmat);
//   }
//   else
//   {
//     vnl_matrix < double > mi = vnl_matrix_inverse<double>(Md.first);
//     MATRIX * tmpmat = MyMatrix::convertVNL2MATRIX(mi,NULL);
//     mri_dst = MRIlinearTransform(mri_in,mri_dst, tmpmat);
//     MatrixFree(&tmpmat);
//     Md.second = 1.0/ p[p.size()-1];
//   }
//   
// 
//   // adjust intensity
//   if (doiscale)
//   {
//     //if (verbose >1) cout << "   - adjusting intensity ( "<< fmd.second << " ) " << endl;
//     MyMRI::MRIvalscale(mri_dst,mri_dst,Md.second);
//   }
// 
//   return mri_dst;
// 
// }

// (mr) moved to RegistrationStep and converted to vnl_matrix
// MATRIX* Registration::constructR(MATRIX* p)
// // Construct restriction matrix (to restrict the affine problem to less parameters)
// // if p->rows == 6 use only rigid
// // if p->rows == 7 use also intensity scale
// // if p->rows == 3 use only trans
// // if p->rows == 4 use only trans + intensity
// {
//   assert(p != NULL);
//   assert((p->rows == 6 || p->rows==7) && p->cols ==1);
// 
//   int adim = 12;
//   if (iscale)
//   {
//     assert(p->rows == 7 || p->rows ==4);
//     adim++;
//   }
//   MATRIX* R = MatrixAlloc(p->rows,adim,MATRIX_REAL);
//   MatrixClear(R);
// 
// 
//   // translation p1,p2,p3 map to m4,m8,m12
//   *MATRIX_RELT(R, 1,  4) = 1.0;
//   *MATRIX_RELT(R, 2,  8) = 1.0;
//   *MATRIX_RELT(R, 3, 12) = 1.0;
// 
//   // iscale (p7 -> m13)
//   if (p->rows ==7) *MATRIX_RELT(R, 7, 13) = 1.0;
//   if (p->rows ==4) *MATRIX_RELT(R, 4, 13) = 1.0;
// 
//   if (p->rows <=4) return R;
// 
//   // rotation derivatives (dm_i/dp_i)
//   double s4 = sin(*MATRIX_RELT(p, 4, 1));
//   double c4 = cos(*MATRIX_RELT(p, 4, 1));
//   double s5 = sin(*MATRIX_RELT(p, 5, 1));
//   double c5 = cos(*MATRIX_RELT(p, 5, 1));
//   double s6 = sin(*MATRIX_RELT(p, 6, 1));
//   double c6 = cos(*MATRIX_RELT(p, 6, 1));
// 
//   *MATRIX_RELT(R, 5,  1) = -s5*c6;
//   *MATRIX_RELT(R, 6,  1) = -c5*s6;
// 
//   *MATRIX_RELT(R, 5,  2) = -s5*s6;
//   *MATRIX_RELT(R, 6,  2) =  c5*c6;
// 
//   *MATRIX_RELT(R, 5,  3) = -c5;
// 
//   *MATRIX_RELT(R, 4,  5) =  c4*s5*c6+s4*s6;
//   *MATRIX_RELT(R, 5,  5) =  s4*c5*c6;
//   *MATRIX_RELT(R, 6,  5) = -s4*s5*s6-c4*c6;
// 
//   *MATRIX_RELT(R, 4,  6) =  c4*s5*s6-s4*c6;
//   *MATRIX_RELT(R, 5,  6) =  s4*c5*s6;
//   *MATRIX_RELT(R, 6,  6) =  s4*s5*c6-c4*s6;
// 
//   *MATRIX_RELT(R, 4,  7) =  c4*c5;
//   *MATRIX_RELT(R, 5,  7) = -s4*s5;
// 
//   *MATRIX_RELT(R, 4,  9) = -s4*s5*c6+c4*s6;
//   *MATRIX_RELT(R, 5,  9) =  c4*c5*c6;
//   *MATRIX_RELT(R, 6,  9) = -c4*s5*s6+s4*c6;
// 
//   *MATRIX_RELT(R, 4, 10) = -s4*s5*s6-c4*c6;
//   *MATRIX_RELT(R, 5, 10) =  c4*c5*s6;
//   *MATRIX_RELT(R, 6, 10) =  c4*s5*c6+s4*s6;
// 
//   *MATRIX_RELT(R, 4, 11) = -s4*c5;
//   *MATRIX_RELT(R, 5, 11) = -c4*s5;
// 
//   return R;
// }
// 

// MATRIX * Registration::rt2mat(MATRIX * r, MATRIX * t, MATRIX *outM)
// // converts rot vector (3x1) and translation vector (3x1)
// // into an affine matrix (homogeneous coord) 4x4
// // if global rtype ==1 r1,r2,r3 are as in robust paper (axis, and length is angle)
// // if global rtype ==2 then r1,r2,r3 are angles around x,y,z axis (order 1zrot,2yrot,3xrot)
// {
//   if (outM == NULL)
//     outM = MatrixAlloc(4, 4, MATRIX_REAL);
// 
//   assert(r->rows == 3 && r->cols == 1);
//   assert(t->rows == 3 && t->cols == 1);
// 
//   MATRIX *rmat;
// 
//   if (rtype == 2)
//   {
// //      MATRIX* rx = MatrixAllocRotation(3,*MATRIX_RELT(r, 1, 1),X_ROTATION);
// //      MATRIX* ry = MatrixAllocRotation(3,*MATRIX_RELT(r, 2, 1),Y_ROTATION);
// //      MATRIX* rz = MatrixAllocRotation(3,*MATRIX_RELT(r, 3, 1),Z_ROTATION);
// //      MATRIX* tmp= MatrixMultiply(rx,ry,NULL);
// //      rmat = MatrixMultiply(tmp,rz,NULL);
// //      MatrixFree(&rx); MatrixFree(&ry); MatrixFree(&rz); MatrixFree(&tmp);
// //      MatrixPrintFmt(stdout,"% 2.8f",rmat);
// //      cout << endl;
// 
// 
//     // first convert rotation to quaternion (clockwise)
//     Quaternion q;
//     q.importZYXAngles(-*MATRIX_RELT(r, 3, 1), -*MATRIX_RELT(r, 2, 1), -*MATRIX_RELT(r, 1, 1));
//     // then to rotation matrix
//     rmat = MyMatrix::getMatrix(q.getRotMatrix3d(),3);
//     //MatrixPrintFmt(stdout,"% 2.8f",rmat2);
// 
//   }
//   else if (rtype ==1)
//   {
// 
//     // first convert rotation to quaternion
//     Quaternion q;
//     q.importRotVec(*MATRIX_RELT(r, 1, 1),*MATRIX_RELT(r, 2, 1),*MATRIX_RELT(r, 3, 1));
//     // then to rotation matrix
//     rmat = MyMatrix::getMatrix(q.getRotMatrix3d(),3);
// 
//   }
//   else assert (1==2);
// 
//   int rr, cc;
//   for (rr=1;rr<=3;rr++)
//   {
//     for (cc=1;cc<=3;cc++) // copy rot-matrix
//       *MATRIX_RELT(outM, rr, cc) = *MATRIX_RELT(rmat, rr, cc);
// 
//     // copy translation into 4th column
//     *MATRIX_RELT(outM, rr, 4) = *MATRIX_RELT(t, rr, 1);
//     // set 4th row to zero
//     *MATRIX_RELT(outM, 4, rr) = 0.0;
//   }
//   //except 4,4
//   *MATRIX_RELT(outM, 4, 4) = 1.0;
// 
//   MatrixFree(&rmat);
// 
//   return outM;
// }
// 
// 
// MATRIX * Registration::p2mat(MATRIX * p6, MATRIX *outM)
// // converts trans-rot vector (6x1)
// // into an affine matrix (homogeneous coord) 4x4
// // if rtype ==2 , then p4,p5,p6 are angles around x,y,z axis
// // else rtype ==1 they are as in robust paper
// {
//   assert(p6->rows == 6 && p6->cols == 1);
//   MATRIX* t = MatrixAlloc(3, 1, MATRIX_REAL);
//   MATRIX* r = MatrixAlloc(3, 1, MATRIX_REAL);
//   for (int rr = 1;rr<=3;rr++)
//   {
//     *MATRIX_RELT(t, rr, 1) = *MATRIX_RELT(p6, rr, 1);
//     *MATRIX_RELT(r, rr, 1) = *MATRIX_RELT(p6, rr+3, 1);
//   }
// 
//   outM = rt2mat(r,t,outM);
//   MatrixFree(&r);
//   MatrixFree(&t);
//   return outM;
// }
// 
// pair < MATRIX*, double > Registration::convertP2MATRIXd(MATRIX* p)
// // rtype : use restriction (if 2) or rigid from robust paper
// // returns registration as 4x4 matrix M, and iscale
// {
// //   cout << " Registration::convertP2Md(MATRIX* p) (p->rows: " << p->rows << " )" << flush;
//   pair < MATRIX*, double> ret(NULL,1.0);
//   MATRIX* pt;
// 
//   if (p->rows == 4 ||p->rows == 7 || p->rows == 13) // iscale
//   {
//     //cout << " has intensity " << endl;
//     // cut off intensity scale
//     
// //    //ret.second = 1.0-*MATRIX_RELT(p, p->rows, 1);    
//     ret.second = 1.0/(1.0+*MATRIX_RELT(p, p->rows, 1));
// //    ret.second = *MATRIX_RELT(p, p->rows, 1);
//     
//     pt=  MatrixAlloc(p->rows -1, 1,MATRIX_REAL);
//     for (int rr = 1; rr<p->rows; rr++)
//       *MATRIX_RELT(pt, rr, 1) = *MATRIX_RELT(p, rr, 1);
//   }
//   else pt = MatrixCopy(p,NULL);
// 
//   if (pt->rows == 12)
//   {
//     ret.first = MyMatrix::aff2mat(pt,NULL);
//   } 
//   else if (pt->rows == 6)
//   {
//     ret.first = p2mat(pt,NULL);
//   }
//   else if (pt->rows ==3)
//   {
//     ret.first = MatrixIdentity(4,NULL);
//     *MATRIX_RELT(ret.first, 1, 4) = *MATRIX_RELT(pt, 1, 1);
//     *MATRIX_RELT(ret.first, 2, 4) = *MATRIX_RELT(pt, 2, 1);
//     *MATRIX_RELT(ret.first, 3, 4) = *MATRIX_RELT(pt, 3, 1);
//   }
//   else
//   {
//     cerr << " parameter neither 3,6 nor 12 : " << pt->rows <<" ??" << endl;
//     assert(1==2);
//   }
// 
//   MatrixFree(&pt);
// //   cout << " -- DONE " << endl;
//   return ret;
// }
// 

/**
 Both MRI are assumed to have conform voxels (same size).
 \param min no dimension should get smaller than min voxels, default 16
 \param max no dimension will be larger than max in both images (can happen in one)
 */
pair<int, int> Registration::getGPLimits(MRI *mriS, MRI *mriT, int min, int max)
{

  if (verbose > 0)
    cout << "   - Get Gaussian Pyramid Limits ( min size: " << min
        << " max size: " << max << " ) " << endl;

  int smallest = mriS->width;
  if (mriT->width < smallest)
    smallest = mriT->width;
  if (mriS->height < smallest)
    smallest = mriS->height;
  if (mriT->height < smallest)
    smallest = mriT->height;
  //bool is2d = true;
  if (mriS->depth != 1 || mriT->depth != 1)
  {
    if (mriS->depth < smallest)
      smallest = mriS->depth;
    if (mriT->depth < smallest)
      smallest = mriT->depth;
    //is2d=true;
  }

  if (smallest < min)
  {
    cout << endl;
    cout << " Input image is smaller than min dimension: " << smallest << " < "
        << min << endl;
    cout << " Input dimensions (after conformalizing) MOV: " << mriS->width
        << " " << mriS->height << " " << mriS->depth << endl;
    cout << " Input dimensions (after conformalizing) DST: " << mriT->width
        << " " << mriT->height << " " << mriT->depth << endl;
    cout << " Specify larger minimum ..." << endl;
    exit(1);
  }

  // how often can I subdivide smallest to be >= min
  int maxsubdiv, minsubdiv;
  int temp = smallest / 2;
  for (maxsubdiv = 0; temp > min; maxsubdiv++)
    temp = temp / 2; // floor

  // if max size irrelevant, done:
  if (max == -1)
    return pair<int, int>(0, maxsubdiv);

  // get smaller value in each dimension:
  int t1 = mriT->width;
  if (mriS->width < t1)
    t1 = mriS->width;
  int t2 = mriT->height;
  if (mriS->height < t2)
    t2 = mriS->height;
  int t3 = mriT->depth;
  if (mriS->depth < t3)
    t3 = mriS->depth;

  //cout << " t1: " << t1 << " t2: " << t2 << " t3: " << t3 << endl;

  // get largest dimension
  temp = t1;
  if (t2 > temp)
    temp = t2;
  if (t3 > temp)
    temp = t3;

  //cout << " largest of smallest : " << temp;
  for (minsubdiv = 0; minsubdiv < maxsubdiv; minsubdiv++)
  {
    if (temp < max)
      break; // done, smaller value in all dimensions below max
    temp = temp / 2;
  }
  //cout << " minsubdiv : "<< minsubdiv << endl;
  //cout << "   then    : " << temp;

  return pair<int, int>(minsubdiv, maxsubdiv);

}

/** Here limits is the min and max iterations (subdivision level),
 meaning:  start highest resolution after min steps
 don't do more than max steps.
 */
vector<MRI*> Registration::buildGPLimits(MRI * mri_in,
    std::pair<int, int> limits)
{
  if (verbose > 0)
    cout << "   - Build Gaussian Pyramid ( Limits min steps: " << limits.first
        << " max steps: " << limits.second << " ) " << endl;

  int n = limits.second - limits.first + 1;
  vector<MRI*> p(n);

  // smoothing kernel:
  MRI *mri_kernel;
  mri_kernel = MRIgaussian1d(1.08, 5);
  //mri_kernel = MRIalloc(5,1,1, MRI_FLOAT);
  MRIFvox(mri_kernel, 0, 0, 0) = 0.0625;
  MRIFvox(mri_kernel, 1, 0, 0) = 0.25;
  MRIFvox(mri_kernel, 2, 0, 0) = 0.375;
  MRIFvox(mri_kernel, 3, 0, 0) = 0.25;
  MRIFvox(mri_kernel, 4, 0, 0) = 0.0625;

  MRI * mri_tmp;
  if (verbose > 1)
    cout << "        dim: " << mri_in->width << " " << mri_in->height << " "
        << mri_in->depth << endl;

  // smooth high res:
  //p[0] = MRIconvolveGaussian(mri_in, NULL, mri_kernel);
  p[0] = MRIcopy(mri_in, NULL);
  mri_tmp = mri_in;

  // subsample until highest resolution is small enough
  int i;
  for (i = 0; i < limits.first; i++)
  {
    //subsample:
    mri_tmp = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel);
    MRIfree(&p[0]);
    //p[0] = MRIdownsample2(mri_tmp,NULL);
    p[0] = MRIdownsample2BSpline(mri_tmp, NULL);
    MRIfree(&mri_tmp);
    mri_tmp = p[0];
  }

  //cout << " w[0]: " << p[0]->width << endl;
  p[0]->outside_val = mri_in->outside_val;
  int j = 1;
  for (; i < limits.second; i++)
  {
    //subsample:
    mri_tmp = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel);
    //p[j] = MRIdownsample2(mri_tmp,NULL);
    p[j] = MRIdownsample2BSpline(mri_tmp, NULL);
    p[j]->outside_val = mri_in->outside_val;
    MRIfree(&mri_tmp);
    mri_tmp = p[j];
    j++;
    //cout << " w[" << i<<"]: " << p[i]->width << endl;
    //cout << "i: " << i << " j: " << j << endl;
    //cout << " rows " << mri_tmp->r_to_i__->rows << " cols " << mri_tmp->r_to_i__->cols <<  endl;

  }

  assert(j==n);
  // check that all fields were filled

  MRIfree(&mri_kernel);
  return p;

}

/**
 \param mri_in input MRI
 \param min no dimension should get smaller than min voxels, default 16
 \param max no dimension will be larger than max
 */
vector<MRI*> Registration::buildGaussianPyramid(MRI * mri_in, int min, int max)
{

  if (verbose > 0)
    cout << "   - Build Gaussian Pyramid ( min: " << min << " max: " << max
        << " ) " << endl;

  // if max not passed allow pyramid to go up to highest resolution:
  if (max == -1)
    max = mri_in->width + mri_in->height + mri_in->depth;

  if (mri_in->width < min || mri_in->height < min || mri_in->depth < min)
  {
    cout << " Input image is smaller than min dimension: " << min << endl;
    cout << " Input dimensions (after conformalizing): " << mri_in->width << " "
        << mri_in->height << " " << mri_in->depth << endl;
    cout << " Specify larger minimum ..." << endl;
    exit(1);
  }

  int n = mri_in->depth; // choose n too large and adjust below

  vector<MRI*> p(n);
  MRI * mri_tmp;
// if (mri_in->type == MRI_UCHAR) cout << " MRI_UCHAR" << endl;
// else cout << " type: " << mri_in->type << endl;

  if (verbose > 1)
    cout << "        dim: " << mri_in->width << " " << mri_in->height << " "
        << mri_in->depth << endl;

  // smoothing kernel:
  MRI *mri_kernel;
  mri_kernel = MRIgaussian1d(1.08, 5);
  //mri_kernel = MRIalloc(5,1,1, MRI_FLOAT);
  MRIFvox(mri_kernel, 0, 0, 0) = 0.0625;
  MRIFvox(mri_kernel, 1, 0, 0) = 0.25;
  MRIFvox(mri_kernel, 2, 0, 0) = 0.375;
  MRIFvox(mri_kernel, 3, 0, 0) = 0.25;
  MRIFvox(mri_kernel, 4, 0, 0) = 0.0625;

  // smooth high res:
  p[0] = MRIconvolveGaussian(mri_in, NULL, mri_kernel);
  mri_tmp = mri_in;

  // if max passed subsample until highest resolution is below max
  if (max != -1)
  {

    // subsample if p[0] too large:
    while (p[0]->width > max || p[0]->height > max || p[0]->depth > max)
    {
      //subsample:
      mri_tmp = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel);
      MRIfree(&p[0]);
      //p[0] = MRIdownsample2(mri_tmp,NULL);
      p[0] = MyMRI::subSample(mri_tmp, NULL, true);
      //p[0] = MRIdownsample2BSpline(mri_tmp, NULL);
      MRIfree(&mri_tmp);
      mri_tmp = p[0];

    }
  }
  cout << " downsample2 " << endl;
  //cout << " w[0]: " << p[0]->width << endl;
  p[0]->outside_val = mri_in->outside_val;
  int i;
  for (i = 1; i < n; i++)
  {
    if (p[i - 1]->width < min || p[i - 1]->height < min
        || p[i - 1]->depth < min)
      break;
    //else subsample:
    mri_tmp = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel);
    //p[i] = MRIdownsample2(mri_tmp,NULL);
    p[i] = MyMRI::subSample(mri_tmp, NULL, true);
    p[i]->outside_val = mri_in->outside_val;
    //p[i] = MRIdownsample2BSpline(mri_tmp, NULL);
    MRIfree(&mri_tmp);
    mri_tmp = p[i];
    //cout << " w[" << i<<"]: " << p[i]->width << endl;
  }
  if (i < n)
    p.resize(i - 1);
  else
    assert(1==2);
  // should never get here as n should be large enough

  MRIfree(&mri_kernel);

  return p;
}

void Registration::freeGaussianPyramid(std::vector<MRI*>& p)
{
  for (uint i = 0; i < p.size(); i++)
    MRIfree(&p[i]);
  p.clear();
}

void Registration::saveGaussianPyramid(std::vector<MRI*>& p,
    const std::string & prefix)
{
  cout << "  Saving Pyramid " << prefix << endl;
  std::string name;
  for (uint i = 0; i < p.size(); i++)
  {
    std::stringstream out;
    out << i;
    name = prefix + "-r" + out.str() + ".mgz";
    MRIwrite(p[i], name.c_str());
  }
}

// ---------------------- Initial Transform using Moments -----------------------------

vnl_matrix_fixed<double, 4, 4> Registration::initializeTransform(MRI *mriS,
    MRI *mriT)
{
  if (verbose > 0)
    cout << "   - computing centroids \n";

// removed: do not trust RAS coordinates
//  it can happen that SRC is outside of Target frame
//  MATRIX* myinit = MRIgetVoxelToVoxelXform(mriS,mriT) ;
//  MRI * mri_tmp = MRIlinearTransform(mriS,NULL,myinit);
//  if (verbose > 1)
//  {
//      cout << endl <<"Vox2Vox (mov->dst): " << endl;
//      MatrixPrintFmt(stdout,"% 2.8f",myinit);
//  }

  //MATRIX* myinit = MatrixIdentity(4,NULL) ;
  vnl_matrix_fixed<double, 4, 4> myinit;
  myinit.set_identity();

  // find centroids:
  centroidS.clear();
  centroidT.clear();
  centroidS = CostFunctions::centroid(mriS);
  centroidT = CostFunctions::centroid(mriT);

  if (verbose > 1)
  {
    cout << "        Centroid S: " << centroidS[0] << ", " << centroidS[1]
        << ", " << centroidS[2] << endl;
    cout << "        Centroid T: " << centroidT[0] << ", " << centroidT[1]
        << ", " << centroidT[2] << endl;
  }

  if (!inittransform)
  {
    // default (new) initialize based on ras coordinates
    // later: allow option to use voxel identity for init?
    if (verbose > 0)
      cout << "   - initialize transform based on RAS\n";
    MATRIX *mv2rS, *mv2rT;
    mv2rS = extract_i_to_r(mriS);
    mv2rT = extract_i_to_r(mriT);
    // first convert to double, then compute vox2vox (using inverse)
    vnl_matrix_fixed<double, 4, 4> sv2r(MyMatrix::convertMATRIX2VNL(mv2rS));
    vnl_matrix_fixed<double, 4, 4> tv2r(MyMatrix::convertMATRIX2VNL(mv2rT));
    MatrixFree(&mv2rS);
    MatrixFree(&mv2rT);

    myinit = vnl_inverse(tv2r) * sv2r;
    //vnl_matlab_print(vcl_cerr,myinit,"myinit",vnl_matlab_print_format_long);cerr << endl;

    // OLD: inaccurate due to float inversion internally:
//     MATRIX* v2v = MRIgetVoxelToVoxelXform(mriS,mriT);
//     myinit = MyMatrix::convertMATRIX2VNL(v2v);
//     MatrixFree(&v2v);
//     vnl_matlab_print(vcl_cerr,myinit,"myinit",vnl_matlab_print_format_long);cerr << endl;

    return myinit;
  }

  if (verbose > 0)
    cout << "   - computing initial transform\n";

  if (initscaling)
  {
    if (symmetry)
    {
      cout
          << " ERROR: initscaling only meaningful with --nosym (else image sizes will agree)!\n";
      exit(1);
    }
    double rw = ((double) mriT->width) / ((double) mriS->width);
    double rh = ((double) mriT->height) / ((double) mriS->height);
    double rd = ((double) mriT->depth) / ((double) mriS->depth);
    //cout << " rw : " << rw << " rh : " << rh  << " rd :  " << rd << endl;
    //cout << " sw : " << mriS->width << " sh : " << mriS->height  << " sd :  " << mriS->depth << endl;
    double s;
    if (mriT->depth > 1 && mriS->depth > 1)
    {
      s = pow(rw * rh * rd, 1.0 / 3.0);
      myinit[0][0] = s;
      myinit[1][1] = s;
      myinit[2][2] = s;
    }
    else //assume 2d
    {
      s = sqrt(rw * rh);
      myinit[0][0] = s;
      myinit[1][1] = s;
    }
    //  if (verbose > 1)   
    cout << "     -- adjusting image sizes ( " << s << " ) \n";
    return myinit;
  }

  //bool initorient = false; // do not use orientation (can be off due to different cropping)
  //bool initorient = true;
  if (initorient)
  {
    if (verbose > 0)
      cout << "     -- trying to use orientation info (--initorient) \n";
    // find orientation:
    //MATRIX * evT = CostFunctions::orientation(mriT);
    //MATRIX * evS = CostFunctions::orientation(mriS);
    vnl_matrix_fixed<double, 3, 3> evT = CostFunctions::orientation(mriT);
    vnl_matrix_fixed<double, 3, 3> evS = CostFunctions::orientation(mriS);

    if (verbose > 1)
    {
      vnl_matlab_print(vcl_cerr,evT,"evT",vnl_matlab_print_format_long);
      vnl_matlab_print(vcl_cerr,evS,"evS",vnl_matlab_print_format_long);
    }

    // adjust orientation (flip) to target
    int fcount = 0;
    for (int c = 0; c < 3; c++)
    {
      double d = 0;
      for (int r = 0; r < 3; r++)
        d += evT[r][c] * evS[r][c];
      if (fabs(d) < 0.1)
      {
        cout << "       WARNING: Over 84 degree difference ( "
            << 360 * acos(fabs(d)) / M_PI << " ) in column " << c
            << ", skipping pre-orientation" << endl;
        initorient = false;
        break;
      }
      if (d < 0) // flip
      {
        if (verbose > 1)
          cout << "       --- flipping column " << c+1 <<" \n";
        fcount++;
        for (int r = 0; r < 3; r++)
          evS[r][c] = -evS[r][c];
      }

    }
    assert(fcount ==0 || fcount==2);
    // make sure det > 0 still

    if (debug)
    {
      cout << endl << " evecs S: " << endl;
      //MatrixPrintFmt(stdout,"% 2.8f",evS);
      cout << evS << endl;

      cout << endl << " evecs T: " << endl;
      //MatrixPrintFmt(stdout,"% 2.8f",evT);
      cout << evT << endl;
    }

    if (initorient) // still use orientation info for intial alignment
    {
      if (verbose > 0)
        cout << "     -- using orientation info\n";
      // make 4x4
      // MATRIX * evS2 = MatrixIdentity(4,NULL);
      // MATRIX * evT2 = MatrixIdentity(4,NULL);
      vnl_matrix_fixed<double, 4, 4> evS2;
      evS2.set_identity();
      vnl_matrix_fixed<double, 4, 4> evT2;
      evT2.set_identity();
      for (int c = 0; c < 3; c++)
        for (int r = 0; r < 3; r++)
        {
          evS2[r][c] = evS[r][c];
          evT2[r][c] = evT[r][c];
        }

      // first move Src centroid to origin
      myinit[0][3] -= centroidS[0];
      myinit[1][3] -= centroidS[1];
      myinit[2][3] -= centroidS[2];

      // Rotate evT * evS^-1
      //MATRIX* evS2i = MatrixInverse(evS2,NULL);
      vnl_matrix_fixed<double, 4, 4> evS2i = vnl_inverse(evS2);
      //assert(evS2i);
      //myinit = MatrixMultiply(evS2i,myinit,myinit);
      //myinit = MatrixMultiply(evT2,myinit,myinit);
      myinit = evT2 * evS2i * myinit;

      // finally move to T centroid
      myinit[0][3] += centroidT[0];
      myinit[1][3] += centroidT[1];
      myinit[2][3] += centroidT[2];

//       if (debug)
//       {
//         MATRIX * mtmp = MyMatrix::convertVNL2MATRIX(myinit,NULL);
//         MRI* mri_tmp = MRIlinearTransform(mriS,NULL,mtmp); // move to mymri using vnl matrix later
//         MatrixFree(&mtmp);
//         MRIwrite(mri_tmp,"init-align-rot.mgz");
//         MRIfree(&mri_tmp);        
//       }
      //MatrixFree(&evS2);
      //MatrixFree(&evS2i);
      //MatrixFree(&evT2);
    }
    //if (evT) MatrixFree(&evT);
    //if (evS) MatrixFree(&evS);

  }

  if (!initorient) // if orientation did not work
  {
    if (verbose > 0)
      cout << "     -- using translation info\n";
    myinit[0][3] += centroidT[0] - centroidS[0];
    myinit[1][3] += centroidT[1] - centroidS[1];
    myinit[2][3] += centroidT[2] - centroidS[2];
//     if (debug)
//     {
//       MATRIX * mtmp = MyMatrix::convertVNL2MATRIX(myinit,NULL);
//       MRI * mri_tmp = MRIlinearTransform(mriS,NULL,mtmp);
//       MatrixFree(&mtmp);
//       MRIwrite(mri_tmp,"initS-align-trans.mgz");
//       MRIfree(&mri_tmp);
//       MRIwrite(mriT,"initT-align-trans.mgz");
//     }
  }

  return myinit;
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
  MATRIX *m_scaling;
  float sx, sy, sz;
  MRI_REGION in_bbox, ref_bbox;

  m_scaling = MatrixIdentity(4, NULL);

  MRIboundingBox(mri_in, 60, &in_bbox);
  MRIboundingBox(mri_ref, 60, &ref_bbox);
  sx = (float) ref_bbox.dx / (float) in_bbox.dx;
  sy = (float) ref_bbox.dy / (float) in_bbox.dy;
  sz = (float) ref_bbox.dz / (float) in_bbox.dz;

  if (sx > MAX_DX)
    sx = MAX_DX;
  if (sx < MIN_DX)
    sx = MIN_DX;
  if (sy > MAX_DY)
    sy = MAX_DY;
  if (sy < MIN_DY)
    sy = MIN_DY;
  if (sz > MAX_DZ)
    sz = MAX_DZ;
  if (sz < MIN_DZ)
    sz = MIN_DZ;
//  if (Gdiag & DIAG_SHOW)
  fprintf(stderr, "initial scaling: (%2.2f, %2.2f, %2.2f) <-- "
      "(%d/%d,%d/%d,%d/%d)\n", sx, sy, sz, ref_bbox.dx, in_bbox.dx, ref_bbox.dy,
      in_bbox.dy, ref_bbox.dz, in_bbox.dz);
  *MATRIX_RELT(m_scaling, 1, 1) = sx;
  *MATRIX_RELT(m_scaling, 2, 2) = sy;
  *MATRIX_RELT(m_scaling, 3, 3) = sz;

#if 0
  *MATRIX_RELT(m_L, 1, 4) = dx;
  *MATRIX_RELT(m_L, 2, 4) = dy;
  *MATRIX_RELT(m_L, 3, 4) = dz;
#endif
  MatrixMultiply(m_scaling, m_L, m_L);
  return (NO_ERROR);
}

/** Important in some cases where input voxel sizes are negative
 based on weird RAS coordinates?
 Also adjusts resample matrices.
 */
bool Registration::flipInputs()
{
  bool flipped = false;
  int x = 1;
  int y = 2;
  int z = 3;
  if (mri_source->xsize < 0)
    x = -1;
  if (mri_source->ysize < 0)
    y = -2;
  if (mri_source->zsize < 0)
    z = -3;
  if (x < 0 || y < 0 || z < 0)
  {
    cout << "   Flipping axes in movable ( " << x << " " << y << " " << z
        << " ) ..." << endl;
    resample = true;
    flipped = true;
    //swap stuff:
    MRI* mri_temp = mri_source;
    mri_source = MRIreorder(mri_temp, NULL, x, y, z);

    MATRIX *v2v = MRIgetVoxelToVoxelXform(mri_source, mri_temp);
    vnl_matrix_fixed<double, 4, 4> Sreorderinv = MyMatrix::convertMATRIX2VNL(
        v2v);
    MatrixFree(&v2v);
    MRIfree(&mri_temp);

    if (debug)
    {
      string fn = getName() + "-mriS-flipped.mgz";
      cout << "   Writing flipped source as " << fn << endl;
      MRIwrite(mri_source, fn.c_str());
    }

    // adjust reslice matrix of source          
    // Rsrc points from resliced/reordered back to original input
    Rsrc = Rsrc * Sreorderinv;

  }

  if (mri_target->xsize < 0)
    x = -1;
  if (mri_target->ysize < 0)
    y = -2;
  if (mri_target->zsize < 0)
    z = -3;
  if (x < 0 || y < 0 || z < 0)
  {
    cout << "   Flipping axes in dest ( " << x << " " << y << " " << z
        << " ) ..." << endl;
    resample = true;
    flipped = true;
    //swap stuff:
    MRI* mri_temp = mri_target;
    mri_target = MRIreorder(mri_temp, NULL, x, y, z);

    MATRIX *v2v = MRIgetVoxelToVoxelXform(mri_target, mri_temp);
    vnl_matrix_fixed<double, 4, 4> Treorderinv = MyMatrix::convertMATRIX2VNL(
        v2v);
    MatrixFree(&v2v);
    MRIfree(&mri_temp);

    if (debug)
    {
      string fn = getName() + "-mriT-flipped.mgz";
      cout << "   Writing flipped dest as " << fn << endl;
      MRIwrite(mri_target, fn.c_str());
    }

    // adjust reslice matrix of source          
    // Rsrc points from resliced/reordered back to original input
    Rtrg = Rtrg * Treorderinv;

  }
  return flipped;
}

/** Potentially changes mri_Source and the resample matrix.
 This ensures that vox2vox rot is small and dimensions agree,
 important for Gaussian pyramid.
 */
bool Registration::reorientSource()
{
  //cout << "Registration::reorientSource ..." << endl;
  assert(mri_source);
  assert(mri_target);
  vnl_matrix_fixed<double, 4, 4> myinit;

  if (!Minit.empty())
  {
    //cout << " get Minit after resampling" << endl;
    myinit = getMinitResampled();
  }
  else
  {
    //cout << " Minit empty, get vox2vox" << endl;
    MATRIX* v2v = MRIgetVoxelToVoxelXform(mri_source, mri_target);
    myinit = MyMatrix::convertMATRIX2VNL(v2v);
    MatrixFree(&v2v);
  }
  //vnl_matlab_print(vcl_cout,myinit,"initv2v",vnl_matlab_print_format_long);cout << endl;
  //cout << " Determinant : " << vnl_determinant(myinit) << endl;

  // swap (and possibly invert) axis, so that rotation gets smaller
  int xd = 1;
  int yd = 1;
  int zd = 1;
  // determine max in each column:
  if (fabs(myinit[1][0]) > fabs(myinit[0][0]))
    xd = 2;
  if (fabs(myinit[2][0]) > fabs(myinit[0][0])
      && fabs(myinit[2][0]) > fabs(myinit[1][0]))
    xd = 3;
  if (fabs(myinit[1][1]) > fabs(myinit[0][1]))
    yd = 2;
  if (fabs(myinit[2][1]) > fabs(myinit[0][1])
      && fabs(myinit[2][1]) > fabs(myinit[1][1]))
    yd = 3;
  if (fabs(myinit[1][2]) > fabs(myinit[0][2]))
    zd = 2;
  if (fabs(myinit[2][2]) > fabs(myinit[0][2])
      && fabs(myinit[2][2]) > fabs(myinit[1][2]))
    zd = 3;
  // sign
  if (myinit[xd - 1][0] < 0.0)
    xd = -xd;
  if (myinit[yd - 1][1] < 0.0)
    yd = -yd;
  if (myinit[zd - 1][2] < 0.0)
    zd = -zd;

  //cout << " xd: " << xd << " yd: " << yd << " zd: " << zd << endl;

  if (xd == 1 && yd == 2 && zd == 3)
    return false; // nothing to swap

  if (abs(xd) * abs(yd) * abs(zd) != 6) // collision
  {
    cout << "  Registration::reorientSource: reorder not clear, solving conflict ..." << endl;
    vnl_matlab_print(vcl_cout,myinit,"v2v",vnl_matlab_print_format_long);
    cout << endl;
    
    // fix cases where 2 want to swap (and the third does not
    // want to stay in his spot, else there would be no conflict)
    // here force the third to stay:
    if (abs(xd) == 2 && abs(yd) ==1)
    {
      zd = 3;
      if (myinit[2][2] < 0.0) zd = -zd;      
    }
    else if (abs(xd) == 3 && abs(zd) ==1)
    {
      yd = 2;
      if (myinit[1][1] < 0.0) yd = -yd;    
    }    
    else if (abs(yd) == 3 && abs(zd) ==2)
    {
      xd = 1;
      if (myinit[0][0] < 0.0) xd = -xd;
    } 
    else
    {
      // there remain cases where (a) one wants to remain and others don't 
      // want to swap with each other, or (b) two want to remain
      // and third wants to swap:  keep all in their spot
      xd = 1;
      if (myinit[0][0] < 0.0) xd = -xd;
      yd = 2;
      if (myinit[1][1] < 0.0) yd = -yd;
      zd = 3;
      if (myinit[2][2] < 0.0) zd = -zd;
      
    }
        
  }

  cout << "   Reordering axes in mov to better fit dst... ( " << xd << " " << yd
      << " " << zd << " )" << endl;
  resample = true;

  //swap stuff:
  MRI* mri_temp = mri_source;
  mri_source = MRIreorder(mri_temp, NULL, xd, yd, zd);

  // get reordermatrix:
  MATRIX *v2v = MRIgetVoxelToVoxelXform(mri_source, mri_temp);
  vnl_matrix_fixed<double, 4, 4> Sreorderinv = MyMatrix::convertMATRIX2VNL(v2v);
  MatrixFree(&v2v);
  MRIfree(&mri_temp);

  if (debug)
  {
    string fn = getName() + "-mriS-reorder.mgz";
    cout << "   Writing reordered source as " << fn << endl;
    MRIwrite(mri_source, fn.c_str());
  }

//    vnl_matlab_print(vcl_cout,Sreorderinv,"reorder",vnl_matlab_print_format_long);cout << endl;
//    cout << " Determinant : " << vnl_determinant(Sreorderinv) << endl;

  // adjust reslice matrix of source          
  // Rsrc points from resliced/reordered back to original input
  Rsrc = Rsrc * Sreorderinv;


  // check if this fixed orientation problems (should have worked):
  double det = 0;
  vnl_matrix_fixed<double, 4, 4> v2vfull2;
  if (!Minit.empty())
  {
    v2vfull2 = getMinitResampled();
  }
  else
  {
    MATRIX *v2vfull = MRIgetVoxelToVoxelXform(mri_source, mri_target);
    v2vfull2 = MyMatrix::convertMATRIX2VNL(v2vfull);
    MatrixFree(&v2vfull); 
  }
  det = vnl_determinant(v2vfull2);
  cout << " Determinant after swap : " << det << endl;  
  if (det < 0 && symmetry) // this should not happen
  {  // cannot run sqrt later if det < 0
    cout << "ERROR: vox2vox det: " << det << " < 0"
         << endl;
     vnl_matlab_print(vcl_cout,v2vfull2,"v2v",vnl_matlab_print_format_long);cout << endl;
    cout << "       Symmetry does not work with flipped volumes."
         << endl;
    cout << "       Make sure volumes are in same voxel orientation." << endl;
    cout << "       Or run with --nosym flag." << endl;
    exit(1);
  }

  return true;

}

void Registration::setTransformation(bool is2d)
{
  if (trans)
    delete (trans);
  if (transonly)
  {
    if (is2d)
      trans = new Transform2dTranslate;
    else
      trans = new Transform3dTranslate;
  }
  else if (rigid)
  {
    if (is2d)
      trans = new Transform2dRigid;
    else
      trans = new Transform3dRigid;
  }
  else if (isoscale)
  {
    if (is2d)
      trans = new Transform2dIsoscale;
    else
      trans = new Transform3dIsoscale;
  }
  else if (affine)
  {
    if (is2d)
      trans = new Transform2dAffine;
    else
      trans = new Transform3dAffine;
  }
  else
  {
    trans = new Transform3dIdentity; // no geometric transform
    if (! iscale) // should never get here as --iscaleonly switches on --iscale
    {
       ErrorExit(ERROR_BADPARM,
          "Registration::setTransformation please specify either --iscale or a transformation model\n");

    }
  }
}

/** Both inputs need to be in the same voxel space.
 If not keeptype, ensure both are float.
 If keeptype, and types disagree, ensure both float (is that still true or is it uchar???)
 */
void Registration::setSourceAndTarget(MRI * s, MRI * t, bool keeptype)
{
  if (verbose > 0)
  {
    string kt = "FALSE";
    if (keeptype)
      kt = "TRUE";
    cout << "Registration::setSourceAndTarget(MRI s, MRI t, keeptype = " << kt
        << " )" << endl;
  }
  if (s->type != t->type)
    keeptype = false;

  if (verbose > 0)
  {
    cout << "   Type Source : " << s->type << "  Type Target : " << t->type;
    if (keeptype)
      cout << "  keeping type " << endl;
    else
      cout << "  ensure both FLOAT (" << MRI_FLOAT << ")" << endl;
  }

  // Init reslice matrices:
  vnl_matrix_fixed<double, 4, 4> Mid;
  Mid.set_identity();
  Rsrc = Mid;
  Rtrg = Mid;
  
  // init source and target mri:
  if (mri_source) MRIfree(&mri_source);
  mri_source = MRIcopy(s,NULL);
  if (mri_target) MRIfree(&mri_target);
  mri_target = MRIcopy(t,NULL);
  
  // reorder axis of srouce to match target orientation
  // flip and reorder axis of source based on RAS alignment or ixform:
  // this ensures that vox2vox rot is small and dimensions agree 
  // important for gaussian pyramid
  //cout << " Sin " << mri_source->width << " " << mri_source->height << " " << mri_source->depth << endl;
  //cout << " Tin " << mri_target->width << " " << mri_target->height << " " << mri_target->depth << endl;
  reorientSource();  
  //cout << " Sout " << mri_source->width << " " << mri_source->height << " " << mri_source->depth << endl;
  //cout << " Tout " << mri_target->width << " " << mri_target->height << " " << mri_target->depth << endl;

  
  // now we will make images isotropic

  // get smallest voxel length
  double mins = fabs(mri_source->xsize);
  if (fabs(mri_source->ysize) < mins)
    mins = fabs(mri_source->ysize);
  if (mri_source->depth > 1 && fabs(mri_source->zsize) < mins)
    mins = fabs(mri_source->zsize);
  double mint = fabs(mri_target->xsize);
  if (fabs(mri_target->ysize) < mint)
    mint = fabs(t->ysize);
  if (t->depth > 1 && fabs(mri_target->zsize) < mint)
    mint = fabs(mri_target->zsize);
  // select the larger of the smallest sides
  double isosize = mins;
  if (mint > mins)
    isosize = mint;
  vector<int> s_dim = MyMRI::findRightSize(mri_source, isosize, false);
  vector<int> t_dim = MyMRI::findRightSize(mri_target, isosize, false);

  //cout <<"   isosize: " << isosize << endl;
  //cout <<"   sdim : " << s_dim[0] << " " << s_dim[1] << " " << s_dim[2] << endl;
  //cout <<"   tdim : " << t_dim[0] << " " << t_dim[1] << " " << t_dim[2] << endl;

  if (symmetry) // adjust dimensions to match
  {
    for (uint i = 0; i < 3; i++)
    {
      if (s_dim[i] < t_dim[i])
        s_dim[i] = t_dim[i];
      else
        t_dim[i] = s_dim[i];
    }
    //cout <<"   sdim2 : " << s_dim[0] << " " << s_dim[1] << " " << s_dim[2] << endl;
    //cout <<"   tdim2 : " << t_dim[0] << " " << t_dim[1] << " " << t_dim[2] << endl;
  }
  

  if (verbose > 0)
  {
    cout << "   Mov: (" << mri_source->xsize << ", " << mri_source->ysize << ", " << mri_source->zsize
        << ")mm  and dim (" << mri_source->width << ", " << mri_source->height << ", " << mri_source->depth
        << ")" << endl;
    cout << "   Dst: (" << mri_target->xsize << ", " << mri_target->ysize << ", " << mri_target->zsize
        << ")mm  and dim (" << mri_target->width << ", " << mri_target->height << ", " << mri_target->depth
        << ")" << endl;

    cout << "   Asserting both images: " << isosize << "mm isotropic " << endl; //and (" << s_dim[0] << ", " << s_dim[1] << ", " << s_dim[2] <<") voxels" <<endl;
  }

  // source
  if (needReslice(mri_source, isosize,  s_dim[0], s_dim[1], s_dim[2], keeptype))
  {
    if (verbose > 0) cout << "    - reslicing Mov ..." << endl;
    pair<MRI*, vnl_matrix_fixed<double, 4, 4> > mm = makeIsotropic(mri_source, NULL,
      isosize, s_dim[0], s_dim[1], s_dim[2], keeptype);
    Rsrc = Rsrc * mm.second;
    assert(mm.first != NULL);
    if (mri_source)
      MRIfree(&mri_source);
    mri_source = mm.first;
  
    if (debug)
    {
      cout << "   Reslice Mov Matrix: " << endl << mm.second << endl;
      string n = name+string("-mriS-resample.mgz");
      cout << "   Writing resampled source as " << n << endl;
      MRIwrite(mri_source,n.c_str());
    }
  }
  else
  {
    if (verbose > 0)
      cout << "    - no Mov reslice necessary" << endl;
    //if (mri_source)
    //  MRIfree(&mri_source);
    //vnl_matrix_fixed<double, 4, 4> Rm;
    //Rm.set_identity();
    //Rsrc = Rm;
    //mri_source = MRIcopy(s,NULL);
  }

  // target
  if (needReslice(mri_target, isosize,  t_dim[0], t_dim[1], t_dim[2], keeptype))
  {
    if (verbose > 0) cout << "    - reslicing Dst ..." << endl;
    pair<MRI*, vnl_matrix_fixed<double, 4, 4> > mm = makeIsotropic(mri_target, NULL,
      isosize, t_dim[0], t_dim[1], t_dim[2], keeptype);
    Rtrg = Rtrg * mm.second;
    assert(mm.first != NULL);
    if (mri_target)
      MRIfree(&mri_target);
    mri_target = mm.first;
    if (debug)
    {
      cout << "   Reslice Dst Matrix: " << endl << mm.second << endl;
      string n = name + string("-mriT-resample.mgz");
      cout << "   Writing resampled target as " << n << endl;
      MRIwrite(mri_target, n.c_str());
    }
  }
  else
  {
    if (verbose > 0)
      cout << "    - no Dst reslice necessary" << endl;
    //if (mri_target)
    //  MRIfree(&mri_target);
    //vnl_matrix_fixed<double, 4, 4> Rm;
    //Rm.set_identity();
    //Rtrg = Rm;
    //mri_target = MRIcopy(t,NULL);
  }
  
  // set outside value to background:
  float bgvals = MyMRI::getBackground(mri_source);
  float bgvalt = MyMRI::getBackground(mri_target);
  //cout << "    - checkBackground: outside_val " << mri_source->outside_val << " src, " << mri_target->outside_val << " trg" << endl;
  //cout << "          suspected background val " << bgvals << " src, " << bgvalt << " trg" << endl;
  //MRIwrite(mri_source,"test.mgz");
  if (bgvals != mri_source->outside_val)
  {
    cout << "    - WARNING: src mri outside_val = " << mri_source->outside_val << "  but suspected background: " << bgvals << endl;     
    cout << "               If background is black or darkgray, ignore this. If white, pass --whitebgmov" << endl;
    //cout << "    - updating outside value..." << endl;
    //mri_source->outside_val = bgvals;
  }
  if (bgvalt != mri_target->outside_val)
  {
    cout << "    - WARNING: trg mri outside_val = " << mri_target->outside_val << "  but suspected background: " << bgvalt << endl;     
    cout << "               If background is black or darkgray, ignore this. If white, pass --whitebgdst" << endl;
    //cout << "    - updating outside value..." << endl;
    //mri_target->outside_val = bgvalt;
  }  

//  // flip and reorder axis of source based on RAS alignment or ixform:
//  // this ensures that vox2vox rot is small and dimensions agree 
//  // important for gaussian pyramid
//  reorientSource();  
//  cout << " S " << mri_source->width << " " << mri_source->height << " " << mri_source->depth << endl;
//  cout << " T " << mri_target->width << " " << mri_target->height << " " << mri_target->depth << endl;

  if (gpS.size() > 0)
    freeGaussianPyramid(gpS);
  centroidS.clear();
  if (gpT.size() > 0)
    freeGaussianPyramid(gpT);
  centroidT.clear();

  // initialize the correct registration type:
  bool is2d = (s->depth == 1 || t->depth == 1);
  setTransformation(is2d);

  if (verbose > 1)
    cout << " DONE setSourceAndTarget " << endl;
}

void Registration::setSource(MRI * s, bool conform, bool keeptype)
//local copy
{
  double vsize = -1.0;
  if (conform)
  {
    vsize = 1.0;
  }

  pair<MRI*, vnl_matrix_fixed<double, 4, 4> > mm = makeIsotropic(s, NULL, vsize,
      -1, -1, -1, keeptype);
  if (mri_source)
    MRIfree(&mri_source);
  mri_source = mm.first;
  Rsrc = mm.second;
  if (debug)
  {
    string n = name + string("-mriS-resample.mgz");
    MRIwrite(mri_source, n.c_str());
  }

  if (gpS.size() > 0)
    freeGaussianPyramid(gpS);
  centroidS.clear();

  // initialize the correct registration type:
  bool is2d = (s->depth == 1);
  setTransformation(is2d);
  
  //cout << "mri_source" << mri_source << endl;
}

void Registration::setTarget(MRI * t, bool conform, bool keeptype)
//local copy
{
  double vsize = -1.0;
  if (conform)
  {
    vsize = 1.0;
  }

  pair<MRI*, vnl_matrix_fixed<double, 4, 4> > mm = makeIsotropic(t, NULL, vsize,
      -1, -1, -1, keeptype);
  if (mri_target)
    MRIfree(&mri_target);
  mri_target = mm.first;
  Rtrg = mm.second;
  if (debug)
  {
    string n = name + string("-mriT-resample.mgz");
    MRIwrite(mri_target, n.c_str());
  }

  if (gpT.size() > 0)
    freeGaussianPyramid(gpT);
  centroidT.clear();
  //cout << "mri_target" << mri_target << endl;

  // initialize the correct registration type:
  bool is2d = (t->depth == 1);
  setTransformation(is2d);
}

/** To set name and nbase (base name without path)
 */
void Registration::setName(const std::string &n)
{
  name = n;
  nbase = n;
  int rf = nbase.rfind("/");
  if (rf != -1)
  {
    nbase = nbase.substr(rf + 1, nbase.length());
  }
}

bool Registration::needReslice(MRI *mri, double vsize, int xdim, int ydim,
    int zdim, bool keeptype)
{
  //cout << "Registration::needReslice(mri , "<< vsize << ", " << xdim << ", " << ydim << ", " << zdim << ", " << keeptype << ")" <<endl;

  int mw = mri->width;
  int mh = mri->height;
  int md = mri->depth;
  double mx = mri->xsize;
  double my = mri->ysize;
  double mz = mri->zsize;
  if (md == 1)
    mz = my; // adjust voxel size for tests below if 2D image

  //cout << mw << " "<< mh <<" "<<md << endl;
  //cout << mx << " "<< my <<" "<<mz << endl;

  // don't change type if keeptype or if already float:
  bool notypeconvert = (keeptype || mri->type == MRI_FLOAT);
//  if (notypeconvert && verbose > 1) cout << "     - no TYPE conversion necessary" << endl;
  // dont change voxel size if 
  // if already conform and no vsize specified
  bool novoxconvert = (vsize < 0 && mx == my && my == mz);
//  if (novoxconvert && verbose > 1) cout << "     - no vsize and allready conform "<< mri->xsize << endl;
  // if conform like vsize and no dims specified:
  double eps = 0.0001;
  bool conformvsize = (fabs(mx - vsize) < eps && fabs(my - vsize) < eps
      && fabs(mz - vsize) < eps);
//  if (conformvsize && verbose > 1) cout << "     - allready conform to "<< vsize << endl;
  novoxconvert = novoxconvert
      || (xdim < 0 && ydim < 0 && zdim < 0 && conformvsize);
//  if (novoxconvert && verbose > 1) cout << "     - no voxel conversion necessary (dimensions not passed)" << endl;
  // or if all values are specified and agree:
  bool dimagree = (xdim == mw && ydim == mh && zdim == md);
//  if (dimagree && verbose > 1) cout << "     - image dimensions agree" << endl;
  novoxconvert = novoxconvert || (conformvsize && dimagree);
//  if (novoxconvert && verbose > 1) cout << "     - no VOXEL conversion necessary" << endl;

  return !(novoxconvert && notypeconvert);
}

/** Makes the voxel size isotropic (and vox2ras map standard)
 */
std::pair<MRI*, vnl_matrix_fixed<double, 4, 4> > Registration::makeIsotropic(
    MRI *mri, MRI *out, double vsize, int xdim, int ydim, int zdim,
    bool keeptype)
{

  vnl_matrix_fixed<double, 4, 4> Rm;
  Rm.set_identity();

  if (!needReslice(mri, vsize, xdim, ydim, zdim, keeptype))
    return std::pair<MRI *, vnl_matrix_fixed<double, 4, 4> >(NULL, Rm);

  out = MRIcopy(mri, out);
//  if (out->depth == 1)
//    out->zsize = out->ysize; // for 2d images adjust z voxels size

  if (verbose > 1)
    cout << "     - will resample image or adjust type" << endl;
  // determine conform size
  double conform_size = vsize;
  if (conform_size < 0) // if not set use minimum:
  {
    conform_size = mri->xsize;
    if (mri->ysize < conform_size)
      conform_size = mri->ysize;
    if (mri->zsize < conform_size)
      conform_size = mri->zsize;
  }
  // get dimensions:
  vector<int> conform_dimensions = MyMRI::findRightSize(out, conform_size, false);
  if (verbose > 1)
    // will be cropped below, if larger than target dimensions
    cout<< "     - conform dim: " << conform_dimensions[0] << " "  << conform_dimensions[1] << " " << conform_dimensions[2] << endl;

  if (xdim > 0)
  {
    if (xdim < conform_dimensions[0])
    {
      cerr
          << "Registration::makeConform specified xdim too small to cover image"
          << endl;
      exit(1);
    }
    else
      conform_dimensions[0] = xdim;
  }
  if (ydim > 0)
  {
    if (ydim < conform_dimensions[1])
    {
      cerr
          << "Registration::makeConform specified ydim too small to cover image"
          << endl;
      exit(1);
    }
    else
      conform_dimensions[1] = ydim;
  }
  if (zdim > 0)
  {
    if (zdim < conform_dimensions[2])
    {
      cerr
          << "Registration::makeConform specified zdim too small to cover image"
          << endl;
      exit(1);
    }
    else
      conform_dimensions[2] = zdim;
  }

  // histogram based methods need uchar later
  // for rob at some point float performed better (needs further testing)
  MRI * temp;
//   if (costfun == ROB || costfun == LS)
//   {
  temp = MRIallocHeader(conform_dimensions[0], conform_dimensions[1],
      conform_dimensions[2], MRI_FLOAT, 1);
  MRIcopyHeader(mri, temp);
  temp->type = MRI_FLOAT;
//   }
//   else
//   {
//     temp = MRIallocHeader(conform_dimensions[0],
//                           conform_dimensions[1],
//                           conform_dimensions[2],
//                           MRI_UCHAR,
//                           1);
//     MRIcopyHeader(mri, temp);
//     temp->type   = MRI_UCHAR;
//   }

  temp->width = conform_dimensions[0];
  temp->height= conform_dimensions[1];
  temp->depth = conform_dimensions[2];
  temp->imnr0 = 1;
  temp->imnr1 = temp->depth;
  temp->thick = conform_size;
  temp->ps    = conform_size;
  temp->xsize = temp->ysize = temp->zsize = conform_size;  
  // if 2D image
  if (mri->depth == 1)
    temp->zsize = mri->zsize;
    
  temp->xstart = -conform_dimensions[0] / 2;
  temp->ystart = -conform_dimensions[1] / 2;
  temp->zstart = -conform_dimensions[2] / 2;
  temp->xend = conform_dimensions[0] / 2;
  temp->yend = conform_dimensions[1] / 2;
  temp->zend = conform_dimensions[2] / 2;
  // keep directional cosines from original (why rotate, just skrews up the dimenisons)
//   temp->x_r = -1.0;
//   temp->x_a =  0.0;
//   temp->x_s =  0.0;
//   temp->y_r =  0.0;
//   temp->y_a =  0.0;
//   temp->y_s = -1.0;
//   temp->z_r =  0.0;
//   temp->z_a =  1.0;
//   temp->z_s =  0.0;

  /* ----- change type if necessary ----- */
  if (mri->type != temp->type && !keeptype)
  {
    int no_scale_flag = FALSE;
    printf("       -- changing data type from %d to %d (noscale = %d)...\n",
        mri->type, temp->type, no_scale_flag);
    MRI * mri2 = MRISeqchangeType(out, temp->type, 0.0, 0.999, no_scale_flag);
    if (mri2 == NULL)
    {
      printf("ERROR: MRISeqchangeType\n");
      exit(1);
    }
    MRIfree(&out);
    out = mri2;
  }

  /* ----- reslice if necessary ----- */
//   if ((mri->xsize != temp->xsize ||
//        mri->ysize != temp->ysize ||
//        mri->zsize != temp->zsize ||
//        mri->width != temp->width ||
//        mri->height != temp->height ||
//        mri->depth != temp->depth ||
//        mri->x_r != temp->x_r ||
//        mri->x_a != temp->x_a ||
//        mri->x_s != temp->x_s ||
//        mri->y_r != temp->y_r ||
//        mri->y_a != temp->y_a ||
//        mri->y_s != temp->y_s ||
//        mri->z_r != temp->z_r ||
//        mri->z_a != temp->z_a ||
//        mri->z_s != temp->z_s ||
//        mri->c_r != temp->c_r ||
//        mri->c_a != temp->c_a ||
//        mri->c_s != temp->c_s) && fixvoxel > 0)
  double eps = 0.00001;
  if (fabs(mri->xsize - temp->xsize) > eps
      || fabs(mri->ysize - temp->ysize) > eps
      || fabs(mri->zsize - temp->zsize) > eps || mri->width != temp->width
      || mri->height != temp->height || mri->depth != temp->depth)
  {
    resample = true;
    // store resample matrix
    MATRIX * mtemp = MRIgetResampleMatrix(out, temp);
    Rm = MyMatrix::convertMATRIX2VNL(mtemp);
    MatrixFree(&mtemp);

    // NOT NEEDED if Cosines are kept from original
//     // determine rotation (do set target width, height, depth, as it might be not conform)
//     vector < int > rlookup(3,-1);
//     int product = 1;
//     for (uint c = 0 ; c<3; c++)
//     {
//       for (uint r = 0 ; r<3; r++)
//        if (fabs(Rm->rptr[r+1][c+1]) > fabs(Rm->rptr[((r+1)%3)+1][((c+1)%3)+1])  && 
//            fabs(Rm->rptr[r+1][c+1]) > fabs(Rm->rptr[((r+2)%3)+1][((c+2)%3)+1]) )
//            {
//              rlookup[c] = r; 
//              continue;
//            }
//        //cout << " c " << c << " from row " << rlookup[c] << endl; 
//       product *= (rlookup[c]+1);   // supposed to be 1*2*3 in the end to ensure we found each direction
//     }
//     
//     if (verbose > 1 || product != 6)
//     {
//        cout << " Resample Matrix : " << endl;
//        MatrixPrintFmt(stdout,"% 2.8f",Rm);
//        for (uint c = 0 ; c<3; c++)
//          cout << " resampled direction " << c << " from orig dir: " << rlookup[c] << endl; 
//     }
//     assert (product == 6);
//     // set target dimensions (according to possible geometry rotation)
//     temp->width  = conform_dimensions[rlookup[0]];
//     temp->height = conform_dimensions[rlookup[1]];
//     temp->depth  = conform_dimensions[rlookup[2]];
//     temp->imnr1  = temp->depth;
//     temp->xstart = -conform_dimensions[rlookup[0]]/2;
//     temp->ystart = -conform_dimensions[rlookup[1]]/2;
//     temp->zstart = -conform_dimensions[rlookup[2]]/2;
//     temp->xend   = conform_dimensions[rlookup[0]]/2;
//     temp->yend   = conform_dimensions[rlookup[1]]/2;
//     temp->zend   = conform_dimensions[rlookup[2]]/2;

    printf(
        "       -- Original : (%g, %g, %g)mm and (%d, %d, %d) voxels.\n",
        mri->xsize, mri->ysize, mri->zsize, mri->width, mri->height,
        mri->depth);
    printf(
        "       -- Resampled: (%g, %g, %g)mm and (%d, %d, %d) voxels.\n",
        temp->xsize, temp->ysize, temp->zsize, temp->width, temp->height,
        temp->depth);

    //int resample_type_val = SAMPLE_TRILINEAR;
    int resample_type_val = SAMPLE_CUBIC_BSPLINE;

    printf("       -- Reslicing using ");
    switch (resample_type_val)
    {
    case SAMPLE_TRILINEAR:
      printf("trilinear interpolation \n");
      break;
    case SAMPLE_NEAREST:
      printf("nearest \n");
      break;
    case SAMPLE_SINC:
      printf("sinc \n");
      break;
    case SAMPLE_CUBIC:
      printf("cubic \n");
      break;
    case SAMPLE_WEIGHTED:
      printf("weighted \n");
      break;
    case SAMPLE_CUBIC_BSPLINE:
      printf("cubic bspline \n");
      break;
    }
    MRI * mri2 = MRIresample(out, temp, resample_type_val);
//    printf("   Output   : (%g, %g, %g) mm size and (%d, %d, %d) voxels.\n",
//           mri2->xsize,mri2->ysize,mri2->zsize, mri2->width,mri2->height,mri2->depth);

    if (mri2 == NULL)
    {
      cerr << "makeIsotropic: MRIresample did not return MRI" << endl;
      exit(1);
    }

    MRIfree(&out);
    out = mri2;
  }
  else
  {
      cerr << "makeIsotropic WARNING: not different enough, won't reslice!" << endl;
  }

  MRIfree(&temp);
  return std::pair<MRI *, vnl_matrix_fixed<double, 4, 4> >(out, Rm);

}

// std::vector < double > Registration::getCentroidS()
// // map centroid back to original space
// {
//   vnl_vector_fixed < double, 4 > ncenter;
//   for (uint ii = 0; ii<3;ii++)
//      ncenter[ii] = centroidS[ii];
//   ncenter[3] = 1.0;
//   ncenter = Rsrc * ncenter;
//   vector < double > cs (3);
//   for (uint ii = 0; ii<3;ii++)
//      cs[ii] = ncenter[ii];
//   
//   return cs;
// }
// 
// std::vector < double > Registration::getCentroidT()
// // map centroid back to original space
// {
//   vnl_vector_fixed < double, 4 > ncenter;
//   for (uint ii = 0; ii<3;ii++)
//      ncenter[ii] = centroidT[ii];
//   ncenter[3] = 1.0;
//   ncenter = Rtrg * ncenter;
//   vector < double > ct (3);
//   for (uint ii = 0; ii<3;ii++)
//      ct[ii] = ncenter[ii];
//   return ct;
// }

/** Maps centroid back to original space
 */
vnl_vector_fixed<double, 4> Registration::getCentroidS()
{
  vnl_vector_fixed<double, 4> ncenter;
  for (uint ii = 0; ii < 3; ii++)
    ncenter[ii] = centroidS[ii];
  ncenter[3] = 1.0;
  ncenter = Rsrc * ncenter;
  return ncenter;
}

/** Maps centroid back to original space
 */
vnl_vector_fixed<double, 4> Registration::getCentroidT()
{
  vnl_vector_fixed<double, 4> ncenter;
  for (uint ii = 0; ii < 3; ii++)
    ncenter[ii] = centroidT[ii];
  ncenter[3] = 1.0;
  ncenter = Rtrg * ncenter;
  return ncenter;
}

/** Map source centroid to target space,
 only makes sense to call this after registration
 */
vnl_vector_fixed<double, 4> Registration::getCentroidSinT()
{
  vnl_vector_fixed<double, 4> ncenter = getCentroidS();
  vnl_matrix_fixed<double, 4, 4> v2v = getFinalVox2Vox();
  return v2v * ncenter;
}

/**
 <pre>
 (Mfinal)
 VOXS ------------> VOXT   (resampled)
 |   |                 |     |
 |   |                 |     |
 (Rsrc)|  RAS               RAS    |(Rtrg)
 |   |                 |     |
 V   |                 |     V
 VOXS ------------> VOXT   (original)
 (finalV2V)
 </pre>
 The matrix Mfinal is computed for the resampled (isotropic) images
 They have the same RAS coords as the original images.
 to compute finalV2V we use the resample matrices: Rsrc and Rtrg
 these were stored during resampling (conformalizing), they map
 from the  resampled/conform image to the original (vox2vox)
 */
vnl_matrix_fixed<double, 4, 4> Registration::getFinalVox2Vox()
{

  if (!resample)
    return Mfinal;

  cout << "Adjusting final transform due to initial resampling (voxel or size changes) ..." << endl;

//    // Make Mfinal to RAS2RAS
//    MATRIX * m = NULL;
//    m = MRIvoxelXformToRasXform (mri_source,mri_target, Mfinal, NULL) ;
// 
//    m = MRIrasXformToVoxelXform (src,trg,m,m);
// //   return m;
//    cout << " FINALVOX2VOX : " << endl;
//     MatrixPrintFmt(stdout,"% 2.8f",m);

  // since resampling did not change RAS coord, this should 
  // be the same as for the original data

  // the resample matrices are the inverse of the 
  // transform that happens when resampling
  vnl_matrix_fixed<double, 4, 4> m2;

  m2 = Rtrg * Mfinal * vnl_inverse(Rsrc);

  //cout << " OWN FINALVOX2VOX : " << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",m2);

  return m2;
}

std::pair<vnl_matrix_fixed<double, 4, 4>, vnl_matrix_fixed<double, 4, 4> > Registration::getHalfWayMaps()
{
  std::pair<vnl_matrix_fixed<double, 4, 4>, vnl_matrix_fixed<double, 4, 4> > md2w(
      mov2weights, dst2weights);

  if (!resample)
    return md2w;

//   if (Rsrc)
//   {   
//     md2w.first = MatrixInverse(Rsrc,md2w.first);
//     md2w.first = MatrixMultiply(mov2weights,md2w.first,md2w.first);
//   }
  md2w.first = mov2weights * vnl_inverse(Rsrc);

//   if (Rtrg)
//   {   
//     md2w.second = MatrixInverse(Rtrg,md2w.second);
//     md2w.second = MatrixMultiply(dst2weights,md2w.second,md2w.second);
//   }
  md2w.second = dst2weights * vnl_inverse(Rtrg);

  return md2w;
}

/** Returns the Minit matrix as vox2vox for the resampled images
 (it is originally the vox2vox of the input images before resampling),
 see also description of Registration::getFinalVox2Vox*/
vnl_matrix<double> Registration::getMinitResampled()
{
  if (Minit.empty())
    return Minit;
  if (!resample)
    return Minit;

  vnl_matrix<double> MIR(4, 4);

  MIR = vnl_inverse(Rtrg) * Minit * Rsrc;

  //cout << " Rscr: " << endl << Rsrc << endl;
  //cout << " Rtrg: " << endl << Rtrg << endl;

  return MIR;
}

/* Performs several checks (e.g. if sqrt produces a reflection, or projection). 
 If testrigid is true, check if sqrt produces reflection or scale, and if close to a rotation.
 Expects sqrt to be passed as M.
 */
bool Registration::checkSqrtM(const vnl_matrix_fixed<double, 4, 4> & M,
    bool testrigid)
{

  if (testrigid)
  {
    vnl_matrix<double> R(3, 3), S(3, 3), A(3, 3), I(3, 3);
    I.set_identity();
    M.extract(A);
    MyMatrix::PolarDecomposition(A, R, S);
    if (S[0][0] < 0.0 || S[1][1] < 0.0 || S[2][2] < 0.0)
      ErrorExit(ERROR_OUT_OF_BOUNDS,
          "Internal Error: Matrix Sqrt produced reflection.\n");
    double eps = 0.000001; // cannot be smaller due to scaling in ras2ras -> vox2vox conversion

    double fnorm1 = (S - I).frobenius_norm();
    if (fnorm1 > eps)
    {
      std::cerr << "Internal Error: " << std::endl;
      std::cerr << " Sqrt of Rotation should not scale ( fnorm(S-I) = "
          << fnorm1 << " )" << std::endl;
      std::cerr << " Debug Info: " << std::endl;
      vnl_matlab_print(vcl_cerr,A,"A",vnl_matlab_print_format_long);
      std::cerr << std::endl;
      vnl_matlab_print(vcl_cerr,R,"R",vnl_matlab_print_format_long);
      std::cerr << std::endl;
      vnl_matlab_print(vcl_cerr,S,"S",vnl_matlab_print_format_long);
      std::cerr << std::endl;
      ErrorExit(ERROR_OUT_OF_BOUNDS,
          "Internal Error: Sqrt of Rotation should not scale.\n");
    }

    double fnorm2 = (A - R).frobenius_norm();
    if (fnorm2 > eps)
    {
      ErrorExit(ERROR_OUT_OF_BOUNDS,
          "Internal Error: Sqrt should be a rotation.\n");
    }

  }
  else //affine 
  {
    vnl_matrix<double> R(3, 3), S(3, 3), A(3, 3);
    vnl_diag_matrix<double> D(3), I(3, 1.0);
    M.extract(A);
    MyMatrix::Polar2Decomposition(A, R, S, D);

    if (D[0] < 0.0 || D[1] < 0.0 || D[2] < 0.0)
      ErrorExit(ERROR_OUT_OF_BOUNDS,
          "Internal Error: Matrix Sqrt produced reflection.\n");

    // actually should not even get close to zero
    double eps = 0.001;
    if (D[0] < eps || D[1] < eps || D[2] < eps)
      ErrorExit(ERROR_OUT_OF_BOUNDS,
          "Internal Error: Matrix Sqrt produced near projection.\n");

  }
  return true;
}

/** Resample both images (or if not symmetric only src) to new space.
 Also return half way transforms.
 Checks members: symmetry and iscale
 Also rigid and Minit (if symmetry == true)
 */
void Registration::mapToNewSpace(const vnl_matrix_fixed<double, 4, 4>& M,
    double iscaleval, MRI * mriS, MRI* mriT, MRI *& mri_Swarp, MRI*& mri_Twarp,
    vnl_matrix_fixed<double, 4, 4>& mh, vnl_matrix_fixed<double, 4, 4>& mhi)
{
  if (verbose > 1)
    std::cout << " Registration::mapToNewSpace " << std::endl;
    
  if (iscaleonly)
  {
    mh.set_identity();
    mhi.set_identity();
    // adjust intensity   
    if (verbose > 1)
      std::cout << "   - adjusting only intensity ( " << iscaleval << " ) "
          << std::endl;
    // ISCALECHANGE:
    double si = sqrt(iscaleval);
    mri_Swarp = MyMRI::MRIvalscale(mriS, mri_Swarp, si);
    mri_Twarp = MyMRI::MRIvalscale(mriT, mri_Twarp, 1.0 / si);
    mri_Swarp->outside_val = mriS->outside_val;
    mri_Twarp->outside_val = mriT->outside_val;
    return;
  }
  else if (symmetry)
  {
    // here symmetrically warp both images SQRT(M)
    // this keeps the problem symmetric
    if (verbose > 1)
      std::cout << "   - resampling MOV and DST (sqrt)" << std::endl;
    // half way voxelxform
    //mh  = MyMatrix::MatrixSqrtAffine(M); // does not seem to work (creates imag results ...)?
//      vnl_matrix_fixed<double , 4, 4> &Mnew =&M;
//       // 2D correction
//       if ( mriS->depth == 1)
//       {
//         assert (mriT->depth ==1);
//         Mnew[2][0] = 0.0;
//         Mnew[0][2] = 0.0;
//         Mnew[2][1] = 0.0;
//         Mnew[1][2] = 0.0;
//         Mnew[2][2] = 1.0;
//         Mnew[2][3] = 0.0;
//       }
    mh = MyMatrix::MatrixSqrt(M); // symmetry slighlty destroyed here? !!

    // check if we moved out of our space:
    checkSqrtM(mh, rigid && Minit.empty()); // if minit was passed, it might be an affine initialization
        // also in old code for affine we checked fmd.first (here M) instead of mh (bug??)

    // do not just assume m = mh*mh, rather m = mh2 * mh
    // for transforming target we need mh2^-1 = mh * m^-1
    vnl_matrix_fixed<double, 4, 4> mi = vnl_inverse(M);
    mhi = mh * mi;
    //vnl_matlab_print(vcl_cerr,mh,"mh",vnl_matlab_print_format_long);std::cerr << std::endl;
    //vnl_matlab_print(vcl_cerr,mhi,"mhi",vnl_matlab_print_format_long);std::cerr << std::endl;

    // map both to mid space:
    if (mri_Swarp)
      MRIfree(&mri_Swarp);
    mri_Swarp = MRIclone(mriS, NULL);
    mri_Swarp = MyMRI::MRIlinearTransform(mriS, mri_Swarp, mh);
    if (mri_Twarp)
      MRIfree(&mri_Twarp);
    mri_Twarp = MRIclone(mriS, NULL); // set source geometry here: bring them to same space !! symmetry slightly destroyed here!!
    mri_Twarp = MyMRI::MRIlinearTransform(mriT, mri_Twarp, mhi);
    //cout << " Ss : " << mri_Swarp->width << " " << mri_Swarp->height << " " << mri_Swarp->depth << endl;
    //cout << " Ts : " << mri_Twarp->width << " " << mri_Twarp->height << " " << mri_Twarp->depth << endl;
    // probably better to:
    //   a) keep inputs at their input size in setSourceTarget
    //   b) compute overlap here (bbox)
    //   c) compute image dimensions based on bbox

    // adjust intensity symmetrically 
    if (iscale)
    {
      if (verbose > 1)
        std::cout << "   - adjusting intensity ( " << iscaleval << " ) "
            << std::endl;
      // ISCALECHANGE:
      double si = sqrt(iscaleval);
      MyMRI::MRIvalscale(mri_Swarp, mri_Swarp, si);
      MyMRI::MRIvalscale(mri_Twarp, mri_Twarp, 1.0 / si);
    }
    mri_Swarp->outside_val = mriS->outside_val;
    mri_Twarp->outside_val = mriT->outside_val;
    return;
  }
  else // resample at target location (using target geometry)
  {
    if (verbose > 1)
      std::cout << "   - resampling MOV to DST " << std::endl;
    if (mri_Swarp)
      MRIfree(&mri_Swarp);
    mri_Swarp = MRIclone(mriT, NULL); // set target geometry here
    mh = M;
//       // 2D correction
//       if ( mriS->depth == 1)
//       {
//         assert (mriT->depth ==1);
//         mh[2][0] = 0.0;
//         mh[0][2] = 0.0;
//         mh[2][1] = 0.0;
//         mh[1][2] = 0.0;
//         mh[2][2] = 1.0;
//         mh[2][3] = 0.0;
//       }
    mri_Swarp = MyMRI::MRIlinearTransform(mriS, mri_Swarp, mh);
    mhi.set_identity();
    if (!mri_Twarp)
      mri_Twarp = MRIcopy(mriT, mri_Twarp);
      
    // adjust intensity of source
    if (iscale)
    {
      if (verbose > 1)
        std::cout << "   - adjusting intensity ( " << iscaleval << " ) "
            << std::endl;
      // ISCALECHANGE:
      MyMRI::MRIvalscale(mri_Swarp, mri_Swarp, iscaleval);
    }
    mri_Swarp->outside_val = mriS->outside_val;
    mri_Twarp->outside_val = mriT->outside_val;
    return;
  }

}

std::pair<vnl_matrix_fixed<double, 4, 4>, double> Registration::convertP2Md(
    const vnl_vector<double>& p) const
{

  std::pair<vnl_matrix_fixed<double, 4, 4>, double> ret;
  ret.second = 0.0;

  int psize = p.size();

  if (iscale)
  {
    //std::cout << " has intensity " << std::endl;
    // last is intensity scale    
    // ISCALECHANGE:
    psize--;
    ret.second = (double) p[psize];
  }

  trans->setParameters(p);
  ret.first = trans->getMatrix();
  return ret;

}

// void Registration::testRotations()
// {
// 
//   int i;
//   std::vector < double > cost(24);
//   
//   for (i=0;i<24;i++)
//   {
//   
//   
//   
//   }
// 
// }
