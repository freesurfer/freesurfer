/**
 * @file RegPowell.cpp
 * @brief A class to compute a robust registration using Powell
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2011/10/07 22:28:51 $
 *    $Revision: 1.13 $
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
#include <iomanip>
#include "fs_vnl/fs_powell.h"
#include "fs_vnl/fs_cost_function.h"
#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_powell.h>

using namespace std;

RegPowell* RegPowell::tocurrent = NULL;
MRI * RegPowell::scf = NULL;
MRI * RegPowell::tcf = NULL;
int RegPowell::pcount = 0;
vnl_matrix_fixed < double, 4, 4> RegPowell::mh1 = vnl_matrix_fixed < double, 4, 4>();
vnl_matrix_fixed < double, 4, 4> RegPowell::mh2 = vnl_matrix_fixed < double, 4, 4>();
int RegPowell::icount = 0;
int RegPowell::subsamp = 1;
bool RegPowell::is2d = false;

class my_cost_function : public fs_cost_function
//class my_cost_function : public vnl_cost_function
{
private:
  double (*mFunction)(const vnl_vector<double> &);
public:
  my_cost_function( double (*function)(const vnl_vector<double> &), int pdim ):fs_cost_function(NULL),mFunction(function)
    {set_number_of_unknowns(pdim);};
//  my_cost_function( double (*function)(const vnl_vector<double> &), int pdim ):vnl_cost_function(pdim),mFunction(function){};
  virtual double f(const vnl_vector<double> & p){return ( *mFunction )(p);};
};



double RegPowell::costFunction(const vnl_vector < double >& p)
{
  //cout << "RegPowell::costFunction " << flush;

  // transform into matrix and iscale double
  static pair < vnl_matrix_fixed < double , 4, 4 > , double > Md;
  if (tocurrent->is2d)
  {
    Md = RegistrationStep<double>::convertP2Md2(p,tocurrent->iscale,tocurrent->rtype);  
  }
  else
    Md = RegistrationStep<double>::convertP2Md(p,tocurrent->rtype);
  Md.second = exp(Md.second); // compute full factor (source to target)
  //cout << endl;
  //cout << " M(p)   : " << endl << Md.first << endl;
  //cout << " psize : " << p.size() << endl;
  //vnl_matlab_print(vcl_cerr,Md.first,"M",vnl_matlab_print_format_long);
  //cout << " iscale : " << Md.second << endl;
  //vnl_matlab_print(vcl_cerr,p,"p",vnl_matlab_print_format_long);
  //vnl_matlab_print(vcl_cerr,Md.first,"M",vnl_matlab_print_format_long);
 
  // new full M = mh2 * cm * mh1
	Md.first = mh2 * Md.first * mh1;
  //cout << " M adj  : " << endl << Md.first << endl;
  //vnl_matlab_print(vcl_cerr,Md.first,"Madj",vnl_matlab_print_format_long);
  //exit(1);
  // maps from half way space (or target if not sym) back to source and to target
  static vnl_matrix_fixed<double,4,4> msi;
  static vnl_matrix_fixed<double,4,4> mti;
  
  if (tocurrent->symmetry)
  {
    // compute new half way maps and here assign to mti (hwspace to target)
    mti = MyMatrix::MatrixSqrt(Md.first);
    //vnl_matlab_print(vcl_cerr,mti,"Mti",vnl_matlab_print_format_long);
    // do not just assume m = mti*mti, rather m = mti * mh2
    // for sampling source we need mh2^-1 = m^-1 * mti
    msi = vnl_inverse(Md.first) * mti;
  }
  else
  {
    // only compute inverse to resample source at target
    msi = vnl_inverse(Md.first);
    mti.set_identity();
  }
  //vnl_matlab_print(vcl_cerr,mi,"minv",vnl_matlab_print_format_long);
  //cout << "tcf: " << MRIgetVoxVal(tcf,36,36,36,0) << " scf : " << MRIgetVoxVal(scf,36,36,36,0) << endl;
  //static vnl_matrix < double > HM;
  //static vnl_vector_fixed < double, 3 > ss;
  //ss[0] = tocurrent->subsamp; ss[1] = tocurrent->subsamp; ss[2] = tocurrent->subsamp;
  ////cout <<  ss[0]<< " " << flush;
  //hist2(HM, msi,mti,scf,tcf,ss);
  //hist2(HM, msi,mti,scf,tcf,tocurrent->subsamp,tocurrent->subsamp,tocurrent->subsamp);
  ////vnl_matlab_print(vcl_cerr,HM,"HM",vnl_matlab_print_format_long);
    
  static JointHisto H;
  H.create(scf,tcf,msi,mti,tocurrent->subsamp, tocurrent->subsamp, tocurrent->subsamp);
  //H.set(HM);
  //H.print();
  //exit(1);
  
  //H.clip(0.000001);
  
  H.smooth(7);
  double dd;
  switch (tocurrent->costfun)
  {
    case MI : dd=-H.computeMI();  break;
    case NMI: dd=-H.computeNMI(); break;
    case ECC: dd=-H.computeECC(); break;
    case NCC: dd=-H.computeNCC(); break;
    case LS:  dd= H.computeLS();  break;
    case SCR: dd=-H.computeSCR(); break;
    default:
      cout << " RegPowell::costFunction ERROR cannot deal with cost function " << tocurrent->costfun << " ! " << endl;
      exit(1);
  }
  
  // report eval:
  if (tocurrent->debug)
  {
    cout.setf(ios::fixed,ios::floatfield);   // floatfield set to fixed
    cout << " e = " << setprecision(14) << dd << "  @  " << flush;
    vnl_matlab_print(vcl_cerr,p,"p",vnl_matlab_print_format_long);
    if (icount ==0)
    {
      MRI* mri_Swarp = MRIclone(tcf,NULL);
      mri_Swarp = MyMRI::MRIlinearTransform(scf,mri_Swarp, vnl_inverse(msi));
      MRI* mri_Twarp = MRIclone(tcf,NULL); // bring them to same space (just use dst geometry)
      mri_Twarp = MyMRI::MRIlinearTransform(tcf,mri_Twarp, vnl_inverse(mti));
      MRIwrite(mri_Swarp,"shw.mgz");
      MRIwrite(mri_Twarp,"thw.mgz");

      
      MRIfree(&mri_Swarp);
      MRIfree(&mri_Twarp);
    }
  }
  icount++;
  
  //exit(1);
  return dd;

//   // compute new half way maps
//   vnl_matrix_fixed < double , 4, 4 > mh  = MyMatrix::MatrixSqrt(Md.first);
//   // do not just assume m = mh*mh, rather m = mh2 * mh
//   // for transforming target we need mh2^-1 = mh * m^-1
//   vnl_matrix_fixed < double , 4, 4 > mhi = mh1 * vnl_inverse(Md.first);
// 	
//   // map both to half way
//   MRI* mri_Swarp = MRIclone(scf,NULL);
//   mri_Swarp = MyMRI::MRIlinearTransform(scf,mri_Swarp, mh);
//   MRI* mri_Twarp = MRIclone(scf,NULL); // bring them to same space (just use src geometry)
//   mri_Twarp = MyMRI::MRIlinearTransform(tcf,mri_Twarp, mhi);
// 
//   // adjust intensity
//   if (tocurrent->iscale)
//   {
//     //cout << "   - adjusting intensity ( "<< fmd.second << " ) " << endl;
// 
//     MyMRI::MRIvalscale(mri_Swarp,mri_Swarp,(1.0+Md.second)*0.5);
//     MyMRI::MRIvalscale(mri_Twarp,mri_Twarp,(1.0+ 1.0/Md.second)*0.5);
//   }

  vnl_matrix_fixed < double , 4, 4 > mh;
  vnl_matrix_fixed < double , 4, 4 > mhi;
  MRI* mri_Swarp  = NULL;
  MRI* mri_Twarp  = NULL;
  tocurrent->mapToNewSpace(Md.first,Md.second,scf,tcf,mri_Swarp,mri_Twarp,mh,mhi);
  if (icount ==0)
  {
    MRIwrite(mri_Swarp,"shw.mgz");
    MRIwrite(mri_Twarp,"thw.mgz");
  }

  // compute error
  float d;
//  if (tocurrent->robust) d = (float) CostFunctions::tukeyBiweight(mri_Swarp,mri_Twarp,tocurrent->sat);
//  else d = (float) CostFunctions::leastSquares(mri_Swarp,mri_Twarp);

  d = (float) CostFunctions::normalizedMutualInformation(mri_Swarp,mri_Twarp);
  //d = (float) CostFunctions::MutualInformation(mri_Swarp,mri_Twarp);

   //cout << d << endl;

  icount++;
  if (icount%100 == 0) cout << "*" << flush;
//    if (icount%20 == 0)
    {
       //cout << endl << " iteration : " << icount << endl;
       vcl_cerr << icount << "  |  " << setprecision(8) << d << "    " << flush; vnl_matlab_print(vcl_cerr,p,"p",vnl_matlab_print_format_long);
//       int ii = icount / 20;
//       std::stringstream out;
//       out << ii;
//       MRIwrite(mri_Swarp,("shw-"+out.str()+".mgz").c_str());
//       MRIwrite(mri_Twarp,("thw-"+out.str()+".mgz").c_str());
    }

  // clean up
  MRIfree(&mri_Swarp);
  MRIfree(&mri_Twarp);

  return d;
}


void RegPowell::computeIterativeRegistration( int nmax,double epsit,MRI * mriS, MRI* mriT, const vnl_matrix < double > & m , double scaleinit)
// computes 4x4 matrix and iscale value in members iscalefinal and Mfinal
// nmax and epsit is ignored 
{
  cout << "RegPowell::computeIterativeRegistration" << endl;
  cout << " sym: " << symmetry << endl; 

  if (!mriS) mriS = mri_source;
  if (!mriT) mriT = mri_target;

  assert (mriS && mriT);

  is2d = false;
  if ( mriS->depth == 1 || mriT->depth == 1)
  {
    if (mriT->depth != mriS->depth){cout << "ERROR: both source and target need to be 2D or 3D" << endl; exit(1);}
    is2d = true;
  }
  
  bool cleanupS = true;
  if (mriS->type != MRI_UCHAR)
  {
    int no_scale_flag = FALSE;
    printf("changing data type from %d to %d (noscale = %d)...\n",
           mriS->type,MRI_UCHAR,no_scale_flag);
    mriS = MRISeqchangeType(mriS, MRI_UCHAR, 0.0, 0.999, no_scale_flag);
    if (mriS == NULL)
    {
      printf("ERROR: MRISeqchangeType\n");
      exit(1);
    }
    //cleanupS = true;
  }
  else
     mriS = MRIcopy(mriS,NULL);
  MyMRI::MRInorm255(mriS,mriS);
  
  bool cleanupT = true;
  if (mriT->type != MRI_UCHAR)
  {
    int no_scale_flag = FALSE;
    printf("changing data type from %d to %d (noscale = %d)...\n",
           mriT->type,MRI_UCHAR,no_scale_flag);
    mriT = MRISeqchangeType(mriT, MRI_UCHAR, 0.0, 0.999, no_scale_flag);
    if (mriT == NULL)
    {
      printf("ERROR: MRISeqchangeType\n");
      exit(1);
    }
    //cleanupT = true;
  }
  else
     mriT = MRIcopy(mriT,NULL);
  MyMRI::MRInorm255(mriT,mriT);
     
     

  tocurrent = this; // so that we can access this from static cost function
  rtype = 2;
  
  pair < vnl_matrix_fixed < double, 4, 4> , double > fmd(vnl_matrix_fixed < double, 4 , 4> () ,0.0);

  // check if mi (inital transform) is passed
  if (!m.empty())
  {
    fmd.first = m;
    if (verbose > 1) cout << "     initial M passed ... "  << endl;
  }
  else if (!Minit.empty())
  {
    fmd.first = getMinitResampled();
    if (verbose > 1) cout << "     initial M from resample ... "  << endl;
  }
  else
  {
    fmd.first = initializeTransform(mriS,mriT) ;
    if (verbose > 1) cout << "     initialize transform M ... "  << endl;
  }
  vnl_matrix < double > initialM = fmd.first;

  // ISCALECHANGE:
  // intensity model: R(s,IS,IT) = exp(-0.5 s) IT - exp(0.5 s) IS
  //                  R'  = -0.5 ( exp(-0.5 s) IT + exp(0.5 s) IS)
  //  thus iscalefinal= exp(0.5 s) * (1/exp(-0.5 s)) = exp(s)
  //     ==> s = log (iscalefinal)
  iscalefinal = 1.0;
  if (scaleinit != 1.0) iscalefinal = scaleinit;
	else iscalefinal = iscaleinit;
  fmd.second = log(iscalefinal);
  

  if (verbose > 1)
  {
    std::cout << "   - initial iscale: " << iscalefinal << std::endl;
    std::cout << "   - initial transform:\n" ;
		std::cout << fmd.first << std::endl;
   // std::cout << "   - initial det: " <<vnl_determinant(fmd.first) << endl ;
  }

  // decompose initial transform to apply in cost function (if symmetric):
  if (symmetry)
  {
    // here  symmetrically warp both images SQRT(M)
    // this keeps the problem symmetric
    mh1 = MyMatrix::MatrixSqrt(fmd.first);
    // do not just assume m = mh*mh, rather m = mh2 * mh
    // for transforming target we need mh2^-1 = mh * m^-1
    vnl_matrix_fixed < double, 4, 4 > mhi = mh1 * vnl_inverse(fmd.first);
	
    //set static
    mh2 = vnl_inverse(mhi); // M = mh2 * mh1
  }
  else
  {
     mh1 = fmd.first;
     mh2.set_identity();
  }
  
//cout << "M:" << endl << fmd.first << endl;
//cout << "m':" << endl << mh2*mh1 << endl;

  // set static pointers to images:
  scf = mriS;
  tcf = mriT;

  // create parameter vector:
  pcount = 3; // transolny
  if (rigid) pcount = 6;
  else pcount = 12;
  if (pcount==3) assert(transonly);
  if (is2d)
  { 
    pcount = 6;
    if (transonly) pcount = 2;
    else if (rigid) pcount = 3;
  }
  
  if (iscale) pcount++;

  subsamp = 1;
  if (subsamplesize <= 0) subsamplesize = 150;
  if (subsamplesize > 0 && mriT->width > subsamplesize && mriT->height > subsamplesize && mriT->depth > subsamplesize)
  {
    cout << "   - subsampling this resolution ... " << endl;
    subsamp = 2;
  }
  // compute Registration
  cout << "   - compute new registration ( " << pcount << " params )" << endl;

//pcount = 3;

//  double tols[6] = { 0.01, 0.01, 0.01, 0.001, 0.001 ,0.001};
  double tols[13] = { 0.02, 0.02, 0.02, 0.001, 0.001 ,0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001,0.001};
  // initial parameters
  vnl_vector < double > p(pcount,0.0);

  if (is2d)
  { 
    tols[2] = tols[1];
    if (pcount >= 6) p[3] = 1.0 ; p[4] = 1.0;
  }

  if (pcount >=12)
  {
     p[6] = 1.0; p[7] = 1.0; p[8] = 1.0;
  }
  if (iscale) p[pcount-1] = iscaleinit;
  // initial steps:
  vnl_matrix < double > xi(pcount,pcount,0.0);
  for (int i = 0;i<pcount;i++)
    xi[i][i] = 20*tols[i];
  
  if (iscale)
  {
    p[pcount-1] = tols[12];
    xi[pcount-1][pcount-1] = 20*tols[12];
  }
  
//p[0] = -0.2;
//p[1] = 0.13;
//p[2] = -0.15;
//p[3] = 0.08;
//p[4] = -0.03;
//p[5] = -0.01;
//p[0] =0.5187 ;  p[1] = -0.167 ; p[2] = -0.308;

//p[0] =0.1458980337503;

//    int nn = 200;
//    vnl_vector < double > sx(nn);
//    vnl_vector < double > sy(nn);
//    for (int ii = 0;ii<nn;ii++)
//    {
//      sx [ii] = (0.5*ii)/nn;
//      p[0] = sx[ii];
//      sy [ii] = costFunction(p);
//    }
//      vnl_matlab_print(vcl_cerr,sx,"sx",vnl_matlab_print_format_long); cout << endl;   
//      vnl_matlab_print(vcl_cerr,sy,"sy",vnl_matlab_print_format_long); cout << endl;   
//  exit(1);

  double min_sse = costFunction(p) ;
  cout << "   - initial cost: " << min_sse << endl;
//exit(1);


  icount = 0;
  cout << " START POWELL" << endl;
  my_cost_function myCostFunction( &costFunction , pcount );
  fs_powell minimizer( &myCostFunction );
  //vnl_powell minimizer( &myCostFunction );
 // double tol = 1e-4; //-8
  //int maxiter = 5;
  //minimizer.set_linmin_xtol(tol);
  //double tol = 1e-4; //
  double tol = 1e-5; //
  minimizer.set_x_tolerance(tol);
  minimizer.set_f_tolerance(tol);
  //minimizer.set_max_function_evals(maxiter);

  //minimizer.set_trace(0);
  //minimizer.set_verbose(1);

  int returnCode = minimizer.minimize( p, &xi );
  //int returnCode = minimizer.minimize( p );

  // exit if failure
  if ( returnCode == vnl_nonlinear_minimizer::FAILED_TOO_MANY_ITERATIONS )
  {
    ErrorExit(ERROR_BADPARM, "powell exceeding maximum iterations.");
  }
  else if ( returnCode == fs_powell::ERROR_FAILURE ||
            returnCode == fs_powell::ERROR_DODGY_INPUT ||
            returnCode == fs_powell::FAILED_FTOL_TOO_SMALL ||
            returnCode == fs_powell::FAILED_XTOL_TOO_SMALL ||
            returnCode == fs_powell::FAILED_GTOL_TOO_SMALL )
  {
    ErrorExit(ERROR_BADPARM, "powell error.");
  }
  else
  {
    // success
    cout << " iter: " << minimizer.get_num_iterations() << endl;
    cout << " final e = " << setprecision(14) << minimizer.get_end_error() << endl;
    vnl_matlab_print(vcl_cerr,p,"p",vnl_matlab_print_format_long); cout << endl;   
  }
  
  //OpenPowell(p, xi, pcount, tol, &iter, &fret, costFunction);
 // OpenPowell2(p, xi, pcount, tol,tol,maxiter, &iter, &fret, costFunction);
  //cout << endl << "best alignment initial powell: " << fret << " (" << iter << " steps)" << endl;

//   int count = 0;
//   float flast;
//   do
//   {
//     count++;
//     // reinitialize powell directions
//     for (int r = 1 ; r <= pcount ; r++)
//     {
//       for (int c = 1 ; c <= pcount ; c++)
//       {
//         xi[r][c] = r == c ? 20*tols[r-1] : 0 ;
//       }
//     }
// 
//     flast = fret ;
//     icount = 0;
//     OpenPowell(p, xi, pcount, tol, &iter, &fret, costFunction);
//     cout << endl << " return from " << count << " powell: " << fret << " (" << iter << " steps)" << endl;
// 
// 
// //    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
//     //   printf("best alignment after powell: %2.3f (%d steps)\n", fret, iter) ;
//     //cout << "best alignment after powell: " << fret << " (" << iter << " steps)" << endl;
// 
//     if ((flast-fret)/flast < tol)
//       break ;
//   }
//   while (fret < flast) ;

//  vnl_vector < double > v(pcount);
//  for (int i = 1;i<=pcount; i++)
//    v[i-1] = p[i];
// // MatrixPrintFmt(stdout,"% 2.8f",v);
//	cout << v << endl;

  //if (fmd.first) MatrixFree(&fmd.first);
  if (is2d)
   fmd = RegistrationStep<double>::convertP2Md2(p,iscale,rtype);  
  else
   fmd = RegistrationStep<double>::convertP2Md(p,rtype);
  //vnl_matlab_print(vcl_cerr,fmd.first,"M",vnl_matlab_print_format_long); cout << endl;   

  if (symmetry)
  {
    // new M = mh2 * cm * mh1
    fmd.first = (mh2 * fmd.first) * mh1;
  }
  else fmd.first = fmd.first * initialM;
  //vnl_matlab_print(vcl_cerr,fmd.first,"M",vnl_matlab_print_format_long); cout << endl;   

  // ISCALECHANGE:
  if (iscale)
  {
    iscalefinal = exp(fmd.second); // compute full factor (source to target)
    //idiff = fabs(cmd.second);
    //std::ostringstream istar;
    //if (idiff <= ieps) istar << " <= " << ieps << "  :-)" ;
    //if (verbose >0 ) std::cout << "     -- intensity log diff: abs(" << cmd.second << ") " << istar.str() << std::endl;
  }

  // adjust half way maps to new midpoint based on final transform
  if (verbose >1) std::cout << "     -- adjusting half-way maps " << std::endl;
  vnl_matrix_fixed < double, 4, 4 > ch = MyMatrix::MatrixSqrt(fmd.first);
  // do not just assume c = ch*ch, rather c = ch2 * ch
  // for transforming target we need ch2^-1 = ch * c^-1
  vnl_matrix_fixed < double, 4, 4 >  ci  = vnl_inverse(fmd.first);
  vnl_matrix_fixed < double, 4, 4 >  chi = ch*ci;
  mov2weights = ch; 
  dst2weights = chi;

  Mfinal = fmd.first;

  cout << endl << " DONE " << endl;
  //cout << endl << "Final Transform:" << endl;
  
  //MatrixPrintFmt(stdout,"% 2.8f",Mfinal);
	//cout << Mfinal << endl;

  //MRI* mri_Swarp = MRIclone(mriT,NULL);
  //mri_Swarp = MyMRI::MRIlinearTransform(mriS,mri_Swarp,Mfinal);
  //MRIwrite(mriS,"mriS.mgz");
  //MRIwrite(mri_Swarp,"mriSwarp.mgz");
  //MRIwrite(mriT,"mriT.mgz");
  //MRIfree(&mri_Swarp);
  
  if (cleanupS) MRIfree(&mriS);
  if (cleanupT) MRIfree(&mriT);

}
