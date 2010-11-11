/**
 * @file  mri_robust_register.cpp
 * @brief Linear registration of two volumes using robust statistics
 *
 * See also "Robust Multiresolution Alignment of MRI Brain Volumes"
 * by Nestares and Heeger (2000)
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/11/11 22:32:17 $
 *    $Revision: 1.43 $
 *
 * Copyright (C) 2008-2012
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


//
// mri_robust_register.cpp
//
// written by Martin Reuter
// Nov. 4th ,2008
//
////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include <vnl/vnl_inverse.h>
#include <vnl/vnl_matrix_fixed.h>

#include "Registration.h"
#include "Regression.h"
#include "RegPowell.h"
#include "CostFunctions.h"
#include "MyMRI.h"
#include "MyMatrix.h"
#include "mri_robust_register.help.h"

// all other software are all in "C"
#ifdef __cplusplus
extern "C"
{
#endif
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"

#ifdef __cplusplus
}
#endif

using namespace std;

#define SAT -1 // leave blank, either passed by user or --satit
//#define SAT 4.685 // this is suggested for gaussian noise
//#define SAT 20
#define SSAMPLE -1

struct Parameters
{
  string mov;
  string dst;
  string lta;
  string maskmov;
  string maskdst;
  string halfmov;
  string halfdst;
  string halfweights;
  string halfmovlta;
  string halfdstlta;
  string weightsout;
  bool   satit;
  bool   satest;
  bool   nomulti;
  bool   conform;
  bool   floattype;
  bool   lta_vox2vox;
  bool   affine;
  bool   iscale;
  bool   transonly;
  string transform;
  bool   leastsquares;
  int    iterate;
  double epsit;
  double sat;
  bool   warp;
  string warpout;
  int    subsamplesize;
  int    debug;
  MRI*   mri_mov;
  MRI*   mri_dst;
  bool   dosatest;
  bool   initorient;
	bool   inittrans;
	int    verbose;
	int    highit;
	bool   doubleprec;
	double wlimit;
	bool   oneminusweights;
	bool   symmetry;
	string iscaleout;
};
static struct Parameters P =
  { "","","","","","","","","","","",false,false,false,false,false,false,false,false,false,"",false,5,0.01,SAT,false,"",SSAMPLE,0,NULL,NULL,false,false,true,1,-1,false,0.16,false,true,""
  };


static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[],Parameters & P) ;
static void initRegistration(Registration & R, Parameters & P) ;

static char vcid[] = "$Id: mri_robust_register.cpp,v 1.43 2010/11/11 22:32:17 mreuter Exp $";
char *Progname = NULL;

//static MORPH_PARMS  parms ;
//static FILE *diag_fp = NULL ;


void conv(MRI * i)
{
  cout << " adsf" << endl;
  Registration R;
  MRI * fmri = MRIalloc(i->width,i->height,i->depth,MRI_FLOAT);
  MRIcopyHeader(i,fmri);
  fmri->type = MRI_FLOAT;
  float f;
  for (int z=0;z<i->depth;z++)
    for (int y=0;y<i->height;y++)
      for (int x=0;x<i->width;x++)
      {
        f = MRIgetVoxVal(i,x,y,z,0);
        // cout << " f " << f << endl;
        MRIsetVoxVal(fmri,x,y,z,0,f);
      }
  cout << "asdfasdf" << endl;
  MRIwrite(fmri,"float-1.mgz");
  MRI * sfmri;
  sfmri = MyMRI::MRIvalscale(fmri,NULL,100);
  MRIwrite(sfmri,"float-100.mgz");
  sfmri = MyMRI::MRIvalscale(fmri,sfmri,1000);
  MRIwrite(sfmri,"float-1000.mgz");
  exit(0);
}

void testRegression()
{

  int n = 200;
  vnl_matrix < double > A(n,1);
  vnl_vector < double > b(n);

  for (int i = 0;i<n; i++)
	{
       A[i][0] = i;
       b[i]    = 4*i;
	}
  for (int i = 0;i<n; i+=5)
	{
    b[i] = 0;
	}
	
	Regression<double> R1(A,b);
  vnl_vector < double >  M1 = R1.getLSEst();
	cout << M1 << endl;
	cout << endl <<endl;

	Regression<double> R2(A,b);
  vnl_vector < double >  M2 = R2.getRobustEst();
	cout << M1 << endl;
	cout << endl <<endl;

  exit(0);

}


int main(int argc, char *argv[])
{
  { // for valgrind, so that everything is freed
  cout << vcid << endl << endl;
//  setenv("SURFER_FRONTDOOR","",1) ;
  // set the environment variable
  // to store mri as chunk in memory:
//  setenv("FS_USE_MRI_CHUNK","",1) ;
  if (getenv("FS_USE_MRI_CHUNK") != NULL)
  {
    cerr << "Error: do not set FS_USE_MRI_CHUNK while it is still buggy!" << endl;
    exit(1);
  }

  // Default initialization
  int nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
//  DiagInit(NULL, NULL, NULL) ;

  if (!parseCommandLine(argc, argv, P))
  {
    printUsage();
    exit(1);
  }

  // Timer
  struct timeb start ;
  int    msec,minutes,seconds;
  TimerStart(&start) ;


  // init registration from Parameters
//  RegPowell R;
  Registration R;
  initRegistration(R,P);
//conv(P.mri_dst);

//  cout << " mean mov : " << CostFunctions::mean(P.mri_mov) << "  mean dst: " << CostFunctions::mean(P.mri_dst) << endl;
//  cout << " sdev mov : " << CostFunctions::sdev(P.mri_mov) << "  sdev dst: " << CostFunctions::sdev(P.mri_dst) << endl;
//  cout << " median mov : " << CostFunctions::median(P.mri_mov) << "  median dst: " << CostFunctions::median(P.mri_dst) << endl;
//  cout << " mad mov : " << CostFunctions::mad(P.mri_mov) << "  mad dst: " << CostFunctions::mad(P.mri_dst) << endl;
// does not work with different image dimensions:
//  cout << " LS difference before: " << CostFunctions::leastSquares(P.mri_mov,P.mri_dst) << endl;
//  cout << " NC difference before: " << CostFunctions::normalizedCorrelation(P.mri_mov,P.mri_dst) << endl;

  // compute Alignment
  //std::pair <MATRIX*, double> Md;
  if (P.satest) R.computeSatEstimate(2,P.iterate,P.epsit);
//  else if (P.satit) Md = R.computeIterativeRegSat(P.iterate,P.epsit);
  else if (P.satit) {  R.findSaturation(); R.computeMultiresRegistration(0,P.iterate,P.epsit); }
  else if (P.nomulti) R.computeIterativeRegistration(P.iterate,P.epsit);
  else R.computeMultiresRegistration(0,P.iterate,P.epsit);

  if (P.satest) // old stuff, can be removed ?
  {
    cout << "run:" << endl;
    cout << " gnuplot " << R.getName() << "-sat.plot ; \\ " << endl;
    cout << " epstopdf " << R.getName() << "-sat.eps " << endl;
    cout << " and view the pdf " << endl << endl;
    msec = TimerStop(&start) ;
    seconds = nint((float)msec/1000.0f) ;
    minutes = seconds / 60 ;
    seconds = seconds % 60 ;
    cout << "registration took "<<minutes<<" minutes and "<<seconds<<" seconds." << endl;
    exit(0);
  }

  //Md.first = MatrixReadTxt("xform.txt",NULL);
  //Md.second = 1;

  // Print results:
  std::pair <MATRIX*, double> Md;
  cout << endl << "Final Transform:" << endl;
  Md.first = MyMatrix::convertVNL2MATRIX(R.getFinalVox2Vox());
	Md.second = R.getFinalIscale();
  MatrixPrintFmt(stdout,"% 2.8f",Md.first);
  if (R.isIscale()) cout << "Intenstiy Scale Factor: " << Md.second << endl;
  cout << endl ;


  // writing transform section here
  cout << "writing output transformation to "<<P.lta <<" ..." << endl;
  char reg[STRLEN];
  strcpy(reg, P.lta.c_str());
  LTA * lta = LTAalloc(1,P.mri_mov);
  if (!P.lta_vox2vox) // do ras to ras
  {
    cout << "converting VOX to RAS and saving RAS2RAS..." << endl ;
    //cout << "VOX2VOX:" << endl ;
    //MatrixPrint(stdout, Md.first) ;
    lta->xforms[0].m_L = MRIvoxelXformToRasXform (P.mri_mov, P.mri_dst, Md.first, lta->xforms[0].m_L) ;
    //cout << "RAS2RAS:" << endl ;
    //MatrixPrint(stdout,lta->xforms[0].m_L) ;
    lta->type = LINEAR_RAS_TO_RAS ;
  }
  else // vox to vox
  {
    cout << "saving VOX2VOX..." << endl ;
    lta->xforms[0].m_L = MatrixCopy(Md.first, lta->xforms[0].m_L) ;
    lta->type = LINEAR_VOX_TO_VOX ;
  }
  // add src and dst info
  getVolGeom(P.mri_mov, &lta->xforms[0].src);
  getVolGeom(P.mri_dst, &lta->xforms[0].dst);
  LTAwriteEx(lta, reg) ;

  if (R.isIscale() && Md.second >0 && P.iscaleout != "")
//  if (R.isIscale() && Md.second >0)
  {
    //string fn;
    //if (P.iscaleout != "") fn = P.iscaleout;
    //else fn = R.getName() + "-intensity.txt";
    //ofstream f(fn.c_str(),ios::out);
    ofstream f(P.iscaleout.c_str(),ios::out);
    f << Md.second;
    f.close();
  }

  //  MatrixWriteTxt("xform.txt",Md.first);
  // end of writing transform

  // here do scaling of intensity values
  if (R.isIscale() && Md.second > 0)
  {
    cout << "Adjusting Intensity of MOV by " << Md.second << endl;
    P.mri_mov = MyMRI::MRIvalscale(P.mri_mov, P.mri_mov, Md.second);
  }

  // maybe warp source to target:
  if (P.warpout != "")
  {
    //cout << "using lta" << endl;
    int nframes = P.mri_mov->nframes;
		if (P.mri_mov->nframes > 1) cout << " WARNING: movable has more than one frame !!! Only warp first ..." << endl;
    P.mri_mov->nframes = 1 ; // only map frame 1
    MRI *mri_aligned = MRIclone(P.mri_dst,NULL);
    mri_aligned = LTAtransform(P.mri_mov,mri_aligned, lta);
    P.mri_mov->nframes = nframes ;
		
		// keep acquisition params:
		MRIcopyPulseParameters(P.mri_mov,mri_aligned);
		
//    sprintf(fname, "%s_after_final_alignment", parms.base_name) ;
//    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
//    sprintf(fname, "%s_target", parms.base_name) ;
//    MRIwriteImageViews(mri_dst, fname, IMAGE_SIZE) ;

//      cout << " mean warp : " << CostFunctions::mean(mri_aligned) << "  mean dst: " << CostFunctions::mean(P.mri_dst) << endl;
//      cout << " sdev warp : " << CostFunctions::sdev(mri_aligned) << "  sdev dst: " << CostFunctions::sdev(P.mri_dst) << endl;
//      cout << " median warp : " << CostFunctions::median(mri_aligned) << "  median dst: " << CostFunctions::median(P.mri_dst) << endl;
//      cout << " mad warp : " << CostFunctions::mad(mri_aligned) << "  mad dst: " << CostFunctions::mad(P.mri_dst) << endl;
// does not work with different image dimensions:
//    cout << " LS difference after: " << CostFunctions::leastSquares(mri_aligned,P.mri_dst) << endl;
//    cout << " NC difference after: " << CostFunctions::normalizedCorrelation(mri_aligned,P.mri_dst) << endl;

    MRIwrite(mri_aligned, P.warpout.c_str()) ;
    MRIfree(&mri_aligned) ;

    cout << endl;
    cout << "To check aligned result, run:" << endl;
    cout << "  tkmedit -f "<< P.dst <<" -aux " << P.warpout << endl;
  }

  // maybe write out weights in target space:
  if (P.weightsout!="")
  {

    MRI * mri_weights = R.getWeights(); // in target space
    if (mri_weights != NULL)
    {
		
		  if (P.oneminusweights) mri_weights = MRIlinearScale(mri_weights,NULL,-1,1,0);
		  MRIwrite(mri_weights,P.weightsout.c_str()) ;
			if (P.oneminusweights) mri_weights = R.getWeights();
		
//       // map to target and use target geometry 
//       std::pair < vnl_matrix_fixed < double, 4, 4>, vnl_matrix_fixed < double, 4, 4> > map2weights = R.getHalfWayMaps();
//       vnl_matrix_fixed < double, 4, 4> hinv = vnl_inverse(map2weights.second);
//       MRI * wtarg = MRIalloc(P.mri_dst->width,P.mri_dst->height,P.mri_dst->depth,MRI_FLOAT);
//       MRIcopyHeader(P.mri_dst,wtarg);
//       wtarg->type = MRI_FLOAT;
//       wtarg = MyMRI::MRIlinearTransform(mri_weights,wtarg, hinv);
//       MRIwrite(wtarg, P.weightsout.c_str()) ;
//       MRIfree(&wtarg);
//       //MatrixFree(&hinv);
      cout << "or even overlay the weights:" <<endl;
      cout << "  tkmedit -f "<< P.dst <<" -aux "<< P.warpout << " -overlay " << P.weightsout <<endl;
    }
    else
      cout << "Warning: no weights have been computed! Maybe you ran with --leastsquares??" << endl;
  }

  // write out images in half way space
  if (P.halfmov != "" || P.halfdst != "" || P.halfweights != "" || P.halfdstlta != "" || P.halfmovlta != "")
  {
    cout << endl;
    cout << "Creating half way data ..." << endl;
    std::pair < vnl_matrix_fixed < double, 4, 4>, vnl_matrix_fixed < double, 4, 4> > maps2weights = R.getHalfWayMaps();
    MRI * mri_weights = R.getWeights();
				
    LTA * m2hwlta = LTAalloc(1,P.mri_mov);
    LTA * d2hwlta = LTAalloc(1,P.mri_dst);
    if (!P.lta_vox2vox) // do ras to ras
    {
      // cout << "converting VOX to RAS and saving RAS2RAS..." << endl ;
      // (use geometry of destination space for half-way)
      m2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform (P.mri_mov, mri_weights, MyMatrix::convertVNL2MATRIX(maps2weights.first), m2hwlta->xforms[0].m_L) ;
      m2hwlta->type = LINEAR_RAS_TO_RAS ;
      d2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform (P.mri_dst, mri_weights, MyMatrix::convertVNL2MATRIX(maps2weights.second), d2hwlta->xforms[0].m_L) ;
      d2hwlta->type = LINEAR_RAS_TO_RAS ;
    }
    else // vox to vox
    {
      // cout << "saving VOX2VOX..." << endl ;
      //m2hwlta->xforms[0].m_L = MatrixCopy(maps2weights.first, m2hwlta->xforms[0].m_L) ;
      m2hwlta->xforms[0].m_L = MyMatrix::convertVNL2MATRIX(maps2weights.first, m2hwlta->xforms[0].m_L) ;
      m2hwlta->type = LINEAR_VOX_TO_VOX ;
      //d2hwlta->xforms[0].m_L = MatrixCopy(maps2weights.second, d2hwlta->xforms[0].m_L) ;
      d2hwlta->xforms[0].m_L = MyMatrix::convertVNL2MATRIX(maps2weights.second, d2hwlta->xforms[0].m_L) ;
      d2hwlta->type = LINEAR_VOX_TO_VOX ;
    }
    // add src and dst info (use mri_weights as target geometry in both cases)
    getVolGeom(P.mri_mov, &m2hwlta->xforms[0].src);
    getVolGeom(mri_weights, &m2hwlta->xforms[0].dst);
    getVolGeom(P.mri_dst, &d2hwlta->xforms[0].src);
    getVolGeom(mri_weights, &d2hwlta->xforms[0].dst);

    // write lta to half way
    if (P.halfmovlta != "") LTAwriteEx(m2hwlta, P.halfmovlta.c_str()) ;
    if (P.halfdstlta != "") LTAwriteEx(d2hwlta, P.halfdstlta.c_str()) ;

    if (P.halfmov != "")
    {
      cout << " creating half-way movable ..." << endl;
      // take dst geometry info from lta:
      MRI* mri_Swarp = LTAtransform(P.mri_mov,NULL, m2hwlta);

      //cout << " MOV       RAS: " << P.mri_mov->c_r << " , " <<	P.mri_mov->c_a << " , " <<	P.mri_mov->c_s << endl;
      //cout << " DST       RAS: " << P.mri_dst->c_r << " , " <<	P.mri_dst->c_a << " , " <<	P.mri_dst->c_s << endl;
      //cout << " weights   RAS: " << mri_weights->c_r << " , " <<	mri_weights->c_a << " , " <<	mri_weights->c_s << endl;
      //cout << " Swarp_old RAS: " << mri_Swarp_old->c_r << " , " <<	mri_Swarp_old->c_a << " , " <<	mri_Swarp_old->c_s << endl;
      //MRI* mri_Swarp = MRIalloc(mri_weights->width, mri_weights->height, mri_weights->depth, P.mri_mov->type);
      //MRIcopyHeader(mri_weights,mri_Swarp);
      //mri_Swarp->type = P.mri_mov->type;
      //LTAtransform(P.mri_mov,mri_Swarp, m2hwlta);
      //cout << " Swarp     RAS: " << mri_Swarp->c_r << " , " <<	mri_Swarp->c_a << " , " <<	mri_Swarp->c_s << endl;
		  MRIcopyPulseParameters(P.mri_mov,mri_Swarp);
      MRIwrite(mri_Swarp,P.halfmov.c_str());

      if (P.debug)
      {
        MRIiterator mw(mri_weights);
        MRIiterator ms(mri_Swarp);
        double meanw1=0, meanw0=0, mean = 0, meanw = 0, countw = 0;
        int countw1=0,countw0=0,count=0;
        for (ms.begin(); !ms.isEnd(); ms++)
        {
          if (fabs(*mw )>0.0001)
          {
            meanw0+= (*ms);
            countw0++;
          }
          if (fabs(*mw-1.0) < 0.0001)
          {
            meanw1+= *ms;
            countw1++;
          }

          mean+= *ms;
          count++;

          meanw+= *ms * *mw;
          countw+= *mw;

          assert(! (mw.isEnd() && !ms.isEnd()));
          mw++;
        }
			
        cout << " mov int means: " << mean/count << " ( " << count << " )  w0: " << meanw0/countw0 << " ( " << countw0 << " ) w1: " << meanw1/countw1 << " ( " << countw1 << " )  weighted: " << meanw/countw<<" ( " << countw << " )" << endl;
      }

      MRIfree(&mri_Swarp);

      //MRIwrite(P.mri_mov,"movable-original.mgz");
      //mri_Swarp = R.makeConform(P.mri_mov,NULL,false,true);
      //MRIwrite(mri_Swarp,"movable-uhar.mgz");
      //MRI * tttemp = MRIclone(mri_Swarp,NULL);
      //tttemp =  MRIlinearTransform(mri_Swarp,tttemp, mh);
      //MRIwrite(tttemp,"movable-uhar-half.mgz");
      //MRIfree(&mri_Swarp);
      //MRIfree(&tttemp);

    }
    if (P.halfdst != "")
    {
      cout << " creating half-way destination ..." << endl;
      MRI* mri_Twarp = LTAtransform(P.mri_dst,NULL, d2hwlta);
		  MRIcopyPulseParameters(P.mri_dst,mri_Twarp);
      MRIwrite(mri_Twarp,P.halfdst.c_str());
      MRI * mri_weights = R.getWeights();
			
      if (P.debug)
      {
        MRIiterator mw(mri_weights);
        MRIiterator ms(mri_Twarp);
        double meanw1=0, meanw0=0, mean = 0, meanw = 0, countw = 0;
        int countw1=0,countw0=0,count=0;
        for (ms.begin(); !ms.isEnd(); ms++)
        {
          if (fabs(*mw )>0.0001)
          {
            meanw0+= (*ms);
            countw0++;
          }
          if (fabs(*mw-1.0) < 0.0001)
          {
            meanw1+= *ms;
            countw1++;
          }

          mean+= *ms;
          count++;

          meanw+= *ms * *mw;
          countw+= *mw;

          assert(! (mw.isEnd() && !ms.isEnd()));
          mw++;
        }
        cout << " mov int means: " << mean/count << " ( " << count << " )  w0: " << meanw0/countw0 << " ( " << countw0 << " ) w1: " << meanw1/countw1 << " ( " << countw1 << " )  weighted: " << meanw/countw<<" ( " << countw << " )" << endl;
      }
			
      MRIfree(&mri_Twarp);
    }
    if (P.halfweights != "")
    {
      //MRI * mri_weights = R.getWeights();
      if (mri_weights != NULL)
      {
        cout << " saving half-way weights ..." << endl;
        MRI* mri_wtemp = LTAtransform(mri_weights,NULL, d2hwlta);
		    if (P.oneminusweights) mri_wtemp = MRIlinearScale(mri_wtemp,mri_wtemp,-1,1,0);
        MRIwrite(mri_wtemp,P.halfweights.c_str());
				MRIfree(&mri_wtemp);
        //MRIwrite(mri_weights,P.halfweights.c_str());
      }
      else
        cout << "Warning: no weights have been computed! Maybe you ran with --leastsquares??" << endl;
    }
  }

  if (P.debug >0)
  {
    cout << endl;
    cout << "To check debug output, run:" << endl;
    std::string name = R.getName();
    cout << "  tkmedit -f " << name << "-mriS-warp.mgz -aux " << name << "-mriT-warp.mgz -overlay " << name << "-mriS-weights.mgz" << endl;
  }

  cout << endl;
  cout << "To check transform, run:" << endl;
  cout << "  tkregister2 --mov "<< P.mov <<" --targ " << P.dst <<" --lta " << P.lta << " --reg " << R.getName() << ".reg" << endl;



  // cleanup
  if (Md.first) MatrixFree(&Md.first) ;
  if (P.mri_mov) MRIfree(&P.mri_mov);
  if (P.mri_dst) MRIfree(&P.mri_dst);

  ///////////////////////////////////////////////////////////////
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  cout << endl << "Registration took "<<minutes<<" minutes and "<<seconds<<" seconds." << endl;
  //if (diag_fp) fclose(diag_fp) ;
  } // for valgrind, so that everything is free
  exit(0) ;
  return(0) ;
}

// int main(int argc, char *argv[])
// {
//
//   char         *fname_src, *fname_dst, *fname_out, fname[STRLEN];
//   MRI          *mri_src, *mri_dst, *mri_tmp;
//
//   int          nargs,ninputs,i,msec,minutes,seconds;
//   struct timeb start ;
//
//   // defaults
//   Progname = argv[0] ;
//
//   Registration R; // sets its own default parameters
//
// //  int ac = argc ;
// //  char **av = argv ;
//   for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
//   {
//     nargs = get_option(argc, argv,R) ;
//     argc -= nargs ;
//     argv += nargs ;
//   }
//
//   if (argc < 4)
//   {
//     printUsage();
//     exit(1);
//   }
//
//   ninputs = argc-3 ;
//   //cout << "reading "<<ninputs<<" input volumes..."<< endl;
//   fname_dst = argv[ninputs+1] ;
//   fname_out = argv[ninputs+2] ;
//   FileNameOnly(fname_out, fname) ;
//   FileNameRemoveExtension(fname, fname) ;
//   strcpy(parms.base_name, fname) ;
//  // Gdiag |= DIAG_WRITE ;
//  // cout << "logging results to "<< parms.base_name <<".log" << endl;
//
//   TimerStart(&start) ;
//   ///////////  read MRI Target //////////////////////////////////////////////////
//   cout << endl << "reading target '"<<fname_dst<<"'..."<< endl;;
//   fflush(stdout) ;
//   mri_dst = MRIread(fname_dst) ;
//   if (mri_dst == NULL)
//     ErrorExit(ERROR_NOFILE, "%s: could not open MRI Target %s.\n",
//               Progname, fname_dst) ;
//
//    //////////////////////////////////////////////////////////////
//   // create a list of MRI volumes
//   cout << "reading "<<ninputs<<" source (movable) volumes..."<< endl;
//   for (i = 0 ; i < ninputs ; i++) {
//     fname_src = argv[i+1] ;
//     cout << "reading source '"<<fname_src<<"'..." << endl;
//     fflush(stdout) ;
//     mri_tmp = MRIread(fname_src) ;
//     if (!mri_tmp)
//       ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
//                 Progname, fname_src) ;
//
// //    TRs[i] = mri_tmp->tr ;
// //    fas[i] = mri_tmp->flip_angle ;
// //    TEs[i] = mri_tmp->te ;
//
// //    if (mask_fname) {
// //      MRI *mri_mask ;
// //
// //      mri_mask = MRIread(mask_fname) ;
// //      if (!mri_mask)
// //        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
// //                  Progname, mask_fname) ;
// //      MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
// //      MRIfree(&mri_mask) ;
// //    }
//     if (i == 0)
//     {
//       mri_src = MRIallocSequence(mri_tmp->width,
//                                 mri_tmp->height,
//                                 mri_tmp->depth,
//                                 mri_tmp->type,
//                                 ninputs) ;
//       MRIcopyHeader(mri_tmp, mri_src) ;
//     }
//     MRIcopyFrame(mri_tmp, mri_src, 0, i) ;
//     MRIfree(&mri_tmp) ;
//   }
//   cout << endl;
//   //////////////////////////////////////////////////////////////
//
// //  if (!FZERO(blur_sigma)) {
// //    MRI *mri_tmp, *mri_kernel ;
//
// //    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
// //    mri_tmp = MRIconvolveGaussian(mri_in, NULL, mri_kernel) ;
// //    MRIfree(&mri_in) ;
// //    mri_in = mri_tmp ;
// // }
//
//     ////////////////////////////////////////////////////
//     // now start working (remember this is vox-to-vox transform)
//   //  parms.lta->xforms[0].m_L = MatrixIdentity(4, NULL) ;
//
//   // real work done here
// //  MATRIX * Minit = MRIgetVoxelToVoxelXform(mri_src,mri_dst) ;
// //  cout << "initial transform:\n" ;
// //  MatrixPrintFmt(stdout,"% 2.8f",Minit);
// //  //std::pair <MATRIX*, double> Md = R.computeIterativeRegistration(10,mri_src,mri_dst,Minit);
// //  std::pair <MATRIX*, double> Md = R.computeMultiresRegistration(mri_src,mri_dst,Minit);
// //  std::pair <MATRIX*, double> Md = R.computeMultiresRegistration(mri_src,mri_dst);
//
//   std::pair <MATRIX*, double> Md;
//   if (nit > 0)  Md = R.computeIterativeRegistration(nit,mri_src,mri_dst);
//   else Md = R.computeMultiresRegistration(mri_src,mri_dst);
//
//
//   cout << "final transform:\n" ;
//   MatrixPrintFmt(stdout,"% 2.8f",Md.first);
//   if (R.isIscale()) cout << " iscale: " << Md.second << endl;
//   cout << endl ;
//
//   /////////////////////////diagnostics/////////////////////////////////
// //  if ((Gdiag & DIAG_WRITE) && (parms.write_iterations != 0)) {
//     MRI *mri_aligned ;
//     int nframes = mri_src->nframes;
//     mri_src->nframes = 1 ;
//     mri_aligned = MRIlinearTransform(mri_src, NULL,Md.first);
//     // here also do scaling of intensity values
//     mri_src->nframes = nframes ;
// //    sprintf(fname, "%s_after_final_alignment", parms.base_name) ;
// //    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
// //    sprintf(fname, "%s_target", parms.base_name) ;
// //    MRIwriteImageViews(mri_dst, fname, IMAGE_SIZE) ;
//     sprintf(fname, "%s.mgz", parms.base_name) ;
//     MRIwrite(mri_aligned, fname) ;
//     MRIfree(&mri_aligned) ;
//     cout << "To check results, run:" << endl;
//     cout << "  tkmedit -f "<< fname_dst <<" -aux " << fname << endl;
// //    cout << " or " << endl;
// //    cout << "  tkmedit -f "<< fname_src <<" -aux " << fname << endl;
//
// //  }
//
//    /////////////////////////////////////////////////////////////////////
//   cout << "writing output transformation to "<<fname_out<<"..." << endl;
//   // writing transform section here
//   // create gca volume for outputting dirction cosines and c_(ras)
//   //mri_dst = MRIallocHeader(gca->width, gca->height, gca->depth, mri_in->type);
//   //GCAcopyDCToMRI(gca, mri_dst);
//   //strcpy(mri_dst->fname,gca_fname); // copy gca name
//   LTA * lta = LTAalloc(1,mri_src);
//
//   if (!stricmp(fname_out+strlen(fname_out)-3, "XFM"))
//   {
//     cout << "converting xform to RAS..." << endl ;
//     cout << "initial:" << endl ;
//
//     MatrixPrint(stdout, Md.first) ;
//     lta->xforms[0].m_L = MRIvoxelXformToRasXform (mri_src, mri_dst, Md.first, lta->xforms[0].m_L) ;
//     cout << "final:" << endl ;
//     MatrixPrint(stdout,lta->xforms[0].m_L) ;
//     lta->type = LINEAR_RAS_TO_RAS ;
//   }
//   else
//   {
//     lta->xforms[0].m_L = MatrixCopy(Md.first, lta->xforms[0].m_L) ;
//     lta->type = LINEAR_VOX_TO_VOX ;
//   }
//
//   // add src and dst info
//     getVolGeom(mri_src, &lta->xforms[0].src);
//     getVolGeom(mri_dst, &lta->xforms[0].dst);
//     LTAwriteEx(lta, fname_out) ;
//
//   ///////////////////////////////////////////// end of writing transform
//
//    if (mri_src)
//     MRIfree(&mri_src) ;
//    if (mri_dst)
//     MRIfree(&mri_dst) ;
//    if (Md.first)
//     MatrixFree(&Md.first) ;
//
//
//   ///////////////////////////////////////////////////////////////
//   msec = TimerStop(&start) ;
//   seconds = nint((float)msec/1000.0f) ;
//   minutes = seconds / 60 ;
//   seconds = seconds % 60 ;
//   cout << "registration took "<<minutes<<" minutes and "<<seconds<<" seconds." << endl;
//   if (diag_fp)
//     fclose(diag_fp) ;
//   exit(0) ;
//   return(0) ;
//
// }

/*----------------------------------------------------------------------
  ----------------------------------------------------------------------*/
static void printUsage(void)
{
//  outputHelp("mri_robust_register");
  outputHelpMemory(mri_robust_register_help_xml,mri_robust_register_help_xml_len);

#ifdef GREGT
  cout << endl << endl;
  cout << "Usage: mri_robust_register <required arguments>" << endl <<endl;

  cout << "Short description: finds translation and rotation (using robust statistics)" << endl << endl;

  cout << "Required arguments" << endl << endl ;
  cout << "  --mov srcvol.mgz       movable/source volume to be registered" << endl;
  cout << "  --dst dstvol.mgz       destination/target volume for the registration" << endl;
  cout << "  --lta regfile.lta      output registration file" << endl;
	cout << "  Either --satit or --sat <real> (if not --leastsquares) for sensitivity" << endl << endl;

  cout << "Optional arguments" << endl << endl;
  cout << "  --warp outvol.mgz      apply final xform to source, write to outvol.mgz" << endl;
  cout << "  --weights wvol.mgz     output weights (in target space) as wvol.mgz" << endl;

  cout << "  --halfmov hm.mgz       outputs half-way mov (mapped to halfway space)" << endl;
  cout << "  --halfdst hd.mgz       outputs half-way dst (mapped to halfway space)" << endl;
  cout << "  --halfweights hw.mgz   outputs half-way weights (mapped to halfway space)" << endl;
  cout << "  --halfmovlta hm.lta    outputs transform from mov to half-way space" << endl;
  cout << "  --halfdstlta hd.lta    outputs transform from dst to half-way space" << endl;

//  cout << "  -A, --affine (testmode)    find 12 parameter affine xform (default is 6-rigid)" << endl;
  cout << "  --iscale               estimate intensity scale factor (default no)" << endl;
  cout << "                            !!Highly recommended for unnormalized images!!" << endl;
  cout << "  --iscaleout <str>      output file for iscale value (only with --iscale)" << endl;
	cout << "                            default: <regfile>-intensitiy.txt (regfile: --lta regifle.lta)" << endl;
  cout << "  --transonly            find 3 parameter translation only" << endl;
  cout << "  --transform lta        use initial transform lta on source ('id'=identity)" << endl;
  cout << "                            default is align center (using moments)" << endl;
  cout << "  --initorient           use moments for orientation init. (default false)" << endl;
  cout << "                            (recommended for stripped brains, but not with" << endl;
  cout << "                             with full head images with different cropping)"<<endl;
	cout << "  --noinit               skip transform init, default: transl. of centers" << endl;
  cout << "  --vox2vox              output VOX2VOX lta file (default is RAS2RAS)" << endl;
  cout << "  --leastsquares         use least squares instead of robust M-estimator" << endl;
  cout << "  --maxit <#>            iterate max # times on each resolution (default "<<P.iterate<<")"  << endl;
  cout << "  --highit <#>           iterate max # times on highest resol. (default "<<P.iterate<<")"  << endl;
  cout << "  --epsit <real>         stop iterations when below <real> (default "<<P.epsit <<")" << endl;
  cout << "  --nomulti              work on highest resolution (no multiscale)" << endl;
  cout << "  --sat <real>           set outlier sensitivity explicitly (e.g. '--sat 4.685' )" << endl;
	cout << "                             higher values mean less sensitivity" << endl;
//	cout << "                             default: automatically determine sat for head scans" << endl;
  cout << "  --satit                auto-detect good sensitivity (for head scans)" << endl;
	cout << "  --wlimit <real>        sets maximal outlier limit in satit (default "<<P.wlimit<<")" << endl;
  cout << "  --subsample <#>        subsample if dim > # on all axes (default no subs.)" << endl;
  cout << "  --doubleprec           double precision (default: float) for intensities (!!memory!!)" << endl;
  cout << "  --maskmov mask.mgz     mask mov/src with mask.mgz" << endl;
  cout << "  --maskdst mask.mgz     mask dst/target with mask.mgz" << endl;
//  cout << "  --conform              conform output volumes 1mm uchar vox (256^3)" << endl; // not implemented
  cout << "  --floattype            use float intensities (default keep input type)" << endl; 
  cout << "  --debug                create debug hw-images (default: no debug files)" << endl;
  cout << "  --verbose              0 quiet, 1 normal (default), 2 detail" << endl;
//  cout << "      --test i mri         perform test number i on mri volume" << endl;

  cout << endl;
  cout << endl;
  cout << "Description:" << endl;
	cout << "This program symmetrically aligns two volumes. It uses a method based on robust statistics to detect outliers and removes them from the registration. This leads to highly accurate registrations even with local changes in the image (e.g. jaw movement). The main purpose is to find the rigid registration (translation, rotation) of longitudinal data, but the method can be used to rigidly align different images. An additional optional intensity scale parameter can be used to adjust for global intensity differences. The extension to affine registration is being tested."<<endl;
  cout << endl;
  cout << "If the registration fails: " << endl;
	cout << "The registration can fail because of several reasons, most likeley due to large intensity differences or non-linear differences in the image. You can try:"<< endl;
	cout << " * Switch on intensity scaling (--iscale)." << endl;
	cout << " * When specifying a manual saturation (--sat) too many voxels might be considered outlier early in the process. You can check this by outputing the weights (--weights ow.mgz) and by looking at them in:" << endl;
	cout << "   > tkmedit -f dst.mgz -aux mov.mgz -overlay ow.mgz " << endl;
	cout << "   If most of the brain is labeled outlier, try to set the saturation to a higher value (eg. --sat 12) or use --satit to automatically determine a good sat value." << endl;
  cout << " * When using automatic saturation estimation (--satit) you can try specifying the sensitivity manually or twiddle around with --wlimit (which is around 0.16 by default). A lower wlimit should reduce the number of outlier voxels." << endl;
  cout << endl;
  cout << " Report bugs to: freesurfer@nmr.mgh.harvard.edu" << endl;
  cout << endl;

#endif

}

/*!
\fn void initRegistration(Registration & R, const Parameters & P)
\brief Initializes a Registration with Parameters (affine, iscale, transonly, leastsquares, sat and trans)
\param R  Registration to be initialized
\param P  Paramters for the initialization
*/
static void initRegistration(Registration & R, Parameters & P)
{
  R.setRigid(!P.affine);
  R.setIscale(P.iscale);
  R.setTransonly(P.transonly);
  R.setRobust(!P.leastsquares);
  R.setSaturation(P.sat);
  R.setVerbose(P.verbose); // set before debug, as debug sets its own verbose level
  R.setDebug(P.debug);
  R.setHighit(P.highit);
  R.setInitTransform(P.inittrans);
  R.setInitOrient(P.initorient);
  R.setDoublePrec(P.doubleprec);
  R.setWLimit(P.wlimit);
  R.setSymmetry(P.symmetry);
  //R.setOutputWeights(P.weights,P.weightsout);


  int pos = P.lta.rfind(".");
  if (pos > 0) R.setName(P.lta.substr(0,pos));
  else  R.setName(P.lta);

  if (P.subsamplesize > 0) R.setSubsamplesize(P.subsamplesize);

//   //////////////////////////////////////////////////////////////
//   // create a list of MRI volumes
//   //cout << "reading "<<ninputs<<" source (movable) volumes..."<< endl;
//   int ninputs = 1; // later we might want to load several frames
//   MRI* mri_src = NULL;
//   MRI* mri_tmp = NULL;
//   if (P.mri_mov) MRIfree(&P.mri_mov);
//   for (int i = 0 ; i < ninputs ; i++)
//   {
//     cout << "reading source '"<<P.mov<<"'..." << endl;
//     fflush(stdout) ;
// 
//     mri_tmp = MRIread(P.mov.c_str()) ;
//     if (!mri_tmp)
//     {
//       ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
//                 Progname, P.mov.c_str()) ;
//       //cerr << Progname << " could not open input volume " << P.mov << endl;
//       //exit(1);
//     }
// 
//     if (i == 0)
//     {
//       mri_src = MRIallocSequence(mri_tmp->width,
//                                  mri_tmp->height,
//                                  mri_tmp->depth,
//                                  mri_tmp->type,
//                                  ninputs) ;
//       MRIcopyHeader(mri_tmp, mri_src) ;
//       P.mri_mov = MRIallocSequence(mri_tmp->width,
//                                    mri_tmp->height,
//                                    mri_tmp->depth,
//                                    mri_tmp->type,
//                                    ninputs) ;
//       MRIcopyHeader(mri_tmp, P.mri_mov) ;
//     }
//     MRIcopyFrame(mri_tmp, P.mri_mov, 0, i) ; // store input in P.mri_mov
// 
//     if (P.maskmov != "") // work only on mri_src to init registration (not P.mri_mov)
//     {
//       MRI *mri_mask = MRIread(P.maskmov.c_str());
//       if (!mri_mask)
//         ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
//                   Progname, P.maskmov.c_str()) ;
//       MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
//       MRIfree(&mri_mask) ;
//     }
// 
//     MRIcopyFrame(mri_tmp, mri_src, 0, i) ;
//     MRIfree(&mri_tmp) ;
//   }
////   R.setSource(mri_src,P.fixvoxel,P.fixtype);
////   MRIfree(&mri_src);

  ///////////  read MRI Source //////////////////////////////////////////////////
  cout << endl;
  cout <<  "reading source '"<<P.mov<<"'..."<< endl ;
  fflush(stdout) ;

  MRI* mri_mov = MRIread(P.mov.c_str()) ;
  if (mri_mov == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not open MRI source %s.\n",
              Progname, P.mov.c_str()) ;
    //cerr << Progname << " could not open MRI Target " << P.mov << endl;
    //exit(1);
  }
	if (mri_mov->nframes != 1)
	{
    ErrorExit(ERROR_NOFILE, "%s: only pass single frame MRI source %s.\n",
              Progname, P.mov.c_str()) ;	
	}
  P.mri_mov = MRIcopy(mri_mov,P.mri_mov); // save dst mri

  if (P.maskmov != "")
  {
    MRI *mri_mask = MRIread(P.maskmov.c_str());
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                Progname, P.maskmov.c_str()) ;
    MRImask(mri_mov, mri_mask, mri_mov, 0, 0) ;
    MRIfree(&mri_mask) ;
  }

  ///////////  read MRI Target //////////////////////////////////////////////////
  cout <<  "reading target '"<<P.dst<<"'..."<< endl ;
  fflush(stdout) ;

  MRI* mri_dst = MRIread(P.dst.c_str()) ;
  if (mri_dst == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not open MRI target %s.\n",
              Progname, P.dst.c_str()) ;
    //cerr << Progname << " could not open MRI Target " << P.dst << endl;
    //exit(1);
  }
	if (mri_dst->nframes != 1)
	{
    ErrorExit(ERROR_NOFILE, "%s: only pass single frame MRI target %s.\n",
              Progname, P.dst.c_str()) ;	
	}
  P.mri_dst = MRIcopy(mri_dst,P.mri_dst); // save dst mri

  if (P.maskdst != "")
  {
    MRI *mri_mask = MRIread(P.maskdst.c_str());
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                Progname, P.maskdst.c_str()) ;
    MRImask(mri_dst, mri_mask, mri_dst, 0, 0) ;
    MRIfree(&mri_mask) ;
  }
	
  // Set initial transform //////////////////////////////////////////////////
  if (P.transform != "")
  {
    cout << endl << "reading initial transform '"<<P.transform<<"'..."<< endl;
    
//     // try to read simple text
//     bool st = true;
//     MATRIX* mi = MatrixAlloc(4,4,MATRIX_REAL);
//     ifstream f;
//     while (1==1) // fake while loop (to be run once)
//     {
//       string sin (P.transform);
//       if (sin == "id" || sin == "identity.nofile")
//       {
//         mi = MatrixIdentity(4,mi);
//         break;
//       }
// 
//       f.open(P.transform.c_str(),ios::in);
//       if (!f)
//       {
//         cerr <<" Read 4x4: could not open initial transform file " << P.transform <<endl;
//         st = false;
//         break;
//       }
// 
//       int row, col;
//       for (row = 1 ; row <= 4 ; row++)
//       {
//         string s;
//         getline(f,s);
//         istringstream s1(s);
//         for (col = 1 ; col <= 4 ; col++)
//         {
//           s1>> mi->rptr[row][col];
//         }
//         if (!s1.good())
//         {
//           st = false;
//           break;
//         }
//       }
//       break; // quit fake while loop
//     }
//     f.close();
// 
// 
//     if (st)
//     {
//       R.setMinit(MyMatrix::convertMATRIX2VNL(mi));
//     }
//     MatrixFree(&mi);
// 
// 
//     if (!st)
//     {
      // try to read other transform
      TRANSFORM * trans = TransformRead(P.transform.c_str());
      LTA* lta =  (LTA *)trans->xform ;
      if (!lta)
        ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",Progname, P.transform.c_str()) ;
      if (! lta->xforms[0].src.valid )
			{
			  cout << " WARNING: no source geometry (RAS) in transform, assuming movable !!!" << endl;
        getVolGeom(mri_mov, &lta->xforms[0].src);
			}
      if (! lta->xforms[0].dst.valid )
			{
			  cout << " WARNING: no target geometry (RAS) in transform, assuming destination !!!" << endl;
        getVolGeom(mri_dst, &lta->xforms[0].dst);
			}
      lta = LTAchangeType(lta,LINEAR_VOX_TO_VOX);
      if (lta->type!=LINEAR_VOX_TO_VOX)
      {
        ErrorExit(ERROR_BADFILE, "%s: must be LINEAR_VOX_TO_VOX (=0), but %d", Progname, P.transform.c_str(), lta->type) ;
      }
      R.setMinit(MyMatrix::convertMATRIX2VNL(lta->xforms[0].m_L));
      //if (P.debug) // apply init transform to input source image directly
      //{
      //  MRI * mri_tmp = LTAtransform(mri_mov,NULL, lta);
	    //  string fn = R.getName() + "-source-init.mgz";
      //  MRIwrite(mri_tmp,fn.c_str());
      //  MRIfree(&mri_tmp);
      //}
//    }
  }
	
	cout << endl;
	
  // now actually set source and target (and possibly reslice):
  R.setSourceAndTarget(mri_mov,mri_dst,!P.floattype);	
  MRIfree(&mri_mov);
  MRIfree(&mri_dst);

}

/*!
\fn int parseNextCommand(int argc, char **argv)
\brief Parses the command-line for next command
\param   argc  number of command line arguments
\param   argv  pointer to a character pointer
\param      P  reference to parameters
\returns       number of used arguments for this command
*/
static int parseNextCommand(int argc, char *argv[], Parameters & P)
{
  int  nargs = 0 ;
  char *option ;

  option = argv[0] + 1 ;                     // remove '-'
  if (option[0] == '-') option = option +1;  // remove second '-'
  StrUpper(option) ;

  //cout << " option: " << option << endl;

  if (!strcmp(option, "MOV") ||  !strcmp(option, "M") )
  {
    P.mov = string(argv[1]);
    nargs = 1;
    cout << "--mov: Using "<< P.mov << " as movable/source volume." << endl;
  }
  else if (!strcmp(option, "DST") || !strcmp(option, "D") )
  {
    P.dst = string(argv[1]);
    nargs = 1;
    cout << "--dst: Using "<< P.dst << " as target volume." << endl;
  }
  else if (!strcmp(option, "LTA")   )
  {
    P.lta = string(argv[1]);
    nargs = 1;
    cout << "--lta: Output transform as "<< P.lta << " . " << endl;
  }
  else if (!strcmp(option, "VOX2VOX")   )
  {
    P.lta_vox2vox = true;
    cout << "--vox2vox: Output transform as VOX2VOX. " << endl;
  }
  else if (!strcmp(option, "AFFINE") || !strcmp(option, "A") )
  {
    P.affine = true;
    cout << "--affine: Enableing affine transform!" << endl;
  }
  else if (!strcmp(option, "ISCALE") || !strcmp(option, "I") )
  {
    P.iscale = true;
    cout << "--iscale: Enableing intensity scaling!" << endl;
  }
  else if (!strcmp(option, "TRANSONLY"))
  {
    P.transonly = true;
    cout << "--transonly: Using only translation!" << endl;
  }
  else if (!strcmp(option, "TRANSFORM") || !strcmp(option, "T") )
  {
    P.transform = string(argv[1]);
    nargs = 1;
    cout << "--transform: Using previously computed initial transform: "<< argv[1] << endl;
  }
  else if (!strcmp(option, "INITORIENT"))
  {
    P.initorient = true;
    cout << "--initorient: Using moments for initial orientation!" << endl;
  }
  else if (!strcmp(option, "NOINIT"))
  {
    P.inittrans = false;
    cout << "--noinit: Skipping init of transform !" << endl;
  }
  else if (!strcmp(option, "LEASTSQUARES") || !strcmp(option, "L")  )
  {
    P.leastsquares = true;
    cout << "--leastsquares: Using standard least squares (non-robust)!" << endl;
  }
  else if (!strcmp(option, "MAXIT")  )
  {
    P.iterate = atoi(argv[1]);
    nargs = 1 ;
    cout << "--maxit: Performing maximal " << P.iterate << " iterations on each resolution" << endl;
  }
  else if (!strcmp(option, "HIGHIT")  )
  {
    P.highit = atoi(argv[1]);
    nargs = 1 ;
    cout << "--highit: Performing maximal " << P.highit << " iterations on highest resolution" << endl;
  }
  else if (!strcmp(option, "EPSIT") )
  {
    P.epsit = atof(argv[1]);
    nargs = 1 ;
    cout << "--epsit: Stop iterations when change is less than " << P.epsit << " . " << endl;
  }
  else if (!strcmp(option, "NOMULTI") )
  {
    P.nomulti = true;
    nargs = 0 ;
    cout << "--nomulti: Will work on highest resolution only (nomulti)!" << endl;
  }
  else if (!strcmp(option, "SAT")  )
  {
    P.sat = atof(argv[1]);
    nargs = 1 ;
    cout << "--sat: Using saturation " << P.sat << " in M-estimator!" << endl;
  }
  else if (!strcmp(option, "WLIMIT")  )
  {
    P.wlimit = atof(argv[1]);
    nargs = 1 ;
    cout << "--wlimit: Using wlimit in satit " << P.wlimit <<  endl;
  }
  else if (!strcmp(option, "SUBSAMPLE") )
  {
    P.subsamplesize = atoi(argv[1]);
    nargs = 1 ;
    if (P.subsamplesize >= 0) cout << "--subsample: Will subsample if size is larger than " << P.subsamplesize << " on all axes!" << endl;
    else cout << "--subsample -1: Will not subsample on any scale!" << endl;
  }
  else if (!strcmp(option, "SATIT") )
  {
    P.satit = true;
    nargs = 0 ;
    cout << "--satit: Will iterate with different SAT to ensure outliers below wlimit!" << endl;
  }
  else if (!strcmp(option, "SATEST") ) // old  remove
  {
    P.satest = true;
    nargs = 0 ;
    cout << "--satest: Will estimate SAT (never really tested, use --satit instead!)" << endl;
  }
  else if (!strcmp(option, "SATEST") ) // never reached???  - old remove
  {
    P.dosatest = true;
    nargs = 0 ;
    cout << "--satest: Trying to estimate SAT value!" << endl;
  }
  else if (!strcmp(option, "DOUBLEPREC") )
  {
    P.doubleprec = true;
    nargs = 0 ;
    cout << "--doubleprec: Will perform algorithm with double precision (higher mem usage)!" << endl;
  }
  else if (!strcmp(option, "DEBUG") )
  {
    P.debug = 1;
    nargs = 0 ;
    cout << "--debug: Will output debug info and files!" << endl;
  }
  else if (!strcmp(option, "VERBOSE") )
  {
    P.verbose = atoi(argv[1]);
    nargs = 1 ;
    cout << "--verbose: Will use verbose level : " << P.verbose << endl;
  }
  else if (!strcmp(option, "WEIGHTS") )
  {
    P.weightsout = string(argv[1]);
    nargs = 1 ;
    cout << "--weights: Will output weights transformed to target space as "<<P.weightsout<<" !" << endl;
  }
  else if (!strcmp(option, "WARP") || !strcmp(option, "W") )
  {
    P.warp = true;
    P.warpout = string(argv[1]);
    nargs = 1 ;
    cout << "--warp: Will save warped source as "<<P.warpout <<" !" << endl;
  }
  else if (!strcmp(option, "HALFMOV") )
  {
    P.halfmov = string(argv[1]);
    nargs = 1 ;
    cout << "--halfmov: Will output final half way MOV !" << endl;
  }
  else if (!strcmp(option, "HALFDST") )
  {
    P.halfdst = string(argv[1]);
    nargs = 1 ;
    cout << "--halfdst: Will output final half way DST !" << endl;
  }
  else if (!strcmp(option, "HALFWEIGHTS") )
  {
    P.halfweights = string(argv[1]);
    nargs = 1 ;
    cout << "--halfweights: Will output half way WEIGHTS from last step to " <<P.halfweights<<" !" << endl;
  }
  else if (!strcmp(option, "HALFMOVLTA") )
  {
    P.halfmovlta = string(argv[1]);
    nargs = 1 ;
    cout << "--halfmovlta: Will output half way transform (mov) " <<P.halfmovlta << " !" << endl;
  }
  else if (!strcmp(option, "HALFDSTLTA") )
  {
    P.halfdstlta = string(argv[1]);
    nargs = 1 ;
    cout << "--halfdstlta: Will output half way transform (dst) " <<P.halfdstlta << " !" << endl;
  }
  else if (!strcmp(option, "MASKMOV") )
  {
    P.maskmov = string(argv[1]);
    nargs = 1 ;
    cout << "--maskmov: Will apply "<<P.maskmov <<" to mask mov/src !" << endl;
  }
  else if (!strcmp(option, "MASKDST") )
  {
    P.maskdst = string(argv[1]);
    nargs = 1 ;
    cout << "--maskdst: Will apply "<<P.maskdst <<" to mask dst/target !" << endl;
  }
  else if (!strcmp(option, "TEST"))
  {
    cout << "--test: TEST-MODE " << endl;
    Registration R;
    R.testRobust(argv[2], atoi(argv[1]));
    nargs = 2 ;
    exit(0);
  }
  else if (!strcmp(option, "CONFORM") )
  {
    P.conform = true;
    nargs = 0 ;
    cout << "--conform: Will conform images to 256^3 and voxels to 1mm!" << endl;
  }
  else if (!strcmp(option, "FLOATTYPE") )
  {
    P.floattype = true;
    nargs = 0 ;
    cout << "--floattype: Use float images internally (independent of input)!" << endl;
  }
  else if (!strcmp(option, "ONEMINUSW") )
  {
    P.oneminusweights = true;
    nargs = 0 ;
    cout << "--oneminusw: Will output 1-weights!" << endl;
  }
  else if (!strcmp(option, "NOSYM") )
  {
    P.symmetry = false;
    nargs = 0 ;
    cout << "--nosym: Will resample source to target (no half-way space)!" << endl;
  }
  else if (!strcmp(option, "ISCALEOUT") )
  {
    P.iscaleout = string(argv[1]);
    nargs = 1 ;
    P.iscale = true;
    cout << "--iscaleout: Will do --iscale and ouput intensity scale to "<<P.iscaleout <<  endl;
  }
  else
  {
    cerr << endl << endl << "ERROR: Option: " << argv[0] << " unknown !! " << endl << endl;
    exit(1);
  }

  fflush(stdout);

  return(nargs) ;
}

/*!
\fn int parseCommandLine(int argc, char **argv)
\brief Parses the command-line
\param   argc  number of command line arguments
\param   argv  pointer to a character pointer
\param      P  reference to parameters
\returns       if all necessary parameters were set
*/
static bool parseCommandLine(int argc, char *argv[], Parameters & P)
{
  int nargs;

  for ( ; argc > 0 && ISOPTION(*argv[0]) ; argc--, argv++)
  {
    nargs = parseNextCommand(argc, argv,P) ;
    argc -= nargs ;
    argv += nargs ;
  }

  bool test1 = ( P.mov != "" && P.dst != "" && P.lta != "" );
	if (!test1)
	{
	  cerr << endl << "Please specify --mov --dst and --lta !  "<< endl;
	}
	bool test2 = ( P.satit || P.sat > 0 || P.leastsquares );
	if (!test2)
	{
	  cerr << endl << "Please specify either --satit or --sat <float> !  "<< endl;
	}
	bool test3 = ( P.iscaleout == "" || P.iscale);
	if (!test3)
	{
	  cerr << endl << "Please spedify --iscale together with --iscaleout to compute and output global intensity scaling! " << endl;
	}
  return (test1 && test2 && test3);
}
