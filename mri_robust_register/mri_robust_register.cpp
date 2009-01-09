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
 *    $Date: 2009/01/09 19:10:08 $
 *    $Revision: 1.2 $
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

#include "Registration.h"

// all other software are all in "C"
#ifdef __cplusplus
extern "C" {
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

//#define SAT 4.685 // this is suggested for gaussian noise
#define SAT 20
#define SSAMPLE -1

struct Parameters
{
  string mov;
  string dst;
  string lta;
  bool   affine;
  bool   iscale; 
  bool   transonly;
  string transform;
  bool   leastsquares;
  int    iterate;
  double sat;
  bool   warp;
  string warpout;
  int    subsamplesize;
  MRI*   mri_src;
  MRI*   mri_dst;
};
static struct Parameters P = { "","","",false,false,false,"",false,-1,SAT,false,"",SSAMPLE,NULL,NULL};


static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[],Parameters & P) ;
static void initRegistration(Registration & R, Parameters & P) ;

static char vcid[] = "$Id: mri_robust_register.cpp,v 1.2 2009/01/09 19:10:08 mreuter Exp $";
char *Progname = NULL;

//static MORPH_PARMS  parms ;
//static FILE *diag_fp = NULL ;


using namespace std;


int main(int argc, char *argv[])
{
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
  Registration R;
  initRegistration(R,P);
  
  // compute Alignment
  std::pair <MATRIX*, double> Md;
  if (P.iterate > 0)  Md = R.computeIterativeRegistration(P.iterate);
  else Md = R.computeMultiresRegistration();

  // Print results:
  cout << "final transform:\n" ;
  MatrixPrintFmt(stdout,"% 2.8f",Md.first);
  if (R.isIscale()) cout << " iscale: " << Md.second << endl;
  cout << endl ;

  // maybe warp source to target:
   char reg[STRLEN];
   strcpy(reg, P.lta.c_str());
   if (P.warp)
   {
//   // construct base name of P.reg
//      char foutbase[STRLEN];
//      FileNameOnly(reg, foutbase) ;
//      FileNameRemoveExtension(foutbase, foutbase) ;
//       string outfn (foutbase);
//       outfn += ".mgz";
      R.warpSource(P.warpout);
      cout << "To check results, run:" << endl;
      cout << "  tkmedit -f "<< P.dst <<" -aux " << P.warpout << endl;
   }


  // writing transform section here
  cout << "writing output transformation to "<<P.lta <<" ..." << endl;
  LTA * lta = LTAalloc(1,P.mri_src);
  
//  if (!stricmp(fname_out+strlen(fname_out)-3, "XFM"))
  if ( (P.lta.substr(P.lta.length() - 3) == "XFM"))
  {
    cout << "converting xform to RAS..." << endl ;
    cout << "initial:" << endl ;
    
    MatrixPrint(stdout, Md.first) ;
    lta->xforms[0].m_L = MRIvoxelXformToRasXform (P.mri_src, P.mri_dst, Md.first, lta->xforms[0].m_L) ;
    cout << "final:" << endl ;
    MatrixPrint(stdout,lta->xforms[0].m_L) ;
    lta->type = LINEAR_RAS_TO_RAS ;
  } 
  else
  {
    lta->xforms[0].m_L = MatrixCopy(Md.first, lta->xforms[0].m_L) ;
    lta->type = LINEAR_VOX_TO_VOX ;
  }

  // add src and dst info
    getVolGeom(P.mri_src, &lta->xforms[0].src);
    getVolGeom(P.mri_dst, &lta->xforms[0].dst);
    LTAwriteEx(lta, reg) ;

  ///////////////////////////////////////////// end of writing transform

   if (Md.first)
    MatrixFree(&Md.first) ;


  ///////////////////////////////////////////////////////////////
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  cout << "registration took "<<minutes<<" minutes and "<<seconds<<" seconds." << endl;
  //if (diag_fp) fclose(diag_fp) ;
  exit(0) ;
  return(0) ; 
}

// int main(int argc, char *argv[])
// {
// //   testRegression();
// 
// //   testRobust(argv[1],atoi(argv[2]));
// //     MRI *mri_kernel ;
// // 
// //     mri_kernel = MRIgaussian1d(1.08, 5) ;
// //     
// //   for (int c=0; c < mri_kernel->width; c++)
// //   for (int r=0; r < mri_kernel->height; r++)
// //   for (int s=0; s < mri_kernel->depth; s++)
// //   for (int f=0; f < mri_kernel->nframes ; f++)
// //   {
// //   	cout << " c: " << c << "  r: " << r << "  s: " << s << "  f: " << f << " = " ;
// // 	cout << MRIgetVoxVal(mri_kernel, c, r, s, f) << endl;
// //   }
// //   exit(1);
// // 
// //   MATRIX * p = MatrixAlloc(6,1,MATRIX_REAL);
// //   *MATRIX_RELT(p, 1, 1) = -100;
// //   *MATRIX_RELT(p, 2, 1) = 2;
// //   *MATRIX_RELT(p, 3, 1) = 800;
// //   *MATRIX_RELT(p, 4, 1) = 5;
// //   *MATRIX_RELT(p, 5, 1) = 4;
// //   *MATRIX_RELT(p, 6, 1) = 3;
// //   MATRIX * m = param2matrixt(p);
// //   MatrixPrint(stdout,m);
// //   MATRIX *m2 = param2matrix(p);
// //   cout << endl;
// //   MatrixPrint(stdout,m2);
// // 
// // //MatrixPrint(stdout, p);
// // //cout << VectorMedian(p) << endl;
// // //double ss = getSigmaMAD(p);
// // //cout << endl << ss << endl;
// //   exit(1);
// //   
// 
// //  char         *in_Sfname, *in_Tfname, *out_fname,fname[STRLEN];
// //  MRI          *mri_inS, *mri_inT, *mri_tmp, *mri_dst;
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
static void printUsage(void) {
  cout << endl << endl;
  cout << "Usage: mri_robust_register <required arguments>" << endl <<endl;

  cout << "Short description: finds translation and rotation (using robust statistics)" << endl << endl;

  cout << "Required arguments" << endl << endl ;
  cout << "  -M, --mov srcvol.mgz       movable/source volume to be registered" << endl;
  cout << "  -D, --dst dstvol.mgz       destination/target volume for the registration" << endl;
  cout << "      --lta regfile.lta      output registration file" << endl << endl;

  cout << "Optional arguments" << endl << endl;
  cout << "  -A, --affine               find 12 parameter affine xform (default is 6-rigid)" << endl;
  cout << "  -I, --iscale               allow also intensity scaling" << endl;
  cout << "      --transonly            find 3 parameter translation" << endl;
  cout << "  -T, --transform lta        use transform lta on source, type 'id' for identity" << endl;
  cout << "                                default is RAS-mov to VOX to RAS-dst" << endl;
  cout << "  -L, --leastsquares         use least squares instead of robust M-estimator" << endl;
  cout << "  -N, --iterate <#>          iterate # times on highest res. (no multiscale)" << endl;
  cout << "      --sat <float>          set saturation for robust estimator (default "<<SAT<<")" << endl;
  cout << "      --subsample <#>        subsample if dim > # on all axes" << endl;
  cout << "                                (default never subsample)" << endl;
//  cout << "                                set to -1 to avoid subsampling on all scales" << endl;
  cout << "  -W, --warp outvol.mgz      apply final xform on source, write to outvol.mgz" << endl;
//  cout << "      --test i mri         perform test number i on mri volume" << endl;

  cout << endl;
  cout << "Mandatory or optional arguments to long options are also mandatory or optional for any corresponding short options." << endl;
  cout << endl;
  
  cout << "Report bugs to: analysis-bugs@nmr.mgh.harvard.edu" << endl;

  
  /*printf("  -dist distance\n");
  printf("  -nomap\n");
  printf("  -flash\n");
  printf("  -mask mask\n");
  printf("  -skull\n");
  printf("  -uns nbrspacing\n");
  printf("  -diag diagfile\n");
  printf("  -debug_voxel x y z\n");
  printf("  -debug_label label\n");
  printf("  -tr TR\n");
  printf("  -te TE\n");
  printf("  -alpha alpha\n");
  printf("  -example T1 seg\n");
  printf("  -samples fname\n");
  printf("  -fsamples fname\n");
  printf("  -nsamples fname\n");
  printf("  -contrast\n");
  printf("  -flash_parms parameter\n");
  printf("  -transonly \n");
  printf("  -write_mean fname\n");
  printf("  -prior min_prior\n");
  printf("  -spacing max_spacing\n");
  printf("  -scales nscales\n");
  printf("  -novar\n");
  printf("  -dt dt\n");
  printf("  -tol tol\n");
  printf("  -center\n");
  printf("  -noscale\n");
  printf("  -noiscale\n");
  printf("  -num num_xforms\n");
  printf("  -area area\n");
  printf("  -nlarea nlarea\n");
  printf("  -levels levels\n");
  printf("  -intensity intensity\n");
  printf("  -reduce nreductions\n");
  printf("  -nsamples nsamples\n");
  printf("  -norm fname\n");
  printf("  -trans max_trans\n");
  printf("  -steps max_angles\n");
  printf("  -l xform long_reg\n");
  printf("  -f controlpoints\n");
  printf("  -d tx ty tz\n");
  printf("  -r rx ry rz\n");
  printf("  -t xform\n");
  printf("  -b blur_sigma\n");
  printf("  -v diagno\n");
  printf("  -s max_angles\n");
  printf("  -max_angle max_angle in radians (def=15 deg)\n");
  printf("  -n niters\n");
  printf("  -w write_iters\n");
  printf("  -p ctl_point_pct : use top pct percent wm points as control points\n");
  printf("  -m momentum\n");*/

 
 cout << endl;

}

// /*!
// \fn void initRegistration(Registration & R, const Parameters & P)
// \brief Converts uchar MRI to float 0..1
// \param mri_uchar MRI volume of type uchar (0..255)
// \returns convered float MRI volume (0..1)
// */
// static MRI* convertFloat01(MRI * mri_uchar)
// {
//   int x,y,z;
//   if (mri_uchar->type != MRI_UCHAR)
//   {
//   	cerr << " convertFloat01 failed: input must be of type UCHAR" << endl;
// 	assert(mri_uchar->type == MRI_UCHAR);
//   }
//   
//   MRI* mri_float = MRIalloc( mri_uchar->width,  mri_uchar->height, mri_uchar->depth,MRI_FLOAT) ;
//   MRIcopyHeader(mri_uchar, mri_float) ;
//   mri_float->type = MRI_FLOAT;
//   
//   for (z = 0 ; z < mri_uchar->depth ; z++)
//   for (y = 0 ; y < mri_uchar->height; y++)
//   for (x = 0 ; x < mri_uchar->width ; x++)
//   {
//      MRIsetVoxVal(mri_float, x,y,z,0, MRIgetVoxVal(mri_uchar,x,y,z,0)/255.0);
//      //if (MRIgetVoxVal(mri_uchar,x,y,z,0) > 150)
//      //cout << " uchar: " << MRIgetVoxVal(mri_uchar,x,y,z,0) << "   float: " << MRIgetVoxVal(mri_float,x,y,z,0) << endl;
// 
//   }
//   return mri_float;
// }

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
   
   if (P.subsamplesize > 0) R.setSubsamplesize(P.subsamplesize);

   // Set initial transform
   if (P.transform != "")
   {
     cout << endl << "reading initial transform '"<<P.transform<<"'..."<< endl;;
       // try to read simple text 
       bool st = true;
       MATRIX* mi = MatrixAlloc(4,4,MATRIX_REAL);
       ifstream f;
       while (1==1) // fake while loop (to be run once)
       {
         string sin (P.transform);
	 if (sin == "id") 
	 {
	   mi = MatrixIdentity(4,mi);
	   break;
	  }
       
	f.open(P.transform.c_str(),ios::in);
	if(!f) 
        {
	   cerr <<" Read 4x4: could not open initial transform file " << P.transform <<endl;
	   st = false;
	   break;
	}

	int row, col;
	  for (row = 1 ; row <= 4 ; row++)
	  {
	    string s;
	    getline(f,s);
	    istringstream s1(s);
	    for (col = 1 ; col <= 4 ; col++)
	    {
		 s1>> mi->rptr[row][col];
	    }
	    if (!s1.good()) 
	    {
	    	st = false;
	    	break;
	    }
	  }
	  break; // quit fake while loop
	}
	f.close();
       
       
       if (st)
       {
         R.setMinit(mi);
       }
       MatrixFree(&mi);
       
       
       if (!st)
       {
       
         // try to read other transform
         char ctrans[STRLEN];
         strcpy(ctrans, P.transform.c_str());
         TRANSFORM * trans = TransformRead(ctrans);
         LTA* lta =  (LTA *)trans->xform ;
         if (!lta)
         ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",Progname, ctrans) ;
          if (lta->type!=LINEAR_VOX_TO_VOX)
          {
            ErrorExit(ERROR_BADFILE, "%s: must be LINEAR_VOX_TO_VOX (=0), but %d", Progname, ctrans, lta->type) ;
	    exit(1);
	  }
        R.setMinit(lta->xforms[0].m_L);
       }  
    } 
    
   //////////////////////////////////////////////////////////////
  // create a list of MRI volumes
  //cout << "reading "<<ninputs<<" source (movable) volumes..."<< endl;
  int ninputs = 1; // later we might want to load several frames
  MRI* mri_src = NULL;
  MRI* mri_tmp = NULL;
  for (int i = 0 ; i < ninputs ; i++)
  {
    cout << "reading source '"<<P.mov<<"'..." << endl;
    fflush(stdout) ;
   char mov[STRLEN];
   strcpy(mov, P.mov.c_str());
    
    mri_tmp = MRIread(mov) ;
    if (!mri_tmp)
    {
      cerr << Progname << " could not open input volume " << P.mov << endl;
      exit(1);
    }

//    TRs[i] = mri_tmp->tr ;
//    fas[i] = mri_tmp->flip_angle ;
//    TEs[i] = mri_tmp->te ;

//    if (mask_fname) {
//      MRI *mri_mask ;
//
//      mri_mask = MRIread(mask_fname) ;
//      if (!mri_mask)
//        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
//                  Progname, mask_fname) ;
//      MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
//      MRIfree(&mri_mask) ;
//    }
    if (i == 0)
    {
      mri_src = MRIallocSequence(mri_tmp->width,
                                mri_tmp->height,
                                mri_tmp->depth,
                                mri_tmp->type,
                                ninputs) ;
      MRIcopyHeader(mri_tmp, mri_src) ;
    }
    MRIcopyFrame(mri_tmp, mri_src, 0, i) ;
    MRIfree(&mri_tmp) ;
  }
//   if (mri_src->type == MRI_UCHAR) 
//   {
//      mri_tmp = mri_src;
//      mri_src = convertFloat01(mri_tmp);
//      MRIfree(&mri_tmp);
//   }
//   assert (mri_src->type == MRI_FLOAT);
  R.setSource(mri_src);
  P.mri_src = mri_src;


  ///////////  read MRI Target //////////////////////////////////////////////////
  cout <<  "reading target '"<<P.dst<<"'..."<< endl ;
  fflush(stdout) ;
   char dst[STRLEN];
   strcpy(dst, P.dst.c_str());
  MRI* mri_dst = MRIread(dst) ;
  if (mri_dst == NULL)
  {
     cerr << Progname << " could not open MRI Target " << P.dst << endl;
     exit(1);
  }
  cout << endl;
//   if (mri_dst->type == MRI_UCHAR) 
//   {
//      mri_tmp = mri_dst;
//      mri_dst = convertFloat01(mri_tmp);
//      MRIfree(&mri_tmp);
//   }
//   assert (mri_dst->type == MRI_FLOAT);
  R.setTarget(mri_dst);
  P.mri_dst = mri_dst;
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
  if (option[0] == '-') option = option +1;  // remove '-' 
  StrUpper(option) ;
  if (!strcmp(option, "MOV") || *option =='M' )
  {
     P.mov = string(argv[1]);
     nargs = 1;
     cout << "Using "<< P.mov << " as movable/source volume." << endl;
  }
  else if (!strcmp(option, "DST") || *option =='D' )
  {
     P.dst = string(argv[1]);
     nargs = 1;
     cout << "Using "<< P.dst << " as target volume." << endl;
  }
  else if (!strcmp(option, "LTA")   )
  {
     P.lta = string(argv[1]);
     nargs = 1;
     cout << "Output transform as "<< P.lta << " . " << endl;
  }
  else if (!strcmp(option, "AFFINE") || *option =='A' )
  {
     P.affine = true;
     cout << "Enableing affine transform!" << endl;
  }
  else if (!strcmp(option, "ISCALE") || *option =='I' )
  {
     P.iscale = true;
     cout << "Enableing intensity scaling!" << endl;
  }
  else if (!strcmp(option, "TRANSONLY"))
  {
     P.transonly = true;
     cout << "Using only translation!" << endl;
  }
  else if (!strcmp(option, "TRANSFORM") || *option =='T' )
  {
     P.transform = string(argv[1]);
     nargs = 1;
     cout << "Using previously computed initial transform: "<< argv[1] << endl;
  }
  else if (!strcmp(option, "LEASTSQUARES") || *option =='L'  )
  {
     P.leastsquares = true;
     cout << "Using standard least squares (non-robust)!" << endl;
  }
  else if (!strcmp(option, "ITERATE") || *option =='N' )
  {
       P.iterate = atoi(argv[1]);
       nargs = 1 ;
       cout << "Performing " << P.iterate << " iterations on highest resolutin (no multires)" << endl;
  }
  else if (!strcmp(option, "SAT")  )
  {
     P.sat = atof(argv[1]);
     nargs = 1 ;
     cout << "Using saturation " << P.sat << " in M-estimator!" << endl;
  }
  else if (!strcmp(option, "SUBSAMPLE") )
  {
     P.subsamplesize = atoi(argv[1]);
     nargs = 1 ;
     if (P.subsamplesize >= 0) cout << "Will subsample if size is larger than " << P.subsamplesize << " on all axes!" << endl;
     else cout << "Will not subsample on any scale!" << endl;
  }
  else if (!strcmp(option, "WARP") || *option =='W' )
  {
     P.warp = true;
     P.warpout = string(argv[1]);
     nargs = 1 ;
     cout << "Will save warped source as "<<P.warpout <<"!" << endl;
  }
  else if (!strcmp(option, "TEST"))
  {
     Registration R;
     R.testRobust(argv[2], atoi(argv[1]));
     nargs = 2 ;
     exit(0);
  }
  else
    cout << "Option: " << argv[0] << " unknown !! " << endl;

  
/*  if (!strcmp(option, "DIST") || !strcmp(option, "DISTANCE")) {
    // seems like not used.
    parms.l_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_dist = %2.2f\n", parms.l_dist) ;
  } else if (!strcmp(option, "NOMAP")) {
    // seems not used
    nomap = 1 ;
  } else if (!stricmp(option, "FLASH")) {
    map_to_flash = 1 ;
    printf("using FLASH forward model to predict intensity values...\n") ;
  } else if (!stricmp(option, "MAX_ANGLE")) {
    MAX_ANGLE = RADIANS(atof(argv[2])) ;
    printf("using %2.2f deg as max angle for rotational search\n",
           DEGREES(MAX_ANGLE)) ;
    nargs = 1 ;
  } else if (!stricmp(option, "BABY")) {
    baby = 1 ;
    printf("using baby brain intensity model\n") ;
  } else if (!strcmp(option, "MASK")) {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  } else if (!strcmp(option, "SKULL")) {
    unknown_nbr_spacing = 5 ;
    printf("aligning to atlas containing skull, "
           "setting unknown_nbr_spacing = %d\n",
           unknown_nbr_spacing) ;
    skull = 1 ;
  } else if (!strcmp(option, "RIGID")) {
    rigid = 1 ;
    printf("constraining transform to be rigid\n") ;
  } else if (!strcmp(option, "UNS")) {
    unknown_nbr_spacing = atoi(argv[2]) ;
    nargs = 1 ;
    printf("aligning to atlas containing skull, "
           "setting unknown_nbr_spacing = %d\n",
           unknown_nbr_spacing) ;
  }
  /////// debug options //////////////////////////////////
  else if (!strcmp(option, "DIAG")) {
    diag_fp = fopen(argv[2], "w") ;
    if (!diag_fp)
      ErrorExit
      (ERROR_NOFILE,
       "%s: could not open diag file %s for writing",
       Progname, argv[2]) ;
    printf("opening diag file %s for writing\n", argv[2]) ;
    nargs = 1 ;
  } else if (!strcmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!strcmp(option, "DEBUG_LABEL")) {
    Ggca_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging label %s (%d)\n",
           cma_label_to_name(Ggca_label), Ggca_label) ;
  }
  ////////// TR, TE, Alpha ////////////////////////////////
  else if (!strcmp(option, "TR")) {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
  } else if (!strcmp(option, "TE")) {
    TE = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TE=%2.1f msec\n", TE) ;
  } else if (!strcmp(option, "ALPHA")) {
    nargs = 1 ;
    alpha = RADIANS(atof(argv[2])) ;
    printf("using alpha=%2.0f degrees\n", DEGREES(alpha)) ;
  } else if (!strcmp(option, "EXAMPLE")) {
    example_T1 = argv[2] ;
    example_segmentation = argv[3] ;
    printf("using %s and %s as example T1 and segmentations respectively.\n",
           example_T1, example_segmentation) ;
    nargs = 2 ;
  }
  /////////////// writing out various samples /////////////////
  else if (!strcmp(option, "SAMPLES")) {
    sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing control points to %s...\n", sample_fname) ;
  } else if (!strcmp(option, "FSAMPLES") || !strcmp(option, "ISAMPLES")) {
    transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing transformed control points to %s...\n",
           transformed_sample_fname) ;
  } else if (!strcmp(option, "NSAMPLES")) {
    normalized_transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing  transformed normalization control points to %s...\n",
           normalized_transformed_sample_fname) ;
  }
  ///////////////////
  else if (!strcmp(option, "CONTRAST")) {
    use_contrast = 1 ;
    printf("using contrast to find labels...\n") ;
  } else if (!strcmp(option, "RENORM")) {
    renormalization_fname = argv[2] ;
    nargs = 1 ;
    printf("renormalizing using predicted intensity values in %s...\n",
           renormalization_fname) ;
  } else if (!strcmp(option, "FLASH_PARMS")) {
    tissue_parms_fname = argv[2] ;
    nargs = 1 ;
    printf("using FLASH forward model and tissue parms in %s to predict"
           " intensity values...\n", tissue_parms_fname) ;
  } else if (!strcmp(option, "TRANSONLY")) {
    translation_only = 1 ;
    printf("only computing translation parameters...\n") ;
  } else if (!strcmp(option, "WRITE_MEAN")) {
    gca_mean_fname = argv[2] ;
    nargs = 1 ;
    printf("writing gca means to %s...\n", gca_mean_fname) ;
  } else if (!strcmp(option, "PRIOR")) {
    min_prior = atof(argv[2]) ;
    nargs = 1 ;
    printf("using prior threshold %2.2f\n", min_prior) ;
  } else if (!strcmp(option, "SPACING")) {
    max_spacing = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using max GCA spacing %d...\n", max_spacing) ;
  } else if (!stricmp(option, "SCALES") || !stricmp(option, "SCALES")) {
    nscales = atoi(argv[2]) ;
    nargs = 1 ;
    printf("finding optimal linear transform over %d scales...\n", nscales);
  } else if (!stricmp(option, "NOVAR")) {
    novar = 1 ;
    printf("not using variance estimates\n") ;
  } else if (!strcmp(option, "DT")) {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("dt = %2.2e\n", parms.dt) ;
  } else if (!strcmp(option, "TOL")) {
    tol = parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("tol = %2.2e\n", parms.tol) ;
  } else if (!strcmp(option, "CENTER")) {
    center = 1 ;
    printf("using GCA centroid as origin of transform\n") ;
  } else if (!strcmp(option, "NOSCALE")) {
    noscale = 1 ;
    printf("disabling scaling...\n") ;
  } else if (!strcmp(option, "NOISCALE")) {
    noiscale = 1 ;
    printf("disabling intensity scaling...\n") ;
  } else if (!strcmp(option, "NUM")) {
    num_xforms = atoi(argv[2]) ;
    nargs = 1 ;
    printf("finding a total of %d linear transforms\n", num_xforms) ;
  } else if (!strcmp(option, "AREA")) {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_area = %2.2f\n", parms.l_area) ;
  } else if (!strcmp(option, "NLAREA")) {
    parms.l_nlarea = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_nlarea = %2.2f\n", parms.l_nlarea) ;
  } else if (!strcmp(option, "LEVELS")) {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("levels = %d\n", parms.levels) ;
  } else if (!strcmp(option, "INTENSITY") || !strcmp(option, "CORR")) {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_intensity = %2.2f\n", parms.l_intensity) ;
  } else if (!stricmp(option, "reduce")) {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    printf("reducing input images %d times before aligning...\n",
           nreductions) ;
  } else if (!stricmp(option, "nsamples")) {
    nsamples = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d samples of GCA...\n", nsamples) ;
  } else if (!stricmp(option, "norm")) {
    norm_fname = argv[2] ;
    nargs = 1 ;
    printf("intensity normalizing and writing to %s...\n",norm_fname);
  } else if (!stricmp(option, "trans")) {
    MAX_TRANS = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting max translation search range to be %2.1f\n", MAX_TRANS) ;
  } else if (!stricmp(option, "steps")) {
    max_angles = atoi(argv[2]) ;
    nargs = 1 ;
    printf("taking %d angular steps...\n", max_angles) ;
  } else switch (*option) {
    case 'L':   // for longitudinal analysis 
    {
      TRANSFORM *reg_transform ;

      xform_name = argv[2] ;
      long_reg_fname = argv[3] ;
      nargs = 2 ;
      printf("reading previously computed atlas xform %s "
             "and applying registration %s\n",
             xform_name, long_reg_fname) ;
      parms.transform = transform = TransformRead(argv[2]) ;
      if (transform == NULL)
        ErrorExit
        (ERROR_NOFILE,
         "%s: could not read transform from %s",
         Progname, argv[2]) ;
      Glta = parms.lta = (LTA *)transform->xform ;
      reg_transform = TransformRead(argv[3]) ;
      if (reg_transform == NULL)
        ErrorExit
        (ERROR_NOFILE,
         "%s: could not read registration from %s",
         Progname, argv[3]) ;
      transform_loaded = 1 ;
      TransformInvert(reg_transform, NULL) ;
      MatrixMultiply(((LTA *)(transform->xform))->xforms[0].m_L,
                     ((LTA *)(reg_transform->xform))->inv_xforms[0].m_L,
                     ((LTA *)(transform->xform))->xforms[0].m_L) ;
      TransformFree(&reg_transform) ;
    }
    break ;
    case 'F':
      ctl_point_fname = argv[2] ;
      nargs = 1 ;
      printf("reading manually defined control points from %s\n",
             ctl_point_fname) ;
      break ;
    case 'D':
      tx = atof(argv[2]) ;
      ty = atof(argv[3]) ;
      tz = atof(argv[4]) ;
      nargs = 3 ;
      break ;
    case 'R':
      rxrot = RADIANS(atof(argv[2])) ;
      ryrot = RADIANS(atof(argv[3])) ;
      rzrot = RADIANS(atof(argv[4])) ;
      nargs = 3 ;
      break ;
    case 'T':
      parms.transform = transform = TransformRead(argv[2]) ;
      Glta = parms.lta = (LTA *)transform->xform ;
#if 0
      parms.lta = LTAreadEx(argv[2]) ; // used to be LTAread()
#endif
      if (!parms.lta)
        ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",
                  Progname, argv[2]) ;
      if (parms.lta->type!=LINEAR_VOX_TO_VOX)
        ErrorExit(ERROR_BADFILE, "%s: must be LINEAR_VOX_TO_VOX (=0), but %d",
                  Progname, argv[2], parms.lta->type) ;
      nargs = 1 ;
      printf("using previously computed transform %s\n", argv[2]) ;
      if (parms.lta->type != LINEAR_VOX_TO_VOX) {
        fprintf(stdout,
                "ERROR: must use LINEAR_VOX_TO_VOX (=0) transform. "
                "The type was %d.\n",
                parms.lta->type);
        exit(1);
      }
      transform_loaded = 1 ;
      break ;
    case 'B':
      blur_sigma = atof(argv[2]) ;
      nargs = 1 ;
      printf("blurring input image with sigma=%2.3f\n", blur_sigma);
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'S':
#if 0
      parms.sigma = atof(argv[2]) ;
      printf("using sigma=%2.3f as upper bound on blurring.\n",
             parms.sigma) ;
      nargs = 1 ;
#else
      MAX_ANGLES = MAX_TRANS_STEPS = max_angles = (float)atoi(argv[2]) ;
      nargs = 1 ;
      printf("examining %2.0f different trans/rot/scale values...\n",
             MAX_ANGLES);
#endif
      break ;
    case '?':
    case 'U':
      printUsage();
      exit(1) ;
      break ;
    case 'N':
      parms.niterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("niterations = %d\n", parms.niterations) ;
      break ;
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'P':
      ctl_point_pct = atof(argv[2]) ;
      nargs = 1 ;
      printf("using top %2.1f%% wm points as control points....\n",
             100.0*ctl_point_pct) ;
      break ;
    case 'M':
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      printf("momentum = %2.2f\n", parms.momentum) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }
    */
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

  return (P.mov != "" && P.dst != "" && P.lta != "");
}
