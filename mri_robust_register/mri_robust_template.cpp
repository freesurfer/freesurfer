/**
 * @file  mri_robust_template.cpp
 * @brief combine multiple volumes by mean or median
 *
 * Creation of robust template of several volumes together with
 * their linear registration
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/08/05 16:30:00 $
 *    $Revision: 1.6.2.4 $
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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include "Registration.h"
#include "Regression.h"
#include "RobustGaussian.h"
#include "CostFunctions.h"

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

#define SAT 4.685 // this is suggested for gaussian noise
//#define SAT 20
#define SSAMPLE -1

// unsigned int printmemusage ()
// {
//   char buf[30];
//   snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
//   FILE* pf = fopen(buf, "r");
//   unsigned size = 0; //       total program size
//   if (pf) {
//     //unsigned resident;//   resident set size
//     //unsigned share;//      shared pages
//     //unsigned text;//       text (code)
//     //unsigned lib;//        library
//     //unsigned data;//       data/stack
//     //unsigned dt;//         dirty pages (unused in Linux 2.6)
//     //fscanf(pf, "%u" /* %u %u %u %u %u"*/, &size/*, &resident, &share, &text, &lib, &data*/);
//     fscanf(pf, "%u" , &size);
//     //DOMSGCAT(MSTATS, std::setprecision(4) << size / (1024.0) << "MB mem used");
//     cout <<  size / (1024.0) << " MB mem used" << endl;
//     fclose(pf);
//   }
// // while (1)
// // {
// //     if ('n' == getchar())
// //        break;
// // }
//   return size;
// }

struct Parameters
{
  vector < std::string > mov;
  string mean;
  MRI * mri_mean;
  vector <string> nltas;
  vector <string> nweights;
  bool   fixvoxel;
  bool   fixtype;
  bool   lta_vox2vox;
  bool   affine;
  bool   iscale;
  bool   transonly;
  bool   leastsquares;
  int    iterate;
  double epsit;
  double sat;
  vector <string> nwarps;
  int    debug;
  int    average;
  vector < MRI* > mri_mov;
  vector < LTA* > ltas;
  vector < MRI* > mri_warps;
  vector < MRI* > mri_weights;
  vector < double > intensities;
  int    inittp;
  bool   noit;
};

static struct Parameters P =
{
  vector < string >(0),"",NULL,vector < string >(0),vector < string >(0),false,false,false,false,false,false,false,5,0.01,SAT,vector < string >(0),0,1,vector < MRI* >(0),vector < LTA* >(0),vector < MRI* >(0),vector < MRI* >(0),vector < double >(0),1,false
};


static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[],Parameters & P) ;
static MRI* averageSet(const vector < MRI * >& set,
                       MRI* mean, int method, double sat);
static MRI* initialAverageSet(const vector < MRI * >& set,
                              MRI* mean, int method, double sat);
static int loadMovables(Parameters & P);
static void initRegistration(Registration & R, Parameters & P) ;

void computeTemplate (Parameters & P);
void initialXforms (Parameters & P, int tpi);
void halfWayTemplate (Parameters & P);

static char vcid[] =
"$Id: mri_robust_template.cpp,v 1.6.2.4 2009/08/05 16:30:00 nicks Exp $";
char *Progname = NULL;

//static MORPH_PARMS  parms ;
//static FILE *diag_fp = NULL ;


using namespace std;

int main(int argc, char *argv[])
{
  cout << vcid << endl;
  // set the environment variable
  // to store mri as chunk in memory:
//  setenv("FS_USE_MRI_CHUNK","",1) ;
  if (getenv("FS_USE_MRI_CHUNK") != NULL)
  {
    cerr <<
      "Error: do not set FS_USE_MRI_CHUNK while it is still buggy!" << endl;
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
//  if (P.outdir[P.outdir.length()-1] != '/') P.outdir += "/";


  // Timer
  struct timeb start ;
  int    msec,minutes,seconds;
  TimerStart(&start) ;
  ///////////////////////////////////////////////////////////////

  int nin = loadMovables(P);
  if (P.ltas.size() == 0) P.ltas.resize(nin,NULL);
  else assert (P.ltas.size() == P.mri_mov.size());
  P.mri_warps.resize(nin,NULL);
  P.intensities.resize(nin,1.0);
  if (!P.leastsquares) P.mri_weights.resize(nin,NULL);

  assert(nin >=2);

  if (nin == 2 && ! P.noit) halfWayTemplate(P);
  else
  {
    // if no initial xforms are given, use initialization:
    if(P.ltas[0] == NULL && P.inittp > 0 && ! P.noit) initialXforms(P,P.inittp);
		if (P.noit) assert(P.ltas[0] != NULL);
    computeTemplate(P);
  }

  cout << "Writing final template: " << P.mean << endl;
  strncpy(P.mri_mean->fname, P.mean.c_str(),STRLEN);
  MRIwrite(P.mri_mean,P.mean.c_str());
  MRIfree(&P.mri_mean);

  // output transforms and warps
  cout << "Writing final transforms (warps etc.)..." << endl;
  for (int i = 0;i<nin;i++)
  {
    if (P.nltas.size() > 0)
    {
      if (P.lta_vox2vox)
      {
        LTAchangeType(P.ltas[i], LINEAR_VOX_TO_VOX);
      }
      else assert(P.ltas[i]->type == LINEAR_RAS_TO_RAS);
      strncpy(P.ltas[i]->xforms[0].dst.fname, P.mean.c_str(),STRLEN);
      LTAwriteEx(P.ltas[i], P.nltas[i].c_str()) ;
			LTAfree(&P.ltas[i]);
    }
    if (P.nwarps.size() > 0)
		{
		  MRIwrite(P.mri_warps[i],P.nwarps[i].c_str()) ;
			MRIfree(&P.mri_warps[i]);
	  }
    if (P.iscale && P.intensities[i] >0 && P.nwarps.size() > 0)
    {
      size_t  pos = P.nwarps[i].rfind(".");    // position of "." in str
      if (pos!=string::npos)
      {
        string fn = P.nwarps[i].substr(0,pos) + "-intensity.txt";
        ofstream f(fn.c_str(),ios::out);
        f << P.intensities[i];
        f.close();
      }
      else
      {
        cerr << " output warp :  " <<
          P.nwarps[i] << " should end with .mgz!" << endl;
        exit (1);
      }

    }

    if (P.mri_weights[i] != NULL && P.nweights.size() > 0)
    {
      MRIwrite(P.mri_weights[i], P.nweights[i].c_str()) ;
      MRIfree(&P.mri_weights[i]);
    }
    MRIfree(&P.mri_mov[i]);
  }


  ///////////////////////////////////////////////////////////////
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  cout << "registration took "<<minutes<<" minutes and "<<
    seconds<<" seconds." << endl;
  //if (diag_fp) fclose(diag_fp) ;

  exit(0) ;
  return(0) ;
}

LTA* VOXmatrix2LTA(MATRIX * m, MRI* src, MRI* dst)
{
  LTA* ret =  LTAalloc(1,src);
  ret->xforms[0].m_L = MRIvoxelXformToRasXform (src,dst,m,NULL) ;
  ret->type = LINEAR_RAS_TO_RAS ;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

LTA* RASmatrix2LTA(MATRIX * m, MRI* src, MRI* dst)
{
  LTA* ret =  LTAalloc(1,src);
  ret->xforms[0].m_L = MatrixCopy(m,ret->xforms[0].m_L) ;
  ret->type = LINEAR_RAS_TO_RAS ;
  getVolGeom(src, &ret->xforms[0].src);
  getVolGeom(dst, &ret->xforms[0].dst);

  return ret;
}

void computeTemplate(Parameters & P)
// returns P.mri_mean, P.ltas, P.mri_warps (P.mri_weights and P.intensities if available)
{
  int nin = (int) P.mri_mov.size();

  string outdir = "./";
  size_t  pos = P.mean.rfind("/");    // position of "." in str
  if (pos!=string::npos)
  {
    outdir = P.mean.substr(0,pos-1);
  }
  
  cout << "Computing first template" << endl;  
  bool havexforms = (P.ltas[0] != NULL);
  if (!havexforms) 
    P.mri_mean = initialAverageSet(P.mri_mov,NULL,P.average,P.sat);
  else // we have initial transforms
  {  

    cout << "warping movs and creating template..." << endl;
    for (int i = 0;i<nin;i++)
    {  
      if (P.mri_warps[i]) MRIfree(&P.mri_warps[i]); // should not happen
      P.mri_warps[i] = MRIclone(P.mri_mov[0],P.mri_warps[i]);
      P.mri_warps[i] = LTAtransform(P.mri_mov[i],P.mri_warps[i], P.ltas[i]);
      if (P.debug)
      {
        ostringstream oss;
        oss << outdir << "tp" << i+1 << "_to_template-it0.mgz";
        MRIwrite(P.mri_warps[i], oss.str().c_str()) ;
      }
    }
    P.mri_mean = averageSet(P.mri_warps, P.mri_mean,P.average,P.sat);
  }
  
  if (P.noit) // no iterations necessary just return with mean (and ltas and warps);
  {
    strncpy(P.mri_mean->fname, P.mean.c_str(),STRLEN);
    return;
  }

  
  strncpy(P.mri_mean->fname,(outdir+"template-it0.mgz").c_str(),STRLEN);
  if (P.debug)
  {
    cout << "debug: saving template-it0.mgz" << endl;
    MRIwrite(P.mri_mean,(outdir+"template-it0.mgz").c_str());
  }

   cout << "template fname: " << P.mri_mean->fname << endl;

  int itmax  = 10;
  double eps = 0.025;

  int itcount = 0;
  double maxchange = 100;

  vector < MATRIX * > transforms(P.mri_mov.size(),NULL);
  if (havexforms) //init transforms
  {
     for (int i=0;i<nin;i++)
     {
        LTAchangeType(P.ltas[i],LINEAR_VOX_TO_VOX); //vox2vox for registration
        transforms[i] = MatrixCopy(P.ltas[i]->xforms[0].m_L, NULL);
     }
	
  }
  
  vector < Registration > Rv(P.mri_mov.size());
  for (int i = 0;i<nin;i++) Rv[i].setSource(P.mri_mov[i],P.fixvoxel,P.fixtype);

  LTA * lastlta = NULL;
  while (itcount < itmax && maxchange > eps)
  {
    itcount++;
    if (itcount > 1) maxchange = 0;
    cout << endl << "Working on global iteration : " << itcount << endl;
    cout << endl;
    //cout << "========================================" << endl;
    //printmemusage();
    //cout << "========================================" << endl << endl;
		
    // register all inputs to mean
    vector < double > dists(nin,1000); // should be larger than maxchange!
    for (int i = 0;i<nin;i++)
    {
      Rv[i].clear();
      initRegistration(Rv[i],P); //set parameter
      Rv[i].setTarget(P.mri_mean,
                      P.fixvoxel,
                      P.fixtype); // gaussian pyramid will be constructed for
                                  // each Rv[i], could be optimized
      ostringstream oss;
      oss << outdir << "tp" << i+1 << "_to_template-it" << itcount;
      Rv[i].setName(oss.str());

      // compute Alignment
      std::pair <MATRIX*, double> Md;
      int iterate = P.iterate;
      double epsit= P.epsit;
      int maxres = 0;
      int subsamp = -1;
      switch (itcount) // simplify first steps:
      {
      case 1:
        maxres = 2;
        break;
      case 2:
        maxres = 2;
        break;
      case 3:
        maxres = 1;
        break;
      case 4:
        maxres = 1;
        break;
      case 5:
        subsamp = 180;
        break;
      }

      if (havexforms)
      {
        //// tried high res only: (does not make sense, we need 2 iter. on high res anyway)
        //Md= Rv[i].computeIterativeRegistration(1,epsit,NULL,NULL,transforms[i],P.intensities[i]);
        //dists[i] = sqrt(Rv[i].AffineTransDistSq(Md.first, transforms[i]));
        if (itcount ==1) subsamp = 180;	
        maxres = 0; //go up to hig-res
      }
      //else dists[i] = 100; // we need to run multires...

      //if (dists[i] > eps)
      Rv[i].setSubsamplesize(subsamp);
      Md = Rv[i].computeMultiresRegistration(maxres,
                                               iterate,
                                               epsit,
                                               NULL,
                                               NULL,
                                               transforms[i],
                                               P.intensities[i]);
      if (transforms[i]) MatrixFree(&transforms[i]);
      transforms[i] = Md.first;
      P.intensities[i] = Md.second;

      // convert Matrix to LTA ras to ras
      if (lastlta)  LTAfree(&lastlta);
      if (P.ltas[i]) lastlta = P.ltas[i];
      P.ltas[i] = VOXmatrix2LTA(Md.first,P.mri_mov[i],P.mri_mean);
      //P.ltas[i] = LTAalloc(1,P.mri_mov[i]);
      //P.ltas[i]->xforms[0].m_L =
      //  MRIvoxelXformToRasXform (P.mri_mov[i],
      //                           P.mri_mean,
      //                           Md.first,
      //                           P.ltas[i]->xforms[0].m_L) ;
      //P.ltas[i]->type = LINEAR_RAS_TO_RAS ;
      //getVolGeom(P.mri_mov[i], &P.ltas[i]->xforms[0].src);
      //getVolGeom(P.mri_mean, &P.ltas[i]->xforms[0].dst);

      // compute maxchange
      if (lastlta)
      {
        LTAchangeType(lastlta,LINEAR_RAS_TO_RAS); //measure dist in RAS coords
        dists[i] = sqrt(Rv[i].AffineTransDistSq(lastlta->xforms[0].m_L,
                                                P.ltas[i]->xforms[0].m_L));
        LTAfree(&lastlta);
        if (dists[i] > maxchange) maxchange = dists[i];
        cout << endl << "tp " << i+1 << " distance: " << dists[i] << endl;
      }

      // here do scaling of intensity values
      if (Rv[i].isIscale() && Md.second > 0)
      {
        cout << "Adjusting Intensity of MOV by " << Md.second << endl;
        P.mri_mov[i] = Rv[i].MRIvalscale(P.mri_mov[i],
                                         P.mri_mov[i], Md.second);
      }

      // warp mov to mean
      cout << "warping mov to template..." << endl;
      if (P.mri_warps[i]) MRIfree(&P.mri_warps[i]);
      P.mri_warps[i] = MRIclone(P.mri_mean,P.mri_warps[i]);
      P.mri_warps[i] = LTAtransform(P.mri_mov[i],P.mri_warps[i], P.ltas[i]);
      
      //cout << " LS difference after: " <<
      //CF.leastSquares(mri_aligned,P.mri_dst) << endl;
      //cout << " NC difference after: " <<
      //CF.normalizedCorrelation(mri_aligned,P.mri_dst) << endl;


      if (P.debug)
      {
        cout << "debug: writing transforms, warps, weights ..." << endl;
        LTAwriteEx(P.ltas[i], (oss.str()+".lta").c_str()) ;

        MRIwrite(P.mri_warps[i], (oss.str()+".mgz").c_str()) ;

        if (Rv[i].isIscale() && Md.second >0)
        {
          string fn = oss.str() + "-intensity.txt";
          ofstream f(fn.c_str(),ios::out);
          f << Md.second;
          f.close();
        }

        //if (P.mri_weights[i]) MRIfree(&P.mri_weights[i]);
        P.mri_weights[i] = Rv[i].getWeights(); // in original half way space
        if (P.mri_weights[i] != NULL)
        {
          std::pair <MATRIX*, MATRIX*> map2weights = Rv[i].getHalfWayMaps();
          MATRIX * hinv = MatrixInverse(map2weights.second,NULL);

          MatrixPrint(stdout,map2weights.first) ;
          MatrixPrint(stdout,map2weights.second) ;
          MatrixPrint(stdout,hinv) ;

          MRI * wtarg = MRIalloc(P.mri_weights[i]->width,
                                 P.mri_weights[i]->height,
                                 P.mri_weights[i]->depth,
                                 MRI_FLOAT);
          MRIcopyHeader(P.mri_weights[i],wtarg);
          MATRIX * v2r  = MRIgetVoxelToRasXform(P.mri_mean);
          MRIsetVoxelToRasXform(wtarg,v2r);
          wtarg->type = MRI_FLOAT;
          wtarg->i_to_r__ = MatrixCopy(P.mri_mean->i_to_r__, wtarg->i_to_r__);
          wtarg->r_to_i__ = MatrixCopy(P.mri_mean->r_to_i__, wtarg->r_to_i__);

          wtarg = MRIlinearTransform(P.mri_weights[i],wtarg, hinv);
          MRIwrite(wtarg, (oss.str()+"-weights.mgz").c_str()) ;
          MRIwrite(P.mri_weights[i], (oss.str()+"-www.mgz").c_str());
          MRIfree(&wtarg);
          MatrixFree(&hinv);
          MatrixFree(&v2r);
        }
      }
    
      Rv[i].clear();
      Rv[i].freeGPT();
      cout << endl << "Finished TP : " << i+1 << endl;
      cout << endl;
      //cout << "========================================" << endl;
      //printmemusage();
      //cout << "========================================" << endl << endl;

    }

    if (dists[0] <= maxchange) // it was computed, so print it:
    {
      cout << "Iteration " << itcount <<" Distances: " << endl;
      for (unsigned int i=0;i<dists.size();i++) cout << dists[i] <<" ";
      cout << endl;
    }

    // create new mean
    cout << "Computing new template " <<itcount << endl;
    P.mri_mean = averageSet(P.mri_warps, P.mri_mean,P.average,P.sat);
    if (P.debug)
    {
      ostringstream oss;
      oss << outdir << "template-it" << itcount << ".mgz";
      cout << "debug: writing template to " << oss.str() << endl;
      MRIwrite(P.mri_mean,oss.str().c_str());
      strncpy(P.mri_mean->fname, oss.str().c_str(),STRLEN);
    }

  } // end while

  strncpy(P.mri_mean->fname, P.mean.c_str(),STRLEN);

  // copy weights if we have them (as the RV->weights will be destroyed):
  if (P.mri_weights[0])
    for (unsigned int i = 0;i<P.mri_weights.size();i++)
      P.mri_weights[i] = MRIcopy(P.mri_weights[i],NULL);

}

pair < MATRIX *, VECTOR * > getRTfromM(MATRIX * M, MATRIX * R, VECTOR * T)
{
   // check dimenstions:
   assert (M->rows == 4);
   assert (M->cols == 4);
   if (R == NULL) R = MatrixAlloc(3,3,M->type);
   if (T == NULL) T = VectorAlloc(3,M->type);
   assert(R->rows ==3);
   assert(R->cols ==3);
   assert(T->rows ==3);
   assert(T->cols ==1);
   
   // check M
   double eps = 0.000001;
   assert (fabs(M->rptr[4][1]) < eps);
   assert (fabs(M->rptr[4][2]) < eps);
   assert (fabs(M->rptr[4][3]) < eps);
   assert (fabs(M->rptr[4][4] - 1.0) < eps);
   
    for (int c=1; c<4;c++)
    {
      for (int r=1; r<4;r++)
      {
        R->rptr[r][c] = M->rptr[r][c];
      }
      T->rptr[c][1] = M->rptr[c][4];
   }
   
   return pair <MATRIX *, VECTOR *> (R,T);
}

MATRIX * getMfromRT(MATRIX * R, VECTOR * T, MATRIX * M)
{
   int type;
   if (R != NULL) type = R->type;
   else if ( T != NULL) type = T->type;
   else assert(R != NULL || T != NULL);
   
   if (M == NULL ) M = MatrixAlloc(4,4,type);
   if (R == NULL) R = MatrixIdentity(3,NULL);
   if (T == NULL) T = MatrixZero(3,1,NULL);
   
   // check dimensions:
   assert (M->rows == 4);
   assert (M->cols == 4);
   assert(R->rows ==3);
   assert(R->cols ==3);
   assert(T->rows ==3);
   assert(T->cols ==1);
   
   
    for (int c=1; c<4;c++)
    {
      for (int r=1; r<4;r++)
      {
        M->rptr[r][c] = R->rptr[r][c];
      }
      M->rptr[c][4] = T->rptr[c][1];
      M->rptr[4][c] = 0;
    }
    M->rptr[4][4] = 1;
   
   return M;
}

void initialXforms(Parameters & P, int tpi)
// will be returned in P.ltas
// tpi 1....n
{
  if (tpi <= 0) return; // no init
  cout << endl << "initializing Xforms (init " << tpi << " ) : " << endl;
  assert(P.ltas[0] == NULL);
  
  int nin = (int) P.mri_mov.size();
	
  tpi--; // tpi 0....n-1
  assert (tpi>=0 && tpi <nin);
	
  // create index set for lookup (tpi should switch places with 0)
  vector < int > index(nin);
  for (int i = 0;i<nin;i++) index[i] = i;
  index[0] = tpi;
  index[tpi] = 0;

  string outdir = "./";
  size_t  pos = P.mean.rfind("/");    // position of "." in str
  if (pos!=string::npos)
  {
    outdir = P.mean.substr(0,pos-1);
  }

  // Register everything to first TP (coarse)
  vector < Registration > Rv(nin);
  vector < std::pair <MATRIX*, double> > Md(nin);
  Md[0].first = MatrixIdentity(4,NULL);
  Md[0].second= 1.0;
  for (int i = 1;i<nin;i++) 
  {
	  int j = index[i];
    cout << endl << "  computing coarse registration of TP "<< j+1 <<" ( "<<P.mov[j]<<" ) wrt to TP "<<tpi+1<<" ( "<<P.mov[tpi]<<" )" << endl;

    ostringstream oss;
    oss << outdir << "tp" << j+1 << "_to_tp" << tpi;
    
    Registration R;
    R.setSource(P.mri_mov[j],P.fixvoxel,P.fixtype);
    initRegistration(R,P); //set parameter
    R.setTarget(P.mri_mov[tpi], P.fixvoxel, P.fixtype);
    R.setName(oss.str());
    // compute Alignment
    int iterate = 5;
    double epsit= 0.05;
    int maxres  = 1;
    Md[i] = R.computeMultiresRegistration(maxres,iterate,epsit);

//    ostringstream oss2;
//    oss2 << "tp" << j+1;    
//    MatrixWrite(Md[i].first,(oss2.str()+".fsmat").c_str(),oss2.str().c_str());
//    Md[i].first = MatrixRead((oss2.str()+".fsmat").c_str());
//    Md[i].second = 1.0;

//     if (P.debug)
//     {
//       LTA * nlta = VOXmatrix2LTA(Md[i].first,P.mri_mov[j],P.mri_mov[tpi]);
//       LTAwriteEx(nlta, (oss.str()+".lta").c_str()) ;
// 
//       MRI* warped = MRIclone(P.mri_mov[tpi],NULL);
//       warped = LTAtransform(P.mri_mov[j],warped, nlta);
// 
//       if (R.isIscale() && Md[i].second >0)
//       {
//         string fn = oss.str() + "-intensity.txt";
//         ofstream f(fn.c_str(),ios::out);
//         f << Md[i].second;
//         f.close();
// 	
//         warped = R.MRIvalscale(warped,warped, Md[i].second);
//       }
//       MRIwrite(warped, (oss.str()+".mgz").c_str()) ;
//       // todo: maybe output weights (not necessary so far) 
//       LTAfree(&nlta);
//       MRIfree(&warped);     
//     }
  }
  
  // find mean center
  vector <MATRIX* > mras(nin);
  mras[0] = MatrixIdentity(4,NULL);
  MATRIX* rot = MatrixIdentity(3,NULL);
  MATRIX* roti = NULL;
  MATRIX* trans = VectorAlloc(3,MATRIX_REAL);
  MATRIX* meant = MatrixZero(3,1,NULL);
  MATRIX* meanr = MatrixIdentity(3,NULL);
  for (int i = 1;i<nin;i++) 
  {
	  int j = index[i];
    cout << "  computing coord of TP "<< j+1 <<" ( "<<P.mov[j]<<" ) wrt to TP "<<tpi<<" ( "<<P.mov[tpi]<<" )" << endl;
    mras[i] = MRIvoxelXformToRasXform (P.mri_mov[j],P.mri_mov[tpi],Md[i].first,mras[i]);
    //MatrixPrintFmt(stdout,"% 2.8f",mras);cout << endl;
    // split into rotation translation
    getRTfromM(mras[i],rot,trans);
   //MatrixPrintFmt(stdout,"% 2.8f",trans); cout << endl;
   //cout << " Mras: " << endl;
   //MatrixPrintFmt(stdout,"% 2.8f",mras); cout << endl;

    // reverse order (first translate, then rotate)
    // rot stays, trans changes: Rx+t = R (x + R^{-1} t)
    roti  = MatrixTranspose(rot,roti); // inverse
    //cout << " roti: " << endl;MatrixPrintFmt(stdout,"% 2.8f",roti); cout << endl;
    
    trans = MatrixMultiply(roti,trans,trans);
    //cout << "transi: - " << endl; MatrixPrintFmt(stdout,"% 2.8f",trans[i]); cout << endl;
    
//     if (P.debug) // output transonly
//     {
//       ostringstream oss;
//       oss << outdir << "tp" << j+1 << "_to_tp"<<tpi<<"-transonly";
//       MATRIX *Mtemp = getMfromRT(NULL,trans,NULL);
//       LTA * nlta = RASmatrix2LTA(Mtemp,P.mri_mov[j],P.mri_mov[tpi]);
//       MRI* warped = LTAtransform(P.mri_mov[j],NULL, nlta);
//       MRIwrite(warped, (oss.str()+".mgz").c_str()) ;
//       LTAfree(&nlta);
//       MRIfree(&warped);
//       MatrixFree(&Mtemp);         
//     }
    meant = MatrixSubtract(meant,trans,meant);
    meanr = MatrixAdd(meanr,roti,meanr);
  }
  
  MatrixFree(&roti);
  MatrixFree(&trans);
  MatrixFree(&rot);

  //average
  meant = MatrixScalarMul(meant,1.0/nin,meant);
  //cout << "meant: " << endl; MatrixPrintFmt(stdout,"% 2.8f",meant); cout << endl;
  meanr = MatrixScalarMul(meanr,1.0/nin,meanr);
  //cout << "meanr: " << endl;MatrixPrintFmt(stdout,"% 2.8f",meanr); cout << endl;
  
  // project meanr back to SO(3) (using polar decomposition)
  VECTOR * vz = VectorAlloc(3,MATRIX_REAL);
  MATRIX * mv = MatrixSVD(meanr,vz,NULL);
  //MatrixPrintFmt(stdout,"% 2.8f",meanr); cout << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",vz); cout << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",mv); cout << endl;
  mv = MatrixTranspose(mv,mv);
  //MatrixPrintFmt(stdout,"% 2.8f",mv); cout << endl;
  meanr = MatrixMultiply(meanr,mv,meanr);
  //cout << " meanr.proj: " << endl;MatrixPrintFmt(stdout,"% 2.8f",meanr); cout << endl;
  MatrixFree(&mv);
  MatrixFree(&vz);

  // construct maps from each image to the mean space
  MATRIX * Mm = getMfromRT(meanr,meant,NULL);
  MATRIX * M = NULL;
  //cout << " Mm: " << endl; MatrixPrintFmt(stdout,"% 2.8f",Mm); cout << endl;
  for (int i = 0;i<nin;i++) 
  {
	  int j = index[i];
    cout << "  computing mean coord of TP "<< j+1 <<" ( "<<P.mov[j]<<" ) " << endl;
    // concat transforms: meant meanr mras
    // from right to left: mras[j] (aligns to T1)
    //                     then do the mean rot and move to mean location

    //MatrixPrintFmt(stdout,"% 2.8f",mras[i]); cout << endl;
    M = MatrixMultiply(Mm,mras[i],M);
    //MatrixPrintFmt(stdout,"% 2.8f",M); cout << endl;
    
    // make lta from M (M is RAS to RAS)
    assert(P.ltas[j] == NULL);
    P.ltas[j] = RASmatrix2LTA(M,P.mri_mov[j],P.mri_mov[tpi]); // use geometry of tpi
    
    if (P.debug)
    {
      ostringstream oss;
      oss << outdir << "tp" << j+1 << "_to_mcoord";
      LTAwriteEx(P.ltas[j], (oss.str()+".lta").c_str()) ;
      MRI * warped = LTAtransform(P.mri_mov[j],NULL, P.ltas[j]);
//      if (R.isIscale() && Md[i].second >0)
//      {
//        string fn = oss.str() + "-intensity.txt";
//        ofstream f(fn.c_str(),ios::out);
//        f << Md[i].second;
//        f.close();
//	
//        warped = R.MRIvalscale(warped,warped, Md[i].second);
//      }
      MRIwrite(warped, (oss.str()+".mgz").c_str()) ;
      MRIfree(&warped);     
    }
    
    // cleanup
    MatrixFree(&mras[i]);
    MatrixFree(&Md[i].first);
  }
  //cleanup
  MatrixFree(&M);
  MatrixFree(&Mm);
  
}

void halfWayTemplate(Parameters & P)
// can be used only with two input images
{
  int nin = (int) P.mri_mov.size();
  assert (nin == 2);

  // register 1 with 2

  Registration R;
  initRegistration(R,P); //set parameter
  R.setSource(P.mri_mov[0],P.fixvoxel,P.fixtype);
  R.setTarget(P.mri_mov[1],P.fixvoxel,P.fixtype);

//   ostringstream oss;
//   oss << outdir << "tp" << i+1 << "_to_template-it" << itcount;
  R.setName(P.mean);

  // use initial transform if given:
  MATRIX *minit = NULL;
  if (P.ltas[0])
  {
    assert (P.ltas[0]->type == LINEAR_VOX_TO_VOX);
    assert (P.ltas[1]->type == LINEAR_VOX_TO_VOX);
    minit = MatrixInverse(P.ltas[1]->xforms[0].m_L,NULL);
    minit = MatrixMultiply(minit,P.ltas[0]->xforms[0].m_L,minit);
    R.setDebug(1);
  }		
  // compute Alignment
  std::pair <MATRIX*, double> Md;
  Md = R.computeMultiresRegistration(0,P.iterate,P.epsit,NULL,NULL,minit);
	if (minit) MatrixFree(&minit);
  P.intensities[0] = Md.second;
  P.intensities[1] = 1;

  cout << "creating half way data ..." << endl;
  std::pair < MATRIX*, MATRIX*> maps2weights = R.getHalfWayMaps();
  LTA * m2hwlta = LTAalloc(1,P.mri_mov[0]);
  LTA * d2hwlta = LTAalloc(1,P.mri_mov[1]);
  if (!P.lta_vox2vox) // do ras to ras
  {
    // cout << "converting VOX to RAS and saving RAS2RAS..." << endl ;
    // (use geometry of destination space for half-way)
    m2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform (P.mri_mov[0],
                                                      P.mri_mov[1],
                                                      maps2weights.first,
                                                      m2hwlta->xforms[0].m_L) ;
    m2hwlta->type = LINEAR_RAS_TO_RAS ;
    d2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform (P.mri_mov[1],
                                                      P.mri_mov[1],
                                                      maps2weights.second,
                                                      d2hwlta->xforms[0].m_L) ;
    d2hwlta->type = LINEAR_RAS_TO_RAS ;
  }
  else // vox to vox
  {
    // cout << "saving VOX2VOX..." << endl ;
    m2hwlta->xforms[0].m_L = MatrixCopy(maps2weights.first,
                                        m2hwlta->xforms[0].m_L) ;
    m2hwlta->type = LINEAR_VOX_TO_VOX ;
    d2hwlta->xforms[0].m_L = MatrixCopy(maps2weights.second,
                                        m2hwlta->xforms[0].m_L) ;
    d2hwlta->type = LINEAR_VOX_TO_VOX ;
  }
  // add src and dst info (use src as target geometry in both cases)
  getVolGeom(P.mri_mov[0], &m2hwlta->xforms[0].src);
  getVolGeom(P.mri_mov[0], &m2hwlta->xforms[0].dst);
  getVolGeom(P.mri_mov[1], &d2hwlta->xforms[0].src);
  getVolGeom(P.mri_mov[0], &d2hwlta->xforms[0].dst);

  P.ltas.resize(2,NULL);
  if (P.ltas[0]) LTAfree(&P.ltas[0]);
  if (P.ltas[1]) LTAfree(&P.ltas[1]);
  P.ltas[0] = m2hwlta;
  P.ltas[1] = d2hwlta;

  cout << " creating half-way template ..." << endl;
  // take dst info from lta:
  if (P.mri_warps[0]) MRIfree(&P.mri_warps[0]);
  if (P.mri_warps[1]) MRIfree(&P.mri_warps[1]);
  P.mri_warps[0] = LTAtransform(P.mri_mov[0],NULL, m2hwlta);
  P.mri_warps[1] = LTAtransform(P.mri_mov[1],NULL, d2hwlta);
  // here do scaling of intensity values
  if (R.isIscale() && Md.second > 0)
  {
    cout << "Adjusting Intensity of MOV by " << Md.second << endl;
    P.mri_warps[0] = R.MRIvalscale(P.mri_warps[0], P.mri_warps[0], Md.second);
  }
  P.mri_mean = averageSet(P.mri_warps, P.mri_mean,P.average,P.sat);


  P.mri_weights[0] = R.getWeights(); // in original half way space
  P.mri_weights[1] = P.mri_weights[0];


  // copy weights if we have them (as the RV->weights will be destroyed):
  if (P.mri_weights[0])
    for (unsigned int i = 0;i<P.mri_weights.size();i++)
      P.mri_weights[i] = MRIcopy(P.mri_weights[i],NULL);
			
  MatrixFree(&Md.first);
}


/*----------------------------------------------------------------------
  ----------------------------------------------------------------------*/
static void printUsage(void)
{
  cout << endl << endl;
  cout << "Usage: mri_robust_template <required arguments>" << endl <<endl;

  cout << "Short description: constructs template and rigidly registers volumes (using robust statistics)" << endl;
  cout << "                   will output the template volume and the transformations (lta)"<< endl << endl;

  cout << "Required arguments" << endl << endl ;
  cout << "  --mov tp1.mgz tp2.mgz ...  movable/source volumes to be registered" << endl;
  cout << "  --template template.mgz    output template volume" << endl;
  cout << endl;
  cout << "Optional arguments" << endl << endl;
//  cout << "      --outdir               output directory (default ./ )" << endl;
  cout << "  --lta tp1.lta tp2.lta ...  output xforms to template (for each input)" << endl;
  cout << "  --warp warp1.mgz ...       map each input to template" << endl;
  cout << "  --weights weights1.mgz ... output weights in target space" << endl;
  cout << "  --average #                construct template from: 0 Mean, 1 Median (default)" << endl;
  cout << "  --inittp #                 use TP# for spacial init (default 1), 0: no init" << endl;

//  cout << "  -A, --affine (testmode)    find 12 parameter affine xform (default is 6-rigid)" << endl;
  cout << "  -I, --iscale               allow also intensity scaling on high-res. (default no)" << endl;
//  cout << "      --transonly            find 3 parameter translation only" << endl;
//  cout << "  -T, --transform lta        use initial transform lta on source ('id'=identity)" << endl;
  cout << "  --ixforms t1.lta t2.lta .. use initial transforms (lta) on source  ('id'=identity)" << endl;
//  cout << "                                default is geometry (RAS2VOX_dst * VOX2RAS_mov)" << endl;
  cout << "  --vox2vox                  output VOX2VOX lta file (default is RAS2RAS)" << endl;
  cout << "  --leastsquares             use least squares instead of robust M-estimator" << endl;
  cout << "  --noit                     do not iterate, just create first template" << endl;
//  cout << "      --maxit <#>            iterate max # times on each resolution (default "<<P.iterate<<")"  << endl;
//  cout << "      --epsit <float>        stop iterations when below <float> (default "<<P.epsit <<")" << endl;
//  cout << "      --nomulti              work on highest resolution (no multiscale)" << endl;
  cout << "  --sat <float>              set saturation for robust estimator (default "<<SAT<<")" << endl;
//  cout << "      --subsample <#>        subsample if dim > # on all axes (default no subs.)" << endl;
//  cout << "      --maskmov mask.mgz     mask mov/src with mask.mgz" << endl;
//  cout << "      --maskdst mask.mgz     mask dst/target with mask.mgz" << endl;
  cout << "  --uchar                    set input type to UCHAR (with intensity scaling)" << endl;
  cout << "  --conform                  conform volumes to 1mm vox (256^3)" << endl;
//  cout << "      --satit                iterate on highest res with different sat" << endl;
  cout << "  --debug                    show debug output (default no debug output)" << endl;
//  cout << "      --test i mri         perform test number i on mri volume" << endl;

  cout << endl;
  cout << "Mandatory or optional arguments to long options are also mandatory or optional for any" << endl;
  cout << " corresponding short options." << endl;
  cout << endl;

  cout << "Report bugs to: analysis-bugs@nmr.mgh.harvard.edu" << endl;

  cout << endl;
}

/*!
  \fn MRI* averageSet(const vector < MRI * >& set, MRI* mean, int method, double sat)
  \brief Averages the movable volumes depending on method
  \param set vector of movable volumes
  \param mean  output of mean volume
  \param method  0 = mean, 1 = median, 2 = tukey biweight (testing)
  \param sat     saturation for tukey biweight
*/
static MRI* averageSet(const vector < MRI * >& set,
                       MRI* mean, int method, double sat)
{
  assert(set.size() > 1);

//    for(uint i = 0;i<set.size();i++)
//    {
//       cout << " TP " << i+1 << endl;
//       cout << " mean   : " << CostFunctions::mean(set[i])   << endl;
//       cout << " sdev   : " << CostFunctions::sdev(set[i])   << endl;
//       cout << " median : " << CostFunctions::median(set[i]) << endl;
//       cout << " mad    : " << CostFunctions::mad(set[i])    << endl;
//     }

  if (method == 0)
  {
    // mean
    cout << "    using mean" << endl;
    for (unsigned int i = 0;i<set.size();i++)
    {
      mean = MRIaverage(set[i],i,mean);
    }
  }
  else if (method ==1)
  {
    cout << " using median " << endl;
    // robust
    int x,y,z,i;
    assert(set.size() > 0);
    double dd[set.size()];
    if (!mean) mean = MRIclone(set[0],NULL);
    //MRI * midx = MRIalloc(set[0]->width,set[0]->height,set[0]->depth, MRI_FLOAT);
    //MRIcopyHeader(mean,midx);
    //midx->type = MRI_FLOAT;
    pair < double, double > mm;
    for (z = 0 ; z < set[0]->depth ; z++)
      for (y = 0 ; y < set[0]->height ; y++)
        for (x = 0 ; x < set[0]->width ; x++)
        {
          for (i=0; i<(int) set.size();i++)
            dd[i] = MRIgetVoxVal(set[i],x,y,z,0);
          mm = RobustGaussian::medianI(dd,(int)set.size());
          MRIsetVoxVal(mean,x,y,z,0,mm.first);
	  //MRIsetVoxVal(midx,x,y,z,0,mm.second);
        }
    //MRIwrite(midx,"midx.mgz");
    //MRIwrite(mean,"m.mgz");
    //assert(1==2);
  }
  else if (method ==2)
  {
    cout << "    using tukey biweight" << endl;
    // robust tukey biweight
    int x,y,z,i;
    MATRIX* b = MatrixAlloc(set.size(),1,MATRIX_REAL);
    pair < MATRIX* , MATRIX* >  muw;
    assert(set.size() > 0);
    if (!mean) mean = MRIclone(set[0],NULL);
    for (z = 0 ; z < set[0]->depth ; z++)
      for (y = 0 ; y < set[0]->height ; y++)
        for (x = 0 ; x < set[0]->width ; x++)
        {
          // cout << " x: " << x << " y: " << y << " z: " <<z << "  size: " << set.size() <<endl;
          for (i=0; i<(int) set.size();i++)
          {
            //cout << "i: " << i << endl;
            b->rptr[i+1][1] = MRIgetVoxVal(set[i],x,y,z,0);
          }
          //cout << endl << "intgensities at this voxel:" << endl; ;
          //MatrixPrintFmt(stdout,"% 2.8f",b);

          Regression R(b);
          muw = R.getRobustEstW(sat);
          //    cout << " tukey mean: " << muw.first->rptr[1][1] << endl;
          MRIsetVoxVal(mean,x,y,z,0,muw.first->rptr[1][1]);
          MatrixFree(&muw.first);
          MatrixFree(&muw.second);
        }
  }
  else if (method ==3) // needs more development (sigma..)
  {
    cout << "    using tukey biweight (alltoghether)" << endl;
    // robust tukey biweight
    int x,y,z,i;
    int n = set[0]->depth * set[0]->height *  set[0]->width ;
    MATRIX* b = MatrixAlloc(set.size()*n,1,MATRIX_REAL);
    MATRIX* A = MatrixAlloc(set.size()*n,n,MATRIX_REAL);
    A = MatrixZero(A->rows,A->cols,A);
    pair < MATRIX* , MATRIX* >  muw;
    assert(set.size() > 0);
    if (!mean) mean = MRIclone(set[0],NULL);
    for (i=0; i<(int) set.size();i++)
    {
      assert(set[i]->width == set[0]->width);
      assert(set[i]->height == set[0]->height);
      assert(set[i]->depth == set[0]->depth);
      int pcount = 0;
      for (z = 0 ; z < set[0]->depth ; z++)
        for (y = 0 ; y < set[0]->height ; y++)
          for (x = 0 ; x < set[0]->width ; x++)
          {
            // cout << " x: " << x << " y: " << y << " z: " <<z << "  size: " << set.size() <<endl;
            //cout << "i: " << i << endl;
            b->rptr[pcount*set.size()+i+1][1] =
              MRIgetVoxVal(set[i],x,y,z,0);
            A->rptr[pcount*set.size()+i+1][pcount+1] = 1;
            pcount++;
          }
    }
    //cout << endl << "intgensities at this voxel:" << endl; ;
    //MatrixPrintFmt(stdout,"% 2.8f",b);

    Regression R(A,b);
    muw = R.getRobustEstW(sat);
    //    cout << " tukey mean: " << muw.first->rptr[1][1] << endl;
    int pcount = 0;
    for (z = 0 ; z < set[0]->depth ; z++)
      for (y = 0 ; y < set[0]->height ; y++)
        for (x = 0 ; x < set[0]->width ; x++)
        {
          pcount++;
          MRIsetVoxVal(mean,x,y,z,0,muw.first->rptr[pcount][1]);
        }
    MatrixFree(&muw.first);
    MatrixFree(&muw.second);

  }
  else
  {

    cerr <<  " averageSet  method " << method << " unknown" << endl;
    exit(1);
  }
//       cout << " Average Vol "  << endl;
//       cout << " mean   : " << CostFunctions::mean(mean) << endl;
//       cout << " sdev   : " << CostFunctions::sdev(mean) << endl;
//       cout << " median : " << CostFunctions::median(mean) << endl;
//       cout << " mad    : " << CostFunctions::mad(mean)<< endl;

  return mean;
}

/*!
  \fn MRI* initialAverageSet(const vector < MRI * >& set, MRI* mean, int method, double sat)
  \brief Averages the movable volumes depending on method
  \param set vector of movable volumes
  \param mean  output of mean volume
  \param method  0 = mean, 1 = median, 2 = tukey biweight (testing)
  \param sat     saturation for tukey biweight
*/
static MRI* initialAverageSet(const vector < MRI * >& set,
                              MRI* mean, int method, double sat)
{
// the initial average can be anyplace as long as it is
// not identical to one of the
// initial images (to avoid introducing a bias).
// it can be the average of all images but that is extremely blurry, and if
// we have only few images, far apart, maybe the algorithm will not converge
// or converge very slowly (has not happened yet).
// therefore we find an initial coarse alignment by using moments

  assert(set.size() > 1);
  int n = (int) set.size();
  vector < double >  centroid(3,0);
  vector < vector < double > >centroids(n);
  for (int i = 0;i<n;i++)
  {
    centroids[i] = CostFunctions::centroid(set[i]);
    centroid[0] += centroids[i][0];
    centroid[1] += centroids[i][1];
    centroid[2] += centroids[i][2];
  }
  centroid[0] /= n;
  centroid[1] /= n;
  centroid[2] /= n;

  // map set to mean of centroids
  MATRIX* Mtrans = MatrixIdentity(4,NULL);
  vector < MRI* > newset(n);
  for (int i = 0;i<n;i++)
  {
    Mtrans->rptr[1][4] = centroid[0]-centroids[i][0];
    Mtrans->rptr[2][4] = centroid[1]-centroids[i][1];
    Mtrans->rptr[3][4] = centroid[2]-centroids[i][2];
    newset[i] = MRIlinearTransform(set[i],NULL,Mtrans);
  }

  return averageSet(newset,mean,method,sat);
}

/*!
  \fn void loadMovables(Parameters & P)
  \brief Loads the movable volumes as specified on command line
  \param P  Paramters for the initialization
*/
static int loadMovables(Parameters & P)
{

  assert (P.mri_mov.size () == 0);
  int n = (int) P.mov.size();
  P.mri_mov.resize(n);

  for (unsigned int i = 0;i<P.mov.size(); i++)
  {
    cout << "reading source '"<<P.mov[i]<<"'..." << endl;

    P.mri_mov[i] = MRIread(P.mov[i].c_str()) ;
    if (!P.mri_mov[i])
    {
      ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
                Progname, P.mov[i].c_str()) ;
    }
  }

  return P.mov.size();
}

/*!
  \fn void initRegistration(Registration & R, const Parameters & P)
  \brief Initializes a Registration with Parameters (affine, iscale, transonly, leastsquares, sat and trans)
  \param R  Registration to be initialized
  \param P  Paramters for the initialization
*/
static void initRegistration(Registration & R, Parameters & P)
{
  // assert(n < (int) P.mov.size());

  R.setRigid(!P.affine);
  R.setIscale(P.iscale);
  R.setTransonly(P.transonly);
  R.setRobust(!P.leastsquares);
  R.setSaturation(P.sat);
  //R.setDebug(P.debug);

//   int pos = P.mov[n].rfind(".");
//   if (pos > 0) R.setName(P.mov[n].substr(0,pos));
//   else  R.setName(P.mov[n]);
//
//  // if (P.subsamplesize > 0) R.setSubsamplesize(P.subsamplesize);
//
//
//   R.setSource(P.mri_mov[n],P.fixvoxel,P.fixtype);
//   R.setTarget(P.mri_mean,P.fixvoxel,P.fixtype);
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

  if (!strcmp(option, "MOV") )
  {
    nargs = 0;
    do
    {
      option = argv[nargs+1];
      if (option[0] != '-')
      {
        nargs++;
        P.mov.push_back(string(argv[nargs]));
        cout << "Using "<< P.mov.back() <<
          " as movable/source volume." << endl;
      }
    }
    while (nargs+1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "    Using " << nargs << " input volumes" << endl;
  }
//  else if (!strcmp(option, "OUTDIR") )
//  {
//     P.outdir = string(argv[1]);
//     nargs = 1;
//     cout << "Using "<< P.outdir << " as output directory." << endl;
//  }
  else if (!strcmp(option, "TEMPLATE") )
  {
    P.mean = string(argv[1]);
    nargs = 1;
    cout << "Using "<< P.mean << " as template output volume." << endl;
  }
  else if (!strcmp(option, "LTA")   )
  {
    nargs = 0;
    do
    {
      option = argv[nargs+1];
      if (option[0] != '-')
      {
        nargs++;
        P.nltas.push_back(string(argv[nargs]));
        //cout << "Using "<< P.nltas.back() << " as LTA." << endl;
      }
    }
    while (nargs+1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "Will output LTA transforms" << endl;
  }
  else if (!strcmp(option, "IXFORMS")   )
  {
    nargs = 0;
    do
    {
      option = argv[nargs+1];
      if (option[0] != '-')
      {
        nargs++;
	
	
        TRANSFORM * trans = TransformRead(argv[nargs]);
        LTA* lta =  (LTA *)trans->xform ;
        if (!lta)
          ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",Progname, argv[nargs]) ;
        lta = LTAchangeType(lta,LINEAR_VOX_TO_VOX);

        P.ltas.push_back(lta);
        //cout << "Using "<< P.nltas.back() << " as LTA." << endl;
      }
    }
    while (nargs+1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "Will use init XFORMS." << endl;
  }
  else if (!strcmp(option, "AVERAGE") )
  {
    P.average = atoi(argv[1]);
    nargs = 1;
    cout << "Using Method: " << P.average <<
      " for template computation." <<  endl;
  }
  else if (!strcmp(option, "VOX2VOX")   )
  {
    P.lta_vox2vox = true;
    cout << "Output transforms as VOX2VOX. " << endl;
  }
  else if (!strcmp(option, "AFFINE") || !strcmp(option, "A") )
  {
    P.affine = true;
    cout << "Enableing affine transform!" << endl;
  }
  else if (!strcmp(option, "ISCALE") || !strcmp(option, "I") )
  {
    P.iscale = true;
    cout << "Enableing intensity scaling!" << endl;
  }
  else if (!strcmp(option, "TRANSONLY"))
  {
    P.transonly = true;
    cout << "Using only translation!" << endl;
  }
  else if (!strcmp(option, "LEASTSQUARES") || !strcmp(option, "L")  )
  {
    P.leastsquares = true;
    cout << "Using standard least squares (non-robust)!" << endl;
  }
  else if (!strcmp(option, "MAXIT")  )
  {
    P.iterate = atoi(argv[1]);
    nargs = 1 ;
    cout << "Performing maximal " << P.iterate <<
      " iterations on each resolution" << endl;
  }
  else if (!strcmp(option, "EPSIT") )
  {
    P.epsit = atof(argv[1]);
    nargs = 1 ;
    cout << "Stop iterations when change is less than " << P.epsit <<
      " . " << endl;
  }
  else if (!strcmp(option, "SAT")  )
  {
    P.sat = atof(argv[1]);
    nargs = 1 ;
    cout << "Using saturation " << P.sat << " in M-estimator!" << endl;
  }
  else if (!strcmp(option, "DEBUG") )
  {
    P.debug = 1;
    nargs = 0 ;
    cout << "Will output debug info and files!" << endl;
  }
  else if (!strcmp(option, "NOIT") )
  {
    P.noit = true;
    nargs = 0 ;
    cout << "Will output only first template (no iterations)!" << endl;
  }
  else if (!strcmp(option, "WEIGHTS") )
  {
    nargs = 0;
    do
    {
      option = argv[nargs+1];
      if (option[0] != '-')
      {
        nargs++;
        P.nweights.push_back(string(argv[nargs]));
        //cout << "Using "<< P.nweights.back() << " as weights volume." << endl;
      }
    }
    while (nargs+1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "Will output weights in target space" << endl;
  }
  else if (!strcmp(option, "WARP") || !strcmp(option, "W") )
  {
    nargs = 0;
    do
    {
      option = argv[nargs+1];
      if (option[0] != '-')
      {
        nargs++;
        P.nwarps.push_back(string(argv[nargs]));
        //cout << "Using "<< P.nwarps.back() << " as weights volume." << endl;
      }
    }
    while (nargs+1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "Will save warped sources !" << endl;
  }
  else if (!strcmp(option, "TEST"))
  {
    cout << " TEST-MODE " << endl;
    Registration R;
    R.testRobust(argv[2], atoi(argv[1]));
    nargs = 2 ;
    exit(0);
  }
  else if (!strcmp(option, "CONFORM") )
  {
    P.fixvoxel = true;
    nargs = 0 ;
    cout << "Will conform images to 256^3 and voxels to 1mm!" << endl;
  }
  else if (!strcmp(option, "UCHAR") )
  {
    P.fixtype = true;
    nargs = 0 ;
    cout << "Changing type to UCHAR (with intesity scaling)!" << endl;
  }
  else if (!strcmp(option, "INITTP") )
  {
    P.inittp = atoi(argv[1]);
    nargs = 1 ;
    if (P.inittp ==0 ) cout<< " No initialization, construct first mean from original TPs" << endl;
    else cout << "Using TP " <<P.inittp <<" as target for initialization" << endl;
  }
  else
  {
    cerr << "Option: " << argv[0] << " unknown !! " << endl;
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

  bool ntest = true;
  if (P.nwarps.size() > 0 && P.mov.size() != P.nwarps.size())
  {
    ntest = false;
    cerr << " No. of filnames for --warp should agree with no. of inputs!"
         << endl;
    exit(1);
  }
  if (P.nltas.size() > 0 && P.mov.size() != P.nltas.size())
  {
    ntest = false;
    cerr << " No. of filnames for --lta should agree with no. of inputs!"
         << endl;
    exit(1);
  }
  if (P.nweights.size() > 0 && P.mov.size() != P.nweights.size())
  {
    ntest = false;
    cerr << " No. of filnames for --weights should agree with no. of inputs!"
         << endl;
    exit(1);
  }

  return (P.mov.size() >= 2 && P.mean != "" && ntest);
}
