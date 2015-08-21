/**
 * @file  mri_coreg.c
 * @brief Computes registration between two volumes
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2015/08/21 21:05:30 $
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


/*!
\file dummy.c
\brief Example c file that can be used as a template.
\author Douglas Greve
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "mrisurf.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "timer.h"
#include "mrimorph.h"
#include "fmriutils.h"
#include "fsenv.h"
#include "matfile.h"
#include "icosahedron.h"
#include "cpputils.h"
#include "numerics.h"
#include "randomfields.h"

double round(double x);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_coreg.c,v 1.13 2015/08/21 21:05:30 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

typedef struct {
  char *mov;
  char *ref;
  char *refmask;
  char *movmask;
  char *outreg;
  char *regdat;
  char *subject;
  int DoCoordDither;
  int DoIntensityDither;
  int dof;
  double params[12];
  int nsep, seplist[10];
  int DoInitCostOnly;
  int DoSmoothing;
  int cras0;
  double ftol,linmintol;
  int nitersmax;
} CMDARGS;

CMDARGS *cmdargs;

MRI *MRIrescaleToUChar(MRI *mri, MRI *ucmri, double sat);
unsigned char *MRItoUCharVect(MRI *mri, RFS *rfs);
MATRIX *MRIgetVoxelToVoxelXformBase(MRI *mri_src, MRI *mri_dst, MATRIX *SrcRAS2DstRAS, MATRIX *SrcVox2DstVox, int base);

double **conv1dmat(double **M, int rows, int cols, double *v, int nv, int dim, 
		   double **C, int *pcrows, int *pcols);
int conv1dmatTest(void);
double *conv1d(double *v1, int n1, double *v2, int n2, double *vc);

double **AllocDoubleMatrix(int rows, int cols);
int FreeDoubleMatrix(double **M, int rows, int cols);
int WriteDoubleMatrix(char *fname, char *fmt, double **M, int rows, int cols);
int PrintDoubleMatrix(FILE *fp, char *fmt, double **M, int rows, int cols);
double *SumVectorDoubleMatrix(double **M, int rows, int cols, int dim, double *sumvect, int *nv);

typedef struct {
  MRI *ref, *mov, *refmask, *movmask;
  int seplist[10],nsep,sep,sepmin;
  double SatPct,refsat, movsat;
  unsigned char *g,*f;
  MATRIX *M,*V2V;
  double reffwhm[3],refgstd[3];
  double movfwhm[3],movgstd[3];
  double histfwhm[2];
  double params[12];
  int nparams;
  double H01d[256*256];
  double **H0;
  double cost;
  int nCostEvaluations;
  double tLastEval;
  double ftol,linmintol;
  float fret;
  int nitersmax,niters;
  int startmin;
  int nhits,nvoxmov;
  double pcthits;
  int DoCoordDither;
  RFS *crfs;
  MRI *cdither;
  int DoIntensityDither;
  RFS *refirfs,*movirfs;
  int DoSmoothing;
} COREG;

double COREGcost(COREG *coreg);
float COREGcostPowell(float *pPowel) ;
int COREGMinPowell();
float MRIgetPercentile(MRI *mri, double Pct, int frame);
int COREGfwhm(MRI *mri, double sep, double fwhm[3]);
int COREGpreproc(COREG *coreg);
LTA *LTAcreate(MRI *src, MRI *dst, MATRIX *T, int type);
int COREGhist(COREG *coreg);
long COREGvolIndex(int ncols, int nrows, int nslices, int c, int r, int s);
double COREGsamp(unsigned char *f, const double c, const double r, const double s, 
		  const int ncols, const int nrows, const int nslices);
double NMICost(double **H, int cols, int rows);
MATRIX *COREGmatrix(double *p0, int np, MATRIX *M);
double *COREGparams9(MATRIX *M9, double *p);
int COREGprint(FILE *fp, COREG *coreg);

COREG *coreg;
FSENV *fsenv;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err,n;
  LTA *lta;
  struct timeb timer;

  TimerStart(&timer);

  cmdargs = (CMDARGS *)calloc(sizeof(CMDARGS),1);
  cmdargs->mov = NULL;
  cmdargs->ref = NULL;
  cmdargs->outreg = NULL;
  cmdargs->regdat = NULL;
  cmdargs->subject = NULL;
  cmdargs->DoCoordDither = 1;
  cmdargs->DoIntensityDither = 1;
  cmdargs->dof = 6;
  for(n=0; n < 12; n++) cmdargs->params[n] = 0;
  for(n=6; n < 9; n++) cmdargs->params[n] = 1;
  cmdargs->nsep = 0;
  cmdargs->DoInitCostOnly = 0;
  cmdargs->DoSmoothing = 1;
  cmdargs->cras0 = 1;
  cmdargs->nitersmax = 4;
  cmdargs->ftol = 10e-8;
  cmdargs->linmintol = .001;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  fsenv = FSENVgetenv();
  check_options();
  if(checkoptsonly) return(0);
  dump_options(stdout);

  coreg = (COREG *) calloc(sizeof(COREG),1);

  printf("Reading in mov %s\n",cmdargs->mov);
  coreg->mov = MRIread(cmdargs->mov);
  if(!coreg->mov) exit(1);
  coreg->nvoxmov = coreg->mov->width * coreg->mov->height * coreg->mov->depth;

  printf("Reading in ref %s\n",cmdargs->ref);  
  coreg->ref = MRIread(cmdargs->ref);
  if(!coreg->ref) exit(1);

  if(cmdargs->refmask){
    printf("Reading in refmask %s\n",cmdargs->refmask);
    coreg->refmask = MRIread(cmdargs->refmask);
    if(!coreg->refmask) exit(1);
    int c,r,s,m;
    for(c=0; c < coreg->ref->width; c++){
      for(r=0; r < coreg->ref->width; r++){
	for(s=0; s < coreg->ref->width; s++){
	  m = MRIgetVoxVal(coreg->refmask,c,r,s,0);
	  if(m < 0.5) MRIsetVoxVal(coreg->ref,c,r,s,0,0.0);
	}
      }
    }
  }

  if(cmdargs->movmask){
    printf("Reading in movmask %s\n",cmdargs->movmask);
    coreg->movmask = MRIread(cmdargs->movmask);
    if(!coreg->movmask) exit(1);
    int c,r,s,m;
    for(c=0; c < coreg->mov->width; c++){
      for(r=0; r < coreg->mov->width; r++){
	for(s=0; s < coreg->mov->width; s++){
	  m = MRIgetVoxVal(coreg->movmask,c,r,s,0);
	  if(m < 0.5) MRIsetVoxVal(coreg->mov,c,r,s,0,0.0);
	}
      }
    }
  }

  if(cmdargs->cras0){
    printf("Setting cras translation parameters to align centers\n");
    MATRIX *Vref, *Vmov, *invVref, *M;
    Vref = MRIxfmCRS2XYZ(coreg->ref, 0);
    Vmov = MRIxfmCRS2XYZ(coreg->mov, 0);
    invVref = MatrixInverse(Vref,NULL);
    M = MatrixMultiplyD(Vmov,invVref,NULL);
    cmdargs->params[0] = M->rptr[1][4];
    cmdargs->params[1] = M->rptr[2][4];
    cmdargs->params[2] = M->rptr[3][4];
    MatrixFree(&Vref);
    MatrixFree(&Vmov);
    MatrixFree(&invVref);
    MatrixFree(&M);
  }

  coreg->histfwhm[0] = 7;
  coreg->histfwhm[1] = 7;
  coreg->nparams = cmdargs->dof;
  coreg->nitersmax = cmdargs->nitersmax;
  coreg->ftol = cmdargs->ftol;
  coreg->linmintol = cmdargs->linmintol;
  coreg->SatPct = 99.99;
  coreg->DoCoordDither = cmdargs->DoCoordDither;
  coreg->DoIntensityDither = cmdargs->DoIntensityDither;
  coreg->refirfs = NULL;
  coreg->movirfs = NULL;
  coreg->DoSmoothing = cmdargs->DoSmoothing;

  if(coreg->DoCoordDither){
    // Creating a dither volume is needed for thread safety
    printf("Creating random numbers for coordinate dithering\n");
    coreg->crfs = RFspecInit(53,NULL);
    coreg->crfs->name = strcpyalloc("uniform");
    coreg->cdither = MRIallocSequence(coreg->ref->width, coreg->ref->height, coreg->ref->depth, MRI_FLOAT,3);
    RFsynth(coreg->cdither, coreg->crfs, NULL);
  } else printf("NOT Creating random numbers for coordinate dithering\n");
  if(coreg->DoIntensityDither){
    printf("Performing intensity dithering\n");
    coreg->refirfs = RFspecInit(53,NULL);
    coreg->refirfs->name = strcpyalloc("uniform");
    coreg->movirfs = RFspecInit(53,NULL);
    coreg->movirfs->name = strcpyalloc("uniform");
  } else printf("NOT Performing intensity dithering\n");
  fflush(stdout);

  // Initial parameters
  for(n=0; n < 12; n++) coreg->params[n] = cmdargs->params[n];
  printf("Initial parameters ");
  for(n=0; n < 12; n++) printf("%7.4lf ",coreg->params[n]);
  printf("\n");

  coreg->nsep = cmdargs->nsep;
  coreg->sepmin = 100;
  for(n=0; n < cmdargs->nsep; n++){
    coreg->seplist[n] = cmdargs->seplist[n];
    if(coreg->seplist[n] < coreg->sepmin) coreg->sepmin = coreg->seplist[n];
  }

  COREGprint(stdout, coreg);

  COREGpreproc(coreg);

  printf("Testing if mov and target overlap\n");
  coreg->sep = 4;
  COREGcost(coreg);
  printf("Init cost %15.10lf\n",coreg->cost);
  printf("nhits = %d out of %d, Percent Overlap: %5.1f\n",coreg->nhits,coreg->nvoxmov,coreg->pcthits);
  if(coreg->nhits == 0){
    printf("ERROR: mov and ref do not overlap\n");
    exit(1);
  }
  printf("Initial  RefRAS-to-MovRAS\n");
  MatrixPrint(stdout,coreg->M);
  printf("Initial  RefVox-to-MovVox\n");
  MatrixPrint(stdout,coreg->V2V);
  if(cmdargs->DoInitCostOnly) exit(0);

  for(n=0; n < coreg->nsep; n++){
    coreg->sep = coreg->seplist[n];
    printf("sep = %d -----------------------------------\n",coreg->sep);
    coreg->startmin = 1;
    COREGMinPowell();
  }

  MATRIX *invV2V;
  invV2V = MatrixInverse(coreg->V2V,NULL);
  lta = LTAcreate(coreg->mov, coreg->ref, invV2V, LINEAR_VOX_TO_VOX);
  if(cmdargs->subject)  sprintf(lta->subject,cmdargs->subject);
  else                  sprintf(lta->subject,"unknown");
  err = LTAwrite(lta,cmdargs->outreg);
  if(err) exit(1);

  if(cmdargs->regdat){
    LTAchangeType(lta,REGISTER_DAT);
    err = LTAwrite(lta,cmdargs->regdat);
    if(err) exit(1);
  }

  printf("Final  RefRAS-to-MovRAS\n");
  MatrixPrint(stdout,coreg->M);
 printf("Final  RefVox-to-MovVox\n");
  MatrixPrint(stdout,coreg->V2V);

  printf("mri_coreg RunTimeSec %4.1f sec\n",TimerStop(&timer)/1000.0);
  printf("mri_coreg done\n\n");

  printf("To check run:\n");
  printf("   tkregister --mov %s --targ %s --reg %s ",cmdargs->mov,cmdargs->ref,cmdargs->outreg);
  if(cmdargs->subject) printf("--s %s --surfs ",cmdargs->subject);
  printf("\n\n");

  exit(0);
}

/* -------------------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--init-cost-only")) cmdargs->DoInitCostOnly = 1;
    else if (!strcasecmp(option, "--no-smooth")) cmdargs->DoSmoothing = 0;
    else if (!strcasecmp(option, "--cras0")) cmdargs->cras0 = 1;
    else if (!strcasecmp(option, "--no-cras0"))  cmdargs->cras0 = 0;
    else if (!strcasecmp(option, "--regheader")) cmdargs->cras0 = 0;

    else if (!strcasecmp(option, "--mov")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->mov = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ref") || !strcasecmp(option, "--targ")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->ref = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ref-mask")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->refmask = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mov-mask")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->movmask = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--reg") || !strcasecmp(option, "--lta")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->outreg = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--regdat")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->regdat = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--s")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->subject = pargv[0];
      nargsused = 1;
    } 
    else if(!strcmp(option, "--sd") || !strcmp(option, "-SDIR")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--dof")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&cmdargs->dof);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--6")) cmdargs->dof = 6;
    else if (!strcasecmp(option, "--9")) cmdargs->dof = 9;
    else if (!strcasecmp(option, "--12")) cmdargs->dof = 12;
    else if (!strcasecmp(option, "--nitersmax")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&cmdargs->nitersmax);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mat2par")) {
      if(nargc < 1) CMDargNErr(option,1);
      LTA *lta;
      lta = LTAread(pargv[0]);
      if(lta==NULL) exit(1);
      LTAchangeType(lta,LINEAR_RAS_TO_RAS);
      LTAinvert(lta,lta);
      COREGparams9(lta->xforms[0].m_L,NULL);
      nargsused = 1;
      exit(0);
    } 
    else if (!strcasecmp(option, "--ftol")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->ftol);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--linmintol")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->linmintol);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--sep")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&cmdargs->seplist[cmdargs->nsep]);
      cmdargs->nsep++;
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--trans")) {
      if(nargc < 3) CMDargNErr(option,3);
      sscanf(pargv[0],"%lf",&cmdargs->params[0]);
      sscanf(pargv[1],"%lf",&cmdargs->params[1]);
      sscanf(pargv[2],"%lf",&cmdargs->params[2]);
      cmdargs->cras0 = 0;
      nargsused = 3;
    } 
    else if (!strcasecmp(option, "--rot")) {
      if(nargc < 3) CMDargNErr(option,3);
      sscanf(pargv[0],"%lf",&cmdargs->params[3]);
      sscanf(pargv[1],"%lf",&cmdargs->params[4]);
      sscanf(pargv[2],"%lf",&cmdargs->params[5]);
      nargsused = 3;
    } 
    else if (!strcasecmp(option, "--scale")) {
      if(nargc < 3) CMDargNErr(option,3);
      sscanf(pargv[0],"%lf",&cmdargs->params[6]);
      sscanf(pargv[1],"%lf",&cmdargs->params[7]);
      sscanf(pargv[2],"%lf",&cmdargs->params[8]);
      nargsused = 3;
    } 
    else if (!strcasecmp(option, "--shear")) {
      if(nargc < 3) CMDargNErr(option,3);
      sscanf(pargv[0],"%lf",&cmdargs->params[9]);
      sscanf(pargv[1],"%lf",&cmdargs->params[10]);
      sscanf(pargv[2],"%lf",&cmdargs->params[11]);
      nargsused = 3;
    } 
    else if(!strcmp(option, "--no-coord-dither")) cmdargs->DoCoordDither = 0;  
    else if(!strcmp(option, "--no-intensity-dither"))   cmdargs->DoIntensityDither = 0;
    else if(!strcmp(option, "--no-dither")){
      cmdargs->DoCoordDither = 0;  
      cmdargs->DoIntensityDither = 0;
    }
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      int nthreads;
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diag")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diag-show"))    Gdiag = (Gdiag & DIAG_SHOW);
    else if (!strcasecmp(option, "--diag-verbose")) Gdiag = (Gdiag & DIAG_VERBOSE);
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* -------------------------------------------------------- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -------------------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --mov movvol : source volume\n");
  printf("   --ref refvol : target volume (can use --targ too)\n");
  printf("   --reg reg.lta : output registration (can use --lta too)\n");
  printf("\n");
  printf("   --s subject \n");
  printf("   --dof DOF : default is %d (also: --6, --9, --12)\n",cmdargs->dof);
  printf("   --ref-mask refmaskvol : mask ref with refmaskvol\n");
  printf("   --mov-mask movmaskvol : mask ref with movmaskvol\n");
  printf("   --threads nthreads\n");
  printf("   --sd SUBJECTS_DIR \n");
  printf("   --regdat reg.dat \n");
  printf("   --no-coord-dither: turn off coordinate dithering\n");
  printf("   --no-intensity-dither: turn off intensity dithering\n");
  printf("   --sep voxsep1 <--sep voxsep2> : set spatial scales (def is 2 and 4)\n");
  printf("   --trans Tx Ty Tz : initial translation in mm (implies --no-cras0)\n");
  printf("   --rot   Rx Ry Rz : initial rotation in deg\n");
  printf("   --scale Sx Sy Sz : initial scale\n");
  printf("   --shear Hxy Hxz Hyz : initial shear\n");
  printf("   --no-cras0 : do not Sett translation parameters to align centers of mov and ref\n");
  printf("   --regheader : same as no-cras0\n");
  printf("   --nitersmax n : default is %d\n",cmdargs->nitersmax);
  printf("   --ftol ftol : default is %5.3le\n",cmdargs->ftol);
  printf("   --linmintol linmintol : default is %5.3le\n",cmdargs->linmintol);
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* -------------------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("This is a program that performs a linear registration between\n");
  printf("two volumes in a way that should more-or-less replicate spm_coreg\n");
  printf("as called by the FreeSurfer program spmregister. \n");
  exit(1) ;
}
/* -------------------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* -------------------------------------------------------- */
static void check_options(void) {
  char tmpstr[3000];
  if(cmdargs->mov == NULL){
    printf("ERROR: must spec --mov\n");
    exit(1);
  }
  if(cmdargs->subject){
    if(cmdargs->ref == NULL) cmdargs->ref = "brainmask.mgz";
    if(!fio_FileExistsReadable(cmdargs->ref)){
      sprintf(tmpstr,"%s/%s/mri/%s",fsenv->SUBJECTS_DIR,cmdargs->subject,cmdargs->ref);
      cmdargs->ref = strcpyalloc(tmpstr);
    }
    if(cmdargs->refmask && !fio_FileExistsReadable(cmdargs->refmask)){
      sprintf(tmpstr,"%s/%s/mri/%s",fsenv->SUBJECTS_DIR,cmdargs->subject,cmdargs->refmask);
      cmdargs->refmask = strcpyalloc(tmpstr);
    }
  }
  if(cmdargs->ref == NULL){
    printf("ERROR: must spec --ref\n");
    exit(1);
  }
  if(cmdargs->outreg == NULL){
    printf("ERROR: must spec --reg\n");
    exit(1);
  }
  if(cmdargs->nsep == 0){
    cmdargs->nsep = 2;
    cmdargs->seplist[0] = 4;
    cmdargs->seplist[1] = 2;
  }
  return;
}
/* -------------------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"cras0    %d\n",cmdargs->cras0);
  return;
}

long COREGvolIndex(int ncols, int nrows, int nslices, int c, int r, int s)
{
  // Index calc must be consistent with MRItoUCharVect(), 0-based, col fastest
  long i;
  i = c + (long)r*ncols + (long)s*nrows*ncols;
  return(i);
}

/*!
  \fn double COREGsamp(unsigned char *f, const double c, const double r, const double s, 
                       const int ncols, const int nrows, const int nslices)
  \brief Trilinear interpolation		       
 */
double COREGsamp(unsigned char *f, const double c, const double r, const double s, 
		 const int ncols, const int nrows, const int nslices)	
{
  int cm,rm,sm,cp,rp,sp;
  double val,cmd,rmd,smd,cpd,rpd,spd;

  cm = floor(c);
  rm = floor(r);
  sm = floor(s);

  cp = ceil(c);
  rp = ceil(r);
  sp = ceil(s);

  cmd = c - cm ;
  rmd = r - rm ;
  smd = s - sm ;
  cpd = (1.0 - cmd) ;
  rpd = (1.0 - rmd) ;
  spd = (1.0 - smd) ;

  val =
    cpd * rpd * spd * f[COREGvolIndex(ncols,nrows,nslices, cm, rm, sm)] +
    cpd * rpd * smd * f[COREGvolIndex(ncols,nrows,nslices, cm, rm, sp)] +
    cpd * rmd * spd * f[COREGvolIndex(ncols,nrows,nslices, cm, rp, sm)] +
    cpd * rmd * smd * f[COREGvolIndex(ncols,nrows,nslices, cm, rp, sp)] +
    cmd * rpd * spd * f[COREGvolIndex(ncols,nrows,nslices, cp, rm, sm)] +
    cmd * rpd * smd * f[COREGvolIndex(ncols,nrows,nslices, cp, rm, sp)] +
    cmd * rmd * spd * f[COREGvolIndex(ncols,nrows,nslices, cp, rp, sm)] +
    cmd * rmd * smd * f[COREGvolIndex(ncols,nrows,nslices, cp, rp, sp)] ;

  return(val);
}


/*!
  \fn int COREGhist(COREG *coreg)
  \brief Compute joint histogram. Somewhat based on spm_hist2.c
 */
int COREGhist(COREG *coreg)
{
  int n,c,r,cref,k,nthreads;
  long nhits;
  double V2V[16],**HH;

  // Pack vox2voxl matrix into an array for speed
  V2V[0] = coreg->V2V->rptr[1][1];
  V2V[1] = coreg->V2V->rptr[2][1];
  V2V[2] = coreg->V2V->rptr[3][1];
  V2V[3] = 0;
  V2V[4] = coreg->V2V->rptr[1][2];
  V2V[5] = coreg->V2V->rptr[2][2];
  V2V[6] = coreg->V2V->rptr[3][2];
  V2V[7] = 0;
  V2V[8]  = coreg->V2V->rptr[1][3];
  V2V[9]  = coreg->V2V->rptr[2][3];
  V2V[10] = coreg->V2V->rptr[3][3];
  V2V[11] = 0;
  V2V[12] = coreg->V2V->rptr[1][4];
  V2V[13] = coreg->V2V->rptr[2][4];
  V2V[14] = coreg->V2V->rptr[3][4];
  V2V[15] = 0;

  nthreads = 1;
  #ifdef _OPENMP
  nthreads = omp_get_max_threads();
  #endif
  HH = (double **)calloc(sizeof(double*),nthreads);
  for(n=0; n < nthreads; n++) 
    HH[n] = (double *)calloc(sizeof(double),256*256);

  nhits = 0;
  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:nhits)
  #endif
  for(cref=0; cref < coreg->ref->width; cref += coreg->sep){
    int rref,sref;
    double dcref,drref,dsref;
    double dcmov,drmov,dsmov;
    double vf, vg;
    int   ivf, ivg;
    double *H;
    int threadno = 1;
    //int iran;

    #ifdef _OPENMP
    threadno = omp_get_thread_num(); 
    #endif
    H = HH[threadno];

    for(rref=0; rref < coreg->ref->height; rref += coreg->sep){
      for(sref=0; sref < coreg->ref->depth; sref += coreg->sep){
	dcref  = cref;
	drref  = rref;
	dsref  = sref;

	if(coreg->DoCoordDither){
	  dcref += MRIgetVoxVal(coreg->cdither,cref,rref,sref,0);
	  drref += MRIgetVoxVal(coreg->cdither,cref,rref,sref,1);
	  dsref += MRIgetVoxVal(coreg->cdither,cref,rref,sref,2);
	}

	dcmov  = V2V[0]*dcref + V2V[4]*drref + V2V[ 8]*dsref +  V2V[12];
	drmov  = V2V[1]*dcref + V2V[5]*drref + V2V[ 9]*dsref +  V2V[13];
	dsmov  = V2V[2]*dcref + V2V[6]*drref + V2V[10]*dsref +  V2V[14];

	if(dcmov < 0 || dcmov > coreg->mov->width-1)  continue;
	if(drmov < 0 || drmov > coreg->mov->height-1) continue;
	if(dsmov < 0 || dsmov > coreg->mov->depth-1)  continue;
	nhits ++;

	vg = COREGsamp(coreg->g, dcref, drref, dsref, coreg->ref->width,coreg->ref->height,coreg->ref->depth);
	vf = COREGsamp(coreg->f, dcmov, drmov, dsmov, coreg->mov->width,coreg->mov->height,coreg->mov->depth);

	ivf = floor(vf);
	ivg = floor(vg+0.5);
	H[ivf+ivg*256] += (1-(vf-ivf));
	if(ivf<255) H[ivf+1+ivg*256] += (vf-ivf);
      }
    }
  }

  // Collect the threads
  for(k=0; k < 256*256; k++){
    coreg->H01d[k] = 0;
    for(n=0; n < nthreads; n++) coreg->H01d[k] += HH[n][k];
  }
  for(n=0; n < nthreads; n++) free(HH[n]);
  free(HH); HH=NULL;

  // Repackage Histogram into a 2D array
  if(!coreg->H0) coreg->H0 = AllocDoubleMatrix(256,256);
  n = 0;
  for(c=0; c < 256; c++){
    for(r=0; r < 256; r++){
      coreg->H0[r][c] = coreg->H01d[n];
      n++;
    }
  }

  // This is good for computing whether and how much the mov and ref overlap
  coreg->nhits = nhits;
  coreg->pcthits = pow(coreg->sep,3)*(double) 100.0*nhits/coreg->nvoxmov;

  return(nhits);
}


/*!
  \fn MRIrescaleToUChar(MRI *mri, MRI *ucmri, double sat)
  \brief rescales to uchar range , but ucmri is actually a FLOAT
 */
MRI *MRIrescaleToUChar(MRI *mri, MRI *ucmri, double sat)
{
  int c,r,s;
  double min, max, v;
  unsigned char ucv;

  if(ucmri == NULL){
    ucmri = MRIalloc(mri->width,mri->height,mri->depth,MRI_FLOAT) ;
    MRIcopyHeader(mri,ucmri);
  }

  max = -10e10;
  min = +10e10;
  for(c=0; c < mri->width; c++){
    for(r=0; r < mri->height; r++){
      for(s=0; s < mri->depth; s++){
	v = MRIgetVoxVal(mri,c,r,s,0);
	if(v > max) max = v;
	if(v < min) min = v;
      }
    }
  }
  //printf("min = %lf max = %lf\n",min,max);
  if(max > sat) max = sat;

  for(c=0; c < mri->width; c++){
    for(r=0; r < mri->height; r++){
      for(s=0; s < mri->depth; s++){
	v = MRIgetVoxVal(mri,c,r,s,0);
	if(v > sat) v = sat;
	ucv = nint(255*(v-min)/(max-min));
	MRIsetVoxVal(ucmri,c,r,s,0,ucv);
      }
    }
  }
  return(ucmri);
}

/*!
  \fn unsigned char *MRItoUCharVect(MRI *mri)
  Converts mri values to a uchar vector.
  Important! Must be consistent withh COREGvolIndex()
*/
unsigned char *MRItoUCharVect(MRI *mri, RFS *rfs)
{
  int c,r,s,nvox;
  unsigned char *a, *pa;
  float val,dval=0.0;

  nvox = mri->width*mri->height*mri->depth;
  a = (unsigned char *) calloc(sizeof(unsigned char),nvox);
  pa = a;
  for(s=0; s < mri->depth; s++){
    for(r=0; r < mri->height; r++){
      for(c=0; c < mri->width; c++){
	if(rfs) dval = RFdrawVal(rfs);
	val = MRIgetVoxVal(mri,c,r,s,0)+dval;
	if(val < 0)   val = 0;
	if(val > 255) val = 255;
	*pa = (unsigned char) nint(val);
	pa++;
      }
    }
  }
  return(a);
}

MATRIX *MRIgetVoxelToVoxelXformBase(MRI *mri_src, MRI *mri_dst, MATRIX *SrcRAS2DstRAS, MATRIX *SrcVox2DstVox, int base)
{
  MATRIX *m_ras2vox_dst, *m_vox2ras_src, *m_vox2ras_dst ;
  m_vox2ras_src = MRIxfmCRS2XYZ(mri_src, base);
  m_vox2ras_dst = MRIxfmCRS2XYZ(mri_dst, base);
  m_ras2vox_dst = MatrixInverse(m_vox2ras_dst,NULL);
  if(SrcRAS2DstRAS){
    SrcVox2DstVox = MatrixMultiplyD(m_ras2vox_dst, SrcRAS2DstRAS, SrcVox2DstVox) ;
    MatrixMultiplyD(SrcVox2DstVox, m_vox2ras_src, SrcVox2DstVox);
  } 
  else  SrcVox2DstVox = MatrixMultiplyD(m_ras2vox_dst, m_vox2ras_src, SrcVox2DstVox) ;

  MatrixFree(&m_vox2ras_dst);
  MatrixFree(&m_vox2ras_src) ;
  MatrixFree(&m_ras2vox_dst) ;
  return(SrcVox2DstVox) ;
}

/*!
  \fn double **AllocDoubleMatrix(int rows, int cols)
  Just a simple routine to allocate a matrix
*/
double **AllocDoubleMatrix(int rows, int cols)
{
  double **M;
  int r;
  M = (double **) calloc(rows,sizeof(double*));
  for(r=0; r < rows; r++)
    M[r] = (double *) calloc(cols,sizeof(double));
  return(M);
}

int FreeDoubleMatrix(double **M, int rows, int cols)
{
  int r;
  for(r=0; r < rows; r++) free(M[r]);
  free(M);
  return(0);
}
int WriteDoubleMatrix(char *fname, char *fmt, double **M, int rows, int cols)
{
  FILE *fp;
  fp = fopen(fname,"w");
  PrintDoubleMatrix(fp, fmt, M, rows, cols);
  fclose(fp);
  return(0);
}
int PrintDoubleMatrix(FILE *fp, char *fmt, double **M, int rows, int cols)
{
  int r,c;
  if(fmt==NULL) fmt = "%20.10lf";
  for(r=0; r < rows; r++){
    for(c=0; c < cols; c++) fprintf(fp,fmt,M[r][c]);
    fprintf(fp,"\n");
  }
  return(0);
}

/*!
  \fn double *SumVectorDoubleMatrix(double **M, int rows, int cols, int dim, double *sumvect, int *nv)
  Compute a vector that is the sum of either the rows or he columns of the given matrix
*/
double *SumVectorDoubleMatrix(double **M, int rows, int cols, int dim, double *sumvect, int *nv)
{
  int r,c;

  if(dim == 1){
    if(!sumvect) sumvect = (double *) calloc(sizeof(double),rows);
    *nv=rows;
    for(r=0; r < rows; r++){
      sumvect[r] = 0;
      for(c=0; c < cols; c++) sumvect[r] += M[r][c];
    }
  }
  else {
    if(!sumvect) sumvect = (double *) calloc(sizeof(double),cols);
    *nv=cols;
    for(c=0; c < cols; c++){
      sumvect[c] = 0;
      for(r=0; r < rows; r++) sumvect[c] += M[r][c];
    }
  }
  return(sumvect);
}


double *conv1d(double *v1, int n1, double *v2, int n2, double *vc)
{
  int n, ntot = n1+n2-1;
  if(vc == NULL) vc = (double*) calloc(ntot,sizeof(double));

  //To OMP this code, must make vc thread safe
  for(n=0; n < ntot; n++){
    int k, kmin, kmax;
    if(n >= n2-1) kmin = n-(n2-1);
    else          kmin = 0;
    if(n < n1-1)  kmax = n;
    else          kmax = n1-1;
    vc[n] = 0;
    for(k=kmin; k <= kmax; k++) vc[n] += v1[k]*v2[n-k];
  }
  //for(n=0; n < ntot; n++)  printf("%3d  %10.8lf\n",n,vc[n]);

  return(vc);
}

double **conv1dmat(double **M, int rows, int cols, double *v, int nv, int dim, 
		   double **C, int *pcrows, int *pccols)
{
  int crows, ccols, c, r, ntot;
  double *vc,*Mv;

  if(dim == 1) {
    // convolve with column vector
    crows = rows+nv-1;
    ccols = cols;
    ntot = crows;
  }
  else{
    // convolve with row vector
    crows = rows;
    ccols = cols+nv-1;
    ntot = ccols;
  }
  //printf("M size %d,%d, nv=%d, C size %d %d\n",rows,cols,nv,crows,cols);
  if(C==NULL) C = AllocDoubleMatrix(crows, ccols);

  if(dim == 1) {
    // convolve v with column vector
    Mv = (double *) calloc(rows,sizeof(double));
    vc = (double *) calloc(ntot,sizeof(double));
    for(c=0; c < cols; c++){
      for(r=0; r < rows; r++) Mv[r] = M[r][c];
      conv1d(Mv, rows, v, nv, vc);
      for(r=0; r < crows; r++) C[r][c] = vc[r];
    }
    free(Mv);
    free(vc);
  }  
  if(dim == 2) {
    // convolve v with row vector
    for(r=0; r < rows; r++) conv1d(M[r], cols, v, nv, C[r]);
  }  

  *pcrows = crows;
  *pccols = ccols;
  return(C);
}

/*!
  \fn int conv1dmatTest(void)
  This is a routine to test conv1dmat and, by extension, conv1d.
  The idea is to run this function to generate several output files,
  then run the code below in matlab.
    m = load('m.dat');
    v = load('v.dat');
    c1 = load('c1.dat');
    c2 = load('c2.dat');
    c1a = conv2(m,v');
    c2a = conv2(c1a,v);
    d1 = max(abs(c1(:)-c1a(:)));
    d2 = max(abs(c2(:)-c2a(:)));
    % d1 and d2 should be very small < 10e-7
 */
int conv1dmatTest(void)
{
  double **M, **C1, **C2, *v;
  int rows=30, cols=200, nv=40;
  int C1rows, C1cols, C2rows, C2cols;
  int c,r;
  FILE *fp;

  M = AllocDoubleMatrix(rows, cols);
  for(r=0; r < rows; r++)
    for(c=0; c < cols; c++) M[r][c] = drand48()-0.5;

  v = (double *) calloc(nv,sizeof(double));
  for(c=0; c < nv; c++) v[c] = drand48()-0.5;  

  fp = fopen("m.dat","w");
  for(r=0; r < rows; r++){
    for(c=0; c < cols; c++) fprintf(fp,"%20.10lf",M[r][c]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  fp = fopen("v.dat","w");
  for(c=0; c < nv; c++) fprintf(fp,"%20.10lf\n",v[c]);
  fclose(fp);

  C1 = conv1dmat(M, rows, cols, v, nv, 2, NULL, &C1rows, &C1cols);
  WriteDoubleMatrix("c1.dat", "%20.10lf", C1, C1rows, C1cols);

  C2 = conv1dmat(C1, C1rows, C1cols, v, nv, 1, NULL,&C2rows, &C2cols);
  WriteDoubleMatrix("c2.dat", "%20.10lf", C2, C2rows, C2cols);

  return(0);
}

double NMICost(double **H, int cols, int rows)
{
  double *s1, *s2;
  int ns1, ns2,n,c,r;
  double den,cost,s1sum,s2sum;

  s1 = SumVectorDoubleMatrix(H, rows, cols, 1, NULL, &ns1);
  s2 = SumVectorDoubleMatrix(H, rows, cols, 2, NULL, &ns2);

  den = 0;
  for(c=0; c < cols; c++){
    for(r=0; r < rows; r++){
      den += (H[r][c]*log2(H[r][c]));
    }
  }

  s1sum = 0;
  for(n=0; n < ns1; n++) s1sum += (s1[n]*log2(s1[n]));
  s2sum = 0;
  for(n=0; n < ns2; n++) s2sum += (s2[n]*log2(s2[n]));
  
  cost = -(s1sum+s2sum)/(den+FLT_EPSILON);
  //printf("%20.15lf %20.15lf %20.15lf \n",s1sum,s2sum,cost);

  free(s1); s1=NULL;
  free(s2); s2=NULL;
  return(cost);
}

/*!
  \fn MATRIX *COREGmatrix(double *p0, int np, MATRIX *M)
  \brief Computes a RAS-to-RAS transformation matrix given
  the parameters. p0[0-2] translation, p0[3-5] rotation
  in degrees, p0[6-8] scale, p[9-11] shear. 
  M = T*R1*R2*R3*SCALE*SHEAR
  Consistent with COREGparams9()
 */
MATRIX *COREGmatrix(double *p0, int np, MATRIX *M)
{
  MATRIX *T, *R1, *R2, *R3, *R, *ZZ, *S;
  double p[12];
  int n;
  for(n=0; n<12; n++) p[n] = 0;
  p[6] = p[7] = p[8] = 1; // scaling
  for(n=0; n<np; n++) p[n] = p0[n];

  //printf("p = [");
  //for(n=0; n<np; n++) printf("%10.3lf ",p[n]);
  //printf("];\n");

  // translations
  T = MatrixIdentity(4,NULL);
  T->rptr[1][4] = p[0];
  T->rptr[2][4] = p[1];
  T->rptr[3][4] = p[2];

  // rotations
  R1 = MatrixIdentity(4,NULL);
  R1->rptr[2][2] = cos(p[3]*M_PI/180);
  R1->rptr[2][3] = sin(p[3]*M_PI/180);
  R1->rptr[3][2] = -sin(p[3]*M_PI/180);
  R1->rptr[3][3] = cos(p[3]*M_PI/180);

  R2 = MatrixIdentity(4,NULL);
  R2->rptr[1][1] = cos(p[4]*M_PI/180);
  R2->rptr[1][3] = sin(p[4]*M_PI/180);
  R2->rptr[3][1] = -sin(p[4]*M_PI/180);
  R2->rptr[3][3] = cos(p[4]*M_PI/180);

  R3 = MatrixIdentity(4,NULL);
  R3->rptr[1][1] = cos(p[5]*M_PI/180);
  R3->rptr[1][2] = sin(p[5]*M_PI/180);
  R3->rptr[2][1] = -sin(p[5]*M_PI/180);
  R3->rptr[2][2] = cos(p[5]*M_PI/180);

  R = MatrixMultiplyD(R1,R2,NULL);
  MatrixMultiplyD(R,R3,R);

  // scale, use ZZ because some idiot #defined Z
  ZZ = MatrixIdentity(4,NULL);
  ZZ->rptr[1][1] = p[6];
  ZZ->rptr[2][2] = p[7];
  ZZ->rptr[3][3] = p[8];
  ZZ->rptr[4][4] = 1;

  // shear
  S = MatrixIdentity(4,NULL);
  S->rptr[1][2] = p[9];
  S->rptr[1][3] = p[10];
  S->rptr[2][3] = p[11];

  // M = T*R*ZZ*S
  M = MatrixMultiplyD(T,R,M);
  MatrixMultiplyD(M,ZZ,M);
  MatrixMultiplyD(M,S,M);
  //MatrixPrint(stdout,M);

  MatrixFree(&T);
  MatrixFree(&R1);
  MatrixFree(&R2);
  MatrixFree(&R3);
  MatrixFree(&R);
  MatrixFree(&ZZ);
  MatrixFree(&S);

  return(M);
}

/*!
  \fn double *COREGparams9(MATRIX *M9, double *p)
  \brief Extracts parameter from a 9 dof transformation matrix.
  This is consistent with COREGmatrix(). Still need to figure
  out how to do 12 dof. Note: p will be alloced to 12.
  Angles are in degrees.
 */
double *COREGparams9(MATRIX *M9, double *p)
{
  double sum;
  int n, c, r;
  MATRIX *R;

  if(p==NULL) p = (double*) calloc(12,sizeof(double));

  // translation
  for(r=0; r < 3; r++) p[r] = M9->rptr[r+1][4];

  // R is the rotation matrix
  R = MatrixAlloc(3,3,MATRIX_REAL);
  for(c=0; c < 3; c++){
    sum = 0;
    for(r=0; r < 3; r++) sum += (M9->rptr[r+1][c+1]*M9->rptr[r+1][c+1]);
    p[c+6] = sqrt(sum); //scale
    for(r=0; r < 3; r++) R->rptr[r+1][c+1] = M9->rptr[r+1][c+1]/sqrt(sum);
  }

  // extract rotation params
  p[3] = atan2(R->rptr[2][3],R->rptr[3][3])*180/M_PI;
  p[4] = atan2(R->rptr[1][3],sqrt(pow(R->rptr[2][3],2) + pow(R->rptr[3][3],2)))*180/M_PI;
  p[5] = atan2(R->rptr[1][2],R->rptr[1][1])*180/M_PI;

  MatrixFree(&R);

  for(n=0; n < 9; n++) printf("%10.8lf ",p[n]);
  printf("\n");

  return(p);
}


double COREGcost(COREG *coreg)
{
  double **H1,**H;
  double *g1, *g2, sum, std1, std2;
  int r,c,n,lim1,lim2,ng1,ng2;
  int H1rows,H1cols,Hrows,Hcols;

  // RefRAS-to-MovRAS
  coreg->M = COREGmatrix(coreg->params, coreg->nparams, coreg->M);

  // AnatVox-to-FuncVox
  coreg->V2V = MRIgetVoxelToVoxelXformBase(coreg->ref,coreg->mov,coreg->M,coreg->V2V,0);

  // Compute joint histogram
  COREGhist(coreg);

  //printf("M  %20.18f %20.18f %20.18f \n",coreg->M->rptr[1][1],coreg->V2V->rptr[1][1],coreg->H0[0][0]);

  // filter for the column vectors
  std1 = coreg->histfwhm[0]/sqrt(log(256.0));
  lim1 = ceil(2*coreg->histfwhm[0]);
  ng1 = 2*lim1+1;
  g1 = (double *) calloc(ng1,sizeof(double));
  sum = 0;
  for(n=-lim1; n <= lim1; n++){
    g1[n+lim1] = exp(-(n*n)/(2*(std1*std1)))/(std1*sqrt(2*M_PI));
    sum += g1[n+lim1]; 
  }
  for(n=0; n < ng1; n++) {
    g1[n] /= sum;
    //printf("%3d  %lf\n",n,g1[n]);
  }

  // filter for the row vectors
  std2 = coreg->histfwhm[1]/sqrt(log(256.0));
  lim2 = ceil(2*coreg->histfwhm[1]);
  ng2 = 2*lim2+1;
  g2 = (double *) calloc(2*lim2+1,sizeof(double));
  sum = 0;
  for(n=-lim2; n <= lim2; n++){
    g2[n+lim2] = exp(-(n*n)/(2*(std2*std2)))/(std2*sqrt(2*M_PI));
    sum += g2[n+lim2]; 
  }
  for(n=0; n < ng2; n++) g2[n] /= sum;

  // Apply filters
  H1 = conv1dmat(coreg->H0, 256, 256, g2, ng2, 2, NULL,&H1rows,&H1cols);
  H  = conv1dmat(H1, H1rows, H1cols, g1, ng1, 1, NULL,&Hrows,&Hcols);

  for(c=0; c < Hcols; c++) for(r=0; r < Hrows; r++) H[r][c] += FLT_EPSILON;

  sum = 0;
  for(c=0; c < Hcols; c++) for(r=0; r < Hrows; r++) sum += H[r][c];
  for(c=0; c < Hcols; c++) for(r=0; r < Hrows; r++) H[r][c] /= sum;

  coreg->cost = NMICost(H, Hcols, Hrows);

  FreeDoubleMatrix(H1,H1rows,H1cols); H1=NULL;
  FreeDoubleMatrix(H,Hrows,Hcols);    H=NULL;
  free(g1); g1=NULL;
  free(g2); g2=NULL;

  return(coreg->cost);
}


/*--------------------------------------------------------------------------*/
float COREGcostPowell(float *pPowel) 
{
  extern COREG *coreg;
  int n,newmin;
  float curcost;
  static float initcost=-1,mincost=-1,ppmin[100];
  FILE *fp;

  for(n=0; n < coreg->nparams; n++) coreg->params[n] = pPowel[n+1];
  
  // compute cost
  curcost = COREGcost(coreg);

  newmin = 0;
  if(coreg->startmin) {
    newmin = 1;
    initcost = curcost;
    mincost = curcost;
    for(n=0; n<coreg->nparams; n++) ppmin[n] = coreg->params[n];
    printf("InitialCost %20.10lf \n",initcost);
    coreg->startmin = 0;
  }

  if(mincost > curcost) {
    newmin = 1;
    mincost = curcost;
    for(n=0; n<coreg->nparams; n++) ppmin[n] = coreg->params[n];
  }

  if(0){
    fp = stdout;
    fprintf(fp,"%4d  ",coreg->nCostEvaluations);
    for(n=0; n<coreg->nparams; n++) fprintf(fp,"%12.8f ",coreg->params[n]);
    fprintf(fp,"  %12.10f %12.5f\n",curcost/initcost,curcost);
    fflush(fp);
  }

  if(newmin){
    printf("#@# %2d %4d  ",coreg->sep,coreg->nCostEvaluations);
    for(n=0; n<coreg->nparams; n++) printf("%7.5f ",ppmin[n]);
    printf("  %9.7f\n",mincost);
    fflush(stdout);
  }

  coreg->nCostEvaluations++;
  return((float)curcost);
}

/*---------------------------------------------------------*/
int COREGMinPowell()
{
  extern COREG *coreg;
  float *pPowel, **xi;
  int    r, c, n,dof;
  struct timeb timer;

  TimerStart(&timer);
  dof = coreg->nparams;

  printf("\n\n---------------------------------\n");
  printf("Init Powel Params dof = %d\n",dof);
  pPowel = vector(1, dof) ;
  for(n=0; n < dof; n++) pPowel[n+1] = coreg->params[n];

  xi = matrix(1, dof, 1, dof) ;
  for (r = 1 ; r <= dof ; r++) {
    for (c = 1 ; c <= dof ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }
  printf("Starting OpenPowel2(), sep = %d\n",coreg->sep);
  OpenPowell2(pPowel, xi, dof, coreg->ftol, coreg->linmintol, coreg->nitersmax, 
	      &coreg->niters, &coreg->fret, COREGcostPowell);
  printf("Powell done niters total = %d\n",coreg->niters);
  printf("OptTimeSec %4.1f sec\n",TimerStop(&timer)/1000.0);
  printf("OptTimeMin %5.2f min\n",(TimerStop(&timer)/1000.0)/60);
  printf("nEvals %d\n",coreg->nCostEvaluations);
  //printf("EvalTimeSec %4.1f sec\n",(TimerStop(&timer)/1000.0)/coreg->nCostEvaluations);
  fflush(stdout);

  printf("Final parameters ");
  for(n=0; n < coreg->nparams; n++){
    coreg->params[n] = pPowel[n+1];
    printf("%12.8f ",coreg->params[n]);
  }
  printf("\n");

  COREGcost(coreg);
  printf("Final cost %20.15lf\n ",coreg->cost);

  free_matrix(xi, 1, dof, 1, dof);
  free_vector(pPowel, 1, dof);
  printf("\n\n---------------------------------\n");
  return(NO_ERROR) ;
}

int COREGpreproc(COREG *coreg)
{
  int n, DoSmooth;
  MRI *mritmp;

  // Rescale and maybe smooth the moveable
  coreg->movsat = MRIgetPercentile(coreg->mov, coreg->SatPct, 0);
  printf("movsat = %6.4lf\n",coreg->movsat);
  mritmp = MRIrescaleToUChar(coreg->mov,NULL,coreg->movsat);
  COREGfwhm(coreg->mov, coreg->sepmin, coreg->movfwhm);
  DoSmooth = 0;
  for(n=0; n < 3; n++){
    coreg->movgstd[n] = coreg->movfwhm[n]/sqrt(log(256.0));
    if(coreg->movgstd[n]>0) DoSmooth = 1;
  }
  printf("mov gstd %6.4f %6.4f %6.4f\n",coreg->movgstd[0],coreg->movgstd[1],coreg->movgstd[2]);
  if(DoSmooth && coreg->DoSmoothing){
    printf("Smoothing mov\n");
    MRIgaussianSmoothNI(mritmp, coreg->movgstd[0], coreg->movgstd[1], coreg->movgstd[2], mritmp);
  } else printf("NOT Smoothing mov\n");
  if(coreg->f) free(coreg->f);
  coreg->f = MRItoUCharVect(mritmp,coreg->movirfs);
  MRIfree(&mritmp);
  fflush(stdout);

  // Rescale and maybe smooth the reference
  coreg->refsat = MRIgetPercentile(coreg->ref, coreg->SatPct, 0);
  printf("refsat = %6.4lf\n",coreg->refsat);
  mritmp = MRIrescaleToUChar(coreg->ref,NULL,coreg->refsat);
  COREGfwhm(coreg->ref, coreg->sepmin, coreg->reffwhm);
  DoSmooth = 0;
  for(n=0; n < 3; n++){
    coreg->refgstd[n] = coreg->reffwhm[n]/sqrt(log(256.0));
    if(coreg->refgstd[n]>0) DoSmooth = 1;
  }
  printf("ref gstd %6.4f %6.4f %6.4f\n",coreg->refgstd[0],coreg->refgstd[1],coreg->refgstd[2]);
  if(DoSmooth && coreg->DoSmoothing){
    printf("Smoothing ref\n");
    MRIgaussianSmoothNI(mritmp, coreg->refgstd[0], coreg->refgstd[1], coreg->refgstd[2], mritmp);
    //MRIwrite(mritmp,"ref.smoothed.mgh");
  } else printf("NOT Smoothing ref\n");
  if(coreg->g) free(coreg->g);
  coreg->g = MRItoUCharVect(mritmp,coreg->refirfs);
  MRIfree(&mritmp);
  fflush(stdout);

  printf("COREGpreproc() done\n");

  return(0);
}

/*!
  \fn int COREGfwhm(MRI *mri, double sep, double fwhm[3])
  Compute fwhm for smoothing to take into account the separation
 */
int COREGfwhm(MRI *mri, double sep, double fwhm[3])
{
  int n;
  double voxsize[3],val;
  voxsize[0] = mri->xsize;
  voxsize[1] = mri->ysize;
  voxsize[2] = mri->zsize;
  for(n=0; n < 3; n++){
    val = sep*sep-voxsize[n]*voxsize[n];
    if(val > 0) fwhm[n] = 1.15*sqrt(val);
    else        fwhm[n] = 0;
    //fwhm[n] = 0;
  }
  return(0);
}

/*!
  \fn float MRIgetPercentile(MRI *mri, double Pct, int frame)
  Compute the value of the given intensity percentile. Pct 0-100.
 */
float MRIgetPercentile(MRI *mri, double Pct, int frame)
{
  float *vect,val;
  int nvox,c,r,s,n;

  nvox = mri->width * mri->height * mri->depth;
  vect = (float *)calloc(nvox,sizeof(float));

  n = 0;
  for(s=0; s < mri->depth; s++){
    for(c=0; c < mri->width; c++){
      for(r=0; r < mri->height; r++){
	vect[n] = MRIgetVoxVal(mri,c,r,s,frame);
	n++;
      }
    }
  }
  qsort((void *) vect, nvox, sizeof(float), compare_floats);

  n = round(nvox*Pct/100.0);
  if(n==0 || n >= nvox-1) val = vect[n];
  else val = (vect[n]+vect[n-1]+vect[n+1])/3.0;
  free(vect); vect=NULL;
  return(val);
}


int COREGprint(FILE *fp, COREG *coreg)
{
  int n;
  fprintf(fp,"Separation list (%d): ",coreg->nsep);
  for(n=0; n < coreg->nsep; n++)  fprintf(fp,"%2d ",coreg->seplist[n]);
  fprintf(fp,"  min = %d\n",coreg->sepmin);
  fprintf(fp,"DoSmoothing %d\n",coreg->DoSmoothing);
  fprintf(fp,"DoCoordDither %d\n",coreg->DoCoordDither);
  fprintf(fp,"DoIntensityDither %d\n",coreg->DoIntensityDither);
  fprintf(fp,"nitersmax %d\n",coreg->nitersmax);
  fprintf(fp,"ftol %5.3le\n",coreg->ftol);
  fprintf(fp,"linmintol %5.3le\n",coreg->linmintol);
  fprintf(fp,"SatPct %lf\n",coreg->SatPct);
  fprintf(fp,"Hist FWHM %lf %lf\n",coreg->histfwhm[0],coreg->histfwhm[1]);
#ifdef _OPENMP
  fprintf(fp,"nthreads %d\n",omp_get_max_threads());
#else
  fprintf(fp,"nthreads %d\n",1);
#endif
  return(0);
}

