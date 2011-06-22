/**
 * @file  mri_ibmc.c
 * @brief Intersection-based motion correction
 *
 * Intersection-based motion correction based on Kim, et al, TMI,
 * 2010. Registers three volumes.
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2011/06/22 14:27:15 $
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

/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <float.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "fsenv.h"
#include "registerio.h"
#include "matrix.h"
#include "mri2.h"
#include "transform.h"
#include "timer.h"
#include "numerics.h"


double round(double);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_ibmc.c,v 1.8 2011/06/22 14:27:15 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;
MATRIX *MRIangles2RotMatB(double *angles, MATRIX *R);
MRI *MRIextractSlice(MRI *vol, MRI *slice, int Sno)
{
  int width, height, w, h, f;
  MATRIX *Sv,*P0, *crs;
  double v;

  if(Sno >= vol->depth){
    printf("ERROR: MRIextractSlice: slice no %d too large (%d)\n",Sno,vol->depth);
    return(NULL);
  }

  width = vol->width;
  height = vol->height;
  if(slice == NULL){
    slice = MRIallocSequence(width,height, 1, MRI_FLOAT, vol->nframes);
    MRIcopyHeader(vol,slice);
  }
  else {
    if(slice->width != width || slice->height != height) {
      printf("ERROR: MRIextractSlice: size mismatch\n");
      return(NULL);
    }
  }

  for(f=0; f < vol->nframes; f++){
    for(w=0; w < width; w++){
      for(h=0; h < height; h++){
	v = MRIgetVoxVal(vol,w,h,Sno,f);
	MRIsetVoxVal(slice,w,h,0,f,v);
      }
    }
  }

  crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[1][3] = Sno;
  crs->rptr[1][4] = 1;
  Sv = MRIxfmCRS2XYZ(vol, 0);
  P0 = MatrixMultiply(Sv,crs,NULL);
  MRIp0ToCRAS(slice, P0->rptr[1][1], P0->rptr[1][2], P0->rptr[1][3]);
  MatrixFree(&crs);
  MatrixFree(&Sv);
  MatrixFree(&P0);
  
  return(slice);
}

#define IBMC_A2B_NONE 0
#define IBMC_A2B_RIGID 1

/*---------------------------------------------------------------*/
typedef struct 
{
  MRI *mri;
  MATRIX *Rs;
  MATRIX *Rs0;
  double beta[6];
  int Vno, Sno;
  MATRIX *D, *iD;
  MATRIX *Ts;
  MATRIX *iRs;
  MATRIX *iRsTs; // inv(Rs)*Ts
} IBMC_SLICE;
/*---------------------------------------------------------------*/
typedef struct 
{
  MRI *mri;
  MATRIX *Rv;
  MATRIX *Tv,*iTv;
  int nalpha;
  int a2bmethod;
  double *alpha;
  int nbeta;
  double *beta;
  IBMC_SLICE **slice;
  int Vno;
} IBMC_VOL;
/*---------------------------------------------------------------*/
typedef struct 
{
  int nvols;
  IBMC_VOL *vols;
  int a2bmethod;
} IBMC;
/*---------------------------------------------------------------*/
IBMC_SLICE *IBMCinitSlice(IBMC_VOL *vol, int Sno)
{
  IBMC_SLICE *slice;

  slice = (IBMC_SLICE *) calloc(1,sizeof(IBMC_SLICE));
  slice->Sno = Sno;
  slice->mri = MRIcrop(vol->mri, 0,0,Sno, vol->mri->width-1,vol->mri->height-1,Sno);
  slice->Ts  = MRIxfmCRS2XYZtkreg(slice->mri);
  slice->D = MatrixIdentity(4,NULL);
  slice->D->rptr[3][4] = -Sno;
  slice->iD   = MatrixInverse(slice->D,NULL);
  slice->Rs0  = MatrixMultiply(slice->Ts,slice->iD,NULL);
  slice->Rs0  = MatrixMultiply(slice->Rs0,vol->iTv,slice->Rs0);
  slice->Rs    = MatrixAlloc(4,4,MATRIX_REAL);
  slice->iRs   = MatrixAlloc(4,4,MATRIX_REAL);
  slice->iRsTs = MatrixAlloc(4,4,MATRIX_REAL);
  return(slice);
}
/*---------------------------------------------------------------*/
int IBMCsetSliceBeta(IBMC_SLICE *slice, const double beta[6])
{
  int k;
  static double angles[3];

  for(k=0; k < 6; k++) slice->beta[k] = beta[k];
  for(k=0; k < 3; k++) angles[k] = beta[k+3]*M_PI/180;
  slice->Rs = MRIangles2RotMatB(angles, slice->Rs);
  for(k=0; k < 3; k++) slice->Rs->rptr[k+1][4] = beta[k];
  slice->iRs = MatrixInverse(slice->Rs,slice->iRs);
  slice->iRsTs = MatrixMultiply(slice->iRs,slice->Ts,slice->iRsTs);
  slice->iRsTs = MatrixMultiply(slice->Rs,slice->Rs0,slice->Rs);
  return(0);
}
/*---------------------------------------------------------------*/
IBMC_VOL *IBMCinitVol(MRI *mri, MATRIX *Rv, int abmethod)
{
  IBMC_VOL *vol;
  int nthslice;

  vol = (IBMC_VOL *)calloc(1,sizeof(IBMC_VOL));
  vol->mri = mri;
  vol->Rv = Rv;
  vol->slice = (IBMC_SLICE **)calloc(mri->depth,sizeof(IBMC_SLICE *));
  vol->Tv  = MRIxfmCRS2XYZtkreg(vol->mri);
  vol->iTv = MatrixInverse(vol->Tv,NULL);

  for(nthslice = 0; nthslice < vol->mri->depth; nthslice++){
    vol->slice[nthslice] = IBMCinitSlice(vol, nthslice);
    //printf("%2d %3d\n",nthslice,vol->slice[0]->mri->depth);
  }
  return(vol);
}

/*---------------------------------------------------------------*/
int IBMCwriteSlices(IBMC_VOL *vol, char *base)
{
  char tmpstr[2000];
  int nthslice;

  for(nthslice = 0; nthslice < vol->mri->depth; nthslice++){
    sprintf(tmpstr,"%s.%03d.nii",base,nthslice);
    //printf("%2d %3d %s\n",nthslice,vol->slice[nthslice]->mri->depth,tmpstr);
    MRIwrite(vol->slice[nthslice]->mri,tmpstr);
  }
  return(0);
}


/*---------------------------------------------------------------*/
int IBMCalpha2Beta(IBMC_VOL *vol)
{
  int n,s,k,nthbeta;

  switch(vol->a2bmethod){
  case IBMC_A2B_NONE:
    for(n=0; n < vol->nalpha; n++) vol->beta[n] = vol->alpha[n];
    break;
  case IBMC_A2B_RIGID:
    nthbeta = 0;
    for(s=0; s < vol->mri->depth; s++){
      for(k=0; k < 6; k++){
	vol->beta[nthbeta] = vol->alpha[k];
	nthbeta ++;
      }
    }
    break;
  }
  return(0);
}

/*---------------------------------------------------------------*/
int IBMCvol2slices(MRI *src, IBMC_VOL *dst)
{
  int nthslice;
  MRI *slc;
  MATRIX *iRvRs=NULL, *RvRs=NULL;
  for(nthslice=0; nthslice < dst->mri->depth; nthslice++){
    RvRs = MatrixMultiply(dst->Rv,dst->slice[nthslice]->Rs,RvRs);
    iRvRs = MatrixInverse(RvRs,iRvRs);
    printf("iRvRs %d ----------------\n",nthslice);
    MatrixPrint(stdout,dst->slice[nthslice]->Rs);
    slc = dst->slice[nthslice]->mri;
    MRIvol2VolTkRegVSM(src, slc, iRvRs, SAMPLE_TRILINEAR, 0, NULL);
  }
  return(0);
}

/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
MRI *targ=NULL, *mov=NULL;
MATRIX *Rtarg=NULL,*Rmov=NULL;
char *AlphaFile=NULL;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,k,nthslice;
  IBMC_VOL *vol;
  double beta[6];

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
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);

  vol = IBMCinitVol(mov, Rmov, 0);

  for(k=0; k<6; k++) beta[k] = 0;
  for(nthslice=0; nthslice < vol->mri->depth; nthslice++)
    IBMCsetSliceBeta(vol->slice[nthslice], beta);
  IBMCvol2slices(targ, vol);
  IBMCwriteSlices(vol, "base");

  beta[0] = 5;
  for(nthslice=0; nthslice < vol->mri->depth; nthslice++)
    IBMCsetSliceBeta(vol->slice[nthslice], beta);
  IBMCvol2slices(targ, vol);
  IBMCwriteSlices(vol, "base2");


  return 0;
}

/*------------------------------------------------------------------*/
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
    else if (!strcasecmp(option, "--debug"))  debug = 1;
    else if (!strcasecmp(option, "--diag")) Gdiag_no = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--targ")) {
      if(nargc < 1) CMDargNErr(option,1);
      printf("Reading targ\n");
      targ = MRIread(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mov")) {
      if(nargc < 2) CMDargNErr(option,2);
      printf("Reading mov\n");
      mov = MRIread(pargv[0]);
      printf("Reading Rmov\n");
      Rmov = regio_read_registermat(pargv[1]);
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--alpha")) {
      if(nargc < 1) CMDargNErr(option,1);
      AlphaFile = pargv[0];
      nargsused = 1;
    } 
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
/*------------------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*------------------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --v1 vol1 <regfile>: template volume \n");
  printf("   --v2 vol2 <regfile>: template volume \n");
  printf("   --v3 vol3 <regfile>: template volume \n");
  printf("\n");
  printf("  --nmax nmax   : max number of powell iterations (def 36)\n");
  printf("  --tol   tol   : powell inter-iteration tolerance on cost\n");
  printf("       This is the fraction of the cost that the difference in \n");
  printf("       successive costs must drop below to stop the optimization.  \n");
  printf("  --tol1d tol1d : tolerance on powell 1d minimizations\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/*------------------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/*------------------------------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/*------------------------------------------------------------------*/
static void check_options(void) {
  return;
}
/*------------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  return;
}
/*-----------------------------------------------------*/
MATRIX *MRIangles2RotMatB(double *angles, MATRIX *R)
{
  double gamma, beta, alpha;
  int r,c;
  MATRIX *R3, *Rx, *Ry, *Rz;

  gamma = angles[0];
  beta  = angles[1];
  alpha = angles[2];

  //printf("angles %g %g %g\n",angles[0],angles[1],angles[2]);

  Rx = MatrixZero(3,3,NULL);
  Rx->rptr[1][1] = +1;
  Rx->rptr[2][2] = +cos(gamma);
  Rx->rptr[2][3] = -sin(gamma);
  Rx->rptr[3][2] = +sin(gamma);
  Rx->rptr[3][3] = +cos(gamma);
  //printf("Rx ----------------\n");
  //MatrixPrint(stdout,Rx);

  Ry = MatrixZero(3,3,NULL);
  Ry->rptr[1][1] = +cos(beta);
  Ry->rptr[1][3] = +sin(beta);
  Ry->rptr[2][2] = 1;
  Ry->rptr[3][1] = -sin(beta);
  Ry->rptr[3][3] = +cos(beta);
  //printf("Ry ----------------\n");
  //MatrixPrint(stdout,Ry);

  Rz = MatrixZero(3,3,NULL);
  Rz->rptr[1][1] = +cos(alpha);
  Rz->rptr[1][2] = -sin(alpha);
  Rz->rptr[2][1] = +sin(alpha);
  Rz->rptr[2][2] = +cos(alpha);
  Rz->rptr[3][3] = +1;
  //printf("Rz ----------------\n");
  //MatrixPrint(stdout,Rz);

  // This will be a 3x3 matrix
  R3 = MatrixMultiply(Rz,Ry,NULL);
  R3 = MatrixMultiply(R3,Rx,R3);

  // Stuff 3x3 into a 4x4 matrix, with (4,4) = 1
  R = MatrixZero(4,4,R);
  for(c=1; c <= 3; c++){
    for(r=1; r <= 3; r++){
      R->rptr[r][c] = R3->rptr[r][c];
    }
  }
  R->rptr[4][4] = 1;

  MatrixFree(&Rx);
  MatrixFree(&Ry);
  MatrixFree(&Rz);
  MatrixFree(&R3);

  //printf("R ----------------\n");
  //MatrixPrint(stdout,R);

  return(R);
}
