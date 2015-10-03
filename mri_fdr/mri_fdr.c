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
 *    $Date: 2015/09/30 21:48:49 $
 *    $Revision: 1.1 $
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
#include "mri_conform.h"

double round(double x);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_fdr.c,v 1.1 2015/09/30 21:48:49 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

typedef struct {
  char *input;
  char *mask;
  double fdr;
  char *output;
  int log10flag, signid, frame;
  char *vthfile;
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

FSENV *fsenv;
MRI *input, *mask=NULL, *output;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err;
  struct timeb timer;
  double vthresh;

  TimerStart(&timer);

  cmdargs = (CMDARGS *)calloc(sizeof(CMDARGS),1);
  cmdargs->input = NULL;
  cmdargs->mask = NULL;
  cmdargs->output = NULL;
  cmdargs->fdr = -1;
  cmdargs->log10flag = 1;
  cmdargs->signid = 0;
  cmdargs->frame = 0;

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

  input = MRIread(cmdargs->input);
  if(!input) exit(1);

  if(cmdargs->mask){
    mask = MRIread(cmdargs->mask);
    if(!mask) exit(1);
  }

  output = NULL;
  if(cmdargs->output){
    output = MRIallocSequence(input->width,input->height,
			      input->depth,MRI_FLOAT,1);
    if(output == NULL) exit(1);
    MRIcopyHeader(input, output); /* does not change dimensions */
    MRIcopyPulseParameters(input, output);
  }

  err = MRIfdr2vwth(input, cmdargs->frame, cmdargs->fdr, cmdargs->signid,
		    cmdargs->log10flag, mask, &vthresh,output);
  printf("\nvoxel-wise-threshold %20.10lf\n",vthresh);
  printf("\n\n");

  if(cmdargs->vthfile){
    FILE *fp;
    fp = fopen(cmdargs->vthfile,"w");
    fprintf(fp,"%20.10lf\n",vthresh);    
    fclose(fp);
  }

  if(cmdargs->output){
    printf("Writing to %s\n",cmdargs->output);
    err = MRIwrite(output,cmdargs->output);
    if(err) exit(1);
  }

  printf("mri_fdr done\n");
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
    else if (!strcasecmp(option, "--pos")) cmdargs->signid = +1;
    else if (!strcasecmp(option, "--neg")) cmdargs->signid = -1;
    else if (!strcasecmp(option, "--abs")) cmdargs->signid =  0;
    else if (!strcasecmp(option, "--no-log10p")) cmdargs->log10flag = 0;
    else if (!strcasecmp(option, "--i")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->input = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--m")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->mask = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--o")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->output = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--thfile")) {
      if(nargc < 1) CMDargNErr(option,1);
      cmdargs->vthfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--fdr")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cmdargs->fdr);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--f")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&cmdargs->frame);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diag")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
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
/* -------------------------------------------------------- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -------------------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --i input : source volume or surface overlay\n");
  printf("   --fdr FDR : value between 0 and 1, typically .05\n");
  printf("\n");
  printf("   --m mask : binary mask\n");
  printf("   --o output : thresholded volume or surface overlay\n");
  printf("   --f frame : use frame as input (default is 0)\n");
  printf("   --pos : only consider positive voxels\n");
  printf("   --neg : only consider negative voxels\n");
  printf("   --abs : consider all voxels regardless of sign (default)\n");
  printf("   --no-log10p : input is raw p-values, not -log10(p)\n");
  printf("   --thfile txtfile : write threshold to text file\n");
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
  printf("This is a program that performs False Discovery Rate\n");
  printf("correction.\n");
  exit(1) ;
}
/* -------------------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* -------------------------------------------------------- */
static void check_options(void) {
  if(cmdargs->input == NULL){
    printf("ERROR: must spec --i\n");
    exit(1);
  }
  if(cmdargs->fdr < 0){
    printf("ERROR: must spec --fdr\n");
    exit(1);
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
  fprintf(fp,"fdr     %lf\n",cmdargs->fdr);
  return;
}
