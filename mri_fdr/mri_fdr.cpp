/**
 * @brief Computes registration between two volumes
 *
 */
/*
 * Original Author: Douglas N. Greve
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
#include "romp_support.h"
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

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

typedef struct {
  int ninputs;
  char *inputlist[100];
  char *masklist[100];
  char *outputlist[100];
  int  defaultframe, framelist[100];
  int log10flag, signid;
  double fdr;
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
MRI *inputlist[100], *masklist[100], *outputlist[100];

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err,n;
  Timer timer;
  double vthresh;

  timer.reset();

  cmdargs = (CMDARGS *)calloc(sizeof(CMDARGS),1);
  cmdargs->ninputs = 0;
  for(n=0; n < 100; n++){
    cmdargs->inputlist[n] = NULL;
    cmdargs->framelist[n] = -1;
    cmdargs->masklist[n] = NULL;
    cmdargs->outputlist[n] = NULL;
  }
  cmdargs->defaultframe = 0;
  cmdargs->fdr = -1;
  cmdargs->log10flag = 1;
  cmdargs->signid = 0;

  nargs = handleVersionOption(argc, argv, "mri_fdr");
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
  for(n=0; n < cmdargs->ninputs; n++){
    if(cmdargs->framelist[n] < 0) cmdargs->framelist[n] = cmdargs->defaultframe;
    inputlist[n] = MRIread(cmdargs->inputlist[n]);
    if(!inputlist[n]) exit(1);
    if(cmdargs->masklist[n]){
      masklist[n] = MRIread(cmdargs->masklist[n]);
      if(!masklist[n]) exit(1);
    }
    if(cmdargs->outputlist[n]){
      outputlist[n] = MRIallocSequence(inputlist[n]->width,inputlist[n]->height,
				       inputlist[n]->depth,MRI_FLOAT,1);
      if(outputlist[n] == NULL) exit(1);
      MRIcopyHeader(inputlist[n], outputlist[n]); /* does not change dimensions */
      MRIcopyPulseParameters(inputlist[n], outputlist[n]);
    }
  }

  err = MRIfdr2vwth(inputlist, cmdargs->ninputs, cmdargs->framelist, cmdargs->fdr, cmdargs->signid,
		    cmdargs->log10flag, masklist, &vthresh,outputlist);
  printf("\nvoxel-wise-threshold %20.10lf\n",vthresh);
  printf("\n\n");

  if(cmdargs->vthfile){
    FILE *fp;
    fp = fopen(cmdargs->vthfile,"w");
    fprintf(fp,"%20.10lf\n",vthresh);    
    fclose(fp);
  }

  for(n=0; n < cmdargs->ninputs; n++){
    if(cmdargs->outputlist[n]){
      printf("Writing to %s\n",cmdargs->outputlist[n]);
      err = MRIwrite(outputlist[n],cmdargs->outputlist[n]);
      if(err) exit(1);
    }
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
      int k = cmdargs->ninputs;
      cmdargs->inputlist[k] = pargv[0];
      nargsused = 1;
      if(CMDnthIsArg(nargc, pargv, 1)){
	if(strcmp(pargv[1],"nomask") != 0) cmdargs->masklist[k] = pargv[1];
        nargsused ++;
	if(CMDnthIsArg(nargc, pargv, 2)){
	  if(strcmp(pargv[2],"nooutput") != 0) cmdargs->outputlist[k] = pargv[2];
	  nargsused ++;
	  if(CMDnthIsArg(nargc, pargv, 3)){
	    sscanf(pargv[3],"%d",&cmdargs->framelist[k]);
	    nargsused ++;
	  } 
	} 
      } 
      cmdargs->ninputs++;
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
      sscanf(pargv[0],"%d",&cmdargs->defaultframe);
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
  printf("   --i input <mask> <output> <frame>: source volume or surface overlay\n");
  printf("     --i input <mask> <output> <frame>: another source volume or surface overlay\n");
  printf("   --fdr FDR : value between 0 and 1, typically .05\n");
  printf("\n");
  printf("   --f defaultframe : use input frame when not specing frame in --i\n");
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
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* -------------------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("This is a program that performs False Discovery Rate correction.\n");
  printf("Basic usage:\n");
  printf("  mri_fdr --fdr .05 --i sig.mgh mask.mgh output.mgh \n");
  printf("To perform FDR on both the left and right hemis simultaneously:\n");
  printf("  mri_fdr --fdr .05 --i lh.sig.mgh lh.mask.mgh lh.output.mgh \n");
  printf("                    --i rh.sig.mgh rh.mask.mgh rh.output.mgh \n");
  printf("Masks are optional\n");
  printf("Outputs are optional (in which case you probably just want the threshold)\n");
  printf("  You can capture the threshold in a text file with --thfile\n");
  printf("The input is assumed to be -log10(p). If that is not the case then --no-log10p\n");
  printf("Output threshold will be -log10(p) unless --no-log10p\n");
  printf("If you want to spec an output but not a mask,  then set the mask file to 'nomask'\n");
  printf("If you want to spec a frame but not an output, then set the output file to 'nooutput'\n");
  printf("\n");
  printf("  Thresholding of Statistical Maps in Functional Neuroimaging Using\n");
  printf("  the False Discovery Rate.  Christopher R. Genovese, Nicole A. Lazar,\n");
  printf("  Thomas E. Nichols (2002).  NeuroImage 15:870-878.\n");
  printf("\n");

  exit(1) ;
}
/* -------------------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* -------------------------------------------------------- */
static void check_options(void) {
  if(cmdargs->ninputs == 0){
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
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"fdr       %lf\n",cmdargs->fdr);
  fprintf(fp,"ninputs   %d\n",cmdargs->ninputs);
  return;
}
