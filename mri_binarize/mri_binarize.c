// $Id: mri_binarize.c,v 1.1 2006/05/31 20:39:53 greve Exp $

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
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "matrix.h"
#include "annotation.h"
#include "fmriutils.h"
#include "cmdargs.h"
#include "fsglm.h"
#include "gsl/gsl_cdf.h"
#include "pdf.h"
#include "fsgdf.h"
#include "timer.h"
#include "matfile.h"
#include "volcluster.h"
#include "surfcluster.h"


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_binarize.c,v 1.1 2006/05/31 20:39:53 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *InVolFile=NULL;
char *OutVolFile=NULL;
char *MergeVolFile=NULL;
double MinThresh, MaxThresh;
int MinThreshSet=0, MaxThreshSet=0;
int BinVal=1;
int BinValNot=0;

MRI *InVol,*OutVol,*MergeVol;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int nargs, c, r, s, nhits;
  double val,outputval;

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
  if(argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if(checkoptsonly) return(0);
  dump_options(stdout);

  InVol = MRIread(InVolFile);
  if(InVol==NULL) exit(1);

  if(MergeVolFile){
    MergeVol = MRIread(MergeVolFile);
    if(MergeVol==NULL) exit(1);
    if(MergeVol->width != InVol->width){
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
    if(MergeVol->height != InVol->height){
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
    if(MergeVol->depth != InVol->depth){
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
  }

  OutVol = MRIalloc(InVol->width,InVol->height,InVol->depth,MRI_INT);
  if(OutVol == NULL) exit(1);
  MRIcopyHeader(InVol, OutVol);

  nhits = 0;
  for(c=0; c < InVol->width; c++){
    for(r=0; r < InVol->height; r++){
      for(s=0; s < InVol->depth; s++){
	val = MRIgetVoxVal(InVol,c,r,s,0);
	if((MinThreshSet && (val < MinThresh)) || 
	   (MaxThreshSet && (val > MinThresh)) ){
	  // Not in the Range
	  if(MergeVol)
	    outputval = MRIgetVoxVal(MergeVol,c,r,s,0);
	  else
	    outputval = BinValNot;
	}
	else {
	  outputval = BinVal;
	  nhits ++;
	}
	MRIsetVoxVal(OutVol,c,r,s,0,outputval);
      }
    }
  }
  printf("Found %d values in range\n",nhits);
  MRIwrite(OutVol,OutVolFile);

  

  printf("mri_binarize done\n");
  exit(0);
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if(argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--i")){
      if(nargc < 1) CMDargNErr(option,1);
      InVolFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      OutVolFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--m")){
      if(nargc < 1) CMDargNErr(option,1);
      MergeVolFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--min")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&MinThresh);
      MinThreshSet = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--max")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&MaxThresh);
      MaxThreshSet = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--binval")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&BinVal);
      nargsused = 1;
    }
    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(CMDsingleDash(option))
	fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --i invol  : input volume \n");
  printf("   --min min  : min thresh (def is -inf)\n");
  printf("   --max max  : max thresh (def is +inf)\n");
  printf("   --o outvol : output volume \n");
  printf("   \n");
  printf("   --binval val : set vox within thresh to val (default is 1) \n");
  printf("   --binvalnot notval : set vox outside range to notval (default is 0) \n");
  printf("   --m mergevol : merge with mergevolume \n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(InVolFile == NULL){
    printf("ERROR: must specify input volume\n");
    exit(1);
  }
  if(OutVolFile == NULL){
    printf("ERROR: must specify output volume\n");
    exit(1);
  }
  if(MinThreshSet == 0 && MaxThreshSet == 0){
    printf("ERROR: must specify minimum and/or maximum threshold\n");
    exit(1);
  }
  if(MaxThreshSet && MinThreshSet && MaxThresh < MinThresh){
    printf("ERROR: max thresh = %g < min thresh = %g\n",
	   MaxThresh,MinThresh);
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"\n");
  fprintf(fp,"input      %s\n",InVolFile);
  fprintf(fp,"output     %s\n",OutVolFile);
  fprintf(fp,"min        %g\n",MinThresh);
  if(MinThreshSet) 
    fprintf(fp,"min        %g\n",MaxThresh);
  else
    fprintf(fp,"min        -infinity\n");
  if(MaxThreshSet) 
    fprintf(fp,"max        %g\n",MaxThresh);
  else
    fprintf(fp,"max        +infinity\n");
  fprintf(fp,"binval        %d\n",BinVal);
  fprintf(fp,"binvalnot     %d\n",BinValNot);
  if(MergeVolFile) 
    fprintf(fp,"merge      %s\n",MergeVolFile);
  return;
}
