// Random field simulator
//
// Things to do:
// 1. how to do two-sided tests --one-tailed --two-tailed
// 2. surf
// 3. smooth
// 4. cluster
// 5. glm
// 6. smooth with ubermask
// 7. label as ubermask
// 8. binarize mask
// 9. merge masks
// 10. invert mask
// 11. permute?
// 12. power?

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
#include "pdf.h"
#include "fsgdf.h"
#include "timer.h"
#include "matfile.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "randomfields.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_mcsim.c,v 1.4 2006/05/11 21:56:17 nicks Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

MRI  *TempVol=NULL;
char *TempVolFile=NULL;

char *subject=NULL, *hemi=NULL, *SUBJECTS_DIR=NULL;
MRIS *surf=NULL;

double VWPThresh = -1; // voxel-wise p-threshold
char  *VWSign = "abs";

char  *UberMaskFile=NULL;
MRI   *UberMask=NULL;

MRI   *Mask=NULL;
char  *MaskFile=NULL;
double MaskThresh=0;
char  *MaskSign = "abs";
int    MaskInvertFlag=0;

double fwhm =-1, gstd;
int nreps=0;

char *OutFile = NULL;

RFS *rfs = NULL;
int SynthSeed = 0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int nargs, nthrep;
  double gmean, gstddev, gmax;
  MRI *mritmp;

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

  if(SynthSeed > 0) RFspecSetSeed(rfs,SynthSeed);

  TempVol = MRIread(TempVolFile);
  if(TempVol == NULL){
    printf("ERROR: reading %s\n",TempVolFile);
    exit(1);
  }

  if(UberMaskFile){
    UberMask = MRIread(UberMaskFile);
    if(UberMask == NULL){
      printf("ERROR: reading %s\n",UberMaskFile);
      exit(1);
    }
  }

  if(TempVol->type != MRI_FLOAT){
    printf("INFO: changing template type to float\n");
    mritmp = MRIchangeType(TempVol,MRI_FLOAT,0,0,0);
    if(mritmp == NULL){
      printf("ERROR: could change type\n");
      exit(1);
    }
    MRIfree(&TempVol);
    TempVol = mritmp;
  }

  for(nthrep = 0; nthrep < nreps; nthrep++){
    RFsynth(TempVol,rfs,UberMask);
    RFglobalStats(TempVol, UberMask, &gmean, &gstddev, &gmax);
    printf("%3d %lf %lf %lf\n",nthrep,gmean,gstddev,gmax);
    MRIwrite(TempVol,"rf.mgh");
  }

  return 0;
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;
  int n;

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
    else if (!strcasecmp(option, "--mask-invert")) MaskInvertFlag = 1;

    else if (!strcasecmp(option, "--temp-vol")){
      if(nargc < 1) CMDargNErr(option,1);
      TempVolFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--temp-surf")){
      if(nargc < 2) CMDargNErr(option,2);
      subject = pargv[0];
      hemi    = pargv[1];
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      if(SUBJECTS_DIR == NULL){
	printf("ERROR: SUBJECTS_DIR not defined in environment\n");
	exit(1);
      }
      nargsused = 2;
    }
    else if (!strcasecmp(option, "--vw-pthresh")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&VWPThresh);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--vw-sign")){
      if(nargc < 1) CMDargNErr(option,1);
      VWSign = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--ubermask")){
      if(nargc < 1) CMDargNErr(option,1);
      UberMaskFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--mask")){
      if(nargc < 1) CMDargNErr(option,1);
      MaskFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--mask-thresh")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&MaskThresh);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--mask-sign")){
      if(nargc < 1) CMDargNErr(option,1);
      MaskSign = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--fwhm")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&fwhm);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--nreps")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nreps);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--field")){
      if(nargc < 1) CMDargNErr(option,1);
      rfs = RFspecInit(0,NULL);
      rfs->name = strcpyalloc(pargv[0]);
      if(RFname2Code(rfs)==-1){
	printf("ERROR: field %s not supported\n",rfs->name);
	exit(1);
      }
      printf("code = %d\n",rfs->code);
      RFnparams(rfs);
      if(nargc < rfs->nparams + 1){
	printf("ERROR: field %s requires %d parameters\n",rfs->name,rfs->nparams);
	exit(1);
      }
      for(n=0; n < rfs->nparams; n++)
	sscanf(pargv[n+1],"%lf",&(rfs->params[n]));
      RFexpectedMeanStddev(rfs);
      nargsused = rfs->nparams+1;
    }
    else if (!strcasecmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      OutFile = pargv[0];
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
  printf("   --temp-vol volfile : template volume \n");
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
  if(TempVolFile == NULL && subject == NULL){
    printf("ERROR: need template volume or surface\n");
    exit(1);
  }
  if(VWPThresh < 0){
    printf("ERROR: --vw-pthresh needs to be set\n");
    exit(1);
  }
  if(rfs==NULL){
    printf("ERROR: need to specify a random field\n");
    exit(1);
  }
  if(OutFile==NULL){
    printf("ERROR: need to specify an output file\n");
    exit(1);
  }
  if(nreps < 1){
    printf("ERROR: need to specify number of simulation repitions\n");
    exit(1);
  }
  if(fwhm < 0){
    printf("ERROR: need to specify FWHM\n");
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

  if(TempVolFile != NULL)
    fprintf(fp,"tempvol     %s\n",TempVolFile);
  if(subject != NULL)
    fprintf(fp,"subject     %s\n",subject);
  if(hemi != NULL)
    fprintf(fp,"hemi     %s\n",hemi);
  if(UberMaskFile != NULL)
    fprintf(fp,"ubermask     %s\n",UberMaskFile);
  if(Mask != NULL){
    fprintf(fp,"mask     %s\n",MaskFile);
    fprintf(fp,"maskthresh  %lf\n",MaskThresh);
    fprintf(fp,"masksign    %s\n",MaskSign);
    fprintf(fp,"maskinvert  %d\n",MaskInvertFlag);    
  }
  fprintf(fp,"vwpthresh  %lf\n",VWPThresh);
  fprintf(fp,"vwsign     %s\n",VWSign);
  fprintf(fp,"fwhm     %lf\n",fwhm);
  fprintf(fp,"nreps    %d\n",nreps);
  fprintf(fp,"OutFile  %s\n",OutFile);
  RFprint(stdout, rfs);

  return;
}
#if 0
/*-------------------------------------------------------------------*/
MRI *RFsynth(MRI *temp, RFS *rfs, MRI *rf, int nframes)
{

}
#endif
