/**
 * @brief Manages the computation of the gray/white stats used to
 * place the surface (currently in MRIScomputeBorderValues())
 * 
 */
/*
 * Original Author: Douglas N Greve (but basically a rewrite of mris_make_surfaces by BF)
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <float.h>
#include <errno.h>

#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "surfgrad.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "annotation.h"
#include "cmdargs.h"
#include "romp_support.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

struct utsname uts;
char *cmdline, cwd[2000];
int debug = 0, checkoptsonly = 0;
int nthreads = 1;

int main(int argc, char *argv[]) ;

const char *Progname = "mris_autodet_gwstats";
AutoDetGWStats adgws;
char *outfile = NULL;
char *involpath=NULL;
char *wmvolpath=NULL;
char *insurfpath=NULL;
char *lhsurfpath=NULL;
char *rhsurfpath=NULL;
const char *subject=NULL,*insurfname = "orig",*involname="brain.finalsurfs.mgz", *wmvolname="wm.mgz";
char tmpstr[2000];
char *SUBJECTS_DIR=NULL;
int hemicode = -1;

/*--------------------------------------------------*/
int main(int argc, char **argv) 
{
  int nargs;
  char *cmdline2, cwd[2000];

  nargs = handleVersionOption(argc, argv, "mris_autodet_gwstats");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  cmdline2 = argv2cmdline(argc,argv);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  Gdiag |= DIAG_SHOW ;

  if(argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if(checkoptsonly) return(0);

  // print out version of this program 
  printf("%s\n",getVersion().c_str());
  printf("\n");
  printf("cd %s\n",cwd);
  printf("setenv SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  printf("%s\n",cmdline2);
  printf("\n");
  fflush(stdout);

  printf("Reading in intensity volume %s\n",involpath);
  adgws.mri_T1 = MRIread(involpath);
  if(adgws.mri_T1==NULL) exit(1);

  printf("Reading in wm volume %s\n",wmvolpath);
  adgws.mri_wm = MRIread(wmvolpath);
  if(adgws.mri_wm==NULL) exit(1);

  if(lhsurfpath){
    printf("Reading in lhsurf %s\n",lhsurfpath);
    adgws.mrisADlh = MRISread(lhsurfpath);
    if(adgws.mrisADlh == NULL) exit(1);
    adgws.hemicode = 1;
  }
  if(rhsurfpath){
    printf("Reading in rhsurf %s\n",rhsurfpath);
    adgws.mrisADrh = MRISread(rhsurfpath);
    if(adgws.mrisADrh == NULL) exit(1);
    adgws.hemicode = 2;
  }
  if(insurfpath){
    printf("Reading in surf %s\n",insurfpath);
    adgws.mrisAD = MRISread(insurfpath);
    if(adgws.mrisAD == NULL) exit(1);
    adgws.hemicode = hemicode;
  }

  int err = adgws.AutoDetectStats();
  if(err) exit(1);
  adgws.Write(outfile);

  printf("#VMPC# mris_autodet_gwstats VmPeak  %d\n",GetVmPeak());
  printf("mris_autodet_gwstats done\n");

  return(0);

}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
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
    else if(!strcasecmp(option, "--version")) print_version() ;
    else if(!strcasecmp(option, "--debug"))   debug = 1;
    else if(!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if(!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if(!strcmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--i")){
      if(nargc < 1) CMDargNErr(option,1);
      involpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--wm")){
      if(nargc < 1) CMDargNErr(option,1);
      wmvolpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--surf")){
      if(nargc < 1) CMDargNErr(option,1);
      insurfpath = pargv[0];
      char *surfbasename = fio_basename(pargv[0],NULL);
      hemicode = -1;
      if(strncmp(surfbasename,"lh",2)==0) hemicode = 1;
      if(strncmp(surfbasename,"rh",2)==0) hemicode = 2;
      nargsused = 1;
    }
    else if(!strcmp(option, "--surfs")){
      if(nargc < 2) CMDargNErr(option,2);
      lhsurfpath = pargv[0];
      rhsurfpath = pargv[1];
      nargsused = 2;
    }
    else if(!strcmp(option, "--lh-surf")){
      if(nargc < 1) CMDargNErr(option,1);
      lhsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--rh-surf")){
      if(nargc < 1) CMDargNErr(option,1);
      rhsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--sd")){
      if(nargc < 1) CMDargNErr(option,1);
      printf("using %s as SUBJECTS_DIR...\n", pargv[0]) ;
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    }
    else if(!strcmp(option, "--s")){
      if(nargc < 1) CMDargNErr(option,1);
      subject = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--min_border_white") || !strcmp(option, "-seg-wlo")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.min_border_white = atof(pargv[0]);
      adgws.min_border_white_set = 1;
      nargsused = 1;
    }
    else if(!strcmp(option, "--max_border_white")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.max_border_white = atof(pargv[0]);
      adgws.max_border_white_set = 1;
      nargsused = 1;
    }
    else if(!strcmp(option, "--min_gray_at_white_border")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.min_gray_at_white_border = atof(pargv[0]);
      adgws.min_gray_at_white_border_set = 1;
      nargsused = 1;
    }
    else if(!strcmp(option, "--max_gray") || !strcmp(option, "-seg-ghi")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.max_gray = atof(pargv[0]);
      adgws.max_gray_set = 1;
      nargsused = 1;
    }
    else if(!strcmp(option, "--max_gray_at_csf_border")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.max_gray_at_csf_border = atof(pargv[0]);
      adgws.max_gray_at_csf_border_set = 1;
      nargsused = 1;
    }
    else if(!strcmp(option, "--min_gray_at_csf_border")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.min_gray_at_csf_border = atof(pargv[0]);
      adgws.min_gray_at_csf_border_set = 1;
      nargsused = 1;
    }
    else if(!strcmp(option, "--max_csf")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.max_csf = atof(pargv[0]);
      adgws.max_csf_set = 1;
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
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
/* --------------------------------------------- */
static void check_options(void) {
  if(outfile == NULL){
    printf("ERROR: no output set\n");
    exit(1);
  }
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(subject != NULL){
    if(involpath == NULL){
      sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,involname);
      involpath = strcpyalloc(tmpstr);
    }
    if(wmvolpath == NULL){
      sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,wmvolname);
      wmvolpath = strcpyalloc(tmpstr);
    }
    if(lhsurfpath == NULL){
      sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,insurfname);
      lhsurfpath = strcpyalloc(tmpstr);
    }
    if(rhsurfpath == NULL){
      sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,insurfname);
      rhsurfpath = strcpyalloc(tmpstr);
    }
  }
  if(involpath==NULL){
    printf("ERROR: no input volume set\n");
    exit(1);
  }
  if(wmvolpath==NULL){
    printf("ERROR: no wm volume set\n");
    exit(1);
  }
  if(lhsurfpath==NULL && rhsurfpath==NULL && insurfpath==NULL){
    printf("ERROR: no surface set\n");
    exit(1);
  }
  if(lhsurfpath!=NULL && insurfpath!=NULL){
    printf("ERROR: cannot spec lhsurf and insurf\n");
    exit(1);
  }
  if(rhsurfpath!=NULL && insurfpath!=NULL){
    printf("ERROR: cannot spec rhsurf and insurf\n");
    exit(1);
  }

  return;
}

/* --------------------------------------------- */
static void print_usage(void) 
{
  printf("\n");
  printf("PROGRAM: mris_autodet_gwstats\n");
  printf("Manages the computation of the gray/white statistics used to place the \n");
  printf("  white and pial surfaces (currently in MRIScomputeBorderValues()\n");
  printf(" --o outputfile : output text file with stats\n");
  printf(" --i  T1wvolume (usually brain.finalsurfs.mgz)\n");
  printf(" --wm wmvolume  (usually wm.mgz)\n");
  printf(" --surf  surf (usually ?h.orig)\n");
  printf(" --surfs lhsurf rhsurf \n");
  printf(" --lh-surf lhsurf \n");
  printf(" --rh-surf rhsurf \n");
  printf(" --s subject : reads in brain.finalsurfs.mgz, wm.mgz, lh.orig and rh.orig\n");
  printf(" --sd SUBJECTS_DIR \n");
  printf(" --min_border_white MinBW : Min border white  (can also use -seg-wlo for compatibility)\n");
  printf(" --max_border_white MaxBW : Max border white\n");
  printf(" --min_gray_at_white_border MinGWB\n");
  printf(" --max_gray MaxG (can also use -seg-ghi for compatibility)\n");
  printf(" --max_gray_at_csf_border MaxGCSFB\n");
  printf(" --min_gray_at_csf_border MinGCSFB\n");
  printf(" --max_csf MaxCSF\n");
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n");
  printf("\n");
  exit(1) ;
}

/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  return;
}

