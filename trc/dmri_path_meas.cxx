/**
 * @file  dmri_path_meas.c
 * @brief Write diffusion measures along path to text file
 *
 * Write diffusion measures along path to text file
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2010/12/16 06:56:46 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2010
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <float.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "error.h"
#include "diag.h"
#include "mri.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "timer.h"

#include "spline.h"

using namespace std;

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]) ;

static char vcid[] = "";
const char *Progname = "dmri_path_meas";

char *inFile = NULL, *outFile = NULL, *maskFile = NULL, *dtBase = NULL,
     fname[PATH_MAX];
MRI *l1, *l2, *l3, *v1, *meas[5];

struct utsname uts;
char *cmdline, cwd[2000];

struct timeb cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd, 2000);

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

  printf("Reading inputs...\n");
  sprintf(fname, "%s_L1.nii.gz", dtBase);
  l1 = MRIread(fname);
  sprintf(fname, "%s_L2.nii.gz", dtBase);
  l2 = MRIread(fname);
  sprintf(fname, "%s_L3.nii.gz", dtBase);
  l3 = MRIread(fname);
  sprintf(fname, "%s_V1.nii.gz", dtBase);
  v1 = MRIread(fname);

  Spline myspline(inFile, maskFile);

  printf("Computing spline...\n");
  TimerStart(&cputimer);
  myspline.InterpolateSpline();

  printf("Computing diffusion measures...\n");

  sprintf(fname, "%s_L1.nii.gz", dtBase);
  meas[0] = MRIread(fname);			// Axial diffusivity
  sprintf(fname, "%s_L2.nii.gz", dtBase);
  meas[1] = MRIread(fname);
  MRIadd(l3, meas[1], meas[1]);
  MRIscalarMul(meas[1], meas[1], .5);		// Radial diffusivity
  sprintf(fname, "%s_MD.nii.gz", dtBase);
  meas[2] = MRIread(fname);			// Mean diffusivity
  sprintf(fname, "%s_FA.nii.gz", dtBase);
  meas[3] = MRIread(fname);			// Fractional anisotropy

  cputime = TimerStop(&cputimer);
  printf("Done in %g sec.\n", cputime/1000.0);

  myspline.WriteValues(meas, 4, outFile);

  printf("dmri_path_meas done\n");
  return(0);
  exit(0);
}

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc, nargsused;
  char **pargv, *option;

  if (argc < 1) usage_exit();

  nargc = argc;
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
    else if (!strcmp(option, "--cpts")) {
      if (nargc < 1) CMDargNErr(option,1);
      inFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--out")) {
      if (nargc < 1) CMDargNErr(option,1);
      outFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--dtbase")) {
      if (nargc < 1) CMDargNErr(option,1);
      dtBase = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}

/* --------------------------------------------- */
static void print_usage(void) 
{
  printf("\n");
  printf("USAGE: ./dmri_path_meas\n");
  printf("\n");
  printf("   --cpts <filename>:   input text file containing control points\n");
  printf("   --mask <filename>:   input mask volume\n");
  printf("   --dtbase <filename>: base name of input dtifit files\n");
  printf("   --out <filename>:    output text file\n");
  printf("\n");
  printf("\n");
  printf("   --debug:     turn on debugging\n");
  printf("   --checkopts: don't run anything, just check options and exit\n");
  printf("   --help:      print out information on how to use this program\n");
  printf("   --version:   print out version and exit\n");
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n");
  printf("...\n");
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
  printf("%s\n", vcid) ;
  exit(1) ;
}

/* --------------------------------------------- */
static void check_options(void) {
  if(!inFile) {
    printf("ERROR: must specify input text file\n");
    exit(1);
  }
  if(!maskFile) {
    printf("ERROR: must specify mask volume\n");
    exit(1);
  }
  if(!outFile) {
    printf("ERROR: must specify output file\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  fprintf(fp, "Control points: %s\n", inFile);
  fprintf(fp, "Mask volume: %s\n", maskFile);
  fprintf(fp, "Input DTIFIT base: %s\n", dtBase);
  fprintf(fp, "Output volume: %s\n", outFile);

  return;
}

