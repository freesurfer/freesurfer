/**
 * @file  dmri_spline.c
 * @brief Interpolate a spline from its control points
 *
 * Interpolate a spline from its control points
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2012/07/18 16:44:56 $
 *    $Revision: 1.6 $
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
#include <limits.h>

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
const char *Progname = "dmri_spline";

char *inFile = NULL, *maskFile = NULL,
     *outVolFile = NULL, *outTextFile = NULL, *outVecBase = NULL;
//char *inFiles[50];

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

  Spline myspline(inFile, maskFile);

  printf("Computing spline...\n");
  TimerStart(&cputimer);

  myspline.InterpolateSpline();

  cputime = TimerStop(&cputimer);
  printf("Done in %g sec.\n", cputime/1000.0);

  if (outVolFile)
    myspline.WriteVolume(outVolFile);

  if (outTextFile)
    myspline.WriteAllPoints(outTextFile);

  if (outVecBase) {
    char fname[PATH_MAX];

    printf("Computing analytical tangent, normal, and curvature...\n");
    TimerStart(&cputimer);

    myspline.ComputeTangent(true);
    myspline.ComputeNormal(true);
    myspline.ComputeCurvature(true);

    cputime = TimerStop(&cputimer);
    printf("Done in %g sec.\n", cputime/1000.0);

    // Write tangent, normal, and curvature (analytical) to text files
    sprintf(fname, "%s_tang.txt", outVecBase);
    myspline.WriteTangent(fname);
    sprintf(fname, "%s_norm.txt", outVecBase);
    myspline.WriteNormal(fname);
    sprintf(fname, "%s_curv.txt", outVecBase);
    myspline.WriteCurvature(fname);

    printf("Computing discrete tangent, normal, and curvature...\n");
    TimerStart(&cputimer);

    myspline.ComputeTangent(false);
    myspline.ComputeNormal(false);
    myspline.ComputeCurvature(false);

    cputime = TimerStop(&cputimer);
    printf("Done in %g sec.\n", cputime/1000.0);

    // Write tangent, normal, and curvature (discrete) to text files
    sprintf(fname, "%s_tang_diff.txt", outVecBase);
    myspline.WriteTangent(fname);
    sprintf(fname, "%s_norm_diff.txt", outVecBase);
    myspline.WriteNormal(fname);
    sprintf(fname, "%s_curv_diff.txt", outVecBase);
    myspline.WriteCurvature(fname);
  }

  printf("dmri_spline done\n");
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
      outVolFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--outpts")) {
      if (nargc < 1) CMDargNErr(option,1);
      outTextFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--outvec")) {
      if (nargc < 1) CMDargNErr(option,1);
      outVecBase = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskFile = fio_fullpath(pargv[0]);
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
static void print_usage(void) 
{
  printf("\n");
  printf("USAGE: ./dmri_spline\n");
  printf("\n");
  printf("   --cpts   <file>: input text file containing control points\n");
  printf("   --mask   <file>: input mask volume\n");
  printf("   --out    <file>: output spline volume\n");
  printf("   --outpts <file>: output text file containing all spline points\n");
  printf("   --outvec <base>: base name of output text files containing the\n");
  printf("                    tangent vectors, normal vectors, and curvatures\n");
  printf("                    of the spline at every point along the spline\n");
  printf("                    (both analytical and finite-difference version)\n");
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
  if(!outVolFile && !outTextFile && !outVecBase) {
    printf("ERROR: must specify at least one type of output file\n");
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
  if (outVolFile)
    fprintf(fp, "Output volume: %s\n", outVolFile);
  if (outTextFile)
    fprintf(fp, "Output text file: %s\n", outTextFile);
  if (outVecBase)
    fprintf(fp, "Output tangent vector file base name: %s\n", outVecBase);

  return;
}

