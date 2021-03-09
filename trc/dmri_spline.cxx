/**
 * @brief Interpolate a spline from its control points
 *
 * Interpolate a spline from its control points
 */
/*
 * Original Author: Anastasia Yendiki
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
static void print_usage(void);
static void usage_exit(void);
static void print_help(void);
static void print_version(void);
static void dump_options();

int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]);

const char *Progname = "dmri_spline";

bool showControls = false;
std::string inFile, outVolFile, outVecBase, maskFile, outTextFile;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;

  nargs = handleVersionOption(argc, argv, "dmri_spline");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd, 2000);

  Progname = argv[0];
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL);
  DiagInit(NULL, NULL, NULL);

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);

  dump_options();

  Spline myspline(inFile.c_str(), maskFile.c_str());

  printf("Computing spline...\n");
  cputimer.reset();

  myspline.InterpolateSpline();

  cputime = cputimer.milliseconds();
  printf("Done in %g sec.\n", cputime/1000.0);

  if (!outVolFile.empty() ) {
    myspline.WriteVolume(outVolFile.c_str(), showControls);
  }

  if (!outTextFile.empty()) {
    myspline.WriteAllPoints(outTextFile.c_str());
  }

  if (!outVecBase.empty()) {
    std::string fname;

    printf("Computing analytical tangent, normal, and curvature...\n");
    cputimer.reset();

    myspline.ComputeTangent(true);
    myspline.ComputeNormal(true);
    myspline.ComputeCurvature(true);

    cputime = cputimer.milliseconds();
    printf("Done in %g sec.\n", cputime/1000.0);

    // Write tangent, normal, and curvature (analytical) to text files
    fname = outVecBase + "_tang.txt";
    myspline.WriteTangent(fname.c_str());
    fname = outVecBase + "_norm.txt";
    myspline.WriteNormal(fname.c_str());
    fname = outVecBase + "_curv.txt";
    myspline.WriteCurvature(fname.c_str());

    printf("Computing discrete tangent, normal, and curvature...\n");
    cputimer.reset();

    myspline.ComputeTangent(false);
    myspline.ComputeNormal(false);
    myspline.ComputeCurvature(false);

    cputime = cputimer.milliseconds();
    printf("Done in %g sec.\n", cputime/1000.0);

    // Write tangent, normal, and curvature (discrete) to text files
    fname = outVecBase + "_tang_diff.txt";
    myspline.WriteTangent(fname.c_str());
    fname = outVecBase + "_norm_diff.txt";
    myspline.WriteNormal(fname.c_str());
    fname = outVecBase + "_curv_diff.txt";
    myspline.WriteCurvature(fname.c_str());
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

    if (!strcasecmp(option, "--help"))  print_help();
    else if (!strcasecmp(option, "--version")) print_version();
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
    else if (!strcmp(option, "--show"))
      showControls = true;
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
static void print_usage(void) {
  cout
  << endl << "USAGE: " << Progname << endl << endl
  << "Basic inputs" << endl
  << "   --cpts <file>:" << endl
  << "     Input text file containing control points" << endl
  << "   --mask <file>:" << endl
  << "     Input mask volume (spline is not allowed to stray off mask)" << endl
  << endl 
  << "Outputs (at least one output type must be specified)" << endl
  << "   --out <file>:" << endl
  << "     Output volume of the interpolated spline" << endl
  << "   --show:" << endl
  << "     Highlight control points in output volume (default: no)" << endl
  << "   --outpts <file>:" << endl
  << "     Output text file containing all interpolated spline points" << endl
  << "   --outvec <base>:" << endl
  << "     Base name of output text files containing tangent vectors," << endl
  << "     normal vectors, and curvatures at every point along the" << endl
  << "     spline (both analytical and finite-difference versions)" << endl
  << endl
  << "Other options" << endl
  << "   --debug:     turn on debugging" << endl
  << "   --checkopts: don't run anything, just check options and exit" << endl
  << "   --help:      print out information on how to use this program" << endl
  << "   --version:   print out version and exit" << endl
  << endl;
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage();

  cout << endl
       << "..." << endl
       << endl;

  exit(1);
}

/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage();
  exit(1);
}

/* --------------------------------------------- */
static void print_version(void) {
  cout << getVersion() << endl;
  exit(1);
}

/* --------------------------------------------- */
static void check_options(void) {
  if(inFile.size() == 0) {
    cout << "ERROR: Must specify input text file" << endl;
    exit(1);
  }
  if(maskFile.size() == 0) {
    cout << "ERROR: Must specify mask volume" << endl;
    exit(1);
  }
  if((outVolFile.size() + outTextFile.size() + outVecBase.size()) == 0) {
    cout << "ERROR: Must specify at least one type of output file" << endl;
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options() {
  cout << endl
       << getVersion() << endl
       << "cwd " << cwd << endl
       << "cmdline " << cmdline << endl
       << "sysname  " << uts.sysname << endl
       << "hostname " << uts.nodename << endl
       << "machine  " << uts.machine << endl
       << "user     " << VERuser() << endl;

  cout << "Control points: " << inFile << endl;
  cout << "Mask volume: " << maskFile << endl;
  if (outVolFile.size() != 0) {
    cout << "Output volume: " << outVolFile << endl
         << "Show controls: " << showControls << endl;
  }
  if (outTextFile.size() != 0) {
    cout << "Output text file: " << outTextFile << endl;
  }
  if (outVecBase.size() != 0 ) {
    cout << "Output tangent vector file base name: " << outVecBase << endl;
  }

  return;
}

