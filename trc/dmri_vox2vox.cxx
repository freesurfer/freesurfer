/**
 * @brief Apply affine and non-linear warp to voxel coordinates in text file
 *
 * Apply affine and non-linear warp to voxel coordinates in text file
 */
/*
 * Original Author: Anastasia Yendiki
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

#include "vial.h"	// Needs to be included first because of CVS libs

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

const char *Progname = "dmri_vox2vox";

int doInvNonlin = 0;
unsigned int nFile = 0;
string inDir, outDir, inRefFile, outRefFile, affineXfmFile, nonlinXfmFile;
vector<string> inFile, outFile;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  string fname;
  vector<float> point(3);
  MRI *inref = 0, *outref = 0;
  AffineReg affinereg;
#ifndef NO_CVS_UP_IN_HERE
  NonlinReg nonlinreg;
#endif

  nargs = handleVersionOption(argc, argv, "dmri_vox2vox");
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

  // Read reference volumes
  inref = MRIread(inRefFile.c_str());
  outref = MRIread(outRefFile.c_str());

  // Read transform files
#ifndef NO_CVS_UP_IN_HERE
  if (!nonlinXfmFile.empty()) {
    if (!affineXfmFile.empty())
      affinereg.ReadXfm(affineXfmFile, inref, outref);
    nonlinreg.ReadXfm(nonlinXfmFile, outref);
  }
  else
#endif
  if (!affineXfmFile.empty())
    affinereg.ReadXfm(affineXfmFile, inref, outref);

  for (unsigned int ifile = 0; ifile < nFile; ifile++) {
    float coord;
    ifstream infile;
    ofstream outfile;
    vector<float> inpts;

    printf("Processing coordinate file %d of %d...\n", ifile+1, nFile);
    cputimer.reset();

    // Read input text file
    fname = inFile[ifile];

    if (!inDir.empty())
      fname = inDir + "/" + fname;

    infile.open(fname, ios::in);
    if (!infile) {
      cout << "ERROR: Could not open " << fname << " for reading" << endl;
      exit(1);
    }

    inpts.clear();
    while (infile >> coord)
      inpts.push_back(coord);

    if (inpts.size() % 3 != 0) {
      cout << "ERROR: File " << fname
           << " must contain triplets of coordinates" << endl;
      exit(1);
    }

    infile.close();

    for (vector<float>::iterator ipt = inpts.begin();
                                 ipt < inpts.end(); ipt += 3) {
      copy(ipt, ipt+3, point.begin());

      // Apply affine transform
      if (!affinereg.IsEmpty())
        affinereg.ApplyXfm(point, point.begin());

#ifndef NO_CVS_UP_IN_HERE
      // Apply nonlinear transform
      if (!nonlinreg.IsEmpty()) {
        if (doInvNonlin)
          nonlinreg.ApplyXfmInv(point, point.begin());
        else
          nonlinreg.ApplyXfm(point, point.begin());
      }
#endif

      copy(point.begin(), point.end(), ipt);
    }

    // Write output text file
    fname = outFile[ifile];

    if (!outDir.empty())
      fname = outDir + "/" + fname;

    outfile.open(fname, ios::out);
    if (!outfile) {
      cout << "ERROR: Could not open " << fname << " for writing" << endl;
      exit(1);
    }

    for (vector<float>::const_iterator ipt = inpts.begin();
                                       ipt < inpts.end(); ipt += 3)
        outfile << ipt[0] << " " << ipt[1] << " " << ipt[2] << endl;

    outfile.close();

    cputime = cputimer.milliseconds();
    printf("Done in %g sec.\n", cputime/1000.0);
  }

  MRIfree(&inref);
  MRIfree(&outref);

  printf("dmri_vox2vox done\n");
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
    else if (!strcmp(option, "--indir")) {
      if (nargc < 1) CMDargNErr(option,1);
      inDir = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--in")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        inFile.push_back(pargv[nargsused]);
        nargsused++;
      }
    } 
    else if (!strcmp(option, "--outdir")) {
      if (nargc < 1) CMDargNErr(option,1);
      outDir = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--out")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        outFile.push_back(pargv[nargsused]);
        nargsused++;
      }
    } 
    else if (!strcmp(option, "--inref")) {
      if (nargc < 1) CMDargNErr(option,1);
      inRefFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--outref")) {
      if (nargc < 1) CMDargNErr(option,1);
      outRefFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--reg")) {
      if (nargc < 1) CMDargNErr(option,1);
      affineXfmFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--regnl")) {
      if (nargc < 1) CMDargNErr(option,1);
      nonlinXfmFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--invnl"))
      doInvNonlin = 1;
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
  printf("USAGE: ./dmri_vox2vox\n");
  printf("\n");
  printf("Basic inputs\n");
  printf("   --in <file> [...]:\n");
  printf("     Input text file(s)\n");
  printf("   --indir <dir>:\n");
  printf("     Input directory (optional)\n");
  printf("     If specified, names of input text files are relative to this\n");
  printf("   --out <file> [...]:\n");
  printf("     Output text file(s), as many as inputs\n");
  printf("   --outdir <dir>:\n");
  printf("     Output directory (optional)\n");
  printf("     If specified, names of output text files are relative to this)\n");
  printf("   --inref <file>:\n");
  printf("     Input reference volume\n");
  printf("   --outref <file>:\n");
  printf("     Output reference volume\n");
  printf("   --reg <file>:\n");
  printf("     Affine registration (.mat), applied first\n");
  printf("   --regnl <file>:\n");
  printf("     Nonlinear registration (.m3z), applied second\n");
  printf("   --invnl:\n");
  printf("     Apply inverse of nonlinear warp (with --regnl, default: no)\n");
  printf("\n");
  printf("Other options\n");
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
  std::cout << getVersion() << std::endl;
  exit(1) ;
}

/* --------------------------------------------- */
static void check_options(void) {
  nFile = inFile.size();

  if(nFile == 0) {
    printf("ERROR: must specify input text file(s)\n");
    exit(1);
  }
  if(outFile.size() != nFile) {
    printf("ERROR: must specify as many output text files as input files\n");
    exit(1);
  }
  if(inRefFile.empty()) {
    printf("ERROR: must specify input reference volume\n");
    exit(1);
  }
  if(outRefFile.empty()) {
    printf("ERROR: must specify output reference volume\n");
    exit(1);
  }
  return;
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

  if (!inDir.empty())
    fprintf(fp, "Input directory: %s\n", inDir.c_str());
  fprintf(fp, "Input files:");
  for (unsigned int ifile = 0; ifile < nFile; ifile++)
    fprintf(fp, " %s", inFile[ifile].c_str());
  fprintf(fp, "\n");
  if (!outDir.empty())
    fprintf(fp, "Output directory: %s\n", outDir.c_str());
  fprintf(fp, "Output files:");
  for (unsigned int ifile = 0; ifile < nFile; ifile++)
    fprintf(fp, " %s", outFile[ifile].c_str());
  fprintf(fp, "\n");
  fprintf(fp, "Input reference: %s\n", inRefFile.c_str());
  fprintf(fp, "Output reference: %s\n", outRefFile.c_str());
  if (!affineXfmFile.empty())
    fprintf(fp, "Affine registration: %s\n", affineXfmFile.c_str());
  if (!nonlinXfmFile.empty()) {
    fprintf(fp, "Nonlinear registration: %s\n", nonlinXfmFile.c_str());
    fprintf(fp, "Invert nonlinear morph: %d\n", doInvNonlin);
  }

  return;
}

