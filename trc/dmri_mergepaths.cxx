/**
 * @brief Merge posterior distributions from multiple paths into a 4D file
 *
 * Merge posterior distributions from multiple paths into a 4D file
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

#include "cma.h"
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

const char *Progname = "dmri_mergepaths";

int nframe = 0;
float dispThresh = 0;
string inDir, outFile, ctabFile;
vector<string> inFile;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  string fname, pname;
  MRI *invol = 0, *outvol = 0;

  nargs = handleVersionOption(argc, argv, "dmri_mergepaths");
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

  cputimer.reset();

  for (int iframe = 0; iframe < nframe; iframe++) {
    float inmax = 0;

    cout << "Merging volume " << iframe+1 << " of " << nframe << "... " << endl;

    // Read input volume
    fname = inFile[iframe];

    if (!inDir.empty())
      fname = inDir + "/" + fname;

    invol = MRIread(fname.c_str());

    // Copy input volume to output volume series
    if (invol) {
      if (!outvol) {
        // Allocate output 4D volume
        outvol = MRIcloneBySpace(invol, invol->type, nframe);

        // Read color table
        outvol->ct = CTABreadASCII(ctabFile.c_str());
      }

      MRIcopyFrame(invol, outvol, 0, iframe);

      if (dispThresh > 0)
        inmax = (float) MRIfindPercentile(invol, .99, 0);	// Robust max
    }

    // Set display threshold for current volume
    outvol->frames[iframe].thresh = dispThresh * inmax;

    // Look up pathway name for current volume
    pname = inFile[iframe];
    pname = pname.substr(0, pname.rfind("/"));
    if (pname.rfind("/") != string::npos)
      pname = pname.substr(pname.rfind("/")+1, string::npos);
    pname = pname.substr(0, pname.find("_"));

    outvol->frames[iframe].label = 0;

    for (int ict = outvol->ct->nentries; ict > 0; ict--) {
      CTE *cte = outvol->ct->entries[ict];

      if (cte != NULL && strstr(pname.c_str(), cte->name)
                      && strlen(pname.c_str()) == strlen(cte->name)) {
        outvol->frames[iframe].label = ict;
        strcpy(outvol->frames[iframe].name, cma_label_to_name(ict));

        continue;
      }
    }

    cout << "Threshold: " << outvol->frames[iframe].thresh
         << " Name: "     << outvol->frames[iframe].name << endl;
  }

  // Write output file
  if (outvol)
    MRIwrite(outvol, outFile.c_str());
  else {
    cout << "ERROR: could not open any of the input files" << endl;
    exit(1);
  }

  cputime = cputimer.milliseconds();
  printf("Done in %g sec.\n", cputime/1000.0);

  printf("dmri_mergepaths done\n");
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
    else if (!strcmp(option, "--out")) {
      if (nargc < 1) CMDargNErr(option,1);
      outFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--ctab")) {
      if (nargc < 1) CMDargNErr(option,1);
      ctabFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&dispThresh);
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
  printf("USAGE: ./dmri_mergepaths\n");
  printf("\n");
  printf("Basic inputs\n");
  printf("   --in <file> [...]:\n");
  printf("     Input volume(s)\n");
  printf("   --indir <dir>:\n");
  printf("     Input directory (optional)\n");
  printf("     If specified, names of input files are relative to this\n");
  printf("   --out <file>:\n");
  printf("     Output volume\n");
  printf("   --ctab <file>:\n");
  printf("     Color table file\n");
  printf("   --thresh <num>:\n");
  printf("     Lower threshold for display (0<=num<=1, as fraction of max)\n");
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
  nframe = inFile.size();

  if(nframe == 0) {
    printf("ERROR: must specify input volume(s)\n");
    exit(1);
  }
  if(outFile.empty()) {
    printf("ERROR: must specify output volume\n");
    exit(1);
  }
  if(ctabFile.empty()) {
    printf("ERROR: must specify color table file\n");
    exit(1);
  }
  if(dispThresh < 0 || dispThresh > 1) {
    printf("ERROR: display threshold must a number between 0 and 1\n");
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
  for (int k = 0; k < nframe; k++)
    fprintf(fp, " %s", inFile[k].c_str());
  fprintf(fp, "\n");
  fprintf(fp, "Output file: %s\n", outFile.c_str());
  fprintf(fp, "Color table file: %s\n", ctabFile.c_str());
  fprintf(fp, "Lower threshold for display: %f\n", dispThresh);

  return;
}

