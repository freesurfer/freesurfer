/**
 * @file  dmri_train.cxx
 * @brief Training of priors for probabilistic global tractography
 *
 * Training of priors for probabilistic global tractography
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2011/02/25 03:52:35 $
 *    $Revision: 1.2 $
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

#include "blood.h"

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
const char *Progname = "dmri_train";

bool useTrunc = false;
int nControl[100], nout = 0, ntrk = 0, nroi1 = 0, nroi2 = 0, nlab = 0, ncpt = 0;
float trainMaskLabel[100];
char *outDir = NULL, *outBase[100], *trainListFile = NULL,
     *trainTrkList[100], *trainRoi1List[100], *trainRoi2List[100],
     *trainAsegFile = NULL, *trainMaskFile = NULL, *testMaskFile = NULL;

struct utsname uts;
char *cmdline, cwd[2000];

struct timeb cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  char fbase[PATH_MAX];

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

  Blood myblood(trainListFile, trainTrkList[0],
                nroi1 ? trainRoi1List[0] : 0, nroi2 ? trainRoi2List[0] : 0,
                trainAsegFile, trainMaskFile, nlab ? trainMaskLabel[0] : 0,
                testMaskFile, useTrunc);

  for (int itrk = 0; itrk < ntrk; itrk++) {
    if (itrk > 0)
      myblood.ReadStreamlines(trainListFile, trainTrkList[itrk],
                              nroi1 ? trainRoi1List[itrk] : 0,
                              nroi2 ? trainRoi2List[itrk] : 0,
                              nlab ? trainMaskLabel[itrk] : 0);

    printf("Processing pathway %d of %d...\n", itrk+1, ntrk);
    TimerStart(&cputimer);

    myblood.ComputePriors();

    for (int icpt = 0; icpt < ncpt; icpt++)
      myblood.SelectControlPoints(nControl[icpt]);

    if (outDir)
      sprintf(fbase, "%s/%s", outDir, outBase[itrk]);
    else
      strcpy(fbase, outBase[itrk]);

    myblood.WriteOutputs(fbase);

    cputime = TimerStop(&cputimer);
    printf("Done in %g sec.\n", cputime/1000.0);
  }

  printf("dmri_train done\n");
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
    else if (!strcmp(option, "--outdir")) {
      if (nargc < 1) CMDargNErr(option,1);
      outDir = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--out")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (strncmp(pargv[nargsused], "--", 2)) {
        outBase[nout] = pargv[nargsused];
        nargsused++;
        nout++;
      }
    }
    else if (!strcmp(option, "--slist")) {
      if (nargc < 1) CMDargNErr(option,1);
      trainListFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--trk")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (strncmp(pargv[nargsused], "--", 2)) {
        trainTrkList[ntrk] = pargv[nargsused];
        nargsused++;
        ntrk++;
      }
    }
    else if (!strcmp(option, "--rois")) {
      if (nargc < 2) CMDargNErr(option,1);
      while (strncmp(pargv[nargsused], "--", 2)) {
        trainRoi1List[nroi1] = pargv[nargsused];
        nargsused++;
        nroi1++;
        trainRoi2List[nroi2] = pargv[nargsused];
        nargsused++;
        nroi2++;
      }
    }
    else if (!strcmp(option, "--seg")) {
      if (nargc < 1) CMDargNErr(option,1);
      trainAsegFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--cmask")) {
      if (nargc < 1) CMDargNErr(option,1);
      trainMaskFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--lmask")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (strncmp(pargv[nargsused], "--", 2)) {
        sscanf(pargv[nargsused],"%f",&trainMaskLabel[nlab]);
        nargsused++;
        nlab++;
      }
    }
    else if (!strcmp(option, "--bmask")) {
      if (nargc < 1) CMDargNErr(option,1);
      testMaskFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--ncpts")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (strncmp(pargv[nargsused], "--", 2)) {
        sscanf(pargv[nargsused],"%u",&nControl[ncpt]);
        nargsused++;
        ncpt++;
      }
    }
    else if (!strcmp(option, "--trunc"))
      useTrunc = true;

    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}

/* --------------------------------------------- */
static void print_usage(void) 
{
  printf("\n");
  printf("USAGE: ./dmri_train\n");
  printf("\n");
  printf("Basic inputs\n");
  printf("   --out <base> [...]:\n");
  printf("     Base name(s) of output(s), one per path\n");
  printf("   --outdir <dir>:\n");
  printf("     Output directory (optional)\n");
  printf("     If specified, base names of outputs are relative to this\n");
  printf("   --slist <file>:\n");
  printf("     Text file with list of training subject directories\n");
  printf("   --trk <file> [...]:\n");
  printf("     Name(s) of input .trk file(s), one per path\n");
  printf("     (Names relative to training subject directory)\n");
  printf("   --rois <file1> <file2> [...]:\n");
  printf("     Optional, names of input tract labeling ROIs, two per path\n");
  printf("     (Names relative to training subject directory)\n");
  printf("   --seg <file>:\n");
  printf("     Name of input aparc+aseg volume\n");
  printf("     (Name relative to training subject directory)\n");
  printf("   --cmask <file>:\n");
  printf("     Name of input cortex mask volume\n");
  printf("   --lmask <id> [...]:\n");
  printf("     Add a label ID from aparc+aseg to cortex mask, one per path\n");
  printf("     (0 doesn't add any label)\n");
  printf("   --bmask <file>:\n");
  printf("     Input brain mask volume for test subject\n");
  printf("   --ncpts <num> [...]:\n");
  printf("     Number of control points for initial spline\n");
  printf("   --trunc:\n");
  printf("     Also save results using all streamlines, truncated or not\n");
  printf("     (Default: Only save results using non-truncated streamlines)\n");
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
  printf("%s\n", vcid) ;
  exit(1) ;
}

/* --------------------------------------------- */
static void check_options(void) {
  if(nout == 0) {
    printf("ERROR: must specify at least one output name\n");
    exit(1);
  }
  if(!trainListFile) {
    printf("ERROR: must specify training subject list file\n");
    exit(1);
  }
  if(ntrk == 0) {
    printf("ERROR: must specify location of at least one streamline file\n");
    exit(1);
  }
  if(ntrk != nout) {
    printf("ERROR: numbers of input .trk files and output names must match\n");
    exit(1);
  }
  if(!trainAsegFile) {
    printf("ERROR: must specify location of aparc+aseg volume\n");
    exit(1);
  }
  if(!trainMaskFile) {
    printf("ERROR: must specify location of cortex mask volume\n");
    exit(1);
  }
  if(nroi1 > 0 && nroi1 != ntrk) {
    printf("ERROR: numbers of input .trk files and start ROIs must match\n");
    exit(1);
  }
  if(nroi2 > 0 && nroi2 != ntrk) {
    printf("ERROR: numbers of input .trk files and end ROIs must match\n");
    exit(1);
  }
  if(nlab > 0 && nlab != ntrk) {
    printf("ERROR: numbers of input .trk files and label IDs must match\n");
    exit(1);
  }
  if(!testMaskFile) {
    printf("ERROR: must specify brain mask volume for output subject\n");
    exit(1);
  }
  if(ncpt == 0) {
    printf("ERROR: must specify number of control points for initial spline\n");
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

  if (outDir)
    fprintf(fp, "Output directory: %s\n", outDir);
  fprintf(fp, "Output base:");
  for (int k = 0; k < nout; k++)
    fprintf(fp, " %s", outBase[k]);
  fprintf(fp, "\n");
  fprintf(fp, "Training subject directory list: %s\n", trainListFile);
  fprintf(fp, "Location of streamline files relative to base:");
  for (int k = 0; k < ntrk; k++)
    fprintf(fp, " %s", trainTrkList[k]);
  fprintf(fp, "\n");
  if (nroi1 > 0) {
    fprintf(fp, "Location of start ROI files relative to base:");
    for (int k = 0; k < nroi1; k++)
      fprintf(fp, " %s", trainRoi1List[k]);
    fprintf(fp, "\n");
    fprintf(fp, "Location of end ROI files relative to base:");
    for (int k = 0; k < nroi2; k++)
      fprintf(fp, " %s", trainRoi2List[k]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "Location of cortex masks relative to base: %s\n", trainMaskFile);
  if (nlab > 0) {
    fprintf(fp, "Label ID's from aparc+aseg to add to cortex mask:");
    for (int k = 0; k < nlab; k++)
      fprintf(fp, " %1.0f ", trainMaskLabel[k]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "Location of aparc+aseg's relative to base: %s\n", trainAsegFile);
  fprintf(fp, "Brain mask for output subject: %s\n", testMaskFile);
  fprintf(fp, "Number of control points for initial spline:");
  for (int k = 0; k < ncpt; k++)
    fprintf(fp, " %d ", nControl[k]);
  fprintf(fp, "\n");
  fprintf(fp, "Use truncated streamlines: %d\n", useTrunc);

  return;
}

