/**
 * @file  dmri_track.c
 * @brief Probabilistic global tractography
 *
 * Probabilistic global tractography
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2010/07/23 22:34:16 $
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

#include "coffin.h"

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
const char *Progname = "dmri_track";

unsigned int nTract = 1, 
             nBurnIn = 5000, nSample = 5000, nKeepSample = 10, nUpdateProp = 40,
             asegPriorType = 0;
char *outDir = NULL, *dwiFile = NULL,
     *gradFile = NULL, *bvalFile = NULL,
     *maskFile = NULL, *bedpostDir = NULL,
     *roiFile1 = NULL, *roiFile2 = NULL,
     *roiMeshFile1 = NULL, *roiMeshFile2 = NULL,
     *roiRefFile1 = NULL, *roiRefFile2 = NULL,
     *xfmFile = NULL, *initFile = NULL,
     *priorFile0 = NULL, *priorFile1 = NULL,
     *asegPriorFile0 = NULL, *asegPriorFile1 = NULL, *asegIdFile = NULL,
     *asegTrainFile = NULL, *pathTrainFile = NULL,
     *asegFile = NULL,
     *stdPropFile = NULL;

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

  Coffin mycoffin(outDir, dwiFile,
                  gradFile, bvalFile,
                  maskFile,
                  bedpostDir, nTract,
                  roiFile1, roiFile2,
                  roiMeshFile1, roiMeshFile2,
                  roiRefFile1, roiRefFile2,
                  xfmFile, initFile,
                  priorFile0, priorFile1,
                  asegPriorType,
                  asegPriorFile0, asegPriorFile1, asegIdFile,
                  asegTrainFile, pathTrainFile,
                  asegFile,
                  nBurnIn, nSample,
                  nKeepSample, nUpdateProp,
                  stdPropFile,
                  debug);

  TimerStart(&cputimer);
  //mycoffin.RunMCMC();
mycoffin.RunMCMC1();
  cputime = TimerStop(&cputimer);
  printf("Done in %g sec.\n", cputime/1000.0);
  mycoffin.WriteOutputs();

  printf("dmri_track done\n");
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
    else if (!strcmp(option, "--dwi")) {
      if (nargc < 1) CMDargNErr(option,1);
      dwiFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--grad")) {
      if (nargc < 1) CMDargNErr(option,1);
      gradFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--bval")) {
      if (nargc < 1) CMDargNErr(option,1);
      bvalFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--bpdir")) {
      if (nargc < 1) CMDargNErr(option,1);
      bedpostDir = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--ntr")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%u",&nTract);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--roi1")) {
      if (nargc < 1) CMDargNErr(option,1);
      roiFile1 = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--roi2")) {
      if (nargc < 1) CMDargNErr(option,1);
      roiFile2 = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--roimesh1")) {
      if (nargc < 1) CMDargNErr(option,1);
      roiMeshFile1 = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--roimesh2")) {
      if (nargc < 1) CMDargNErr(option,1);
      roiMeshFile2 = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--roiref1")) {
      if (nargc < 1) CMDargNErr(option,1);
      roiRefFile1 = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--roiref2")) {
      if (nargc < 1) CMDargNErr(option,1);
      roiRefFile2 = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--dwi2roi")) {
      if (nargc < 1) CMDargNErr(option,1);
      xfmFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--init")) {
      if (nargc < 1) CMDargNErr(option,1);
      initFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--sdp")) {
      if (nargc < 1) CMDargNErr(option,1);
      stdPropFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--prior")) {
      if (nargc < 2) CMDargNErr(option,2);
      priorFile0 = fio_fullpath(pargv[0]);
      priorFile1 = fio_fullpath(pargv[1]);
      nargsused = 2;
    }
    else if (!strcmp(option, "--segtype")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%u",&asegPriorType);
      nargsused = 1;
    }
    else if (!strcmp(option, "--segprior")) {
      if (nargc < 3) CMDargNErr(option,2);
      asegPriorFile0 = fio_fullpath(pargv[0]);
      asegPriorFile1 = fio_fullpath(pargv[1]);
      asegIdFile     = fio_fullpath(pargv[2]);
      nargsused = 3;
    }
    else if (!strcmp(option, "--segtrain")) {
      if (nargc < 2) CMDargNErr(option,1);
      asegTrainFile = fio_fullpath(pargv[0]);
      pathTrainFile = fio_fullpath(pargv[1]);
      nargsused = 2;
    }
    else if (!strcmp(option, "--seg")) {
      if (nargc < 1) CMDargNErr(option,1);
      asegFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--nb")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%u",&nBurnIn);
      nargsused = 1;
    }
    else if (!strcmp(option, "--ns")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%u",&nSample);
      nargsused = 1;
    }
    else if (!strcmp(option, "--nk")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%u",&nKeepSample);
      nargsused = 1;
    }
    else if (!strcmp(option, "--nu")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%u",&nUpdateProp);
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
  printf("USAGE: ./dmri_track\n");
  printf("\n");
  printf("Basic inputs\n");
  printf("   --outdir <dir>: output directory\n");
  printf("   --dwi <file>: DWI volume series\n");
  printf("   --grad <file>: text file of diffusion gradients\n");
  printf("   --bval <file>: text file of diffusion b-values\n");
  printf("   --mask <file>: mask volume\n");
  printf("   --bpdir <dir>: BEDPOST directory\n");
  printf("   --ntr <number>: max number of tracts per voxel (default 1)\n");
  printf("   --roi1 <file>: end ROI 1 (volume or label)\n");
  printf("   --roi2 <file>: end ROI 2 (volume or label)\n");
  printf("   --roimesh1 <file>: mesh for end ROI 1 (if label)\n");
  printf("   --roimesh2 <file>: mesh for end ROI 2 (if label)\n");
  printf("   --roiref1 <file>: reference volume for end ROI 1 (if label)\n");
  printf("   --roiref2 <file>: reference volume for end ROI 2 (if label)\n");
  printf("   --dwi2roi <file>: DWI-to-ROI .mat registration file\n");
  printf("   --init <file>: text file of initialization control points\n");
  printf("\n");
  printf("Prior inputs\n");
  printf("   --prior <file0 file1>: spatial path prior\n");
  printf("     (negative log-likelihoods for H0 and H1)\n");
  printf("   --segtype <number>: 1 = local aseg prior\n");
  printf("                       2 = local neighborhood aseg prior\n");
  printf("                       3 = nearest-neighbor aseg prior\n");
  printf("                       4 = distance aseg prior\n");
  printf("   --segprior <file0 file1 fileID>: aseg priors\n");
  printf("     (negative log-likelihoods for H0 and H1, list of labels)\n");
  printf("   --segtrain <pathfile asegfile>: local aseg prior\n");
  printf("     (training sets of segmentation maps and paths)\n");
  printf("   --seg <file>: segmentation map of test subject\n");
  printf("\n");
  printf("MCMC options\n");
  printf("   --nb <number>: number of burn-in samples (default 5000)\n");
  printf("   --ns <number>: number of post-burn-in samples (default 5000)\n");
  printf("   --nk <number>: keep every nk-th sample (default 10)\n");
  printf("   --nu <number>: update proposal every nu-th sample (default 40)\n");
  printf("   --sdp <file>: text file with initial proposal SD's for control point perturbations (default SD=1 for all control points)\n");
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
  if(!outDir) {
    printf("ERROR: must specify output directory\n");
    exit(1);
  }
  if(!dwiFile) {
    printf("ERROR: must specify DWI volume series\n");
    exit(1);
  }
  if(!gradFile) {
    printf("ERROR: must specify gradient text file\n");
    exit(1);
  }
  if(!bvalFile) {
    printf("ERROR: must specify b-value text file\n");
    exit(1);
  }
  if(!maskFile) {
    printf("ERROR: must specify mask volume\n");
    exit(1);
  }
  if(!bedpostDir) {
    printf("ERROR: must specify BEDPOST directory\n");
    exit(1);
  }
  if(!roiFile1 || !roiFile2) {
    printf("ERROR: must specify both end ROIs\n");
    exit(1);
  }
  if(strstr(roiFile1, ".label") && (!roiMeshFile1 || !roiRefFile1)) {
    printf("ERROR: must specify mesh and reference volume for label ROI 1\n");
    exit(1);
  }
  if(strstr(roiFile2, ".label") && (!roiMeshFile2 || !roiRefFile2)) {
    printf("ERROR: must specify mesh and reference volume for label ROI 2\n");
    exit(1);
  }
  if(!initFile) {
    printf("ERROR: must specify control point initialization file\n");
    exit(1);
  }
  if(!asegFile && (asegTrainFile || asegPriorFile0)) {
    printf("ERROR: must specify seg. map file with aseg prior\n");
    exit(1);
  }
  if(!asegPriorType && asegFile) {
    printf("ERROR: must specify type of aseg prior (1-4)\n");
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

  fprintf(fp, "Output directory: %s\n", outDir);
  fprintf(fp, "DWIs: %s\n", dwiFile);
  fprintf(fp, "Gradients: %s\n", gradFile);
  fprintf(fp, "B-values: %s\n", bvalFile);
  fprintf(fp, "Mask: %s\n", maskFile);
  fprintf(fp, "BEDPOST directory: %s\n", bedpostDir);
  fprintf(fp, "Max number of tracts per voxel: %u\n", nTract);
  fprintf(fp, "End ROI 1: %s\n", roiFile1);
  fprintf(fp, "End ROI 2: %s\n", roiFile2);
  if (roiMeshFile1) fprintf(fp, "End ROI 1 mesh: %s\n", roiMeshFile1);
  if (roiMeshFile2) fprintf(fp, "End ROI 2 mesh: %s\n", roiMeshFile2);
  if (roiRefFile1) fprintf(fp, "End ROI 1 reference volume: %s\n", roiRefFile1);
  if (roiRefFile2) fprintf(fp, "End ROI 2 reference volume: %s\n", roiRefFile2);
  if (xfmFile) fprintf(fp, "DWI-to-ROI: %s\n", xfmFile);
  fprintf(fp, "Initial control point file: %s\n", initFile);
  if (priorFile0) {
    fprintf(fp, "Spatial prior (off path): %s\n", priorFile0);
    fprintf(fp, "Spatial prior (on path): %s\n", priorFile1);
  }
  if (asegPriorFile0) {
    fprintf(fp, "Aseg prior (off path): %s\n", asegPriorFile0);
    fprintf(fp, "Aseg prior (on path): %s\n", asegPriorFile1);
    fprintf(fp, "Aseg label ID list: %s\n", asegIdFile);
  }
  if (asegTrainFile) {
    fprintf(fp, "Seg. map training set: %s\n", asegTrainFile);
    fprintf(fp, "Path training set: %s\n", pathTrainFile);
  }
  if (asegFile) {
    fprintf(fp, "Segmentation map: %s\n", asegFile);
    fprintf(fp, "Type of aseg prior: %u\n", asegPriorType);
  }
  fprintf(fp, "Number of burn-in samples: %u\n", nBurnIn);
  fprintf(fp, "Number of post-burn-in samples: %u\n", nSample);
  fprintf(fp, "Keep every: %u-th sample\n", nKeepSample);
  fprintf(fp, "Update proposal every: %u-th sample\n", nUpdateProp);
  if (stdPropFile) fprintf(fp, "Initial proposal SD file: %s\n", stdPropFile);

  return;
}

