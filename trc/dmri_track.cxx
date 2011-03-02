/**
 * @file  dmri_track.c
 * @brief Probabilistic global tractography
 *
 * Probabilistic global tractography
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:41 $
 *    $Revision: 1.4 $
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

#include "coffin.h"

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

static char vcid[] = "";
const char *Progname = "dmri_track";

unsigned int nTract = 1, 
             nBurnIn = 5000, nSample = 5000, nKeepSample = 10, nUpdateProp = 40,
             asegPriorType = 0;
int nout = 0, ninit = 0, nroi1 = 0, nroi2 = 0, nlab1 = 0, nlab2 = 0,
    nmesh1 = 0, nmesh2 = 0, nref1 = 0, nref2 = 0,
    npri = 0, nseg = 0, ntra = 0, nnei = 0, nloc = 0, nsdp = 0;
float fminPath = 0;
char *outDir[100], *dwiFile = NULL,
     *gradFile = NULL, *bvalFile = NULL,
     *maskFile = NULL, *bedpostDir = NULL,
     *initFile[100],
     *roiFile1[100], *roiFile2[100],
     *roiMeshFile1[100], *roiMeshFile2[100],
     *roiRefFile1[100], *roiRefFile2[100],
     *priorFile0[100], *priorFile1[100],
     *asegPriorFile0[100], *asegPriorFile1[100],
     *asegIdFile[100],
     *asegTrainFile[100], *pathTrainFile[100],
     *neighPriorFile[100], *neighIdFile[100],
     *localPriorFile[100], *localIdFile[100],
     *asegFile = NULL,
     *affineXfmFile = NULL, *nonlinXfmFile = NULL,
     *stdPropFile[100];

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

  Coffin mycoffin(outDir[0], dwiFile,
                  gradFile, bvalFile,
                  maskFile, bedpostDir,
                  nTract, fminPath,
                  initFile[0],
                  roiFile1[0], roiFile2[0],
                  nlab1 ? roiMeshFile1[0] : 0, nlab2 ? roiMeshFile2[0] : 0,
                  nlab1 ? roiRefFile1[0] : 0, nlab2 ? roiRefFile2[0] : 0,
                  npri ? priorFile0[0] : 0, npri ? priorFile1[0] : 0,
                  asegPriorType,
                  nseg ? asegPriorFile0[0] : 0, nseg ? asegPriorFile1[0] : 0,
                  nseg ? asegIdFile[0] : 0,
                  ntra ? asegTrainFile[0] : 0, ntra ? pathTrainFile[0] : 0,
                  nnei ? neighPriorFile[0] : 0, nnei ? neighIdFile[0] : 0,
                  nloc ? localPriorFile[0] : 0, nloc ? localIdFile[0] : 0,
                  asegFile,
                  affineXfmFile, nonlinXfmFile,
                  nBurnIn, nSample,
                  nKeepSample, nUpdateProp,
                  nsdp ? stdPropFile[0] : 0,
                  debug);

  for (int k = 0; k < nout; k++) {
    if (k > 0) {
      mycoffin.SetOutputDir(outDir[k]);
      mycoffin.SetPathway(initFile[k],
                  roiFile1[k], roiFile2[k],
                  nlab1 ? roiMeshFile1[k] : 0, nlab2 ? roiMeshFile2[k] : 0,
                  nlab1 ? roiRefFile1[k] : 0, nlab2 ? roiRefFile2[k] : 0,
                  npri ? priorFile0[k] : 0, npri ? priorFile1[k] : 0,
                  asegPriorType,
                  nseg ? asegPriorFile0[k] : 0, nseg ? asegPriorFile1[k] : 0,
                  nseg ? asegIdFile[k] : 0,
                  ntra ? asegTrainFile[k] : 0, ntra ? pathTrainFile[k] : 0,
                  nnei ? neighPriorFile[k] : 0, nnei ? neighIdFile[k] : 0,
                  nloc ? localPriorFile[k] : 0, nloc ? localIdFile[k] : 0);
      mycoffin.SetMCMCParameters(nBurnIn, nSample,
                  nKeepSample, nUpdateProp,
                  nsdp ? stdPropFile[k] : 0);
    }

    printf("Processing pathway %d of %d...\n", k+1, nout);
    TimerStart(&cputimer);

    //if (mycoffin.RunMCMC())
if (mycoffin.RunMCMC1())
      mycoffin.WriteOutputs();
    else
      cout << "ERROR: Pathway reconstruction failed" << endl;

    cputime = TimerStop(&cputimer);
    printf("Done in %g sec.\n", cputime/1000.0);
  }

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
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        outDir[nout] = pargv[nargsused];
        nargsused++;
        nout++;
      }
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
    else if (!strcmp(option, "--fmin")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&fminPath);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--roi1")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        roiFile1[nroi1] = fio_fullpath(pargv[nargsused]);
        if(strstr(roiFile1[nroi1], ".label"))
          nlab1++;
        nargsused++;
        nroi1++;
      }
    } 
    else if (!strcmp(option, "--roi2")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        roiFile2[nroi2] = fio_fullpath(pargv[nargsused]);
        if(strstr(roiFile2[nroi2], ".label"))
          nlab2++;
        nargsused++;
        nroi2++;
      }
    } 
    else if (!strcmp(option, "--roimesh1")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        roiMeshFile1[nmesh1] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        nmesh1++;
      }
    } 
    else if (!strcmp(option, "--roimesh2")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        roiMeshFile2[nmesh2] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        nmesh2++;
      }
      nmesh2 += nargsused;
    } 
    else if (!strcmp(option, "--roiref1")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        roiRefFile1[nref1] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        nref1++;
      }
    } 
    else if (!strcmp(option, "--roiref2")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        roiRefFile2[nref2] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        nref2++;
      }
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
    else if (!strcmp(option, "--init")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        initFile[ninit] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        ninit++;
      }
    }
    else if (!strcmp(option, "--sdp")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        stdPropFile[nsdp] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        nsdp++;
      }
    }
    else if (!strcmp(option, "--prior")) {
      if (nargc < 2) CMDargNErr(option,2);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        priorFile0[npri] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        priorFile1[npri] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        npri++;
      }
    }
    else if (!strcmp(option, "--segtype")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%u",&asegPriorType);
      nargsused = 1;
    }
    else if (!strcmp(option, "--segprior")) {
      if (nargc < 3) CMDargNErr(option,3);
      while (strncmp(pargv[nargsused], "--", 2)) {
        asegPriorFile0[nseg] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        asegPriorFile1[nseg] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        asegIdFile[nseg]     = fio_fullpath(pargv[nargsused]);
        nargsused++;
        nseg++;
      }
    }
    else if (!strcmp(option, "--segtrain")) {
      if (nargc < 2) CMDargNErr(option,2);
      nargsused = 0;
      while (strncmp(pargv[nargsused], "--", 2)) {
        asegTrainFile[ntra] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        pathTrainFile[ntra] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        ntra++;
      }
    }
    else if (!strcmp(option, "--nprior")) {
      if (nargc < 2) CMDargNErr(option,2);
      while (strncmp(pargv[nargsused], "--", 2)) {
        neighPriorFile[nnei] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        neighIdFile[nnei]    = fio_fullpath(pargv[nargsused]);
        nargsused++;
        nnei++;
      }
    }
    else if (!strcmp(option, "--lprior")) {
      if (nargc < 2) CMDargNErr(option,2);
      while (strncmp(pargv[nargsused], "--", 2)) {
        localPriorFile[nloc] = fio_fullpath(pargv[nargsused]);
        nargsused++;
        localIdFile[nloc]    = fio_fullpath(pargv[nargsused]);
        nargsused++;
        nloc++;
      }
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
  printf("   --outdir <dir> [...]: output directory (one per path)\n");
  printf("   --dwi <file>: DWI volume series\n");
  printf("   --grad <file>: text file of diffusion gradients\n");
  printf("   --bval <file>: text file of diffusion b-values\n");
  printf("   --mask <file>: mask volume\n");
  printf("   --bpdir <dir>: BEDPOST directory\n");
  printf("   --ntr <number>: max number of tracts per voxel (default 1)\n");
  printf("   --fmin <number>: tract volume fraction threshold (default 0)\n");
  printf("   --init <file> [...]: text file of initialization control points (one per path)\n");
  printf("   --roi1 <file> [...]: end ROI 1 (volume or label, one per path)\n");
  printf("   --roi2 <file> [...]: end ROI 2 (volume or label, one per path)\n");
  printf("   --roimesh1 <file> [...]: mesh for end ROI 1 (for label ROIs)\n");
  printf("   --roimesh2 <file> [...]: mesh for end ROI 2 (for label ROIs)\n");
  printf("   --roiref1 <file> [...]: reference volume for end ROI 1 (for label ROIs)\n");
  printf("   --roiref2 <file> [...]: reference volume for end ROI 2 (for label ROIs)\n");
  printf("\n");
  printf("Prior inputs\n");
  printf("   --prior <file0 file1> [...]: spatial path prior\n");
  printf("     (negative log-likelihoods for H0 and H1, one pair per path)\n");
  printf("   --segtype <number>: 1 = local aseg prior\n");
  printf("                       2 = local neighborhood aseg prior\n");
  printf("                       3 = nearest-neighbor aseg prior\n");
  printf("                       4 = distance aseg prior\n");
  printf("   --segprior <file0 file1 fileID> [...]: aseg priors\n");
  printf("     (negative log-likelihoods for H0 and H1 and list of labels, one triplet per path)\n");
  printf("   --segtrain <pathfile asegfile> [...]: local aseg prior\n");
  printf("     (training sets of segmentation maps and paths, one pair per path)\n");
  printf("   --nprior <priorfile idfile> [...]: near neighbor aseg priors\n");
  printf("     (negative log-likelihood and list of labels, one pair per path)\n");
  printf("   --lprior <priorfile idfile> [...]: local neighbor aseg priors\n");
  printf("     (negative log-likelihood and list of labels, one pair per path)\n");
  printf("   --seg <file>: segmentation map of test subject\n");
  printf("   --reg <file>:  DWI-to-atlas affine registration (.mat)\n");
  printf("   --regnl <file>: DWI-to-atlas nonlinear registration (.tm3d)\n");
  printf("\n");
  printf("MCMC options\n");
  printf("   --nb <number>: number of burn-in samples (default 5000)\n");
  printf("   --ns <number>: number of post-burn-in samples (default 5000)\n");
  printf("   --nk <number>: keep every nk-th sample (default 10)\n");
  printf("   --nu <number>: update proposal every nu-th sample (default 40)\n");
  printf("   --sdp <file> [...]: text file with initial proposal SD's for control point perturbations (one per path or default SD=1 for all control points)\n");
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
  if(ninit != nout) {
    printf("ERROR: must specify as many control point initialization files as outputs\n");
    exit(1);
  }
  if(nroi1 != nout) {
    printf("ERROR: must specify as many end ROI 1's as outputs\n");
    exit(1);
  }
  if(nroi2 != nout) {
    printf("ERROR: must specify as many end ROI 2's as outputs\n");
    exit(1);
  }
  if(nmesh1 != nlab1) {
    printf("ERROR: must specify as many meshes as labels for ROI 1\n");
    exit(1);
  }
  if(nref1 != nlab1) {
    printf("ERROR: must specify as many reference volumes as labels for ROI 1\n");
    exit(1);
  }
  if(nmesh2 != nlab2) {
    printf("ERROR: must specify as many meshes as labels for ROI 2\n");
    exit(1);
  }
  if(nref2 != nlab2) {
    printf("ERROR: must specify as many reference volumes as labels for ROI 2\n");
    exit(1);
  }
  if(npri > 0 && npri != nout) {
    printf("ERROR: must specify as many spatial prior pairs as outputs\n");
    exit(1);
  }
  if(nseg > 0 && nseg != nout) {
    printf("ERROR: must specify as many aseg prior triplets as outputs\n");
    exit(1);
  }
  if(ntra > 0 && ntra != nout) {
    printf("ERROR: must specify as many training paths and aseg's as outputs\n");
    exit(1);
  }
  if(nnei > 0 && nnei != nout) {
    printf("ERROR: must specify as many neighbor aseg prior pairs as outputs\n");
    exit(1);
  }
  if(nloc > 0 && nloc != nout) {
    printf("ERROR: must specify as many local aseg prior pairs as outputs\n");
    exit(1);
  }
  if(!asegFile && (nseg > 0 || ntra > 0 || nnei > 0 || nloc > 0)) {
    printf("ERROR: must specify segmentation map file with aseg prior\n");
    exit(1);
  }
  if(nsdp > 0 && nsdp != nout) {
    printf("ERROR: must specify as many control point proposal SD files as outputs\n");
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

  fprintf(fp, "Output directory:");
  for (int k = 0; k < nout; k++)
    fprintf(fp, " %s", outDir[k]);
  fprintf(fp, "\n");
  fprintf(fp, "DWIs: %s\n", dwiFile);
  fprintf(fp, "Gradients: %s\n", gradFile);
  fprintf(fp, "B-values: %s\n", bvalFile);
  fprintf(fp, "Mask: %s\n", maskFile);
  fprintf(fp, "BEDPOST directory: %s\n", bedpostDir);
  fprintf(fp, "Max number of tracts per voxel: %u\n", nTract);
  fprintf(fp, "Tract volume fraction threshold: %g\n", fminPath);
  fprintf(fp, "Initial control point file:");
  for (int k = 0; k < ninit; k++)
    fprintf(fp, " %s", initFile[k]);
  fprintf(fp, "\n");
  fprintf(fp, "End ROI 1:");
  for (int k = 0; k < nroi1; k++)
    fprintf(fp, " %s", roiFile1[k]);
  fprintf(fp, "\n");
  if (nlab1 > 0) {
    fprintf(fp, "End ROI 1 mesh:");
    for (int k = 0; k < nlab1; k++)
      fprintf(fp, " %s", roiMeshFile1[k]);
    fprintf(fp, "\n");
    fprintf(fp, "End ROI 1 reference volume:");
    for (int k = 0; k < nlab1; k++)
      fprintf(fp, " %s", roiRefFile1[k]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "End ROI 2:");
  for (int k = 0; k < nroi2; k++)
    fprintf(fp, " %s", roiFile2[k]);
  fprintf(fp, "\n");
  if (nlab2 > 0) {
    fprintf(fp, "End ROI 2 mesh:");
    for (int k = 0; k < nlab2; k++)
      fprintf(fp, " %s", roiMeshFile2[k]);
    fprintf(fp, "\n");
    fprintf(fp, "End ROI 2 reference volume:");
    for (int k = 0; k < nlab2; k++)
      fprintf(fp, " %s", roiRefFile2[k]);
    fprintf(fp, "\n");
  }
  if (npri > 0) {
    fprintf(fp, "Spatial prior (off path):");
    for (int k = 0; k < npri; k++)
      fprintf(fp, " %s", priorFile0[k]);
    fprintf(fp, "\n");
    fprintf(fp, "Spatial prior (on path):");
    for (int k = 0; k < npri; k++)
      fprintf(fp, " %s", priorFile1[k]);
    fprintf(fp, "\n");
  }
  if (nseg > 0) {
    fprintf(fp, "Aseg prior (off path):");
    for (int k = 0; k < nseg; k++)
      fprintf(fp, " %s", asegPriorFile0[k]);
    fprintf(fp, "\n");
    fprintf(fp, "Aseg prior (on path):");
    for (int k = 0; k < nseg; k++)
      fprintf(fp, " %s", asegPriorFile1[k]);
    fprintf(fp, "\n");
    fprintf(fp, "Aseg label ID list:");
    for (int k = 0; k < nseg; k++)
      fprintf(fp, " %s", asegIdFile[k]);
    fprintf(fp, "\n");
  }
  if (ntra > 0) {
    fprintf(fp, "Seg. map training set:");
    for (int k = 0; k < ntra; k++)
      fprintf(fp, " %s", asegTrainFile[k]);
    fprintf(fp, "\n");
    fprintf(fp, "Path training set:");
    for (int k = 0; k < ntra; k++)
      fprintf(fp, " %s", pathTrainFile[k]);
    fprintf(fp, "\n");
  }
  if (nnei > 0) {
    fprintf(fp, "Neighbor aseg prior:");
    for (int k = 0; k < nnei; k++)
      fprintf(fp, " %s", neighPriorFile[k]);
    fprintf(fp, "\n");
    fprintf(fp, "Neighbor aseg label ID list:");
    for (int k = 0; k < nnei; k++)
      fprintf(fp, " %s", neighIdFile[k]);
    fprintf(fp, "\n");
  }
  if (nloc > 0) {
    fprintf(fp, "Local aseg prior:");
    for (int k = 0; k < nloc; k++)
      fprintf(fp, " %s", localPriorFile[k]);
    fprintf(fp, "\n");
    fprintf(fp, "Local aseg label ID list:");
    for (int k = 0; k < nloc; k++)
      fprintf(fp, " %s", localIdFile[k]);
    fprintf(fp, "\n");
  }
  if (asegFile) {
    fprintf(fp, "Segmentation map: %s\n", asegFile);
    fprintf(fp, "Type of aseg prior: %u\n", asegPriorType);
  }
  if (affineXfmFile)
    fprintf(fp, "DWI-to-atlas affine registration: %s\n", affineXfmFile);
  if (nonlinXfmFile)
    fprintf(fp, "DWI-to-atlas nonlinear registration: %s\n", nonlinXfmFile);
  fprintf(fp, "Number of burn-in samples: %u\n", nBurnIn);
  fprintf(fp, "Number of post-burn-in samples: %u\n", nSample);
  fprintf(fp, "Keep every: %u-th sample\n", nKeepSample);
  fprintf(fp, "Update proposal every: %u-th sample\n", nUpdateProp);
  if (nsdp > 0) {
    fprintf(fp, "Initial proposal SD file:");
    for (int k = 0; k < nsdp; k++)
      fprintf(fp, " %s", stdPropFile[k]);
    fprintf(fp, "\n");
  }

  return;
}

