/**
 * @brief Probabilistic global tractography
 *
 * Probabilistic global tractography
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
static void print_usage(void);
static void usage_exit(void);
static void print_help(void);
static void print_version(void);
static void dump_options();

int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]);

const char *Progname = "dmri_paths";

unsigned int nlab1 = 0, nlab2 = 0;
unsigned int nTract = 1, 
             nBurnIn = 5000, nSample = 5000, nKeepSample = 10, nUpdateProp = 40,
             localPriorSet = 15, neighPriorSet = 14;
float fminPath = 0;
std::string dwiFile, gradFile, bvalFile, maskFile, bedpostDir,
  baseXfmFile, baseMaskFile, affineXfmFile, nonlinXfmFile;
vector<std::string> outDir, inDirList, initFile, roiFile1, roiFile2,
               roiMeshFile1, roiMeshFile2, roiRefFile1, roiRefFile2,
               xyzPriorFile0, xyzPriorFile1, tangPriorFile, curvPriorFile,
               neighPriorFile, neighIdFile, localPriorFile, localIdFile,
               asegList, stdPropFile;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  bool doxyzprior = true,
       dotangprior = true,
       docurvprior = true,
       doneighprior = true,
       dolocalprior = true,
       dopropinit = true;
  int nargs, cputime, ilab1 = 0, ilab2 = 0;

  nargs = handleVersionOption(argc, argv, "dmri_paths");
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

  srand(6875);
  srand48(6875);

  if (xyzPriorFile0.empty())  doxyzprior = false;
  if (tangPriorFile.empty())  dotangprior = false;
  if (curvPriorFile.empty())  docurvprior = false;
  if (neighPriorFile.empty()) doneighprior = false;
  if (localPriorFile.empty()) dolocalprior = false;
  if (stdPropFile.empty())    dopropinit = false;

  Coffin mycoffin(outDir[0], inDirList, dwiFile,
                  gradFile, bvalFile,
                  maskFile, bedpostDir,
                  nTract, fminPath,
                  baseXfmFile, baseMaskFile,
                  initFile[0],
                  roiFile1[0], roiFile2[0],
                  strstr(roiFile1[0].c_str(), ".label") ? roiMeshFile1[ilab1] : std::string(),
                  strstr(roiFile2[0].c_str(), ".label") ? roiMeshFile2[ilab2] : std::string(),
                  strstr(roiFile1[0].c_str(), ".label") ? roiRefFile1[ilab1] : std::string(),
                  strstr(roiFile2[0].c_str(), ".label") ? roiRefFile2[ilab2] : std::string(),
                  doxyzprior ? xyzPriorFile0[0] : std::string(),
                  doxyzprior ? xyzPriorFile1[0] : std::string(),
                  dotangprior ? tangPriorFile[0] : std::string(),
                  docurvprior ? curvPriorFile[0] : std::string(),
                  doneighprior ? neighPriorFile[0] : std::string(),
                  doneighprior ? neighIdFile[0] : std::string(),
                  doneighprior ? neighPriorSet : 0,
                  dolocalprior ? localPriorFile[0] : std::string(),
                  dolocalprior ? localIdFile[0] : std::string(),
                  dolocalprior ? localPriorSet : 0,
                  asegList,
                  affineXfmFile, nonlinXfmFile,
                  nBurnIn, nSample, nKeepSample, nUpdateProp,
                  dopropinit ? stdPropFile[0] : std::string(),
                  debug);

  if (strstr(roiFile1[0].c_str(), ".label")) ilab1++;
  if (strstr(roiFile2[0].c_str(), ".label")) ilab2++;

  for (unsigned int iout = 0; iout < outDir.size(); iout++) {
    if (iout > 0) {
      mycoffin.SetOutputDir(outDir[iout]);
      mycoffin.SetPathway(initFile[iout],
			  roiFile1[iout], roiFile2[iout],
			  strstr(roiFile1[iout].c_str(), ".label") ? roiMeshFile1[ilab1] : std::string(),
			  strstr(roiFile2[iout].c_str(), ".label") ? roiMeshFile2[ilab2] : std::string(),
			  strstr(roiFile1[iout].c_str(), ".label") ? roiRefFile1[ilab1] : std::string(),
			  strstr(roiFile2[iout].c_str(), ".label") ? roiRefFile2[ilab2] : std::string(),
			  doxyzprior ? xyzPriorFile0[iout] : std::string(),
			  doxyzprior ? xyzPriorFile1[iout] : std::string(),
			  dotangprior ? tangPriorFile[iout] : std::string(),
			  docurvprior ? curvPriorFile[iout] : std::string(),
			  doneighprior ? neighPriorFile[iout] : std::string(),
			  doneighprior ? neighIdFile[iout] : std::string(),
			  dolocalprior ? localPriorFile[iout] : std::string(),
			  dolocalprior ? localIdFile[iout] : std::string());
      mycoffin.SetMcmcParameters(nBurnIn, nSample, nKeepSample, nUpdateProp,
				 dopropinit ? stdPropFile[iout] : std::string());

      if (strstr(roiFile1.at(iout).c_str(), ".label")) ilab1++;
      if (strstr(roiFile2.at(iout).c_str(), ".label")) ilab2++;
    }

    cout << "Processing pathway " << iout+1 << " of " << outDir.size() << "..."
         << endl;
    cputimer.reset();

    //if (mycoffin.RunMcmcFull())
    if (mycoffin.RunMcmcSingle()) {
      mycoffin.WriteOutputs();
    } else {
      cout << "ERROR: Pathway reconstruction failed" << endl;
    }

    cputime = cputimer.milliseconds();
    printf("Done in %g sec.\n", cputime/1000.0);
  }

  printf("dmri_paths done\n");
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
    else if (!strcmp(option, "--outdir")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        outDir.push_back(pargv[nargsused]);
        nargsused++;
      }
    } 
    else if (!strcmp(option, "--indir")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        inDirList.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    } 
    else if (!strcmp(option, "--dwi")) {
      if (nargc < 1) CMDargNErr(option,1);
      dwiFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--grad")) {
      if (nargc < 1) CMDargNErr(option,1);
      gradFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--bval")) {
      if (nargc < 1) CMDargNErr(option,1);
      bvalFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--bpdir")) {
      if (nargc < 1) CMDargNErr(option,1);
      bedpostDir = pargv[0];
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
    else if (!strcmp(option, "--basereg")) {
      if (nargc < 1) CMDargNErr(option,1);
      baseXfmFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--basemask")) {
      if (nargc < 1) CMDargNErr(option,1);
      baseMaskFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--roi1")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        roiFile1.push_back(fio_fullpath(pargv[nargsused]));
        if(strstr(roiFile1.back().c_str(), ".label"))
          nlab1++;
        nargsused++;
      }
    } 
    else if (!strcmp(option, "--roi2")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        roiFile2.push_back(fio_fullpath(pargv[nargsused]));
        if(strstr(roiFile2.back().c_str(), ".label"))
          nlab2++;
        nargsused++;
      }
    } 
    else if (!strcmp(option, "--roimesh1")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        roiMeshFile1.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    } 
    else if (!strcmp(option, "--roimesh2")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        roiMeshFile2.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    } 
    else if (!strcmp(option, "--roiref1")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        roiRefFile1.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    } 
    else if (!strcmp(option, "--roiref2")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        roiRefFile2.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
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
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        initFile.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--sdp")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        stdPropFile.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--prior")) {
      if (nargc < 2) CMDargNErr(option,2);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        xyzPriorFile0.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
        xyzPriorFile1.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--nprior")) {
      if (nargc < 2) CMDargNErr(option,2);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        neighPriorFile.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
        neighIdFile.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--nset")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%u",&neighPriorSet);
      nargsused = 1;
    }
    else if (!strcmp(option, "--lprior")) {
      if (nargc < 2) CMDargNErr(option,2);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        localPriorFile.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
        localIdFile.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--lset")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%u",&localPriorSet);
      nargsused = 1;
    }
    else if (!strcmp(option, "--seg")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        asegList.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--tprior")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        tangPriorFile.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--cprior")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        curvPriorFile.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
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
static void print_usage(void) {
  cout
  << endl << "USAGE: " << Progname << endl << endl
  << "Basic inputs (native DWI space)" << endl
  << "   --indir <dir> [...]:" << endl
  << "     Input subject directory (optional)" << endl
  << "     If specified, names of all basic inputs are relative to this" << endl
  << "     Specify multiple input directories for longitudinal data" << endl
  << "   --outdir <dir> [...]:" << endl
  << "     Output directory (one per path)" << endl
  << "   --dwi <file>:" << endl
  << "     DWI volume series" << endl
  << "   --grad <file>:" << endl
  << "     Text file of diffusion gradients" << endl
  << "   --bval <file>:" << endl
  << "     Text file of diffusion b-values" << endl
  << "   --mask <file>:" << endl
  << "     Mask volume" << endl
  << "   --bpdir <dir>:" << endl
  << "     BEDPOST directory" << endl
  << "   --ntr <num>:" << endl
  << "     Max number of tracts per voxel (default 1)" << endl
  << "   --fmin <num>:" << endl
  << "     Tract volume fraction threshold (default 0)" << endl
  << "   --basereg <file> [...]:" << endl
  << "     Base-to-DWI registration, needed for longitudinal data only" << endl
  << "     (.mat, as many as input directories)" << endl
  << endl
  << "Longitudinal inputs (base template space)" << endl
  << "   --basemask <file>:" << endl
  << "     Base template mask volume" << endl
  << endl
  << "End ROIs (atlas space)" << endl
  << "   --roi1 <file> [...]:" << endl
  << "     End ROI 1 (volume or label, one per path)" << endl
  << "   --roi2 <file> [...]:" << endl
  << "     End ROI 2 (volume or label, one per path)" << endl
  << "   --roimesh1 <file> [...]:" << endl
  << "     Mesh for end ROI 1 (for label ROIs)" << endl
  << "   --roimesh2 <file> [...]:" << endl
  << "     Mesh for end ROI 2 (for label ROIs)" << endl
  << "   --roiref1 <file> [...]:" << endl
  << "     Reference volume for end ROI 1 (for label ROIs)" << endl
  << "   --roiref2 <file> [...]:" << endl
  << "     Reference volume for end ROI 2 (for label ROIs)" << endl
  << endl
  << "Prior-related inputs (atlas space)" << endl
  << "   --prior <file0 file1> [...]:" << endl
  << "     Spatial path priors (negative log-likelihoods off and" << endl
  << "     on the path, one pair per path)" << endl
  << "   --nprior <priorfile idfile> [...]:" << endl
  << "     Near-neighbor label priors (negative log-likelihood" << endl
  << "     and list of labels, one pair per path)" << endl
  << "   --nset <num>:" << endl
  << "     Subset of near-neighbor label priors (default all)" << endl
  << "   --lprior <priorfile idfile> [...]:" << endl
  << "     Local-neighbor label priors (negative log-likelihood " << endl
  << "     and list of labels, one pair per path)" << endl
  << "   --lset <num>:" << endl
  << "     Subset of local-neighbor label priors (default all)" << endl
  << "   --seg <file> [...]:" << endl
  << "     Segmentation map of test subject" << endl
  << "     Specify multiple segmentation maps for longitudinal data" << endl
  << "   --tprior <file> [...]:" << endl
  << "     Path tangent vector priors (negative log-likelihood," << endl
  << "     one per path)" << endl
  << "   --cprior <file> [...]:" << endl
  << "     Path curvature priors (negative log-likelihood," << endl
  << "     one per path)" << endl
  << "   --reg <file>:" << endl
  << "     DWI-to-atlas affine registration (.mat)" << endl
  << "   --regnl <file>:" << endl
  << "     DWI-to-atlas nonlinear registration (.m3z)" << endl
  << endl
  << "MCMC options (native diffusion or base template space)" << endl
  << "   --init <file> [...]:" << endl
  << "     Text file of initial control points (one per path)" << endl
  << "   --nb <num>:" << endl
  << "     Number of burn-in samples (default 5000)" << endl
  << "   --ns <num>:" << endl
  << "     Number of post-burn-in samples (default 5000)" << endl
  << "   --nk <num>:" << endl
  << "     Keep every nk-th sample (default 10)" << endl
  << "   --nu <num>:" << endl
  << "     Update proposal every nu-th sample (default 40)" << endl
  << "   --sdp <file> [...]:" << endl
  << "     Text file with initial proposal standard deviations" << endl
  << "     for control point perturbations (one per path or" << endl
  << "     default SD=1 for all control points and all paths)" << endl
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
  if (outDir.empty()) {
    cout << "ERROR: Must specify output directory" << endl;
    exit(1);
  }
  if (dwiFile.empty()) {
    cout << "ERROR: Must specify DWI volume series" << endl;
    exit(1);
  }
  if (gradFile.empty()) {
    cout << "ERROR: Must specify gradient text file" << endl;
    exit(1);
  }
  if (bvalFile.empty()) {
    cout << "ERROR: Must specify b-value text file" << endl;
    exit(1);
  }
  if (maskFile.empty()) {
    cout << "ERROR: Must specify mask volume" << endl;
    exit(1);
  }
  if (bedpostDir.empty()) {
    cout << "ERROR: Must specify BEDPOST directory" << endl;
    exit(1);
  }
  if (initFile.size() != outDir.size()) {
    cout << "ERROR: Must specify as many control point initialization files"
         << " as outputs" << endl;
    exit(1);
  }
  if (roiFile1.size() != outDir.size()) {
    cout << "ERROR: Must specify as many end ROI 1's as outputs" << endl;
    exit(1);
  }
  if (roiFile2.size() != outDir.size()) {
    cout << "ERROR: Must specify as many end ROI 2's as outputs" << endl;
    exit(1);
  }
  if (roiMeshFile1.size() != nlab1) {
    cout << "ERROR: Must specify as many meshes as labels for ROI 1" << endl;
    exit(1);
  }
  if (roiRefFile1.size() != nlab1) {
    cout << "ERROR: Must specify as many reference volumes as labels for ROI 1"          << endl;
    exit(1);
  }
  if (roiMeshFile2.size() != nlab2) {
    cout << "ERROR: Must specify as many meshes as labels for ROI 2" << endl;
    exit(1);
  }
  if (roiRefFile2.size() != nlab2) {
    cout << "ERROR: Must specify as many reference volumes as labels for ROI 2"          << endl;
    exit(1);
  }
  if (!xyzPriorFile0.empty() && xyzPriorFile0.size() != outDir.size()) {
    cout << "ERROR: Must specify as many spatial prior pairs as outputs"                 << endl;
    exit(1);
  }
  if (!neighPriorFile.empty() && neighPriorFile.size() != outDir.size()) {
    cout << "ERROR: Must specify as many neighbor aseg prior pairs as outputs"
         << endl;
    exit(1);
  }
  if (!localPriorFile.empty() && localPriorFile.size() != outDir.size()) {
    cout << "ERROR: Must specify as many local aseg prior pairs as outputs"
         << endl;
    exit(1);
  }
  if (!localPriorFile.empty() &&
      (localPriorSet != 1 && localPriorSet != 7 && localPriorSet != 15)) {
    cout << "ERROR: invalid set of local aseg priors" << endl;
    exit(1);
  }
  if (!neighPriorFile.empty() &&
      (neighPriorSet != 6 && neighPriorSet != 14)) {
    cout << "ERROR: invalid set of neighbor aseg priors" << endl;
    exit(1);
  }
  if (asegList.empty() &&
      (!neighPriorFile.empty() || !localPriorFile.empty())) {
    cout << "ERROR: Must specify segmentation map file with aseg prior" << endl;
    exit(1);
  }
  if (!stdPropFile.empty() && stdPropFile.size() != outDir.size()) {
    cout << "ERROR: Must specify as many control point proposal"
         << " standard deviation files as outputs" << endl;
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options() {
  vector<std::string>::const_iterator istr;

  cout << endl
       << getVersion() << endl
       << "cwd " << cwd << endl
       << "cmdline " << cmdline << endl
       << "sysname  " << uts.sysname << endl
       << "hostname " << uts.nodename << endl
       << "machine  " << uts.machine << endl
       << "user     " << VERuser() << endl;

  cout << "Output directory:";
  for (istr = outDir.begin(); istr < outDir.end(); istr++)
    cout << " " << *istr;
  cout << endl;

  if (!inDirList.empty()) {
    cout << "Input directory:";
    for (istr = inDirList.begin(); istr < inDirList.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }

  cout << "DWIs: " << dwiFile << endl
       << "Gradients: " << gradFile << endl
       << "B-values: " << bvalFile << endl
       << "Mask: " << maskFile << endl
       << "BEDPOST directory: " << bedpostDir << endl
       << "Max number of tracts per voxel: " << nTract << endl
       << "Tract volume fraction threshold: " << fminPath << endl;

  cout << "Initial control point file:";
  for (istr = initFile.begin(); istr < initFile.end(); istr++) {
    cout << " " << *istr;
  }
  cout << endl;

  cout << "End ROI 1:";
  for (istr = roiFile1.begin(); istr < roiFile1.end(); istr++) {
    cout << " " << *istr;
  }
  cout << endl;

  if (nlab1 > 0) {
    cout << "End ROI 1 mesh:";
    for (istr = roiMeshFile1.begin(); istr < roiMeshFile1.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;

    cout << "End ROI 1 reference volume:";
    for (istr = roiRefFile1.begin(); istr < roiRefFile1.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }

  cout << "End ROI 2:";
  for (istr = roiFile2.begin(); istr < roiFile2.end(); istr++) {
    cout << " " << *istr;
  }
  cout << endl;

  if (nlab2 > 0) {
    cout << "End ROI 2 mesh:";
    for (istr = roiMeshFile2.begin(); istr < roiMeshFile2.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;

    cout << "End ROI 2 reference volume:";
    for (istr = roiRefFile2.begin(); istr < roiRefFile2.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }

  if (!xyzPriorFile0.empty()) {
    cout << "Spatial prior (off path):";
    for (istr = xyzPriorFile0.begin(); istr < xyzPriorFile0.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;

    cout << "Spatial prior (on path):";
    for (istr = xyzPriorFile1.begin(); istr < xyzPriorFile1.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }

  if (!neighPriorFile.empty()) {
    cout << "Neighbor aseg prior:";
    for (istr = neighPriorFile.begin(); istr < neighPriorFile.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;

    cout << "Neighbor aseg label ID list:";
    for (istr = neighIdFile.begin(); istr < neighIdFile.end(); istr++) {
      cout << " " << *istr;
    }
    cout <<  endl;

    cout << "Neighbor aseg prior set: " << neighPriorSet << endl;
  }

  if (!localPriorFile.empty()) {
    cout << "Local aseg prior:";
    for (istr = localPriorFile.begin(); istr < localPriorFile.end(); istr++) {
      cout << " " << *istr;
    }
    cout  << endl;

    cout << "Local aseg label ID list:";
    for (istr = localIdFile.begin(); istr < localIdFile.end(); istr++) {
      cout << " " << *istr;
    }
    cout  << endl;

    cout << "Local aseg prior set: " << localPriorSet << endl;
  }

  if (!asegList.empty()) {
    cout << "Segmentation map: ";
    for (istr = asegList.begin(); istr < asegList.end(); istr++) {
      cout << " " << *istr;
    }
    cout  << endl;
  }

  if (!affineXfmFile.empty()) {
    cout << "DWI-to-atlas affine registration: " << affineXfmFile << endl;
  }

  if (!nonlinXfmFile.empty()) {
    cout << "DWI-to-atlas nonlinear registration: " << nonlinXfmFile << endl;
  }

  cout << "Number of burn-in samples: " << nBurnIn << endl
       << "Number of post-burn-in samples: " << nSample << endl
       << "Keep every: " << nKeepSample << "-th sample" << endl
       << "Update proposal every: " << nUpdateProp << "-th sample" << endl;

  if (!stdPropFile.empty()) {
    cout << "Initial proposal SD file:";
    for (istr = stdPropFile.begin(); istr < stdPropFile.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }

  return;
}

