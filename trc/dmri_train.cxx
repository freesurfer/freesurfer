/**
 * @brief Training of priors for probabilistic global tractography
 *
 * Training of priors for probabilistic global tractography
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

#include "blood.h"

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

const char *Progname = "dmri_train";

bool useTrunc = false, useAnatomy = false, useShape = false, excludeStr = false;
int numStrMax = INT_MAX;
vector<float> trainMaskLabel;
vector< vector<int> > nControl;
vector<string> outTrkList, outPriorBase,
               trainTrkList, trainRoi1List, trainRoi2List,
               testMaskList, testFaList, testBaseXfmList;
string outPriorDir, outTestDir,
       trainListFile, trainAsegFile, trainMaskFile,
       testAffineXfmFile, testNonlinXfmFile,
       testNonlinRefFile, testBaseMaskFile;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  string excfile, fbase;

  nargs = handleVersionOption(argc, argv, "dmri_train");
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

  if (excludeStr) {
    excfile = outPriorBase[0] + "_cpts_all.bad.txt";

    if (!outPriorDir.empty())
      excfile = outPriorDir + "/" + excfile;
  }

  if (trainRoi1List.empty()) {
    trainRoi1List.resize(trainTrkList.size());
    fill(trainRoi1List.begin(), trainRoi1List.end(), string());
  }

  if (trainRoi2List.empty()) {
    trainRoi2List.resize(trainTrkList.size());
    fill(trainRoi2List.begin(), trainRoi2List.end(), string());
  }

  if (trainMaskLabel.empty()) {
    trainMaskLabel.resize(trainTrkList.size());
    fill(trainMaskLabel.begin(), trainMaskLabel.end(), 0.0f);
  }

  if (nControl.empty())
    nControl.push_back(vector<int>());

  Blood myblood(trainListFile, trainTrkList[0],
                trainRoi1List[0], trainRoi2List[0],
                trainAsegFile, trainMaskFile, trainMaskLabel[0],
                excfile,
                testMaskList, testFaList,
                testAffineXfmFile, testNonlinXfmFile, testNonlinRefFile,
                testBaseXfmList, testBaseMaskFile,
                useTrunc, useAnatomy, useShape, nControl[0], numStrMax,
                debug);

  for (unsigned int itrk = 0; itrk < trainTrkList.size(); itrk++) {
    if (itrk > 0) {
      if (excludeStr) {
        excfile = outPriorBase[itrk] + "_cpts_all.bad.txt";

        if (!outPriorDir.empty())
          excfile = outPriorDir + "/" + excfile;
      }

      if (nControl.size() > 1)		// Variable number of control points
        myblood.SetNumControls(nControl[itrk]);

      myblood.ReadStreamlines(trainListFile, trainTrkList[itrk],
                              trainRoi1List[itrk], trainRoi2List[itrk],
                              trainMaskLabel[itrk], excfile);
    }

    cout << "Processing pathway " << itrk+1 << " of " << trainTrkList.size()
         << "..." << endl;
    cputimer.reset();

    if (!outTrkList.empty()) {		// Prep training streamlines
      myblood.PrepStreamlines();

      myblood.WriteStreamlines(trainListFile, outTrkList[itrk]);
    }
    else {				// Compute priors for test subject
      myblood.ComputePriors();

      fbase = outPriorBase[itrk];

      if (!outPriorDir.empty())
        fbase = outPriorDir + "/" + fbase;

      if (!outTestDir.empty()) {
        string ftbase = outTestDir + "/" + outPriorBase[itrk];

        myblood.WriteOutputs(fbase.c_str(), ftbase.c_str());
      }
      else
        myblood.WriteOutputs(fbase.c_str());

      cputime = cputimer.milliseconds();
      cout << "Done in " << cputime/1000.0 << " sec." << endl;
    }
  }

  cout << "dmri_train done" << endl;
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
      outPriorDir = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--cptdir")) {
      if (nargc < 1) CMDargNErr(option,1);
      outTestDir = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--out")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        outPriorBase.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--outtrk")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        outTrkList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--slist")) {
      if (nargc < 1) CMDargNErr(option,1);
      trainListFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--trk")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        trainTrkList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--rois")) {
      if (nargc < 2) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        trainRoi1List.push_back(pargv[nargsused]);
        nargsused++;
        trainRoi2List.push_back(pargv[nargsused]);
        nargsused++;
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
      float labid;
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        sscanf(pargv[nargsused], "%f", &labid);
        trainMaskLabel.push_back(labid);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--bmask")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        testMaskList.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--fa")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        testFaList.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--reg")) {
      if (nargc < 1) CMDargNErr(option,1);
      testAffineXfmFile = fio_fullpath(pargv[nargsused]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--regnl")) {
      if (nargc < 1) CMDargNErr(option,1);
      testNonlinXfmFile = fio_fullpath(pargv[nargsused]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--refnl")) {
      if (nargc < 1) CMDargNErr(option,1);
      testNonlinRefFile = fio_fullpath(pargv[nargsused]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--basereg")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        testBaseXfmList.push_back(fio_fullpath(pargv[nargsused]));
        nargsused++;
      }
    }
    else if (!strcmp(option, "--baseref")) {
      if (nargc < 1) CMDargNErr(option,1);
      testBaseMaskFile = fio_fullpath(pargv[nargsused]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--ncpts")) {
      if (nargc < 1) CMDargNErr(option,1);
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        int ncpts;
        vector<int> ncptlist;

        sscanf(pargv[nargsused], "%d", &ncpts);
        ncptlist.push_back(ncpts);
        nControl.push_back(ncptlist);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--max")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &numStrMax);
      nargsused = 1;
    }
    else if (!strcmp(option, "--trunc"))
      useTrunc = true;
    else if (!strcmp(option, "--aprior"))
      useAnatomy = true;
    else if (!strcmp(option, "--sprior"))
      useShape = true;
    else if (!strcmp(option, "--xstr"))
      excludeStr = true;
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
  << "Basic inputs (all volumes in common space)" << endl
  << "   --slist <file>:" << endl
  << "     Text file with list of training subject directories" << endl
  << "   --trk <file> [...]:" << endl
  << "     Name(s) of input .trk file(s), one per path" << endl
  << "     (Names relative to training subject directory)" << endl
  << "   --rois <file1> <file2> [...]:" << endl
  << "     Optional, names of input tract labeling ROIs, two per path" << endl
  << "     (Names relative to training subject directory)" << endl
  << "   --seg <file>:" << endl
  << "     Name of input aparc+aseg volume" << endl
  << "     (Name relative to training subject directory)" << endl
  << "   --cmask <file>:" << endl
  << "     Name of input cortex mask volume" << endl
  << "   --lmask <id> [...]:" << endl
  << "     Add a label ID from aparc+aseg to cortex mask, one per path" << endl
  << "     (0 doesn't add any label)" << endl
  << endl
  << "Inputs used to prep training data (all volumes in common space)" << endl
  << "   --bmask <file> [...]:" << endl
  << "     Input brain mask volume(s)" << endl
  << "   --outtrk <file> [...]:" << endl
  << "     Name(s) of output, pre-sorted .trk file(s), one per path" << endl
  << "     (Names relative to training subject directory)" << endl
  << endl
  << "Inputs used to compute priors (all volumes in common space)" << endl
  << "   --bmask <file> [...]:" << endl
  << "     Input brain mask volume(s) for test subject" << endl
  << "   --fa <file> [...]:" << endl
  << "     Input FA volume(s) for test subject (optional)" << endl
  << "   --reg <file>:" << endl
  << "     Affine registration from atlas to base space (optional)" << endl
  << "   --regnl <file>:" << endl
  << "     Nonlinear registration from atlas to base space (optional)" << endl
  << "   --refnl <file>:" << endl
  << "     Nonlinear registration source reference volume (optional)" << endl
  << "   --basereg <file> [...]:" << endl
  << "     Affine registration(s) from base to FA volume(s) (optional)" << endl
  << "   --baseref <file> [...]:" << endl
  << "     Base space reference volume (optional)" << endl
  << "   --ncpts <num> [...]:" << endl
  << "     Number of control points for initial spline, one per path" << endl
  << "     or one for all paths" << endl
  << "   --max <num>:" << endl
  << "     Maximum number of training streamlines to keep per path" << endl
  << "   --xstr:" << endl
  << "     Exclude previously chosen center streamline(s) (Default: No)" << endl
  << "   --aprior:" << endl
  << "     Compute priors on underlying anatomy (Default: No)" << endl
  << "   --sprior:" << endl
  << "     Compute priors on shape (Default: No)" << endl
  << "   --trunc:" << endl
  << "     Use all training streamlines, truncated or not" << endl
  << "     (Default: Only save results using non-truncated streamlines)" << endl
  << "   --out <base> [...]:" << endl
  << "     Base name(s) of output(s) for test subject, one per path" << endl
  << "   --outdir <dir>:" << endl
  << "     Output directory (optional)" << endl
  << "   --cptdir <dir>:" << endl
  << "     Output directory for control points in test subject's space" << endl
  << "     (optional, requires registration files)" << endl
  << "     If specified, base names of outputs are relative to this" << endl
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
  if (outTrkList.empty() && outPriorBase.empty()) {
    cout << "ERROR: Must specify at least one output name" << endl;
    exit(1);
  }
  if (trainListFile.empty()) {
    cout << "ERROR: Must specify training subject list file" << endl;
    exit(1);
  }
  if (trainTrkList.empty()) {
    cout << "ERROR: Must specify location of at least one streamline file"
         << endl;
    exit(1);
  }
  if (!outPriorBase.empty() && (trainTrkList.size() != outPriorBase.size())) {
    cout << "ERROR: Numbers of input .trk files and output names must match"
         << endl;
    exit(1);
  }
  if (!outTrkList.empty() && (trainTrkList.size() != outTrkList.size())) {
    cout << "ERROR: Numbers of input and output .trk files must match"
         << endl;
    exit(1);
  }
  if (trainAsegFile.empty()) {
    cout << "ERROR: Must specify location of aparc+aseg volume" << endl;
    exit(1);
  }
  if (trainMaskFile.empty()) {
    cout << "ERROR: Must specify location of cortex mask volume" << endl;
    exit(1);
  }
  if (!trainRoi1List.empty() && trainRoi1List.size() != trainTrkList.size()) {
    cout << "ERROR: Numbers of input .trk files and start ROIs must match"
         << endl;
    exit(1);
  }
  if (!trainRoi2List.empty() && trainRoi2List.size() != trainTrkList.size()) {
    cout << "ERROR: Numbers of input .trk files and end ROIs must match"
         << endl;
    exit(1);
  }
  if (!trainMaskLabel.empty() && trainMaskLabel.size() != trainTrkList.size()) {
    cout << "ERROR: Numbers of input .trk files and mask label IDs must match"
         << endl;
    exit(1);
  }
  if (!outPriorBase.empty()) {
    if (!testFaList.empty() && (testFaList.size() != testMaskList.size())) {
      cout << "ERROR: Must specify as many FA volumes as brain masks" << endl;
      exit(1);
    }
    if (nControl.empty()) {
      cout << "ERROR: Must specify number of control points for initial spline"
           << endl;
      exit(1);
    }
    if (nControl.size() > 1 && nControl.size() != trainTrkList.size()) {
      cout << "ERROR: Must specify number of control points for each .trk file"
           << "ERROR: or a single number of control points for all .trk files"
           << endl;
      exit(1);
    }
    if (!testNonlinXfmFile.empty() && testNonlinRefFile.empty()) {
      cout << "ERROR: Must specify source reference volume for nonlinear warp"
           << endl;
      exit(1);
    }
    if (!testBaseXfmList.empty() && testBaseMaskFile.empty()) {
      cout << "ERROR: Must specify reference volume for base space" << endl;
      exit(1);
    }
    if (!testBaseXfmList.empty() &&
        (testBaseXfmList.size() != testFaList.size())) {
      cout << "ERROR: Must specify as many base registrations as FA volumes"
           << endl;
      exit(1);
    }
  }
  return;
}

/* --------------------------------------------- */
static void dump_options() {
  vector<string>::const_iterator istr;

  cout << endl
       << getVersion() << endl
       << "cwd " << cwd << endl
       << "cmdline " << cmdline << endl
       << "sysname  " << uts.sysname << endl
       << "hostname " << uts.nodename << endl
       << "machine  " << uts.machine << endl
       << "user     " << VERuser() << endl;

  if (!outPriorDir.empty())
    cout << "Output directory: " << outPriorDir << endl;

  if (!outTestDir.empty())
    cout << "Output directory in test subject's space: " << outTestDir << endl;

  if (!outPriorBase.size()) {
    cout << "Output base for anatomical priors:";
    for (istr = outPriorBase.begin(); istr < outPriorBase.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }

  cout << "Training subject directory list: " << trainListFile << endl;

  if (!outTrkList.empty()) {
    cout << "Location of output, pre-sorted streamline files relative to base:";
    for (istr = outTrkList.begin(); istr < outTrkList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }

  cout << "Location of streamline files relative to base:";
  for (istr = trainTrkList.begin(); istr < trainTrkList.end(); istr++)
    cout << " " << *istr;
  cout << endl;

  if (!trainRoi1List.empty()) {
    cout << "Location of start ROI files relative to base:";
    for (istr = trainRoi1List.begin(); istr < trainRoi1List.end(); istr++)
      cout << " " << *istr;
    cout << endl;

    cout << "Location of end ROI files relative to base:";
    for (istr = trainRoi2List.begin(); istr < trainRoi2List.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }

  cout << "Location of cortex masks relative to base: " << trainMaskFile
       << endl;

  if (!trainMaskLabel.empty()) {
    cout << "Label ID's from aparc+aseg to add to cortex mask:";
    for (vector<float>::const_iterator ilab = trainMaskLabel.begin();
                                       ilab < trainMaskLabel.end(); ilab++)
      cout << " " << (int) *ilab;
    cout << endl;
  }

  cout << "Location of aparc+aseg's relative to base: " << trainAsegFile
       << endl;

  cout << "Brain mask for output subject:";
  for (vector<string>::const_iterator ifile = testMaskList.begin();
                                      ifile < testMaskList.end(); ifile++)
    cout << " " << *ifile;
  cout << endl;

  if (!testFaList.empty()) {
    cout << "FA map for output subject:";
    for (vector<string>::const_iterator ifile = testFaList.begin();
                                        ifile < testFaList.end(); ifile++)
      cout << " " << *ifile;
    cout << endl;
  }

  if (!testAffineXfmFile.empty())
    cout << "Affine registration from atlas to base for output subject: "
         << testAffineXfmFile << endl;

  if (!testNonlinXfmFile.empty())
    cout << "Nonlinear registration from atlas to base for output subject: "
         << testNonlinXfmFile << endl;

  if (!testNonlinRefFile.empty())
    cout << "Nonlinear registration source reference for output subject: "
         << testNonlinRefFile << endl;

  if (!testBaseXfmList.empty()) {
    cout << "Affine registration from base to FA map for output subject:";
    for (vector<string>::const_iterator ifile = testBaseXfmList.begin();
                                        ifile < testBaseXfmList.end(); ifile++)
      cout << " " << *ifile;
    cout << endl;
  }

  if (!testBaseMaskFile.empty())
    cout << "Base mask for output subject: " << testBaseMaskFile << endl;

  cout << "Number of control points for initial spline:";
  for (vector< vector<int> >::const_iterator inlist = nControl.begin();
                                             inlist < nControl.end(); inlist++)
    for (vector<int>::const_iterator inum = inlist->begin();
                                     inum < inlist->end(); inum++)
      cout << " " << *inum;
  cout << endl;

  if (numStrMax < INT_MAX)
    cout << "Maximum number of training streamlines per path: "
         << numStrMax << endl;

  cout << "Exclude previously chosen center streamlines: " << excludeStr << endl
       << "Use truncated streamlines: " << useTrunc << endl
       << "Compute priors on underlying anatomy: " << useAnatomy << endl
       << "Compute priors on shape: " << useShape << endl;

  return;
}

