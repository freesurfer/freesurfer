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

bool useTrunc = false, excludeStr = false;
vector<float> trainMaskLabel;
vector<vector<int>> nControl;
vector<std::string> outBase, trainTrkList, trainRoi1List, trainRoi2List,
               testMaskList, testFaList, testBaseXfmList;
std::string outDir, outTestDir, trainListFile,
  trainAsegFile, trainMaskFile, testAffineXfmFile,
  testNonlinXfmFile, testNonlinRefFile, testBaseMaskFile;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  std::string excfile, fbase;

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
    if (!outDir.empty()) {
      excfile = outDir + '/' + outBase.at(0) + "_cpts_all.bad.txt";
    } else {
      excfile = outBase.at(0) + "_cpts_all.bad.txt";
    }
  }
  
  Blood myblood(trainListFile, trainTrkList[0],
                trainRoi1List.size() ? trainRoi1List[0] : std::string(),
                trainRoi2List.size() ? trainRoi2List[0] : std::string(),
                trainAsegFile, trainMaskFile,
                trainMaskLabel.size() ? trainMaskLabel[0] : 0.0f,
                excludeStr ? excfile : std::string(),
                testMaskList, testFaList,
                testAffineXfmFile, testNonlinXfmFile, testNonlinRefFile,
                testBaseXfmList, testBaseMaskFile,
                useTrunc, nControl[0],
                debug);

  for (unsigned int itrk = 0; itrk < trainTrkList.size(); itrk++) {
    if (itrk > 0) {
      if (excludeStr) {
        if (!outDir.empty()) {
	  excfile = outDir + '/' + outBase.at(itrk) + "_cpts_all.bad.txt";
        } else {
	  excfile = outBase.at(itrk) + "_cpts_all.bad.txt";
	}
      }

      if (nControl.size() > 1)		// Variable number of controls
        myblood.SetNumControls(nControl[itrk]);

      myblood.ReadStreamlines(trainListFile, trainTrkList[itrk],
                              trainRoi1List.size() ? trainRoi1List[itrk] : 0,
                              trainRoi2List.size() ? trainRoi2List[itrk] : 0,
                              trainMaskLabel.size() ? trainMaskLabel[itrk] : 0,
                              excludeStr ? excfile : 0);
    }

    cout << "Processing pathway " << itrk+1 << " of " << trainTrkList.size()
         << "..." << endl;
    cputimer.reset();

    myblood.ComputePriors();

    if (!outDir.empty()) {
      fbase = outDir + '/' + outBase.at(itrk);
    } else {
      fbase = outBase.at(itrk);
    }

    if (!outTestDir.empty()) {
      std::string ftbase;

      ftbase = outTestDir + outBase.at(itrk);
      myblood.WriteOutputs(fbase.c_str(), ftbase.c_str());
    }
    else {
      myblood.WriteOutputs(fbase.c_str());
    }

    cputime = cputimer.milliseconds();
    cout << "Done in " << cputime/1000.0 << " sec." << endl;
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
      outDir = fio_fullpath(pargv[0]);
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
        outBase.push_back(pargv[nargsused]);
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
    else if (!strcmp(option, "--trunc"))
      useTrunc = true;
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
  << "Basic inputs (all must be in common space)" << endl
  << "   --out <base> [...]:" << endl
  << "     Base name(s) of output(s), one per path" << endl
  << "   --outdir <dir>:" << endl
  << "     Output directory (optional)" << endl
  << "   --cptdir <dir>:" << endl
  << "     Output directory for control points in test subject's space" << endl
  << "     (optional, requires registration files)" << endl
  << "     If specified, base names of outputs are relative to this" << endl
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
  << "   --xstr:" << endl
  << "     Exclude previously chosen center streamline(s) (Default: No)" << endl
  << "   --trunc:" << endl
  << "     Also save results using all streamlines, truncated or not" << endl
  << "     (Default: Only save results using non-truncated streamlines)" << endl
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
  if (outBase.empty()) {
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
  if (trainTrkList.size() != outBase.size()) {
    cout << "ERROR: Numbers of input .trk files and output names must match"
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
  if (testMaskList.empty()) {
    cout << "ERROR: Must specify brain mask volume for output subject" << endl;
    exit(1);
  }
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
  if ((!testNonlinXfmFile.empty()) && testNonlinRefFile.empty()) {
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

  if (!outDir.empty()) {
    cout << "Output directory: " << outDir << endl;
  }

  if (!outTestDir.empty()) {
    cout << "Output directory in test subject's space: " << outTestDir << endl;
  }

  cout << "Output base:";
  for (istr = outBase.begin(); istr < outBase.end(); istr++) {
    cout << " " << *istr;
  }
  cout << endl;

  cout << "Training subject directory list: " << trainListFile << endl;

  cout << "Location of streamline files relative to base:";
  for (istr = trainTrkList.begin(); istr < trainTrkList.end(); istr++) {
    cout << " " << *istr;
  }
  cout << endl;

  if (!trainRoi1List.empty()) {
    cout << "Location of start ROI files relative to base:";
    for (istr = trainRoi1List.begin(); istr < trainRoi1List.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;

    cout << "Location of end ROI files relative to base:";
    for (istr = trainRoi2List.begin(); istr < trainRoi2List.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }

  cout << "Location of cortex masks relative to base: " << trainMaskFile
       << endl;

  if (!trainMaskLabel.empty()) {
    cout << "Label ID's from aparc+aseg to add to cortex mask:";
    for (vector<float>::const_iterator ilab = trainMaskLabel.begin();
	 ilab < trainMaskLabel.end(); ilab++) {
      cout << " " << (int) *ilab;
    }
    cout << endl;
  }

  cout << "Location of aparc+aseg's relative to base: " << trainAsegFile
       << endl;

  cout << "Brain mask for output subject:";
  for (auto ifile = testMaskList.begin();
       ifile < testMaskList.end(); ifile++) {
    cout << " " << *ifile;
  }
  cout << endl;

  if (!testFaList.empty()) {
    cout << "FA map for output subject:";
    for (auto ifile = testFaList.begin();
	 ifile < testFaList.end(); ifile++) {
      cout << " " << *ifile;
    }
    cout << endl;
  }

  if (!testAffineXfmFile.empty()) {
    cout << "Affine registration from atlas to base for output subject: "
         << testAffineXfmFile << endl;
  }

  if (!testNonlinXfmFile.empty()) {
    cout << "Nonlinear registration from atlas to base for output subject: "
         << testNonlinXfmFile << endl;
  }

  if (!testNonlinRefFile.empty()) {
    cout << "Nonlinear registration source reference for output subject: "
         << testNonlinRefFile << endl;
  }

  if (!testBaseXfmList.empty()) {
    cout << "Affine registration from base to FA map for output subject:";
    for (auto ifile = testBaseXfmList.begin();
	 ifile < testBaseXfmList.end(); ifile++) {
      cout << " " << *ifile;
    }
    cout << endl;
  }

  if (!testBaseMaskFile.c_str()) {
    cout << "Base mask for output subject: " << testBaseMaskFile << endl;
  }

  cout << "Number of control points for initial spline:";
  for (vector< vector<int> >::const_iterator inlist = nControl.begin();
       inlist < nControl.end(); inlist++) {
    for (vector<int>::const_iterator inum = inlist->begin();
	 inum < inlist->end(); inum++) {
      cout << " " << *inum;
    }
  }
  cout << endl;

  cout << "Exclude previously chosen center streamlines: " << excludeStr << endl
       << "Use truncated streamlines: " << useTrunc << endl;

  return;
}

