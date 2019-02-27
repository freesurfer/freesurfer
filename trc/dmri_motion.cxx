/**
 * @file  dmri_motion.cxx
 * @brief Compute measures of head motion in DWIs
 *
 * Compute measures of head motion in DWIs
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2014/05/27 14:49:34 $
 *    $Revision: 1.6 $
 *
 * Copyright © 2013 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <limits.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <set>
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
static void dump_options();

int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]) ;

static char vcid[] = "";
const char *Progname = "dmri_motion";

float T = 100, D = .001;

char *inMatFile = NULL, *inDwiFile = NULL, *inBvalFile = NULL,
     *outFile = NULL, *outFrameFile = NULL;

MRI *dwi;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  float travg = 0, roavg = 0, score = 0, pbad = 0;
  vector<int> nbadframe;
  vector<float> trframe, roframe, scoreframe;
  string matline;
  ofstream outfile;

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

  dump_options();

  cputimer.reset();

  if (inMatFile) {		// Estimate between-volume motion
    int nframe = 0;
    vector<float> xform, tr0(3,0), ro0(3,0), trtot(3,0), rotot(3,0);
    ifstream infile;

    cout << "Loading volume-to-baseline affine transformations" << endl;
    infile.open(inMatFile, ios::in);
    if (!infile) {
      cout << "ERROR: Could not open " << inMatFile << " for reading" << endl;
      exit(1);
    }

    cout << "Computing between-volume head motion measures" << endl;

    while (getline(infile, matline))
      if (~matline.empty() && ~isalpha(matline[0])) {
        float xval;
        istringstream matstr(matline);

        while (matstr >> xval)
          xform.push_back(xval);

        if (xform.size() == 16) {
          vector<float> tr, ro;
          AffineReg reg(xform);

          // Decompose affine registration matrix into its parameters
          reg.DecomposeXfm();

          // Translations with respect to first frame
          tr = reg.GetTranslate();

          // Frame-to-frame translations
          for (int k = 0; k < 3; k++)
            trframe.push_back(tr[k] - tr0[k]);

          // Cumulative frame-to-frame translations
          for (int k = 0; k < 3; k++)
            trtot[k] += fabs(*(trframe.end()-3+k));

          copy(tr.begin(), tr.end(), tr0.begin());

          // Rotations with respect to first frame
          ro = reg.GetRotate();

          // Frame-to-frame rotations
          for (int k = 0; k < 3; k++)
            roframe.push_back(ro[k] - ro0[k]);

          // Cumulative frame-to-frame rotations
          for (int k = 0; k < 3; k++)
            rotot[k] += fabs(*(roframe.end()-3+k));

          copy(ro.begin(), ro.end(), ro0.begin());

          xform.clear();

          nframe++;
        }
      }

    infile.close();

    cout << "INFO: Processed transforms for " << nframe << " volumes" << endl;

    travg = sqrt(trtot[0]*trtot[0] + trtot[1]*trtot[1] + trtot[2]*trtot[2])
          / nframe;
    roavg = (rotot[0] + rotot[1] + rotot[2]) / nframe;
  }

  if (inBvalFile) {		// Estimate within-volume motion
    int nx, ny, nz, nd, nxy, nslice = 0, nbad = 0;
    float minvox, b;
    vector<int> r, r1;
    vector<int>::const_iterator ir1;
    vector<float> bvals;
    vector<float>::const_iterator ibval;
    ifstream infile;
 
    // Read DWI volume series
    cout << "Loading DWI volume series from " << inDwiFile << endl;
    dwi = MRIread(inDwiFile);
    if (!dwi) {
      cout << "ERROR: Could not read " << inDwiFile << endl;
      exit(1);
    }

    nx = dwi->width;
    ny = dwi->height;
    nz = dwi->depth;
    nd = dwi->nframes;

    nxy = nx*ny;
    minvox = 0.05*nxy;

    nbadframe.resize(nd);
    fill(nbadframe.begin(), nbadframe.end(), 0);
    scoreframe.resize(nd);
    fill(scoreframe.begin(), scoreframe.end(), 0.0);

    // Read b-value table
    cout << "Loading b-value table from " << inBvalFile << endl;
    infile.open(inBvalFile, ios::in);
    if (!infile) {
      cout << "ERROR: Could not open " << inBvalFile << " for reading" << endl;
      exit(1);
    }

    while (infile >> b)
      bvals.push_back(b);

    infile.close();

    cout << "Computing within-volume head motion measures" << endl;

    // Find unique b-values
    set<float> blist(bvals.begin(), bvals.end());

    // Compare frames acquired with each b-value separately
    for (set<float>::const_iterator ib = blist.begin(); ib != blist.end();
                                                        ib++) {
      const float thresh = T * exp(-(*ib)*D);
      vector<int>::iterator inbad = nbadframe.begin();
      vector<float>::iterator iscore = scoreframe.begin();

      r1.clear();
      ibval = bvals.begin();

      for (int id = 0; id < nd; id++) {
        if (*ibval == *ib) {
          r.clear();

          // Find number of voxels above threshold for each slice in this frame
          for (int iz = 0; iz < nz; iz++) {
            int count = 0;

            for (int iy = 0; iy < ny; iy++)
              for (int ix = 0; ix < nx; ix++)
                if (MRIgetVoxVal(dwi, ix, iy, iz, id) > thresh)
                  count++;

            r.push_back(count);
          }

          if (r1.empty())		// First frame with this b-value
            r1.insert(r1.begin(), r.begin(), r.end());

          // Motion score (from Benner et al MRM 2011)
          ir1 = r1.begin();
          for (vector<int>::const_iterator ir = r.begin(); ir < r.end(); ir++) {
            if (*ir >= minvox) {	// Do not count empty slices
              const float S = 2 - *ir / (0.7 * (*ir1));

              nslice++;

              if (S > 1) {
                (*inbad)++;
                (*iscore) += S;
              }
            }

            ir1++;
          }

          nbad += *inbad;
          score += *iscore;

          // Average motion score of bad slices in this frame
          if (*inbad > 0)
            *iscore /= *inbad;
          else
            *iscore = 1;
        }

        ibval++;
        inbad++;
        iscore++;
      }
    }

    // Percentage of bad slices among all non-empty slices
    if (nslice > 0)
      pbad = nbad / (float) nslice * 100;

    // Average motion score of bad slices
    if (nbad > 0)
      score /= nbad;
    else
      score = 1;
  }

  // Write overall measures to file
  outfile.open(outFile, ios::out);
  outfile << "AvgTranslation AvgRotation PercentBadSlices AvgDropoutScore"
          << endl
          << travg << " " << roavg << " " << pbad << " " << score << endl;
  outfile.close();

  // Write frame-by-frame measures to file
  if (outFrameFile) {
    vector<float>::const_iterator itr = trframe.begin(),
                                  iro = roframe.begin(),
                                  iscore = scoreframe.begin();

    if (trframe.empty()) {
      trframe.resize(nbadframe.size() * 3);
      fill(trframe.begin(), trframe.end(), 0.0);
    }
    else if (trframe.size() != nbadframe.size() * 3) {
      cout << "ERROR: inconsistent number of frames for "
           << "between-volume (" << trframe.size()/3 << ") and "
           << "within-volume (" << nbadframe.size() << ") measures" << endl;
      exit(1);
    }

    if (roframe.empty()) {
      roframe.resize(nbadframe.size() * 3);
      fill(roframe.begin(), roframe.end(), 0.0);
    }
    else if (roframe.size() != nbadframe.size() * 3) {
      cout << "ERROR: inconsistent number of frames for "
           << "between-volume (" << roframe.size()/3 << ") and "
           << "within-volume (" << nbadframe.size() << ") measures" << endl;
      exit(1);
    }

    if (nbadframe.empty()) {
      nbadframe.resize(trframe.size() / 3);
      fill(nbadframe.begin(), nbadframe.end(), 0);
    }
    else if (nbadframe.size() != trframe.size() / 3) {
      cout << "ERROR: inconsistent number of frames for "
           << "between-volume (" << trframe.size()/3 << ") and "
           << "within-volume (" << nbadframe.size() << ") measures" << endl;
      exit(1);
    }

    if (scoreframe.empty()) {
      scoreframe.resize(trframe.size() / 3);
      fill(scoreframe.begin(), scoreframe.end(), 0.0);
    }
    else if (scoreframe.size() != trframe.size() / 3) {
      cout << "ERROR: inconsistent number of frames for "
           << "between-volume (" << trframe.size()/3 << ") and "
           << "within-volume (" << scoreframe.size() << ") measures" << endl;
      exit(1);
    }

    outfile.open(outFrameFile, ios::out);
    outfile << "TranslationX TranslationY TranslationZ "
            << "RotationX RotationY RotationZ "
            << "PercentBadSlices AvgDropoutScore" << endl;

    for (vector<int>::const_iterator inbad = nbadframe.begin();
                                     inbad < nbadframe.end(); inbad++) {
    
      outfile << itr[0] << " " << itr[1] << " " << itr[2] << " "
              << iro[0] << " " << iro[1] << " " << iro[2] << " "
              << *inbad << " " << *iscore << endl;

      itr += 3;
      iro += 3;
      iscore++;
    }

    outfile.close();
  }

  cputime = cputimer.milliseconds();
  cout << "Done in " << cputime/1000.0 << " sec." << endl;

  cout << "dmri_motion done" << endl;
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
    else if (!strcmp(option, "--dwi")) {
      if (nargc < 1) CMDargNErr(option,1);
      inDwiFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--bval")) {
      if (nargc < 1) CMDargNErr(option,1);
      inBvalFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--mat")) {
      if (nargc < 1) CMDargNErr(option,1);
      inMatFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--T")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%f", &T);
      nargsused = 1;
    }
    else if (!strcmp(option, "--D")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%f", &D);
      nargsused = 1;
    }
    else if (!strcmp(option, "--out")) {
      if (nargc < 1) CMDargNErr(option,1);
      outFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--outf")) {
      if (nargc < 1) CMDargNErr(option,1);
      outFrameFile = fio_fullpath(pargv[0]);
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
  cout
  << endl << "USAGE: " << Progname << endl << endl
  << "Required arguments" << endl
  << "   --out <file>:" << endl
  << "     Output text file of motion measures" << endl
  << endl
  << "Optional arguments" << endl
  << "   --outf <file>:" << endl
  << "     Output text file of frame-by-frame motion measures" << endl
  << endl
  << "Arguments needed for between-volume motion measures" << endl
  << "   --mat <file>:" << endl
  << "     Input text file of volume-to-baseline affine transformations" << endl
  << endl
  << "Arguments needed for within-volume motion measures" << endl
  << "(see Benner et al MRM 2011):" << endl
  << "   --dwi <file>:" << endl
  << "     Input DWI volume series, unprocessed" << endl
  << "   --bval <file>:" << endl
  << "     Input b-value table" << endl
  << "   --T <num>:" << endl
  << "     Low-b image intensity threshold (default: 100)" << endl
  << "   --D <num>:" << endl
  << "     Nominal diffusivity (default: .001)" << endl
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
  cout << vcid << endl;
  exit(1);
}

/* --------------------------------------------- */
static void check_options(void) {
  if (!outFile) {
    cout << "ERROR: must specify output file" << endl;
    exit(1);
  }
  if ((inDwiFile && !inBvalFile) || (!inDwiFile && inBvalFile)) {
    cout << "ERROR: must specify both DWI and b-value files" << endl;
    exit(1);
  }
  if (!inBvalFile && !inMatFile) {
    cout << "ERROR: must specify inputs for between-volume and/or "
         << "within-volume motion measures" << endl;
    exit(1);
  }
  if (D < 0) {
    cout << "ERROR: diffusivity must be positive" << endl;
    exit(1);
  }
  return;
}

static void dump_options() {
  cout << endl
       << vcid << endl
       << "cwd " << cwd << endl
       << "cmdline " << cmdline << endl
       << "sysname  " << uts.sysname << endl
       << "hostname " << uts.nodename << endl
       << "machine  " << uts.machine << endl
       << "user     " << VERuser() << endl;

  cout << "Output motion measure file: " << outFile << endl;

  if (outFrameFile)
    cout << "Output frame-by-frame motion measure file: " << outFrameFile
         << endl;

  if (inMatFile)
    cout << "Input transform file: " << inMatFile << endl;

  if (inBvalFile) {
    cout << "Input DWI file: " << inDwiFile << endl;
    cout << "Input b-value table: " << inBvalFile << endl;
    cout << "Low-b image intensity threshold: " << T << endl;
    cout << "Nominal diffusivity: " << D << endl;
  }

  return;
}

