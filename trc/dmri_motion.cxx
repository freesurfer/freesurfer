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
 *    $Date: 2013/10/25 03:12:59 $
 *    $Revision: 1.1 $
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
static void dump_options(FILE *fp);

int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]) ;

static char vcid[] = "";
const char *Progname = "dmri_motion";

float T = 100, D = .001;

char *inDwiFile = NULL, *inMatFile = NULL, *inBvalFile = NULL,
     *outFile = NULL;

MRI *dwi;

struct utsname uts;
char *cmdline, cwd[2000];

struct timeb cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime, nx, ny, nz, nd;
  float score = 0, pbad = 0;
  vector<unsigned int> lowbs;
  vector<float> trtot(3,0), rotot(3,0), bvals;
  string matline;
  ifstream infile;
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

  dump_options(stdout);

  TimerStart(&cputimer);

  // Read DWIs
  cout << "Reading DWI volume series..." << endl;
  dwi = MRIread(inDwiFile);
  nx = dwi->width;
  ny = dwi->height;
  nz = dwi->depth;
  nd = dwi->nframes;

  if (inMatFile) {		// Estimate between-volume motion
    int nframe = 0;
    vector<float> xform, tr0(3,0), ro0(3,0);

    cout << "Computing between-volume head motion measures..." << endl;

    // Read volume-to-baseline affine transformations
    infile.open(inMatFile, ios::in);
    if (!infile) {
      cout << "ERROR: Could not open " << inMatFile << " for reading" << endl;
      exit(1);
    }

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

          // Cumulative frame-to-frame translations
          for (int k = 0; k < 3; k++)
            trtot[k] += fabs(tr[k] - tr0[k]);

          copy(tr.begin(), tr.end(), tr0.begin());

          // Rotations with respect to first frame
          ro = reg.GetRotate();

          // Cumulative frame-to-frame rotations
          for (int k = 0; k < 3; k++)
            rotot[k] += fabs(ro[k] - ro0[k]);

          copy(ro.begin(), ro.end(), ro0.begin());

          xform.clear();

          nframe++;
        }
      }

    infile.close();

    if (nframe != nd) {
      cout << "ERROR: Number of volumes in " << inDwiFile
           << " and in " << inMatFile << " do not match" << endl;
      exit(1);
    }
  }

  if (inBvalFile) {		// Estimate within-volume motion
    const int nxy = nx*ny;
    int nslice = 0, nbad = 0;
    float b;
    const float minvox = 0.05*nxy;
    vector<int> r, r1;
    vector<int>::const_iterator ir1;
    vector<float>::const_iterator ibval;
 
    cout << "Computing within-volume head motion measures..." << endl;

    // Read b-value table
    infile.open(inBvalFile, ios::in);
    if (!infile) {
      cout << "ERROR: Could not open " << inBvalFile << " for reading" << endl;
      exit(1);
    }

    while (infile >> b)
      bvals.push_back(b);

    infile.close();

    // Find unique b-values
    set<float> blist(bvals.begin(), bvals.end());

    // Compare frames acquired with each b-value separately
    for (set<float>::const_iterator ib = blist.begin(); ib != blist.end();
                                                        ib++) {
      const float thresh = T * exp(-(*ib)*D);
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
                nbad++;
                score += S;
              }
            }

            ir1++;
          }
        }

        ibval++;
      }
    }

    // Percentage of bad slices among all non-empty slices
    pbad = nbad / (float) nslice * 100;

    // Average motion score of bad slices
    if (nbad > 0)
      score /= nbad;
    else
      score = 1;
  }

  outfile.open(outFile, ios::out);
  outfile << "AvgTranslation AvgRotation PercentBadSlices AvgDropoutScore"
          << endl
          << sqrt(trtot[0]*trtot[0] + trtot[1]*trtot[1] + trtot[2]*trtot[2])/nd
          << " "
          << (rotot[0] + rotot[1] + rotot[2])/nd
          << " "
          << pbad
          << " "
          << score << endl;
  outfile.close();

  cputime = TimerStop(&cputimer);
  printf("Done in %g sec.\n", cputime/1000.0);

  printf("dmri_motion done\n");
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
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}

/* --------------------------------------------- */
static void print_usage(void) 
{
  printf("\n");
  printf("USAGE: ./dmri_motion\n");
  printf("\n");
  printf("   --dwi <file>:\n");
  printf("     Input DWI volume series (unprocessed)\n");
  printf("   --bval <file>:\n");
  printf("     Input b-value table\n");
  printf("   --mat <file>:\n");
  printf("     Input text file of volume-to-baseline affine transformations\n");
  printf("   --out <file>:\n");
  printf("     Output text file of motion measures\n");
  printf("\n");
  printf("Optional arguments (see Benner et al MRM 2011):\n");
  printf("   --T <num>: Low-b image intensity threshold (default: 100)\n");
  printf("   --D <num>: Nominal diffusivity (default: .001)\n");
  printf("\n");
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
  if (!inDwiFile) {
    printf("ERROR: must specify input DWI file\n");
    exit(1);
  }
  if (!inBvalFile && !inMatFile) {
    printf("ERROR: must specify b-value file and/or affine transform file\n");
    exit(1);
  }
  if (!outFile) {
    printf("ERROR: must specify output file\n");
    exit(1);
  }
  if (D < 0) {
    printf("ERROR: diffusivity must be positive\n");
    exit(1);
  }
  return;
}

static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  fprintf(fp, "Input DWI file: %s\n", inDwiFile);
  fprintf(fp, "Input b-value table: %s\n", inBvalFile);
  fprintf(fp, "Input transform file: %s\n", inMatFile);
  fprintf(fp, "Output motion measure file: %s\n", outFile);

  return;
}

