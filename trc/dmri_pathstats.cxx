/**
 * @file  dmri_pathstats.c
 * @brief Compute measures on probabilistic or deterministic tractography paths
 *
 * Compute measures on probabilistic or deterministic tractography paths
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2011/02/27 07:42:19 $
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
 *
 */

#include "blood.h"
#include "spline.h"

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
const char *Progname = "dmri_pathstats";

char *inTrkFile = NULL, *inRoi1File = NULL, *inRoi2File = NULL,
     *inTrcDir = NULL, *dtBase = NULL,
     *outFile = NULL, *outVoxFile = NULL, *outStrFile = NULL,
     fname[PATH_MAX];

MRI *l1, *l2, *l3, *v1;

struct utsname uts;
char *cmdline, cwd[2000];

struct timeb cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  vector<MRI *> meas;
  ofstream fout;

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

  printf("Computing statistics on %s...\n", inTrcDir?inTrcDir:inTrkFile);
  TimerStart(&cputimer);

  if (dtBase) {
    sprintf(fname, "%s_L1.nii.gz", dtBase);
    l1 = MRIread(fname);
    sprintf(fname, "%s_L2.nii.gz", dtBase);
    l2 = MRIread(fname);
    sprintf(fname, "%s_L3.nii.gz", dtBase);
    l3 = MRIread(fname);
    sprintf(fname, "%s_V1.nii.gz", dtBase);
    v1 = MRIread(fname);

    sprintf(fname, "%s_L1.nii.gz", dtBase);
    meas.push_back(MRIread(fname));		// Axial diffusivity
    sprintf(fname, "%s_L2.nii.gz", dtBase);
    meas.push_back(MRIread(fname));
    MRIadd(l3, meas[1], meas[1]);
    MRIscalarMul(meas[1], meas[1], .5);		// Radial diffusivity
    sprintf(fname, "%s_MD.nii.gz", dtBase);
    meas.push_back(MRIread(fname));		// Mean diffusivity
    sprintf(fname, "%s_FA.nii.gz", dtBase);
    meas.push_back(MRIread(fname));		// Fractional anisotropy
  }

  if (outFile) {
    fout.open(outFile, ios::out);

    fout << "Count" << "\t"
         << "Volume" << "\t"
         << "Length_Min" << "\t"
         << "Length_Max" << "\t"
         << "Length_Avg" << "\t"
         << "Length_Center";

    if (dtBase)
      fout << "\t"
           << "AD_Avg" << "\t"
           << "AD_Avg_Weight" << "\t"
           << "AD_Avg_Center" << "\t"
           << "RD_Avg" << "\t"
           << "RD_Avg_Weight" << "\t"
           << "RD_Avg_Center" << "\t"
           << "MD_Avg" << "\t"
           << "MD_Avg_Weight" << "\t"
           << "MD_Avg_Center" << "\t"
           << "FA_Avg" << "\t"
           << "FA_Avg_Weight" << "\t"
           << "FA_Avg_Center";

    fout << endl;
  }

  if (outVoxFile) {
    ofstream fvox(outVoxFile, ios::out);

    fvox << "x" << "\t"
         << "y" << "\t"
         << "z" << "\t"
         << "AD" << "\t"
         << "RD" << "\t"
         << "MD" << "\t"
         << "FA" << endl;
  }

  if (inTrcDir > 0) {			// Probabilistic paths
    int len, nx, ny, nz, nvox = 0;
    float wtot = 0, thresh = 0;
    char mname[PATH_MAX];
    vector<int> lengths;
    vector<float> avg(meas.size(), 0), wavg(meas.size(), 0);
    vector<float>::iterator iavg, iwavg;
    MRI *post;
    ifstream infile;

    // Read lengths of path samples
    sprintf(fname, "%s/LENGTH_1_1.txt", inTrcDir);
    infile.open(fname, ios::in);
    if (!infile) {
      cout << "ERROR: Could not open " << fname << " for reading" << endl;
      exit(1);
    }

    while (infile >> len)
      lengths.push_back(len);

    infile.close();

    // Sum path lengths
    len = 0;
    for (vector<int>::const_iterator ilen = lengths.begin();
                                     ilen < lengths.end(); ilen++)
      len += *ilen;

    // Read path posterior distribution
    sprintf(fname, "%s/Fsamples_1_1.nii.gz", inTrcDir);
    post = MRIread(fname);
    nx = post->width;
    ny = post->height;
    nz = post->depth;

    // Find maximum value of posterior distribution
    for (int iz = 0; iz < nz; iz++)
      for (int iy = 0; iy < ny; iy++)
        for (int ix = 0; ix < nx; ix++) {
          const float h = MRIgetVoxVal(post, ix, iy, iz, 0);

          if (h > thresh)
            thresh = h;
        }

    // Set threshold at 20% of maximum
    thresh *= .2;

    // Compute average and weighted average of measures on thresholded posterior
    for (int iz = 0; iz < nz; iz++)
      for (int iy = 0; iy < ny; iy++)
        for (int ix = 0; ix < nx; ix++) {
          const float h = MRIgetVoxVal(post, ix, iy, iz, 0);

          if (h > thresh) {
            iavg = avg.begin();
            iwavg = wavg.begin();

            for (vector<MRI *>::const_iterator ivol = meas.begin();
                                               ivol < meas.end(); ivol++) {
              *iavg += MRIgetVoxVal(*ivol, ix, iy, iz, 0);
              *iwavg += h * MRIgetVoxVal(*ivol, ix, iy, iz, 0);

              iavg++;
              iwavg++;
            }

            nvox++;
            wtot += h;
          }
        }

    for (iavg = avg.begin(); iavg < avg.end(); iavg++)
      *iavg /= nvox;

    for (iwavg = wavg.begin(); iwavg < wavg.end(); iwavg++)
      *iwavg /= wtot;

    // Read control points of MAP path sample
    sprintf(fname, "%s/CONTROLS_1_1.txt", inTrcDir);
    sprintf(mname, "%s/Fsamples_1_1.nii.gz", inTrcDir);
    Spline myspline(fname, mname);
    myspline.InterpolateSpline();

    // Overall measures
    if (outFile) {
      fout << lengths.size() << "\t"
           << nvox << "\t"
           << *min_element(lengths.begin(), lengths.end()) << "\t"
           << *max_element(lengths.begin(), lengths.end()) << "\t"
           << len / (float) lengths.size() << "\t"
           << (myspline.GetAllPointsEnd() - myspline.GetAllPointsBegin()) / 3;

      if (dtBase) {
        vector<float> avgmap = myspline.ComputeAvg(meas);

        for (unsigned int k = 0; k < meas.size(); k++)
          fout << "\t" << avg[k]
               << "\t" << wavg[k]
               << "\t" << avgmap[k];
      }
    }

    // Measures by voxel on MAP streamline
    if (outVoxFile)
      myspline.WriteValues(meas, outVoxFile);
  }
  else {				// Deterministic paths
    // Read .trk file
    Blood myblood(inTrkFile, inRoi1File, inRoi2File);

    myblood.ComputeHistogram();
    myblood.MatchStreamlineEnds();
    myblood.FindCenterStreamline();

    // Overall measures
    if (outFile) {
      fout << myblood.GetNumStr() << "\t"
           << myblood.GetVolume() << "\t"
           << myblood.GetLengthMin() << "\t"
           << myblood.GetLengthMax() << "\t"
           << myblood.GetLengthAvg() << "\t"
           << myblood.GetLengthCenter();

      if (dtBase) {
        vector<float> avgpath   = myblood.ComputeAvgPath(meas),
                      wavgpath  = myblood.ComputeWeightAvgPath(meas),
                      avgcenter = myblood.ComputeAvgCenter(meas);

        for (unsigned int k = 0; k < meas.size(); k++)
          fout << "\t" << avgpath[k]
               << "\t" << wavgpath[k]
               << "\t" << avgcenter[k];
      }
    }

    // Measures by voxel on center streamline
    if (outVoxFile)
      myblood.WriteValuesCenter(meas, outVoxFile);

    // Save center streamline
    if (outStrFile)
      myblood.WriteCenterStreamline(outStrFile, inTrkFile);
  }

  if (outFile)
    fout.close();

  cputime = TimerStop(&cputimer);
  printf("Done in %g sec.\n", cputime/1000.0);

  printf("dmri_pathstats done\n");
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
    else if (!strcmp(option, "--intrk")) {
      if (nargc < 1) CMDargNErr(option,1);
      inTrkFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--rois")) {
      if (nargc < 2) CMDargNErr(option,2);
      inRoi1File = fio_fullpath(pargv[0]);
      inRoi2File = fio_fullpath(pargv[1]);
      nargsused = 2;
    }
    else if (!strcmp(option, "--intrc")) {
      if (nargc < 1) CMDargNErr(option,1);
      inTrcDir = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--dtbase")) {
      if (nargc < 1) CMDargNErr(option,1);
      dtBase = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--out")) {
      if (nargc < 1) CMDargNErr(option,1);
      outFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--outvox")) {
      if (nargc < 1) CMDargNErr(option,1);
      outVoxFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--outstr")) {
      if (nargc < 1) CMDargNErr(option,1);
      outStrFile = fio_fullpath(pargv[0]);
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
  printf("USAGE: ./dmri_pathstats\n");
  printf("\n");
  printf("   --intrk <file>:\n");
  printf("     Input trackvis .trk file\n");
  printf("   --rois <file1> <file2>:\n");
  printf("     Input labeling ROIs for .trk file (optional)\n");
  printf("   --intrc <file>:\n");
  printf("     Input tracula directory\n");
  printf("   --dtbase <file>:\n");
  printf("     Base name of input dtifit files (optional)\n");
  printf("   --out <file>:\n");
  printf("     Output text file for overall path measures\n");
  printf("   --outvox <file>:\n");
  printf("     Output text file for voxel-by-voxel measures along path (optional)\n");
  printf("   --outstr <file>:\n");
  printf("     Output volume of center streamline (optional)\n");
  printf("\n");
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
  if(inTrcDir && inTrkFile) {
    printf("ERROR: cannot specify both .trk file and tracula directory\n");
    exit(1);
  }
  if(!inTrcDir && !inTrkFile) {
    printf("ERROR: must specify input .trk file or tracula directory\n");
    exit(1);
  }
  if(!outFile && !outVoxFile && !outStrFile) {
    printf("ERROR: must specify at least one type of output file\n");
    exit(1);
  }
  if(outVoxFile && !dtBase) {
    printf("ERROR: must specify dtifit base name for voxel-by-voxel output\n");
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

  if (inTrkFile)
    fprintf(fp, "Input .trk file: %s\n", inTrkFile);
  if (inRoi1File) {
    fprintf(fp, "Input end ROI 1: %s\n", inRoi1File);
    fprintf(fp, "Input end ROI 2: %s\n", inRoi2File);
  }
  if (inTrcDir)
    fprintf(fp, "Input tracula directory: %s\n", inTrcDir);
  if (dtBase)
    fprintf(fp, "Input DTI fit base: %s\n", dtBase);
  if (outFile)
    fprintf(fp, "Output file for overall measures: %s\n", outFile);
  if (outVoxFile)
    fprintf(fp, "Output file for voxel-by-voxel measures: %s\n", outVoxFile);
  if (outStrFile)
    fprintf(fp, "Output center streamline volume: %s\n", outStrFile);

  return;
}

