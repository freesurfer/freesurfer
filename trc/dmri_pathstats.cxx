/**
 * @file  dmri_pathstats.cxx
 * @brief Compute measures on probabilistic or deterministic tractography paths
 *
 * Compute measures on probabilistic or deterministic tractography paths
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2013/02/12 06:19:26 $
 *    $Revision: 1.8 $
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
static void WriteHeader(char *OutFile);

int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]) ;

static char vcid[] = "";
const char *Progname = "dmri_pathstats";

char *inTrkFile = NULL, *inRoi1File = NULL, *inRoi2File = NULL,
     *inTrcDir = NULL, *dtBase = NULL,
     *outFile = NULL, *outVoxFile = NULL,
     *outStrFile = NULL, *outEndBase = NULL, *refVolFile = NULL,
     fname[PATH_MAX];

MRI *l1, *l2, *l3, *v1;

struct utsname uts;
char *cmdline, cwd[2000], subjName[100], pathName[100] ;

struct timeb cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime, count, volume, lenmin, lenmax, lencent;
  float lenavg;
  vector<float> avg, wavg, cavg;
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

  if (outVoxFile) {
    WriteHeader(outVoxFile);

    ofstream fvox(outVoxFile, ios::app);
    fvox << "# pathway start" << endl;
    fvox << "x y z AD RD MD FA" << endl;
    fvox.close();
  }

  if (inTrcDir > 0) {			// Probabilistic paths
    int len, nx, ny, nz, nvox = 0;
    float wtot = 0, thresh = 0;
    char mname[PATH_MAX];
    vector<int> lengths;
    vector<float>::iterator iavg, iwavg;
    MRI *post;
    ifstream infile;

    // Read lengths of path samples
    sprintf(fname, "%s/length.samples.txt", inTrcDir);
    infile.open(fname, ios::in);
    if (!infile) {
      cout << "ERROR: Could not open " << fname << " for reading" << endl;
      exit(1);
    }

    // Sum path lengths
    lenavg = 0;
    while (infile >> len) {
      lengths.push_back(len);
      lenavg += len;
    }

    infile.close();

    // Read path posterior distribution
    sprintf(fname, "%s/path.pd.nii.gz", inTrcDir);
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
    avg.resize(meas.size());
    fill(avg.begin(), avg.end(), 0.0);
    wavg.resize(meas.size());
    fill(wavg.begin(), wavg.end(), 0.0);

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

    if (nvox > 0)
      for (iavg = avg.begin(); iavg < avg.end(); iavg++)
        *iavg /= nvox;

    if (wtot > 0)
      for (iwavg = wavg.begin(); iwavg < wavg.end(); iwavg++)
        *iwavg /= wtot;

    // Read control points of MAP path sample
    sprintf(fname, "%s/cpts.map.txt", inTrcDir);
    sprintf(mname, "%s/path.pd.nii.gz", inTrcDir);
    Spline myspline(fname, mname);
    myspline.InterpolateSpline();

    // Overall measures
    count   = lengths.size();
    volume  = nvox;
    lenmin  = *min_element(lengths.begin(), lengths.end());
    lenmax  = *max_element(lengths.begin(), lengths.end());
    lenavg  = ( (lenavg > 0) ? (lenavg / (float) lengths.size()) : 0 );
    lencent = (myspline.GetAllPointsEnd() - myspline.GetAllPointsBegin()) / 3;

    if (dtBase)
      cavg = myspline.ComputeAvg(meas);

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
    count   = myblood.GetNumStr();
    volume  = myblood.GetVolume();
    lenmin  = myblood.GetLengthMin();
    lenmax  = myblood.GetLengthMax();
    lenavg  = myblood.GetLengthAvg();
    lencent = myblood.GetLengthCenter();

    if (dtBase) {
      avg  = myblood.ComputeAvgPath(meas);
      wavg = myblood.ComputeWeightAvgPath(meas);
      cavg = myblood.ComputeAvgCenter(meas);
    }

    // Measures by voxel on center streamline
    if (outVoxFile)
      myblood.WriteValuesCenter(meas, outVoxFile);

    // Save center streamline
    if (outStrFile)
      myblood.WriteCenterStreamline(outStrFile, inTrkFile);

    // Save streamline end points
    if (outEndBase) {
      MRI *refvol;

      if (refVolFile)
        refvol = MRIread(refVolFile);
      else
        refvol = l1;

      myblood.WriteEndPoints(outEndBase, refvol);
    }
  }

  if (outFile) {
    WriteHeader(outFile);

    fout.open(outFile, ios::app);

    fout << "Count "      << count   << endl
         << "Volume "     << volume  << endl
         << "Len_Min "    << lenmin  << endl
         << "Len_Max "    << lenmax  << endl
         << "Len_Avg "    << lenavg  << endl
         << "Len_Center " << lencent << endl;

    if (dtBase)
      fout << "AD_Avg "        << avg[0]  << endl
           << "AD_Avg_Weight " << wavg[0] << endl
           << "AD_Avg_Center " << cavg[0] << endl
           << "RD_Avg "        << avg[1]  << endl
           << "RD_Avg_Weight " << wavg[1] << endl
           << "RD_Avg_Center " << cavg[1] << endl
           << "MD_Avg "        << avg[2]  << endl
           << "MD_Avg_Weight " << wavg[2] << endl
           << "MD_Avg_Center " << cavg[2] << endl
           << "FA_Avg "        << avg[3]  << endl
           << "FA_Avg_Weight " << wavg[3] << endl
           << "FA_Avg_Center " << cavg[3] << endl;

    fout.close();
  }

  if (outVoxFile) {
    ofstream fvox(outVoxFile, ios::app);
    fvox << "# pathway end" << endl;
    fvox.close();
  }

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
    else if (!strcmp(option, "--path")) {
      if (nargc < 1) CMDargNErr(option,1);
      strcpy(pathName, pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--subj")) {
      if (nargc < 1) CMDargNErr(option,1);
      strcpy(subjName, pargv[0]);
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
    else if (!strcmp(option, "--outend")) {
      if (nargc < 1) CMDargNErr(option,1);
      outEndBase = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--ref")) {
      if (nargc < 1) CMDargNErr(option,1);
      refVolFile = fio_fullpath(pargv[0]);
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
  printf("   --path <name>:\n");
  printf("     Name of pathway (optional, written to output files)\n");
  printf("   --subj <name>:\n");
  printf("     Name of subject (optional, written to output files)\n");
  printf("   --out <file>:\n");
  printf("     Output text file for overall path measures\n");
  printf("   --outvox <file>:\n");
  printf("     Output text file for voxel-by-voxel measures along path (optional)\n");
  printf("   --outstr <file>:\n");
  printf("     Output .trk file of center streamline (optional)\n");
  printf("   --outend <base>:\n");
  printf("     Base name of output volumes of streamline ends (optional)\n");
  printf("   --ref <file>:\n");
  printf("     Reference volume (needed only if using --outend without --dtbase)\n");
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
  if(!outFile && !outVoxFile && !outStrFile && !outEndBase) {
    printf("ERROR: must specify at least one type of output\n");
    exit(1);
  }
  if(outVoxFile && !dtBase) {
    printf("ERROR: must specify dtifit base name for voxel-by-voxel output\n");
    exit(1);
  }
  if(outStrFile && !inTrkFile) {
    printf("ERROR: must specify input .trk file to use --outstr\n");
    exit(1);
  }
  if(outEndBase && !inTrkFile) {
    printf("ERROR: must specify input .trk file to use --outend\n");
    exit(1);
  }
  if(outEndBase && !refVolFile && !dtBase) {
    printf("ERROR: must specify reference volume to use --outend\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void WriteHeader(char *OutFile) {
  ofstream fout(OutFile, ios::out);

  fout << "# Title Pathway Statistics" << endl
       << "#"                          << endl
       << "# generating_program "      << Progname << endl
       << "# cvs_version "             << vcid << endl
       << "# cmdline "                 << cmdline << endl
       << "# sysname "                 << uts.sysname << endl
       << "# hostname "                << uts.nodename << endl
       << "# machine "                 << uts.machine << endl
       << "# user "                    << VERuser() << endl
       << "# anatomy_type pathway"     << endl
       << "#"                          << endl
       << "# subjectname "             << subjName << endl
       << "# pathwayname "             << pathName << endl
       << "#"                          << endl;

  fout.close();
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
//  if (pathName)
    fprintf(fp, "Pathway name: %s\n", pathName);
//  if (subjName)
    fprintf(fp, "Subject name: %s\n", subjName);
  if (outFile)
    fprintf(fp, "Output file for overall measures: %s\n", outFile);
  if (outVoxFile)
    fprintf(fp, "Output file for voxel-by-voxel measures: %s\n", outVoxFile);
  if (outStrFile)
    fprintf(fp, "Output center streamline .trk file: %s\n", outStrFile);
  if (outEndBase)
    fprintf(fp, "Base name of output end point volumes: %s\n", outEndBase);
  if (refVolFile)
    fprintf(fp, "Reference for output end point volumes: %s\n", refVolFile);

  return;
}

