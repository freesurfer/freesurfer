/**
 * @file  dmri_group.cxx
 * @brief Combine path measures from multiple subjects
 *
 * Combine path measures from multiple subjects
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2013/11/08 20:11:46 $
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
static void print_usage(void);
static void usage_exit(void);
static void print_help(void);
static void print_version(void);
static void dump_options(FILE *fp);

int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]);

static char vcid[] = "";
const char *Progname = "dmri_group";

int nSection = 0;

char *inListFile = NULL, *outRefFile = NULL, *outBase = NULL;

struct utsname uts;
char *cmdline, cwd[2000];

struct timeb cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  unsigned int nmeas, npt, lenmax, ntot;
  float distmin, darc,
        arcmin = numeric_limits<float>::infinity(),
        arcmax = 0.0;
  string listline, filename;
  vector<float>::const_iterator icenter;
  vector< vector<float> >::const_iterator itemplate, iallm;
  vector< vector<float> >::iterator ialla;
  vector< vector<unsigned int> >::const_iterator iallk;
  vector<unsigned int> lengths;
  vector<string> subjlist, measlist;
  vector< vector<unsigned int> > allknots;
  vector< vector<float> > allarc, allpaths, allmeas, allmeasint, allmeassec;
  ifstream listfile;
  ofstream pathfile;
  MATRIX *outv2r;
  MRI *outref = 0;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
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

  dump_options(stdout);

  TimerStart(&cputimer);

  // Read output reference volume
  if (outRefFile) {
    cout << "Loading output reference volume from " << outRefFile << endl;
    outref = MRIread(outRefFile);
    if (!outref) {
      cout << "ERROR: Could not read " << outRefFile << endl;
      exit(1);
    }

    outv2r = MRIgetVoxelToRasXform(outref);
  }

  // Read list of inputs
  cout << "Loading list of inputs from " << inListFile << endl;
  listfile.open(inListFile, ios::in);
  if (!listfile) {
    cout << "ERROR: Could not open " << inListFile << " for reading" << endl;
    exit(1);
  }

  while (getline(listfile, listline)) {
    string measline, subjid;
    vector<unsigned int> knots;
    vector<float> arc, path, meas;
    vector<string> inputs;
    ifstream infile;
    istringstream liststr(listline);
    MRI *inref = 0;
    AffineReg affinereg;
#ifndef NO_CVS_UP_IN_HERE
    NonlinReg nonlinreg;
#endif

    while (liststr >> filename)
      inputs.push_back(filename);

    if (inputs.empty())
      continue;

    // First input on each line is the path directory
    filename = inputs[0] + "/pathstats.byvoxel.txt";

    cout << "Loading measures along the path from " << filename << endl;
    infile.open(filename.c_str(), ios::in);
    if (!infile) {
      cout << "WARN: Could not open " << filename << endl
           << "WARN: Skipping ahead" << endl;
      continue;
    }

    // Second input on each line is an input reference volume
    if (inputs.size() > 1) {
      // Read reference volumes
      cout << "Loading input reference volume from " << inputs[1] << endl;
      inref = MRIread(inputs[1].c_str());
      if (!inref) {
        cout << "ERROR: Could not read " << inputs[1] << endl;
        exit(1);
      }
    }

    // Third input on each line is an affine transform
    // Fourth input on each line is a nonlinear transform
#ifndef NO_CVS_UP_IN_HERE
    if (inputs.size() > 3) {
      affinereg.ReadXfm(inputs[2].c_str(), inref, 0);
      nonlinreg.ReadXfm(inputs[3].c_str(), outref);
    }
    else
#endif
    if (inputs.size() > 2)
      affinereg.ReadXfm(inputs[2].c_str(), inref, outref);

    // Read measures along the path
    while (getline(infile, measline)) {
      if (measline.substr(0,1).compare("#") == 0) {		// Comment line
        string word;
        istringstream linestr(measline);

        linestr >> word;
        linestr >> word;

        if (word.compare("subjectname") == 0)
          linestr >> subjid;
      }
      else if (measline.substr(0,1).compare("x") == 0) {	// Header line
        string word;
        istringstream linestr(measline);

        // The first three column headers are "x", "y", "z"
        linestr >> word;
        linestr >> word;
        linestr >> word;

        // The remaining column headers are the names of the measures
        measlist.clear();
        while (linestr >> word)
          measlist.push_back(word);
      }
      else {							// Value line
        float val;
        vector<float> point, valpoint;
        istringstream linestr(measline);

        // The first three columns are the x, y, z coordinates of this point
        linestr >> val;
        point.push_back(val);
        linestr >> val;
        point.push_back(val);
        linestr >> val;
        point.push_back(val);

        // Apply affine transform
        if (!affinereg.IsEmpty())
          affinereg.ApplyXfm(point, point.begin());

#ifndef NO_CVS_UP_IN_HERE
        // Apply nonlinear transform
        if (!nonlinreg.IsEmpty())
          nonlinreg.ApplyXfm(point, point.begin());
#endif

        // Add length of new path segment to arc lengths
        if (path.empty())
          arc.push_back(0.0);
        else {
          const float dx = *(path.end()-3) - point[0],
                      dy = *(path.end()-2) - point[1],
                      dz = *(path.end()-1) - point[2],
                      seglen = sqrt(dx*dx + dy*dy + dz*dz);

          arc.push_back(*(arc.end()-1) + seglen);
        }

        // Coordinates of new path point
        path.insert(path.end(), point.begin(), point.end());

        // The remaining columns are the values of the measures at this point
        while (linestr >> val)
          meas.push_back(val);
      }
    }

    infile.close();

    if (inref)
      MRIfree(&inref);

    if (arc.empty())
      continue;

    subjlist.push_back(subjid);
    lengths.push_back(arc.size());
    allarc.push_back(arc);
    allpaths.push_back(path);
    allmeas.push_back(meas);

    // Start point of the path
    knots.push_back(0);

    // Points that are a quarter of the way or half way along the path
    for (int k = 1; k < 4; k++) {
      const float arcpart = *(arc.end()-1) * 0.25 * k;
      float dmin = numeric_limits<float>::infinity();
      vector<float>::const_iterator imin = arc.begin();

      for (vector<float>::const_iterator iarc = arc.begin(); iarc < arc.end();
                                                             iarc++) {
        const float darc = fabs(arcpart - *iarc);

	if (darc < dmin) {
          imin = iarc;
          dmin = darc;
        }
        else
          break;
      }

      knots.push_back(imin - arc.begin());
    }

    // End point of the path
    knots.push_back(arc.size() - 1);

    allknots.push_back(knots);
  }

  nmeas = measlist.size();

  // Find the most representative path
  iallk = allknots.begin();
  itemplate = allpaths.begin();
  icenter = itemplate->begin() + (allknots[0][2] * 3);
  distmin = numeric_limits<float>::infinity();

  for (vector< vector<float> >::const_iterator iallp = allpaths.begin();
                                               iallp < allpaths.end();
                                               iallp++) {
    float dist = 0.0;
    vector< vector<unsigned int> >::const_iterator jallk = allknots.begin();

    for (vector< vector<float> >::const_iterator jallp = allpaths.begin();
                                                 jallp < allpaths.end();
                                                 jallp++) {
      if (jallp != iallp) {
        vector<unsigned int >::const_iterator iknot = iallk->begin();

        for (vector<unsigned int>::const_iterator jknot = jallk->begin();
                                                  jknot < jallk->end();
                                                  jknot++) {
          const unsigned int ioffset = (*iknot) * 3,
                             joffset = (*jknot) * 3;
          const float dx = iallp->at(ioffset)     - jallp->at(joffset),
                      dy = iallp->at(ioffset + 1) - jallp->at(joffset + 1),
                      dz = iallp->at(ioffset + 2) - jallp->at(joffset + 2);

          dist += sqrt(dx*dx + dy*dy + dz*dz);

          iknot++;
        }
      }

      jallk++;
    }

    if (dist < distmin) {
      itemplate = iallp;			// Most representative path
      icenter = iallp->begin() + (iallk->at(2) * 3);	// Mid-point on path
      distmin = dist;
    }

    iallk++;
  }

  // Write points of most representative path to file as RAS coords
  filename = string(outBase) + ".median.txt";

  cout << "Writing median path to " << filename << endl;
  pathfile.open(filename.c_str(), ios::out);

  pathfile << "#!ascii label" << endl
           << itemplate->size() / 3 << endl;

  npt = 1;

  for (vector<float>::const_iterator ipt = itemplate->begin();
                                     ipt < itemplate->end(); ipt += 3) {
    pathfile << npt;

    for (int k = 1; k < 4; k++)		// Transform from voxel to RAS coords
      pathfile << " " << ipt[0] * outv2r->rptr[k][1] +
                         ipt[1] * outv2r->rptr[k][2] +
                         ipt[2] * outv2r->rptr[k][3] +
                                  outv2r->rptr[k][4];

    pathfile << " 0" << endl;

    npt++;
  }

  pathfile.close();

  // Reparameterize the arc length on each path
  ialla = allarc.begin();

  for (vector< vector<float> >::const_iterator iallp = allpaths.begin();
                                               iallp < allpaths.end();
                                               iallp++) {
    float dmin = numeric_limits<float>::infinity(),
          arc0 = 0, arcm;
    vector<float>::const_iterator iarc = ialla->begin();

    // Find the closest point to the mid-point of the most representative path
    for (vector<float>::const_iterator ipath = iallp->begin();
                                       ipath < iallp->end(); ipath += 3) {
      const float dx = ipath[0] - icenter[0],
                  dy = ipath[1] - icenter[1],
                  dz = ipath[2] - icenter[2],
                  dist = dx*dx + dy*dy + dz*dz;
//...

      if (dist < dmin) {
        arc0 = *iarc;
        dmin = dist;
      }

      iarc++;
    }

    // Make this point the origin of the arc length for this path
    for (vector<float>::iterator iarcnew = ialla->begin();
                                 iarcnew < ialla->end(); iarcnew++)
      *iarcnew -= arc0;

    arcm = *min_element(ialla->begin(), ialla->end());

    if (arcm < arcmin)
      arcmin = arcm;

    arcm = *max_element(ialla->begin(), ialla->end());

    if (arcm > arcmax)
      arcmax = arcm;

    ialla++;
  }

  // Interpolate measures at the same arc lengths on every path
  lenmax = *max_element(lengths.begin(), lengths.end());
  darc = (arcmax - arcmin) / lenmax;

  iallm = allmeas.begin();

  for (vector< vector<float> >::const_iterator ialla = allarc.begin();
                                               ialla < allarc.end(); ialla++) {
    vector<float> meas;

    for (float larc = arcmin + darc; larc <= arcmax; larc += darc) {
      float slope;
      vector<float>::const_iterator iarc = ialla->begin(),
                                    imeas1, imeas0;

      if (*iarc > larc) {		// No points in this segment, skip ahead
        for (int k = (int) nmeas; k > 0; k--)
          meas.push_back(numeric_limits<float>::infinity());

        continue;
      }

      while (*iarc < larc && iarc < ialla->end())
        iarc++;

      if (iarc == ialla->end()) {	// No points in this segment, skip ahead
        for (int k = (int) nmeas; k > 0; k--)
          meas.push_back(numeric_limits<float>::infinity());

        continue;
      }

      // Linear interpolation
      slope = (larc - *(iarc-1)) / (*iarc - *(iarc-1));

      imeas1 = iallm->begin() + nmeas * (iarc - ialla->begin());
      imeas0 = imeas1 - nmeas;

      // Interpolate values of each measure
      for (int k = (int) nmeas; k > 0; k--) {
        meas.push_back(*imeas0 + (*imeas1 - *imeas0) * slope);

        imeas1++;
        imeas0++;
      }
    }

    allmeasint.push_back(meas);

    iallm++;
  }

  // Write output files
  ntot = allmeasint[0].size();

  for (vector<string>::const_iterator imeas = measlist.begin();
                                      imeas < measlist.end(); imeas++) {
    string outname = string(outBase) + "." + *imeas + ".txt";
    ofstream outfile;

    cout << "Writing group table to " << outname << endl;
    outfile.open(outname.c_str(), ios::out);

    // Write subject names
    for (vector<string>::const_iterator isubj = subjlist.begin();
                                        isubj < subjlist.end(); isubj++)
      outfile << *isubj << " ";

    outfile << endl;

    // Write interpolated values of this measure
    for (unsigned ipt = imeas - measlist.begin(); ipt < ntot; ipt += nmeas) {
      for (iallm = allmeasint.begin(); iallm < allmeasint.end(); iallm++) {
        vector<float>::const_iterator ival = iallm->begin() + ipt;

        if (*ival == numeric_limits<float>::infinity())
          outfile << "NaN ";
        else
          outfile << *ival << " ";
      }

      outfile << endl;
    }

    outfile.close();
  }

  // Average measures over sections along the path
  if (nSection > 0) {
    char nsec[PATH_MAX];

    sprintf(nsec, "%dsec", nSection);

    darc = (arcmax - arcmin) / nSection;

    iallm = allmeas.begin();

    for (vector< vector<float> >::const_iterator ialla = allarc.begin();
                                                 ialla < allarc.end();
                                                 ialla++) {
      float larc = arcmin + darc;
      vector<float>::const_iterator iarc = ialla->begin(),
                                    imeas = iallm->begin();
      vector<float> meas;

      while (larc <= arcmax) {
        int nsamp = 0;
        vector<float> avg(nmeas,0);

        while (*iarc < larc && iarc < ialla->end()) {
          for (vector<float>::iterator iavg = avg.begin(); iavg < avg.end();
                                                           iavg++) {
            *iavg += *imeas;

            imeas++;
          }

          iarc++;
          nsamp++;
        }

        if (nsamp > 0)
          for (vector<float>::iterator iavg = avg.begin(); iavg < avg.end();
                                                            iavg++)
            *iavg /= nsamp;
        else					// No points in this section
          for (vector<float>::iterator iavg = avg.begin(); iavg < avg.end();
                                                           iavg++)
            *iavg = numeric_limits<float>::infinity();

        meas.insert(meas.end(), avg.begin(), avg.end());

        larc += darc;
      }

      allmeassec.push_back(meas);

      iallm++;
    }

    // Write output files
    ntot = nSection * nmeas;

    for (vector<string>::const_iterator imeas = measlist.begin();
                                        imeas < measlist.end(); imeas++) {
      string outname = string(outBase) + "." + *imeas + "." + nsec + ".txt";
      ofstream outfile;

      cout << "Writing group table to " << outname << endl;
      outfile.open(outname.c_str(), ios::out);

      // Write subject names
      for (vector<string>::const_iterator isubj = subjlist.begin();
                                          isubj < subjlist.end(); isubj++)
        outfile << *isubj << " ";

      outfile << endl;

      // Write section averages of values of this measure
      for (unsigned ipt = imeas - measlist.begin(); ipt < ntot; ipt += nmeas) {
        for (iallm = allmeassec.begin(); iallm < allmeassec.end(); iallm++) {
          vector<float>::const_iterator ival = iallm->begin() + ipt;

          if (*ival == numeric_limits<float>::infinity())
            outfile << "NaN ";
          else
            outfile << *ival << " ";
        }

        outfile << endl;
      }

      outfile.close();
    }
  }

  if (outref) {
    MRIfree(&outref);
    MatrixFree(&outv2r);
  }

  cputime = TimerStop(&cputimer);
  cout << "Done in " << cputime/1000.0 << " sec." << endl;

  cout << "dmri_group done" << endl;
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
    else if (!strcmp(option, "--list")) {
      if (nargc < 1) CMDargNErr(option,1);
      inListFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--ref")) {
      if (nargc < 1) CMDargNErr(option,1);
      outRefFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--out")) {
      if (nargc < 1) CMDargNErr(option,1);
      outBase = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--sec")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &nSection);
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
  cout
  << endl << "USAGE: " << Progname << endl << endl
  << "Required arguments" << endl
  << "   --list <file>:" << endl
  << "     Text file with list of individual inputs" << endl
  << "   --ref <file>:" << endl
  << "     Reference volume for output path" << endl
  << "   --out <base>:" << endl
  << "     Base name of output text files" << endl
  << "   --sec <num>:" << endl
  << "     Divide the pathway into a number of sections and output " << endl
  << "     average measures for each section (optional)" << endl
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
  if (!outBase) {
    cout << "ERROR: must specify base name for output files" << endl;
    exit(1);
  }
  if (!inListFile) {
    cout << "ERROR: must specify input list file" << endl;
    exit(1);
  }
  return;
}

static void dump_options(FILE *fp) {
  cout << endl
       << vcid << endl
       << "cwd " << cwd << endl
       << "cmdline " << cmdline << endl
       << "sysname  " << uts.sysname << endl
       << "hostname " << uts.nodename << endl
       << "machine  " << uts.machine << endl
       << "user     " << VERuser() << endl;

  cout << "Base name of output files: " << outBase << endl;
  cout << "Text file with list of individual inputs: " << inListFile << endl;

  return;
}

