/**
 * @brief Combine path measures from multiple subjects
 *
 * Combine path measures from multiple subjects
 */
/*
 * Original Author: Anastasia Yendiki
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

const char *Progname = "dmri_group";

int nSection = 0;

string inListFile, outRefFile, outBase;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  unsigned int nmeas, npt, narc, nsubjmin, ntot;
  float distmin, darc, arcmin, arcmax, arc1m, arc2m, arc1s, arc2s,
        lthresh1, lthresh2, uthresh1, uthresh2;
  string listline, filename;
  vector<bool>::const_iterator iout;
  vector<unsigned int>::const_iterator insubj;
  vector<float>::const_iterator icenter, iarc1, iarc2;
  vector< vector<unsigned int> >::const_iterator iallk;
  vector< vector<float> >::const_iterator itemplate, iallm, iallp;
  vector< vector<float> >::iterator ialla;
  vector<bool> isout;
  vector<unsigned int> lengths, nsubj;
  vector<float> meanpath, arcend1, arcend2;
  vector<string> subjlist, measlist;
  vector< vector<unsigned int> > allknots;
  vector< vector<float> > allarc, allpaths, allmeas, allmeasint, allmeassec;
  ifstream listfile;
  ofstream pathfile, pathrasfile;
  MATRIX *outv2r;
  MRI *outref = 0, *outvol = 0;

  nargs = handleVersionOption(argc, argv, "dmri_group");
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

  cputimer.reset();

  // Read output reference volume
  if (!outRefFile.empty()) {
    cout << "Loading output reference volume from " << outRefFile << endl;
    outref = MRIread(outRefFile.c_str());
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
      affinereg.ReadXfm(inputs[2], inref, outref);
      nonlinreg.ReadXfm(inputs[3], outref);
    }
    else
#endif
    if (inputs.size() > 2)
      affinereg.ReadXfm(inputs[2], inref, outref);

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

  // Choose sampling interval for measures along the path
  darc = *(allarc[itemplate - allpaths.begin()].end() - 1)
       / lengths[itemplate - allpaths.begin()];

  cout << "INFO: Sampling interval along path is ";
  if (outref)
    cout << darc * outref->xsize << " mm" << endl;
  else
    cout << darc << " voxels" << endl;

if (0) {
  // Write points of most representative path to file as RAS coords
  filename = outBase + ".median.txt";

  cout << "Writing median path to " << filename << endl;
  pathrasfile.open(filename.c_str(), ios::out);

  pathrasfile << "#!ascii label, vox2ras=scanner" << endl
              << itemplate->size() / 3 << endl;

  npt = 1;

  for (vector<float>::const_iterator ipt = itemplate->begin();
                                     ipt < itemplate->end(); ipt += 3) {
    pathrasfile << npt;

    for (int k = 1; k < 4; k++)		// Transform from voxel to RAS coords
      pathrasfile << " " << ipt[0] * outv2r->rptr[k][1] +
                            ipt[1] * outv2r->rptr[k][2] +
                            ipt[2] * outv2r->rptr[k][3] +
                                     outv2r->rptr[k][4];

    pathrasfile << " 0" << endl;

    npt++;
  }

  pathrasfile.close();
}

  // Reparameterize the arc length on each path
  ialla = allarc.begin();

  for (vector< vector<float> >::const_iterator iallp = allpaths.begin();
                                               iallp < allpaths.end();
                                               iallp++) {
    float dmin = numeric_limits<float>::infinity(),
          arc0 = 0;
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

    arcend1.push_back(*min_element(ialla->begin(), ialla->end()));
    arcend2.push_back(*max_element(ialla->begin(), ialla->end()));

    ialla++;
  }

  // Find outlier paths based on arc length parameterization
  arc1m = accumulate(arcend1.begin(), arcend1.end(), 0.0) / arcend1.size();
  arc2m = accumulate(arcend2.begin(), arcend2.end(), 0.0) / arcend2.size();

  arc1s = sqrt(inner_product(arcend1.begin(), arcend1.end(),
                             arcend1.begin(), 0.0) / arcend1.size()
               - arc1m * arc1m);
  arc2s = sqrt(inner_product(arcend2.begin(), arcend2.end(),
                             arcend2.begin(), 0.0) / arcend2.size()
               - arc2m * arc2m);

  lthresh1 = arc1m - 3*arc1s;
  lthresh2 = arc2m - 3*arc2s;

  uthresh1 = arc1m + 3*arc1s;
  uthresh2 = arc2m + 3*arc2s;

  isout.resize(allarc.size());
  fill(isout.begin(), isout.end(), false);

  iarc1 = arcend1.begin();
  iarc2 = arcend2.begin();

  for (vector<bool>::iterator iout = isout.begin(); iout < isout.end();
                                                    iout++) {
    if (*iarc1 < lthresh1 || *iarc1 > uthresh1 ||
        *iarc2 < lthresh2 || *iarc2 > uthresh2) {
     *iout = true;
     cout << "Found outlier path: " << subjlist[iout - isout.begin()] << endl;
    }

    iarc1++;
    iarc2++;
  }

  // Interpolate measures at the same arc lengths on every path
  arcmin = *min_element(arcend1.begin(), arcend1.end());
  arcmax = *max_element(arcend2.begin(), arcend2.end());

  narc = (unsigned int) floor((arcmax - arcmin) / darc);

  nsubj.resize(narc);
  fill(nsubj.begin(), nsubj.end(), 0);

  meanpath.resize(narc * 3);
  fill(meanpath.begin(), meanpath.end(), 0.0);

  iallm = allmeas.begin();
  iallp = allpaths.begin();
  iout = isout.begin();

  for (vector< vector<float> >::const_iterator ialla = allarc.begin();
                                               ialla < allarc.end(); ialla++) {
    float larc = arcmin + darc;
    vector<unsigned int>::iterator insubj = nsubj.begin();
    vector<float>::iterator imean = meanpath.begin();
    vector<float> meas;

    for (unsigned int ilen = 0; ilen < narc; ilen++) {
      float slope;
      vector<float>::const_iterator iarc = ialla->begin(),
                                    imeas1, imeas0;

      if (*iarc > larc)			// No points in this segment, skip ahead
        for (int k = (int) nmeas; k > 0; k--)
          meas.push_back(numeric_limits<float>::infinity());
      else {
        while (*iarc < larc && iarc < ialla->end())
          iarc++;

        if (iarc == ialla->end()) 	// No points in this segment, skip ahead
          for (int k = (int) nmeas; k > 0; k--)
            meas.push_back(numeric_limits<float>::infinity());
        else {
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

          if (! *iout) {
            vector<float>::const_iterator ipt = iallp->begin()
                                              + 3 * (iarc - ialla->begin());

            // Increment number of samples in this segment
            (*insubj)++;

            // Add point towards mean path
            for (int k = 0; k < 3; k++)
              imean[k] += ipt[k];
          }
        }
      }

      larc += darc;
      insubj++;
      imean += 3;
    }

    allmeasint.push_back(meas);

    iallm++;
    iallp++;
    iout++;
  }

  // Minimum number of subjects that must contribute to a position on the path
  nsubjmin = (unsigned int) ceil(.2 * (float) subjlist.size());

  // Remove positions from start of path that don't have enough samples
  while (*nsubj.begin() < nsubjmin) {
    nsubj.erase(nsubj.begin());

    meanpath.erase(meanpath.begin(), meanpath.begin()+3);

    for (vector< vector<float> >::iterator iallm = allmeasint.begin();
                                           iallm < allmeasint.end(); iallm++)
      iallm->erase(iallm->begin(), iallm->begin() + nmeas);
  }

  // Remove positions from end of path that don't have enough samples
  while (*(nsubj.end()-1) < nsubjmin) {
    nsubj.erase(nsubj.end()-1);

    meanpath.erase(meanpath.end()-3, meanpath.end());

    for (vector< vector<float> >::iterator iallm = allmeasint.begin();
                                           iallm < allmeasint.end(); iallm++)
      iallm->erase(iallm->end() - nmeas, iallm->end());
  }

  // Divide sums of points by number of samples to get mean path points
  insubj = nsubj.begin();

  for (vector<float>::iterator ipt = meanpath.begin(); ipt < meanpath.end();
                                                       ipt += 3) {
    for (int k = 0; k < 3; k++)
      ipt[k] /= *insubj;

    insubj++;
  }

  // Write points of mean path to file as voxel and RAS coords
  filename = outBase + ".coords.mean.txt";

  cout << "Writing mean path voxel coords to " << filename << endl;
  pathfile.open(filename.c_str(), ios::out);

  filename = outBase + ".path.mean.txt";

  cout << "Writing mean path RAS coords to " << filename << endl;
  pathrasfile.open(filename.c_str(), ios::out);

  pathrasfile << "#!ascii label, vox2ras=scanner" << endl
              << meanpath.size() / 3 << endl;

  npt = 1;

  for (vector<float>::iterator ipt = meanpath.begin(); ipt < meanpath.end();
                                                       ipt += 3) {
    // Write voxel coordinates
    pathfile << ipt[0] << " " << ipt[1] << " " << ipt[2] << endl;

    // Write RAS coordinates (in freeview waypoint file format)
    pathrasfile << npt;

    for (int k = 1; k < 4; k++)	// Transform from voxel to RAS coords
      pathrasfile << " " << ipt[0] * outv2r->rptr[k][1] +
                            ipt[1] * outv2r->rptr[k][2] +
                            ipt[2] * outv2r->rptr[k][3] +
                                     outv2r->rptr[k][4];

    pathrasfile << " 0" << endl;

    npt++;
  }

  npt--;

  pathfile.close();
  pathrasfile.close();

  // Write output files
  outvol = MRIallocSequence(npt, 1, 1, MRI_FLOAT, subjlist.size());
  if (outref)
    MRIcopyHeader(outref, outvol);

  ntot = allmeasint[0].size();

  for (vector<string>::const_iterator imeas = measlist.begin();
                                      imeas < measlist.end(); imeas++) {
    int jpt = 0;
    string outname = outBase + "." + *imeas + ".txt";
    ofstream outfile;

    cout << "Writing group table to " << outname << endl;
    outfile.open(outname.c_str(), ios::out);

    // Write subject names
    for (vector<string>::const_iterator isubj = subjlist.begin();
                                        isubj < subjlist.end(); isubj++)
      outfile << *isubj << " ";

    outfile << endl;

    MRIclear(outvol);

    // Write interpolated values of this measure
    for (unsigned ipt = imeas - measlist.begin(); ipt < ntot; ipt += nmeas) {
      int jsubj = 0;

      for (iallm = allmeasint.begin(); iallm < allmeasint.end(); iallm++) {
        vector<float>::const_iterator ival = iallm->begin() + ipt;

        if (*ival == numeric_limits<float>::infinity())
          outfile << "NaN ";
        else
          outfile << *ival << " ";

        MRIsetVoxVal(outvol, jpt, 0, 0, jsubj, *ival);
        jsubj++;
      }

      outfile << endl;
      jpt++;
    }

    outfile.close();

    outname = outBase + "." + *imeas + ".nii.gz";
    cout << "Writing group table to " << outname << endl;
    MRIwrite(outvol, outname.c_str());
  }

  // Average measures over sections along the path
  if (nSection > 0) {
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
      string outname = outBase + "." + *imeas + "."
                                       + to_string(nSection) + "sec.txt";

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

  MRIfree(&outvol);
  if (outref) {
    MRIfree(&outref);
    MatrixFree(&outv2r);
  }

  cputime = cputimer.milliseconds();
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
  << "   --list <file>:" << endl
  << "     Text file with list of individual inputs" << endl
  << "   --ref <file>:" << endl
  << "     Reference volume for output path" << endl
  << "   --out <base>:" << endl
  << "     Base name of output stats files" << endl
  << endl
  << "Optional arguments" << endl
  << "   --sec <num>:" << endl
  << "     Divide the pathway into a number of sections and output " << endl
  << "     average measures for each section" << endl
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
    cout << "ERROR: must specify base name for output files" << endl;
    exit(1);
  }
  if (inListFile.empty()) {
    cout << "ERROR: must specify input list file" << endl;
    exit(1);
  }
  return;
}

static void dump_options(FILE *fp) {
  cout << endl
       << getVersion() << endl
       << "cwd " << cwd << endl
       << "cmdline " << cmdline << endl
       << "sysname  " << uts.sysname << endl
       << "hostname " << uts.nodename << endl
       << "machine  " << uts.machine << endl
       << "user     " << VERuser() << endl;

  cout << "Base name of output files: " << outBase << endl;
  cout << "Text file with list of individual inputs: " << inListFile << endl;
  if (!outRefFile.empty())
    cout << "Reference volume for output path: " << outRefFile << endl;
  if (nSection > 0)
    cout << "Number of path sections: " << nSection << endl;

  return;
}

