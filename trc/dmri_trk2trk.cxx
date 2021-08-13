/**
 * @brief Transform streamlines in .trk file
 *
 * Transform streamlines in .trk file:
 * Apply affine or non-linear warp, inclusion or exclusion masks,
 * convert to text file or volume file, etc.
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

#include "vial.h"	// Needs to be included first because of CVS libs
#include "spline.h"
#include "TrackIO.h"

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
#include <sstream>
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

const char *Progname = "dmri_trk2trk";

bool doMerge = false;
int doInvXfm = 0, doFill = 0, doMean = 0, doNearMean = 0,
    doNth = 0, strNum = -1, doEvery = 0, everyNum = -1, doSmooth = 0,
    lengthMin = -1, lengthMax = -1;
unsigned int nTract = 0;
string inDir, outDir, inRefFile, outRefFile, affineXfmFile, nonlinXfmFile;
vector<string> inTrkList, inAscList, outTrkList, outAscList, outVolList,
               overList, overnames,
               incMaskList, excMaskList, incTermMaskList, excTermMaskList;
vector< vector<float> > overvals;
vector<MRI *>  incMask, excMask, incTermMask, excTermMask;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  char outorient[4];
  string fname;
  vector<float> point(3), fillstep(3, 0);
  vector< vector<float> > streamlines, overlays, properties;
  MATRIX *outv2r;
  MRI *inref = 0, *outref = 0, *outvol = 0;
  AffineReg affinereg;
#ifndef NO_CVS_UP_IN_HERE
  NonlinReg nonlinreg;
#endif

  nargs = handleVersionOption(argc, argv, "dmri_trk2trk");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd, 2000);

  Progname = argv[0];
  argc--;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);

  dump_options(stdout);

  // Read reference volumes
  if (!inRefFile.empty())
    inref = MRIread(inRefFile.c_str());

  if (!outRefFile.empty())
    outref = MRIread(outRefFile.c_str());

  if (!outVolList.empty())
    outvol = MRIclone(outref, NULL);

  // Output space orientation information
  if (outref) {
    outv2r = MRIgetVoxelToRasXform(outref);
    MRIdircosToOrientationString(outref, outorient);
  }

  // Read transform files
  if (!affineXfmFile.empty()) {
    if (doInvXfm) {
      affinereg.ReadXfm(affineXfmFile, outref, inref);

      if (affinereg.IsInvEmpty()) {
        cout << "ERROR: For --inv, use LTA format" << endl;
        exit(1);
      }
    }
    else
      affinereg.ReadXfm(affineXfmFile, inref, outref);
  }

#ifndef NO_CVS_UP_IN_HERE
  if (!nonlinXfmFile.empty()) {
    if (doInvXfm)
      nonlinreg.ReadXfm(nonlinXfmFile, inref);
    else
      nonlinreg.ReadXfm(nonlinXfmFile, outref);
  }
#endif

  // Read scalar overlay vectors
  for (vector<string>::const_iterator iover = overList.begin();
                                      iover < overList.end(); iover++) {
    int ndim = 0;
    vector<float> overlay;
    MRI *invol = MRIread((*iover).c_str());
    string oname(*iover);

    // Check dimensions of input overlay volume
    if (invol->width   > 1)	ndim++;
    if (invol->height  > 1)	ndim++;
    if (invol->depth   > 1)	ndim++;
    if (invol->nframes > 1)	ndim++;

    if (ndim > 1)
      cout << "WARN: Overlay " << *iover << " has more than one "
           << "non-singular dimension - treating as vector" << endl;

    // Read overlay values
    for (int it = 0; it < invol->nframes; it++)
      for (int iz = 0; iz < invol->depth; iz++)
        for (int iy = 0; iy < invol->height; iy++)
          for (int ix = 0; ix < invol->width; ix++)
            overlay.push_back(MRIgetVoxVal(invol, ix, iy, iz, it));

    overvals.push_back(overlay);

    // Use file name as overlay name
    if (oname.rfind("/") != string::npos)
      oname = oname.substr(oname.rfind("/")+1, string::npos);

    overnames.push_back(oname);

    MRIfree(&invol);
  }

  // Read inclusion masks
  for (vector<string>::const_iterator imask = incMaskList.begin();
                                      imask < incMaskList.end(); imask++)
    incMask.push_back(MRIread((*imask).c_str()));

  // Read exclusion masks
  for (vector<string>::const_iterator imask = excMaskList.begin();
                                      imask < excMaskList.end(); imask++)
    excMask.push_back(MRIread((*imask).c_str()));

  // Read terminal inclusion masks
  for (vector<string>::const_iterator imask = incTermMaskList.begin();
                                      imask < incTermMaskList.end(); imask++)
    incTermMask.push_back(MRIread((*imask).c_str()));

  // Read terminal exclusion masks
  for (vector<string>::const_iterator imask = excTermMaskList.begin();
                                      imask < excTermMaskList.end(); imask++)
    excTermMask.push_back(MRIread((*imask).c_str()));

  for (unsigned int itract = 0; itract < nTract; itract++) {
    int npts, nstr = 0;
    CTrackReader trkreader;
    TRACK_HEADER trkheadin;

    cout << "Processing input file " << itract+1 << " of " << nTract
         << "..." << endl;
    cputimer.reset();

    if (!doMerge) {
      streamlines.clear();
      overlays.clear();
      properties.clear();
    }

    if (doEvery)
      strNum = 0;

    if (!inTrkList.empty()) {		// Read streamlines from .trk file
      fname = inTrkList[itract];

      if (!inDir.empty())
        fname = inDir + "/" + fname;

      if (!trkreader.Open(fname.c_str(), &trkheadin)) {
        cout << "ERROR: Cannot open input file " << fname << endl;
        cout << "ERROR: " << trkreader.GetLastErrorMessage() << endl;
        exit(1);
      }

      while (trkreader.GetNextPointCount(&npts)) {
        const int veclen  = npts * 3,
                  overlen = npts * trkheadin.n_scalars;
        float *iraw, *rawpts  = new float[veclen],
                     *scalars = new float[overlen],
                     *props   = new float[trkheadin.n_properties];
        vector<float> newpts(veclen);

        // Read a streamline from input file
        trkreader.GetNextTrackData(npts, rawpts, scalars, props);

        if ( ((doNth || doEvery) && nstr != strNum) ||
             (lengthMin > -1 && npts <= lengthMin) ||
             (lengthMax > -1 && npts >= lengthMax) ) {
          delete[] rawpts;
          delete[] scalars;
          delete[] props;
          nstr++;
          continue;
        }

        if (doEvery && nstr == strNum)
          strNum += everyNum;

        iraw = rawpts;

        // Divide by input voxel size and make 0-based to get voxel coords
        for (vector<float>::iterator ipt = newpts.begin(); ipt < newpts.end();
                                                           ipt += 3)
          for (int k = 0; k < 3; k++) {
            ipt[k] = *iraw / trkheadin.voxel_size[k] - .5;
            iraw++;
          }

        delete[] rawpts;
        streamlines.push_back(newpts);

        // Store scalar overlays of input streamlines
        if (overlen > 0) {
          vector<float> newscalars(overlen);

          copy(scalars, scalars+overlen, newscalars.begin());

          delete[] scalars;
          overlays.push_back(newscalars);
        }

        // Store properties of input streamlines
        if (trkheadin.n_properties > 0) {
          vector<float> newprops(trkheadin.n_properties);

          copy(props, props+trkheadin.n_properties, newprops.begin());

          delete[] props;
          properties.push_back(newprops);
        }

        nstr++;
      }
    }
    else if (!inAscList.empty()) {	// Read streamlines from text file
      string ptline;
      ifstream infile;
      vector<float> newpts;

      fname = inAscList[itract];

      if (!inDir.empty())
        fname = inDir + "/" + fname;

      infile.open(fname, ios::in);
      if (!infile) {
        cout << "Error: Cannot open input file " << fname << endl;
        exit(1);
      }

      while (getline(infile, ptline)) {
        float val;
        istringstream linestr(ptline);
        vector<float> point;

        while (linestr >> val)
          point.push_back(val);

        if (point.empty()) {		// Empty line marks end of streamline
          if ( (!(doNth || doEvery) || nstr == strNum) &&
               (lengthMin == -1 || (int) newpts.size()/3 > lengthMin) &&
               (lengthMax == -1 || (int) newpts.size()/3 < lengthMax) )
            streamlines.push_back(newpts);

          if (doEvery && nstr == strNum)
            strNum += everyNum;

          newpts.clear();
          nstr++;
        }
        else if (point.size() != 3) {
          cout << "ERROR: Unexpected number of entries in a line ("
               << point.size() << ") in file " << fname << endl;
          exit(1);
        }
        else
          newpts.insert(newpts.end(), point.begin(), point.end());
      }

      infile.close();
    }

    if (doMerge && itract < nTract-1)
      continue;

    nstr = streamlines.size();

    // Apply transformations
    for (int kstr = nstr-1; kstr >= 0; kstr--) {
      vector<float> newpts;

      for (vector<float>::iterator ipt = streamlines[kstr].begin();
                                   ipt < streamlines[kstr].end(); ipt += 3) {
        copy(ipt, ipt+3, point.begin());

        if (doInvXfm) {
#ifndef NO_CVS_UP_IN_HERE
          // Apply inverse of nonlinear transform
          if (!nonlinreg.IsEmpty())
            nonlinreg.ApplyXfmInv(point, point.begin());
#endif

          // Apply inverse of affine transform
          if (!affinereg.IsEmpty())
            affinereg.ApplyXfmInv(point, point.begin());
        }
        else {
          // Apply affine transform
          if (!affinereg.IsEmpty())
            affinereg.ApplyXfm(point, point.begin());

#ifndef NO_CVS_UP_IN_HERE
          // Apply nonlinear transform
          if (!nonlinreg.IsEmpty())
            nonlinreg.ApplyXfm(point, point.begin());
#endif
        }

        copy(point.begin(), point.end(), ipt);
      }

      for (vector<float>::const_iterator ipt = streamlines[kstr].begin();
                                         ipt < streamlines[kstr].end();
                                         ipt += 3) {
        float dmax = 1;		// This will not remove duplicate points

        if (doFill && ipt < streamlines[kstr].end()-3) {
          // Calculate step for filling in gap between points
          // Gaps could result when mapping to a higher-resolution space
          for (int k = 0; k < 3; k++) {
            float dist = ipt[k+3] - ipt[k];

            fillstep[k] = dist;
            dist = fabs(dist);

            if (dist > dmax)
              dmax = dist;
          }

          if (dmax > 0)
            for (int k = 0; k < 3; k++)
              fillstep[k] /= dmax;
        }

        copy(ipt, ipt+3, point.begin());

        for (int istep = (int) round(dmax); istep > 0; istep--) {
          newpts.insert(newpts.end(), point.begin(), point.end());

          for (int k = 0; k < 3; k++)
            point[k] += fillstep[k];
        }
      }

      streamlines[kstr].resize(newpts.size());
      copy(newpts.begin(), newpts.end(), streamlines[kstr].begin());
    }

    // Apply inclusion/exclusion masks
    nstr = streamlines.size();

    for (int kstr = nstr-1; kstr >= 0; kstr--) {
      bool dokeep = true;

      // There must be an endpoint that intersects each terminal inclusion mask
      for (vector<MRI *>::const_iterator imask = incTermMask.begin();
                                         imask < incTermMask.end(); imask++) {
        vector<float>::const_iterator ipt = streamlines[kstr].begin();
        int ix = (int) round(ipt[0]),
            iy = (int) round(ipt[1]),
            iz = (int) round(ipt[2]);

        if (ix < 0)                   ix = 0;
        if (ix >= (*imask)->width)    ix = (*imask)->width-1;
        if (iy < 0)                   iy = 0;
        if (iy >= (*imask)->height)   iy = (*imask)->height-1;
        if (iz < 0)                   iz = 0;
        if (iz >= (*imask)->depth)    iz = (*imask)->depth-1;

        dokeep = false;

        if (MRIgetVoxVal(*imask, ix, iy, iz, 0) > 0)
           dokeep = true;
        else {
          ipt = streamlines[kstr].end() - 3;
          ix = (int) round(ipt[0]);
          iy = (int) round(ipt[1]);
          iz = (int) round(ipt[2]);

          if (ix < 0)                   ix = 0;
          if (ix >= (*imask)->width)    ix = (*imask)->width-1;
          if (iy < 0)                   iy = 0;
          if (iy >= (*imask)->height)   iy = (*imask)->height-1;
          if (iz < 0)                   iz = 0;
          if (iz >= (*imask)->depth)    iz = (*imask)->depth-1;

          if (MRIgetVoxVal(*imask, ix, iy, iz, 0) > 0)
           dokeep = true;
        }

        if (!dokeep)
          break;
      }

      if (!dokeep) {
        streamlines.erase(streamlines.begin() + kstr);
        if (!overlays.empty())
          overlays.erase(overlays.begin() + kstr);
        if (!properties.empty())
          properties.erase(properties.begin() + kstr);
        continue;
      }

      // There must be at least one point that intersects each inclusion mask
      for (vector<MRI *>::const_iterator imask = incMask.begin();
                                         imask < incMask.end(); imask++) {
        dokeep = false;

        for (vector<float>::const_iterator ipt = streamlines[kstr].begin();
                                           ipt < streamlines[kstr].end();
                                           ipt += 3) {
          int ix = (int) round(ipt[0]),
              iy = (int) round(ipt[1]),
              iz = (int) round(ipt[2]);

          if (ix < 0)                   ix = 0;
          if (ix >= (*imask)->width)    ix = (*imask)->width-1;
          if (iy < 0)                   iy = 0;
          if (iy >= (*imask)->height)   iy = (*imask)->height-1;
          if (iz < 0)                   iz = 0;
          if (iz >= (*imask)->depth)    iz = (*imask)->depth-1;

          if (MRIgetVoxVal(*imask, ix, iy, iz, 0) > 0) {
             dokeep = true;
             break;
          }
        }

        if (!dokeep)
          break;
      }

      if (!dokeep) {
        streamlines.erase(streamlines.begin() + kstr);
        if (!overlays.empty())
          overlays.erase(overlays.begin() + kstr);
        if (!properties.empty())
          properties.erase(properties.begin() + kstr);
        continue;
      }

      // There must be no endpoint that intersects any terminal exclusion mask
      for (vector<MRI *>::const_iterator imask = excTermMask.begin();
                                         imask < excTermMask.end(); imask++) {
        vector<float>::const_iterator ipt = streamlines[kstr].begin();
        int ix = (int) round(ipt[0]),
            iy = (int) round(ipt[1]),
            iz = (int) round(ipt[2]);

        if (ix < 0)                   ix = 0;
        if (ix >= (*imask)->width)    ix = (*imask)->width-1;
        if (iy < 0)                   iy = 0;
        if (iy >= (*imask)->height)   iy = (*imask)->height-1;
        if (iz < 0)                   iz = 0;
        if (iz >= (*imask)->depth)    iz = (*imask)->depth-1;

        if (MRIgetVoxVal(*imask, ix, iy, iz, 0) > 0) 
          dokeep = false;
        else {
          ipt = streamlines[kstr].end() - 3;
          ix = (int) round(ipt[0]),
          iy = (int) round(ipt[1]),
          iz = (int) round(ipt[2]);

          if (ix < 0)                   ix = 0;
          if (ix >= (*imask)->width)    ix = (*imask)->width-1;
          if (iy < 0)                   iy = 0;
          if (iy >= (*imask)->height)   iy = (*imask)->height-1;
          if (iz < 0)                   iz = 0;
          if (iz >= (*imask)->depth)    iz = (*imask)->depth-1;

          if (MRIgetVoxVal(*imask, ix, iy, iz, 0) > 0) 
            dokeep = false;
        }

        if (!dokeep)
          break;
      }

      if (!dokeep) {
        streamlines.erase(streamlines.begin() + kstr);
        if (!overlays.empty())
          overlays.erase(overlays.begin() + kstr);
        if (!properties.empty())
          properties.erase(properties.begin() + kstr);
        continue;
      }

      // There must be no point that intersects any exclusion mask
      for (vector<float>::const_iterator ipt = streamlines[kstr].begin();
                                         ipt < streamlines[kstr].end();
                                         ipt += 3) {
        for (vector<MRI *>::const_iterator imask = excMask.begin();
                                           imask < excMask.end(); imask++) {
          int ix = (int) round(ipt[0]),
              iy = (int) round(ipt[1]),
              iz = (int) round(ipt[2]);

          if (ix < 0)                   ix = 0;
          if (ix >= (*imask)->width)    ix = (*imask)->width-1;
          if (iy < 0)                   iy = 0;
          if (iy >= (*imask)->height)   iy = (*imask)->height-1;
          if (iz < 0)                   iz = 0;
          if (iz >= (*imask)->depth)    iz = (*imask)->depth-1;

          if (MRIgetVoxVal(*imask, ix, iy, iz, 0) > 0) {
            dokeep = false;
            break;
          }
        }

        if (!dokeep)
          break;
      }

      if (!dokeep) {
        streamlines.erase(streamlines.begin() + kstr);
        if (!overlays.empty())
          overlays.erase(overlays.begin() + kstr);
        if (!properties.empty())
          properties.erase(properties.begin() + kstr);
      }
    }

    // Find mean of streamlines
    if (doMean && !streamlines.empty()) {
      unsigned int nsteps;
      vector<float>::const_iterator imean0;
      StreamSet bundle;

      bundle.SetLengths(streamlines);
      nsteps = bundle.SetNumStepsAvgLength();
      bundle.ComputeSteps();
      bundle.ComputeMeanStreamline(streamlines);
      imean0 = bundle.GetStreamlineMean();

      // Keep only the mean streamline for writing to disk
      streamlines.clear();
      streamlines.resize(1);
      streamlines[0].insert(streamlines[0].begin(), imean0, imean0 + nsteps*3);

      overlays.clear();
      properties.clear();
    }

    // Find streamline nearest to the mean
    if (doNearMean && !streamlines.empty()) {
      unsigned int kstrmean;
      StreamSet bundle;

      bundle.SetLengths(streamlines);
      bundle.SetNumStepsMinLength();
      bundle.ComputeSteps();
      bundle.ComputeMeanStdStreamline(streamlines);
      kstrmean = bundle.FindMeanNearestStreamline(streamlines);

      // Keep only the chosen streamline for writing to disk
      streamlines.erase(streamlines.begin(), streamlines.begin() + kstrmean);
      streamlines.erase(streamlines.begin() + 1, streamlines.end());

      if (!overlays.empty()) {
        overlays.erase(overlays.begin(), overlays.begin() + kstrmean);
        overlays.erase(overlays.begin() + 1, overlays.end());
      }

      if (!properties.empty()) {
        properties.erase(properties.begin(), properties.begin() + kstrmean);
        properties.erase(properties.begin() + 1, properties.end());
      }
    }

    // Smooth streamlines (assumes integer coordinates for now)
    if (doSmooth) {
      vector<int> strint;
      vector<float> strsmooth;
      vector<int>::iterator iptint;

      for (vector< vector<float> >::iterator istr = streamlines.begin();
                                             istr < streamlines.end(); istr++) {
        strint.resize(istr->size());
        strsmooth.resize(istr->size());

        iptint = strint.begin();
        
        for (vector<float>::const_iterator ipt = istr->begin();
                                           ipt < istr->end(); ipt ++) {
          *iptint = (int) round(*ipt);
          iptint++;
        }

        CurveSmooth(strsmooth, strint);

        copy(strsmooth.begin(), strsmooth.end(), istr->begin());
      }
    }

    // Assign a scalar overlay value to each point on each streamline
    if (!overList.empty()) {
      unsigned int nscalar = overlays.empty() ? 0 :
                             (unsigned int) trkheadin.n_scalars;
      vector<float>::const_iterator imean0;
      StreamSet bundle;

      bundle.SetLengths(streamlines);

      if (overlays.empty())
        overlays.resize(streamlines.size());

      for (vector< vector<float> >::const_iterator ival = overvals.begin();
                                                   ival < overvals.end();
                                                   ival++) {
        unsigned int nsteps = ival->size();

        // Compute mean streamline with as many points as the overlay vector
        if (bundle.GetNumSteps() != nsteps) {
          bundle.SetNumSteps(nsteps);
          bundle.ComputeSteps();
          bundle.ComputeMeanStreamline(streamlines);
          imean0 = bundle.GetStreamlineMean();
        }

        for (vector< vector<float> >::const_iterator istr = streamlines.begin();
                                                     istr < streamlines.end();
                                                     istr++) {
          for (vector<float>::const_iterator ipt = istr->begin();
                                             ipt < istr->end(); ipt += 3) {
            unsigned int knear = 0;
            float dmin = numeric_limits<float>::infinity();
            vector<float>::const_iterator imean = imean0;
            vector<float>::iterator iover_pt;
            vector< vector<float> >::iterator iover_str;

            // Find nearest point on the mean streamline
            for (unsigned int k = 0; k < nsteps; k++) {
              const float dx = ipt[0] - imean[0],
                          dy = ipt[1] - imean[1],
                          dz = ipt[2] - imean[2],
                          dist2 = dx*dx + dy*dy + dz*dz;

              if (dist2 < dmin) {
                dmin = dist2;
                knear = k;
              }

              imean += 3;
            }

            // Insert the corresponding scalar from the overlay vector
            iover_str = overlays.begin() + (istr - streamlines.begin());
            iover_pt = iover_str->begin() + (ipt - istr->begin()) / 3 * nscalar;

            iover_str->insert(iover_pt + nscalar, *(ival->begin() + knear));
          }
        }

        nscalar++;
      }
    }

    // Write transformed streamlines to volume
    if (!outVolList.empty()) {
      if (doMerge)
        fname = outVolList[0];
      else
        fname = outVolList[itract];

      if (!outDir.empty())
        fname = outDir + "/" + fname;

      MRIclear(outvol);

      for (vector< vector<float> >::const_iterator istr = streamlines.begin();
                                                   istr < streamlines.end();
                                                   istr++)
        for (vector<float>::const_iterator ipt = istr->begin();
                                           ipt < istr->end(); ipt += 3) {
          int ix = (int) round(ipt[0]),
              iy = (int) round(ipt[1]),
              iz = (int) round(ipt[2]);

          if (ix < 0)			ix = 0;
          if (ix >= outvol->width)	ix = outvol->width-1;
          if (iy < 0)			iy = 0;
          if (iy >= outvol->height)	iy = outvol->height-1;
          if (iz < 0)			iz = 0;
          if (iz >= outvol->depth)	iz = outvol->depth-1;

          MRIsetVoxVal(outvol, ix, iy, iz, 0,
                       MRIgetVoxVal(outvol, ix, iy, iz, 0) + 1);
        }

      MRIwrite(outvol, fname.c_str());
    }

    // Write transformed streamlines to text file
    if (!outAscList.empty()) {
      ofstream outfile;

      if (doMerge)
        fname = outAscList[0];
      else
        fname = outAscList[itract];

      if (!outDir.empty())
        fname = outDir + "/" + fname;

      outfile.open(fname, ios::out);
      if (!outfile) {
        cout << "ERROR: Could not open " << fname << " for writing" << endl;
        exit(1);
      }

      for (vector< vector<float> >::const_iterator istr = streamlines.begin();
                                                   istr < streamlines.end();
                                                   istr++) {
        for (vector<float>::const_iterator ipt = istr->begin();
                                           ipt < istr->end(); ipt += 3)
          outfile << (int) round(ipt[0]) << " "
                  << (int) round(ipt[1]) << " "
                  << (int) round(ipt[2]) << endl;

        outfile << endl;
      }

      outfile.close();
    }

    // Write transformed streamlines to .trk file
    if (!outTrkList.empty()) {
      CTrackWriter trkwriter;
      TRACK_HEADER trkheadout;
      vector< vector<float> >::iterator iover = overlays.begin(),
                                        iprop = properties.begin();

      // Set output .trk header
      if (inTrkList.empty()) {
        trkheadout.Initialize();

        trkheadout.origin[0] = 0;
        trkheadout.origin[1] = 0;
        trkheadout.origin[2] = 0;
      }
      else
        trkheadout = trkheadin;

      if (outref) {
        trkheadout.voxel_size[0] = outref->xsize;
        trkheadout.voxel_size[1] = outref->ysize;
        trkheadout.voxel_size[2] = outref->zsize;

        trkheadout.dim[0] = outref->width;
        trkheadout.dim[1] = outref->height;
        trkheadout.dim[2] = outref->depth;

        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            trkheadout.vox_to_ras[i][j] = outv2r->rptr[i+1][j+1];

        strcpy(trkheadout.voxel_order, outorient);

        // Find patient-to-scanner coordinate transform:
        // Take x and y vectors from vox2RAS matrix, convert to LPS,
        // divide by voxel size
        trkheadout.image_orientation_patient[0] = 
          - trkheadout.vox_to_ras[0][0] / trkheadout.voxel_size[0];
        trkheadout.image_orientation_patient[1] = 
          - trkheadout.vox_to_ras[1][0] / trkheadout.voxel_size[0];
        trkheadout.image_orientation_patient[2] = 
            trkheadout.vox_to_ras[2][0] / trkheadout.voxel_size[0];
        trkheadout.image_orientation_patient[3] = 
          - trkheadout.vox_to_ras[0][1] / trkheadout.voxel_size[1];
        trkheadout.image_orientation_patient[4] = 
          - trkheadout.vox_to_ras[1][1] / trkheadout.voxel_size[1];
        trkheadout.image_orientation_patient[5] = 
            trkheadout.vox_to_ras[2][1] / trkheadout.voxel_size[1];
      }

      trkheadout.n_count = (int) streamlines.size();

      // In case I have cleared the old overlays/properties (if using --mean)
      if (overlays.empty())
        trkheadout.n_scalars = 0;

      if (properties.empty())
        trkheadout.n_properties = 0;

      // Add names of new scalar overlays, if any
      for (vector<string>::const_iterator iname = overnames.begin();
                                          iname < overnames.end(); iname++) {
        if (trkheadout.n_scalars == 20) {
          cout << "ERROR: Exceeded max number of 20 scalar overlays" << endl;
          exit(1);
        }

        strcpy(trkheadout.scalar_name[trkheadout.n_scalars], (*iname).c_str());

        trkheadout.n_scalars++;
      }

      // Open output .trk file
      if (doMerge)
        fname = outTrkList[0];
      else
        fname = outTrkList[itract];

      if (!outDir.empty())
        fname = outDir + "/" + fname;

      if (!trkwriter.Initialize(fname.c_str(), trkheadout)) {
        cout << "ERROR: Cannot open output file " << fname << endl;
        cout << "ERROR: " << trkwriter.GetLastErrorMessage() << endl;
        exit(1);
      }

      for (vector< vector<float> >::iterator istr = streamlines.begin();
                                             istr < streamlines.end(); istr++) {
        // Make .5-based and multiply back by output voxel size
        for (vector<float>::iterator ipt = istr->begin(); ipt < istr->end();
                                                          ipt += 3)
          for (int k = 0; k < 3; k++)
            ipt[k] = (ipt[k] + .5) * trkheadout.voxel_size[k];

        if (overlays.empty()) {
          if (properties.empty())
            trkwriter.WriteNextTrack(istr->size()/3, &(istr->at(0)));
          else {
            // Transfer properties from input to output streamlines
            trkwriter.WriteNextTrack(istr->size()/3, &(istr->at(0)),
                                     NULL, &(iprop->at(0)));
            iprop++;
          }
        }
        else {
          if (properties.empty()) {
            // Transfer scalar overlays from input to output streamlines
            trkwriter.WriteNextTrack(istr->size()/3, &(istr->at(0)),
                                     &(iover->at(0)), NULL);
            iover++;
          }
          else {
            // Transfer overlays & properties from input to output streamlines
            trkwriter.WriteNextTrack(istr->size()/3, &(istr->at(0)),
                                     &(iover->at(0)), &(iprop->at(0)));
            iover++;
            iprop++;
          }
        }
      }

      trkwriter.Close();
    }

    cputime = cputimer.milliseconds();
    cout << "Done in " << cputime/1000.0 << " sec." << endl;
  }

  if (inref)
    MRIfree(&inref);
  if (outref) {
    MatrixFree(&outv2r);
    MRIfree(&outref);
  }
  if (!outVolList.empty())
    MRIfree(&outvol);
  for (vector<MRI *>::iterator imask = incMask.begin();
                               imask < incMask.end(); imask++)
    MRIfree(&(*imask));
  for (vector<MRI *>::iterator imask = excMask.begin();
                               imask < excMask.end(); imask++)
    MRIfree(&(*imask));
  for (vector<MRI *>::iterator imask = incTermMask.begin();
                               imask < incTermMask.end(); imask++)
    MRIfree(&(*imask));
  for (vector<MRI *>::iterator imask = excTermMask.begin();
                               imask < excTermMask.end(); imask++)
    MRIfree(&(*imask));

  cout << "dmri_trk2trk done" << endl;
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
    else if (!strcmp(option, "--indir")) {
      if (nargc < 1) CMDargNErr(option,1);
      inDir = fio_fullpath(pargv[0]);
      nargsused = 1;
    }
    else if (!strcmp(option, "--in")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        inTrkList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--inasc")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        inAscList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--outdir")) {
      if (nargc < 1) CMDargNErr(option,1);
      outDir = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--out")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        outTrkList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--outasc")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        outAscList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--outvol")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        outVolList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--inref")) {
      if (nargc < 1) CMDargNErr(option,1);
      inRefFile = fio_fullpath(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--outref")) {
      if (nargc < 1) CMDargNErr(option,1);
      outRefFile = fio_fullpath(pargv[0]);
      nargsused = 1;
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
    else if (!strcasecmp(option, "--inv"))
      doInvXfm = 1;
    else if (!strcasecmp(option, "--fill"))
      doFill = 1;
    else if (!strcmp(option, "--over")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        overList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--imask")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        incMaskList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--emask")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        excMaskList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--itmask")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        incTermMaskList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcmp(option, "--etmask")) {
      if (nargc < 1) CMDargNErr(option,1);
      nargsused = 0;
      while (nargsused < nargc && strncmp(pargv[nargsused], "--", 2)) {
        excTermMaskList.push_back(pargv[nargsused]);
        nargsused++;
      }
    }
    else if (!strcasecmp(option, "--lmin")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &lengthMin);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--lmax")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &lengthMax);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--mean"))
      doMean = 1;
    else if (!strcasecmp(option, "--nearmean"))
      doNearMean = 1;
    else if (!strcasecmp(option, "--nth")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &strNum);
      nargsused = 1;
      doNth = 1;
    }
    else if (!strcasecmp(option, "--every")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &everyNum);
      nargsused = 1;
      doEvery = 1;
    }
    else if (!strcasecmp(option, "--smooth"))
      doSmooth = 1;
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
  << "Basic inputs" << endl
  << "   --in <file> [...]:" << endl
  << "     Input .trk file(s)" << endl
  << "   --inasc <file> [...]:" << endl
  << "     Input ASCII plain text file(s), as an alternative to .trk" << endl
  << "     Each line in the text file must have the voxel coordinates of a single point of a streamline, with empty lines between streamlines)" << endl
  << "   --indir <dir>:" << endl
  << "     Input directory (optional)" << endl
  << "     If specified, names of input .trk files are relative to this" << endl
  << "   --out <file> [...]:" << endl
  << "     Output .trk file(s), as many as inputs (or 1 to merge inputs)" << endl
  << "   --outasc <file> [...]:" << endl
  << "     Output ASCII plain text file(s), as many as inputs (or 1 to merge inputs)" << endl
  << "   --outvol <file> [...]:" << endl
  << "     Output volume(s), as many as inputs (or 1 to merge inputs)" << endl
  << "   --outdir <dir>:" << endl
  << "     Output directory (optional)" << endl
  << "     If specified, names of output .trk files and volumes are relative to this)" << endl
  << "   --inref <file>:" << endl
  << "     Input reference volume (needed for --reg/--regnl)" << endl
  << "   --outref <file>:" << endl
  << "     Output reference volume (needed for --reg/--regnl/--outvol)" << endl
  << "   --reg <file>:" << endl
  << "     Affine registration file (.lta or .mat), applied first" << endl
  << "   --regnl <file>:" << endl
  << "     Nonlinear registration file (.m3z), applied second" << endl
  << "   --inv:" << endl
  << "     Apply inverse of registration (default: no)" << endl
  << "   --fill:" << endl
  << "     Fill gaps b/w mapped points by linear interpolation" << endl
  << "     (Default: don't fill)" << endl
  << "   --over <file> [...]:" << endl
  << "     Scalar overlay 1D volume(s), applied to all input .trk files" << endl
  << "   --imask <file> [...]:" << endl
  << "     Inclusion mask(s), applied to all input .trk files" << endl
  << "   --emask <file> [...]:" << endl
  << "     Exclusion mask(s), applied to all input .trk files" << endl
  << "   --itmask <file> [...]:" << endl
  << "     Terminal inclusion mask(s), applied to all input .trk files" << endl
  << "   --etmask <file> [...]:" << endl
  << "     Terminal exclusion mask(s), applied to all input .trk files" << endl
  << "   --lmin <num>:" << endl
  << "     Only save streamlines with length greater than this number" << endl
  << "   --lmax <num>:" << endl
  << "     Only save streamlines with length smaller than this number" << endl
  << "   --mean:" << endl
  << "     Only save the mean of the streamlines (Default: save all)" << endl
  << "   --nearmean:" << endl
  << "     Only save the streamline nearest to the mean (Default: save all)" << endl
  << "   --nth <num>:" << endl
  << "     Only save the n-th (0-based) streamline (Default: save all)" << endl
  << "   --every <num>:" << endl
  << "     Only save every n-th streamline (Default: save all)" << endl
  << "   --smooth:" << endl
  << "     Smooth streamlines (default: no)" << endl
  << endl
  << "Other options" << endl
  << "   --debug:     turn on debugging" << endl
  << "   --checkopts: don't run anything, just check options and exit" << endl
  << "   --help:      print out information on how to use this program" << endl
  << "   --version:   print out version and exit" << endl
  << endl
  << "Order of operations (all optional):" << endl
  << "   1. Keep the n-th streamline or every n-th streamline" << endl
  << "   2. Apply streamline length threshold(s)" << endl
  << "   3. Apply affine transform" << endl
  << "   4. Apply non-linear transform" << endl
  << "   5. Apply inclusion mask(s)" << endl
  << "   6. Apply exclusion mask(s)" << endl
  << "   7. Find mean streamline" << endl
  << "   8. Smooth streamline(s)" << endl
  << endl;
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage();

  cout << endl
       << "..." << endl
       << endl;

  exit(1) ;
}

/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage();
  exit(1) ;
}

/* --------------------------------------------- */
static void print_version(void) {
  cout << getVersion() << endl;
  exit(1) ;
}

/* --------------------------------------------- */
static void check_options(void) {
  if (inTrkList.empty() && inAscList.empty()) {
    cout << "ERROR: must specify input .trk or text file(s)" << endl;
    exit(1);
  }
  if (!inTrkList.empty() && !inAscList.empty()) {
    cout << "ERROR: cannot use both --in and --inasc at the same time" << endl;
    exit(1);
  }

  nTract = inTrkList.size() + inAscList.size();

  if (nTract > 1 && (outTrkList.size() == 1 || outAscList.size() == 1 
                                            || outVolList.size() == 1))
    doMerge = true;

  if (outTrkList.empty() && outAscList.empty() && outVolList.empty()) {
    cout << "ERROR: must specify output .trk or text or volume file(s)" << endl;
    exit(1);
  }
  if (outTrkList.size() > 1 && outTrkList.size() != nTract) {
    cout << "ERROR: number of output .trk files (" << outTrkList.size()
         << ") does not match number of input files (" << nTract << ")" << endl;
    exit(1);
  }
  if (outAscList.size() > 1 && outAscList.size() != nTract) {
    cout << "ERROR: number of output text files (" << outAscList.size()
         << ") does not match number of input files (" << nTract << ")" << endl;
    exit(1);
  }
  if (outVolList.size() > 1 && outVolList.size() != nTract) {
    cout << "ERROR: number of output volumes (" << outVolList.size()
         << ") does not match number of input files (" << nTract << ")" << endl;
    exit(1);
  }
  if (inRefFile.empty() && (!affineXfmFile.empty() || !nonlinXfmFile.empty())) {
    cout << "ERROR: must specify input reference volume when using "
         << "--reg or --regnl" << endl;
    exit(1);
  }
  if (outRefFile.empty() && (!affineXfmFile.empty() || !nonlinXfmFile.empty()
                                                    || !outVolList.empty())) {
    cout << "ERROR: must specify output reference volume when using "
         << "--reg, --regnl, or --outvol" << endl;
    exit(1);
  }
  if (doMean && doNearMean) {
    cout << "ERROR: cannot use both --mean and --nearmean" << endl;
    exit(1);
  }
  if (doMean && doNth) {
    cout << "ERROR: cannot use both --mean and --nth" << endl;
    exit(1);
  }
  if (doEvery && doNth) {
    cout << "ERROR: cannot use both --every and --nth" << endl;
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  cout << endl
       << getVersion() << endl
       << "cwd " << cwd << endl
       << "cmdline " << cmdline << endl
       << "sysname  " << uts.sysname << endl
       << "hostname " << uts.nodename << endl
       << "machine  " << uts.machine << endl
       << "user     " << VERuser() << endl;

  if (!inDir.empty())
    cout << "Input directory: " << inDir << endl;
  if (!inTrkList.empty()) {
    cout << "Input .trk files:";
    for (vector<string>::const_iterator istr = inTrkList.begin();
                                        istr < inTrkList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (!inAscList.empty()) {
    cout << "Input text files:";
    for (vector<string>::const_iterator istr = inAscList.begin();
                                        istr < inAscList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (!outDir.empty())
    cout << "Output directory: " << outDir << endl;
  if (!outTrkList.empty()) {
    cout << "Output .trk files:";
    for (vector<string>::const_iterator istr = outTrkList.begin();
                                        istr < outTrkList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (!outAscList.empty()) {
    cout << "Output text files:";
    for (vector<string>::const_iterator istr = outAscList.begin();
                                        istr < outAscList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (!outVolList.empty()) {
    cout << "Output volumes:";
    for (vector<string>::const_iterator istr = outVolList.begin();
                                        istr < outVolList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (!overList.empty()) {
    cout << "Scalar overlay volumes:";
    for (vector<string>::const_iterator istr = overList.begin();
                                        istr < overList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (!incMaskList.empty()) {
    cout << "Inclusion mask volumes:";
    for (vector<string>::const_iterator istr = incMaskList.begin();
                                        istr < incMaskList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (!excMaskList.empty()) {
    cout << "Exclusion mask volumes:";
    for (vector<string>::const_iterator istr = excMaskList.begin();
                                        istr < excMaskList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (!incTermMaskList.empty()) {
    cout << "Terminal inclusion mask volumes:";
    for (vector<string>::const_iterator istr = incTermMaskList.begin();
                                        istr < incTermMaskList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (!excTermMaskList.empty()) {
    cout << "Terminal exclusion mask volumes:";
    for (vector<string>::const_iterator istr = excTermMaskList.begin();
                                        istr < excTermMaskList.end(); istr++)
      cout << " " << *istr;
    cout << endl;
  }
  if (lengthMin > -1)
    cout << "Lower length threshold: " << lengthMin << endl;
  if (lengthMax > -1)
    cout << "Upper length threshold: " << lengthMax << endl;
  if (!inRefFile.empty())
    cout << "Input reference: " << inRefFile << endl;
  if (!outRefFile.empty())
    cout << "Output reference: " << outRefFile << endl;
  if (!affineXfmFile.empty())
    cout << "Affine registration: " << affineXfmFile << endl;
  if (!nonlinXfmFile.empty())
    cout << "Nonlinear registration: " << nonlinXfmFile << endl;
  if (!(affineXfmFile.empty() && nonlinXfmFile.empty()))
    cout << "Invert registration: " << doInvXfm << endl;
  cout << "Fill gaps between points: " << doFill << endl;
  if (doNth)
    cout << "Keeping n-th streamline only: " << strNum << endl;
  else {
    if (doEvery)
      cout << "Keeping every n-th streamline: " << everyNum << endl;
    if (doMean)
      cout << "Keeping mean of streamlines" << endl;
    else if (doNearMean)
      cout << "Keeping the streamline nearest to the mean" << endl;
  }
  if (doSmooth)
    cout << "Smoothing streamlines" << endl;
  if (doMerge)
    cout << "Merging multiple inputs into a single output" << endl;

  return;
}

