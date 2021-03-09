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
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

int doInvNonlin = 0, doFill = 0, doMean = 0, doNth = 0, strNum = -1,
    lengthMin = -1, lengthMax = -1;
unsigned int nTract = 0;
std::string inDir, outDir, inRefFile, outRefFile, affineXfmFile, nonlinXfmFile;
vector<std::string> inTrkList, inAscList, outTrkList, outAscList, outVolList,
               incMaskList, excMaskList;
vector<MRI *>  incMask, excMask;

struct utsname uts;
char *cmdline, cwd[2000];

Timer cputimer;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs, cputime;
  char outorient[4];
  std::string fname;
  vector<float> point(3), step(3, 0);
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
  inref = MRIread(inRefFile.c_str());
  outref = MRIread(outRefFile.c_str());

  if (!outVolList.empty()) {
    outvol = MRIclone(outref, NULL);
  }

  // Output space orientation information
  outv2r = MRIgetVoxelToRasXform(outref);
  MRIdircosToOrientationString(outref, outorient);

  // Read transform files
#ifndef NO_CVS_UP_IN_HERE
  if (!nonlinXfmFile.empty()) {
    if (!affineXfmFile.empty()) {
      affinereg.ReadXfm(affineXfmFile.c_str(), inref, 0);
    }
    nonlinreg.ReadXfm(nonlinXfmFile.c_str(), outref);
  } else {
#endif
    if (!affineXfmFile.empty()) {
      affinereg.ReadXfm(affineXfmFile.c_str(), inref, outref);
    }
#ifndef NO_CVS_UP_IN_HERE
  }
#endif

  // Read inclusion masks
  for (auto imask = incMaskList.begin();
       imask < incMaskList.end(); imask++) {
    incMask.push_back(MRIread((*imask).c_str()));
  }

  // Read exclusion masks
  for (auto imask = excMaskList.begin();
       imask < excMaskList.end(); imask++) {
    excMask.push_back(MRIread((*imask).c_str()));
  }

  for (unsigned int itract = 0; itract < nTract; itract++) {
    int npts, nstr = 0;
    CTrackReader trkreader;
    TRACK_HEADER trkheadin;
    vector< vector<float> > streamlines;

    cout << "Processing input file " << itract+1 << " of " << nTract
         << "..." << endl;
    cputimer.reset();

    if (!inTrkList.empty()) {		// Read streamlines from .trk file
      if (!inDir.empty()) {
	fname = inDir + "/" + inTrkList.at(itract);
      } else {
	fname = inTrkList.at(itract);
      }

      if (!trkreader.Open(fname.c_str(), &trkheadin)) {
        cout << "ERROR: Cannot open input file " << fname << endl;
        cout << "ERROR: " << trkreader.GetLastErrorMessage() << endl;
        exit(1);
      }

      while (trkreader.GetNextPointCount(&npts)) {
        const int veclen = npts*3;
        float *iraw, *rawpts = new float[veclen];
        vector<float> newpts(veclen);

        // Read a streamline from input file
        trkreader.GetNextTrackData(npts, rawpts);

        if ( (doNth && nstr != strNum) ||
             (lengthMin > -1 && npts <= lengthMin) ||
             (lengthMax > -1 && npts >= lengthMax) ) {
          delete[] rawpts;
          nstr++;
          continue;
        }

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

        nstr++;
      }
    }
    else if (!inAscList.empty()) {	// Read streamlines from text file
      string ptline;
      ifstream infile;
      vector<float> newpts;

      if (!inDir.empty()) {
	fname = inDir + '/' + inAscList.at(itract);
      } else {
	fname = inAscList.at(itract);
      }

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
          if ( (!doNth || nstr == strNum) &&
               (lengthMin == -1 || (int) newpts.size()/3 > lengthMin) &&
               (lengthMax == -1 || (int) newpts.size()/3 < lengthMax) )
            streamlines.push_back(newpts);

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

    nstr = streamlines.size();

    for (int kstr = nstr-1; kstr >= 0; kstr--) {
      vector<float> newpts;

      for (vector<float>::iterator ipt = streamlines[kstr].begin();
                                   ipt < streamlines[kstr].end(); ipt += 3) {
        copy(ipt, ipt+3, point.begin());

        // Apply affine transform
        if (!affinereg.IsEmpty())
          affinereg.ApplyXfm(point, point.begin());

#ifndef NO_CVS_UP_IN_HERE
        // Apply nonlinear transform
        if (!nonlinreg.IsEmpty()) {
          if (doInvNonlin)
            nonlinreg.ApplyXfmInv(point, point.begin());
          else
            nonlinreg.ApplyXfm(point, point.begin());
        }
#endif

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

            step[k] = dist;
            dist = fabs(dist);

            if (dist > dmax)
              dmax = dist;
          }

          if (dmax > 0)
            for (int k = 0; k < 3; k++)
              step[k] /= dmax;
        }

        copy(ipt, ipt+3, point.begin());

        for (int istep = (int) round(dmax); istep > 0; istep--) {
          newpts.insert(newpts.end(), point.begin(), point.end());

          for (int k = 0; k < 3; k++)
            point[k] += step[k];
        }
      }

      streamlines[kstr].resize(newpts.size());
      copy(newpts.begin(), newpts.end(), streamlines[kstr].begin());
    }

    // Apply inclusion masks
    nstr = streamlines.size();

    for (int kstr = nstr-1; kstr >= 0; kstr--) {
      bool dokeep = true;

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

      if (!dokeep)
        streamlines.erase(streamlines.begin() + kstr);
    }

    if (doMean && !streamlines.empty()) {
      unsigned int nstr = streamlines.size(), lmin, kmax, nstrout, kstrmean = 0;
      float dmin = numeric_limits<float>::infinity();
      vector<bool> isout(nstr);
      vector<unsigned int> lengths(nstr);
      vector<float> steps(nstr), strmean, strstd, strU, strL;
      vector<bool>::iterator iout;
      vector<unsigned int>::iterator ilen;
      vector<float>::iterator istep, imean, istd, iupper, ilower;

      // Find the minimum streamline length
      ilen = lengths.begin();

      for (vector< vector<float> >::const_iterator istr = streamlines.begin();
                                                   istr < streamlines.end();
                                                   istr++) {
        *ilen = istr->size()/3;
        ilen++;
      }

      lmin = *min_element(lengths.begin(), lengths.end());
      kmax = lmin - 1;

      // Set the step size for each streamline to have equal number of steps
      istep = steps.begin();

      for (ilen = lengths.begin(); ilen < lengths.end(); ilen++) {
        *istep = (*ilen - 1) / (float) kmax;
        istep++;
      }

      lengths.clear();

      // Compute mean and standard deviation of point coordinates at each step
      strmean.resize(lmin*3);
      fill(strmean.begin(), strmean.end(), 0.0);
      strstd.resize(strmean.size());
      fill(strstd.begin(), strstd.end(), 0.0);

      istep = steps.begin();

      for (vector< vector<float> >::const_iterator istr = streamlines.begin();
                                                   istr < streamlines.end();
                                                   istr++) {
        imean = strmean.begin();
        istd = strstd.begin();

        for (unsigned int kpt = 0; kpt < kmax; kpt++) {
          const unsigned int idx = (unsigned int) round(kpt * (*istep));
          vector<float>::const_iterator ipt = istr->begin() + idx * 3;

          if (ipt > istr->end() - 3)
            ipt = istr->end() - 3;

          for (int k = 0; k < 3; k++) {
            imean[k] += (float) ipt[k];
            istd[k]  += (float) ipt[k] * ipt[k];
          }

          imean += 3;
          istd += 3;
        }
 
        istep++;
      }

      // Compute upper and lower limit for flagging a point as an outlier
      strU.resize(strmean.size());
      strL.resize(strmean.size());

      istd = strstd.begin();
      iupper = strU.begin();
      ilower = strL.begin();

      for (imean = strmean.begin(); imean < strmean.end(); imean++) {
        float dout;

        *imean /= nstr;

        if (nstr > 1)
          *istd = sqrt((*istd - nstr * (*imean) * (*imean)) / (nstr-1));
        else
          *istd = 0;

        if (imean == strmean.begin() || imean == strmean.end() - 1)
          dout = *istd;
        else
          dout = 2 * (*istd);

        *iupper = *imean + dout;
        *ilower = *imean - dout;

        istd++;
        iupper++;
        ilower++;
      }

      // Flag streamlines with at least one outlier point
      fill(isout.begin(), isout.end(), false);

      iout = isout.begin();
      istep = steps.begin();

      for (vector< vector<float> >::const_iterator istr = streamlines.begin();
                                                   istr < streamlines.end();
                                                   istr++) {
        iupper = strU.begin();
        ilower = strL.begin();

        for (unsigned int kpt = 0; kpt < kmax; kpt++) {
          const unsigned int idx = (unsigned int) round(kpt * (*istep));
          vector<float>::const_iterator ipt = istr->begin() + idx * 3;

          for (int k = 0; k < 3; k++) {
            if (*ipt > *iupper || *ipt < *ilower) {
              *iout = true;
              break;
            }
        
            ipt++;
            iupper++;
            ilower++;
          }

          if (*iout)
            break;
        }

        iout++;
        istep++;
      }

      nstrout = count(isout.begin(), isout.end(), true);

      cout << "INFO: Found " << nstrout
           << " streamlines with at least one outlier coordinate"
           << endl;

      if (nstrout == nstr) {
        cout << "INFO: Turning off outlier checks" << endl;

        fill(isout.begin(), isout.end(), false);
      }

      // Find the non-outlier streamline that is closest to the mean streamline
      iout = isout.begin();
      istep = steps.begin();

      for (vector< vector<float> >::const_iterator istr = streamlines.begin();
                                                   istr < streamlines.end();
                                                   istr++) {
        if (!*iout) {
          float dist = 0;

          imean = strmean.begin();
          istd = strstd.begin();

          for (unsigned int kpt = 0; kpt < kmax; kpt++) {
            const unsigned int idx = (unsigned int) round(kpt * (*istep));
            vector<float>::const_iterator ipt = istr->begin() + idx * 3;

            const float dx = ipt[0] - imean[0],
                        dy = ipt[1] - imean[1],
                        dz = ipt[2] - imean[2];

            dist += sqrt(dx*dx + dy*dy + dz*dz);

            imean += 3;
            istd += 3;
          }

          if (dist < dmin) {
            dmin = dist;
            kstrmean = istr - streamlines.begin();
          }
        }

        iout++;
        istep++;
      }

      cout << "INFO: Streamline closest to the mean is " << kstrmean << endl;

      // Keep only the chosen streamline for writing to disk
      streamlines.erase(streamlines.begin(), streamlines.begin() + kstrmean);
      streamlines.erase(streamlines.begin() + 1, streamlines.end());
    }

    // Write transformed streamlines to volume
    if (!outVolList.empty()) {
      if (!outDir.empty()) {
	fname = outDir + "/" + outVolList.at(itract);
      } else {
	fname = outVolList.at(itract);
      }

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

      if (!outDir.empty()) {
	fname = outDir + "/" + outAscList.at(itract);
      } else {
	fname = outAscList.at(itract);
      }

      outfile.open(fname, ios::out);
      if (!outfile) {
        cout << "ERROR: Could not open " << fname << " for writing" << endl;
        exit(1);
      }

      for (auto istr = streamlines.begin();
	   istr < streamlines.end();
	   istr++) {
        for (auto ipt = istr->begin();
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

      // Set output .trk header
      if (inTrkList.empty()) {
        trkheadout.Initialize();

        trkheadout.origin[0] = 0;
        trkheadout.origin[1] = 0;
        trkheadout.origin[2] = 0;
      } else {
        trkheadout = trkheadin;
      }

      trkheadout.voxel_size[0] = outref->xsize;
      trkheadout.voxel_size[1] = outref->ysize;
      trkheadout.voxel_size[2] = outref->zsize;

      trkheadout.dim[0] = outref->width;
      trkheadout.dim[1] = outref->height;
      trkheadout.dim[2] = outref->depth;

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          trkheadout.vox_to_ras[i][j] = outv2r->rptr[i+1][j+1];
	}
      }

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

      trkheadout.n_count = (int) streamlines.size();

      // Open output .trk file
      if (!outDir.empty()) {
	fname = outDir + "/" + outTrkList.at(itract);
      } else {
	fname = outTrkList.at(itract);
      }

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

        trkwriter.WriteNextTrack(istr->size()/3, &(istr->at(0)));
      }

      trkwriter.Close();
    }

    cputime = cputimer.milliseconds();
    cout << "Done in " << cputime/1000.0 << " sec." << endl;
  }

  MatrixFree(&outv2r);
  MRIfree(&inref);
  MRIfree(&outref);
  if (!outVolList.empty())
    MRIfree(&outvol);
  for (vector<MRI *>::iterator imask = incMask.begin();
                               imask < incMask.end(); imask++)
    MRIfree(&(*imask));
  for (vector<MRI *>::iterator imask = excMask.begin();
                               imask < excMask.end(); imask++)
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
    else if (!strcasecmp(option, "--invnl"))
      doInvNonlin = 1;
    else if (!strcasecmp(option, "--fill"))
      doFill = 1;
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
    else if (!strcasecmp(option, "--nth")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &strNum);
      nargsused = 1;
      doNth = 1;
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
  << "     Output .trk file(s), as many as inputs" << endl
  << "   --outasc <file> [...]:" << endl
  << "     Output ASCII plain text file(s), as many as inputs" << endl
  << "   --outvol <file> [...]:" << endl
  << "     Output volume(s), as many as inputs" << endl
  << "   --outdir <dir>:" << endl
  << "     Output directory (optional)" << endl
  << "     If specified, names of output .trk files and volumes are relative to this)" << endl
  << "   --inref <file>:" << endl
  << "     Input reference volume" << endl
  << "   --outref <file>:" << endl
  << "     Output reference volume" << endl
  << "   --reg <file>:" << endl
  << "     Affine registration (.mat), applied first" << endl
  << "   --regnl <file>:" << endl
  << "     Nonlinear registration (.m3z), applied second" << endl
  << "   --invnl:" << endl
  << "     Apply inverse of nonlinear warp (with --regnl, default: no)" << endl
  << "   --fill:" << endl
  << "     Fill gaps b/w mapped points by linear interpolation" << endl
  << "     (Default: don't fill)" << endl
  << "   --imask <file> [...]:" << endl
  << "     Inclusion mask(s), applied to all input .trk files" << endl
  << "   --emask <file> [...]:" << endl
  << "     Exclusion mask(s), applied to all input .trk files" << endl
  << "   --lmin <num>:" << endl
  << "     Only save streamlines with length greater than this number" << endl
  << "   --lmax <num>:" << endl
  << "     Only save streamlines with length smaller than this number" << endl
  << "   --mean:" << endl
  << "     Only save the mean streamline (Default: save all)" << endl
  << "   --nth <num>:" << endl
  << "     Only save the n-th (0-based) streamline (Default: save all)" << endl
  << endl
  << "Other options" << endl
  << "   --debug:     turn on debugging" << endl
  << "   --checkopts: don't run anything, just check options and exit" << endl
  << "   --help:      print out information on how to use this program" << endl
  << "   --version:   print out version and exit" << endl
  << endl
  << "Order of operations (all optional):" << endl
  << "   1. Keep n-th streamline only" << endl
  << "   2. Apply streamline length threshold(s)" << endl
  << "   3. Apply affine transform" << endl
  << "   4. Apply non-linear transform" << endl
  << "   5. Apply inclusion mask(s)" << endl
  << "   6. Apply exclusion mask(s)" << endl
  << "   7. Find mean streamline" << endl
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

  if (outTrkList.empty() && outAscList.empty() && outVolList.empty()) {
    cout << "ERROR: must specify output .trk or text or volume file(s)" << endl;
    exit(1);
  }
  if (!outTrkList.empty() && outTrkList.size() != nTract) {
    cout << "ERROR: must specify as many output .trk files as input files"
         << endl;
    exit(1);
  }
  if (!outAscList.empty() && outAscList.size() != nTract) {
    cout << "ERROR: must specify as many output text files as input .trk files"
         << endl;
    exit(1);
  }
  if (!outVolList.empty() && outVolList.size() != nTract) {
    cout << "ERROR: must specify as many output volumes as input .trk files"
         << endl;
    exit(1);
  }
  if (inRefFile.empty()) {
    cout << "ERROR: must specify input reference volume" << endl;
    exit(1);
  }
  if (outRefFile.empty()) {
    cout << "ERROR: must specify output reference volume" << endl;
    exit(1);
  }
  if (doMean && doNth) {
    cout << "ERROR: cannot use both --mean and --nth" << endl;
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

  if (!inDir.empty()) {
    cout << "Input directory: " << inDir << endl;
  }
  if (!inTrkList.empty()) {
    cout << "Input .trk files:";
    for (auto istr = inTrkList.begin(); istr < inTrkList.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }
  if (!inAscList.empty()) {
    cout << "Input text files:";
    for (auto istr = inAscList.begin(); istr < inAscList.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }
  if (!outDir.empty()) {
    cout << "Output directory: " << outDir << endl;
  }
  if (!outTrkList.empty()) {
    cout << "Output .trk files:";
    for (auto istr = outTrkList.begin(); istr < outTrkList.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }
  if (!outAscList.empty()) {
    cout << "Output text files:";
    for (auto istr = outAscList.begin(); istr < outAscList.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }
  if (!outVolList.empty()) {
    cout << "Output volumes:";
    for (auto istr = outVolList.begin(); istr < outVolList.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }
  if (!incMaskList.empty()) {
    cout << "Inclusion mask volumes:";
    for (auto istr = incMaskList.begin(); istr < incMaskList.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }
  if (!excMaskList.empty()) {
    cout << "Exclusion mask volumes:";
    for (auto istr = excMaskList.begin(); istr < excMaskList.end(); istr++) {
      cout << " " << *istr;
    }
    cout << endl;
  }
  if (lengthMin > -1)
    cout << "Lower length threshold: " << lengthMin << endl;
  if (lengthMax > -1)
    cout << "Upper length threshold: " << lengthMax << endl;
  cout << "Input reference: " << inRefFile << endl;
  cout << "Output reference: " << outRefFile << endl;
  if (!affineXfmFile.empty()) {
    cout << "Affine registration: " << affineXfmFile << endl;
  }
  if (!nonlinXfmFile.empty()) {
    cout << "Nonlinear registration: " << nonlinXfmFile << endl;
    cout << "Invert nonlinear morph: " << doInvNonlin << endl;
  }
  cout << "Fill gaps between points: " << doFill << endl;
  if (doMean) {
    cout << "Saving mean streamline" << endl;
  } else if (doNth) {
    cout << "Saving single streamline: " << strNum << endl;
  } else {
    cout << "Saving all streamlines" << endl;
  }

  return;
}

