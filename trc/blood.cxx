/**
 * @file  blood.cxx
 * @brief Training data and methods that probabilistic tractography feeds on
 *
 * Training data and methods that probabilistic tractography feeds on
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
 *
 * Copyright (C) 2010,
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

#include <blood.h>

using namespace std;

const int Blood::mDistThresh = 2;
const unsigned int Blood::mCurvOffset = 2;
const float Blood::mLengthCutoff = 0.05,
            Blood::mLengthRatio = 3.0;

Blood::Blood(const char *TrainListFile, const char *TrainTrkFile,
             const char *TrainAsegFile, const char *TrainMaskFile,
             const char *TestMaskFile) {
  int dirs[45] = { 0,  0,  0,
                   1,  0,  0,
                  -1,  0,  0,
                   0,  1,  0,
                   0, -1,  0,
                   0,  0,  1,
                   0,  0, -1,
                   1,  1,  1,
                  -1,  1,  1,
                   1, -1,  1,
                  -1, -1,  1,
                   1,  1, -1,
                  -1,  1, -1,
                   1, -1, -1,
                  -1, -1, -1 };

  // Directions in which to look for neighboring labels
  mNumLocal = 15;
  mDirLocal.insert(mDirLocal.begin(), dirs, dirs+45);

  mNumNear = 14;
  mDirNear.insert(mDirNear.begin(), dirs+3, dirs+45);
  
  cout << "Loading brain mask of output subject from " << TestMaskFile << endl;
  mTestMask = MRIread(TestMaskFile);

  // Size of volumes in common space
  mNx = mTestMask->width;
  mNy = mTestMask->height;
  mNz = mTestMask->depth;

  // Read inputs for first pathway
  ReadStreamlines(TrainListFile, TrainTrkFile);
  ReadAnatomy(TrainListFile, TrainAsegFile, TrainMaskFile);
}

Blood::Blood(const char *TrainTrkFile) {
  mTestMask = 0;
  ReadStreamlines(TrainTrkFile);
}

Blood::~Blood() {
  vector<MRI *>::iterator ivol;

  for (ivol = mAseg.begin(); ivol != mAseg.end(); ivol++)
    MRIfree(&(*ivol));

  for (ivol = mMask.begin(); ivol != mMask.end(); ivol++)
    MRIfree(&(*ivol));

  MRIfree(&mTestMask);
}

//
// Read streamlines of training subjects
//
void Blood::ReadStreamlines(const char *TrainListFile,
                            const char *TrainTrkFile) {
  string subjdir;
  ifstream infile(TrainListFile, ios::in);

  if (!infile) {
    cout << "ERROR: Could not open " << TrainListFile << endl;
    exit(1);
  }

  // Clear all variables that depend on pathway
  mStreamlines.clear();
  mLengths.clear();
  mNumLines.clear();
  mIsInEnd1.clear();
  mIsInEnd2.clear();
  mCenterStreamline.clear();
  mControlPoints.clear();
  mControlStd.clear();
  mTangent.clear();
  mTangentArc.clear();
  mTangentMean.clear();
  mTangentStd.clear();
  mCurvature.clear();
  mCurvatureArc.clear();
  mCurvatureMean.clear();
  mCurvatureStd.clear();
  mAsegLocal.clear();
  mAsegLocalArc.clear();
  mHistoLocal.clear();
  mPriorLocal.clear();
  mIdsLocal.clear();
  mAsegNear.clear();
  mAsegNearArc.clear();
  mHistoNear.clear();
  mPriorNear.clear();
  mIdsNear.clear();
  mAsegDist.clear();
  mAsegDistArc.clear();
  mAsegDistMean.clear();
  mAsegDistStd.clear();

  while (infile >> subjdir) {
    CTrackReader trkreader;
    TRACK_HEADER trkheader;
    string fname = subjdir + "/" + TrainTrkFile;
    int npts, nlines = 0;

    if (!trkreader.Open(fname.c_str(), &trkheader)) {
      cout << "WARN: Could not open " << fname << " for reading" << endl;
      cout << "WARN: Error was: " << trkreader.GetLastErrorMessage() << endl;
      cout << "WARN: Skipping to next subject" << endl;
      mNumLines.push_back(0);
      continue;
    }

    cout << "Loading streamlines from " << fname << endl;

    while (trkreader.GetNextPointCount(&npts)) {
      bool isinmask = true;
      vector<int> pts;
      float *iraw, *rawpts = new float[npts*3];

      trkreader.GetNextTrackData(npts, rawpts);

      // Divide by voxel size and round to get voxel coords
      iraw = rawpts; 
      for (int ipt = npts; ipt > 0; ipt--)
        for (int k = 0; k < 3; k++) {
          *iraw = round(*iraw / trkheader.voxel_size[k]);
          iraw++;
        }

      // Remove duplicate points and check that streamline stays within mask
      iraw = rawpts; 
      for (int ipt = npts; ipt > 1; ipt--)
        if ( (iraw[0] != iraw[3]) || (iraw[1] != iraw[4]) ||
                                     (iraw[2] != iraw[5]) ) {
          for (int k = 0; k < 3; k++) {
            pts.push_back((int) *iraw);
            iraw++;
          }

          if (!IsInMask(pts.end()-3)) {
            isinmask = false;
            break;
          }
        }
        else
          iraw += 3;

      // Don't save streamline if it strayed off the test subject's mask
      if (!isinmask) {
        delete[] rawpts;
        continue;
      }

      // Final point
      for (int k = 0; k < 3; k++) {
        pts.push_back((int) *iraw);
        iraw++;
      }

      // Don't save streamline if it strayed off the test subject's mask
      if (!IsInMask(pts.end()-3)) {
        delete[] rawpts;
        continue;
      }

      mStreamlines.push_back(pts);
      mLengths.push_back(pts.size() / 3);
      nlines++;

      delete[] rawpts;
    }

    mNumLines.push_back(nlines);

    trkreader.Close();
  }

  mNumTrain = mNumLines.size();

  if (!mMask.empty() && ((int) mMask.size() != mNumTrain)) {
    cout << "ERROR: Number of training paths and cortex masks do not match"
         << endl;
    exit(1);
  }
}

//
// Read streamlines from single file
//
void Blood::ReadStreamlines(const char *TrainTrkFile) {
  CTrackReader trkreader;
  TRACK_HEADER trkheader;
  int npts, nlines = 0;

  // Clear all variables that depend on pathway
  mStreamlines.clear();
  mLengths.clear();
  mNumLines.clear();
  mIsInEnd1.clear();
  mIsInEnd2.clear();
  mCenterStreamline.clear();
  mControlPoints.clear();
  mControlStd.clear();
  mTangent.clear();
  mTangentArc.clear();
  mTangentMean.clear();
  mTangentStd.clear();
  mCurvature.clear();
  mCurvatureArc.clear();
  mCurvatureMean.clear();
  mCurvatureStd.clear();
  mAsegLocal.clear();
  mAsegLocalArc.clear();
  mHistoLocal.clear();
  mPriorLocal.clear();
  mIdsLocal.clear();
  mAsegNear.clear();
  mAsegNearArc.clear();
  mHistoNear.clear();
  mPriorNear.clear();
  mIdsNear.clear();
  mAsegDist.clear();
  mAsegDistArc.clear();
  mAsegDistMean.clear();
  mAsegDistStd.clear();

  if (!trkreader.Open(TrainTrkFile, &trkheader)) {
    cout << "ERROR: Could not open " << TrainTrkFile << " for reading" << endl;
    exit(1);
  }

  cout << "Loading streamlines from " << TrainTrkFile << endl;

  while (trkreader.GetNextPointCount(&npts)) {
    bool isinmask = true;
    vector<int> pts;
    float *iraw, *rawpts = new float[npts*3];

    trkreader.GetNextTrackData(npts, rawpts);

    // Divide by voxel size and round to get voxel coords
    iraw = rawpts; 
    for (int ipt = npts; ipt > 0; ipt--)
      for (int k = 0; k < 3; k++) {
        *iraw = round(*iraw / trkheader.voxel_size[k]);
        iraw++;
      }

    // Remove duplicate points and check that streamline stays within mask
    iraw = rawpts; 
    for (int ipt = npts; ipt > 1; ipt--)
      if ( (iraw[0] != iraw[3]) || (iraw[1] != iraw[4]) ||
                                   (iraw[2] != iraw[5]) ) {
        for (int k = 0; k < 3; k++) {
          pts.push_back((int) *iraw);
          iraw++;
        }

        if (!IsInMask(pts.end()-3)) {
          isinmask = false;
          break;
        }
      }
      else
        iraw += 3;

    // Don't save streamline if it strayed off the test subject's mask
    if (!isinmask) {
      delete[] rawpts;
      continue;
    }

    // Final point
    for (int k = 0; k < 3; k++) {
      pts.push_back((int) *iraw);
      iraw++;
    }

    // Don't save streamline if it strayed off the test subject's mask
    if (!IsInMask(pts.end()-3)) {
      delete[] rawpts;
      continue;
    }

    mStreamlines.push_back(pts);
    mLengths.push_back(pts.size() / 3);
    nlines++;

    delete[] rawpts;
  }

  mNumLines.push_back(nlines);

  trkreader.Close();
}

//
// Check that a point is inside the test subject's brain mask
//
bool Blood::IsInMask(vector<int>::const_iterator Point) {
  if (mTestMask)
    return (Point[0] > -1) && (Point[0] < mNx) &&
           (Point[1] > -1) && (Point[1] < mNy) &&
           (Point[2] > -1) && (Point[2] < mNz) &&
           (MRIgetVoxVal(mTestMask, Point[0], Point[1], Point[2], 0) > 0);
  else
    return true;
}

//
// Read segmentations and cortical masks of training subjects
//
void Blood::ReadAnatomy(const char *TrainListFile, const char *TrainAsegFile,
                                                   const char *TrainMaskFile) {
  string subjdir;
  vector<MRI *>::iterator ivol;
  ifstream infile(TrainListFile, ios::in);

  if (!infile) {
    cout << "ERROR: Could not open " << TrainListFile << endl;
    exit(1);
  }

  for (ivol = mAseg.begin(); ivol != mAseg.end(); ivol++)
    MRIfree(&(*ivol));
  mAseg.clear();  

  for (ivol = mMask.begin(); ivol != mMask.end(); ivol++)
    MRIfree(&(*ivol));
  mMask.clear();  

  // Clear all variables that depend on underlying anatomy
  mAsegLocal.clear();
  mAsegLocalArc.clear();
  mHistoLocal.clear();
  mPriorLocal.clear();
  mIdsLocal.clear();
  mAsegNear.clear();
  mAsegNearArc.clear();
  mHistoNear.clear();
  mPriorNear.clear();
  mIdsNear.clear();
  mAsegDist.clear();
  mAsegDistArc.clear();
  mAsegDistMean.clear();
  mAsegDistStd.clear();

  while (infile >> subjdir) {
    string fname;

    fname = subjdir + "/" + TrainAsegFile;
    cout << "Loading segmentation from " << fname << endl;
    mAseg.push_back(MRIread(fname.c_str()));

    fname = subjdir + "/" + TrainMaskFile;
    cout << "Loading cortex mask from " << fname << endl;
    mMask.push_back(MRIread(fname.c_str()));
  }

  mNumTrain = mMask.size();

  if (!mNumLines.empty() && ((int) mNumLines.size() != mNumTrain)) {
    cout << "ERROR: Number of training paths and training masks do not match"
         << endl;
    exit(1);
  }
}

//
// Do main computation for priors
//
void Blood::ComputePriors() {
  cout << "Found a total of " << mStreamlines.size() << " streamlines" << endl;

  cout << "Matching streamline ends" << endl;
  MatchStreamlineEnds();
  cout << "Found " << mNumMask << " streamlines with ends in mask"  << endl;

  SetArcSegments();
  cout << "Splitting streamlines in " << mNumArc << " segments"  << endl;

  cout << "Computing prior on underlying anatomy" << endl;
  ComputeAnatomyPrior();

  cout << "Computing prior on curvature" << endl;
  ComputeCurvaturePrior();

  cout << "Finding center streamline" << endl;
  FindCenterStreamline();
}

//
// Select initial control points from center streamline
//
void Blood::SelectControlPoints(int NumControls) {
  cout << "Selecting " << NumControls << " points on center streamline" << endl;
  FindPointsOnStreamline(mCenterStreamline, NumControls);
}

//
// Make sure all streamline start/end points are consistent
// and save those that are within (a small distance of) the cortical mask
//
void Blood::MatchStreamlineEnds() {
  int imin;
  vector<bool>::iterator ivalid1, ivalid2;
  vector<int>::const_iterator ishort1, ishort2;
  vector< vector<int> >::iterator istr;
  vector<MRI *>::const_iterator imask;

  mIsInEnd1.resize(mStreamlines.size());
  mIsInEnd2.resize(mStreamlines.size());

  // Find the shortest streamline
  imin = min_element(mLengths.begin(), mLengths.end()) - mLengths.begin();
  ishort1 = mStreamlines[imin].begin();
  ishort2 = mStreamlines[imin].end()-3;

  // Assign streamline starts/ends to match start/end of shortest streamline
  // (As a first pass, this will align most but not all streamlines correctly)
  for (istr = mStreamlines.begin(); istr != mStreamlines.end(); istr++)
    if ( pow(istr->at(0) - ishort1[0], 2) +
         pow(istr->at(1) - ishort1[1], 2) +
         pow(istr->at(2) - ishort1[2], 2) >
         pow(istr->at(0) - ishort2[0], 2) +
         pow(istr->at(1) - ishort2[1], 2) +
         pow(istr->at(2) - ishort2[2], 2) )
      FlipStreamline(istr);

  // Do a second pass, aligning all streamlines with the majority
  mNumMask = 0;
  imask = mMask.begin();
  istr = mStreamlines.begin();
  ivalid1 = mIsInEnd1.begin();
  ivalid2 = mIsInEnd2.begin();

  for (vector<int>::const_iterator inum = mNumLines.begin();
                                   inum != mNumLines.end(); inum++) {
    for (int k = *inum; k > 0; k--) {
      double dist1 = 0, dist2 = 0;
      vector<int>::const_iterator iend1 = istr->begin(),
                                  iend2 = istr->end() - 3;
      vector< vector<int> >::const_iterator jstr;

      for (jstr = mStreamlines.begin(); jstr != mStreamlines.end(); jstr++)
        if (jstr != istr) {
          vector<int>::const_iterator jend1 = jstr->begin(),
                                      jend2 = jstr->end() - 3;

          // Is the start of path i closer to the start or the end of path j?
          dist1 += sqrt( pow(iend1[0] - jend2[0], 2) +
                         pow(iend1[1] - jend2[1], 2) +
                         pow(iend1[2] - jend2[2], 2) )
                 - sqrt( pow(iend1[0] - jend1[0], 2) +
                         pow(iend1[1] - jend1[1], 2) +
                         pow(iend1[2] - jend1[2], 2) );

          // Is the end of path i closer to the start or the end of path j?
          dist2 += sqrt( pow(iend2[0] - jend1[0], 2) +
                         pow(iend2[1] - jend1[1], 2) +
                         pow(iend2[2] - jend1[2], 2) )
                 - sqrt( pow(iend2[0] - jend2[0], 2) +
                         pow(iend2[1] - jend2[1], 2) +
                         pow(iend2[2] - jend2[2], 2) );
        }

      if ( (dist1 < 0 && dist2 < 0) ||
           (dist1 < 0 && dist2 > 0 && dist2 < -dist1) ||
           (dist1 > 0 && dist2 < 0 && dist1 < -dist2) )
        FlipStreamline(istr);

      if (dist1 < 0 && dist2 > 0)
        *ivalid1 = false;
      else
        *ivalid1 = IsValidEnd1(istr, *imask);

      if (dist1 > 0 && dist2 < 0)
        *ivalid2 = false;
      else
        *ivalid2 = IsValidEnd2(istr, *imask);

      if (*ivalid1 && *ivalid2)
        mNumMask++;

      istr++;
      ivalid1++;
      ivalid2++;
    }

    imask++;
  }
}

//
// Invert the order of points in a streamline
//
void Blood::FlipStreamline(vector< vector<int> >::iterator Streamline) {
  vector<int>::iterator itop = Streamline->begin(),
                        ibottom = Streamline->end() - 3;

  while (itop < ibottom) {
    for (int k = 0; k < 3; k++) {
      const int tmp = itop[k];
      itop[k] = ibottom[k];
      ibottom[k] = tmp;
    }

    itop += 3;
    ibottom -= 3;
  }
}

//
// Check if a streamline start point is within (a small distance of) a mask
//
bool Blood::IsValidEnd1(vector< vector<int> >::iterator Streamline, MRI *Mask) {
  vector<int>::iterator itop = Streamline->begin();

  if (MRIgetVoxVal(Mask, itop[0], itop[1], itop[2], 0) > 0)
    return true;
  else {
    // Extend the streamline by a few voxels if that gets it inside the mask
    int newpt[3] = { itop[0], itop[1], itop[2] };
    int diff[3]  = { itop[0] - itop[3],
                     itop[1] - itop[4],
                     itop[2] - itop[5] };
    vector<int> extend;

    for (int d = mDistThresh; d > 0; d--) {
      for (int k = 0; k < 3; k++)
        newpt[k] += diff[k];

      extend.insert(extend.begin(), newpt, newpt+3);

      if (MRIgetVoxVal(Mask, newpt[0], newpt[1], newpt[2], 0) > 0) {
        Streamline->insert(Streamline->begin(), extend.begin(), extend.end());
        mLengths[Streamline - mStreamlines.begin()] = Streamline->size() / 3;
        return true; 
      }
    }
  }

  return false;
}

//
// Check if a streamline end point is within (a small distance of) a mask
//
bool Blood::IsValidEnd2(vector< vector<int> >::iterator Streamline, MRI *Mask) {
  vector<int>::iterator ibottom = Streamline->end() - 3;

  if (MRIgetVoxVal(Mask, ibottom[0], ibottom[1], ibottom[2], 0) > 0)
    return true;
  else {
    // Extend the streamline by a few voxels if that gets it inside the mask
    int newpt[3] = { ibottom[0], ibottom[1], ibottom[2] };
    int diff[3] = { ibottom[0] - ibottom[-3],
                    ibottom[1] - ibottom[-2],
                    ibottom[2] - ibottom[-1] };
    vector<int> extend;

    for (int d = mDistThresh; d > 0; d--) {
      for (int k = 0; k < 3; k++)
        newpt[k] += diff[k];

      extend.insert(extend.end(), newpt, newpt+3);

      if (MRIgetVoxVal(Mask, newpt[0], newpt[1], newpt[2], 0) > 0) {
        Streamline->insert(Streamline->end(), extend.begin(), extend.end());
        mLengths[Streamline - mStreamlines.begin()] = Streamline->size() / 3;
        return true;
      }
    }
  }

  return false;
}

//
// Determine number of arc length segments for priors
//
void Blood::SetArcSegments() {
  int cutoff, count = 0;
  vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                               ivalid2 = mIsInEnd2.begin();
  vector<int>::const_iterator ilen = mLengths.begin();
  vector<int> lengths;

  // Find lengths of streamlines with valid end points
  for (vector<int>::const_iterator ilen = mLengths.begin();
                                   ilen != mLengths.end(); ilen++) {
    if (*ivalid1 && *ivalid2)
      lengths.push_back(*ilen);

    ivalid1++;
    ivalid2++;
  }

  // Find shortest streamline length above a cutoff for outliers
  cutoff = (int) round(lengths.size() * mLengthCutoff);
  
  sort(lengths.begin(), lengths.end());

  ilen = lengths.begin();
  while (count < cutoff) {
    ilen++;
    count++;
  }

  // Set the number of arc segments to a fraction of the shortest length
  mNumArc = (int) round(*ilen / mLengthRatio);
/*
  mNumArc = (int) round(*min_element(lengths.begin(), lengths.end())
                        / mLengthRatio);
*/
}

//
// Compute prior on underlying anatomy by streamline arc length
//
void Blood::ComputeAnatomyPrior() {
  vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                               ivalid2 = mIsInEnd2.begin();
  vector<int>::const_iterator ilen = mLengths.begin();
  vector< vector<int> >::const_iterator istr = mStreamlines.begin();
  vector<MRI *>::const_iterator iaseg = mAseg.begin();

  mAsegLocalArc.resize(mNumLocal * mNumArc);
  mAsegNearArc.resize(mNumNear * mNumArc);
  mAsegDistArc.resize(mNumNear * mNumArc);

  // Find anatomical labels around each point on each streamline
  for (vector<int>::const_iterator inum = mNumLines.begin();
                                   inum != mNumLines.end(); inum++) {
    for (int k = *inum; k > 0; k--) {
      if (*ivalid1 && *ivalid2) {	 // TODO: Use other streamlines, too?
        unsigned int ilocal = 0, inear = 0;
        const double darc = mNumArc / (double) *ilen;
        double larc = 0;
        vector<unsigned int> seglocal, segnear;	//TODO: this is not needed
        vector<int> segdist;

        for (vector<int>::const_iterator ipt = istr->begin();
                                         ipt != istr->end(); ipt += 3) {
          const int ix0 = ipt[0], iy0 = ipt[1], iz0 = ipt[2];
          const float seg0 = MRIgetVoxVal(*iaseg, ix0, iy0, iz0, 0);

          // Save local neighbor labels
          for (vector<int>::const_iterator idir = mDirLocal.begin();
                                           idir != mDirLocal.end(); idir += 3) {
            const int ix = ix0 + idir[0],
                      iy = iy0 + idir[1],
                      iz = iz0 + idir[2];

            seglocal.push_back((unsigned int) MRIgetVoxVal(*iaseg,
                               ((ix > -1 && ix < mNx) ? ix : ix0),
                               ((iy > -1 && iy < mNy) ? iy : iy0),
                               ((iz > -1 && iz < mNz) ? iz : iz0),
                               0));
            mAsegLocalArc[ilocal].push_back((unsigned int) MRIgetVoxVal(*iaseg, 
                                            ((ix > -1 && ix < mNx) ? ix : ix0),
                                            ((iy > -1 && iy < mNy) ? iy : iy0),
                                            ((iz > -1 && iz < mNz) ? iz : iz0),
                                            0));
            ilocal++;
          }

          // Save nearest neighbor labels
          for (vector<int>::const_iterator idir = mDirNear.begin();
                                           idir != mDirNear.end(); idir += 3) {
            int dist = 0, ix = ix0 + idir[0],
                          iy = iy0 + idir[1],
                          iz = iz0 + idir[2];
            float seg = seg0;

            while ((ix > -1) && (ix < mNx) &&
                   (iy > -1) && (iy < mNy) && 
                   (iz > -1) && (iz < mNz) && (seg == seg0)) { 
              seg = MRIgetVoxVal(*iaseg, ix, iy, iz, 0);
              dist++;

              ix += idir[0];
              iy += idir[1];
              iz += idir[2];
            }

            segnear.push_back((unsigned int) seg);
            segdist.push_back(dist);
            mAsegNearArc[inear].push_back((unsigned int) seg);
            mAsegDistArc[inear].push_back(dist);

            inear++;
          }

          larc += darc;

          if (larc > 1)		// Move to the next segment
            larc -= 1;
          else {		// Stay in the same segment
            ilocal -= mNumLocal;
            inear -= mNumNear;
          }
        }

        mAsegLocal.push_back(seglocal);
        mAsegNear.push_back(segnear);
        mAsegDist.push_back(segdist);
      }

      ivalid1++;
      ivalid2++;
      ilen++;
      istr++;
    }

    iaseg++;
  }

  // Compute priors on neighboring anatomical labels by arc length
  for (vector< vector<unsigned int> >::const_iterator
       iseg = mAsegLocalArc.begin(); iseg != mAsegLocalArc.end(); iseg++) {
    set<unsigned int> idlist(iseg->begin(), iseg->end());
    vector<int> histo(idlist.size());
    vector<float> prior(idlist.size() + 1);
    vector<int>::iterator ihisto = histo.begin();
    vector<float>::iterator iprior = prior.begin();
    const float denom = iseg->size() + idlist.size() + 1;

    for (set<unsigned int>::const_iterator iid = idlist.begin();
                                           iid != idlist.end(); iid++) {
      const unsigned int nmatch = count(iseg->begin(), iseg->end(), *iid);

      *ihisto = (int) nmatch;
      *iprior = -log((nmatch + 1) / denom);

      ihisto++;
      iprior++;
    }

    *iprior = -log(1 / denom);

    mIdsLocal.push_back(idlist);
    mHistoLocal.push_back(histo);
    mPriorLocal.push_back(prior);
  }

  vector< vector<int> >::const_iterator idist = mAsegDistArc.begin();

  for (vector< vector<unsigned int> >::const_iterator
       iseg = mAsegNearArc.begin(); iseg != mAsegNearArc.end(); iseg++) {
    set<unsigned int> idlist(iseg->begin(), iseg->end());
    vector<int> histo(idlist.size());
    vector<float> prior(idlist.size() + 1),
                  dmean(idlist.size(), 0),
                  dstd(idlist.size(), 0);
    vector<int>::iterator ihisto = histo.begin();
    vector<float>::iterator iprior = prior.begin(),
                            idmean = dmean.begin(),
                            idstd = dstd.begin();
    const float denom = iseg->size() + idlist.size() + 1;

    for (set<unsigned int>::const_iterator iid = idlist.begin();
                                           iid != idlist.end(); iid++) {
      const unsigned int nmatch = count(iseg->begin(), iseg->end(), *iid);
      vector<unsigned int>::const_iterator
                         imatch = find(iseg->begin(), iseg->end(), *iid);

      *ihisto = (int) nmatch;
      *iprior = -log((nmatch + 1) / denom);

      while (imatch < iseg->end()) {
        const int d = idist->at(imatch - iseg->begin());

        *idmean += (float) d;
        *idstd  += (float) d*d;

        imatch = find(imatch+1, iseg->end(), *iid);
      }

      *idmean /= nmatch;

      if (nmatch > 1)
        *idstd = sqrt((*idstd - nmatch * (*idmean) * (*idmean)) / (nmatch-1));
      else
        *idstd = 0;

      ihisto++;
      iprior++;
      idmean++;
      idstd++;
    }

    *iprior = -log(1 / denom);

    mIdsNear.push_back(idlist);
    mHistoNear.push_back(histo);
    mPriorNear.push_back(prior);
    mAsegDistMean.push_back(dmean);
    mAsegDistStd.push_back(dstd);

    idist++;
  }
}

//
// Compute prior on tangent vector and curvature by streamline arc length
//
void Blood::ComputeCurvaturePrior() {
  const unsigned int offset = 3 * mCurvOffset;
  float curv;
  vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                               ivalid2 = mIsInEnd2.begin();
  vector<int>::const_iterator ilen = mLengths.begin();
  vector<float> tang(3), norm(3), tangmean(3), tangstd(3);

  mTangentArc.resize(mNumArc);
  mCurvatureArc.resize(mNumArc);

  // Compute tangent vector and curvature along streamlines
  for (vector< vector<int> >::const_iterator istr = mStreamlines.begin();
                                           istr != mStreamlines.end(); istr++) {
    if (*ivalid1 && *ivalid2 && (*ilen >= 2 * (int) mCurvOffset + 1)) {
 					// TODO: Use other streamlines, too?
      unsigned int iarc = 0;
      const double darc = mNumArc / (double) *ilen;
      double larc = 0;
      vector<float> strtang, strcurv;	//TODO: this is not needed

      for (vector<int>::const_iterator ipt = istr->begin() + offset;
                                       ipt != istr->end() - offset;) {
        // Quadratic approximation of tangent and normal vectors from 3 points
        for (int k = 0; k < 3; k++) {
          tang[k] = (*(ipt + offset) - *(ipt - offset)) / 2.0;
          norm[k] = *(ipt + offset) - 2 * (*ipt) - *(ipt - offset);
          ipt++;
        }

        strtang.insert(strtang.end(), tang.begin(), tang.end());
        mTangentArc[iarc].insert(mTangentArc[iarc].end(),
                                 tang.begin(), tang.end());

        // Curvature = |r' x r''| / |r'|^3
        curv = 
          sqrt( ( pow(tang[1] * norm[2] - tang[2] * norm[1], 2) +
                  pow(tang[2] * norm[0] - tang[0] * norm[2], 2) +
                  pow(tang[0] * norm[1] - tang[1] * norm[0], 2) ) /
                pow(pow(tang[0], 2) + pow(tang[1], 2) + pow(tang[2], 2), 3) );

        strcurv.push_back(curv);
        mCurvatureArc[iarc].push_back(curv);

        larc += darc;

        if (larc > 1) {		// Move to the next segment
          larc -= 1;
          iarc++;
        }
      }

      mTangent.push_back(strtang);
      mCurvature.push_back(strcurv);
    }

    ivalid1++;
    ivalid2++;
    ilen++;
  }

  // Compute mean and variance of tangent vector by arc length
  for (vector< vector<float> >::const_iterator
       itang = mTangentArc.begin(); itang != mTangentArc.end(); itang++) {
    unsigned int nsamp = itang->size() / 3;

    fill(tangmean.begin(), tangmean.end(), 0.0);
    fill(tangstd.begin(), tangstd.end(), 0.0);

    for (vector<float>::const_iterator isamp = itang->begin(); 
                                       isamp != itang->end(); isamp += 3)
      for (int k = 0; k < 3; k++) {
        tangmean[k] += isamp[k];
        tangstd[k]  += isamp[k] * isamp[k];
      }

    for (int k = 0; k < 3; k++) {
      tangmean[k] /= nsamp;

      if (nsamp > 1)
        tangstd[k] = sqrt((tangstd[k] - nsamp * tangmean[k] * tangmean[k])
                           / (nsamp-1));
      else
        tangstd[k] = 0;
    }

    mTangentMean.insert(mTangentMean.end(), tangmean.begin(), tangmean.end());
    mTangentStd.insert(mTangentStd.end(), tangstd.begin(), tangstd.end());
  }

  // Compute mean and variance of curvature by arc length
  for (vector< vector<float> >::const_iterator
       icurv = mCurvatureArc.begin(); icurv != mCurvatureArc.end(); icurv++) {
    unsigned int nsamp = icurv->size();
    float curvmean = 0, curvstd = 0;

    for (vector<float>::const_iterator isamp = icurv->begin(); 
                                       isamp != icurv->end(); isamp++) {
      curvmean += *isamp;
      curvstd  += (*isamp) * (*isamp);
    }

    curvmean /= nsamp;

    if (nsamp > 1)
      curvstd = sqrt((curvstd - nsamp * curvmean * curvmean) / (nsamp-1));
    else
      curvstd = 0;

    mCurvatureMean.push_back(curvmean);
    mCurvatureStd.push_back(curvstd);
  }
}

//
// Find central streamline among streamlines with valid end points
//
void Blood::FindCenterStreamline() {
  double hdmin = numeric_limits<double>::infinity();
  vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                               ivalid2 = mIsInEnd2.begin();
  vector< vector<int> >::const_iterator icenter;

  for (vector< vector<int> >::const_iterator istr = mStreamlines.begin();
                                           istr != mStreamlines.end(); istr++) {
    if (*ivalid1 && *ivalid2) {
      double hdtot = 0;
      vector<bool>::const_iterator jvalid1 = mIsInEnd1.begin(),
                                   jvalid2 = mIsInEnd2.begin();
      vector<int>::const_iterator jlen = mLengths.begin();

      for (vector< vector<int> >::const_iterator jstr = mStreamlines.begin();
                                           jstr != mStreamlines.end(); jstr++) {
        double hd = 0;

        if (*jvalid1 && *jvalid2 && (jstr != istr)) {
          for (vector<int>::const_iterator jpt = jstr->begin();
                                           jpt != jstr->end(); jpt += 3) {
            double dmin = 1e6;

            for (vector<int>::const_iterator ipt = istr->begin();
                                             ipt != istr->end(); ipt += 3) {
              const double dist = sqrt( pow(ipt[0] - jpt[0], 2) +
                                        pow(ipt[1] - jpt[1], 2) +
                                        pow(ipt[2] - jpt[2], 2) );

              if (dist < dmin)
                dmin = dist;
            }

            hd += dmin;
          }

          hd /= (*jlen);
        }

        hdtot += hd;

        jvalid1++;
        jvalid2++;
        jlen++;
      }

      if (hdtot < hdmin) {
        hdmin = hdtot;
        icenter = istr;
      }
    }

    ivalid1++;
    ivalid2++;
  }

  mCenterStreamline = *icenter;
}

//
// Select a representative subset of uniformly spaced points on a streamline
//
void Blood::FindPointsOnStreamline(vector<int> &Streamline, int NumPoints) {
  bool sign1, sign2;
  const int nptot = Streamline.size() / 3;
  int kmax, nseg = 1, npthave = 1;
  double lentot = 0;
  vector<bool>::iterator iturn;
  vector<int>::const_iterator ipt;
  vector<int>::iterator inpt;
  vector<double>::iterator ilen, inptdec;
  vector<double>::const_iterator idist;
  vector<bool> isturnpt(nptot, true);
  vector<int> sum(3, 0), nptseg, cpts, diff(3, 0), diffmin(3, 0);
  vector<float> cstd;
  vector<double> sumsq(3, 0), ptdist(nptot-1, 0), lenseg, nptsegdec;

  if (NumPoints > (int) Streamline.size() / 3) {
    cout << "ERROR: Selected streamline has fewer than " << NumPoints
         << " points" << endl;
    exit(1);
  }

if (0) {
  // Find the dimension in which the points have the greatest variance
  for (ipt = Streamline.begin(); ipt != Streamline.end(); ipt += 3)
    for (int k = 0; k < 3; k++) {
      sum[k] += ipt[k];
      sumsq[k] += pow(ipt[k], 2);
    }

  for (int k = 0; k < 3; k++)
    sumsq[k] -= (pow(sum[k], 2) / nptot);

  kmax = max_element(sumsq.begin(), sumsq.end()) - sumsq.begin();
}
else {
  // Find the dimension that the points traverse the most
  vector<int> sumdiff(3, 0);
  for (ipt = Streamline.begin() + 3; ipt != Streamline.end(); ipt += 3)
    for (int k = 0; k < 3; k++)
      sumdiff[k] += abs(ipt[k] - ipt[k-3]);
  kmax = max_element(sumdiff.begin(), sumdiff.end()) - sumdiff.begin();
}

  // Find turning points along the dimension of greatest variance
  ipt = Streamline.begin() + kmax;
  sign2 = (*(ipt+3) > *ipt);
  ipt += 3;
  for (iturn = isturnpt.begin() + 1; iturn != isturnpt.end() - 1; iturn++) {
    sign1 = sign2;
    sign2 = (*(ipt+3) > *ipt);

    // Look for changes in direction and make sure adjacent points aren't turns
    if (sign1 == sign2 || *(iturn-1))
      *iturn = false;
    else
      nseg++;

    ipt += 3;
  }

  lenseg.resize(nseg);
  fill(lenseg.begin(), lenseg.end(), 0);
  nptseg.resize(nseg);
  fill(nptseg.begin(), nptseg.end(), 0);
  nptsegdec.resize(nseg);
  fill(nptsegdec.begin(), nptsegdec.end(), 0);

  ipt = Streamline.begin();
  iturn = isturnpt.begin() + 1;
  ilen = lenseg.begin();
  inpt = nptseg.begin();

  for (vector<double>::iterator idist = ptdist.begin();
                                idist != ptdist.end(); idist++) {
    // Compute distances between consecutive points
    for (int k = 0; k < 3; k++) {
      const int dpt = *(ipt+3) - *ipt;
      *idist += dpt*dpt;
      ipt++;
    }

    *idist = sqrt(*idist);

    // Sum to find total length and lengths of segments between turn points
    lentot += *idist;
    *ilen += *idist;
    (*inpt)++;

    if (*iturn) {
      ilen++;
      inpt++;
    }
    iturn++;
  }

  // Remove turn points if there are more than points requested
  while (nseg >= NumPoints) {
    int imerge,
        idel = min_element(lenseg.begin(), lenseg.end()) - lenseg.begin();

    if (idel == 0)
      imerge = 1;
    else if (idel == nseg-1)
      imerge = nseg-2;
    else
      imerge = (lenseg[idel-1] < lenseg[idel+1]) ? idel-1 : idel+1;

    lenseg[imerge] += lenseg[idel];
    lenseg.erase(lenseg.begin() + idel);

    nptseg[imerge] += nptseg[idel];
    nptseg.erase(nptseg.begin() + idel);

    iturn = isturnpt.begin() + 1;
    if (idel == 0)
      while (! *iturn)
        iturn++;
    while (idel > 0) {
      iturn++;
      if (*iturn) idel--;
    }
    *iturn = false;
   
    nseg--;
  }

  // Determine how many points to choose over each segment between turn points
  // based on length of segments
  inpt = nptseg.begin();
  inptdec = nptsegdec.begin();

  for (vector<double>::const_iterator ilen = lenseg.begin();
                                      ilen != lenseg.end(); ilen++) {
    const double n = (NumPoints-nseg-1) * (*ilen) / lentot + 1;

    if (n < *inpt) {
      *inpt = (int) round(n);
      *inptdec = n - *inpt;
    }

    npthave += *inpt;

    inpt++;
    inptdec++;
  }

  // Check if total selected points are more than needed due to rounding up
  while (npthave > NumPoints) {
    const int ichange = min_element(nptsegdec.begin(), nptsegdec.end())
                      - nptsegdec.begin();

    nptseg[ichange]--;
    nptsegdec[ichange] = 0;
    npthave--;
  }

  // Check if total selected points are fewer than needed due to rounding down
  while (npthave < NumPoints) {
    const int ichange = max_element(nptsegdec.begin(), nptsegdec.end())
                      - nptsegdec.begin();

    nptseg[ichange]++;
    nptsegdec[ichange] = 0;
    npthave++;
  }

  // Find uniformly spaced points over each segment between turn points
  ipt = Streamline.begin();
  idist = ptdist.begin();
  ilen = lenseg.begin();

  cpts.insert(cpts.end(), ipt, ipt+3);

  for (vector<int>::const_iterator inpt = nptseg.begin();
                                   inpt != nptseg.end(); inpt++) {
    const double inc = *ilen / *inpt;
    double targetlen = 0, runlen = 0, runlen0 = 0;

    for (int k = *inpt; k > 0; k--) {
      targetlen += inc;

      while (runlen < targetlen && idist < ptdist.end()) {
        runlen0 = runlen;
        runlen += *idist;
        ipt += 3;
        idist++;
      }

      if (fabs(targetlen - runlen0) < fabs(targetlen - runlen)) {
        runlen = runlen0;
        ipt -= 3;
        idist--;
      }

      cpts.insert(cpts.end(), ipt, ipt+3);
    }

    ilen++;
  }

  mControlPoints.push_back(cpts);

  // Find spread of streamlines around each control point
  for (vector<int>::const_iterator icpt = cpts.begin();
                                   icpt != cpts.end(); icpt +=3) {
    vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                                 ivalid2 = mIsInEnd2.begin();

    fill(sum.begin(), sum.end(), 0);
    fill(sumsq.begin(), sumsq.end(), 0);

    for (vector< vector<int> >::const_iterator istr = mStreamlines.begin();
                                           istr != mStreamlines.end(); istr++) {
      if (*ivalid1 && *ivalid2) {
        double dmin = 1e6;

        for (vector<int>::const_iterator ipt = istr->begin();
                                         ipt != istr->end(); ipt += 3) {
          double dist = 0;

          for (int k = 0; k < 3; k++) {
            diff[k] = ipt[k] - icpt[k];
            dist += pow(diff[k], 2);
          }

          if (dist < dmin && dist > 0) {
            dmin = dist;
            copy(diff.begin(), diff.end(), diffmin.begin());
          }
        }

        for (int k = 0; k < 3; k++) {
          sum[k] += diffmin[k];
          sumsq[k] += pow(diffmin[k], 2);
        }
      }

      ivalid1++;
      ivalid2++;
    }

    if (mNumMask > 2)
      for (int k = 0; k < 3; k++)
        sumsq[k] = sqrt((sumsq[k] - sum[k] * sum[k] / (mNumMask-1))
                        / (mNumMask-2));
    else
      copy(sum.begin(), sum.end(), sumsq.begin());

    cstd.insert(cstd.end(), sumsq.begin(), sumsq.end());
  }

  mControlStd.push_back(cstd);
}

//
// Write histograms, priors, end ROIs, and initialization control points
//
void Blood::WriteOutputs(const char *OutBase) {
  const int denom = mNumTrain + 2;
  char fname[PATH_MAX];
  vector<bool>::const_iterator ivalid1, ivalid2;
  vector< vector<int> >::const_iterator istr, ihisto;
  vector< vector<float> >::const_iterator iprior, idmean, idstd;
  vector< set<unsigned int> >::const_iterator iids;
  MRI *out1 = MRIclone(mTestMask, NULL);
  MRI *out2 = MRIclone(mTestMask, NULL);
  MRI *out3 = MRIclone(mTestMask, NULL);
  ofstream outfile;

  cout << "Writing output files to " << OutBase << "_*" << endl;

  // Write streamline-wise histogram to volume
  sprintf(fname, "%s_histo_str.nii.gz", OutBase);

  for (istr = mStreamlines.begin(); istr != mStreamlines.end(); istr++)
    for (vector<int>::const_iterator ipt = istr->begin();
                                     ipt != istr->end(); ipt += 3)
      MRIsetVoxVal(out1, ipt[0], ipt[1], ipt[2], 0, 
                   MRIgetVoxVal(out1, ipt[0], ipt[1], ipt[2], 0) + 1);

  MRIwrite(out1, fname);

  // Write total number of streamlines to text file
  sprintf(fname, "%s_histo_nstr.txt", OutBase);
  ofstream nstrfile(fname, ios::out);

  if (!nstrfile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  nstrfile << mStreamlines.size() << endl;

  // Write subject-wise histogram to volume
  sprintf(fname, "%s_histo.nii.gz", OutBase);
  MRIclear(out1);
  istr = mStreamlines.begin();

  for (vector<int>::const_iterator inum = mNumLines.begin();
                                   inum != mNumLines.end(); inum++) {
    MRIclear(out2);

    for (int k = *inum; k > 0; k--) {
      for (vector<int>::const_iterator ipt = istr->begin();
                                       ipt != istr->end(); ipt += 3)
        MRIsetVoxVal(out2, ipt[0], ipt[1], ipt[2], 0, 1);
      istr++;
    }

    MRIadd(out1, out2, out1);
  }

  MRIwrite(out1, fname);

  // Convert histogram to negative log-likelihood
  MRIclear(out2);

  for (int iz = 0; iz < mNz; iz++)
    for (int iy = 0; iy < mNy; iy++)
      for (int ix = 0; ix < mNx; ix++) {
        const float ratio = (MRIgetVoxVal(out1, ix, iy, iz, 0) + 1) / denom;
        MRIsetVoxVal(out2, ix, iy, iz, 0, -log(ratio));
        MRIsetVoxVal(out3, ix, iy, iz, 0, -log(1-ratio));
      }

  sprintf(fname, "%s_logprior_1.nii.gz", OutBase);
  MRIwrite(out2, fname);

  sprintf(fname, "%s_logprior_0.nii.gz", OutBase);
  MRIwrite(out3, fname);

  // Save anatomical label IDs found in training set, histograms and priors
  iids = mIdsLocal.begin();
  ihisto = mHistoLocal.begin();
  iprior = mPriorLocal.begin();

  for (vector<int>::const_iterator idir = mDirLocal.begin();
                                   idir < mDirLocal.end(); idir += 3) {
    const int idx = idir[0], idy = idir[1], idz = idir[2];

    sprintf(fname, "%s_fsids_%d_%d_%d.txt", OutBase, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (vector< set<unsigned int> >::const_iterator iarc = iids; 
         iarc < mIdsLocal.end(); iarc += mNumLocal) {
      for (set<unsigned int>::const_iterator ival = iarc->begin();
                                             ival != iarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;
    }

    outfile.close();

    sprintf(fname, "%s_fshisto_%d_%d_%d.txt", OutBase, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (vector< vector<int> >::const_iterator iarc = ihisto;
         iarc < mHistoLocal.end(); iarc += mNumLocal) {
      for (vector<int>::const_iterator ival = iarc->begin();
                                       ival != iarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;
    }

    outfile.close();

    sprintf(fname, "%s_fsprior_%d_%d_%d.txt", OutBase, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (vector< vector<float> >::const_iterator iarc = iprior;
         iarc < mPriorLocal.end(); iarc += mNumLocal) {
      for (vector<float>::const_iterator ival = iarc->begin();
                                         ival != iarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;
    }

    outfile.close();

    iids++;
    ihisto++;
    iprior++;
  }

  iids = mIdsNear.begin();
  ihisto = mHistoNear.begin();
  iprior = mPriorNear.begin();
  idmean = mAsegDistMean.begin();
  idstd = mAsegDistStd.begin();

  for (vector<int>::const_iterator idir = mDirNear.begin();
                                   idir < mDirNear.end(); idir += 3) {
    const int idx = idir[0], idy = idir[1], idz = idir[2];

    sprintf(fname, "%s_fsnnids_%d_%d_%d.txt", OutBase, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (vector< set<unsigned int> >::const_iterator iarc = iids; 
         iarc < mIdsNear.end(); iarc += mNumNear) {
      for (set<unsigned int>::const_iterator ival = iarc->begin();
                                             ival != iarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;
    }

    outfile.close();

    sprintf(fname, "%s_fsnnhisto_%d_%d_%d.txt", OutBase, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (vector< vector<int> >::const_iterator iarc = ihisto;
         iarc < mHistoNear.end(); iarc += mNumNear) {
      for (vector<int>::const_iterator ival = iarc->begin();
                                       ival != iarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;
    }

    outfile.close();

    sprintf(fname, "%s_fsnnprior_%d_%d_%d.txt", OutBase, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (vector< vector<float> >::const_iterator iarc = iprior;
         iarc < mPriorNear.end(); iarc += mNumNear) {
      for (vector<float>::const_iterator ival = iarc->begin();
                                         ival != iarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;
    }

    outfile.close();

    sprintf(fname, "%s_fsnndmean_%d_%d_%d.txt", OutBase, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (vector< vector<float> >::const_iterator iarc = idmean;
         iarc < mAsegDistMean.end(); iarc += mNumNear) {
      for (vector<float>::const_iterator ival = iarc->begin();
                                         ival != iarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;
    }

    outfile.close();

    sprintf(fname, "%s_fsnndstd_%d_%d_%d.txt", OutBase, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (vector< vector<float> >::const_iterator iarc = idstd;
         iarc < mAsegDistStd.end(); iarc += mNumNear) {
      for (vector<float>::const_iterator ival = iarc->begin();
                                         ival != iarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;
    }

    outfile.close();

    iids++;
    ihisto++;
    iprior++;
    idmean++;
    idstd++;
  }

  // Save tangent vector distribution in training set
  sprintf(fname, "%s_tangmean.txt", OutBase);
  outfile.open(fname, ios::out);

  for (vector<float>::const_iterator imean = mTangentMean.begin();
                                     imean != mTangentMean.end(); imean += 3)
    outfile << imean[0] << " " << imean[1] << " " << imean[2] << endl;

  outfile.close();

  sprintf(fname, "%s_tangstd.txt", OutBase);
  outfile.open(fname, ios::out);

  for (vector<float>::const_iterator istd = mTangentStd.begin();
                                     istd != mTangentStd.end(); istd += 3)
    outfile << istd[0] << " " << istd[1] << " " << istd[2] << endl;

  outfile.close();

  // Save curvature distribution in training set
  sprintf(fname, "%s_curvmean.txt", OutBase);
  outfile.open(fname, ios::out);

  for (vector<float>::const_iterator imean = mCurvatureMean.begin();
                                     imean != mCurvatureMean.end(); imean++)
    outfile << *imean << endl;

  outfile.close();

  sprintf(fname, "%s_curvstd.txt", OutBase);
  outfile.open(fname, ios::out);

  for (vector<float>::const_iterator istd = mCurvatureStd.begin();
                                     istd != mCurvatureStd.end(); istd++)
    outfile << *istd << endl;

  outfile.close();

  // Write end ROIs to volumes
  MRIclear(out1);
  MRIclear(out2);
  ivalid1 = mIsInEnd1.begin();
  ivalid2 = mIsInEnd2.begin();

  for (istr = mStreamlines.begin(); istr != mStreamlines.end(); istr++) {
    if (*ivalid1) {
      vector<int>::const_iterator iend1 = istr->begin();
    
      MRIsetVoxVal(out1, iend1[0], iend1[1], iend1[2], 0, 
                   MRIgetVoxVal(out1, iend1[0], iend1[1], iend1[2], 0) + 1);
    }

    if (*ivalid2) {
      vector<int>::const_iterator iend2 = istr->end() - 3;
    
      MRIsetVoxVal(out2, iend2[0], iend2[1], iend2[2], 0, 
                   MRIgetVoxVal(out2, iend2[0], iend2[1], iend2[2], 0) + 1);
    }

    ivalid1++;
    ivalid2++;
  }

  sprintf(fname, "%s_end1.nii.gz", OutBase);
  MRIwrite(out1, fname);

  sprintf(fname, "%s_end2.nii.gz", OutBase);
  MRIwrite(out2, fname);

  // Write central streamline to text file and volume
  sprintf(fname, "%s_cpts_all.txt", OutBase);
  ofstream centfile(fname, ios::out);

  if (!centfile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  MRIclear(out1);

  for (vector<int>::const_iterator ipt = mCenterStreamline.begin();
                                   ipt < mCenterStreamline.end(); ipt += 3) {
    centfile << ipt[0] << " " << ipt[1] << " " << ipt[2] << endl;
    MRIsetVoxVal(out1, ipt[0], ipt[1], ipt[2], 0, 1);
  }

  sprintf(fname, "%s_cpts_all.nii.gz", OutBase);
  MRIwrite(out1, fname);

  // Write each set of control points chosen from central streamline
  // to text file and volume
  for (istr = mControlPoints.begin(); istr != mControlPoints.end(); istr++) {
    const unsigned int ncpt = istr->size() / 3;

    sprintf(fname, "%s_cpts_%d.txt", OutBase, ncpt);
    ofstream cptfile(fname, ios::out);

    if (!cptfile) {
      cout << "ERROR: Could not open " << fname << " for writing" << endl;
      exit(1);
    }

    for (vector<int>::const_iterator icpt = istr->begin();
                                     icpt < istr->end(); icpt += 3)
      cptfile << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;

    Spline spline(*istr, mTestMask); 
    spline.InterpolateSpline();

    sprintf(fname, "%s_cpts_%d.nii.gz", OutBase, ncpt);
    spline.WriteVolume(fname, true);
  }

  // Write spread of streamlines around each control point to text file
  for (vector< vector<float> >::const_iterator istd = mControlStd.begin();
                                            istd != mControlStd.end(); istd++) {
    const unsigned int ncpt = istd->size() / 3;

    sprintf(fname, "%s_cpts_%d_std.txt", OutBase, ncpt);
    ofstream stdfile(fname, ios::out);

    if (!stdfile) {
      cout << "ERROR: Could not open " << fname << " for writing" << endl;
      exit(1);
    }

    for (vector<float>::const_iterator ival = istd->begin();
                                       ival < istd->end(); ival += 3)
      stdfile << ival[0] << " " << ival[1] << " " << ival[2] << endl;
  }

  MRIfree(&out1);
  MRIfree(&out2);
  MRIfree(&out3);
}

void Blood::PrintStreamline(int SubjIndex, int LineIndex) {
  vector<int>::const_iterator inum = mNumLines.begin();
  vector< vector<int> >::const_iterator istr = mStreamlines.begin();

  if ((SubjIndex < 0) || (SubjIndex >= mNumTrain))
    cout << "WARN: Subject index " << SubjIndex << " is out of range" << endl;

  for (int k = SubjIndex; k > 0; k--) {
    istr += (*inum);
    inum++;
  }

  if ((LineIndex < 0) || (LineIndex >= (*inum)))
    cout << "WARN: Line index " << LineIndex << " is out of range" << endl;

  istr += LineIndex;

  for (vector<int>::const_iterator ipt = istr->begin();
                                   ipt != istr->end(); ipt += 3)
    cout << ipt[0] << " " << ipt[1] << " " << ipt[2] << endl;
}

