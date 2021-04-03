/**
 * @brief Training data and methods that probabilistic tractography feeds on
 *
 * Training data and methods that probabilistic tractography feeds on
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

#include <blood.h>

using namespace std;

const int Blood::mDistThresh = 4,
          Blood::mEndDilation = 2;
const unsigned int Blood::mDiffStep = 3;
const float Blood::mLengthCutoff = 0.05,
            Blood::mLengthRatio = 3.0,
            Blood::mHausStepRatio = 0.05,
            Blood::mControlStepRatio = 0.8,
            Blood::mTangentBinSize = 1/3.0,	// 0.1,
            Blood::mCurvatureBinSize = 0.01;	// 0.002;

Blood::Blood(const std::string TrainListFile, const std::string TrainTrkFile,
             const std::string TrainRoi1File, const std::string TrainRoi2File,
             const std::string TrainAsegFile, const std::string TrainMaskFile,
             float TrainMaskLabel, const std::string ExcludeFile,
             const vector<std::string> &TestMaskList,
             const vector<std::string> &TestFaList,
             const std::string TestAffineXfmFile,
             const std::string TestNonlinXfmFile,
             const std::string TestNonlinRefFile,
             const vector<std::string> &TestBaseXfmList,
             const std::string TestBaseMaskFile,
             bool UseTruncated, vector<int> &NumControls,
             bool Debug) :
             mDebug(Debug),
             mUseTruncated(UseTruncated),
             mMaskLabel(TrainMaskLabel) {
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
  MRI *testvol;

  // Directions in which to look for neighboring labels
  mNumLocal = 15;
  mDirLocal.insert(mDirLocal.begin(), dirs, dirs+45);

  mNumNear = 14;
  mDirNear.insert(mDirNear.begin(), dirs+3, dirs+45);
  
  // Read brain mask
  for (auto ifile = TestMaskList.begin();
       ifile < TestMaskList.end(); ifile++) {
    cout << "Loading brain mask of output subject from " << *ifile << endl;
    testvol = MRIread((*ifile).c_str());
    if (!testvol) {
      cout << "ERROR: Could not read " << *ifile << endl;
      exit(1);
    }
    mTestMask.push_back(testvol);
  }

  // Size of volumes in common space
  mNx = mTestMask[0]->width;
  mNy = mTestMask[0]->height;
  mNz = mTestMask[0]->depth;

  // Resolution of volumes in common space
  mDx = mTestMask[0]->xsize;

  // Read FA map
  for (auto ifile = TestFaList.begin();
       ifile < TestFaList.end(); ifile++) {
    cout << "Loading FA map of output subject from " << *ifile << endl;
    testvol = MRIread((*ifile).c_str());
    if (!testvol) {
      cout << "ERROR: Could not read " << *ifile << endl;
      exit(1);
    }
    mTestFa.push_back(testvol);
  }

  // Read base reference volume, if any
  if (!TestBaseMaskFile.empty()) {
    cout << "Loading base mask of output subject from " << TestBaseMaskFile
         << endl;
    mTestBaseMask = MRIread(TestBaseMaskFile.c_str());
    if (!mTestBaseMask) {
      cout << "ERROR: Could not read " << TestBaseMaskFile << endl;
      exit(1);
    }
  } else if (!mTestFa.empty()) {
    mTestBaseMask = mTestFa[0];
  } else {
    mTestBaseMask = mTestMask[0];
  }

  // Read atlas-to-base registration files
#ifndef NO_CVS_UP_IN_HERE
  if (!TestNonlinXfmFile.empty()) {
    mTestNonlinReg.ReadXfm(TestNonlinXfmFile.c_str(), mTestMask[0]);

    if (!TestAffineXfmFile.empty()) {
      MRI *refvol;
      
      cout << "Loading non-linear registration source for output subject from "
           << TestNonlinRefFile << endl;
      refvol = MRIread(TestNonlinRefFile.c_str());
      if (!refvol) {
        cout << "ERROR: Could not read " << TestNonlinRefFile << endl;
        exit(1);
      }

      mTestAffineReg.ReadXfm(TestAffineXfmFile.c_str(), refvol, mTestBaseMask);
    }
  } else {
#endif
    if (!TestAffineXfmFile.empty()) {
      mTestAffineReg.ReadXfm(TestAffineXfmFile.c_str(), mTestMask[0], mTestBaseMask);
    }
#ifndef NO_CVS_UP_IN_HERE
  }
#endif

  // Read base-to-DWI registration files
  for (auto ifile = TestBaseXfmList.begin();
       ifile < TestBaseXfmList.end(); ifile++) {
    const unsigned int iframe = ifile - TestBaseXfmList.begin();
    AffineReg basereg;

    basereg.ReadXfm((*ifile).c_str(), mTestBaseMask, mTestFa[iframe]);
    mTestBaseReg.push_back(basereg);
  }

  // Set number(s) of control points to fit for first pathway
  SetNumControls(NumControls);

  // Read training data for first pathway
  ReadStreamlines(TrainListFile, TrainTrkFile, TrainRoi1File, TrainRoi2File,
                  TrainMaskLabel, ExcludeFile);

  // Read training subjects' anatomy
  ReadAnatomy(TrainListFile, TrainAsegFile, TrainMaskFile);

  // Allocate space for histograms
  mHistoStr  = MRIclone(mTestMask[0], NULL);
  mHistoSubj = MRIclone(mTestMask[0], NULL);
}

Blood::Blood(const std::string TrainTrkFile,
             const std::string TrainRoi1File, const std::string TrainRoi2File,
             bool Debug) : 
             mDebug(Debug),
             mUseTruncated(false),
             mNx(0), mNy(0), mNz(0) {
  // Read single input streamline file
  ReadStreamlines(0, TrainTrkFile, TrainRoi1File, TrainRoi2File, 0, 0);

  if (!TrainRoi1File.c_str()) {
    // Allocate space for histograms
    mHistoStr  = MRIclone(mRoi1[0], NULL);
    mHistoSubj = MRIclone(mRoi1[0], NULL);
  }
  else {
    mHistoStr  = 0;
    mHistoSubj = 0;
  }
}

Blood::~Blood() {
  vector<MRI *>::iterator ivol;

  for (ivol = mRoi1.begin(); ivol != mRoi1.end(); ivol++)
    MRIfree(&(*ivol));

  for (ivol = mRoi2.begin(); ivol != mRoi2.end(); ivol++)
    MRIfree(&(*ivol));

  for (ivol = mAseg.begin(); ivol != mAseg.end(); ivol++)
    MRIfree(&(*ivol));

  for (ivol = mMask.begin(); ivol != mMask.end(); ivol++)
    MRIfree(&(*ivol));

  for (ivol = mTestMask.begin(); ivol != mTestMask.end(); ivol++)
    MRIfree(&(*ivol));

  for (ivol = mTestFa.begin(); ivol != mTestFa.end(); ivol++)
    MRIfree(&(*ivol));

  if (mHistoStr)
    MRIfree(&mHistoStr);

  if (mHistoSubj)
    MRIfree(&mHistoSubj);
}

//
// Set the number(s) of control points to fit to the center streamline
//
void Blood::SetNumControls(vector<int> &NumControls) {
  mNumControls.resize(NumControls.size());
  copy(NumControls.begin(), NumControls.end(), mNumControls.begin());

  mControlPoints.clear();
  mControlPointsTest.clear();
  mControlStd.clear();
  mControlStdTest.clear();
}

//
// Read streamlines of training subjects
// (and the ROIs that were used to label them)
//
void Blood::ReadStreamlines(const std::string TrainListFile,
                            const std::string TrainTrkFile,
                            const std::string TrainRoi1File,
                            const std::string TrainRoi2File,
                            float TrainMaskLabel,
                            const std::string ExcludeFile) {
  int nrejmask = 0, nrejrev = 0;
  vector<string> dirlist;
  vector<MRI *>::iterator ivol;

  if (!TrainListFile.empty()) {		// Read multiple inputs from a list
    string dirname;
    ifstream listfile(TrainListFile, ios::in);

    if (!listfile) {
      cout << "ERROR: Could not open " << TrainListFile << endl;
      exit(1);
    }

    while (listfile >> dirname)
      dirlist.push_back(dirname + "/");
  }
  else				// Read a single input
    dirlist.push_back("");

  mMaskLabel = TrainMaskLabel;

  // Clear all variables that depend on pathway
  mStreamlines.clear();
  mLengths.clear();
  mTruncatedLengths.clear();
  mNumLines.clear();
  mIsInEnd1.clear();
  mIsInEnd2.clear();
  mMidPoints.clear();
  mMeanEnd1.clear();
  mMeanEnd2.clear();
  mMeanMid.clear();
  mVarEnd1.clear();
  mVarEnd2.clear();
  mVarMid.clear();
  mExcludedStreamlines.clear();
  mCenterStreamline.clear();
  mControlPoints.clear();
  mControlPointsTest.clear();
  mControlStd.clear();
  mControlStdTest.clear();
  mHistoLocal.clear();
  mHistoLocalAll.clear();
  mPriorLocal.clear();
  mPriorLocalAll.clear();
  mIdsLocal.clear();
  mIdsLocalAll.clear();
  mHistoNear.clear();
  mHistoNearAll.clear();
  mPriorNear.clear();
  mPriorNearAll.clear();
  mIdsNear.clear();
  mIdsNearAll.clear();
  mAsegDistMean.clear();
  mAsegDistMeanAll.clear();
  mAsegDistStd.clear();
  mAsegDistStdAll.clear();
  mHistoTangent.clear();
  mHistoTangentAll.clear();
  mPriorTangent.clear();
  mPriorTangentAll.clear();
  mTangentXByArc.clear();
  mTangentXByArcAll.clear();
  mTangentYByArc.clear();
  mTangentYByArcAll.clear();
  mHistoCurvature.clear();
  mHistoCurvatureAll.clear();
  mPriorCurvature.clear();
  mPriorCurvatureAll.clear();
  mCurvatureByArc.clear();
  mCurvatureByArcAll.clear();

  for (ivol = mRoi1.begin(); ivol != mRoi1.end(); ivol++)
    MRIfree(&(*ivol));
  mRoi1.clear();

  for (ivol = mRoi2.begin(); ivol != mRoi2.end(); ivol++)
    MRIfree(&(*ivol));
  mRoi2.clear();

  mVolume = 0;
  mNumStr = 0;
  mLengthMin = 0;
  mLengthMax = 0;
  mLengthAvg = 0;
  mNumStrEnds = 0;
  mLengthMinEnds = 0;
  mLengthMaxEnds = 0;
  mLengthAvgEnds = 0;

  for (vector<string>::const_iterator idir = dirlist.begin();
                                      idir != dirlist.end(); idir++) {
    int npts, nlines = 0;
    CTrackReader trkreader;
    TRACK_HEADER trkheader;
    string fname;

    if (!TrainRoi1File.empty()) {
      fname = *idir + TrainRoi1File;
      cout << "Loading streamline start ROI from " << fname << endl;
      mRoi1.push_back(MRIread(fname.c_str()));
      if (! *(mRoi1.end()-1)) {
        cout << "ERROR: Could not read " << fname << endl;
        exit(1);
      }
    }

    if (!TrainRoi2File.empty()) {
      fname = *idir + TrainRoi2File;
      cout << "Loading streamline end ROI from " << fname << endl;
      mRoi2.push_back(MRIread(fname.c_str()));
      if (! *(mRoi2.end()-1)) {
        cout << "ERROR: Could not read " << fname << endl;
        exit(1);
      }
    }

    if (mTestMask.empty()) {	// Get size of reference volume from ROIs
      if (!mRoi1.empty()) {
        mNx = mRoi1[0]->width;
        mNy = mRoi1[0]->height;
        mNz = mRoi1[0]->depth;
      }
      else if (!mRoi2.empty()) {
        mNx = mRoi2[0]->width;
        mNy = mRoi2[0]->height;
        mNz = mRoi2[0]->depth;
      }
    }

    fname = *idir + TrainTrkFile;
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
      int nptsmask = npts,
          xroi1 = 0, xroi2 = 0, xroi12 = 0, xroi21 = 0;
      float forwback1 = 0, forw1 = 0,
            forwback2 = 0, forw2 = 0;
      vector<int> pts;
      vector<float> dir1(3), dir2(3);
      float *rawpts = new float[npts*3], *rawptsmask = rawpts, *iraw;

      // Read a streamline from input file
      trkreader.GetNextTrackData(npts, rawpts);

      // Divide by voxel size, make 0-based, and round to get voxel coords
      iraw = rawpts; 
      for (int ipt = npts; ipt > 0; ipt--)
        for (int k = 0; k < 3; k++) {
          *iraw = round(*iraw / trkheader.voxel_size[k] - .5);
          iraw++;
        }

      // Disregard first points if they are off the mask
      while (!IsInMask(rawptsmask) && nptsmask > 0) {
        rawptsmask += 3;
        nptsmask--;
      }

      // Disregard last points if they are off the mask
      iraw = rawpts + (npts-1)*3; 
      while (!IsInMask(iraw) && nptsmask > 0) {
        iraw -= 3;
        nptsmask--;
      }

      if (nptsmask == 0) {			// No points in mask
        delete[] rawpts;
        continue;
      }

      iraw = rawptsmask; 
      for (int ipt = nptsmask; ipt > 1; ipt--)	// Remove duplicate points
        if ( (iraw[0] != iraw[3]) || (iraw[1] != iraw[4]) ||
                                     (iraw[2] != iraw[5]) ) {
          for (int k = 0; k < 3; k++) {
            pts.push_back((int) *iraw);
            iraw++;
          }

          // Check that streamline stays within mask
          if (!IsInMask(pts.end()-3)) {
            isinmask = false;
            nrejmask++;
            break;
          }

          // Check directions in which streamline traverses start ROI
          if (!mRoi1.empty()) {
            if (IsInRoi(pts.end()-3, *(mRoi1.end()-1))) {
              if (xroi1 == 0) {		// Entering ROI for the first time
                xroi1 = 1;

                if (xroi2 > 0)		// Have already been to the other ROI
                  xroi21 = 1;

                iraw -= 3;
                forw1 = 0;
                if (iraw == rawptsmask)
                  iraw += 3;
                else
                  for (int k = 0; k < 3; k++) {		// Direction of entering
                    const float diff = *iraw - *(iraw-3);

                    dir1[k] = diff;
                    forw1 += diff*diff;
                    iraw++;
                  } 
              }

              if (xroi1 == 2) {		// Re-entering ROI
                if (xroi12 > 0) {	// Went back and forth between ROIs
                  xroi1 = 3;
                  nrejrev++;
                  break;
                }

                xroi1 = 1;

                iraw -= 3;
                forwback1 = 0;
                for (int k = 0; k < 3; k++) {		// Direction of entering
                  const float diff = *iraw - *(iraw-3);

                  forwback1 += dir1[k]*diff;
                  iraw++;
                }

                // Check if entering back the way I exited
                if (forwback1 < 0) {		// |angle| > 90
                  xroi1 = 3;
                  nrejrev++;
                  break;
                }
              }
            }
            else
              if (xroi1 == 1) {		// Exiting ROI
                xroi1 = 2;

                if (forw1 == 0) {	// Never entered, started inside ROI
                  iraw -= 3;
                  forw1 = 0;
                  for (int k = 0; k < 3; k++) { 	// Direction of exiting
                    const float diff = *iraw - *(iraw-3);

                    dir1[k] = diff;
                    forw1 += diff*diff;
                    iraw++;
                  }
                }
                else {
                  iraw -= 3;
                  forwback1 = 0;
                  for (int k = 0; k < 3; k++) {		// Direction of exiting
                    const float diff = *iraw - *(iraw-3);

                    forwback1 += dir1[k]*diff;
                    iraw++;
                  } 

                  // Check if exiting back the way I entered
                  if (forwback1 < 0) {		// |angle| > 90
                    xroi1 = 3;
                    nrejrev++;
                    break;
                  }
                }
              }
          }

          // Check directions in which streamline traverses end ROI
          if (!mRoi2.empty()) {
            if (IsInRoi(pts.end()-3, *(mRoi2.end()-1))) {
              if (xroi2 == 0) {		// Entering ROI for the first time
                xroi2 = 1;

                if (xroi1 > 0)		// Have already been to the other ROI
                  xroi12 = 1;

                iraw -= 3;
                forw2 = 0;
                if (iraw == rawptsmask)
                  iraw += 3;
                else
                  for (int k = 0; k < 3; k++) {		// Direction of entering
                    const float diff = *iraw - *(iraw-3);

                    dir2[k] = diff;
                    forw2 += diff*diff;
                    iraw++;
                  } 
              }

              if (xroi2 == 2) {		// Re-entering ROI
                if (xroi21 > 0) {	// Went back and forth between ROIs
                  xroi2 = 3;
                  nrejrev++;
                  break;
                }

                xroi2 = 1;

                iraw -= 3;
                forwback2 = 0;
                for (int k = 0; k < 3; k++) {         // Direction of entering
                  const float diff = *iraw - *(iraw-3);

                  forwback2 += dir2[k]*diff;
                  iraw++;
                }

                // Check if entering back the way I exited
                if (forwback2 < 0) {		// |angle| > 90
                  xroi2 = 3;
                  nrejrev++;
                  break;
                }
              }
            }
            else
              if (xroi2 == 1) {		// Exiting ROI
                xroi2 = 2;

                if (forw2 == 0) {	// Never entered, started inside ROI
                  iraw -= 3;
                  forw2 = 0;
                  for (int k = 0; k < 3; k++) { 	// Direction of exiting
                    const float diff = *iraw - *(iraw-3);

                    dir2[k] = diff;
                    forw2 += diff*diff;
                    iraw++;
                  }
                }
                else {
                  iraw -= 3;
                  forwback2 = 0;
                  for (int k = 0; k < 3; k++) {		// Direction of exiting
                    const float diff = *iraw - *(iraw-3);

                    forwback2 += dir2[k]*diff;
                    iraw++;
                  } 

                  // Check if exiting back the way I entered
                  if (forwback2 < 0) {		// |angle| > 90
                    xroi2 = 3;
                    nrejrev++;
                    break;
                  }
                }
              }
          }
        }
        else
          iraw += 3;

      // Don't save streamline if it strayed off the test subject's mask
      // or if it traversed the same ROI in opposite directions
      if (!isinmask || xroi1 == 3 || xroi2 == 3) {
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

  if (!mAseg.empty() && ((int) mAseg.size() != mNumTrain)) {
    cout << "ERROR: Number of training paths and segmentations does not match"
         << endl;
    exit(1);
  }

  cout << "INFO: Rejected " << nrejmask
       << " streamlines for straying off mask" << endl;
  cout << "INFO: Rejected " << nrejrev
       << " streamlines for reversing direction" << endl;

  RemoveLengthOutliers();

  ComputeStats();

  if (!ExcludeFile.empty()) {
    ReadExcludedStreamlines(ExcludeFile);
  }
}

//
// Read list of streamlines to be excluded from search for center streamline
//
void Blood::ReadExcludedStreamlines(const std::string ExcludeFile) {
  string excline;
  ifstream excfile;

  mExcludedStreamlines.clear();

  cout << "Loading list of excluded streamlines from " << ExcludeFile << endl;
  excfile.open(ExcludeFile, ios::in);

  if (!excfile) {
    cout << "WARN: Could not open " << ExcludeFile << endl;
  }
  else {
    getline(excfile, excline);

    while (excline.compare("exclude") == 0) {
      vector<int> points;

      while (getline(excfile, excline) && excline.compare("exclude") != 0) {
        float coord;
        istringstream excstr(excline);

        while (excstr >> coord)
          points.push_back((int) round(coord));
      }

      if (points.size() % 3 != 0) {
        cout << "ERROR: File " << ExcludeFile
             << " must contain triplets of coordinates" << endl;
        exit(1);
      }

      mExcludedStreamlines.push_back(points);
    }

    excfile.close();
  }

  cout << "INFO: Found " << mExcludedStreamlines.size() 
       << " streamlines to be excluded" << endl;
}

//
// Get info on input streamlines:
// - Total number
// - Average, min, and max length
//
void Blood::ComputeStats() {
  // Number of streamlines
  mNumStr = mStreamlines.size();

  // Streamline length
  if (mNumStr == 0) {
    mLengthMin = 0;
    mLengthMax = 0;
    mLengthAvg = 0;
  }
  else {
    int lsum = 0;

    for (vector<int>::const_iterator ilen = mLengths.begin();
                                     ilen < mLengths.end(); ilen++)
      lsum += *ilen;

    mLengthAvg = lsum / (float) mNumStr;

    mLengthMin = *min_element(mLengths.begin(), mLengths.end());
    mLengthMax = *max_element(mLengths.begin(), mLengths.end());
  }

  cout << "INFO: Have " << mNumStr
       << " total streamlines (min/mean/max length: "
       << mLengthMin << "/" << round(mLengthAvg) << "/"
       << mLengthMax << ")" << endl;

  if (mMask.empty() && mRoi1.empty() && mRoi2.empty()) {
    mNumStrEnds = mNumStr;
    mLengthMinEnds = mLengthMin;
    mLengthMaxEnds = mLengthMax;
    mLengthAvgEnds = mLengthAvg;
  }
}

//
// Get info on input streamlines that don't have truncated end points:
// - Total number
// - Average, min, and max length
//
void Blood::ComputeStatsEnds() {
  if (!mIsInEnd1.empty() && !mIsInEnd2.empty()) {
    int lsum = 0;
    vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                                 ivalid2 = mIsInEnd2.begin();

    mNumStrEnds = 0;
    mLengthMinEnds = mLengthMax;
    mLengthMaxEnds = 0;

    for (vector<int>::const_iterator ilen = mLengths.begin();
                                     ilen < mLengths.end(); ilen++) {
      if (*ivalid1 && *ivalid2) {	// Neither end is truncated
        // Number of streamlines
        mNumStrEnds++;

        // Streamline length
        lsum += *ilen;

        if (*ilen < mLengthMinEnds)
          mLengthMinEnds = *ilen;

        if (*ilen > mLengthMaxEnds)
          mLengthMaxEnds = *ilen;
      }

      ivalid1++;
      ivalid2++;
    }

    if (mNumStrEnds == 0) {
      mLengthMinEnds = 0;
      mLengthAvgEnds = 0;
    }
    else
      mLengthAvgEnds = lsum / (float) mNumStrEnds;

    cout << "INFO: Have " << mNumStrEnds
         << " non-truncated streamlines (min/mean/max length: "
         << mLengthMinEnds << "/" << round(mLengthAvgEnds) << "/"
         << mLengthMaxEnds << ")" << endl;
  }
}

//
// Compute center of mass of end points to constrain center streamline selection
//
void Blood::ComputeEndPointCoM() {
  vector<int>::iterator imidpts;
  vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                               ivalid2 = mIsInEnd2.begin();
  vector<MRI *>::const_iterator iroi1 = mRoi1.begin(),
                                iroi2 = mRoi2.begin();
  vector< vector<int> >::const_iterator istr = mStreamlines.begin();
  vector<int> sum1(3, 0), sum2(3, 0), summ(3, 0),
              sumsq1(3, 0), sumsq2(3, 0), sumsqm(3, 0);

  mMidPoints.resize(mNumStrEnds);
  imidpts = mMidPoints.begin();

  mMeanEnd1.resize(3);
  mMeanEnd2.resize(3);
  mMeanMid.resize(3);
  mVarEnd1.resize(3);
  mVarEnd2.resize(3);
  mVarMid.resize(3);

  for (vector<int>::const_iterator inum = mNumLines.begin();
                                   inum != mNumLines.end(); inum++) {
    for (int k = *inum; k > 0; k--) {
      if (*ivalid1 && *ivalid2) {
        vector<int>::const_iterator itop    = istr->begin(),
                                    ibottom = istr->end() - 3,
                                    ipt1, ipt2, imiddle;

        // Distance from the streamline start to the start ROI
        for (ipt1 = istr->begin(); ipt1 < istr->end(); ipt1 += 3)
          if (MRIgetVoxVal(*iroi1, ipt1[0], ipt1[1], ipt1[2], 0) > 0)
            break;

        // Distance from the streamline start to the end ROI
        for (ipt2 = ipt1; ipt2 < istr->end(); ipt2 += 3)
          if (MRIgetVoxVal(*iroi2, ipt2[0], ipt2[1], ipt2[2], 0) > 0)
            break;

        // Midpoint of streamline between start and end ROI
        *imidpts = (int) (ipt1 - itop) +
                   (int) round(float(ipt2 - ipt1) / 6) * 3;
        imiddle = itop + *imidpts;

        for (int k = 0; k < 3; k++) {
          sum1[k]   += itop[k];
          sumsq1[k] += itop[k]*itop[k];
          sum2[k]   += ibottom[k];
          sumsq2[k] += ibottom[k]*ibottom[k];
          summ[k]   += imiddle[k];
          sumsqm[k] += imiddle[k]*imiddle[k];
        }

        imidpts++;
      }

      istr++;
      ivalid1++;
      ivalid2++;
    }

    iroi1++;
    iroi2++;
  }

  if (mNumStrEnds == 1) {
    for (int k = 0; k < 3; k++) {
      mMeanEnd1[k] = (float) sum1[k];
      mMeanEnd2[k] = (float) sum2[k];
      mMeanMid[k]  = (float) summ[k];
    }

    fill(mVarEnd1.begin(), mVarEnd1.end(), 0.0);
    fill(mVarEnd2.begin(), mVarEnd2.end(), 0.0);
    fill(mVarMid.begin(), mVarMid.end(), 0.0);
  }
  else
    for (int k = 0; k < 3; k++) {
      mMeanEnd1[k] = sum1[k] / float(mNumStrEnds);
      mVarEnd1[k] = (sumsq1[k] - mNumStrEnds * mMeanEnd1[k] * mMeanEnd1[k])
                    / (mNumStrEnds-1);
      mMeanEnd2[k] = sum2[k] / float(mNumStrEnds);
      mVarEnd2[k] = (sumsq2[k] - mNumStrEnds * mMeanEnd2[k] * mMeanEnd2[k])
                    / (mNumStrEnds-1);
      mMeanMid[k] = summ[k] / float(mNumStrEnds);
      mVarMid[k] = (sumsqm[k] - mNumStrEnds * mMeanMid[k] * mMeanMid[k])
                   / (mNumStrEnds-1);
    }

  cout << "INFO: Center of mass of start points: ("
       << round(mMeanEnd1[0]) << "+/-" << round(sqrt(mVarEnd1[0])) << ", "
       << round(mMeanEnd1[1]) << "+/-" << round(sqrt(mVarEnd1[1])) << ", "
       << round(mMeanEnd1[2]) << "+/-" << round(sqrt(mVarEnd1[2])) << ")"
       << endl;
  cout << "INFO: Center of mass of midpoints: ("
       << round(mMeanMid[0]) << "+/-" << round(sqrt(mVarMid[0])) << ", "
       << round(mMeanMid[1]) << "+/-" << round(sqrt(mVarMid[1])) << ", "
       << round(mMeanMid[2]) << "+/-" << round(sqrt(mVarMid[2])) << ")"
       << endl;
  cout << "INFO: Center of mass of end points: ("
       << round(mMeanEnd2[0]) << "+/-" << round(sqrt(mVarEnd2[0])) << ", "
       << round(mMeanEnd2[1]) << "+/-" << round(sqrt(mVarEnd2[1])) << ", "
       << round(mMeanEnd2[2]) << "+/-" << round(sqrt(mVarEnd2[2])) << ")" 
       << endl;
}

//
// Check that a point is inside the test subject's brain mask
//
bool Blood::IsInMask(vector<int>::const_iterator Point) {
  if (mTestMask.empty() && mRoi1.empty() && mRoi2.empty())
    return true;
  else {
    bool isinmask = (Point[0] > -1) && (Point[0] < mNx) &&
                    (Point[1] > -1) && (Point[1] < mNy) &&
                    (Point[2] > -1) && (Point[2] < mNz);

    for (vector<MRI *>::const_iterator imask = mTestMask.begin();
                                       imask < mTestMask.end(); imask++)
      isinmask = isinmask &&
                 (MRIgetVoxVal(*imask, Point[0], Point[1], Point[2], 0) > 0);

    return isinmask;
  }
}

bool Blood::IsInMask(float *Point) {
  if (mTestMask.empty() && mRoi1.empty() && mRoi2.empty())
    return true;
  else {
    const int ix = (int) Point[0],
              iy = (int) Point[1],
              iz = (int) Point[2];
    bool isinmask = (ix > -1) && (ix < mNx) &&
                    (iy > -1) && (iy < mNy) &&
                    (iz > -1) && (iz < mNz);

    for (vector<MRI *>::const_iterator imask = mTestMask.begin();
                                       imask < mTestMask.end(); imask++)
      isinmask = isinmask && (MRIgetVoxVal(*imask, ix, iy, iz, 0) > 0);

    return isinmask;
  }
}

//
// Check that a point is inside one of the end ROIs
//
bool Blood::IsInRoi(vector<int>::const_iterator Point, MRI *Roi) {
  return (MRIgetVoxVal(Roi, Point[0], Point[1], Point[2], 0) > 0);
}

//
// Check that a point is inside a training subject's cortical mask
//
bool Blood::IsInCortex(vector<int>::const_iterator Point,
                       MRI *Mask, MRI *Aseg) {
  return (Point[0] > -1) && (Point[0] < mNx) &&
         (Point[1] > -1) && (Point[1] < mNy) &&
         (Point[2] > -1) && (Point[2] < mNz) &&
         (MRIgetVoxVal(Mask, Point[0], Point[1], Point[2], 0) > 0 ||
          (mMaskLabel > 0 &&
           MRIgetVoxVal(Aseg, Point[0], Point[1], Point[2], 0) == mMaskLabel));
}

//
// Read segmentations and cortical masks of training subjects
//
void Blood::ReadAnatomy(const std::string TrainListFile, const std::string TrainAsegFile,
			const std::string TrainMaskFile) {
  vector<string> dirlist;
  vector<MRI *>::iterator ivol;

  if (!TrainListFile.empty()) {		// Read multiple inputs from a list
    string dirname;
    ifstream listfile(TrainListFile, ios::in);

    if (!listfile) {
      cout << "ERROR: Could not open " << TrainListFile << endl;
      exit(1);
    }

    while (listfile >> dirname)
      dirlist.push_back(dirname + "/");
  }
  else				// Read a single input
    dirlist.push_back("");

  for (ivol = mAseg.begin(); ivol != mAseg.end(); ivol++)
    MRIfree(&(*ivol));
  mAseg.clear();  

  for (ivol = mMask.begin(); ivol != mMask.end(); ivol++)
    MRIfree(&(*ivol));
  mMask.clear();  

  // Clear all variables that depend on underlying anatomy
  mHistoLocal.clear();
  mHistoLocalAll.clear();
  mPriorLocal.clear();
  mPriorLocalAll.clear();
  mIdsLocal.clear();
  mIdsLocalAll.clear();
  mHistoNear.clear();
  mHistoNearAll.clear();
  mPriorNear.clear();
  mPriorNearAll.clear();
  mIdsNear.clear();
  mIdsNearAll.clear();
  mAsegDistMean.clear();
  mAsegDistMeanAll.clear();
  mAsegDistStd.clear();
  mAsegDistStdAll.clear();

  for (vector<string>::const_iterator idir = dirlist.begin();
                                      idir != dirlist.end(); idir++) {
    string fname;

    fname = *idir + TrainAsegFile;
    cout << "Loading segmentation from " << fname << endl;
    mAseg.push_back(MRIread(fname.c_str()));
    if (! *(mAseg.end()-1)) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }

    if (!TrainMaskFile.empty()) {
      fname = *idir + TrainMaskFile;
      cout << "Loading cortex mask from " << fname << endl;
      mMask.push_back(MRIread(fname.c_str()));
      if (! *(mMask.end()-1)) {
        cout << "ERROR: Could not read " << fname << endl;
        exit(1);
      }
    }
  }

  mNumTrain = mAseg.size();

  if (!mNumLines.empty() && ((int) mNumLines.size() != mNumTrain)) {
    cout << "ERROR: Number of training paths and segmentations does not match"
         << endl;
    exit(1);
  }
}

//
// Remove very short and very long streamlines
//
void Blood::RemoveLengthOutliers() {
  const int nlen = mLengths.size();
  int nrejlen = 0;

  if (nlen > 2) {
    const int lmin = *min_element(mLengths.begin(), mLengths.end()),
              lmax = *max_element(mLengths.begin(), mLengths.end());
    int lsum = 0, l2sum = 0, hthresh, llow = lmin, lhigh = lmax;
    float lmean, lstd;
    vector<int> lhisto(lmax-lmin+1, 0);
    vector<int>::const_iterator ihisto;

    // Find mean and standard deviation of lengths
    for (vector<int>::const_iterator ilen = mLengths.begin();
                                     ilen < mLengths.end(); ilen++) {
      lsum += *ilen;
      l2sum += (*ilen) * (*ilen);
    }

    lmean = lsum / (float) nlen;
    lstd = sqrt((l2sum - nlen*lmean*lmean) / (nlen-1));

    // Calculate histogram of lengths
    for (vector<int>::const_iterator ilen = mLengths.begin();
                                     ilen < mLengths.end(); ilen++)
      lhisto.at(*ilen - lmin)++;

    // How many streamlines are too few?
    hthresh = (int) ceil(0.03 * (*max_element(lhisto.begin(), lhisto.end())));

    // Find gap in lower half of histogram
    ihisto = lhisto.begin() + ((int) floor(lmean) - lmin);

    while (ihisto > lhisto.begin()) {
      if (*ihisto < hthresh) {
        int ngap = 1;
        vector<int>::const_iterator ithresh = ihisto;

        while (*(ihisto--) < hthresh && ihisto > lhisto.begin())
          ngap++;

        if (ngap > lstd) {
          llow = lmin + (ithresh + 1 - lhisto.begin());
          break;
        }
      }

      ihisto--;
    }
  
    // Find gap in upper half of histogram
    ihisto = lhisto.begin() + ((int) ceil(lmean) - lmin);

    while (ihisto < lhisto.end()) {
      if (*ihisto < hthresh) {
        int ngap = 1;
        vector<int>::const_iterator ithresh = ihisto;

        while (*(ihisto++) < hthresh && ihisto < lhisto.end())
          ngap++;

        if (ngap > lstd) {
          lhigh = lmin + (ithresh - 1 - lhisto.begin());
          break;
        }
      }

      ihisto++;
    }
  
    // Remove outlier streamlines
    vector<int>::iterator ilen = mLengths.begin();
    for (vector<int>::iterator inum = mNumLines.begin();
                               inum != mNumLines.end(); inum++)
      for (int k = *inum; k > 0; k--)
        if (*ilen < llow || *ilen > lhigh) {
          ilen = mLengths.erase(ilen);
          mStreamlines.erase(mStreamlines.begin() + (ilen - mLengths.begin()));
          (*inum)--;
          nrejlen++;
        }
        else
          ilen++;
  }

  cout << "INFO: Rejected " << nrejlen
       << " streamlines as length outliers" << endl;
}

//
// Do main computation for priors
//
void Blood::ComputePriors() {
  int maxtries;

  cout << "Matching streamline ends" << endl;
  MatchStreamlineEnds();

  SetArcSegments();

  cout << "Computing path histograms" << endl;
  ComputeHistogram();

  cout << "Computing prior on underlying anatomy "
       << "(non-truncated streamlines only)" << endl;
  ComputeAnatomyPrior(false);

  cout << "Computing prior on curvature "
       << "(non-truncated streamlines only)" << endl;
  ComputeShapePrior(false);

  if (mUseTruncated) {
    cout << "Computing prior on underlying anatomy "
         << "(all streamlines)" << endl;
    ComputeAnatomyPrior(true);

    cout << "Computing prior on curvature "
         << "(all streamlines)" << endl;
    ComputeShapePrior(true);
  }

  maxtries = mNumStrEnds - mExcludedStreamlines.size();

  for (int itry = 1; itry <= maxtries; itry++) {
    bool retry = false;

    mCenterStreamline.clear();
    mControlPoints.clear();
    mControlPointsTest.clear();
    mControlStd.clear();
    mControlStdTest.clear();

    cout << "Finding center streamline" << endl;
    FindCenterStreamline();

    for (vector<int>::const_iterator incpt = mNumControls.begin();
                                     incpt < mNumControls.end(); incpt++) {
      cout << "Selecting " << *incpt << " points on center streamline" << endl;
      //FindPointsOnStreamline(mCenterStreamline, *incpt);
      //if (!FindPointsOnStreamlineComb(mCenterStreamline, *incpt)) {
      if (!FindPointsOnStreamlineLS(mCenterStreamline, *incpt)) {
        cout << "WARN: Could not find satisfactory control point fit - try "
             << itry << endl;
        retry = true;

        if (itry < maxtries)
          break;
      }
    }

    if (retry)			// Try to find a better-behaved streamline
      mExcludedStreamlines.push_back(mCenterStreamline);
    else
      break;
  }
}

//
// Make sure all streamline start/end points are consistent
// and identify those that are within (a small distance of) the cortical mask
//
void Blood::MatchStreamlineEnds() {
  vector<bool>::iterator ivalid1, ivalid2;
  vector<int>::iterator itrlen;
  vector< vector<int> >::iterator istr;
  vector<MRI *>::const_iterator imask, iaseg, iroi1, iroi2;

  mIsInEnd1.resize(mStreamlines.size());
  mIsInEnd2.resize(mStreamlines.size());

  if (mRoi1.empty() || mRoi2.empty()) {		// Don't have labeling ROIs
    int imin;
    vector<int>::const_iterator ishort1, ishort2;

    fill(mIsInEnd1.begin(), mIsInEnd1.end(), true);
    fill(mIsInEnd2.begin(), mIsInEnd2.end(), true);

    // Find the shortest streamline
    imin = min_element(mLengths.begin(), mLengths.end()) - mLengths.begin();
    ishort1 = mStreamlines[imin].begin();
    ishort2 = mStreamlines[imin].end()-3;

    // Assign streamline starts/ends to match start/end of shortest streamline
    // (As a first pass, this will align most but not all streamlines correctly)
    for (istr = mStreamlines.begin(); istr != mStreamlines.end(); istr++) {
      const int dx1 = istr->at(0) - ishort1[0],
                dy1 = istr->at(1) - ishort1[1],
                dz1 = istr->at(2) - ishort1[2],
                dx2 = istr->at(0) - ishort2[0],
                dy2 = istr->at(1) - ishort2[1],
                dz2 = istr->at(2) - ishort2[2];

      if (dx1*dx1 + dy1*dy1 + dz1*dz1 >
          dx2*dx2 + dy2*dy2 + dz2*dz2)
        FlipStreamline(istr);
    }

    imask = mMask.begin();
    iaseg = mAseg.begin();
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

        // Do a second pass, aligning all streamlines with the majority
        for (jstr = mStreamlines.begin(); jstr != mStreamlines.end(); jstr++)
          if (jstr != istr) {
            vector<int>::const_iterator jend1 = jstr->begin(),
                                        jend2 = jstr->end() - 3;
            const int dx12 = iend1[0] - jend2[0],
                      dy12 = iend1[1] - jend2[1],
                      dz12 = iend1[2] - jend2[2],
                      dx11 = iend1[0] - jend1[0],
                      dy11 = iend1[1] - jend1[1],
                      dz11 = iend1[2] - jend1[2],
                      dx21 = iend2[0] - jend1[0],
                      dy21 = iend2[1] - jend1[1],
                      dz21 = iend2[2] - jend1[2],
                      dx22 = iend2[0] - jend2[0],
                      dy22 = iend2[1] - jend2[1],
                      dz22 = iend2[2] - jend2[2];

            // Is the start of path i closer to the start or the end of path j?
            dist1 += sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12)
                   - sqrt(dx11*dx11 + dy11*dy11 + dz11*dz11);

            // Is the end of path i closer to the start or the end of path j?
            dist2 += sqrt(dx21*dx21 + dy21*dy21 + dz21*dz21)
                   - sqrt(dx22*dx22 + dy22*dy22 + dz22*dz22);
          }

        if ( (dist1 < 0 && dist2 < 0) ||
             (dist1 < 0 && dist2 > 0 && dist2 < -dist1) ||
             (dist1 > 0 && dist2 < 0 && dist1 < -dist2) )
          FlipStreamline(istr);

        if (dist1 < 0 && dist2 > 0)
          *ivalid1 = false;
        else
          if (!mMask.empty())
            *ivalid1 = IsEnd1InMask(istr, *imask, *iaseg);

        if (dist1 > 0 && dist2 < 0)
          *ivalid2 = false;
        else
          if (!mMask.empty())
            *ivalid2 = IsEnd2InMask(istr, *imask, *iaseg);

        istr++;
        ivalid1++;
        ivalid2++;
      }

      if (!mMask.empty()) {
        imask++;
      }
      iaseg++;
    }
  }
  else {					// Have labeling ROIs
    fill(mIsInEnd1.begin(), mIsInEnd1.end(), false);
    fill(mIsInEnd2.begin(), mIsInEnd2.end(), false);

    iroi1 = mRoi1.begin();
    iroi2 = mRoi2.begin();
    istr = mStreamlines.begin();
    ivalid1 = mIsInEnd1.begin();
    ivalid2 = mIsInEnd2.begin();

    // Assign streamline starts/ends based on proximity to start/end ROIs
    // (As a first pass, this will align streamlines that go through both ROIs)
    for (vector<int>::const_iterator inum = mNumLines.begin();
                                     inum != mNumLines.end(); inum++) {
      for (int k = *inum; k > 0; k--) {
        int dist1 = 1, dist2 = 1;

        // Distance from the streamline start to the start ROI
        for (vector<int>::const_iterator ipt = istr->begin(); ipt < istr->end();
                                                              ipt += 3) {
          if (MRIgetVoxVal(*iroi1, ipt[0], ipt[1], ipt[2], 0) > 0) {
            *ivalid1 = true;
            break;
          }

          dist1++;
        }

        // Distance from the streamline start to the end ROI
        for (vector<int>::const_iterator ipt = istr->begin(); ipt < istr->end();
                                                              ipt += 3) {
          if (MRIgetVoxVal(*iroi2, ipt[0], ipt[1], ipt[2], 0) > 0) {
            *ivalid2 = true;
            break;
          }

          dist2++;
        }

        if (*ivalid1 && *ivalid2 && (dist1 > dist2))
          FlipStreamline(istr);

        istr++;
        ivalid1++;
        ivalid2++;
      }

      iroi1++;
      iroi2++;
    }

    ivalid1 = mIsInEnd1.begin();
    ivalid2 = mIsInEnd2.begin();

    // Align truncated streamlines to the ones that go through both ROIs
    for (istr = mStreamlines.begin(); istr != mStreamlines.end(); istr++) {
      int dist1 = 0, dist2 = 0;
      vector<int>::const_iterator iend1 = istr->begin(),
                                  iend2 = istr->end() - 3;
      vector<bool>::const_iterator jvalid1 = mIsInEnd1.begin(),
                                   jvalid2 = mIsInEnd2.begin();
      vector< vector<int> >::const_iterator jstr;

      if (*ivalid1 && !*ivalid2) {	// Streamline with truncated end point
        for (jstr = mStreamlines.begin(); jstr != mStreamlines.end(); jstr++) {
          if (*jvalid1 && *jvalid2) {
            vector<int>::const_iterator jend1 = jstr->begin();
            const int dx1 = iend1[0] - jend1[0],
                      dy1 = iend1[1] - jend1[1],
                      dz1 = iend1[2] - jend1[2],
                      dx2 = iend2[0] - jend1[0],
                      dy2 = iend2[1] - jend1[1],
                      dz2 = iend2[2] - jend1[2];

            // Distance from the start of path i to the start of path j
            dist1 += dx1*dx1 + dy1*dy1 + dz1*dz1;

            // Distance from the end of path i to the start of path j
            dist2 += dx2*dx2 + dy2*dy2 + dz2*dz2;
          }

          jvalid1++;
          jvalid2++;
        }

        if (dist1 > dist2)
          FlipStreamline(istr);
      }
      else if (!*ivalid1 && *ivalid2) {	// Streamline with truncated start point
        for (jstr = mStreamlines.begin(); jstr != mStreamlines.end(); jstr++) {
          if (*jvalid1 && *jvalid2) {
            vector<int>::const_iterator jend2 = jstr->end() - 3;
            const int dx1 = iend1[0] - jend2[0],
                      dy1 = iend1[1] - jend2[1],
                      dz1 = iend1[2] - jend2[2],
                      dx2 = iend2[0] - jend2[0],
                      dy2 = iend2[1] - jend2[1],
                      dz2 = iend2[2] - jend2[2];

            // Distance from the start of path i to the end of path j
            dist1 += dx1*dx1 + dy1*dy1 + dz1*dz1;

            // Distance from the end of path i to the end of path j
            dist2 += dx2*dx2 + dy2*dy2 + dz2*dz2;
          }

          jvalid1++;
          jvalid2++;
        }

        if (dist2 > dist1)
          FlipStreamline(istr);
      }

      ivalid1++;
      ivalid2++;
    }

    // Check if end points are (within a distance of) the cortical mask
    if (!mMask.empty()) {
      imask = mMask.begin();
      iaseg = mAseg.begin();
      istr = mStreamlines.begin();
      ivalid1 = mIsInEnd1.begin();
      ivalid2 = mIsInEnd2.begin();

      for (vector<int>::const_iterator inum = mNumLines.begin();
                                       inum != mNumLines.end(); inum++) {
        for (int k = *inum; k > 0; k--) {
          *ivalid1 = *ivalid1 && IsEnd1InMask(istr, *imask, *iaseg);
          *ivalid2 = *ivalid2 && IsEnd2InMask(istr, *imask, *iaseg);

          istr++;
          ivalid1++;
          ivalid2++;
        }

        imask++;
        iaseg++;
      }
    }
  }

  ComputeStatsEnds();

  ComputeEndPointCoM();

  // Map each truncated streamline to its nearest streamline that has
  // both start and end point in mask
  mTruncatedLengths.resize(mLengths.size());
  fill(mTruncatedLengths.begin(), mTruncatedLengths.end(), 0);

  if (mUseTruncated) {
    const int lag = max(1, (int) round(mHausStepRatio * mLengthAvgEnds)) * 3;

    ivalid1 = mIsInEnd1.begin();
    ivalid2 = mIsInEnd2.begin();
    itrlen  = mTruncatedLengths.begin();

    for (istr = mStreamlines.begin(); istr < mStreamlines.end(); istr++) {
      if ((*ivalid1 && !*ivalid2) || (!*ivalid1 && *ivalid2)) {
        double hdmin = numeric_limits<double>::infinity();
        vector<bool>::iterator jvalid1 = mIsInEnd1.begin(),
                               jvalid2 = mIsInEnd2.begin();
        vector< vector<int> >::const_iterator jstr, jstrnear;

        // Find nearest streamline that has both start and end point in mask
        for (jstr = mStreamlines.begin(); jstr < mStreamlines.end(); jstr++) {
          if (*jvalid1 && *jvalid2) {
            double hd = 0;

            for (vector<int>::const_iterator ipt = istr->begin();
                                             ipt < istr->end(); ipt += 3) {
              int dmin = 1000000;

              for (vector<int>::const_iterator jpt = jstr->begin();
                                               jpt < jstr->end(); jpt += lag) {
                const int dx = ipt[0] - jpt[0],
                          dy = ipt[1] - jpt[1],
                          dz = ipt[2] - jpt[2],
                          dist = dx*dx + dy*dy + dz*dz;

                if (dist < dmin)
                  dmin = dist;
              }

              hd += sqrt(dmin);
            }

            if (hd < hdmin) {
              hdmin = hd;
              jstrnear = jstr;
            }
          }

          jvalid1++;
          jvalid2++;
        }

        if (*ivalid1) {
          int dmin = 1000000;
          vector<int>::const_iterator jptnear, iend2 = istr->end() - 3;

          // Find point on whole streamline nearest to truncated end point
          for (vector<int>::const_iterator jpt = jstrnear->begin();
                                           jpt < jstrnear->end(); jpt += 3) {
            const int dx = iend2[0] - jpt[0],
                      dy = iend2[1] - jpt[1],
                      dz = iend2[2] - jpt[2],
                      dist = dx*dx + dy*dy + dz*dz;

            if (dist < dmin) {
              dmin = dist;
              jptnear = jpt;
            }
          }

          // Make educated guess about how much of streamline has been truncated
          *itrlen = (jstrnear->end() - 3 - jptnear) / 3;
        }

        if (*ivalid2) {
          int dmin = 1000000;
          vector<int>::const_iterator jptnear, iend1 = istr->begin();

          // Find point on whole streamline nearest to truncated start point
          for (vector<int>::const_iterator jpt = jstrnear->begin();
                                           jpt < jstrnear->end(); jpt += 3) {
            const int dx = iend1[0] - jpt[0],
                      dy = iend1[1] - jpt[1],
                      dz = iend1[2] - jpt[2],
                      dist = dx*dx + dy*dy + dz*dz;

            if (dist < dmin) {
              dmin = dist;
              jptnear = jpt;
            }
          }

          // Make educated guess about how much of streamline has been truncated
          *itrlen = (jptnear - jstrnear->begin()) / 3;
        }
      }

      ivalid1++;
      ivalid2++;
      itrlen++;
    }
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
// Check if a streamline has a valid (not truncated) start point
// by checking if the start point is within (a small distance of) cortex mask
//
bool Blood::IsEnd1InMask(vector< vector<int> >::iterator Streamline,
                         MRI *Mask, MRI *Aseg) {
  vector<int>::iterator itop = Streamline->begin();

  if (MRIgetVoxVal(Mask, itop[0], itop[1], itop[2], 0) > 0 ||
      (mMaskLabel > 0 &&
      MRIgetVoxVal(Aseg, itop[0], itop[1], itop[2], 0) == mMaskLabel))
    return true;

  if (itop < Streamline->end() - 3) {	// If more than one voxel in streamline
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

      if (!IsInMask(extend.begin()))	// Must also be in test subject's mask
        break;

      if (MRIgetVoxVal(Mask, newpt[0], newpt[1], newpt[2], 0) > 0 ||
          (mMaskLabel > 0 &&
          MRIgetVoxVal(Aseg, newpt[0], newpt[1], newpt[2], 0) == mMaskLabel)) {
        Streamline->insert(Streamline->begin(), extend.begin(), extend.end());
        mLengths[Streamline - mStreamlines.begin()] = Streamline->size() / 3;
        return true; 
      }
    }
  }

  return false;
}

//
// Check if a streamline has a valid (not truncated) end point
// by checking if the end point is within (a small distance of) cortex mask
//
bool Blood::IsEnd2InMask(vector< vector<int> >::iterator Streamline,
                         MRI *Mask, MRI *Aseg) {
  vector<int>::iterator ibottom = Streamline->end() - 3;

  if (MRIgetVoxVal(Mask, ibottom[0], ibottom[1], ibottom[2], 0) > 0 ||
      (mMaskLabel > 0 &&
      MRIgetVoxVal(Aseg, ibottom[0], ibottom[1], ibottom[2], 0) == mMaskLabel))
    return true;

  if (ibottom > Streamline->begin()) {	// If more than one voxel in streamline
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

      if (!IsInMask(extend.end()-3))	// Must also be in test subject's mask
        break;

      if (MRIgetVoxVal(Mask, newpt[0], newpt[1], newpt[2], 0) > 0 ||
          (mMaskLabel > 0 &&
          MRIgetVoxVal(Aseg, newpt[0], newpt[1], newpt[2], 0) == mMaskLabel)) {
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
  mNumArc = max(1, (int) round(*ilen / mLengthRatio));
/*
  mNumArc = max(1, (int) round(*min_element(lengths.begin(), lengths.end())
                        / mLengthRatio));
*/

  cout << "INFO: Split streamlines into " << mNumArc << " segments"  << endl;
}

//
// Compute spatial path histogram
//
void Blood::ComputeHistogram() {
  vector< vector<int> >::const_iterator istr;
  MRI *tmp;

  if (!mHistoStr || !mHistoSubj) return;

  // Compute streamline-wise histogram
  MRIclear(mHistoStr);
  mVolume = 0;

  for (istr = mStreamlines.begin(); istr != mStreamlines.end(); istr++)
    for (vector<int>::const_iterator ipt = istr->begin();
                                     ipt != istr->end(); ipt += 3) {
      const float nhits = MRIgetVoxVal(mHistoStr, ipt[0], ipt[1], ipt[2], 0);

      if (nhits < 1)
        mVolume++;

      MRIsetVoxVal(mHistoStr, ipt[0], ipt[1], ipt[2], 0, nhits + 1);
    }

  cout << "INFO: Total streamline volume is " << mVolume << " voxels" << endl;

  // Compute subject-wise histogram
  MRIclear(mHistoSubj);
  tmp = MRIclone(mHistoSubj, NULL);
  istr = mStreamlines.begin();

  for (vector<int>::const_iterator inum = mNumLines.begin();
                                   inum != mNumLines.end(); inum++) {
    MRIclear(tmp);

    for (int k = *inum; k > 0; k--) {
      for (vector<int>::const_iterator ipt = istr->begin();
                                       ipt != istr->end(); ipt += 3)
        MRIsetVoxVal(tmp, ipt[0], ipt[1], ipt[2], 0, 1);
      istr++;
    }

    MRIadd(mHistoSubj, tmp, mHistoSubj);
  }

  MRIfree(&tmp);
}

//
// Compute prior on underlying anatomy by streamline arc length
//
void Blood::ComputeAnatomyPrior(bool UseTruncated) {
  vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                               ivalid2 = mIsInEnd2.begin();
  vector<int>::const_iterator ilen = mLengths.begin(),
                              itrlen = mTruncatedLengths.begin();
  vector< vector<int> >::const_iterator istr = mStreamlines.begin();
  vector<MRI *>::const_iterator iaseg = mAseg.begin();
  vector< vector<int> > distbyarc;
  vector< vector<unsigned int> > localbyarc, nearbyarc;

  localbyarc.resize(mNumLocal * mNumArc);
  nearbyarc.resize(mNumNear * mNumArc);
  distbyarc.resize(mNumNear * mNumArc);

  // Find anatomical labels around each point on each streamline
  for (vector<int>::const_iterator inum = mNumLines.begin();
                                   inum != mNumLines.end(); inum++) {
    for (int k = *inum; k > 0; k--) {
      if ( (UseTruncated && (*ivalid1 || *ivalid2)) ||
           (*ivalid1 && *ivalid2) ) {
        unsigned int ilocal, inear;
        const double darc = mNumArc / (double) (*ilen + *itrlen);
        double larc;

        if (*ivalid1) {
          larc = 0;
          ilocal = 0;
          inear = 0;
        }
        else {				// Skip ahead to truncated start point
          double intpart;
          larc = modf((*itrlen) * darc, &intpart);	// decimal part
          ilocal = (unsigned int) intpart;		// integer part
          inear = ilocal;
        }

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
            const float seg = MRIgetVoxVal(*iaseg, 
                                           ((ix > -1 && ix < mNx) ? ix : ix0),
                                           ((iy > -1 && iy < mNy) ? iy : iy0),
                                           ((iz > -1 && iz < mNz) ? iz : iz0),
                                           0);

            localbyarc[ilocal].push_back((unsigned int) seg);

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

            nearbyarc[inear].push_back((unsigned int) seg);
            distbyarc[inear].push_back(dist);

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
      }

      ivalid1++;
      ivalid2++;
      ilen++;
      itrlen++;
      istr++;
    }

    iaseg++;
  }

  if (UseTruncated) {
    mIdsLocalAll.clear();
    mHistoLocalAll.clear();
    mPriorLocalAll.clear();
    mIdsNearAll.clear();
    mHistoNearAll.clear();
    mPriorNearAll.clear();
    mAsegDistMeanAll.clear();
    mAsegDistStdAll.clear();
  }
  else {
    mIdsLocal.clear();
    mHistoLocal.clear();
    mPriorLocal.clear();
    mIdsNear.clear();
    mHistoNear.clear();
    mPriorNear.clear();
    mAsegDistMean.clear();
    mAsegDistStd.clear();
  }

  // Compute priors on neighboring anatomical labels by arc length
  for (vector< vector<unsigned int> >::const_iterator
       iseg = localbyarc.begin(); iseg != localbyarc.end(); iseg++) {
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

    if (UseTruncated) {
      mIdsLocalAll.push_back(idlist);
      mHistoLocalAll.push_back(histo);
      mPriorLocalAll.push_back(prior);
    }
    else {
      mIdsLocal.push_back(idlist);
      mHistoLocal.push_back(histo);
      mPriorLocal.push_back(prior);
    }
  }

  vector< vector<int> >::const_iterator idist = distbyarc.begin();

  for (vector< vector<unsigned int> >::const_iterator
       iseg = nearbyarc.begin(); iseg != nearbyarc.end(); iseg++) {
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

    if (UseTruncated) {
      mIdsNearAll.push_back(idlist);
      mHistoNearAll.push_back(histo);
      mPriorNearAll.push_back(prior);
      mAsegDistMeanAll.push_back(dmean);
      mAsegDistStdAll.push_back(dstd);
    }
    else {
      mIdsNear.push_back(idlist);
      mHistoNear.push_back(histo);
      mPriorNear.push_back(prior);
      mAsegDistMean.push_back(dmean);
      mAsegDistStd.push_back(dstd);
    }

    idist++;
  }
}

//
// Compute prior on tangent vector and curvature by streamline arc length
//
void Blood::ComputeShapePrior(bool UseTruncated) {
  const int nbin = (int) ceil(2 / mTangentBinSize);
  float tangx, tangy, curv;
  vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                               ivalid2 = mIsInEnd2.begin();
  vector<int>::const_iterator ilen = mLengths.begin(),
                              itrlen = mTruncatedLengths.begin();
  vector<float>::const_iterator id1, id2;
  vector< vector<float> >::const_iterator itangx, itangy, icurv;
  vector<int> tanghisto(nbin*nbin), curvhisto;
  vector<float> strsmooth, diff1, diff2,
                tangprior(tanghisto.size()), curvprior;
  vector< vector<float> > tangxbyarc, tangybyarc, curvbyarc;

  tangxbyarc.resize(mNumArc);
  tangybyarc.resize(mNumArc);
  curvbyarc.resize(mNumArc);

  // Compute tangent vector and curvature at each point on each streamline
  for (vector< vector<int> >::const_iterator istr = mStreamlines.begin();
                                            istr < mStreamlines.end(); istr++) {
    if ( ( (UseTruncated && (*ivalid1 || *ivalid2)) ||
           (*ivalid1 && *ivalid2) ) && *ilen >= (int) mDiffStep ) {
      unsigned int iarc;
      const double darc = mNumArc / (double) (*ilen + *itrlen);
      double larc;

      strsmooth.resize(istr->size());
      diff1.resize(istr->size());
      diff2.resize(istr->size());

      // Smooth discrete point coordinates
      CurveSmooth(strsmooth, *istr);

      // Approximate first derivative by smoothed finite differences
      // of the point coordinates
      CurveFiniteDifferences(diff1, strsmooth, mDiffStep);

      // Approximate second derivative by smoothed finite differences
      // of the first derivative
      CurveFiniteDifferences(diff2, diff1, mDiffStep);

      id2 = diff2.begin();

      if (*ivalid1) {
        larc = 0;
        iarc = 0;
      }
      else {				// Skip ahead to truncated start point
        double intpart;
        larc = modf((*itrlen) * darc, &intpart);	// decimal part
        iarc = (unsigned int) intpart;			// integer part
      }

      for (id1 = diff1.begin(); id1 < diff1.end(); id1 += 3) {
        const float nrm = sqrt(id1[0]*id1[0] + id1[1]*id1[1] + id1[2]*id1[2]);

        if (nrm > 0) {
          // Tangent vector
          if (id1[2] < 0) {	// Use only z>0 octants, direction is irrelevant
            tangx = - id1[0] / nrm;
            tangy = - id1[1] / nrm;
          }
          else {
            tangx = id1[0] / nrm;
            tangy = id1[1] / nrm;
          }

          // Curvature = |r' x r''| / |r'|^3
          curv = 
            sqrt( pow(id1[1] * id2[2] - id1[2] * id2[1], 2) +
                  pow(id1[2] * id2[0] - id1[0] * id2[2], 2) +
                  pow(id1[0] * id2[1] - id1[1] * id2[0], 2) ) / pow(nrm, 3);
        }
        else {
          tangx = 0;
          tangy = 0;
          curv = 0;
        }

        tangxbyarc[iarc].push_back(tangx);
        tangybyarc[iarc].push_back(tangy);
        curvbyarc[iarc].push_back(curv);

        id2 += 3;

        larc += darc;

        if (larc > 1) {		// Move to the next segment
          larc -= 1;
          iarc++;
        }
      }
    }

    ivalid1++;
    ivalid2++;
    ilen++;
    itrlen++;
  }

  if (UseTruncated) {
    mHistoTangentAll.clear();
    mPriorTangentAll.clear();
    mHistoCurvatureAll.clear();
    mPriorCurvatureAll.clear();
  }
  else {
    mHistoTangent.clear();
    mPriorTangent.clear();
    mHistoCurvature.clear();
    mPriorCurvature.clear();
  }

  // Compute priors on tangent vector and curvature by arc length
  itangx = tangxbyarc.begin();
  itangy = tangybyarc.begin();

  for (icurv = curvbyarc.begin(); icurv < curvbyarc.end(); icurv++) {
    int nclass;
    const unsigned int nsamp = icurv->size();
    float denom;
    const float curvmax = *max_element(icurv->begin(), icurv->end());
    vector<int>::const_iterator ihisto;
    vector<float>::const_iterator isampx = itangx->begin();

    // Tangent vector histogram (defined over [-1, 1]^2)
    fill(tanghisto.begin(), tanghisto.end(), 0);

    for (vector<float>::const_iterator isampy = itangy->begin();
                                       isampy < itangy->end(); isampy++) {
      const int ix = (int) floor((*isampx + 1) / mTangentBinSize),
                iy = (int) floor((*isampy + 1) / mTangentBinSize);

      tanghisto.at( ((iy<nbin) ? iy : (nbin-1)) * nbin +
                    ((ix<nbin) ? ix : (nbin-1)) )++;

      isampx++;
    }

    // Tangent vector prior
    nclass = tanghisto.size() - count(tanghisto.begin(), tanghisto.end(), 0);
    denom = (float) (nsamp + nclass + 1);

    ihisto = tanghisto.begin();

    for (vector<float>::iterator iprior = tangprior.begin();
                                 iprior < tangprior.end(); iprior++) {
      *iprior = -log((*ihisto + 1) / denom);
      ihisto++;
    }

    // Curvature histogram (defined over [0, curvmax])
    curvhisto.resize((int) ceil(curvmax / mCurvatureBinSize) + 1);
    fill(curvhisto.begin(), curvhisto.end(), 0);

    for (vector<float>::const_iterator isamp = icurv->begin();
                                       isamp < icurv->end(); isamp++)
      curvhisto.at((int) floor(*isamp / mCurvatureBinSize))++;

    // Curvature prior
    nclass = curvhisto.size() - count(curvhisto.begin(), curvhisto.end(), 0);
    denom = (float) (nsamp + nclass + 1);

    curvprior.resize(curvhisto.size());
    ihisto = curvhisto.begin();

    for (vector<float>::iterator iprior = curvprior.begin();
                                 iprior < curvprior.end(); iprior++) {
      *iprior = -log((*ihisto + 1) / denom);
      ihisto++;
    }

    if (UseTruncated) {
      mHistoTangentAll.push_back(tanghisto);
      mPriorTangentAll.push_back(tangprior);
      mHistoCurvatureAll.push_back(curvhisto);
      mPriorCurvatureAll.push_back(curvprior);
    }
    else {
      mHistoTangent.push_back(tanghisto);
      mPriorTangent.push_back(tangprior);
      mHistoCurvature.push_back(curvhisto);
      mPriorCurvature.push_back(curvprior);
    }

    itangx++;
    itangy++;
  }

  // In debug mode, save all samples of the tangent vector and curvature
  // by arc length
  if (mDebug) {
    if (UseTruncated) {
      mTangentXByArcAll.resize(mNumArc);
      copy(tangxbyarc.begin(), tangxbyarc.end(), mTangentXByArcAll.begin());
      mTangentYByArcAll.resize(mNumArc);
      copy(tangybyarc.begin(), tangybyarc.end(), mTangentYByArcAll.begin());
      mCurvatureByArcAll.resize(mNumArc);
      copy(curvbyarc.begin(), curvbyarc.end(), mCurvatureByArcAll.begin());
    }
    else {
      mTangentXByArc.resize(mNumArc);
      copy(tangxbyarc.begin(), tangxbyarc.end(), mTangentXByArc.begin());
      mTangentYByArc.resize(mNumArc);
      copy(tangybyarc.begin(), tangybyarc.end(), mTangentYByArc.begin());
      mCurvatureByArc.resize(mNumArc);
      copy(curvbyarc.begin(), curvbyarc.end(), mCurvatureByArc.begin());
    }
  }
}

//
// Find central streamline among streamlines with valid end points
//
void Blood::FindCenterStreamline(bool CheckOverlap, bool CheckDeviation,
                                                    bool CheckFa) {
  const int lag = max(1, (int) round(mHausStepRatio * mLengthAvgEnds)) * 3;
  double hdmin = numeric_limits<double>::infinity();
  vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                               ivalid2 = mIsInEnd2.begin();
  vector<int>::const_iterator imidpts = mMidPoints.begin();
  vector< vector<int> >::const_iterator icenter;

  if (mStreamlines.empty() || mNumStrEnds == 0)
    return;

  if (mNumStrEnds == 1)		// Only one candidate streamline
    for (vector< vector<int> >::const_iterator istr = mStreamlines.begin();
                                               istr < mStreamlines.end();
                                               istr++) {
      if (*ivalid1 && *ivalid2) {
        mCenterStreamline = *istr;

        cout << "INFO: Length of center streamline is "
             << mCenterStreamline.size()/3 << " voxels" << endl;

        return;
      }

      ivalid1++;
      ivalid2++;
    }

  cout << "INFO: Step is " << lag/3 << " voxels" << endl;

  for (vector< vector<int> >::const_iterator istr = mStreamlines.begin();
                                             istr < mStreamlines.end();
                                             istr++) {
    if (*ivalid1 && *ivalid2) {
      bool isexcluded = false,
           okhist = true, okfa = true,
           okend1 = true, okend2 = true, okmid = true;
      double hdtot = 0;
      vector<bool>::const_iterator jvalid1 = mIsInEnd1.begin(),
                                   jvalid2 = mIsInEnd2.begin();
      vector<int>::const_iterator jlen = mLengths.begin();

      // Check if this is one of the excluded streamlines, if any
      for (vector< vector<int> >::const_iterator
                                  ixstr = mExcludedStreamlines.begin();
                                  ixstr < mExcludedStreamlines.end(); ixstr++)
        if (equal(istr->begin(), istr->end(), ixstr->begin()))
          isexcluded = true;

      if (!isexcluded && mNumTrain > 1) {	// No checks for single subject
        if (CheckOverlap) {		// Check overlap with histogram
          int nzeros = 0;

          for (vector<int>::const_iterator ipt = istr->begin();
                                           ipt < istr->end(); ipt += 3) {
            const float h = MRIgetVoxVal(mHistoStr, ipt[0], ipt[1], ipt[2], 0);

            if (h < 4.0)		// Little overlap with other streamlines
              nzeros++;
          }

          if (nzeros >= (int) (.1 * istr->size()/3))
            okhist = false;
        }

        if (CheckFa) {		// Check overlap with test subject's FA
          vector<int> nfzeros(mTestFa.size(), 0), basept(3), dwipt(3);

          for (vector<int>::const_iterator ipt = istr->begin();
                                           ipt < istr->end(); ipt += 3) {
            vector<int>::iterator infzeros;
            vector<int>::const_iterator iptbase, iptdwi;

            // Map point to base space if needed
            if (!mTestAffineReg.IsEmpty() || !mTestNonlinReg.IsEmpty()) {
              if (!MapPointToBase(basept.begin(), ipt)) {
                for (infzeros = nfzeros.begin(); infzeros < nfzeros.end();
                                                 infzeros++) {
                  (*infzeros)++;

                  if (*infzeros > 3) {
                    okfa = false;
                    break;
                  }
                }

                if (!okfa)
                  break;

                continue;
              }

              iptbase = basept.begin();
            }
            else
              iptbase = ipt;

            infzeros = nfzeros.begin();

            for (vector<MRI *>::const_iterator ifa = mTestFa.begin();
                                               ifa < mTestFa.end(); ifa++) {
              // Map point to native space if needed
              if (!mTestBaseReg.empty()) {
                if (!MapPointToNative(dwipt.begin(), iptbase,
                                      ifa - mTestFa.begin())) {
                  (*infzeros)++;

                  if (*infzeros > 3) {
                    okfa = false;
                    break;
                  }

                  infzeros++;
                  continue;
                }

                iptdwi = dwipt.begin();
              }
              else
                iptdwi = iptbase;

              // Check anisotropy at this point
              if (MRIgetVoxVal(*ifa, iptdwi[0], iptdwi[1], iptdwi[2], 0)
                                                                      < 0.1) {
                (*infzeros)++;

                if (*infzeros > 3) {
                  okfa = false;
                  break;
                }
              }
              else
                *infzeros = 0;

              infzeros++;
            }

            if (!okfa)
              break;
          }
        }

        if (CheckDeviation) {	// Check endpoint deviation from center of mass
          vector<int>::const_iterator itop    = istr->begin(),
                                      ibottom = istr->end() - 3,
                                      imiddle = itop + *imidpts;

          for (int k = 0; k < 3; k++) {
            const float dist = itop[k] - mMeanEnd1[k];
            okend1 = okend1 && (dist*dist < mVarEnd1[k]);
          }

          for (int k = 0; k < 3; k++) {
            const float dist = ibottom[k] - mMeanEnd2[k];
            okend2 = okend2 && (dist*dist < mVarEnd2[k]);
          }

          for (int k = 0; k < 3; k++) {
            const float dist = imiddle[k] - mMeanMid[k];
            okmid = okmid && (dist*dist < mVarMid[k]);
          }
        }
      }

      if (!isexcluded && okhist && okfa && okend1 && okend2 && okmid) {
        for (vector< vector<int> >::const_iterator
             jstr = mStreamlines.begin(); jstr < mStreamlines.end(); jstr++) {
          double hd = 0;

          if (*jvalid1 && *jvalid2 && (jstr != istr)) {
            for (vector<int>::const_iterator jpt = jstr->begin();
                                             jpt < jstr->end(); jpt += lag) {
              int dmin = 1000000;

              for (vector<int>::const_iterator ipt = istr->begin();
                                               ipt < istr->end(); ipt += lag) {
                const int dx = ipt[0] - jpt[0],
                          dy = ipt[1] - jpt[1],
                          dz = ipt[2] - jpt[2],
                          dist = dx*dx + dy*dy + dz*dz;

                if (dist < dmin)
                  dmin = dist;
              }

              hd += sqrt(dmin);
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

      imidpts++;
    }

    ivalid1++;
    ivalid2++;
  }

  if (hdmin == numeric_limits<double>::infinity()) {
    // In case checks caused failure
    if (CheckDeviation) {
      cout << "WARN: Turning off deviation check for center streamline" << endl;
      FindCenterStreamline(CheckOverlap, false, CheckFa);
    }
    else if (CheckFa) {
      cout << "WARN: Turning off FA check for center streamline" << endl;
      FindCenterStreamline(CheckOverlap, false, false);
    }
    else if (CheckOverlap) {
      cout << "WARN: Turning off overlap check for center streamline" << endl;
      FindCenterStreamline(false, false, false);
    } else {
      cout << "WARN: All checks already off. Exiting." << endl;
      return;
    }
  }
  else {
    mCenterStreamline = *icenter;

    cout << "INFO: Length of center streamline is "
         << mCenterStreamline.size()/3 << " voxels" << endl;
  }
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
  vector<int> sum(3, 0), nptseg, cpts(NumPoints*3);
  vector<double> sumsq(3, 0), ptdist(nptot-1, 0), lenseg, nptsegdec;
  vector<vector<int>::const_iterator> cptopt;

  Spline spline(Streamline, mTestMask[0]); 

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
      sumsq[k] += (ipt[k]*ipt[k]);
    }

  for (int k = 0; k < 3; k++)
    sumsq[k] -= (sum[k]*sum[k] / (double) nptot);

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

  cptopt.push_back(ipt);

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

      cptopt.push_back(ipt);
    }

    ilen++;
  }

  for (int iter = 0; iter < 100; iter++) {
    vector<int>::iterator icpt = cpts.begin();
    vector<float> overlap(NumPoints-1, 0);

    // Fit spline to current control points
    for (vector<vector<int>::const_iterator>::const_iterator
         iopt = cptopt.begin(); iopt < cptopt.end(); iopt++) {
      copy(*iopt, (*iopt)+3, icpt);
      icpt += 3;
    }

    spline.SetControlPoints(cpts);
    if (!spline.InterpolateSpline())
     {}

    // Find overlap of fitted spline segments between controls with histogram
    ipt = spline.GetAllPointsBegin() + 3;
    icpt = cpts.begin() + 3;

    for (vector<float>::iterator iover = overlap.begin();
                                 iover < overlap.end(); iover++) {
      int N = 0;

      while (!(ipt[0] == icpt[0] && ipt[1] == icpt[1] && ipt[2] == icpt[2])) {
        *iover += MRIgetVoxVal(mHistoStr, ipt[0], ipt[1], ipt[2], 0);
        N++;
        ipt += 3;
      }

      *iover /= N;
      ipt += 3;
      icpt += 3;
    }

    // Move a control point towards the segment with the least overlap
    int iworst = min_element(overlap.begin(), overlap.end()) - overlap.begin();
    int imove;
    bool moveback;

    if (iworst == 0) {
      imove = 1;
      moveback = true;
    }
    else if (iworst == NumPoints-2) {
      imove = NumPoints-2;
      moveback = false;
    }
    else {
      if (overlap[iworst-1] > overlap[iworst+1]) {
        imove = iworst;
        moveback = false;
      }
      else {
        imove = iworst+1;
        moveback = true;
      }
    }

    if (moveback && (cptopt[imove] - 3 != cptopt[imove-1]))
      cptopt[imove] -= (int) round((cptopt[imove] - cptopt[imove-1])/6) * 3;
    else if (!moveback && (cptopt[imove] + 3 != cptopt[imove+1]))
      cptopt[imove] += (int) round((cptopt[imove+1] - cptopt[imove])/6) * 3;
  }

  mControlPoints.push_back(cpts);

  ComputeStreamlineSpread(cpts);
}

//
// Find spread of streamlines around each control point
// (in base space, if registration has been specified)
//
void Blood::ComputeStreamlineSpread(vector<int> &ControlPoints) {
  vector<int> point(3, 0), dmin(ControlPoints.size()/3),
                           ptmin(ControlPoints.size(), 0),
                           sum(ControlPoints.size(), 0),
                           sumsq(ControlPoints.size(), 0);
  vector<float> cstd(ControlPoints.size(), 0);

  if (mNumStrEnds > 1) {
    bool isinbase = true;
    vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                                 ivalid2 = mIsInEnd2.begin();
    vector<int>::const_iterator isum   = sum.begin(),
                                isumsq = sumsq.begin();

    for (vector< vector<int> >::const_iterator istr = mStreamlines.begin();
                                               istr < mStreamlines.end();
                                               istr++) {
      if (*ivalid1 && *ivalid2) {
        vector<int>::iterator isum   = sum.begin(),
                              isumsq = sumsq.begin();

        // Find closest point on streamline to each control point in base space
        fill(dmin.begin(), dmin.end(), 1000000);

        for (vector<int>::const_iterator ipt = istr->begin(); ipt < istr->end();
                                                              ipt += 3) {
          vector<int>::iterator idmin  = dmin.begin(),
                                iptmin = ptmin.begin();

          // Map point to base space if registration has been specified
          if (!mTestAffineReg.IsEmpty() || !mTestNonlinReg.IsEmpty())
            isinbase = MapPointToBase(point.begin(), ipt);

          if (!isinbase)
            continue;

          // Find point's distance from each control point in base space
          for (vector<int>::const_iterator icpt = ControlPoints.begin();
                                           icpt < ControlPoints.end();
                                           icpt += 3) {
            int dist = 0;

            for (int k = 0; k < 3; k++) {
              const int diff = point[k] - icpt[k];
              dist += diff*diff;
            }

            if (dist < *idmin) {
              *idmin = dist;
              copy(point.begin(), point.end(), iptmin);
            }

            idmin++;
            iptmin += 3;
          }
        }

        // Compute variance of nearest streamline points for each control point
        for (vector<int>::const_iterator iptmin = ptmin.begin();
                                         iptmin < ptmin.end(); iptmin++) {
          *isum += *iptmin;
          *isumsq += (*iptmin) * (*iptmin);

          isum++;
          isumsq++;
        }
      }

      ivalid1++;
      ivalid2++;
    }

    for (vector<float>::iterator istd = cstd.begin(); istd < cstd.end();
                                                      istd++) {
      *istd = sqrt((*isumsq - (*isum) * (*isum / float(mNumStrEnds)))
                   / (mNumStrEnds-1));
      isum++;
      isumsq++;
    }
  }

  if (mTestAffineReg.IsEmpty() && mTestNonlinReg.IsEmpty())
    mControlStd.push_back(cstd);
  else
    mControlStdTest.push_back(cstd);
}

//
// Select control points on a streamline
// to maximize overlap of the fitted spline with the streamline histogram
//
bool Blood::FindPointsOnStreamlineComb(vector<int> &Streamline, int NumPoints) {
  bool success = true;
  const int strlen = (int) Streamline.size() / 3,
            strdiv = (int) round((strlen-1) / (NumPoints-1)),
            lag = max(1, min((int) round(mControlStepRatio * NumPoints / mDx),
                             strdiv)) * 3;
  double hdmin = numeric_limits<double>::infinity();
  vector<int> cpts(NumPoints*3);
  vector<int>::const_iterator ipt;
  vector<vector<int>::const_iterator> cptopt(NumPoints);
  Spline spline(Streamline, mTestMask[0]); 

  if (NumPoints > strlen) {
    cout << "ERROR: Selected streamline has fewer than " << NumPoints
         << " points" << endl;
    exit(1);
  }

  // If all else fails, I'll use equidistant points
  ipt = Streamline.begin();
  for (vector<int>::iterator icpt = cpts.begin(); icpt < cpts.end() - 3;
                                                  icpt += 3) {
    copy(ipt, ipt+3, icpt);
    ipt += strdiv*3;
  }
  copy(Streamline.end()-3, Streamline.end(), cpts.end()-3);

  cout << "INFO: Step is " << lag/3 << " voxels" << endl;

  // Keep first and last control points fixed 
  cptopt[0] = Streamline.begin();
  cptopt[NumPoints-1] = Streamline.end() - 3;

  // Try combinations of intermediate control points
  for (cptopt[1] = Streamline.begin() + lag;
       cptopt[1] <= Streamline.end() - 3 - (NumPoints-2) * lag;
       cptopt[1] += lag)
    TryControlPoint(hdmin, 2, lag, cpts, cptopt, spline, Streamline);

  if (hdmin == numeric_limits<double>::infinity()) {
    cout << "WARN: Defaulting to equidistant control points" << endl;
    success = false;
  }

  cout << "INFO: Selected control points are" << endl;
  for (vector<int>::const_iterator icpt = cpts.begin(); icpt < cpts.end();
                                                        icpt += 3)
    cout << " " << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;

  cout << "INFO: Distances between consecutive points are";
  for (vector<int>::const_iterator icpt = cpts.begin() + 3; icpt < cpts.end();
                                                            icpt += 3) {
    const int dx = icpt[0] - icpt[-3],
              dy = icpt[1] - icpt[-2],
              dz = icpt[2] - icpt[-1];

    cout << " " << round(sqrt(dx*dx + dy*dy + dz*dz));
  }
  cout << endl;

  mControlPoints.push_back(cpts);

  // Map control points to base space if needed
  if (!mTestAffineReg.IsEmpty() || !mTestNonlinReg.IsEmpty()) {
    vector<int> cptsout(cpts.size()), point(3);
    vector<int>::iterator icptout = cptsout.begin();

    for (vector<int>::const_iterator icpt = cpts.begin(); icpt < cpts.end();
                                                          icpt += 3) {
      MapPointToBase(point.begin(), icpt);
      copy(point.begin(), point.end(), icptout);

      icptout += 3;
    }

    cout << "INFO: Selected control points in test subject's space are" << endl;
    for (vector<int>::const_iterator icpt = cptsout.begin();
                                     icpt < cptsout.end(); icpt += 3)
      cout << " " << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;

    cout << "INFO: Distances between consecutive points "
         << "in test subject's space are";
    for (vector<int>::const_iterator icpt = cptsout.begin() + 3;
                                     icpt < cptsout.end(); icpt += 3) {
      const int dx = icpt[0] - icpt[-3],
                dy = icpt[1] - icpt[-2],
                dz = icpt[2] - icpt[-1];

      cout << " " << round(sqrt(dx*dx + dy*dy + dz*dz));
    }
    cout << endl;

    mControlPointsTest.push_back(cptsout);

    ComputeStreamlineSpread(cptsout);
  }
  else
    ComputeStreamlineSpread(cpts);

  return success;
}

//
// Pick combinations of intermediate control points recursively
//
void Blood::TryControlPoint(double &HausDistMin,
                            int IndexPoint,
                            int SearchLag,
                            vector<int> &ControlPointsOpt,
                            vector<vector<int>::const_iterator> &ControlPoints,
                            Spline &TrySpline,
                            vector<int> &Streamline) {
  const int ncpts = ControlPoints.size();

  if (IndexPoint < ncpts-2)
    for (ControlPoints[IndexPoint] = ControlPoints[IndexPoint-1] + SearchLag;
         ControlPoints[IndexPoint] <= Streamline.end() - 3
                                      - (ncpts-IndexPoint-1) * SearchLag;
         ControlPoints[IndexPoint] += SearchLag)
      TryControlPoint(HausDistMin, IndexPoint + 1, SearchLag,
                      ControlPointsOpt, ControlPoints, TrySpline, Streamline);
  else
    for (ControlPoints[IndexPoint] = ControlPoints[IndexPoint-1] + SearchLag;
         ControlPoints[IndexPoint] < Streamline.end() - SearchLag;
         ControlPoints[IndexPoint] += SearchLag) {
      bool okfa = true;
      int splen, dmin = 1000000, dminmax = 0, nhzeros = 0;
      double hd = 0.0;
      vector<int> cpts, nfzeros(mTestFa.size(), 0), point(3);
      Spline basespline, *checkspline;

      // Fit spline to current control points
      for (vector<vector<int>::const_iterator>::const_iterator
           icpt = ControlPoints.begin(); icpt < ControlPoints.end(); icpt++)
        cpts.insert(cpts.end(), *icpt, (*icpt)+3);

      TrySpline.SetControlPoints(cpts);
      if (!TrySpline.InterpolateSpline())
        continue;

      splen = (TrySpline.GetAllPointsEnd() - TrySpline.GetAllPointsBegin()) / 3;

      // Find Hausdorff distance of true streamline from fitted spline
      for (vector<int>::const_iterator ipt = TrySpline.GetAllPointsBegin();
                                       ipt < TrySpline.GetAllPointsEnd();
                                       ipt += 3) {
        const float h = MRIgetVoxVal(mHistoStr, ipt[0], ipt[1], ipt[2], 0);
        dmin = 1000000;

        // Point distance from true streamline
        for (vector<int>::const_iterator iptrue = Streamline.begin();
                                         iptrue < Streamline.end();
                                         iptrue += 3) {
          const int dx = ipt[0] - iptrue[0],
                    dy = ipt[1] - iptrue[1],
                    dz = ipt[2] - iptrue[2],
                    dist = dx*dx + dy*dy + dz*dz;

          if (dist < dmin)
            dmin = dist;
        }

        if (dmin > 25)			// Point is far from true streamline
          break;

        if (dmin > dminmax)
          dminmax = dmin;

        if (h < 4.0)			// Point is off histogram
          nhzeros++;
      }

      // Don't allow spline if any point is too far from true streamline
      // or if more than 10% of all points are off the histogram
      if (dmin > 25 || nhzeros > (int) (.1 * splen))
        continue;

      // Map control points and fit spline in base space if needed
      if (!mTestAffineReg.IsEmpty() || !mTestNonlinReg.IsEmpty()) {
        vector<int> basecpts(cpts.size());
        vector<int>::iterator icptout = basecpts.begin();

        for (vector<int>::const_iterator icpt = cpts.begin(); icpt < cpts.end();
                                                              icpt += 3) {
          okfa = MapPointToBase(point.begin(), icpt);

          if (!okfa)
            break;

          copy(point.begin(), point.end(), icptout);

          icptout += 3;
        }

        if (!okfa)
          continue;

        basespline.SetMask(mTestBaseMask);
        basespline.SetControlPoints(basecpts);
        okfa = basespline.InterpolateSpline();

        if (!okfa)
          continue;

        checkspline = &basespline;
      }
      else
        checkspline = &TrySpline;

      // Check anisotropy along spline
      for (vector<int>::const_iterator ipt = checkspline->GetAllPointsBegin();
                                       ipt < checkspline->GetAllPointsEnd();
                                       ipt += 3) {
        vector<int>::iterator infzeros = nfzeros.begin();
        vector<int>::const_iterator iptdwi;

        for (vector<MRI *>::const_iterator ifa = mTestFa.begin();
                                           ifa < mTestFa.end(); ifa++) {
          // Map point to native space if needed
          if (!mTestBaseReg.empty()) {
            if (!MapPointToNative(point.begin(), ipt, ifa - mTestFa.begin())) {
              (*infzeros)++;		// Point is outside anisotropy map

              if (*infzeros > 3) {
                okfa = false;
                break;
              }

              infzeros++;
              continue;
            }

            iptdwi = point.begin();
          }
          else
            iptdwi = ipt;

          // Check anisotropy at this point
          if (MRIgetVoxVal(*ifa, iptdwi[0], iptdwi[1], iptdwi[2], 0) < 0.1) {
            (*infzeros)++;		// Point is in low anisotropy area

            if (*infzeros > 3) {
              okfa = false;
              break;
            }
          }
          else
            *infzeros = 0;

          infzeros++;
        }

        if (!okfa)
          break;
      }

      // Don't allow spline if more than 3 contiguous points are in low FA
      if (!okfa)
        continue;

      hd = sqrt(dminmax);

      if (hd < HausDistMin) {
        copy(cpts.begin(), cpts.end(), ControlPointsOpt.begin());
        HausDistMin = hd;
      }
    }
}

//
// Least-squares fit of spline control points to streamline
//
bool Blood::FindPointsOnStreamlineLS(vector<int> &Streamline, int NumPoints) {
  bool success;
  const int strlen = (int) Streamline.size() / 3;
  vector<int> cpts(NumPoints*3), cptsout;
  vector<int>::const_iterator ipt;
  Spline spline(NumPoints, mTestMask[0]); 

  if (NumPoints > strlen) {
    cout << "ERROR: Selected streamline has fewer than " << NumPoints
         << " points" << endl;
    exit(1);
  }

  // Least-squares fitting of spline control points
  success = spline.FitControlPoints(Streamline);

  // Don't allow spline if any point is too far from true streamline
  // or if more than 10% of all points are off the histogram
  if (success) {
    const int splen = (spline.GetAllPointsEnd() -
                       spline.GetAllPointsBegin()) / 3;
    int nhzeros = 0;

    // Find Hausdorff distance of true streamline from fitted spline
    for (vector<int>::const_iterator ipt = spline.GetAllPointsBegin();
                                     ipt < spline.GetAllPointsEnd();
                                     ipt += 3) {
      const float h = MRIgetVoxVal(mHistoStr, ipt[0], ipt[1], ipt[2], 0);
      int dmin = 1000000;

      // Point distance from true streamline
      for (vector<int>::const_iterator iptrue = Streamline.begin();
                                       iptrue < Streamline.end();
                                       iptrue += 3) {
        const int dx = ipt[0] - iptrue[0],
                  dy = ipt[1] - iptrue[1],
                  dz = ipt[2] - iptrue[2],
                  dist = dx*dx + dy*dy + dz*dz;

        if (dist < dmin)
          dmin = dist;
      }

      if (dmin > 25) {			// Point is far from true streamline
        success = false;
        break;
      }

      if (h < 4.0)			// Point is off histogram
        nhzeros++;
    }

    if (nhzeros > (int) (.1 * splen))
      success = false;
  }

  // Don't allow spline if more than 3 contiguous points are in low FA
  if (success) {
    vector<int> nfzeros(mTestFa.size(), 0), point(3);
    Spline basespline, *checkspline;

    // Map control points and fit spline in base space if needed
    if (!mTestAffineReg.IsEmpty() || !mTestNonlinReg.IsEmpty()) {
      vector<int>::iterator icptout;

      cptsout.resize(cpts.size());
      icptout = cptsout.begin();

      for (vector<int>::const_iterator icpt = spline.GetControlPointsBegin();
                                       icpt < spline.GetControlPointsEnd();
                                       icpt += 3) {
        success = MapPointToBase(point.begin(), icpt);

        if (!success)
          break;

        copy(point.begin(), point.end(), icptout);

        icptout += 3;
      }

      if (success) {
        basespline.SetMask(mTestBaseMask); 
        basespline.SetControlPoints(cptsout);
        success = basespline.InterpolateSpline();
      }

      checkspline = &basespline;
    }
    else
      checkspline = &spline;

    // Check anisotropy along the fitted spline
    if (success)
      for (vector<int>::const_iterator ipt = checkspline->GetAllPointsBegin();
                                       ipt < checkspline->GetAllPointsEnd();
                                       ipt += 3) {
        vector<int>::iterator infzeros = nfzeros.begin();
        vector<int>::const_iterator iptdwi;

        for (vector<MRI *>::const_iterator ifa = mTestFa.begin();
                                           ifa < mTestFa.end(); ifa++) {
          // Map point to native space if needed
          if (!mTestBaseReg.empty()) {
            if (!MapPointToNative(point.begin(), ipt, ifa - mTestFa.begin())) {
              (*infzeros)++;		// Point is outside anisotropy map

              if (*infzeros > 3) {
                success = false;
                break;
              }

              infzeros++;
              continue;
            }

            iptdwi = point.begin();
          }
          else
            iptdwi = ipt;

          // Check anisotropy at this point
          if (MRIgetVoxVal(*ifa, iptdwi[0], iptdwi[1], iptdwi[2], 0) < 0.1) {
            (*infzeros)++;		// Point is in low anisotropy area

            if (*infzeros > 3) {
              success = false;
              break;
            }
          }
          else
            *infzeros = 0;

          infzeros++;
        }

        if (!success)
          break;
      }
  }

  if (success) {
    ipt = spline.GetControlPointsBegin();
    for (vector<int>::iterator icpt = cpts.begin(); icpt < cpts.end(); icpt++) {
      *icpt = *ipt;
      ipt++;
    }
  }
  else {		// If fitting fails, I'll use equidistant points
    const int strdiv = (int) round((strlen-1) / (NumPoints-1));

    cout << "WARN: Defaulting to equidistant control points" << endl;

    ipt = Streamline.begin();
    for (vector<int>::iterator icpt = cpts.begin(); icpt < cpts.end() - 3;
                                                    icpt += 3) {
      copy(ipt, ipt+3, icpt);
      ipt += strdiv*3;
    }
    copy(Streamline.end()-3, Streamline.end(), cpts.end()-3);

    // Map control points to base space if needed
    if (!mTestAffineReg.IsEmpty() || !mTestNonlinReg.IsEmpty()) {
      vector<int> point(3);
      vector<int>::iterator icptout;

      cptsout.resize(cpts.size());
      icptout = cptsout.begin();

      for (vector<int>::const_iterator icpt = cpts.begin(); icpt < cpts.end();
                                                            icpt += 3) {
        MapPointToBase(point.begin(), icpt);
        copy(point.begin(), point.end(), icptout);

        icptout += 3;
      }
    }
  }

  cout << "INFO: Selected control points are" << endl;
  for (vector<int>::const_iterator icpt = cpts.begin(); icpt < cpts.end();
                                                        icpt += 3)
    cout << " " << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;

  cout << "INFO: Distances between consecutive points are";
  for (vector<int>::const_iterator icpt = cpts.begin() + 3; icpt < cpts.end();
                                                            icpt += 3) {
    const int dx = icpt[0] - icpt[-3],
              dy = icpt[1] - icpt[-2],
              dz = icpt[2] - icpt[-1];

    cout << " " << round(sqrt(dx*dx + dy*dy + dz*dz));
  }
  cout << endl;

  mControlPoints.push_back(cpts);

  if (!cptsout.empty()) {
    cout << "INFO: Selected control points in test subject's space are" << endl;
    for (vector<int>::const_iterator icpt = cptsout.begin();
                                     icpt < cptsout.end(); icpt += 3)
      cout << " " << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;

    cout << "INFO: Distances between consecutive points "
         << "in test subject's space are";
    for (vector<int>::const_iterator icpt = cptsout.begin() + 3;
                                     icpt < cptsout.end(); icpt += 3) {
      const int dx = icpt[0] - icpt[-3],
                dy = icpt[1] - icpt[-2],
                dz = icpt[2] - icpt[-1];

      cout << " " << round(sqrt(dx*dx + dy*dy + dz*dz));
    }
    cout << endl;

    mControlPointsTest.push_back(cptsout);

    ComputeStreamlineSpread(cptsout);
  }
  else
    ComputeStreamlineSpread(cpts);

  return success;
}

//
// Map point coordinates from atlas space to base space
//
bool Blood::MapPointToBase(vector<int>::iterator OutPoint,
                           vector<int>::const_iterator InPoint) {
  vector<float> point(InPoint, InPoint+3);

#ifndef NO_CVS_UP_IN_HERE
  if (!mTestNonlinReg.IsEmpty())
    mTestNonlinReg.ApplyXfmInv(point, point.begin());	// Inverse warp is used!
#endif
  if (!mTestAffineReg.IsEmpty())
    mTestAffineReg.ApplyXfm(point, point.begin());

  for (int k = 0; k < 3; k++)
    OutPoint[k] = (int) round(point[k]);

  return OutPoint[0] > 0 && OutPoint[0] < mTestBaseMask->width  &&
         OutPoint[1] > 0 && OutPoint[1] < mTestBaseMask->height &&
         OutPoint[2] > 0 && OutPoint[2] < mTestBaseMask->depth;
}

//
// Map point coordinates from base space to diffusion space
//
bool Blood::MapPointToNative(vector<int>::iterator OutPoint,
                             vector<int>::const_iterator InPoint,
                             unsigned int FrameIndex) {
  vector<float> point(InPoint, InPoint+3);

  mTestBaseReg[FrameIndex].ApplyXfm(point, point.begin());

  for (int k = 0; k < 3; k++)
    OutPoint[k] = (int) round(point[k]);

  return OutPoint[0] > 0 && OutPoint[0] < mTestFa[FrameIndex]->width  &&
         OutPoint[1] > 0 && OutPoint[1] < mTestFa[FrameIndex]->height &&
         OutPoint[2] > 0 && OutPoint[2] < mTestFa[FrameIndex]->depth;
}

//
// Write histograms, priors, end ROIs, and initialization control points
//
void Blood::WriteOutputs(const char *OutTrainBase, const char *OutTestBase) {
  const int nstr2 = (int) mStreamlines.size() + 2,
            nsubj2 = mNumTrain + 2;
  char fname[PATH_MAX];
  MRI *out1 = MRIclone(mTestMask[0], NULL);
  MRI *out2 = MRIclone(mTestMask[0], NULL);
  ofstream outfile;

  cout << "Writing output files to " << OutTrainBase << "_*" << endl;

  // Write streamline-wise histogram to volume
  sprintf(fname, "%s_histo_str.nii.gz", OutTrainBase);
  MRIwrite(mHistoStr, fname);

  // Convert streamline-wise histogram to negative log-likelihood
  for (int iz = 0; iz < mNz; iz++)
    for (int iy = 0; iy < mNy; iy++)
      for (int ix = 0; ix < mNx; ix++) {
        const float ratio = (MRIgetVoxVal(mHistoStr, ix, iy, iz, 0) + 1)
                            / nstr2;
        MRIsetVoxVal(out1, ix, iy, iz, 0, -log(ratio));
        MRIsetVoxVal(out2, ix, iy, iz, 0, -log(1-ratio));
      }

  sprintf(fname, "%s_logprior_str_1.nii.gz", OutTrainBase);
  MRIwrite(out1, fname);

  sprintf(fname, "%s_logprior_str_0.nii.gz", OutTrainBase);
  MRIwrite(out2, fname);

  // Write total number of streamlines to text file
  sprintf(fname, "%s_histo_nstr.txt", OutTrainBase);
  ofstream nstrfile(fname, ios::out);

  if (!nstrfile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  nstrfile << mStreamlines.size() << endl;

  MRIclear(out1);
  MRIclear(out2);

  // Write subject-wise histogram to volume
  sprintf(fname, "%s_histo.nii.gz", OutTrainBase);
  MRIwrite(mHistoSubj, fname);

  // Convert subject-wise histogram to negative log-likelihood
  for (int iz = 0; iz < mNz; iz++)
    for (int iy = 0; iy < mNy; iy++)
      for (int ix = 0; ix < mNx; ix++) {
        const float ratio = (MRIgetVoxVal(mHistoSubj, ix, iy, iz, 0) + 1)
                            / nsubj2;
        MRIsetVoxVal(out1, ix, iy, iz, 0, -log(ratio));
        MRIsetVoxVal(out2, ix, iy, iz, 0, -log(1-ratio));
      }

  sprintf(fname, "%s_logprior_1.nii.gz", OutTrainBase);
  MRIwrite(out1, fname);

  sprintf(fname, "%s_logprior_0.nii.gz", OutTrainBase);
  MRIwrite(out2, fname);

  // Save anatomical priors using only non-truncated streamlines
  WritePriors(OutTrainBase, false);

  // Save anatomical priors using all streamlines, truncated or not
  if (mUseTruncated)
    WritePriors(OutTrainBase, true);

  // Write streamline end points to volumes
  WriteEndPoints(OutTrainBase, mTestMask[0]);

  // Write central streamline to text file and volume
  sprintf(fname, "%s_cpts_all.txt", OutTrainBase);
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

  sprintf(fname, "%s_cpts_all.nii.gz", OutTrainBase);
  MRIwrite(out1, fname);

  MRIfree(&out1);
  MRIfree(&out2);

  // Write each set of control points chosen from central streamline
  // to text file and volume
  for (vector< vector<int> >::const_iterator istr = mControlPoints.begin();
                                             istr < mControlPoints.end();
                                             istr++) {
    const unsigned int ncpt = istr->size() / 3;

    sprintf(fname, "%s_cpts_%d.txt", OutTrainBase, ncpt);
    ofstream cptfile(fname, ios::out);

    if (!cptfile) {
      cout << "ERROR: Could not open " << fname << " for writing" << endl;
      exit(1);
    }

    for (vector<int>::const_iterator icpt = istr->begin();
                                     icpt < istr->end(); icpt += 3)
      cptfile << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;

    Spline spline(*istr, mTestMask[0]); 
    spline.InterpolateSpline();

    sprintf(fname, "%s_cpts_%d.nii.gz", OutTrainBase, ncpt);
    spline.WriteVolume(fname, true);
  }

  // Write spread of streamlines around each control point to text file
  for (vector< vector<float> >::const_iterator istd = mControlStd.begin();
                                               istd < mControlStd.end();
                                               istd++) {
    const unsigned int ncpt = istd->size() / 3;

    sprintf(fname, "%s_cpts_%d_std.txt", OutTrainBase, ncpt);
    ofstream stdfile(fname, ios::out);

    if (!stdfile) {
      cout << "ERROR: Could not open " << fname << " for writing" << endl;
      exit(1);
    }

    for (vector<float>::const_iterator ival = istd->begin();
                                       ival < istd->end(); ival += 3)
      stdfile << ival[0] << " " << ival[1] << " " << ival[2] << endl;
  }

  //
  // Also write initialization files in the space of the test subject,
  // if registration was specified
  //
  if (OutTestBase) {
    // Write each set of control points chosen from central streamline
    // to text file and volume
    for (vector< vector<int> >::const_iterator
         istr = mControlPointsTest.begin(); istr < mControlPointsTest.end();
                                            istr++) {
      const unsigned int ncpt = istr->size() / 3;

      sprintf(fname, "%s_cpts_%d.txt", OutTestBase, ncpt);
      ofstream cptfile(fname, ios::out);

      if (!cptfile) {
        cout << "ERROR: Could not open " << fname << " for writing" << endl;
        exit(1);
      }

      for (vector<int>::const_iterator icpt = istr->begin();
                                       icpt < istr->end(); icpt += 3)
        cptfile << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;

      Spline spline(*istr, mTestBaseMask); 
      spline.InterpolateSpline();

      sprintf(fname, "%s_cpts_%d.nii.gz", OutTestBase, ncpt);
      spline.WriteVolume(fname, true);
    }

    // Write spread of streamlines around each control point to text file
    for (vector< vector<float> >::const_iterator
         istd = mControlStdTest.begin(); istd < mControlStdTest.end(); istd++) {
      const unsigned int ncpt = istd->size() / 3;

      sprintf(fname, "%s_cpts_%d_std.txt", OutTestBase, ncpt);
      ofstream stdfile(fname, ios::out);

      if (!stdfile) {
        cout << "ERROR: Could not open " << fname << " for writing" << endl;
        exit(1);
      }

      for (vector<float>::const_iterator ival = istd->begin();
                                         ival < istd->end(); ival += 3)
        stdfile << ival[0] << " " << ival[1] << " " << ival[2] << endl;
    }
  }
}

//
// Save prior information on anatomy and curvature to text files
//
void Blood::WritePriors(const char *OutBase, bool UseTruncated) {
  char fname[PATH_MAX], pfix[5];
  vector<float>::const_iterator ithisto, itprior, ichisto, icprior;
  vector< vector<int> >::const_iterator ihisto;
  vector< vector<float> >::const_iterator iprior, idmean, idstd;
  vector< set<unsigned int> >::const_iterator iids;
  ofstream outfile;

  // Save anatomical label IDs found in training set, histograms and priors
  if (UseTruncated) {
    strcpy(pfix, "_all");
    iids = mIdsLocalAll.begin();
    ihisto = mHistoLocalAll.begin();
    iprior = mPriorLocalAll.begin();
  }
  else {
    strcpy(pfix, "");
    iids = mIdsLocal.begin();
    ihisto = mHistoLocal.begin();
    iprior = mPriorLocal.begin();
  }

  for (vector<int>::const_iterator idir = mDirLocal.begin();
                                   idir < mDirLocal.end(); idir += 3) {
    const int idx = idir[0], idy = idir[1], idz = idir[2];
    vector< set<unsigned int> >::const_iterator iidsarc = iids; 
    vector< vector<int> >::const_iterator ihistoarc = ihisto;
    vector< vector<float> >::const_iterator ipriorarc = iprior;

    sprintf(fname, "%s_fsids%s_%d_%d_%d.txt", OutBase, pfix, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (set<unsigned int>::const_iterator ival = iidsarc->begin();
                                             ival != iidsarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;

      iidsarc += mNumLocal;
    }

    outfile.close();

    sprintf(fname, "%s_fshisto%s_%d_%d_%d.txt", OutBase, pfix, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (vector<int>::const_iterator ival = ihistoarc->begin();
                                       ival != ihistoarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;

      ihistoarc += mNumLocal;
    }

    outfile.close();

    sprintf(fname, "%s_fsprior%s_%d_%d_%d.txt", OutBase, pfix, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (vector<float>::const_iterator ival = ipriorarc->begin();
                                         ival != ipriorarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;

      ipriorarc += mNumLocal;
    }

    outfile.close();

    iids++;
    ihisto++;
    iprior++;
  }

  if (UseTruncated) {
    iids = mIdsNearAll.begin();
    ihisto = mHistoNearAll.begin();
    iprior = mPriorNearAll.begin();
    idmean = mAsegDistMeanAll.begin();
    idstd = mAsegDistStdAll.begin();
  }
  else {
    iids = mIdsNear.begin();
    ihisto = mHistoNear.begin();
    iprior = mPriorNear.begin();
    idmean = mAsegDistMean.begin();
    idstd = mAsegDistStd.begin();
  }

  for (vector<int>::const_iterator idir = mDirNear.begin();
                                   idir < mDirNear.end(); idir += 3) {
    const int idx = idir[0], idy = idir[1], idz = idir[2];
    vector< set<unsigned int> >::const_iterator iidsarc = iids; 
    vector< vector<int> >::const_iterator ihistoarc = ihisto;
    vector< vector<float> >::const_iterator ipriorarc = iprior;
    vector< vector<float> >::const_iterator idmeanarc = idmean;
    vector< vector<float> >::const_iterator idstdarc = idstd;

    sprintf(fname, "%s_fsnnids%s_%d_%d_%d.txt", OutBase, pfix, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (set<unsigned int>::const_iterator ival = iidsarc->begin();
                                             ival != iidsarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;

      iidsarc += mNumNear;
    }

    outfile.close();

    sprintf(fname, "%s_fsnnhisto%s_%d_%d_%d.txt", OutBase, pfix, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (vector<int>::const_iterator ival = ihistoarc->begin();
                                       ival != ihistoarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;

      ihistoarc += mNumNear;
    }

    outfile.close();

    sprintf(fname, "%s_fsnnprior%s_%d_%d_%d.txt", OutBase, pfix, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (vector<float>::const_iterator ival = ipriorarc->begin();
                                         ival != ipriorarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;

      ipriorarc += mNumNear;
    }

    outfile.close();

    sprintf(fname, "%s_fsnndmean%s_%d_%d_%d.txt", OutBase, pfix, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (vector<float>::const_iterator ival = idmeanarc->begin();
                                         ival != idmeanarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;

      idmeanarc += mNumNear;
    }

    outfile.close();

    sprintf(fname, "%s_fsnndstd%s_%d_%d_%d.txt", OutBase, pfix, idx, idy, idz);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (vector<float>::const_iterator ival = idstdarc->begin();
                                         ival != idstdarc->end(); ival++)
        outfile << *ival << " ";
      outfile << endl;

      idstdarc += mNumNear;
    }

    outfile.close();

    iids++;
    ihisto++;
    iprior++;
    idmean++;
    idstd++;
  }

  // Save histograms and priors on tangent vector and curvature
  if (UseTruncated) {
    ihisto = mHistoTangentAll.begin();
    iprior = mPriorTangentAll.begin();
  }
  else {
    ihisto = mHistoTangent.begin();
    iprior = mPriorTangent.begin();
  }

  sprintf(fname, "%s_tanghisto%s.txt", OutBase, pfix);
  outfile.open(fname, ios::out);

  for (int k = mNumArc; k > 0; k--) {
    for (vector<int>::const_iterator ival = ihisto->begin();
                                     ival < ihisto->end(); ival++)
      outfile << *ival << " ";
    outfile << endl;

    ihisto++;
  }

  outfile.close();

  sprintf(fname, "%s_tangprior%s.txt", OutBase, pfix);
  outfile.open(fname, ios::out);

  for (int k = mNumArc; k > 0; k--) {
    for (vector<float>::const_iterator ival = iprior->begin();
                                       ival < iprior->end(); ival++)
      outfile << *ival << " ";
    outfile << endl;

    iprior++;
  }

  outfile.close();

  if (UseTruncated) {
    ihisto = mHistoCurvatureAll.begin();
    iprior = mPriorCurvatureAll.begin();
  }
  else {
    ihisto = mHistoCurvature.begin();
    iprior = mPriorCurvature.begin();
  }

  sprintf(fname, "%s_curvhisto%s.txt", OutBase, pfix);
  outfile.open(fname, ios::out);

  for (int k = mNumArc; k > 0; k--) {
    for (vector<int>::const_iterator ival = ihisto->begin();
                                     ival < ihisto->end(); ival++)
      outfile << *ival << " ";
    outfile << endl;

    ihisto++;
  }

  outfile.close();

  sprintf(fname, "%s_curvprior%s.txt", OutBase, pfix);
  outfile.open(fname, ios::out);

  for (int k = mNumArc; k > 0; k--) {
    for (vector<float>::const_iterator ival = iprior->begin();
                                       ival < iprior->end(); ival++)
      outfile << *ival << " ";
    outfile << endl;

    iprior++;
  }

  outfile.close();

  // In debug mode, save all samples of the tangent vector and curvature
  // found in the training data by arc length
  if (mDebug) {
    vector<float>::const_iterator isamp;
    vector< vector<float> >::const_iterator itangx, itangy, icurv;

    if (UseTruncated) {
      itangx = mTangentXByArcAll.begin();
      itangy = mTangentYByArcAll.begin();
      icurv = mCurvatureByArcAll.begin();
    }
    else {
      itangx = mTangentXByArc.begin();
      itangy = mTangentYByArc.begin();
      icurv = mCurvatureByArc.begin();
    }

    sprintf(fname, "%s_tangx%s.txt", OutBase, pfix);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (isamp = itangx->begin(); isamp < itangx->end(); isamp++)
        outfile << *isamp << " ";
      outfile << endl;

      itangx++;
    }

    outfile.close();

    sprintf(fname, "%s_tangy%s.txt", OutBase, pfix);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (isamp = itangy->begin(); isamp < itangy->end(); isamp++)
        outfile << *isamp << " ";
      outfile << endl;

      itangy++;
    }

    outfile.close();

    sprintf(fname, "%s_curv%s.txt", OutBase, pfix);
    outfile.open(fname, ios::out);

    for (int k = mNumArc; k > 0; k--) {
      for (isamp = icurv->begin(); isamp < icurv->end(); isamp++)
        outfile << *isamp << " ";
      outfile << endl;

      icurv++;
    }

    outfile.close();
  }
}

//
// Save central streamline to .trk file
//
void Blood::WriteCenterStreamline(const std::string CenterTrkFile,
                                  const std::string RefTrkFile) {
  float *icent, *centpts = new float[mCenterStreamline.size()];
  CTrackReader trkreader;
  CTrackWriter trkwriter;
  TRACK_HEADER trkheadin, trkheadout;

  // Open reference .trk file
  if (!trkreader.Open(RefTrkFile.c_str(), &trkheadin)) {
    cout << "ERROR: Cannot open input " << RefTrkFile << endl;
    cout << "ERROR: " << trkreader.GetLastErrorMessage() << endl;
    exit(1);
  }

  // Set output .trk header
  trkheadout = trkheadin;
  trkheadout.n_count = 1;	// Single streamline

  // Open output .trk file
  if (!trkwriter.Initialize(CenterTrkFile.c_str(), trkheadout)) {
    cout << "ERROR: Cannot open output " << CenterTrkFile << endl;
    cout << "ERROR: " << trkwriter.GetLastErrorMessage() << endl;
    exit(1);
  }

  // Make .5-based and multiply back by output voxel size
  icent = centpts;
  for (vector<int>::const_iterator ipt = mCenterStreamline.begin();
                                   ipt < mCenterStreamline.end(); ipt += 3)
    for (int k = 0; k < 3; k++) {
      *icent = (ipt[k] + .5) * trkheadout.voxel_size[k];
      icent++;
    }

  // Write center streamline to output .trk file
  trkwriter.WriteNextTrack(mCenterStreamline.size()/3, centpts);

  trkwriter.Close();
}

//
// Save streamline end points to volumes
//
void Blood::WriteEndPoints(const std::string OutBase, MRI *RefVol) {
  std::string fname;
  vector<bool>::const_iterator ivalid1 = mIsInEnd1.begin(),
                               ivalid2 = mIsInEnd2.begin();
  vector< vector<int> >::const_iterator istr;
  MRI *out1 = MRIclone(RefVol, NULL);
  MRI *out2 = MRIclone(RefVol, NULL);

  // Write end ROIs to volumes
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

  fname = OutBase + "_end1.nii.gz";
  MRIwrite(out1, fname.c_str());

  fname = OutBase + "_end2.nii.gz";
  MRIwrite(out2, fname.c_str());

  // Write dilated end ROIs to volumes
  if (!mMask.empty()) {
    vector<MRI *>::const_iterator imask = mMask.begin(),
                                  iaseg = mAseg.begin();
    vector<int> dilpt(3);

    MRIclear(out1);
    MRIclear(out2);

    istr = mStreamlines.begin();
    ivalid1 = mIsInEnd1.begin();
    ivalid2 = mIsInEnd2.begin();

    for (vector<int>::const_iterator inum = mNumLines.begin();
                                     inum != mNumLines.end(); inum++) {
      for (int k = *inum; k > 0; k--) {
        if (*ivalid1) {
          vector<int>::const_iterator iend1 = istr->begin();

          for (int iz = - mEndDilation; iz <= mEndDilation; iz++)
            for (int iy = - mEndDilation; iy <= mEndDilation; iy++)
              for (int ix = - mEndDilation; ix <= mEndDilation; ix++) {
                dilpt[0] = iend1[0] + ix;
                dilpt[1] = iend1[1] + iy;
                dilpt[2] = iend1[2] + iz;

                if (IsInCortex(dilpt.begin(), *imask, *iaseg))
                  MRIsetVoxVal(out1, dilpt[0], dilpt[1], dilpt[2], 0,
                    MRIgetVoxVal(out1, dilpt[0], dilpt[1], dilpt[2], 0) + 1);
              }
        }

        if (*ivalid2) {
          vector<int>::const_iterator iend2 = istr->end() - 3;

          for (int iz = - mEndDilation; iz <= mEndDilation; iz++)
            for (int iy = - mEndDilation; iy <= mEndDilation; iy++)
              for (int ix = - mEndDilation; ix <= mEndDilation; ix++) {
                dilpt[0] = iend2[0] + ix;
                dilpt[1] = iend2[1] + iy;
                dilpt[2] = iend2[2] + iz;

                if (IsInCortex(dilpt.begin(), *imask, *iaseg))
                  MRIsetVoxVal(out2, dilpt[0], dilpt[1], dilpt[2], 0,
                    MRIgetVoxVal(out2, dilpt[0], dilpt[1], dilpt[2], 0) + 1);
              }
        }

        istr++;
        ivalid1++;
        ivalid2++;
      }

      imask++;
      iaseg++;
    }

    fname = OutBase + "_end1_dil.nii.gz";
    MRIwrite(out1, fname.c_str());

    fname = OutBase + "_end2_dil.nii.gz";
    MRIwrite(out2, fname.c_str());
  }
}

//
// Print coordinates of all points of a streamline
//
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

//
// Average values of input volumes over all streamline voxels
//
vector<float> Blood::ComputeAvgPath(vector<MRI *> &ValueVolumes) {
  int nvox = 0;
  vector<float> avg(ValueVolumes.size(), 0);
  vector<float>::iterator iavg;

  if (mHistoStr) {
    for (int iz = 0; iz < mNz; iz++)
      for (int iy = 0; iy < mNy; iy++)
        for (int ix = 0; ix < mNx; ix++) {
          const float h = MRIgetVoxVal(mHistoStr, ix, iy, iz, 0);

          if (h > 0) {
            iavg = avg.begin();

            for (vector<MRI *>::const_iterator ivol = ValueVolumes.begin();
                                               ivol < ValueVolumes.end();
                                               ivol++) {
              *iavg += MRIgetVoxVal(*ivol, ix, iy, iz, 0);
              iavg++;
            }

            nvox++;
          }
        }

    if (nvox > 0)
      for (iavg = avg.begin(); iavg < avg.end(); iavg++)
        *iavg /= nvox;
  }

  return avg;
}

//
// Weighted average values of input volumes over all streamline voxels
//
vector<float> Blood::ComputeWeightAvgPath(vector<MRI *> &ValueVolumes) {
  float wtot = 0;
  vector<float> avg(ValueVolumes.size(), 0);
  vector<float>::iterator iavg;

  if (mHistoStr) {
    for (int iz = 0; iz < mNz; iz++)
      for (int iy = 0; iy < mNy; iy++)
        for (int ix = 0; ix < mNx; ix++) {
          const float h = MRIgetVoxVal(mHistoStr, ix, iy, iz, 0);

          if (h > 0) {
            iavg = avg.begin();

            for (vector<MRI *>::const_iterator ivol = ValueVolumes.begin();
                                               ivol < ValueVolumes.end();
                                               ivol++) {
              *iavg += h * MRIgetVoxVal(*ivol, ix, iy, iz, 0);
              iavg++;
            }

            wtot += h;
          }
        }

    if (wtot > 0)
      for (iavg = avg.begin(); iavg < avg.end(); iavg++)
        *iavg /= wtot;
  }

  return avg;
}

//
// Average values of input volumes over center streamline
//
vector<float> Blood::ComputeAvgCenter(vector<MRI *> &ValueVolumes) {
  int nvox = (int) mCenterStreamline.size()/3;
  vector<float> avg(ValueVolumes.size(), 0);
  vector<float>::iterator iavg;

  for (vector<int>::const_iterator ipt = mCenterStreamline.begin();
                                   ipt < mCenterStreamline.end(); ipt += 3) {
    iavg = avg.begin();

    for (vector<MRI *>::const_iterator ivol = ValueVolumes.begin();
                                       ivol < ValueVolumes.end(); ivol++) {
      *iavg += MRIgetVoxVal(*ivol, ipt[0], ipt[1], ipt[2], 0);
      iavg++;
    }
  }

  if (nvox > 0)
    for (iavg = avg.begin(); iavg < avg.end(); iavg++)
      *iavg /= nvox;

  return avg;
}

//
// Write values of input volumes point-wise along streamlines
//
void Blood::WriteValuesPointwise(vector<MRI *> &ValueVolumes,
                                 const std::string TextFile) {
  vector<float> valsum(ValueVolumes.size());
  vector<float>::iterator ivalsum;
  ofstream outfile(TextFile, ios::app);
  if (!outfile) {
    cout << "ERROR: Could not open " << TextFile << " for writing" << endl;
    exit(1);
  }

  cout << "Writing point-wise values along streamlines to " << TextFile << endl;

  for (vector<int>::const_iterator ipt = mCenterStreamline.begin();
                                   ipt < mCenterStreamline.end(); ipt += 3) {
    int nsamp = 0;

    // Write coordinates of this point
    outfile << ipt[0] << " " << ipt[1] << " " << ipt[2];

    // Write value of each input volume at this point
    for (vector<MRI *>::const_iterator ivol = ValueVolumes.begin();
                                       ivol < ValueVolumes.end(); ivol++)
      outfile << " " << MRIgetVoxVal(*ivol, ipt[0], ipt[1], ipt[2], 0);

    // Find closest point on each streamline
    fill(valsum.begin(), valsum.end(), 0.0);

    for (vector< vector<int> >::const_iterator ipath = mStreamlines.begin();
                                               ipath < mStreamlines.end();
                                               ipath++) {
      int dmin = 1000000;
      vector<int>::const_iterator iptmin = ipath->begin();

      for (vector<int>::const_iterator ipathpt = ipath->begin();
                                       ipathpt < ipath->end();
                                       ipathpt += 3) {
        int dist = 0;

        for (int k = 0; k < 3; k++) {
          const int diff = ipathpt[k] - ipt[k];
          dist += diff*diff;
        }

        if (dist < dmin) {
          dmin = dist;
          iptmin = ipathpt;
        }
      }

      nsamp++;

      ivalsum = valsum.begin();

      for (vector<MRI *>::const_iterator ivol = ValueVolumes.begin();
                                         ivol < ValueVolumes.end(); ivol++) {
        *ivalsum += MRIgetVoxVal(*ivol, iptmin[0], iptmin[1], iptmin[2], 0);
        ivalsum++;
      }
    }

    // Write average value from each input volume around this point
    ivalsum = valsum.begin();

    for (vector<MRI *>::const_iterator ivol = ValueVolumes.begin();
                                       ivol < ValueVolumes.end(); ivol++) {
      outfile << " " << *ivalsum / nsamp;
      ivalsum++;
    }

    outfile << endl;
  }
}

int Blood::GetVolume() { return mVolume; }

int Blood::GetNumStr() { return mNumStr; }

int Blood::GetLengthMin() { return mLengthMin; }

int Blood::GetLengthMax() { return mLengthMax; }

float Blood::GetLengthAvg() { return mLengthAvg; }

int Blood::GetNumStrEnds() { return mNumStrEnds; }

int Blood::GetLengthMinEnds() { return mLengthMinEnds; }

int Blood::GetLengthMaxEnds() { return mLengthMaxEnds; }

float Blood::GetLengthAvgEnds() { return mLengthAvgEnds; }

int Blood::GetLengthCenter() { return mCenterStreamline.size()/3; }

