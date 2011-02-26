/**
 * @file  blood.h
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
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef BLOOD_H
#define BLOOD_H

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <limits.h>
#include <math.h>
#include "mri.h"

#include "TrackIO.h"
#include "spline.h"

class Blood {
  public:
    Blood(const char *TrainListFile, const char *TrainTrkFile,
          const char *TrainRoi1File, const char *TrainRoi2File,
          const char *TrainAsegFile, const char *TrainMaskFile,
          float TrainMaskLabel, const char *TestMaskFile,
          bool UseTruncated);
    Blood(const char *TrainTrkFile);
    ~Blood();
    void ReadStreamlines(const char *TrainListFile, const char *TrainTrkFile,
                         const char *TrainRoi1File, const char *TrainRoi2File,
                         float TrainMaskLabel);
    void ReadStreamlines(const char *TrainTrkFile);
    void ReadAnatomy(const char *TrainListFile, const char *TrainAsegFile,
                                                const char *TrainMaskFile);
    void FindCenterStreamline();
    void ComputePriors();
    void SelectControlPoints(int NumControls);
    void WriteOutputs(const char *OutBase);
    void WriteCenterStreamline(const char *CenterTrkFile,
                               const char *RefTrkFile);
    void PrintStreamline(int SubjIndex, int LineIndex);
    int GetNumStr();
    int GetLengthMin();
    int GetLengthMax();
    float GetLengthAvg();
    int GetNumStrEnds();
    int GetLengthMinEnds();
    int GetLengthMaxEnds();
    float GetLengthAvgEnds();

  private:
    static const int mDistThresh;
    static const unsigned int mCurvOffset;
    static const float mLengthCutoff, mLengthRatio,
                       mHausStepRatio, mControlStepRatio;

    const bool mUseTruncated;
    int mNx, mNy, mNz, mNumTrain,
        mNumStr, mLengthMin, mLengthMax,
        mNumStrEnds, mLengthMinEnds, mLengthMaxEnds,
        mNumArc, mNumLocal, mNumNear;
    float mMaskLabel, mDx, mLengthAvg, mLengthAvgEnds;
    std::vector<bool> mIsInEnd1, mIsInEnd2;
    std::vector<int> mNumLines, mLengths, mTruncatedLengths, mCenterStreamline,
                     mDirLocal, mDirNear;
    std::vector<float> mTangentMean, mTangentStd,
                       mTangentMeanAll, mTangentStdAll,
                       mCurvatureMean, mCurvatureStd,
                       mCurvatureMeanAll, mCurvatureStdAll;
    std::vector< std::vector<int> > mStreamlines, mControlPoints,
                                    mHistoLocal, mHistoNear,
                                    mHistoLocalAll, mHistoNearAll;
    std::vector< std::vector<float> > mControlStd,
                                      mPriorLocal, mPriorNear,
                                      mPriorLocalAll, mPriorNearAll,
                                      mAsegDistMean, mAsegDistStd,
                                      mAsegDistMeanAll, mAsegDistStdAll;
    std::vector< std::set<unsigned int> > mIdsLocal, mIdsNear, //[{6,7}xmNumArc]
                                          mIdsLocalAll, mIdsNearAll;
    std::vector<MRI *> mRoi1, mRoi2, mAseg, mMask;
    MRI *mTestMask, *mHistoStr, *mHistoSubj;

    void ComputeStats();
    void ComputeStatsEnds();
    bool IsInMask(std::vector<int>::const_iterator Point);
    bool IsInMask(float *Point);
    bool IsInRoi(std::vector<int>::const_iterator Point, MRI *Roi);
    void MatchStreamlineEnds();
    void SetArcSegments();
    void ComputeHistogram();
    void ComputeAnatomyPrior(bool UseTruncated);
    void ComputeCurvaturePrior(bool UseTruncated);
    void FindPointsOnStreamline(std::vector<int> &Streamline, int NumPoints);
    void FindPointsOnStreamlineComb(std::vector<int> &Streamline,
                                    int NumPoints);
    void TryControlPoint(float &OverlapMax,
                         int IndexPoint,
                         int SearchLag,
                         std::vector<int> &ControlPointsMax,
                         std::vector<std::vector<int>::const_iterator>
                           &ControlPoints,
                         Spline &TrySpline,
                         std::vector<int> &Streamline);
    void ComputeStreamlineSpread(std::vector<int> &ControlPoints);
    void FlipStreamline(std::vector< std::vector<int> >::iterator Streamline);
    bool IsEnd1InMask(std::vector< std::vector<int> >::iterator Streamline,
                      MRI *Mask, MRI *Aseg);
    bool IsEnd2InMask(std::vector< std::vector<int> >::iterator Streamline,
                      MRI *Mask, MRI *Aseg);
    void WritePriors(const char *OutBase, bool UseTruncated);
};

#endif

