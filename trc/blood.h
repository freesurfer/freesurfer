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
          const char *TrainAsegFile, const char *TrainMaskFile,
          const char *TestMaskFile);
    Blood(const char *TrainTrkFile);
    ~Blood();
    void ReadStreamlines(const char *TrainListFile, const char *TrainTrkFile);
    void ReadStreamlines(const char *TrainTrkFile);
    void ReadAnatomy(const char *TrainListFile, const char *TrainAsegFile,
                                                const char *TrainMaskFile);
    void ComputePriors();
    void SelectControlPoints(int NumControls);
    void WriteOutputs(const char *OutBase);
    void PrintStreamline(int SubjIndex, int LineIndex);

  private:
    static const int mDistThresh;
    static const unsigned int mCurvOffset;
    static const float mLengthCutoff, mLengthRatio;

    int mNx, mNy, mNz, mNumTrain, mNumMask, mNumArc, mNumLocal, mNumNear;
    std::vector<bool> mIsInEnd1, mIsInEnd2;
    std::vector<int> mNumLines, mLengths, mCenterStreamline,
                     mDirLocal, mDirNear;
    std::vector<float> mTangentMean, mTangentStd,
                       mCurvatureMean, mCurvatureStd;
    std::vector< std::vector<int> > mStreamlines, mControlPoints,
                                    mHistoLocal, mHistoNear,
                                    mAsegDist, mAsegDistArc;
    std::vector< std::vector<unsigned int> > mAsegLocal, mAsegNear,
                                             mAsegLocalArc, mAsegNearArc;
    std::vector< std::vector<float> > mControlStd,
                                      mPriorLocal, mPriorNear,
                                      mAsegDistMean, mAsegDistStd,
                                      mTangent, mTangentArc,
                                      mCurvature, mCurvatureArc;
    std::vector< std::set<unsigned int> > mIdsLocal, mIdsNear; //[{6,7}xmNumArc]
    std::vector<MRI *> mAseg, mMask;
    MRI *mTestMask;

    bool IsInMask(std::vector<int>::const_iterator Point);
    void MatchStreamlineEnds();
    void SetArcSegments();
    void ComputeAnatomyPrior();
    void ComputeCurvaturePrior();
    void FindCenterStreamline();
    void FindPointsOnStreamline(std::vector<int> &Streamline, int NumPoints);
    void FlipStreamline(std::vector< std::vector<int> >::iterator Streamline);
    bool IsValidEnd1(std::vector< std::vector<int> >::iterator Streamline,
                     MRI *Mask);
    bool IsValidEnd2(std::vector< std::vector<int> >::iterator Streamline,
                     MRI *Mask);
};

#endif

