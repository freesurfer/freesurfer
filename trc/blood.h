/**
 * @brief Training data and methods that probabilistic tractography feeds on
 *
 * Training data and methods that probabilistic tractography feeds on
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

#ifndef BLOOD_H
#define BLOOD_H

#include "vial.h"	// Needs to be included first because of CVS libs

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
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
          float TrainMaskLabel, const char *ExcludeFile,
          const std::vector<char *> &TestMaskList,
          const std::vector<char *> &TestFaList,
          const char *TestAffineXfmFile,
          const char *TestNonlinXfmFile,
          const char *TestNonlinRefFile,
          const std::vector<char *> &TestBaseXfmList,
          const char *TestBaseMaskFile,
          bool UseTruncated, bool UseAnatomy, bool UseShape,
          std::vector<int> &NumControls,
          int NumStrMax=INT_MAX, bool Debug=false);
    Blood(const char *TrainTrkFile,
          const char *TrainRoi1File, const char *TrainRoi2File,
          bool Debug=false);
    Blood(const char *TrainTrkFile,
          const char *RefVolFile,
          bool Debug=false);
    ~Blood();
    void SetNumControls(std::vector<int> &NumControls);
    void ReadStreamlines(const char *TrainListFile, const char *TrainTrkFile,
                         const char *TrainRoi1File, const char *TrainRoi2File,
                         float TrainMaskLabel, const char *ExcludeFile);
    void ReadAnatomy(const char *TrainListFile, const char *TrainAsegFile,
                                                const char *TrainMaskFile);
    void RemoveLengthOutliers();
    void PrepStreamlines();
    void MatchStreamlineEnds();
    void SubsampleStreamlines();
    void MatchTruncatedStreamlines();
    void ComputeHistogram();
    void ComputePriors();
    void FindCenterStreamline(bool CheckOverlap=true, bool CheckDeviation=true,
                                                      bool CheckFa=true);
    void WriteOutputs(const char *OutTrainBase, const char *OutTestBase=NULL);
    void WriteCenterStreamline(const char *CenterTrkFile,
                               const char *RefTrkFile);
    void WriteEndPoints(const char *OutBase, MRI *RefVol);
    void WriteStreamlines(const char *TrainListFile, const char *TrainTrkFile);
    void PrintStreamline(int SubjIndex, int LineIndex);
    std::vector<float> ComputeAvgPath(std::vector<MRI *> &ValueVolumes);
    std::vector<float> ComputeWeightAvgPath(std::vector<MRI *> &ValueVolumes);
    std::vector<float> ComputeAvgCenter(std::vector<MRI *> &ValueVolumes);
    void WriteValuesPointwise(std::vector<MRI *> &ValueVolumes,
                              const char *TextFile);
    int GetVolume();
    int GetNumStr();
    int GetLengthMin();
    int GetLengthMax();
    float GetLengthAvg();
    int GetNumStrEnds();
    int GetLengthMinEnds();
    int GetLengthMaxEnds();
    float GetLengthAvgEnds();
    int GetLengthCenter();

  private:
    static const int mDistThresh, mEndDilation;
    static const unsigned int mDiffStep;
    static const float mLengthCutoff, mLengthRatio,
                       mHausStepRatio, mControlStepRatio,
                       mTangentBinSize, mCurvatureBinSize;

    const bool mDebug, mUseTruncated, mUseAnatomy, mUseShape;
    const int mNumStrMax;
    bool mHavePresorted;
    int mNx, mNy, mNz, mNumTrain, mVolume,
        mNumStr, mLengthMin, mLengthMax,
        mNumStrEnds, mLengthMinEnds, mLengthMaxEnds,
        mNumArc, mNumLocal, mNumNear;
    float mMaskLabel, mDx, mLengthAvg, mLengthAvgEnds;
    std::vector<bool> mIsInEnd1, mIsInEnd2, mIsOutHist, mIsOutDev, mIsOutFa;
    std::vector<int> mNumLines, mLengths, mMidPoints, mTruncatedLengths,
                     mDistanceRank,
                     mCenterStreamline, mDirLocal, mDirNear, mNumControls;
    std::vector<float> mMeanEnd1, mMeanEnd2, mMeanMid,
                       mVarEnd1, mVarEnd2, mVarMid;
    std::vector< std::vector<int> > mStreamlines,
                                    mControlPoints, mControlPointsTest,
                                    mExcludedStreamlines,
                                    mHistoLocal, mHistoNear,
                                    mHistoLocalAll, mHistoNearAll,
                                    mHistoTangent, mHistoCurvature,
                                    mHistoTangentAll, mHistoCurvatureAll;
    std::vector< std::vector<float> > mControlStd, mControlStdTest,
                                      mPriorLocal, mPriorNear,
                                      mPriorLocalAll, mPriorNearAll,
                                      mAsegDistMean, mAsegDistStd,
                                      mAsegDistMeanAll, mAsegDistStdAll,
                                      mPriorTangent, mPriorCurvature,
                                      mPriorTangentAll, mPriorCurvatureAll,
                                      mTangentXByArc, mTangentYByArc,
                                      mCurvatureByArc,
                                      mTangentXByArcAll, mTangentYByArcAll,
                                      mCurvatureByArcAll;
    std::vector< std::set<unsigned int> > mIdsLocal, mIdsNear, //[{6,7}xmNumArc]
                                          mIdsLocalAll, mIdsNearAll;
    std::vector<MRI *> mRoi1, mRoi2, mAseg, mMask, mTestMask, mTestFa;
    AffineReg mTestAffineReg;
#ifndef NO_CVS_UP_IN_HERE
    NonlinReg mTestNonlinReg;
#endif
    std::vector<AffineReg> mTestBaseReg;
    MRI *mHistoStr, *mHistoSubj, *mTestBaseMask;

    void AllocateHistogram(MRI *RefVol);
    void ReadExcludedStreamlines(const char *ExcludeFile);
    void ComputeStats();
    void ComputeStatsEnds();
    void ComputeEndPointCoM();
    bool IsInMask(std::vector<int>::const_iterator Point);
    bool IsInMask(float *Point);
    bool IsInRoi(std::vector<int>::const_iterator Point, MRI *Roi);
    bool IsInCortex(std::vector<int>::const_iterator Point,
                    MRI *Mask, MRI *Aseg);
    void FindOutlierStreamlines(bool CheckOverlap, bool CheckDeviation,
                                                   bool CheckFa);
    void RankStreamlineDistance();
    void UpdateDistanceRank();
    void SetArcSegments();
    void ComputeAnatomyPrior();
    void ComputeShapePrior();
    void FindPointsOnStreamline(std::vector<int> &Streamline, int NumPoints);
    bool FindPointsOnStreamlineLS(std::vector<int> &Streamline, int NumPoints);
    bool FindPointsOnStreamlineComb(std::vector<int> &Streamline,
                                    int NumPoints);
    void TryControlPoint(double &HausDistMin,
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
    bool MapPointToBase(std::vector<int>::iterator OutPoint,
                        std::vector<int>::const_iterator InPoint);
    bool MapPointToNative(std::vector<int>::iterator OutPoint,
                          std::vector<int>::const_iterator InPoint,
                          unsigned int FrameIndex);
    void WriteAnatomyPriors(const char *OutBase);
    void WriteShapePriors(const char *OutBase);
};

#endif

