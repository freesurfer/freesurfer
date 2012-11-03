/**
 * @file  coffin.h
 * @brief Main container of tractography data and methods
 *
 * Main container of tractography data and methods
 */
/*
 * Original Author: Anastasia Yendiki
 * CVS Revision Info:
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

#ifndef COFFIN_H
#define COFFIN_H

#include "vial.h"	// Needs to be included first because of CVS libs

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <math.h>
#include <limits.h>
#include "utils.h"
#include "pdf.h"
#include "mri.h"
#include "bite.h"
#include "spline.h"

class Coffin {
  public:
    Coffin(const char *OutDir, const std::vector<char *> InDir,
           const char *DwiFile,
           const char *GradientFile, const char *BvalueFile,
           const char *MaskFile, const char *BedpostDir,
           const int NumTract, const float FminPath,
           const char *InitFile,
           const char *RoiFile1, const char *RoiFile2,
           const char *RoiMeshFile1, const char *RoiMeshFile2,
           const char *RoiRefFile1, const char *RoiRefFile2,
           const char *XyzPriorFile0, const char *XyzPriorFile1,
           const char *TangPriorFile, const char *CurvPriorFile,
           const char *NeighPriorFile, const char *NeighIdFile,
           const int NeighPriorSet,
           const char *LocalPriorFile, const char *LocalIdFile,
           const int LocalPriorSet,
           const char *AsegFile,
           const char *AffineXfmFile, const char *NonlinXfmFile,
           const int NumBurnIn, const int NumSample,
           const int KeepSampleNth, const int UpdatePropNth,
           const char *PropStdFile,
           const bool Debug=false);
    ~Coffin();
    void SetOutputDir(const char *OutDir);
    void SetPathway(const char *InitFile,
                    const char *RoiFile1, const char *RoiFile2,
                    const char *RoiMeshFile1, const char *RoiMeshFile2,
                    const char *RoiRefFile1, const char *RoiRefFile2,
                    const char *XyzPriorFile0, const char *XyzPriorFile1,
                    const char *TangPriorFile, const char *CurvPriorFile,
                    const char *NeighPriorFile, const char *NeighIdFile,
                    const char *LocalPriorFile, const char *LocalIdFile);
    void SetMCMCParameters(const int NumBurnIn, const int NumSample,
                           const int KeepSampleNth, const int UpdatePropNth,
                           const char *PropStdFile);
    bool RunMCMC();
    bool RunMCMC1();
    void WriteOutputs();

  private:
    static const unsigned int mDiffStep;
    static const float mTangentBinSize, mCurvatureBinSize;
    bool mRejectSpline, mRejectPosterior,
         mRejectF, mAcceptF, mRejectTheta, mAcceptTheta;
    const bool mDebug;
    int mNx, mNy, mNz, mNxy, mNumVox, mNumControl,
        mNxAtlas, mNyAtlas, mNzAtlas, mNumArc,
        mPriorSetLocal, mPriorSetNear,
        mNumBurnIn, mNumSample, mKeepSampleNth, mUpdatePropNth;
    double mLikelihoodOnPath, mPriorOnPath, mPosteriorOnPath,
           mLikelihoodOnPathNew, mPriorOnPathNew, mPosteriorOnPathNew,
           mPosteriorOnPathMap,
           mLikelihoodOffPath, mPriorOffPath, mPosteriorOffPath,
           mLikelihoodOffPathNew, mPriorOffPathNew, mPosteriorOffPathNew,
           mXyzPriorOnPath, mXyzPriorOnPathNew,
           mXyzPriorOffPath, mXyzPriorOffPathNew,
           mAnatomicalPrior, mAnatomicalPriorNew,
           mShapePrior, mShapePriorNew;
    std::string mInfoGeneral, mInfoPathway, mInfoMCMC;
    std::vector<bool> mRejectControl;			// [mNumControl]
    std::vector<int> mAcceptCount, mRejectCount,	// [mNumControl]
                     mControlPoints, mControlPointsNew, mControlPointsMap,
                     mPathPoints, mPathPointsNew, mPathPointsMap,
                     mControlPointSamples, mPathLengthSamples,
                     mDirLocal, mDirNear;
    std::vector<float> mProposalStdInit, mProposalStd,	// [mNumControl x 3]
                       mControlPointJumps,		// [mNumControl x 3]
                       mAcceptSpan, mRejectSpan,	// [mNumControl x 3]
                       mPathPhi, mPathPhiNew,
                       mPathTheta, mPathThetaNew,
                       mPathTangent, mPathTangentNew,
                       mPathCurvature, mPathCurvatureNew,
                       mPosteriorSamples;
    std::vector<Bite> mData;				// [mNumVox]
    std::vector<Bite *>mDataMask;			// [mNx x mNy x mNz]
    std::vector< std::vector<int> > mPathPointSamples;
    std::vector< std::vector<unsigned int> > mIdsLocal, mIdsNear;
    std::vector< std::vector<float> > mPriorTangent,	// [mNumArc]
                                      mPriorCurvature,	// [mNumArc]
                                      mPriorLocal, 	// [mNumArc x 7]
                                      mPriorNear;	// [mNumArc x 6]
    char *mOutDir;
    MRI *mMask, *mRoi1, *mRoi2, 
        *mXyzPrior0, *mXyzPrior1,
        *mAseg,
        *mPathSamples;
    std::ofstream mLog;
    Spline mSpline;
    AffineReg mAffineReg;
#ifndef NO_CVS_UP_IN_HERE
    NonlinReg mNonlinReg;
#endif

    void ReadControlPoints(const char *ControlPointFile);
    void ReadProposalStds(const char *PropStdFile);
    bool InitializeMCMC();
    bool JumpMCMC();
    bool JumpMCMC1(int ControlIndex);
    bool ProposePath();
    bool ProposePath1(int ControlIndex);
    void ProposeDiffusionParameters();
    bool AcceptPath();
    double ComputeXyzPriorOffPath(std::vector<int> &PathAtlasPoints);
    double ComputeXyzPriorOnPath(std::vector<int> &PathAtlasPoints);
    double ComputeAnatomicalPrior(std::vector<int> &PathAtlasPoints);
    double ComputeShapePrior(std::vector<int> &PathAtlasPoints);
    void UpdatePath();
    void UpdateAcceptanceRate();
    void UpdateRejectionRate();
    void UpdateAcceptRejectRate1();
    void UpdateProposalStd();
    void SavePathPosterior(bool IsPathAccepted);
    void SavePath();
    bool IsInMask(std::vector<int>::const_iterator Point);
    bool IsInRoi(std::vector<int>::const_iterator Point, MRI *Roi);
    bool IsZigZag(std::vector<int> &ControlPoints,
                  std::vector<int>::const_iterator FirstPerturbedPoint,
                  std::vector<int>::const_iterator LastPerturbedPoint);
    void MapPointToAtlas(std::vector<int>::iterator OutPoint,
                         std::vector<int>::const_iterator InPoint);
    void LogObjective();
    void LogObjectiveNaN();
    void FindPathMAP();
};

#endif

