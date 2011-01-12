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

#ifndef COFFIN_H
#define COFFIN_H

#include "vial.h"	// Needs to be included first because of CVS libs

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
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
    Coffin(const char *OutDir, const char *DwiFile,
           const char *GradientFile, const char *BvalueFile,
           const char *MaskFile, const char *BedpostDir,
           const int NumTract, const float FminPath,
           const char *InitFile,
           const char *RoiFile1, const char *RoiFile2,
           const char *RoiMeshFile1, const char *RoiMeshFile2,
           const char *RoiRefFile1, const char *RoiRefFile2,
           const char *PriorFile0, const char *PriorFile1,
           const int AsegPriorType,
           const char *AsegPriorFile0, const char *AsegPriorFile1,
           const char *AsegIdFile, 
           const char *AsegTrainFile, const char *PathTrainFile,
           const char *NeighPriorFile, const char *NeighIdFile,
           const char *LocalPriorFile, const char *LocalIdFile,
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
                    const char *PriorFile0, const char *PriorFile1,
                    int AsegPriorType,
                    const char *AsegPriorFile0, const char *AsegPriorFile1,
                    const char *AsegIdFile,
                    const char *AsegTrainFile, const char *PathTrainFile,
                    const char *NeighPriorFile, const char *NeighIdFile,
                    const char *LocalPriorFile, const char *LocalIdFile);
    void SetMCMCParameters(const int NumBurnIn, const int NumSample,
                           const int KeepSampleNth, const int UpdatePropNth,
                           const char *PropStdFile);
    bool RunMCMC();
    bool RunMCMC1();
    void WriteOutputs();

  private:
    bool mRejectSpline, mRejectPosterior,
         mRejectF, mAcceptF, mRejectTheta, mAcceptTheta;
    const bool mDebug;
    int mNx, mNy, mNz, mNxy, mNumVox, mNumControl,
        mNxAtlas, mNyAtlas, mNzAtlas, mNumArc, mNumLocal, mNumNear,
        mNumBurnIn, mNumSample, mKeepSampleNth, mUpdatePropNth;
    double mLikelihoodOnPath, mPriorOnPath, mPosteriorOnPath,
           mLikelihoodOnPathNew, mPriorOnPathNew, mPosteriorOnPathNew,
           mPosteriorOnPathMap,
           mLikelihoodOffPath, mPriorOffPath, mPosteriorOffPath,
           mLikelihoodOffPathNew, mPriorOffPathNew, mPosteriorOffPathNew,
           mAnatomicalPrior, mAnatomicalPriorNew;
    std::string mInfoGeneral, mInfoPathway, mInfoMCMC;
    std::vector<bool> mRejectControl;				// [mNumControl]
    std::vector<int> mAcceptCount, mRejectCount;		// [mNumControl]
    std::vector<int> mControlPoints, mControlPointsNew, mControlPointsMap,
                     mPathPoints, mPathPointsNew, mPathPointsMap,
                     mControlPointSamples, mPathLengthSamples,
                     mDirLocal, mDirNear;
    std::vector<float> mProposalStdInit, mProposalStd, mControlPointJumps,
                       mAcceptSpan, mRejectSpan;	// [mNumControl x 3]
    std::vector<float> mLikelihoodOnPathSamples, mPriorOnPathSamples,
                       mPosteriorOnPathSamples;
    std::vector<float> mPathPhi, mPathPhiNew, mPathTheta, mPathThetaNew;
    std::vector<Bite> mData;				// [mNumVox]
    std::vector<Bite *>mDataMask;			// [mNx x mNy x mNz]
    std::vector< std::vector<unsigned int> > mIdsLocal, mIdsNear;
    std::vector< std::vector<float> > mPriorLocal, mPriorNear; //[mNumArcx{6,7}]
    char *mOutDir;
    MRI *mMask, *mRoi1, *mRoi2, 
        *mPathPrior0, *mPathPrior1, 
        *mAsegTrain, *mPathTrain, *mAseg,
        *mPathSamples;
    std::vector<MRI *> mSegPrior0, mSegPrior1;
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
    double ComputeAnatomicalPrior(std::vector<int> &PathPoints);
    void UpdatePath();
    void UpdateAcceptanceRate();
    void UpdateRejectionRate();
    void UpdateAcceptRejectRate1();
    void UpdateProposalStd();
    void SavePath();
    bool IsInMask(std::vector<int>::const_iterator Point);
    bool IsInRoi(std::vector<int>::const_iterator Point, MRI *Roi);
    void MapToAtlas(std::vector<int> &OutPoint,
                    std::vector<int>::const_iterator InPoint);
    void LogInputs();
};

#endif

