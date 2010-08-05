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

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "utils.h"
#include "pdf.h"
#include "mri.h"
#include "bite.h"
#include "spline.h"

class Coffin {
  public:
    Coffin(const char *OutDir, const char *DwiFile,
           const char *GradientFile, const char *BvalueFile,
           const char *MaskFile,
           const char *BedpostDir, const unsigned int NumTract,
           const char *RoiFile1, const char *RoiFile2,
           const char *RoiMeshFile1, const char *RoiMeshFile2,
           const char *RoiRefFile1, const char *RoiRefFile2,
           const char *XfmFile, const char *InitFile,
           const char *PriorFile0, const char *PriorFile1,
           const unsigned int AsegPriorType,
           const char *AsegPriorFile0, const char *AsegPriorFile1,
           const char *AsegIdFile, 
           const char *AsegTrainFile, const char *PathTrainFile,
           const char *AsegFile,
           const unsigned int NumBurnIn, const unsigned int NumSample,
           const unsigned int KeepSampleNth, const unsigned int UpdatePropNth,
           const char *PropStdFile,
           const bool Debug=false);
    ~Coffin();
    void RunMCMC();
    void RunMCMC1();
    void WriteOutputs();

  private:
    bool mRejectSpline, mRejectPosterior,
         mRejectF, mAcceptF, mRejectTheta, mAcceptTheta;
    const bool mDebug;
    unsigned int mNx, mNy, mNz, mNxy, mNumVox, mNumControl;
    const unsigned int mNumBurnIn, mNumSample, mKeepSampleNth, mUpdatePropNth;
    double mLikelihoodOnPath, mPriorOnPath, mPosteriorOnPath,
           mLikelihoodOnPathNew, mPriorOnPathNew, mPosteriorOnPathNew,
           mPosteriorOnPathMap,
           mLikelihoodOffPath, mPriorOffPath, mPosteriorOffPath,
           mLikelihoodOffPathNew, mPriorOffPathNew, mPosteriorOffPathNew;
    std::vector<bool> mRejectControl;				// [mNumControl]
    std::vector<unsigned int> mAcceptCount, mRejectCount;	// [mNumControl]
    std::vector<int> mControlPoints, mControlPointsNew, mControlPointsMap,
                     mPathPoints, mPathPointsNew, mPathPointsMap,
                     mControlPointSamples, mPathLengthSamples;
    std::vector<float> mDwiToRoi;			// [4 x 4]
    std::vector<float> mDwiVoxelSize, mRoiVoxelSize;	// [3]
    std::vector<float> mProposalStdInit, mProposalStd, mControlPointJumps,
                       mAcceptSpan, mRejectSpan;	// [mNumControl x 3]
    std::vector<float> mLikelihoodOnPathSamples, mPriorOnPathSamples,
                       mPosteriorOnPathSamples;
    std::vector<float> mPathPhi, mPathPhiNew, mPathTheta, mPathThetaNew;
    std::vector<Bite> mData;				// [mNumVox]
    std::vector<Bite *>mDataMask;			// [mNx x mNy x mNz]
    char *mOutDir;
    MRI *mMask, *mRoi1, *mRoi2, *mPathSamples;
    Spline mSpline;

    void ReadDwiToRoi(const char *XfmFile);
    void ReadControlPoints(const char *ControlPointFile);
    void ReadProposalStds(const char *PropStdFile);
    void InitializeMCMC();
    bool JumpMCMC();
    bool JumpMCMC1(unsigned int ControlIndex);
    bool ProposePath();
    bool ProposePath1(unsigned int ControlIndex);
    void ProposeDiffusionParameters();
    bool AcceptPath();
    void UpdatePath();
    void UpdateAcceptanceRate();
    void UpdateRejectionRate();
    void UpdateAcceptRejectRate1();
    void UpdateProposalStd();
    void SavePath();
    bool IsInMask(std::vector<int>::const_iterator Point);
    bool IsInRoi(std::vector<int>::const_iterator Point, MRI *Roi);
    void ApplyAffine(std::vector<int> &OutPoint,
                     std::vector<int>::const_iterator InPoint,
                     std::vector<float> &In2Out,
                     std::vector<float> &OutVoxelSize,
                     std::vector<float> &InVoxelSize);
};

#endif

