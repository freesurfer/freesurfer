/**
 * @brief Main container of tractography data and methods
 *
 * Main container of tractography data and methods
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

#include "TrackIO.h"

class Aeon {		// One point in time
  public:
    Aeon();
    ~Aeon();
    static void SetBaseMask(MRI *BaseMask);
    static void SavePathPriors(std::vector<float> &Priors);
    static void SaveBasePath(std::vector<int> &PathPoints);
    static void SetPathMap(unsigned int PathIndex);
    void ReadData(const std::string RootDir, const std::string DwiFile,
                  const std::string GradientFile, const std::string BvalueFile,
                  const std::string MaskFile, const std::string BedpostDir,
                  const int NumTract, const float FminPath,
                  const std::string BaseXfmFile);
    MRI *GetMask() const;
    MRI *GetBaseMask() const;
    float GetDx() const;
    float GetDy() const;
    float GetDz() const;
    unsigned int GetNumSample() const;
    void FreeMask();
    void SetOutputDir(const std::string OutDir);
    const string &GetOutputDir() const;
    void ClearPath();
    bool MapPathFromBase(Spline &BaseSpline);
    void FindDuplicatePathPoints(std::vector<bool> &IsDuplicate);
    void RemovePathPoints(std::vector<bool> &DoRemove, unsigned int NewSize=0);
    void ProposeDiffusionParameters();
    bool ComputePathDataFit();
    int FindErrorSegment(Spline &BaseSpline);
    void UpdatePath();
    void SavePathDataFit(bool IsPathAccepted);
    void SavePath();
    void WriteOutputs();
    unsigned int GetNumFZerosNew() const;
    unsigned int GetNumFZeros() const;
    bool RejectF() const;
    bool AcceptF() const;
    bool RejectTheta() const;
    bool AcceptTheta() const;
    const string &GetLog() const;
    unsigned int GetPathLengthNew() const;
    unsigned int GetPathLength() const;
    double GetLikelihoodOffPathNew() const;
    double GetLikelihoodOffPath() const;
    double GetLikelihoodOnPathNew() const;
    double GetLikelihoodOnPath() const;
    double GetPriorOffPathNew() const;
    double GetPriorOffPath() const;
    double GetPriorOnPathNew() const;
    double GetPriorOnPath() const;
    double GetPosteriorOffPathNew() const;
    double GetPosteriorOffPath() const;
    double GetPosteriorOnPathNew() const;
    double GetPosteriorOnPath() const;
    double GetDataFitNew() const;
    double GetDataFit() const;

  private:
    static const unsigned int mDiffStep;
    static int mMaxAPosterioriPath;
    static unsigned int mMaxAPosterioriPath0;
    static std::vector<float> mPriorSamples;
    static std::vector< std::vector<int> > mBasePathPointSamples;
    static MRI *mBaseMask;

    bool mRejectF, mAcceptF, mRejectTheta, mAcceptTheta;
    int mNx, mNy, mNz, mNxy, mNumVox;
    unsigned int mPathLength, mPathLengthNew;
    double mLikelihoodOnPath, mPriorOnPath, mPosteriorOnPath,
           mLikelihoodOnPathNew, mPriorOnPathNew, mPosteriorOnPathNew,
           mLikelihoodOffPath, mPriorOffPath, mPosteriorOffPath,
           mLikelihoodOffPathNew, mPriorOffPathNew, mPosteriorOffPathNew;
    MRI *mMask;
    string mRootDir, mOutDir, mLog;
    std::vector<int> mPathPoints, mPathPointsNew, mErrorPoint;
    std::vector<float> mPathPhi, mPathPhiNew,
                       mPathTheta, mPathThetaNew,
                       mDataFitSamples;
    std::vector< std::vector<int> > mPathPointSamples;
    std::vector<Bite> mData;				// [mNumVox]
    std::vector<Bite *>mDataMask;			// [mNx x mNy x mNz]
    AffineReg mBaseReg;

    bool IsInMask(std::vector<int>::const_iterator Point);
    void ComputePathLengths(std::vector<int> &PathLengths,
                            std::vector< std::vector<int> > &PathSamples);
    void ComputePathHisto(MRI *HistoVol,
                          std::vector< std::vector<int> > &PathSamples);
    int FindMaxAPosterioriPath(std::vector< std::vector<int> > &PathSamples,
                               std::vector<int> &PathLengths, MRI *PathHisto);
};

class Coffin {		// The main container
 public:
  Coffin(const std::string OutDir, std::vector<std::string> InDirList,
	 const std::string DwiFile,
	 const std::string GradientFile, const std::string BvalueFile,
	 const std::string MaskFile, const std::string BedpostDir,
	 const int NumTract, const float FminPath,
	 const std::string BaseXfmFile, const std::string BaseMaskFile,
	 const std::string InitFile,
	 const std::string RoiFile1, const std::string RoiFile2,
	 const std::string RoiMeshFile1, const std::string RoiMeshFile2,
	 const std::string RoiRefFile1, const std::string RoiRefFile2,
	 const std::string XyzPriorFile0, const std::string XyzPriorFile1,
	 const std::string angPriorFile, const std::string CurvPriorFile,
	 const std::string NeighPriorFile, const std::string NeighIdFile,
	 const int NeighPriorSet,
	 const std::string LocalPriorFile, const std::string LocalIdFile,
	 const int LocalPriorSet,
	 const std::vector<std::string> AsegList,
	 const std::string AffineXfmFile, const std::string NonlinXfmFile,
	 const int NumBurnIn, const int NumSample,
	 const int KeepSampleNth, const int UpdatePropNth,
	 const std::string PropStdFile,
	 const bool Debug=false);
    ~Coffin();
    void SetOutputDir(const std::string OutDir);
    void SetPathway(const std::string InitFile,
                    const std::string RoiFile1, const std::string RoiFile2,
                    const std::string RoiMeshFile1, const std::string RoiMeshFile2,
                    const std::string RoiRefFile1, const std::string RoiRefFile2,
                    const std::string XyzPriorFile0, const std::string XyzPriorFile1,
                    const std::string TangPriorFile, const std::string CurvPriorFile,
                    const std::string NeighPriorFile, const std::string NeighIdFile,
                    const std::string LocalPriorFile, const std::string LocalIdFile);
    void SetMcmcParameters(const int NumBurnIn, const int NumSample,
                           const int KeepSampleNth, const int UpdatePropNth,
                           const std::string PropStdFile);
    bool RunMcmcFull();
    bool RunMcmcSingle();
    void WriteOutputs();

  private:
    static const unsigned int mMaxTryMask, mMaxTryWhite, mDiffStep;
    static const float mTangentBinSize, mCurvatureBinSize;
    bool mRejectSpline, mRejectPosterior,
         mRejectF, mAcceptF, mRejectTheta, mAcceptTheta;
    const bool mDebug;
    int mNx, mNy, mNz, mNxy, mNumControl,
        mNxAtlas, mNyAtlas, mNzAtlas, mNumArc,
        mPriorSetLocal, mPriorSetNear,
        mNumBurnIn, mNumSample, mKeepSampleNth, mUpdatePropNth;
    double mDataPosteriorOnPath, mDataPosteriorOnPathNew,
           mDataPosteriorOffPath, mDataPosteriorOffPathNew,
           mXyzPriorOnPath, mXyzPriorOnPathNew,
           mXyzPriorOffPath, mXyzPriorOffPathNew,
           mAnatomicalPrior, mAnatomicalPriorNew,
           mShapePrior, mShapePriorNew,
           mPosteriorOnPath, mPosteriorOnPathNew, mPosteriorOnPathMap,
           mPosteriorOffPath, mPosteriorOffPathNew;
    std::string mOutDir, mInfoGeneral, mInfoPathway, mInfoMcmc;
    std::vector<bool> mRejectControl;			// [mNumControl]
    std::vector<int> mAcceptCount, mRejectCount,	// [mNumControl]
                     mControlPoints, mControlPointsNew,
                     mPathPoints, mPathPointsNew,
                     mDirLocal, mDirNear;
    std::vector<float> mResolution,			// [3]
                       mProposalStdInit, mProposalStd,	// [mNumControl x 3]
                       mControlPointJumps,		// [mNumControl x 3]
                       mAcceptSpan, mRejectSpan;	// [mNumControl x 3]
    std::vector< std::vector<int> > mAtlasCoords;
    std::vector< std::vector<unsigned int> > mIdsLocal, mIdsNear;
    std::vector< std::vector<float> > mPriorTangent,	// [mNumArc]
                                      mPriorCurvature,	// [mNumArc]
                                      mPriorLocal, 	// [mNumArc x 7]
                                      mPriorNear;	// [mNumArc x 6]
    MRI *mMask, *mRoi1, *mRoi2, *mXyzPrior0, *mXyzPrior1;
    std::ofstream mLog;
    Spline mSpline;
    AffineReg mAffineReg;
#ifndef NO_CVS_UP_IN_HERE
    NonlinReg mNonlinReg;
#endif
    std::vector<MRI *> mAseg;
    std::vector<Aeon> mDwi;

    void ReadControlPoints(const std::string ControlPointFile);
    void ReadProposalStds(const std::string PropStdFile);
    bool InitializeMcmc();
    bool InitializeFixOffMask(int FailSegment);
    bool InitializeFixOffWhite(int FailSegment);
    int FindErrorSegment();
    bool JumpMcmcFull();
    bool JumpMcmcSingle(int ControlIndex);
    bool ProposePathFull();
    bool ProposePathSingle(int ControlIndex);
    void ProposeDiffusionParameters();
    bool AcceptPath(bool UsePriorOnly=false);
    double ComputeXyzPriorOffPath(std::vector<int> &PathAtlasPoints);
    double ComputeXyzPriorOnPath(std::vector<int> &PathAtlasPoints);
    double ComputeAnatomicalPrior(std::vector<int> &PathAtlasPoints);
    double ComputeShapePrior(std::vector<int> &PathAtlasPoints);
    void UpdatePath();
    void UpdateAcceptanceRateFull();
    void UpdateRejectionRateFull();
    void UpdateAcceptRejectRateSingle();
    void UpdateProposalStd();
    void SavePathPosterior(bool IsPathAccepted);
    void SavePath();
    void RemoveDuplicatePathPoints();
    bool IsInMask(std::vector<int>::const_iterator Point);
    bool IsInRoi(std::vector<int>::const_iterator Point, MRI *Roi);
    bool IsZigZag(std::vector<int> &ControlPoints,
                  std::vector<int>::const_iterator FirstPerturbedPoint,
                  std::vector<int>::const_iterator LastPerturbedPoint);
    void MapPointToAtlas(std::vector<int>::iterator OutPoint,
                         std::vector<int>::const_iterator InPoint);
    void LogObjective();
    void LogObjectiveNaN(const unsigned int NumValiData);
};

#endif
