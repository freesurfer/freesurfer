/**
 * @brief Holds utilities for probabilistic tractography
 *
 * Holds utilities for probabilistic tractography
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

#ifndef VIAL_H
#define VIAL_H

//define NO_CVS_UP_IN_HERE

#ifndef NO_CVS_UP_IN_HERE
// Needed for CVS - these must be included first or else they don't compile
#include <cmath>
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include <mpi.h>
#include "../fem_elastic/morph.h"
#include "../fem_elastic/surf_utils.h"
#include "../fem_elastic/morph_utils.h"
#include "gcamorph.h"
#endif

#include <vector>
#include <iostream>
#include <fstream>
#include "mri.h"

class AffineReg {
  public:
    AffineReg();
    AffineReg(std::vector<float> &InToOut);
    ~AffineReg();
    bool IsEmpty();
    bool IsInvEmpty();
    void ReadXfm(const std::string XfmFile, const MRI *InRefVol, 
                                            const MRI *OutRefVol);
    void ApplyXfm(std::vector<float> &OutPoint,
                  std::vector<float>::const_iterator InPoint);
    void ApplyXfmInv(std::vector<float> &OutPoint,
                     std::vector<float>::const_iterator InPoint);
    void DecomposeXfm();
    void PrintScale();
    void PrintShear();
    void PrintRotate();
    std::vector<float>::const_iterator GetTranslate();
    std::vector<float>::const_iterator GetRotate();
    std::vector<float>::const_iterator GetShear();
    std::vector<float>::const_iterator GetScale();

  private:
    std::vector<float> mInToOut, mOutToIn,			// [4 x 4]
                       mInVoxelSize, mOutVoxelSize,		// [3]
                       mTranslate, mRotate, mScale, mShear;	// [3]
};

#ifndef NO_CVS_UP_IN_HERE
class NonlinReg {
  public:
    NonlinReg();
    ~NonlinReg();
    bool IsEmpty();
    void ReadXfm(const std::string XfmFile, MRI *RefVol);
    void ApplyXfm(std::vector<float> &OutPoint,
                  std::vector<float>::const_iterator InPoint);
    void ApplyXfmInv(std::vector<float> &OutPoint,
                     std::vector<float>::const_iterator InPoint);

  private:
    GCAM *mMorph;
};
#endif

class StreamSet {
  public:
    StreamSet();
    ~StreamSet();
    void SetLengths(std::vector< std::vector<float> > &Streamlines);
    void SetNumSteps(const unsigned int NumSteps);
    unsigned int SetNumStepsAvgLength();
    unsigned int SetNumStepsMinLength();
    void ComputeSteps();
    void ComputeMeanStreamline(std::vector< std::vector<float> > &Streamlines);
    void ComputeMeanStdStreamline(std::vector< std::vector<float> >
                                                                 &Streamlines);
    unsigned int FindMeanNearestStreamline(std::vector< std::vector<float> >
                                                                 &Streamlines);
    unsigned int GetNumSteps();
    std::vector<float>::const_iterator GetStreamlineMean();
    std::vector<float>::const_iterator GetStreamlineStd();

  private:
    unsigned int mLengthMax, mLengthMin, mNumSteps, mMeanNearestStr;
    float mLengthAvg;
    std::vector<bool> mIsOut;
    std::vector<unsigned int> mLengths;
    std::vector<float> mSteps, mStreamlineMean, mStreamlineStd;
};

#endif

