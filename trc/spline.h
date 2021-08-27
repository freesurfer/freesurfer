/**
 * @brief A Catmull-Rom spline
 *
 * A Catmull-Rom spline
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

#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <limits>
#include <math.h>
#include "mri.h"

void CurveFiniteDifferences(std::vector<float> &DiffPoints,
                            const std::vector<float> &CurvePoints,
                            const unsigned int DiffStep);

void CurveSmooth(std::vector<float> &SmoothPoints,
                 const std::vector<int> &DiscretePoints);

std::vector<int> CurveFill(const std::vector<int> &InPoints);

class Spline {
  public:
    Spline(const std::string ControlPointFile, const std::string MaskFile);
    Spline(const std::vector<int> &ControlPoints, MRI *Mask);
    Spline(const int NumControl, MRI *Mask);
    Spline();
    ~Spline();
    bool IsDegenerate();
    bool InterpolateSpline();
    bool FitControlPoints(const std::vector<int> &InputPoints);
    unsigned int PointToSegment(unsigned int PointIndex);
    void ComputeTangent(const bool DoAnalytical=true);
    void ComputeNormal(const bool DoAnalytical=true);
    void ComputeCurvature(const bool DoAnalytical=true);
    void ReadControlPoints(const std::string ControlPointFile);
    void ReadMask(const std::string MaskFile);
    void SetControlPoints(const std::vector<int> &ControlPoints);
    void SetMask(MRI *Mask);
    void WriteVolume(const std::string VolumeFile,
                     const bool ShowControls=false);
    void WriteAllPoints(const std::string TextFile);
    void WriteTangent(const std::string TextFile);
    void WriteNormal(const std::string TextFile);
    void WriteCurvature(const std::string TextFile);
    void WriteValues(std::vector<MRI *> &ValueVolumes,
                     const std::string TextFile);
    std::vector<float> ComputeAvg(std::vector<MRI *> &ValueVolumes);
    void PrintControlPoints();
    void PrintAllPoints();
    void PrintTangent();
    void PrintNormal();
    void PrintCurvature();
    std::vector<int>::const_iterator GetControlPointsBegin();
    std::vector<int>::const_iterator GetControlPointsEnd();
    std::vector<int>::const_iterator GetAllPointsBegin();
    std::vector<int>::const_iterator GetAllPointsEnd();
    std::vector<float>::const_iterator GetTangentBegin();
    std::vector<float>::const_iterator GetTangentEnd();
    std::vector<float>::const_iterator GetNormalBegin();
    std::vector<float>::const_iterator GetNormalEnd();
    std::vector<float>::const_iterator GetCurvatureBegin();
    std::vector<float>::const_iterator GetCurvatureEnd();

  private:
    int mNumControl;
    std::vector<int> mControlPoints, mAllPoints;
    std::vector<float> mArcLength,
                       mDerivative1, mDerivative2,
                       mFiniteDifference1, mFiniteDifference2,
                       mTangent, mNormal, mCurvature;
    MRI *mMask, *mVolume;

    void CatmullRomInterpolate(std::vector<int> &InterpPoint,
                               const float t,
                               std::vector<int>::const_iterator ControlPoint1,
                               std::vector<int>::const_iterator ControlPoint2,
                               std::vector<int>::const_iterator ControlPoint3,
                               std::vector<int>::const_iterator ControlPoint4);
    void CatmullRomDerivative1(std::vector<float> &InterpDerivative,
                               const float t,
                               std::vector<int>::const_iterator ControlPoint1,
                               std::vector<int>::const_iterator ControlPoint2,
                               std::vector<int>::const_iterator ControlPoint3,
                               std::vector<int>::const_iterator ControlPoint4);
    void CatmullRomDerivative2(std::vector<float> &InterpDerivative,
                               const float t,
                               std::vector<int>::const_iterator ControlPoint1,
                               std::vector<int>::const_iterator ControlPoint2,
                               std::vector<int>::const_iterator ControlPoint3,
                               std::vector<int>::const_iterator ControlPoint4);
    void CatmullRomFit(const std::vector<int> &InputPoints);
    void CatmullRomFit(const std::vector<int> &InputPoints,
                       const std::vector<float> &ArcLengthParameter);
    bool IsInMask(std::vector<int>::const_iterator Point);
};

#endif

