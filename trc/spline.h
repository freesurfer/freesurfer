/**
 * @file  spline.h
 * @brief A Catmull-Rom spline
 *
 * A Catmull-Rom spline
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

#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "mri.h"

class Spline {
  public:
    Spline(const char *ControlPointFile, const char *MaskFile);
    Spline(const std::vector<int> &ControlPoints, MRI *Mask);
    ~Spline();
    bool IsDegenerate();
    bool InterpolateSpline();
    void ComputeTangent();
    void ComputeNormal();
    void ComputeCurvature();
    void ReadControlPoints(const char *ControlPointFile);
    void ReadMask(const char *MaskFile);
    void SetControlPoints(const std::vector<int> &ControlPoints);
    void WriteVolume(const char *VolumeFile, const bool ShowControls=false);
    void WriteValues(std::vector<MRI *> &ValueVolumes, const char *TextFile);
    std::vector<float> ComputeAvg(std::vector<MRI *> &ValueVolumes);
    void PrintControlPoints();
    void PrintAllPoints();
    void PrintTangent();
    void PrintNormal();
    void PrintCurvature();
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
    std::vector<int> mControlPoints;
    std::vector<int> mAllPoints;
    std::vector<float> mArcLength;
    std::vector<float> mTangent;
    std::vector<float> mNormal;
    std::vector<float> mCurvature;
    MRI *mMask, *mVolume;

    void CatmullRomInterp(std::vector<int> &InterpPoint,
                          const float t,
                          std::vector<int>::const_iterator ControlPoint1,
                          std::vector<int>::const_iterator ControlPoint2,
                          std::vector<int>::const_iterator ControlPoint3,
                          std::vector<int>::const_iterator ControlPoint4);
    void CatmullRomTangent(std::vector<float> &InterpTangent,
                           const float t,
                           std::vector<int>::const_iterator ControlPoint1,
                           std::vector<int>::const_iterator ControlPoint2,
                           std::vector<int>::const_iterator ControlPoint3,
                           std::vector<int>::const_iterator ControlPoint4);
    void CatmullRomNormal(std::vector<float> &InterpNormal,
                          const float t,
                          std::vector<int>::const_iterator ControlPoint1,
                          std::vector<int>::const_iterator ControlPoint2,
                          std::vector<int>::const_iterator ControlPoint3,
                          std::vector<int>::const_iterator ControlPoint4);
    bool IsInMask(std::vector<int>::const_iterator Point);
};

#endif

