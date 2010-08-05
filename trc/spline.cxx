/**
 * @file  spline.cxx
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

#include <spline.h>

using namespace std;

Spline::Spline(const char *ControlPointFile, const char *MaskFile) {
  ReadControlPoints(ControlPointFile);
  ReadMask(MaskFile);
  mVolume = MRIclone(mMask, NULL);
  mAllPoints.clear();
  mTangent.clear();
  mArcLength.clear();
}

Spline::Spline(const vector<int> &ControlPoints, MRI *Mask) {
  SetControlPoints(ControlPoints);
  mMask = MRIcopy(Mask, NULL);
  mVolume = MRIclone(mMask, NULL);
  mAllPoints.clear();
  mTangent.clear();
  mArcLength.clear();
}

Spline::~Spline() {
  mControlPoints.clear();
  mAllPoints.clear();
  mTangent.clear();
  mArcLength.clear();
  MRIfree(&mMask);
  MRIfree(&mVolume);
}

bool Spline::IsDegenerate() {
  vector<int>::const_iterator icpt;

  for (icpt = mControlPoints.begin(); icpt != mControlPoints.end()-3; icpt += 3)
    if ((icpt[0] == icpt[3]) && (icpt[1] == icpt[4]) && (icpt[2] ==icpt[5]))
      return true;

  return false;
}

bool Spline::InterpolateSpline() {
  const unsigned int ncpts = mNumControl-1;
  vector<int>::const_iterator icpt = mControlPoints.begin();
  vector<int> newpoint(3);

  if (IsDegenerate()) {
    cout << "ERROR: Degenerate spline segment" << endl;
    PrintControlPoints();
    return false;
  }

  mAllPoints.clear();
  mArcLength.clear();
  MRIclear(mVolume);

  for (unsigned int kcpt = 1; kcpt <= ncpts; kcpt++) {
    float t = 0, dt = 0;

    // Append current control point to spline
    if (!IsInMask(icpt))
      return false;
    mAllPoints.insert(mAllPoints.end(), icpt, icpt+3);
    mArcLength.push_back(0);
    MRIsetVoxVal(mVolume, icpt[0], icpt[1], icpt[2], 0, 1);

    // Initialize arc length step size
    for (unsigned int k=0; k<3; k++)
      dt += pow(icpt[k] - icpt[k+3], 2);
    dt = 1.0 / sqrt(dt);

    while (t < 1) {
      float tmin = 0, tmax = 1;
      bool incstep, decstep;
      vector<int>::const_iterator ipt = mAllPoints.end()-3;

      do {
        // Interpolate new point
        CatmullRomInterp(newpoint, t+dt, (kcpt==1)?icpt:(icpt-3), 
                                         icpt, icpt+3,
                                         (kcpt==ncpts)?(icpt+3):(icpt+6));

        // Check that the new point is adjacent to the previous point
        incstep = true;
        decstep = false;

        for (unsigned int k=0; k<3; k++)
          switch(abs(ipt[k] - newpoint[k])) {
            case 0:
              break;
            case 1:
              incstep = false;
              break;
            default:
              incstep = false;
              decstep = true;
              break;
          }

        // Adjust arc length step size if neccessary
        if (incstep) {
          tmin = dt;
          dt = (dt + tmax)/2;
        }
        else if (decstep) {
          tmax = dt;
          dt = (tmin + dt)/2;
        }
      } while (incstep || decstep);

      t += dt;

      // Check if the next control point has been reached
      if ((newpoint[0] == icpt[3]) && (newpoint[1] == icpt[4])
                                   && (newpoint[2] == icpt[5]))
        break;

      // Append interpolated point to spline
      if (!IsInMask(newpoint.begin()))
        return false;
      mAllPoints.insert(mAllPoints.end(), newpoint.begin(), newpoint.end());
      mArcLength.push_back(t);
      MRIsetVoxVal(mVolume, newpoint[0], newpoint[1], newpoint[2], 0, 1);
    }

    icpt += 3;
  }

  // Append final control point to spline
  if (!IsInMask(icpt))
    return false;
  mAllPoints.insert(mAllPoints.end(), icpt, icpt+3);
  mArcLength.push_back(0);
  MRIsetVoxVal(mVolume, icpt[0], icpt[1], icpt[2], 0, 1);

  return true;
}

void Spline::ComputeTangent() {
  unsigned int kcpt = 1;
  const unsigned int ncpts = mNumControl-1;
  vector<int>::const_iterator icpt = mControlPoints.begin();
  vector<int>::const_iterator ipt;
  vector<float>::const_iterator iarc = mArcLength.begin();
  vector<float> newtangent(3);

  mTangent.clear();

  for (ipt = mAllPoints.begin(); ipt != mAllPoints.end()-3; ipt += 3) {
    CatmullRomTangent(newtangent, *iarc, (kcpt==1)?icpt:(icpt-3),
                                         icpt, icpt+3,
                                         (kcpt==ncpts)?(icpt+3):(icpt+6));
    mTangent.insert(mTangent.end(), newtangent.begin(), newtangent.end());

    iarc++;
    if (*iarc == 0) {	// Have reached next control point
      kcpt++;
      icpt += 3;
    }
  }

  // Tangent at final control point
  CatmullRomTangent(newtangent, *iarc, icpt-3, icpt, icpt, icpt);
  mTangent.insert(mTangent.end(), newtangent.begin(), newtangent.end());
}

void Spline::ComputeCurvature() {
  unsigned int kcpt = 1;
  const unsigned int ncpts = mNumControl-1;
  vector<int>::const_iterator icpt = mControlPoints.begin();
  vector<int>::const_iterator ipt;
  vector<float>::const_iterator iarc = mArcLength.begin();
  vector<float> newcurv(3);

  mCurvature.clear();

  for (ipt = mAllPoints.begin(); ipt != mAllPoints.end()-3; ipt += 3) {
    CatmullRomCurvature(newcurv, *iarc, (kcpt==1)?icpt:(icpt-3),
                                         icpt, icpt+3,
                                         (kcpt==ncpts)?(icpt+3):(icpt+6));
    mCurvature.insert(mCurvature.end(), newcurv.begin(), newcurv.end());

    iarc++;
    if (*iarc == 0) {	// Have reached next control point
      kcpt++;
      icpt += 3;
    }
  }

  // Curvature at final control point
  CatmullRomCurvature(newcurv, *iarc, icpt-3, icpt, icpt, icpt);
  mCurvature.insert(mCurvature.end(), newcurv.begin(), newcurv.end());
}

void Spline::ReadControlPoints(const char *ControlPointFile) {
  float coord;
  ifstream infile(ControlPointFile, ios::in);

  if (!infile) {
    cout << "ERROR: Could not open " << ControlPointFile << endl;
    exit(1);
  }

  cout << "Loading spline control points from " << ControlPointFile << endl;
  mControlPoints.clear();
  while (infile >> coord)
    mControlPoints.push_back((int) round(coord));

  if (mControlPoints.size() % 3 != 0) {
    cout << "ERROR: File " << ControlPointFile
         << " must contain triplets of coordinates" << endl;
    exit(1);
  }

  mNumControl = mControlPoints.size()/3;
}

void Spline::ReadMask(const char *MaskFile) {
  cout << "Loading spline mask from " << MaskFile << endl;
  mMask = MRIread(MaskFile);
}

void Spline::SetControlPoints(const std::vector<int> &ControlPoints) {
  mControlPoints.resize(ControlPoints.size());
  copy(ControlPoints.begin(), ControlPoints.end(), mControlPoints.begin());

  mNumControl = mControlPoints.size()/3;
}

void Spline::WriteVolume(const char *VolumeFile, const bool ShowControls) {
  if (ShowControls)
    for (vector<int>::const_iterator icpt = mControlPoints.begin();
                                     icpt != mControlPoints.end(); icpt += 3)
      MRIsetVoxVal(mVolume, icpt[0], icpt[1], icpt[2], 0, 2);

  cout << "Writing spline volume to " << VolumeFile << endl;
  MRIwrite(mVolume, VolumeFile);
}

void Spline::PrintControlPoints() {
  vector<int>::const_iterator icpt;

  for (icpt = mControlPoints.begin(); icpt != mControlPoints.end(); icpt += 3)
    cout << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;
}

void Spline::PrintAllPoints() {
  vector<int>::const_iterator ipt;

  for (ipt = mAllPoints.begin(); ipt != mAllPoints.end(); ipt += 3)
    cout << ipt[0] << " " << ipt[1] << " " << ipt[2] << endl;
}

void Spline::PrintTangent() {
  vector<float>::const_iterator ipt;

  for (ipt = mTangent.begin(); ipt != mTangent.end(); ipt += 3)
    cout << ipt[0] << " " << ipt[1] << " " << ipt[2] << endl;
}

void Spline::PrintCurvature() {
  vector<float>::const_iterator ipt;

  for (ipt = mCurvature.begin(); ipt != mCurvature.end(); ipt += 3)
    cout << ipt[0] << " " << ipt[1] << " " << ipt[2] << endl;
}

vector<int>::const_iterator Spline::GetAllPointsBegin() {
  return mAllPoints.begin();
}

vector<int>::const_iterator Spline::GetAllPointsEnd() {
  return mAllPoints.end();
}

vector<float>::const_iterator Spline::GetTangentBegin() {
  return mTangent.begin();
}

vector<float>::const_iterator Spline::GetTangentEnd() {
  return mTangent.end();
}

vector<float>::const_iterator Spline::GetCurvatureBegin() {
  return mCurvature.begin();
}

vector<float>::const_iterator Spline::GetCurvatureEnd() {
  return mCurvature.end();
}

void Spline::CatmullRomInterp(vector<int> &InterpPoint,
                              const float t,
                              vector<int>::const_iterator ControlPoint1,
                              vector<int>::const_iterator ControlPoint2,
                              vector<int>::const_iterator ControlPoint3,
                              vector<int>::const_iterator ControlPoint4) {
  const float t2 = t*t,
              t3 = t2*t,
              a = -.5*(t3 + t) + t2,
              b = 1.5*t3 - 2.5*t2 + 1,
              c = -1.5*t3 + 2*t2 + .5*t,
              d = .5*(t3 - t2);

  for (unsigned int k = 0; k < 3; k++)
    InterpPoint[k] = (int) round(a * ControlPoint1[k] + b * ControlPoint2[k] +
                                 c * ControlPoint3[k] + d * ControlPoint4[k]);
}

void Spline::CatmullRomTangent(vector<float> &InterpTangent,
                               const float t,
                               vector<int>::const_iterator ControlPoint1,
                               vector<int>::const_iterator ControlPoint2,
                               vector<int>::const_iterator ControlPoint3,
                               vector<int>::const_iterator ControlPoint4) {
  const float t2 = t*t,
              a = -1.5*t2 + 2*t - .5,
              b = 4.5*t2 - 5*t,
              c = -4.5*t2 + 4*t + .5,
              d = 1.5*t2 - t;

  for (unsigned int k = 0; k < 3; k++)
    InterpTangent[k] = a * ControlPoint1[k] + b * ControlPoint2[k] +
                       c * ControlPoint3[k] + d * ControlPoint4[k];
}

void Spline::CatmullRomCurvature(vector<float> &InterpCurvature,
                                 const float t,
                                 vector<int>::const_iterator ControlPoint1,
                                 vector<int>::const_iterator ControlPoint2,
                                 vector<int>::const_iterator ControlPoint3,
                                 vector<int>::const_iterator ControlPoint4) {
  const float a = -3*t + 2,
              b = 9*t - 5,
              c = -9*t + 4,
              d = 3*t - 1;

  for (unsigned int k = 0; k < 3; k++)
    InterpCurvature[k] = a * ControlPoint1[k] + b * ControlPoint2[k] +
                         c * ControlPoint3[k] + d * ControlPoint4[k];
}

bool Spline::IsInMask(vector<int>::const_iterator Point) {
  // Check that point is inside mask and has not already been traversed
  return (Point[0] > -1) && (Point[0] < mMask->width) && 
         (Point[1] > -1) && (Point[1] < mMask->height) && 
         (Point[2] > -1) && (Point[2] < mMask->depth) && 
         (MRIgetVoxVal(mMask, Point[0], Point[1], Point[2], 0) > 0) &&
         (MRIgetVoxVal(mVolume, Point[0], Point[1], Point[2], 0) == 0);
}

