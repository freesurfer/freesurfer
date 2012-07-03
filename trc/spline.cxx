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

#include <spline.h>

using namespace std;

Spline::Spline(const char *ControlPointFile, const char *MaskFile) {
  ReadControlPoints(ControlPointFile);
  ReadMask(MaskFile);
  mVolume = MRIclone(mMask, NULL);
  mAllPoints.clear();
  mArcLength.clear();
  mDerivative1.clear();
  mDerivative2.clear();
  mFiniteDifference1.clear();
  mFiniteDifference2.clear();
  mTangent.clear();
  mNormal.clear();
  mCurvature.clear();
}

Spline::Spline(const vector<int> &ControlPoints, MRI *Mask) {
  SetControlPoints(ControlPoints);
  mMask = MRIcopy(Mask, NULL);
  mVolume = MRIclone(mMask, NULL);
  mAllPoints.clear();
  mArcLength.clear();
  mDerivative1.clear();
  mDerivative2.clear();
  mFiniteDifference1.clear();
  mFiniteDifference2.clear();
  mTangent.clear();
  mNormal.clear();
  mCurvature.clear();
}

Spline::Spline(const int NumControl, MRI *Mask) : mNumControl(NumControl) {
  mMask = MRIcopy(Mask, NULL);
  mVolume = MRIclone(mMask, NULL);
  mControlPoints.clear();
  mAllPoints.clear();
  mArcLength.clear();
  mDerivative1.clear();
  mDerivative2.clear();
  mFiniteDifference1.clear();
  mFiniteDifference2.clear();
  mTangent.clear();
  mNormal.clear();
  mCurvature.clear();
}

Spline::~Spline() {
  mControlPoints.clear();
  mAllPoints.clear();
  mArcLength.clear();
  mDerivative1.clear();
  mDerivative2.clear();
  mFiniteDifference1.clear();
  mFiniteDifference2.clear();
  mTangent.clear();
  mNormal.clear();
  mCurvature.clear();
  MRIfree(&mMask);
  MRIfree(&mVolume);
}

//
// Check if two consecutive control points overlap
//
bool Spline::IsDegenerate() {
  vector<int>::const_iterator icpt;

  for (icpt = mControlPoints.begin(); icpt != mControlPoints.end()-3; icpt += 3)
    if ((icpt[0] == icpt[3]) && (icpt[1] == icpt[4]) && (icpt[2] ==icpt[5]))
      return true;

  return false;
}

//
// Interpolate spline given its control points
//
bool Spline::InterpolateSpline() {
  const int ncpts = mNumControl-1;
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

  for (int kcpt = 1; kcpt <= ncpts; kcpt++) {
    float t = 0, dt = 0, newt;

    // Append current control point to spline
    if (!IsInMask(icpt))
      return false;
    mAllPoints.insert(mAllPoints.end(), icpt, icpt+3);
    mArcLength.push_back(0);
    MRIsetVoxVal(mVolume, icpt[0], icpt[1], icpt[2], 0, 1);

    // Initialize arc length step size
    for (int k=0; k<3; k++)
      dt += pow(icpt[k] - icpt[k+3], 2);
    dt = 1.0 / sqrt(dt);

    while (t < 1) {
      float tmin = 0, tmax = 1;
      bool incstep, decstep;
      vector<int>::const_iterator ipt = mAllPoints.end()-3;

      do {
        newt = t + dt;

        if (newt > 1)  {
          incstep = false;
          decstep = true;
        }
        else {
          // Interpolate new point
          CatmullRomInterpolate(newpoint, newt,
                                (kcpt==1)?icpt:(icpt-3), 
                                icpt, icpt+3,
                                (kcpt==ncpts)?(icpt+3):(icpt+6));

          // Check that the new point is adjacent to the previous point
          incstep = true;
          decstep = false;

          for (int k=0; k<3; k++)
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

      t = newt;

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

//
// Fit spline control points to a curve
//
bool Spline::FitControlPoints(const vector<int> &InputPoints) {
  CatmullRomFit(InputPoints);

  return InterpolateSpline();
}

//
// Compute tangent vectors along spline (either from analytical derivatives
// or from finite difference approximation)
//
void Spline::ComputeTangent(const bool DoAnalytical) {
  vector<float>::const_iterator id1;

  if (DoAnalytical) {	// Compute analytically using parametric form
    int kcpt = 1;
    const int ncpts = mNumControl-1;
    vector<int>::const_iterator icpt = mControlPoints.begin();
    vector<int>::const_iterator ipt;
    vector<float>::const_iterator iarc = mArcLength.begin();
    vector<float> newderiv(3);

    mDerivative1.clear();

    for (ipt = mAllPoints.begin(); ipt != mAllPoints.end()-3; ipt += 3) {
      CatmullRomDerivative1(newderiv, *iarc, (kcpt==1)?icpt:(icpt-3),
                                             icpt, icpt+3,
                                             (kcpt==ncpts)?(icpt+3):(icpt+6));
      mDerivative1.insert(mDerivative1.end(), newderiv.begin(), newderiv.end());

      iarc++;
      if (*iarc == 0) {	// Have reached next control point
        kcpt++;
        icpt += 3;
      }
    }

    // First derivative at final control point
    CatmullRomDerivative1(newderiv, *iarc, icpt-3, icpt, icpt, icpt);
    mDerivative1.insert(mDerivative1.end(), newderiv.begin(), newderiv.end());

    id1 = mDerivative1.begin();
  }
  else {		// Approximate using finite differences
    const unsigned int dt = 3,
                       dt2 = 2 * dt,
                       offset = 3 * dt;
    vector<int>::const_iterator ipt;
    vector<float>::const_iterator iptsm;
    vector<float>::iterator ismooth, idiff;
    vector<float> vecsmooth(mAllPoints.size());

    mFiniteDifference1.resize(mAllPoints.size());

    // Smooth discrete point coordinates
    ipt = mAllPoints.begin();
    ismooth = vecsmooth.begin();

    for ( ; ipt < mAllPoints.begin() + 3; ipt++) {
      *ismooth = (float) *ipt;
      ismooth++;
    }

    for ( ; ipt < mAllPoints.end() - 3; ipt++) {
      *ismooth = (*(ipt-3) + 2 * (*ipt) + *(ipt+3)) / 4.0;
      ismooth++;
    }
  
    for ( ; ipt < mAllPoints.end(); ipt++) {
      *ismooth = (float) *ipt;
      ismooth++;
    }

    // Approximate derivatives by 2nd-order central finite differences,
    // everywhere except at the start and end of the spline, where 1st-order
    // forward and backward finite differences are used, respectively
    iptsm = vecsmooth.begin();
    idiff = mFiniteDifference1.begin();

    for ( ; iptsm < vecsmooth.begin() + offset; iptsm++) {
      *idiff = (*(iptsm + offset) - *iptsm) / (float) dt;
      idiff++;
    }

    for ( ; iptsm < vecsmooth.end() - offset; iptsm++) {
      *idiff = (*(iptsm + offset) - *(iptsm - offset)) / (float) dt2;
      idiff++;
    }

    for ( ; iptsm < vecsmooth.end(); iptsm++) {
      *idiff = (*iptsm - *(iptsm - offset)) / (float) dt;
      idiff++;
    }

    // Smooth finite differences
    idiff = mFiniteDifference1.begin();
    ismooth = vecsmooth.begin();

    for ( ; idiff < mFiniteDifference1.begin() + 3; idiff++) {
      *ismooth = (*idiff + *(idiff+3)) / 2.0;
      ismooth++;
    }

    for ( ; idiff < mFiniteDifference1.end() - 3; idiff++) {
      *ismooth = (*(idiff-3) + *idiff + *(idiff+3)) / 3.0;
      ismooth++;
    }
  
    for ( ; idiff < mFiniteDifference1.end(); idiff++) {
      *ismooth = (*(idiff-3) + *idiff) / 2.0;
      ismooth++;
    }
  
    copy(vecsmooth.begin(), vecsmooth.end(), mFiniteDifference1.begin());

    id1 = mFiniteDifference1.begin();
  }

  // Compute tangent vectors from first derivative
  mTangent.resize(mAllPoints.size());

  for (vector<float>::iterator itang = mTangent.begin(); itang < mTangent.end();
                                                         itang += 3) {
    const float nrm = sqrt(id1[0]*id1[0] + id1[1]*id1[1] + id1[2]*id1[2]);

    if (nrm > 0) 			// Normalize
      for (int k = 0; k < 3; k++)
        itang[k] = id1[k] / nrm;
    else
      fill(itang, itang+3, 0.0);

    id1 += 3;
  }
}

//
// Compute normal vectors along spline (either from analytical derivatives
// or from finite difference approximation)
//
void Spline::ComputeNormal(const bool DoAnalytical) {
  vector<float>::const_iterator id2, itang = mTangent.begin();

  if (DoAnalytical) {	// Compute analytically using parametric form
    int kcpt = 1;
    const int ncpts = mNumControl-1;
    vector<int>::const_iterator icpt = mControlPoints.begin();
    vector<int>::const_iterator ipt;
    vector<float>::const_iterator iarc = mArcLength.begin();
    vector<float> newderiv(3);

    mDerivative2.clear();

    for (ipt = mAllPoints.begin(); ipt != mAllPoints.end()-3; ipt += 3) {
      CatmullRomDerivative2(newderiv, *iarc, (kcpt==1)?icpt:(icpt-3),
                                             icpt, icpt+3,
                                             (kcpt==ncpts)?(icpt+3):(icpt+6));
      mDerivative2.insert(mDerivative2.end(), newderiv.begin(), newderiv.end());

      iarc++;
      if (*iarc == 0) {	// Have reached next control point
        kcpt++;
        icpt += 3;
      }
    }

    // Second derivative at final control point
    CatmullRomDerivative2(newderiv, *iarc, icpt-3, icpt, icpt, icpt);
    mDerivative2.insert(mDerivative2.end(), newderiv.begin(), newderiv.end());

    id2 = mDerivative2.begin();
  }
  else {		// Approximate using finite differences
    const unsigned int dt = 3,
                       dt2 = 2 * dt,
                       offset = 3 * dt;
    vector<float>::const_iterator ipt;
    vector<float>::iterator ismooth, idiff;
    vector<float> vecsmooth(mAllPoints.size());

    mFiniteDifference2.resize(vecsmooth.size());

    // Approximate derivatives by 2nd-order central finite differences,
    // everywhere except at the start and end of the spline, where 1st-order
    // forward and backward finite differences are used, respectively
    ipt = mFiniteDifference1.begin();
    idiff = mFiniteDifference2.begin();

    for ( ; ipt < mFiniteDifference1.begin() + offset; ipt++) {
      *idiff = (*(ipt + offset) - *ipt) / (float) dt;
      idiff++;
    }

    for ( ; ipt < mFiniteDifference1.end() - offset; ipt++) {
      *idiff = (*(ipt + offset) - *(ipt - offset)) / (float) dt2;
      *idiff++;
    }

    for ( ; ipt < mFiniteDifference1.end(); ipt++) {
      *idiff = (*ipt - *(ipt - offset)) / (float) dt;
      idiff++;
    }

    // Smooth finite differences
    idiff = mFiniteDifference2.begin();
    ismooth = vecsmooth.begin();

    for ( ; idiff < mFiniteDifference2.begin() + 3; idiff++) {
      *ismooth = (*idiff + *(idiff+3)) / 2.0;
      ismooth++;
    }

    for ( ; idiff < mFiniteDifference2.end() - 3; idiff++) {
      *ismooth = (*(idiff-3) + *idiff + *(idiff+3)) / 3.0;
      ismooth++;
    }
  
    for ( ; idiff < mFiniteDifference2.end(); idiff++) {
      *ismooth = (*(idiff-3) + *idiff) / 2.0;
      ismooth++;
    }
  
    copy(vecsmooth.begin(), vecsmooth.end(), mFiniteDifference2.begin());

    id2 = mFiniteDifference2.begin();
  }

  // Compute normal vectors from first and second derivative
  mNormal.resize(mAllPoints.size());

  for (vector<float>::iterator inorm = mNormal.begin(); inorm < mNormal.end();
                                                        inorm += 3) {
    const float dot = id2[0]*itang[0] + id2[1]*itang[1] + id2[2]*itang[2];
    float nrm = 0;

    for (int k = 0; k < 3; k++) {
      inorm[k] = id2[k] - dot * itang[k];
      nrm += inorm[k]*inorm[k];
    }

    nrm = sqrt(nrm);

    if (nrm > 0) 			// Normalize
      for (int k = 0; k < 3; k++)
        inorm[k] /= nrm;

    id2 += 3;
    itang += 3;
  }
}

//
// Compute curvature along spline (either from analytical derivatives
// or from finite difference approximation)
//
void Spline::ComputeCurvature(const bool DoAnalytical) {
  vector<float>::const_iterator id1, id2;

  mCurvature.resize(mAllPoints.size() / 3);

  if (DoAnalytical) {	// Compute analytically using parametric form
    id1 = mDerivative1.begin();
    id2 = mDerivative2.begin();
  }
  else {		// Approximate using finite differences
    id1 = mFiniteDifference1.begin();
    id2 = mFiniteDifference2.begin();
  }

  // Curvature = |r' x r''| / |r'|^3
  for (vector<float>::iterator icurv = mCurvature.begin();
                               icurv < mCurvature.end(); icurv++) {
    *icurv = sqrt( ( pow(id1[1] * id2[2] - id1[2] * id2[1], 2) +
                     pow(id1[2] * id2[0] - id1[0] * id2[2], 2) +
                     pow(id1[0] * id2[1] - id1[1] * id2[0], 2) ) / 
                   pow(pow(id1[0], 2) + pow(id1[1], 2) + pow(id1[2], 2), 3) );

    id1 += 3;
    id2 += 3;
  }
}

//
// Read control points from file
//
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

//
// Read mask volume from file
//
void Spline::ReadMask(const char *MaskFile) {
  cout << "Loading spline mask from " << MaskFile << endl;
  mMask = MRIread(MaskFile);
}

//
// Set control points
//
void Spline::SetControlPoints(const std::vector<int> &ControlPoints) {
  mControlPoints.resize(ControlPoints.size());
  copy(ControlPoints.begin(), ControlPoints.end(), mControlPoints.begin());

  mNumControl = mControlPoints.size()/3;
}

//
// Write spline to volume
//
void Spline::WriteVolume(const char *VolumeFile, const bool ShowControls) {
  if (ShowControls)
    for (vector<int>::const_iterator icpt = mControlPoints.begin();
                                     icpt != mControlPoints.end(); icpt += 3)
      MRIsetVoxVal(mVolume, icpt[0], icpt[1], icpt[2], 0, 2);

  cout << "Writing spline volume to " << VolumeFile << endl;
  MRIwrite(mVolume, VolumeFile);
}

//
// Write spline coordinates to file
//
void Spline::WriteAllPoints(const char *TextFile) {
  ofstream outfile(TextFile, ios::out);
  if (!outfile) {
    cout << "ERROR: Could not open " << TextFile << " for writing" << endl;
    exit(1);
  }

  cout << "Writing coordinates along spline to " << TextFile << endl;

  for (vector<int>::const_iterator ipt = mAllPoints.begin();
                                   ipt < mAllPoints.end(); ipt += 3)
    outfile << ipt[0] << " " << ipt[1] << " " << ipt[2] << endl;
}

//
// Write tangent vectors along spline to file
//
void Spline::WriteTangent(const char *TextFile) {
  ofstream outfile(TextFile, ios::out);
  if (!outfile) {
    cout << "ERROR: Could not open " << TextFile << " for writing" << endl;
    exit(1);
  }

  cout << "Writing tangent vectors along spline to " << TextFile << endl;

  for (vector<float>::const_iterator itang = mTangent.begin();
                                     itang < mTangent.end(); itang += 3)
    outfile << itang[0] << " " << itang[1] << " " << itang[2] << endl;
}

//
// Write normal vectors along spline to file
//
void Spline::WriteNormal(const char *TextFile) {
  ofstream outfile(TextFile, ios::out);
  if (!outfile) {
    cout << "ERROR: Could not open " << TextFile << " for writing" << endl;
    exit(1);
  }

  cout << "Writing normal vectors along spline to " << TextFile << endl;

  for (vector<float>::const_iterator inorm = mNormal.begin();
                                     inorm < mNormal.end(); inorm += 3)
    outfile << inorm[0] << " " << inorm[1] << " " << inorm[2] << endl;
}

//
// Write curvatures along spline to file
//
void Spline::WriteCurvature(const char *TextFile) {
  ofstream outfile(TextFile, ios::out);
  if (!outfile) {
    cout << "ERROR: Could not open " << TextFile << " for writing" << endl;
    exit(1);
  }

  cout << "Writing curvatures along spline to " << TextFile << endl;

  for (vector<float>::const_iterator icurv = mCurvature.begin();
                                     icurv < mCurvature.end(); icurv++)
    outfile << *icurv << endl;
}

//
// Write the intensity values of each of a set of input volumes along the spline
//
void Spline::WriteValues(vector<MRI *> &ValueVolumes, const char *TextFile) {
  ofstream outfile(TextFile, ios::app);
  if (!outfile) {
    cout << "ERROR: Could not open " << TextFile << " for writing" << endl;
    exit(1);
  }

  cout << "Writing values along spline to " << TextFile << endl;

  for (vector<int>::const_iterator ipt = mAllPoints.begin();
                                   ipt < mAllPoints.end(); ipt += 3) {
    outfile << ipt[0] << " " << ipt[1] << " " << ipt[2];

    for (vector<MRI *>::const_iterator ivol = ValueVolumes.begin();
                                       ivol < ValueVolumes.end(); ivol++)
      outfile << " " << MRIgetVoxVal(*ivol, ipt[0], ipt[1], ipt[2], 0);

    outfile << endl;
  }
}

//
// Compute average intensity of each of a set of input volumes along the spline
//
vector<float> Spline::ComputeAvg(vector<MRI *> &ValueVolumes) {
  int nvox = (int) mAllPoints.size()/3;
  vector<float> avg(ValueVolumes.size(), 0);
  vector<float>::iterator iavg;

  for (vector<int>::const_iterator ipt = mAllPoints.begin();
                                   ipt < mAllPoints.end(); ipt += 3) {
    iavg = avg.begin();

    for (vector<MRI *>::const_iterator ivol = ValueVolumes.begin();
                                       ivol < ValueVolumes.end(); ivol++) {
      *iavg += MRIgetVoxVal(*ivol, ipt[0], ipt[1], ipt[2], 0);
      iavg++;
    }
  }

  for (iavg = avg.begin(); iavg < avg.end(); iavg++)
    *iavg /= nvox;

  return avg;
}

//
// Print control point coordinates
//
void Spline::PrintControlPoints() {
  for (vector<int>::const_iterator icpt = mControlPoints.begin();
                                   icpt != mControlPoints.end(); icpt += 3)
    cout << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;
}

//
// Print spline coordinates
//
void Spline::PrintAllPoints() {
  for (vector<int>::const_iterator ipt = mAllPoints.begin();
                                   ipt != mAllPoints.end(); ipt += 3)
    cout << ipt[0] << " " << ipt[1] << " " << ipt[2] << endl;
}

//
// Print tangent vectors along spline
//
void Spline::PrintTangent() {
  for (vector<float>::const_iterator itang = mTangent.begin();
                                     itang != mTangent.end(); itang += 3)
    cout << itang[0] << " " << itang[1] << " " << itang[2] << endl;
}

//
// Print normal vectors along spline
//
void Spline::PrintNormal() {
  for (vector<float>::const_iterator inorm = mNormal.begin();
                                     inorm != mNormal.end(); inorm += 3)
    cout << inorm[0] << " " << inorm[1] << " " << inorm[2] << endl;
}

//
// Print curvatures along spline
//
void Spline::PrintCurvature() {
  for (vector<float>::const_iterator icurv = mCurvature.begin();
                                     icurv != mCurvature.end(); icurv++)
    cout << *icurv << endl;
}

//
// Return pointer to start of control points
//
vector<int>::const_iterator Spline::GetControlPointsBegin() {
  return mControlPoints.begin();
}

//
// Return pointer to end of control points
//
vector<int>::const_iterator Spline::GetControlPointsEnd() {
  return mControlPoints.end();
}

//
// Return pointer to beginning of control points
//
vector<int>::const_iterator Spline::GetAllPointsBegin() {
  return mAllPoints.begin();
}

//
// Return pointer to end of control points
//
vector<int>::const_iterator Spline::GetAllPointsEnd() {
  return mAllPoints.end();
}

//
// Return pointer to beginning of tangent vectors
//
vector<float>::const_iterator Spline::GetTangentBegin() {
  return mTangent.begin();
}

//
// Return pointer to end of tangent vectors
//
vector<float>::const_iterator Spline::GetTangentEnd() {
  return mTangent.end();
}

//
// Return pointer to beginning of normal vectors
//
vector<float>::const_iterator Spline::GetNormalBegin() {
  return mNormal.begin();
}

//
// Return pointer to end of normal vectors
//
vector<float>::const_iterator Spline::GetNormalEnd() {
  return mNormal.end();
}

//
// Return pointer to beginning of curvatures
//
vector<float>::const_iterator Spline::GetCurvatureBegin() {
  return mCurvature.begin();
}

//
// Return pointer to end of curvatures
//
vector<float>::const_iterator Spline::GetCurvatureEnd() {
  return mCurvature.end();
}

//
// Intepolate Catmull-Rom spline from control points
//
void Spline::CatmullRomInterpolate(vector<int> &InterpPoint,
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

  for (int k = 0; k < 3; k++)
    InterpPoint[k] = (int) round(a * ControlPoint1[k] + b * ControlPoint2[k] +
                                 c * ControlPoint3[k] + d * ControlPoint4[k]);
}

//
// Compute Catmull-Rom spline first derivative vector from control points
//
void Spline::CatmullRomDerivative1(vector<float> &InterpDerivative,
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

  for (int k = 0; k < 3; k++)
    InterpDerivative[k] = a * ControlPoint1[k] + b * ControlPoint2[k] +
                          c * ControlPoint3[k] + d * ControlPoint4[k];
}

//
// Compute Catmull-Rom spline second derivative vector from control points
//
void Spline::CatmullRomDerivative2(vector<float> &InterpDerivative,
                                   const float t,
                                   vector<int>::const_iterator ControlPoint1,
                                   vector<int>::const_iterator ControlPoint2,
                                   vector<int>::const_iterator ControlPoint3,
                                   vector<int>::const_iterator ControlPoint4) {
  const float a = -3*t + 2,
              b = 9*t - 5,
              c = -9*t + 4,
              d = 3*t - 1;

  for (int k = 0; k < 3; k++)
    InterpDerivative[k] = a * ControlPoint1[k] + b * ControlPoint2[k] +
                          c * ControlPoint3[k] + d * ControlPoint4[k];
}

//
// Fit Catmull-Rom spline to a curve
//
void Spline::CatmullRomFit(const vector<int> &InputPoints) {
  const int npts = (int) InputPoints.size() / 3;
  const float smax = (float) (mNumControl - 1),
              ds = smax / (npts - 1);
  float s = ds; 
  vector<int>::const_iterator ipt = InputPoints.begin() + 3;
  vector<int>::iterator icpt;
  MATRIX *A = NULL, *Ap = NULL, *y = NULL, *x = NULL;

  // Fit all but the first and last control points,
  // which will be set equal to the first and last of the input points
  A  = MatrixAlloc(npts-2, mNumControl-2, MATRIX_REAL);
  Ap = MatrixAlloc(mNumControl-2, npts-2, MATRIX_REAL);
  y  = MatrixAlloc(npts-2, 3, MATRIX_REAL);
  x  = MatrixAlloc(mNumControl-2, 3, MATRIX_REAL);

  for (int irow = 1; irow < npts-1; irow++) {
    const int idx = (int) floor(s);
    const float t = s - floor(s),
                t2 = t*t,
                t3 = t2*t,
                a = -.5*(t3 + t) + t2,
                b = 1.5*t3 - 2.5*t2 + 1,
                c = -1.5*t3 + 2*t2 + .5*t,
                d = .5*(t3 - t2);

    if (idx == 0) {
      const float ab = a + b;
      vector<int>::const_iterator ipt1 = InputPoints.begin();

      for (int k = 1; k < 4; k++) {
        y->rptr[irow][k] = *ipt - ab * (*ipt1);
        ipt++;
        ipt1++;
      }

      A->rptr[irow][1] = c;
      A->rptr[irow][2] = d;
    }
    else if (idx == 1) {
      vector<int>::const_iterator ipt1 = InputPoints.begin();

      for (int k = 1; k < 4; k++) {
        y->rptr[irow][k] = *ipt - a * (*ipt1);
        ipt++;
        ipt1++;
      }

      A->rptr[irow][1] = b;
      A->rptr[irow][2] = c;
      A->rptr[irow][3] = d;
    }
    else if (idx == mNumControl - 3) {
      vector<int>::const_iterator iptn = InputPoints.end() - 3;

      for (int k = 1; k < 4; k++) {
        y->rptr[irow][k] = *ipt - d * (*iptn);
        ipt++;
        iptn++;
      }

      A->rptr[irow][mNumControl - 4] = a;
      A->rptr[irow][mNumControl - 3] = b;
      A->rptr[irow][mNumControl - 2] = c;
    }
    else if (idx == mNumControl - 2) {
      const float cd = c + d;
      vector<int>::const_iterator iptn = InputPoints.end() - 3;

      for (int k = 1; k < 4; k++) {
        y->rptr[irow][k] = *ipt - cd * (*iptn);
        ipt++;
        iptn++;
      }

      A->rptr[irow][mNumControl - 3] = a;
      A->rptr[irow][mNumControl - 2] = b;
    }
    else {
      for (int k = 1; k < 4; k++) {
        y->rptr[irow][k] = *ipt;
        ipt++;
      }

      A->rptr[irow][idx-1] = a;
      A->rptr[irow][idx]   = b;
      A->rptr[irow][idx+1] = c;
      A->rptr[irow][idx+2] = d;
    }

    s += ds;
  }

  // Find least-squares fit of all control points but the first and last
  MatrixPseudoInverse(A, Ap);
  MatrixMultiply(Ap, y, x);

  // Copy fitted control points and append the first and last point
  mControlPoints.resize(mNumControl * 3);

  icpt = mControlPoints.begin();
  ipt = InputPoints.begin();

  for (ipt = InputPoints.begin(); ipt < InputPoints.begin() + 3; ipt++) {
    *icpt = *ipt;
    *icpt++;
  }

  for (int irow = 1; irow < mNumControl-1; irow++)
    for (int k = 1; k < 4; k++) {
      *icpt = (int) round(x->rptr[irow][k]);
      icpt++;
    }

  for (ipt = InputPoints.end() - 3; ipt < InputPoints.end(); ipt++) {
    *icpt = *ipt;
    *icpt++;
  }
}

//
// Check if a point is in mask
//
bool Spline::IsInMask(vector<int>::const_iterator Point) {
  // Check that point is inside mask and has not already been traversed
  return (Point[0] > -1) && (Point[0] < mMask->width) && 
         (Point[1] > -1) && (Point[1] < mMask->height) && 
         (Point[2] > -1) && (Point[2] < mMask->depth) && 
         (MRIgetVoxVal(mMask, Point[0], Point[1], Point[2], 0) > 0) &&
         (MRIgetVoxVal(mVolume, Point[0], Point[1], Point[2], 0) == 0);
}

