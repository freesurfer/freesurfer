/**
 * @brief A Catmull-Rom spline
 *
 * A Catmull-Rom spline
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

#include <spline.h>

using namespace std;

//
// Compute finite differences of a curve
// This is used for computing tangent vectors and curvatures of streamlines
//
void CurveFiniteDifferences(vector<float> &DiffPoints,
                            const vector<float> &CurvePoints,
                            const unsigned int DiffStep) {
  const unsigned int ds2 = 2 * DiffStep,
                     offset = 3 * DiffStep;
  vector<float>::const_iterator ipt;
  vector<float>::iterator iout;
  vector<float> diff(CurvePoints.size());

  if (DiffPoints.size() != CurvePoints.size()) {
    cout << "ERROR: Input and output vector sizes do not match" << endl;
    exit(1);
  }

  // Approximate derivatives by 2nd-order central finite differences,
  // everywhere except at the start and end of the spline, where 1st-order
  // forward and backward finite differences are used, respectively
  ipt = CurvePoints.begin();
  iout = diff.begin();

  for ( ; ipt < CurvePoints.begin() + offset; ipt++) {
    *iout = (*(ipt + offset) - *ipt) / (float) DiffStep;
    iout++;
  }

  for ( ; ipt < CurvePoints.end() - offset; ipt++) {
    *iout = (*(ipt + offset) - *(ipt - offset)) / (float) ds2;
    iout++;
  }

  for ( ; ipt < CurvePoints.end(); ipt++) {
    *iout = (*ipt - *(ipt - offset)) / (float) DiffStep;
    iout++;
  }

  // Smooth finite differences
  ipt = diff.begin();
  iout = DiffPoints.begin();

  for ( ; ipt < diff.begin() + 3; ipt++) {
    *iout = (*ipt + *(ipt+3)) / 2.0;
    iout++;
  }

  for ( ; ipt < diff.end() - 3; ipt++) {
    *iout = (*(ipt-3) + *ipt + *(ipt+3)) / 3.0;
    iout++;
  }

  for ( ; ipt < diff.end(); ipt++) {
    *iout = (*(ipt-3) + *ipt) / 2.0;
    iout++;
  }
}

//
// Smooth a curve
// This is used before computing finite differences on a streamline
//
void CurveSmooth(vector<float> &SmoothPoints,
                 const vector<int> &DiscretePoints) {
  vector<int>::const_iterator ipt = DiscretePoints.begin();
  vector<float>::iterator ismooth = SmoothPoints.begin();

  if (SmoothPoints.size() != DiscretePoints.size()) {
    cout << "ERROR: Input and output vector sizes do not match" << endl;
    exit(1);
  }

  for ( ; ipt < DiscretePoints.begin() + 3; ipt++) {
    *ismooth = (float) *ipt;
    ismooth++;
  }

  for ( ; ipt < DiscretePoints.end() - 3; ipt++) {
    *ismooth = (*(ipt-3) + 2 * (*ipt) + *(ipt+3)) / 4.0;
    ismooth++;
  }

  for ( ; ipt < DiscretePoints.end(); ipt++) {
    *ismooth = (float) *ipt;
    ismooth++;
  }
}

//
// Fill gaps between the points of a curve
// Gaps could result after mapping a curve's points to a higher-resolution space
//
vector<int> CurveFill(const vector<int> &InPoints) {
  vector<int> outpts;
  vector<float> point(3), step(3,0);

  for (vector<int>::const_iterator ipt = InPoints.begin(); ipt < InPoints.end();
                                                           ipt += 3) {
    float dmax = 1;		// This will not remove duplicate points

    if (ipt < InPoints.end() - 3) {
      // Calculate step for filling in gap between points
      for (int k = 0; k < 3; k++) {
        float dist = ipt[k+3] - ipt[k];

        step[k] = dist;
        dist = fabs(dist);

        if (dist > dmax)
          dmax = dist;
      }

      if (dmax > 0)
        for (int k = 0; k < 3; k++)
          step[k] /= dmax;
    }

    for (int k = 0; k < 3; k++)
      point[k] = (float) ipt[k];

    for (int istep = (int) round(dmax); istep > 0; istep--)
      for (int k = 0; k < 3; k++) {
        outpts.push_back((int) round(point[k]));
        point[k] += step[k];
      }
  }

  return outpts;
}

//
// Catmull-Rom cubic spline class
//
Spline::Spline(const string ControlPointFile, const string MaskFile) {
  ReadControlPoints(ControlPointFile);
  mMask = 0;
  mVolume = 0;
  ReadMask(MaskFile);
}

Spline::Spline(const vector<int> &ControlPoints, MRI *Mask) {
  SetControlPoints(ControlPoints);
  mMask = MRIcopy(Mask, NULL);
  mVolume = MRIclone(mMask, NULL);
}

Spline::Spline(const int NumControl, MRI *Mask) : mNumControl(NumControl) {
  mMask = MRIcopy(Mask, NULL);
  mVolume = MRIclone(mMask, NULL);
}

Spline::Spline() : mNumControl(0), mMask(0), mVolume(0) { }

Spline::~Spline() {
  if (mMask)   MRIfree(&mMask);
  if (mVolume) MRIfree(&mVolume);
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
// Uses "dominant" points on the curve (cf. Park and Lee, Computer-Aided Design 
// 39:439-451, 2007)
//
bool Spline::FitControlPoints(const vector<int> &InputPoints) {
  bool success;
  const int npts = (int) InputPoints.size() / 3;
  const float seglen = (npts - 1) / (float) (mNumControl - 1),
              ds = 1 / seglen;
  const int lentot = npts - 1,
            mingap = (int) ceil(seglen / 2);
  float s = 0;
  vector<bool> isfinal(mNumControl, false);
  vector<int> peaksegs;
  vector<float> peakcurvs, arcparam(npts);
  vector< vector<int>::const_iterator > peakpts, dompts(mNumControl);
  vector<float>::iterator iarc = arcparam.begin();

  cout << "INFO: Minimum allowable distance between dominant points is "
       << mingap << endl;

  mAllPoints.resize(InputPoints.size());
  copy(InputPoints.begin(), InputPoints.end(), mAllPoints.begin());

  ComputeTangent(false);
  ComputeNormal(false);
  ComputeCurvature(false);

  const float curvtot = accumulate(mCurvature.begin() + 1,
                                   mCurvature.end() - 1, 0.0),
              mincurv = accumulate(mCurvature.begin() + mingap,
                                   mCurvature.end() - mingap, 0.0) /
                        (npts - 2*mingap);

  cout << "INFO: Minimum allowable curvature at a dominant point is "
       << mincurv << endl;

  // Split curve into segments of equal length
  // We will pick only one dominant point in each of these segments,
  // to ensure that control points end up fairly spread out along the length
  // of the streamline
  for (vector< vector<int>::const_iterator>::iterator idom = dompts.begin();
                                                      idom < dompts.end() - 1;
                                                      idom++) {
    *idom = mAllPoints.begin() + 3 * (int) ceil(s);
    s += seglen;
  }

  *(dompts.end() - 1) = mAllPoints.end() - 3;		// Last point

  // Find points where local peaks in curvature occur
  cout << "INFO: Points where local peaks in curvature occur are" << endl;

  for (vector<float>::const_iterator icurv = mCurvature.begin() + mingap;
                                     icurv < mCurvature.end() - mingap;
                                     icurv++) {
    if (*icurv > mincurv && *icurv > *(icurv-1) && *icurv > *(icurv+1)) {
      const int kpt = icurv - mCurvature.begin(),
                idx = (int) floor(ds * kpt);
      vector<int>::const_iterator ipt = mAllPoints.begin() + 3*kpt;

      // Save info for this candidate point
      peaksegs.push_back(idx);
      peakcurvs.push_back(*icurv);
      peakpts.push_back(ipt);

      cout << " " << *ipt << " " << *(ipt+1) << " " << *(ipt+2)
           << " (curv = " << *icurv << ")" << endl;
    }
  }

  *isfinal.begin() = true;
  *(isfinal.end() - 1) = true;

  // Go down the list of candidate points to pick out dominant points
  while (!peakcurvs.empty()) {
    // Find candidate point with maximum curvature
    const int ipeak = max_element(peakcurvs.begin(), peakcurvs.end()) - 
                      peakcurvs.begin(),
              iseg = peaksegs.at(ipeak);
    const bool isLfinal = isfinal.at(iseg),
               isRfinal = isfinal.at(iseg+1);
    vector<int>::iterator imatch;
    vector<int>::const_iterator ipt  = peakpts.at(ipeak),
                                iptL = dompts.at(iseg),
                                iptR = dompts.at(iseg+1);
    const int distL = (ipt - iptL) / 3,
              distR = (iptR - ipt) / 3;

    // Find which end of the segment the candidate point is closest to
    if (distL < distR && !isLfinal && !(isRfinal && distR <= mingap)) {
      // Move left end of segment to this dominant point
      dompts.at(iseg) = ipt;
      isfinal.at(iseg) = true;
    }
    else if (!isRfinal && !(isLfinal && distL <= mingap)) {
      // Move right end of segment to this dominant point
      dompts.at(iseg+1) = ipt;
      isfinal.at(iseg+1) = true;
    }

    // Remove point from list of candidate points
    peaksegs.erase(peaksegs.begin() + ipeak);
    peakcurvs.erase(peakcurvs.begin() + ipeak);
    peakpts.erase(peakpts.begin() + ipeak);

    // Also remove any other candidate points that are in the same segment
    imatch = find(peaksegs.begin(), peaksegs.end(), iseg);

    while (imatch != peaksegs.end()) {
      int irem = imatch - peaksegs.begin(); 

      peaksegs.erase(peaksegs.begin() + irem);
      peakcurvs.erase(peakcurvs.begin() + irem);
      peakpts.erase(peakpts.begin() + irem);

      imatch = find(peaksegs.begin(), peaksegs.end(), iseg);
    }
  }

  cout << "INFO: Intermediate dominant points are" << endl;
  for (vector< vector<int>::const_iterator>::iterator idom = dompts.begin();
                                                      idom < dompts.end();
                                                      idom++)
    cout << " " << *(*idom) << " " << *(*idom+1) << " " << *(*idom+2) << endl;

  // Find segment ends that have not been moved yet
  for (vector<bool>::const_iterator ifinal = isfinal.begin();
                                    ifinal < isfinal.end(); ifinal++)
    if (! *ifinal) {
      const int kdom = ifinal - isfinal.begin();
      int lenL, lenR;
      float curvL, curvR,
            dshape, dshapemin = numeric_limits<double>::infinity();
      vector<int>::const_iterator iptL = dompts.at(kdom-1),
                                  iptR = dompts.at(kdom+1);
      vector<float>::const_iterator icurv = mCurvature.begin() +
                                            (iptL - mAllPoints.begin()) / 3;

      // Find point in the segment between the two neighboring dominant points
      // to split the segment into 2 subsegments with minimum difference in
      // "shape index" (to balance the subsegments' total curvature and length)
      lenL = mingap;
      lenR = (iptR-iptL)/3 - mingap;

      curvL = accumulate(icurv + 1, icurv + mingap - 1, 0.0);
      curvR = accumulate(icurv + mingap + 1, icurv + (iptR-iptL)/3, 0.0);

      iptL += 3*mingap;
      iptR -= 3*mingap;
      icurv += mingap;

      for (vector<int>::const_iterator ipt = iptL; ipt <= iptR; ipt += 3) {
        dshape = fabs( (lenL - lenR) / (float) lentot +
                       (curvL - curvR) / curvtot );

        if (dshape < dshapemin) {
          dompts.at(kdom) = ipt;
          dshape = dshapemin;
        }

        lenL++;
        lenR--;

        curvL += *icurv;
        curvR -= *(icurv+1);

        icurv++;
      }
    }

  cout << "INFO: Final dominant points are" << endl;
  for (vector< vector<int>::const_iterator>::iterator idom = dompts.begin();
                                                      idom < dompts.end();
                                                      idom++)
    cout << " " << *(*idom) << " " << *(*idom+1) << " " << *(*idom+2) << endl;

  // Parameterize the curve based on these dominant points
  for (int kdom = 1; kdom < mNumControl; kdom++) {
    const int seglen = (int) (dompts.at(kdom) - dompts.at(kdom-1)) / 3;
    const float ds = 1.0 / seglen;
    float s = kdom - 1.0;

    for (int kpt = seglen; kpt > 0; kpt--) {
      *iarc = s;
      s += ds;
      iarc++;
    }
  }

  *iarc = (float) (mNumControl - 1);

  CatmullRomFit(InputPoints, arcparam);

  success = InterpolateSpline();

  // If spline interpolation after fitting to these dominant points fails,
  // default to equidistant dominant points
  if (!success) {
    cout << "WARN: Defaulting to equidistant dominant points" << endl;

    CatmullRomFit(InputPoints);

    success = InterpolateSpline();
  }

  return success;
}

//
// Find which spline segment a point belongs to
//
unsigned int Spline::PointToSegment(unsigned int PointIndex) {
  unsigned int iseg = 0;

  if (PointIndex > mArcLength.size()-1) {
    cout << "ERROR: Specified point index (" << PointIndex
         << ") is outside acceptable range (0-" << mArcLength.size()-1
         << ")" << endl;
    exit(1);
  }

  for (vector<float>::const_iterator iarc = mArcLength.begin() + PointIndex;
                                     iarc > mArcLength.begin(); iarc--)
    if (*iarc < *(iarc-1))
      iseg++;

  return iseg;
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
    const unsigned int diffstep = 3;
    vector<float> smoothpts(mAllPoints.size());

    mFiniteDifference1.resize(mAllPoints.size());

    // Smooth discrete point coordinates
    CurveSmooth(smoothpts, mAllPoints);

    // Approximate first derivative by smoothed finite differences 
    // of the point coordinates
    CurveFiniteDifferences(mFiniteDifference1, smoothpts, diffstep);

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
    const unsigned int diffstep = 3;

    mFiniteDifference2.resize(mFiniteDifference1.size());

    // Approximate second derivative by smoothed finite differences 
    // of the first derivative
    CurveFiniteDifferences(mFiniteDifference2, mFiniteDifference1, diffstep);

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
void Spline::ReadControlPoints(const string ControlPointFile) {
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
void Spline::ReadMask(const string MaskFile) {
  if (mMask)   MRIfree(&mMask);
  if (mVolume) MRIfree(&mVolume);

  cout << "Loading spline mask from " << MaskFile << endl;
  mMask = MRIread(MaskFile.c_str());
  mVolume = MRIclone(mMask, NULL);
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
// Set mask volume
//
void Spline::SetMask(MRI *Mask) {
  if (mMask)   MRIfree(&mMask);
  if (mVolume) MRIfree(&mVolume);

  mMask = MRIcopy(Mask, NULL);
  mVolume = MRIclone(mMask, NULL);
}

//
// Write spline to volume
//
void Spline::WriteVolume(const string VolumeFile, const bool ShowControls) {
  if (ShowControls)
    for (vector<int>::const_iterator icpt = mControlPoints.begin();
                                     icpt != mControlPoints.end(); icpt += 3)
      MRIsetVoxVal(mVolume, icpt[0], icpt[1], icpt[2], 0, 2);

  cout << "Writing spline volume to " << VolumeFile << endl;
  MRIwrite(mVolume, VolumeFile.c_str());
}

//
// Write spline coordinates to file
//
void Spline::WriteAllPoints(const string TextFile) {
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
void Spline::WriteTangent(const string TextFile) {
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
void Spline::WriteNormal(const string TextFile) {
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
void Spline::WriteCurvature(const string TextFile) {
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
void Spline::WriteValues(vector<MRI *> &ValueVolumes, const string TextFile) {
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
// Least-squares fit of a Catmull-Rom spline to a set points,
// given the coordinates of the points but no user-specified parameterization
// (default to parameterizing the curve based on equidistant dominant points)
//
void Spline::CatmullRomFit(const vector<int> &InputPoints) {
  const int npts = (int) InputPoints.size() / 3;
  const float smax = (float) (mNumControl - 1),
              ds = smax / (npts - 1);
  float s = 0.0;
  vector<float> arcparam(npts);

  for (vector<float>::iterator iarc = arcparam.begin(); iarc < arcparam.end();
                                                        iarc++) {
    *iarc = s;
    s += ds;
  }

  CatmullRomFit(InputPoints, arcparam);
}

//
// Least-squares fit of a Catmull-Rom spline to a set of points,
// given the coordinates of the points and the respective arc length
// parameter values, which should be reals in the interval [0, mNumControl-1]
//
void Spline::CatmullRomFit(const vector<int> &InputPoints,
                           const vector<float> &ArcLengthParameter) {
  const int npts = (int) InputPoints.size() / 3;
  vector<int>::const_iterator ipt = InputPoints.begin() + 3;
  vector<int>::iterator icpt;
  vector<float>::const_iterator iarc = ArcLengthParameter.begin() + 1;
  MATRIX *A = NULL, *Ap = NULL, *y = NULL, *x = NULL;

  // Fit all but the first and last control points,
  // which will be set equal to the first and last of the input points
  A  = MatrixAlloc(npts-2, mNumControl-2, MATRIX_REAL);
  Ap = MatrixAlloc(mNumControl-2, npts-2, MATRIX_REAL);
  y  = MatrixAlloc(npts-2, 3, MATRIX_REAL);
  x  = MatrixAlloc(mNumControl-2, 3, MATRIX_REAL);

  for (int irow = 1; irow < npts-1; irow++) {
    const int idx = (int) floor(*iarc);
    const float t = *iarc - floor(*iarc),
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

    iarc++;
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

