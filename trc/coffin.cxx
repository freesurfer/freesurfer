/**
 * @brief Container of tractography data and methods
 *
 * Container of tractography data and methods
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


#include <coffin.h>

#include <sstream>
#include <iomanip>

using namespace std;

const unsigned int Aeon::mDiffStep = 3;
int Aeon::mMaxAPosterioriPath;
unsigned int Aeon::mMaxAPosterioriPath0;
vector<float> Aeon::mPriorSamples;
vector< vector<int> > Aeon::mBasePathPointSamples;
MRI *Aeon::mBaseMask;

const unsigned int Coffin::mMaxTryMask = 100,
                   Coffin::mMaxTryWhite = 10,
                   Coffin::mDiffStep = 3;
const float Coffin::mTangentBinSize = 1/3.0,	// 0.1,
            Coffin::mCurvatureBinSize = 0.01;	// 0.002;

//
// A single point in time
//
Aeon::Aeon() {
  mNx = mNy = mNz = mNxy = mNumVox = 0;
  mMask = 0;
  ClearPath();
}

Aeon::~Aeon() {
}

//
// Set the common base template mask for all time points
//
void Aeon::SetBaseMask(MRI *BaseMask) { mBaseMask = BaseMask; }

//
// Save a sample of the atlas-based path priors, common among all time points
//
void Aeon::SavePathPriors(vector<float> &Priors) {
  mPriorSamples.insert(mPriorSamples.end(), Priors.begin(), Priors.end());
}

//
// Save a path sample in base space, common among all time points
//
void Aeon::SaveBasePath(vector<int> &PathPoints) {
  mBasePathPointSamples.push_back(PathPoints);
}

//
// Set a path sample as the MAP path
//
void Aeon::SetPathMap(unsigned int PathIndex) {
  mMaxAPosterioriPath0 = PathIndex;
}

//
// Read data specific to a single time point
//
void Aeon::ReadData(const std::string RootDir, const std::string DwiFile,
                    const std::string GradientFile, const std::string BvalueFile,
                    const std::string MaskFile, const std::string BedpostDir,
                    const int NumTract, const float FminPath,
                    const std::string BaseXfmFile) {
  string dwifile, gradfile, bvalfile, maskfile, bpdir;
  std::string fname;
  MRI *dwi, *phi[NumTract], *theta[NumTract], *f[NumTract],
      *v0[NumTract], *f0[NumTract], *d0;

  if (!RootDir.empty()) {
    mRootDir = RootDir + "/";
  }

  dwifile    = mRootDir + DwiFile;
  gradfile   = mRootDir + GradientFile;
  bvalfile   = mRootDir + BvalueFile;
  maskfile   = mRootDir + MaskFile;
  bpdir      = mRootDir + BedpostDir;

  // Read diffusion-weighted images
  cout << "Loading DWIs from " << dwifile << endl;
  dwi = MRIread(dwifile.c_str());
  if (!dwi) {
    cout << "ERROR: Could not read " << dwifile << endl;
    exit(1);
  }

  // Size of diffusion-weighted images
  mNx = dwi->width;
  mNy = dwi->height;
  mNz = dwi->depth;
  mNxy = mNx * mNy;

  // Read mask
  cout << "Loading mask from " << maskfile << endl;
  mMask = MRIread(maskfile.c_str());
  if (!mMask) {
    cout << "ERROR: Could not read " << maskfile << endl;
    exit(1);
  }

  // Read parameter samples from BEDPOST directory
  cout << "Loading BEDPOST parameter samples from " << bpdir << endl;
  for (int itract = 0; itract < NumTract; itract++) {
    fname = bpdir + "/merged_ph" + std::to_string(itract+1) + "samples.nii.gz";
    phi[itract] = MRIread(fname.c_str());
    if (!phi[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
    fname = bpdir + "/merged_th" + std::to_string(itract+1) + "samples.nii.gz";
    theta[itract] = MRIread(fname.c_str());
    if (!theta[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
    fname = bpdir + "/merged_f" + std::to_string(itract+1) + "samples.nii.gz";
    f[itract] = MRIread(fname.c_str());
    if (!f[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
    fname = bpdir + "/dyads" + std::to_string(itract+1) + ".nii.gz";
    v0[itract] = MRIread(fname.c_str());
    if (!v0[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
    fname = bpdir + "/mean_f" + std::to_string(itract+1) + "samples.nii.gz";
    f0[itract] = MRIread(fname.c_str());
    if (!f0[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
  }

  fname = bpdir + "/mean_dsamples.nii.gz";
  d0 = MRIread(fname.c_str());
  if (!d0) {
    cout << "ERROR: Could not read " << fname << endl;
    exit(1);
  }

  // Initialize voxel-wise diffusion model
  Bite::SetStatic(gradfile.c_str(), bvalfile.c_str(),
                  NumTract, phi[0]->nframes, FminPath);

  if (Bite::GetNumDir() != dwi->nframes) {
    cout << "ERROR: Dimensions of " << bvalfile << " and " << dwifile
         << " do not match" << endl;
    exit(1);
  }

  cout << "INFO: Found "
       << Bite::GetNumB0() << " baseline images (b = "
       << Bite::GetLowBvalue() << ") out of a total of "
       << Bite::GetNumDir() << " frames" << endl;

  mData.clear();
  for (int iz = 0; iz < mNz; iz++)
    for (int iy = 0; iy < mNy; iy++)
      for (int ix = 0; ix < mNx; ix++)
        if (MRIgetVoxVal(mMask, ix, iy, iz, 0)) {
          Bite data = Bite(dwi, phi, theta, f, v0, f0, d0, ix, iy, iz);
          mData.push_back(data);
        }

  mDataMask.clear();
  mNumVox = 0;
  for (int iz = 0; iz < mNz; iz++)
    for (int iy = 0; iy < mNy; iy++)
      for (int ix = 0; ix < mNx; ix++)
        if (MRIgetVoxVal(mMask, ix, iy, iz, 0)) {
          mDataMask.push_back(&mData[mNumVox]);
          mNumVox++;
        }
        else
          mDataMask.push_back(0);

  cout << "INFO: Found " << mNumVox << " voxels in brain mask" << endl;

  // Free temporary variables
  MRIfree(&dwi);

  for (int itract = 0; itract < NumTract; itract++) {
    MRIfree(&phi[itract]);
    MRIfree(&theta[itract]);
    MRIfree(&f[itract]);
    MRIfree(&v0[itract]);
    MRIfree(&f0[itract]);
  }

  MRIfree(&d0);

  // Read transform from base template space to native DWI space
  // (only used for longitudinal data)
  if (!BaseXfmFile.empty()) {
    string regfile = mRootDir + BaseXfmFile;
    if (!mBaseMask) {
      cout << "ERROR: Cannot load base-to-DWI transform without base mask"
           << endl;
      exit(1);
    }
    mBaseReg.ReadXfm(regfile.c_str(), mBaseMask, mMask);
  }
}

//
// Return a pointer to this time point's mask
//
MRI *Aeon::GetMask() const { return mMask; }

//
// Return a pointer to this time point's base template mask
// (null if running in cross-sectional mode)
//
MRI *Aeon::GetBaseMask() const { return mBaseMask; }

//
// Return this time point's spatial resolution
//
float Aeon::GetDx() const { return mMask->xsize; }
float Aeon::GetDy() const { return mMask->ysize; }
float Aeon::GetDz() const { return mMask->zsize; }

//
// Return this time point's current number of saved path samples
//
unsigned int Aeon::GetNumSample() const { return mPathPointSamples.size(); }

//
// Free the space allocation to this time point's mask for clean-up
//
void Aeon::FreeMask() {
  MRIfree(&mMask);
}

//
// Set this time point's output directory for the current pathway
//
void Aeon::SetOutputDir(const std::string OutDir) {
  string cmdline("mkdir -p ");

  mOutDir = mRootDir + OutDir;
  cmdline += mOutDir; 

  cout << "Creating output directory " << mOutDir << endl;
  if (system(cmdline.c_str()) != 0) {
    cout << "ERROR: Could not create directory " << mOutDir << endl;
    exit(1);
  }
}

//
// Return this time point's output directory for the current pathway
//
const string &Aeon::GetOutputDir() const { return mOutDir; }

//
// Clear all path-related variables
//
void Aeon::ClearPath() {
  // Path-related variables that are common among all time points
  mMaxAPosterioriPath = -1;
  mMaxAPosterioriPath0 = 0;
  mPriorSamples.clear();
  mBasePathPointSamples.clear();

  // Path-related variables that are specific to this time point
  mPathPoints.clear();
  mPathPointsNew.clear();
  mPathPointSamples.clear();
  mPathPhi.clear();
  mPathPhiNew.clear();
  mPathTheta.clear();
  mPathThetaNew.clear();

  mDataFitSamples.clear();

  mRejectF = false;
  mAcceptF = false;
  mRejectTheta = false;
  mAcceptTheta = false;
  mLog.clear();
  mErrorPoint.clear();

  mPathLength = 0;
  mPathLengthNew = 0;
  mLikelihoodOnPath = 0;
  mLikelihoodOnPathNew = 0;
  mPriorOnPath = 0;
  mPriorOnPathNew = 0;
  mPosteriorOnPath = 0;
  mPosteriorOnPathNew = 0;
  mLikelihoodOffPath = 0;
  mLikelihoodOffPathNew = 0;
  mPriorOffPath = 0;
  mPriorOffPathNew = 0;
  mPosteriorOffPath = 0;
  mPosteriorOffPathNew = 0;
}

//
// Map proposed path from the base space to this time point's native space
// and calculate orientation angles along the path in the native space
//
bool Aeon::MapPathFromBase(Spline &BaseSpline) {
  vector<float>::const_iterator tangbegin, tangend;
  vector<float>::iterator iphi, itheta;
  vector<float> diff1;

  if (mBaseReg.IsEmpty()) {	// Single time point, there is no base
    // Copy spline points
    mPathPointsNew.resize(BaseSpline.GetAllPointsEnd() -
                          BaseSpline.GetAllPointsBegin());
    copy(BaseSpline.GetAllPointsBegin(), BaseSpline.GetAllPointsEnd(), 
         mPathPointsNew.begin());

    // Calculate tangent vectors along spline
    BaseSpline.ComputeTangent();

    tangbegin = BaseSpline.GetTangentBegin();
    tangend   = BaseSpline.GetTangentEnd();
  }
  else {			// Multiple time points, must map path from base
    vector<int> point(3);
    vector<float> pointf(3), pathsmooth;

    // Map spline points from base to native DWI space, making sure there are
    // no duplicate points due to the higher resolution of the base space
    mPathPointsNew.clear();
    mErrorPoint.clear();

    for (vector<int>::const_iterator iptbase = BaseSpline.GetAllPointsBegin();
                                     iptbase < BaseSpline.GetAllPointsEnd();
                                     iptbase += 3) {

      for (int k = 0; k < 3; k++)
        pointf[k] = (float) iptbase[k];

      mBaseReg.ApplyXfm(pointf, pointf.begin());

      for (int k = 0; k < 3; k++)
        point[k] = (int) round(pointf[k]);

      // Keep point if it's in the mask (a point in the base mask may not be in 
      // an individual time point's mask due to interpolation/rounding errors)
      if (!IsInMask(point.begin())) {
        mErrorPoint.insert(mErrorPoint.begin(), point.begin(), point.end());
        return false;
      }
      else
        mPathPointsNew.insert(mPathPointsNew.end(), point.begin(), point.end());
    }

    // Smooth discrete point coordinates
    pathsmooth.resize(mPathPointsNew.size());
    CurveSmooth(pathsmooth, mPathPointsNew);

    // Approximate first derivative by smoothed finite differences
    // of the point coordinates
    diff1.resize(mPathPointsNew.size());
    CurveFiniteDifferences(diff1, pathsmooth, mDiffStep);

    tangbegin = diff1.begin();
    tangend   = diff1.end();
  }

  // Find path length in the native space
  mPathLengthNew = mPathPointsNew.size() / 3;

  // Calculate orientation angles along path
  mPathPhiNew.resize(mPathPointsNew.size());
  iphi = mPathPhiNew.begin();
  mPathThetaNew.resize(mPathPointsNew.size());
  itheta = mPathThetaNew.begin();

  for (vector<float>::const_iterator itang = tangbegin; itang < tangend;
                                                        itang += 3) {
    *iphi = atan2(itang[1], itang[0]);
    *itheta = acos(itang[2] / sqrt(itang[0]*itang[0] + itang[1]*itang[1]
                                                       + itang[2]*itang[2]));

    iphi++;
    itheta++;
  }

  return true;
}

//
// Find duplicate consecutive points along the current path
//
void Aeon::FindDuplicatePathPoints(vector<bool> &IsDuplicate) {
  vector<bool>::iterator idup = IsDuplicate.begin();
  vector<int>::const_iterator ipt = mPathPoints.begin();

  while (ipt < mPathPoints.end()) {
    vector<int>::const_iterator iptnext = ipt + 3;

    idup++;

    while (iptnext < mPathPoints.end() && *idup) {
      iptnext += 3;
      idup++;
    }

    if (iptnext < mPathPoints.end() &&
        ipt[0] == iptnext[0] && ipt[1] == iptnext[1] && ipt[2] == iptnext[2]) 
      *idup = true;

    ipt = iptnext;

    while (ipt < mPathPoints.end() && *idup) {
      ipt += 3;
      idup++;
    }
  }
}

//
// Remove specified points from current path
//
void Aeon::RemovePathPoints(vector<bool> &DoRemove, unsigned int NewSize) {
  vector<int>::const_iterator ipt;
  vector<int>::iterator iptnew;
  vector<int> newpath;

  if (NewSize > 0)
    newpath.resize(NewSize);
  else {			// If new size was not pre-computed, compute it
    unsigned int newsize = 0;

    for (vector<bool>::const_iterator irem = DoRemove.begin();
                                      irem < DoRemove.end(); irem++)
      if (! *irem)
        newsize += 3;

    newpath.resize(newsize);
  }

  ipt = mPathPoints.begin();
  iptnew = newpath.begin();

  for (vector<bool>::const_iterator irem = DoRemove.begin();
                                    irem < DoRemove.end(); irem++) {
    if (! *irem) {
      copy(ipt, ipt+3, iptnew);
      iptnew += 3;
    }

    ipt += 3;
  }

  mPathPoints.resize(newpath.size());
  copy(newpath.begin(), newpath.end(), mPathPoints.begin());

  mPathLength = mPathPoints.size() / 3;
}

//
// Propose diffusion parameters by sampling from their marginal posteriors
// for this time point along the proposed and current path
//
void Aeon::ProposeDiffusionParameters() {
  vector<int>::const_iterator ipt;

  // Sample parameters on proposed path
  for (ipt = mPathPointsNew.begin(); ipt < mPathPointsNew.end(); ipt += 3) {
    Bite *ivox = mDataMask[ipt[0] + ipt[1]*mNx + ipt[2]*mNxy];
    ivox->SampleParameters();
  }

  // Sample parameters on current path
  for (ipt = mPathPoints.begin(); ipt < mPathPoints.end(); ipt += 3) {
    Bite *ivox = mDataMask[ipt[0] + ipt[1]*mNx + ipt[2]*mNxy];
    ivox->SampleParameters();
  }
}

//
// Compute data-fit terms of the objective function for this time point
// along the proposed and current path
//
bool Aeon::ComputePathDataFit() {
  vector<float>::const_iterator iphi, itheta;

  mRejectF = false;
  mAcceptF = false;
  mRejectTheta = false;
  mAcceptTheta = false;
  mLog.clear();
  mErrorPoint.clear();

  // Compute data-fit terms on proposed path
  mLikelihoodOnPathNew = 0;
  mPriorOnPathNew = 0;
  mLikelihoodOffPathNew = 0;
  mPriorOffPathNew = 0;
  iphi = mPathPhiNew.begin();
  itheta = mPathThetaNew.begin();

  for (vector<int>::iterator ipt = mPathPointsNew.begin();
                             ipt < mPathPointsNew.end(); ipt += 3) {
    Bite *ivox = mDataMask[ipt[0] + ipt[1]*mNx + ipt[2]*mNxy];

    ivox->ComputeLikelihoodOffPath();
    ivox->ComputeLikelihoodOnPath(*iphi, *itheta);
    if (ivox->IsFZero()) {
      ostringstream msg;
      msg << "Reject due to f=0 at "
          << ipt[0] << " " << ipt[1] << " " << ipt[2];
      mLog = msg.str();
      mErrorPoint.insert(mErrorPoint.begin(), ipt, ipt+3);

      mRejectF = true;

      return false;
    }
    if (ivox->IsThetaZero()) {
      ostringstream msg;
      msg << "Accept due to theta=0 at "
          << ipt[0] << " " << ipt[1] << " " << ipt[2];
      mLog = msg.str();
      mErrorPoint.insert(mErrorPoint.begin(), ipt, ipt+3);

      mAcceptTheta = true;

      return false;
    }
    ivox->ComputePriorOffPath();
    ivox->ComputePriorOnPath();

    mLikelihoodOnPathNew += ivox->GetLikelihoodOnPath();
    mPriorOnPathNew += ivox->GetPriorOnPath();

    mLikelihoodOffPathNew += ivox->GetLikelihoodOffPath();
    mPriorOffPathNew += ivox->GetPriorOffPath();

    iphi++;
    itheta++;
  }

  mPosteriorOnPathNew  = mLikelihoodOnPathNew  + mPriorOnPathNew;
  mPosteriorOffPathNew = mLikelihoodOffPathNew + mPriorOffPathNew;

  // Compute data-fit terms on current path
  mLikelihoodOnPath = 0;
  mPriorOnPath = 0;
  mLikelihoodOffPath = 0;
  mPriorOffPath = 0;
  iphi = mPathPhi.begin();
  itheta = mPathTheta.begin();

  for (vector<int>::iterator ipt = mPathPoints.begin();
                             ipt < mPathPoints.end(); ipt += 3) {
    Bite *ivox = mDataMask[ipt[0] + ipt[1]*mNx + ipt[2]*mNxy];

    ivox->ComputeLikelihoodOffPath();
    ivox->ComputeLikelihoodOnPath(*iphi, *itheta);
    if (ivox->IsFZero()) {
      ostringstream msg;
      msg << "Accept due to f=0 at "
          << ipt[0] << " " << ipt[1] << " " << ipt[2];
      mLog = msg.str();
      mErrorPoint.insert(mErrorPoint.begin(), ipt, ipt+3);

      mAcceptF = true;

      return false;
    }
    if (ivox->IsThetaZero()) {
      ostringstream msg;
      msg << "Reject due to theta=0 at "
          << ipt[0] << " " << ipt[1] << " " << ipt[2];
      mLog = msg.str();
      mErrorPoint.insert(mErrorPoint.begin(), ipt, ipt+3);

      mRejectTheta = true;

      return false;
    }
    ivox->ComputePriorOffPath();
    ivox->ComputePriorOnPath();

    mLikelihoodOnPath += ivox->GetLikelihoodOnPath();
    mPriorOnPath += ivox->GetPriorOnPath();

    mLikelihoodOffPath += ivox->GetLikelihoodOffPath();
    mPriorOffPath += ivox->GetPriorOffPath();

    iphi++;
    itheta++;
  }

  mPosteriorOnPath  = mLikelihoodOnPath  + mPriorOnPath;
  mPosteriorOffPath = mLikelihoodOffPath + mPriorOffPath;

  return true;
}

//
// Find spline segment corresponding to the point along the path
// where an error occured
//
int Aeon::FindErrorSegment(Spline &BaseSpline) {
  int iseg = -1;

  if (mErrorPoint.empty())
    return iseg;

  // Find where along the spline the error occured
  if (mBaseReg.IsEmpty()) {	// Single time point, there is no base
    for (vector<int>::const_iterator iptbase = BaseSpline.GetAllPointsBegin();
                                     iptbase < BaseSpline.GetAllPointsEnd();
                                     iptbase += 3)
      if (equal(mErrorPoint.begin(), mErrorPoint.end(), iptbase)) {
        const unsigned int ierr = (iptbase - BaseSpline.GetAllPointsBegin())/3; 
        iseg = (int) BaseSpline.PointToSegment(ierr);
        break;
      }
  }
  else {			// Multiple time points, must map path from base
    vector<int> point(3);
    vector<float> pointf(3);

    for (vector<int>::const_iterator iptbase = BaseSpline.GetAllPointsBegin();
                                     iptbase < BaseSpline.GetAllPointsEnd();
                                     iptbase += 3) {
      for (int k = 0; k < 3; k++)
        pointf[k] = (float) iptbase[k];

      mBaseReg.ApplyXfm(pointf, pointf.begin());

      for (int k = 0; k < 3; k++)
        point[k] = (int) round(pointf[k]);

      if (equal(mErrorPoint.begin(), mErrorPoint.end(), point.begin())) {
        const unsigned int ierr = (iptbase - BaseSpline.GetAllPointsBegin())/3; 
        iseg = (int) BaseSpline.PointToSegment(ierr);
        break;
      }
    }
  }

  return iseg;
}

//
// Copy newly accepted path over current path for this time point
//
void Aeon::UpdatePath() {
  mPathPoints.resize(mPathPointsNew.size());
  copy(mPathPointsNew.begin(), mPathPointsNew.end(), mPathPoints.begin());

  mPathPhi.resize(mPathPhiNew.size());
  copy(mPathPhiNew.begin(), mPathPhiNew.end(), mPathPhi.begin());

  mPathTheta.resize(mPathThetaNew.size());
  copy(mPathThetaNew.begin(), mPathThetaNew.end(), mPathTheta.begin());

  mPathLength        = mPathLengthNew;
  mLikelihoodOnPath  = mLikelihoodOnPathNew; 
  mLikelihoodOffPath = mLikelihoodOffPathNew; 
  mPriorOnPath       = mPriorOnPathNew; 
  mPriorOffPath      = mPriorOffPathNew; 
  mPosteriorOnPath   = mPosteriorOnPathNew; 
  mPosteriorOffPath  = mPosteriorOffPathNew; 
}

//
// Save data-fit terms of accepted and rejected path for this time point
//
void Aeon::SavePathDataFit(bool IsPathAccepted) {
  if (IsPathAccepted) {			// The newly sampled path was accepted
    mDataFitSamples.push_back(GetDataFitNew());
    mDataFitSamples.push_back(GetDataFit());
  }
  else {				// The newly sampled path was rejected
    mDataFitSamples.push_back(GetDataFit());
    mDataFitSamples.push_back(GetDataFitNew());
  }
}

//
// Save current path as an MCMC sample for this time point
//
void Aeon::SavePath() {
  mPathPointSamples.push_back(mPathPoints);
}

//
// Write output files for this time point
//
void Aeon::WriteOutputs() {
  char outorient[4];
  MATRIX *outv2r;
  CTrackWriter trkwriter;
  TRACK_HEADER trkheadout;
  vector<int> lengths(mPathPointSamples.size()), cptsmap, cptsmap0;
  vector<int>::iterator ilen;
  vector<int>::const_iterator iptbase;
  vector<float> fpath;
  vector<float>::iterator ifpt;
  vector<float>::const_iterator ipr;
  vector< vector<int> >::const_iterator pathmap, basepathmap;
  string fname;
  ofstream mapfile;
  MRI *pdvol;

  // Find maximum a posteriori path, if it hasn't been found yet:
  // Case where this is the first of multiple time points
  if (!mBasePathPointSamples.empty() && mMaxAPosterioriPath < 0) {
    pdvol = MRIclone(mBaseMask, NULL);

    ComputePathHisto(pdvol, mBasePathPointSamples);
    ComputePathLengths(lengths, mBasePathPointSamples);

    mMaxAPosterioriPath = FindMaxAPosterioriPath(mBasePathPointSamples, 
                                                 lengths, pdvol);

    MRIfree(&pdvol);
    fill(lengths.begin(), lengths.end(), 0);
  }

  // Save a .trk file with path samples:
  // Output space orientation information
  outv2r = MRIgetVoxelToRasXform(mMask);
  MRIdircosToOrientationString(mMask, outorient);

  // Set output .trk header
  trkheadout.Initialize();

  trkheadout.origin[0] = 0;
  trkheadout.origin[1] = 0;
  trkheadout.origin[2] = 0;

  trkheadout.voxel_size[0] = mMask->xsize;
  trkheadout.voxel_size[1] = mMask->ysize;
  trkheadout.voxel_size[2] = mMask->zsize;

  trkheadout.dim[0] = mMask->width;
  trkheadout.dim[1] = mMask->height;
  trkheadout.dim[2] = mMask->depth;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      trkheadout.vox_to_ras[i][j] = outv2r->rptr[i+1][j+1];

  strcpy(trkheadout.voxel_order, outorient);

  // Find patient-to-scanner coordinate transform:
  // Take x and y vectors from vox2RAS matrix, convert to LPS,
  // divide by voxel size
  trkheadout.image_orientation_patient[0] =
    - trkheadout.vox_to_ras[0][0] / trkheadout.voxel_size[0];
  trkheadout.image_orientation_patient[1] =
    - trkheadout.vox_to_ras[1][0] / trkheadout.voxel_size[0];
  trkheadout.image_orientation_patient[2] =
      trkheadout.vox_to_ras[2][0] / trkheadout.voxel_size[0];
  trkheadout.image_orientation_patient[3] =
    - trkheadout.vox_to_ras[0][1] / trkheadout.voxel_size[1];
  trkheadout.image_orientation_patient[4] =
    - trkheadout.vox_to_ras[1][1] / trkheadout.voxel_size[1];
  trkheadout.image_orientation_patient[5] =
      trkheadout.vox_to_ras[2][1] / trkheadout.voxel_size[1];

  trkheadout.n_count = (int) mPathPointSamples.size();

  // Open output .trk file
  fname = mOutDir + "/path.pd.trk";
  if (!trkwriter.Initialize(fname.c_str(), trkheadout)) {
    cout << "ERROR: Cannot open output file " << fname << endl;
    cout << "ERROR: " << trkwriter.GetLastErrorMessage() << endl;
    exit(1);
  }

  for (vector< vector<int> >::const_iterator ipath = mPathPointSamples.begin();
                                             ipath < mPathPointSamples.end();
                                             ipath++) {
    fpath.resize(ipath->size());
    ifpt = fpath.begin();

    // Make point coordinates .5-based and multiply by voxel size
    for (vector<int>::const_iterator ipt = ipath->begin();
                                     ipt < ipath->end(); ipt += 3) {
      for (int k = 0; k < 3; k++)
        ifpt[k] = (ipt[k] + .5) * trkheadout.voxel_size[k];

      ifpt += 3;
    }

    // Write path to .trk file
    trkwriter.WriteNextTrack(fpath.size()/3, &fpath[0]);
  }

  // Close output .trk file
  trkwriter.Close();
  MatrixFree(&outv2r);

  // Save volume of path samples
  pdvol = MRIclone(mMask, NULL);

  ComputePathHisto(pdvol, mPathPointSamples);

  fname = mOutDir + "/path.pd.nii.gz";
  MRIwrite(pdvol, fname.c_str());

  // Save length of path samples
  ComputePathLengths(lengths, mPathPointSamples);

  fname = mOutDir + "/length.samples.txt";
  ofstream lenfile(fname.c_str(), ios::out);

  if (!lenfile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  for (vector<int>::const_iterator ilen = lengths.begin(); ilen < lengths.end();
                                                           ilen++)
    lenfile << *ilen << endl;

  // Find maximum a posteriori path, if it hasn't been found yet:
  // Case where this is the only time point
  if (mMaxAPosterioriPath < 0)
    mMaxAPosterioriPath = FindMaxAPosterioriPath(mPathPointSamples, 
                                                 lengths, pdvol);

  pathmap = mPathPointSamples.begin() + mMaxAPosterioriPath;

  if (!mBasePathPointSamples.empty()) {
    basepathmap = mBasePathPointSamples.begin() + mMaxAPosterioriPath;
    iptbase = basepathmap->begin();
  }

  // Save maximum a posteriori path coordinates
  fname = mOutDir + "/path.map.txt";
  mapfile.open(fname.c_str(), ios::out);

  if (!mapfile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  for (vector<int>::const_iterator ipt = pathmap->begin();
                                   ipt < pathmap->end(); ipt += 3) {
    mapfile << ipt[0] << " " << ipt[1] << " " << ipt[2];

    if (!mBasePathPointSamples.empty()) {
      mapfile << " " << iptbase[0] << " " << iptbase[1] << " " << iptbase[2];
      iptbase += 3;
    }

    mapfile << endl;
  }

  mapfile.close();

  // Save maximum a posteriori path as a volume
  MRIclear(pdvol);
  for (vector<int>::const_iterator ipt = pathmap->begin();
                                   ipt < pathmap->end(); ipt += 3) {
    const int ix = ipt[0], iy = ipt[1], iz = ipt[2];
    MRIsetVoxVal(pdvol, ix, iy, iz, 0, 1);
  }

  fname = mOutDir + "/path.map.nii.gz";
  MRIwrite(pdvol, fname.c_str());

if (0) {		// !OLD METHOD! Remove this eventually
  // Save maximum a posteriori path as a volume
  pathmap = mPathPointSamples.begin() + mMaxAPosterioriPath0;

  MRIclear(pdvol);
  for (vector<int>::const_iterator ipt = pathmap->begin();
                                   ipt < pathmap->end(); ipt += 3) {
    const int ix = ipt[0], iy = ipt[1], iz = ipt[2];
    MRIsetVoxVal(pdvol, ix, iy, iz, 0, 1);
  }

  fname = mOutDir + "/path.map.nii.gz";
  MRIwrite(pdvol, fname.c_str());
}

  // Save volumes of path start point samples
  MRIclear(pdvol);
  for (vector< vector<int> >::const_iterator ipath = mPathPointSamples.begin();
                                             ipath < mPathPointSamples.end();
                                             ipath++) {
    const int ix = *(ipath->begin()),
              iy = *(ipath->begin()+1),
              iz = *(ipath->begin()+2);

    MRIsetVoxVal(pdvol, ix, iy, iz, 0, MRIgetVoxVal(pdvol, ix, iy, iz, 0) + 1);
  }

  fname = mOutDir + "/endpt1.pd.nii.gz";
  MRIwrite(pdvol, fname.c_str());

  // Save volumes of path end point samples
  MRIclear(pdvol);
  for (vector< vector<int> >::const_iterator ipath = mPathPointSamples.begin();
                                             ipath < mPathPointSamples.end();
                                             ipath++) {
    const int ix = *(ipath->end()-3),
              iy = *(ipath->end()-2),
              iz = *(ipath->end()-1);

    MRIsetVoxVal(pdvol, ix, iy, iz, 0, MRIgetVoxVal(pdvol, ix, iy, iz, 0) + 1);
  }

  fname = mOutDir + "/endpt2.pd.nii.gz";
  MRIwrite(pdvol, fname.c_str());

  // Save data-fit and prior terms of objective function on accepted and
  // rejected path samples
  fname = mOutDir + "/pd.samples.txt";
  ofstream pdfile(fname.c_str(), ios::out);

  if (!pdfile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  pdfile << "DataFit1 XyzPrior1 AnatPrior1 ShapePrior1 "
         << "DataFit0 XyzPrior0 AnatPrior0 ShapePrior0" << endl;
  ipr = mPriorSamples.begin();
  for (vector<float>::const_iterator idf = mDataFitSamples.begin();
                                     idf < mDataFitSamples.end(); idf += 2) {
    pdfile << idf[0] << " " << ipr[0] << " " << ipr[1] << " " << ipr[2] << " "
           << idf[1] << " " << ipr[3] << " " << ipr[4] << " " << ipr[5] << endl;

    ipr += 6;
  }

  MRIfree(&pdvol);
}

//
// Find number of points with f=0 on proposed path
//
unsigned int Aeon::GetNumFZerosNew() const {
  unsigned int nzeros = 0;

  for (vector<int>::const_iterator ipt = mPathPointsNew.begin();
                                   ipt < mPathPointsNew.end(); ipt += 3) {
    Bite *ivox = mDataMask[ipt[0] + ipt[1]*mNx + ipt[2]*mNxy];

    if (ivox->IsAllFZero())
      nzeros++;
  }

  return nzeros;
}

//
// Find number of points with f=0 on current path
//
unsigned int Aeon::GetNumFZeros() const {
  unsigned int nzeros = 0;

  for (vector<int>::const_iterator ipt = mPathPoints.begin();
                                   ipt < mPathPoints.end(); ipt += 3) {
    Bite *ivox = mDataMask[ipt[0] + ipt[1]*mNx + ipt[2]*mNxy];

    if (ivox->IsAllFZero())
      nzeros++;
  }

  return nzeros;
}

//
// Return path-related variables
//
bool Aeon::RejectF() const { return mRejectF; }
bool Aeon::AcceptF() const { return mAcceptF; }

bool Aeon::RejectTheta() const { return mRejectTheta; }
bool Aeon::AcceptTheta() const { return mAcceptTheta; }

const string &Aeon::GetLog() const { return mLog; }

unsigned int Aeon::GetPathLengthNew() const { return mPathLengthNew; }
unsigned int Aeon::GetPathLength() const { return mPathLength; }

double Aeon::GetLikelihoodOffPathNew() const {
  return mLikelihoodOffPathNew/mPathLengthNew;
}
double Aeon::GetLikelihoodOffPath() const {
  return mLikelihoodOffPath/mPathLength;
}

double Aeon::GetLikelihoodOnPathNew() const {
  return mLikelihoodOnPathNew/mPathLengthNew;
}
double Aeon::GetLikelihoodOnPath() const {
  return mLikelihoodOnPath/mPathLength;
}

double Aeon::GetPriorOnPathNew() const {
  return mPriorOnPathNew/mPathLengthNew;
}
double Aeon::GetPriorOnPath() const {
  return mPriorOnPath/mPathLength;
}

double Aeon::GetPriorOffPathNew() const {
  return mPriorOffPathNew/mPathLengthNew;
}
double Aeon::GetPriorOffPath() const {
  return mPriorOffPath/mPathLength;
}

double Aeon::GetPosteriorOffPathNew() const {
  return mPosteriorOffPathNew/mPathLengthNew;
}
double Aeon::GetPosteriorOffPath() const {
  return mPosteriorOffPath/mPathLength;
}

double Aeon::GetPosteriorOnPathNew() const {
  return mPosteriorOnPathNew/mPathLengthNew;
}
double Aeon::GetPosteriorOnPath() const {
  return mPosteriorOnPath/mPathLength;
}

double Aeon::GetDataFitNew() const {
  return (mPosteriorOnPathNew - mPosteriorOffPathNew) / mPathLengthNew;
}
double Aeon::GetDataFit() const {
  return (mPosteriorOnPath - mPosteriorOffPath) / mPathLength;
}

//
// The main container
//
Coffin::Coffin(const std::string OutDir, vector<std::string> InDirList,
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
               const std::string TangPriorFile, const std::string CurvPriorFile,
               const std::string NeighPriorFile, const std::string NeighIdFile,
               const int NeighPriorSet, 
               const std::string LocalPriorFile, const std::string LocalIdFile,
               const int LocalPriorSet, 
               const vector<std::string> AsegList,
               const std::string AffineXfmFile, const std::string NonlinXfmFile,
               const int NumBurnIn, const int NumSample,
               const int KeepSampleNth, const int UpdatePropNth,
               const std::string PropStdFile,
               const bool Debug) :
               mDebug(Debug),
               mPriorSetLocal(LocalPriorSet), mPriorSetNear(NeighPriorSet),
               mMask(0), mRoi1(0), mRoi2(0),
               mXyzPrior0(0), mXyzPrior1(0) {
  vector<std::string >::const_iterator idir;
  MRI *atlasref;
  ostringstream infostr;

  // Save input info for logging
  if (!InDirList.empty()) {
    infostr << "Input directory: ";
    for (idir = InDirList.begin(); idir < InDirList.end(); idir++)
      infostr << " " << *idir;
    infostr << endl;
  }
  infostr  << "DWIs: " << DwiFile << endl
           << "Gradients: " << GradientFile << endl
           << "B-values: " << BvalueFile << endl
           << "Mask: " << MaskFile << endl
           << "BEDPOST directory: " << BedpostDir << endl
           << "Max number of tracts per voxel: " << NumTract << endl
           << "Tract volume fraction threshold: " << FminPath << endl;
  if (!BaseXfmFile.empty())
    infostr << "Base-to-DWI affine registration: " << BaseXfmFile << endl;
  if (!BaseMaskFile.empty())
    infostr << "Base mask: " << BaseMaskFile << endl;
  if (!AffineXfmFile.empty())
    infostr << (!BaseMaskFile.empty()?"Base":"DWI")
            << "-to-atlas affine registration: " << AffineXfmFile << endl;
  if (!NonlinXfmFile.empty())
    infostr << (!BaseMaskFile.empty()?"Base":"DWI")
            << "-to-atlas nonlinear registration: " << NonlinXfmFile << endl;
  if (!AsegList.empty()) {
    infostr << "Segmentation map: ";
    for (auto ifile = AsegList.begin(); ifile < AsegList.end(); ifile++) {
      infostr << " " << *ifile;
    }
    infostr << endl;
  }
  mInfoGeneral = infostr.str();

  // Read base template mask for longitudinal data
  if (!BaseMaskFile.empty()) {
    cout << "Loading base mask from " << BaseMaskFile << endl;
    mMask = MRIread(BaseMaskFile.c_str());
    if (!mMask) {
      cout << "ERROR: Could not read " << BaseMaskFile << endl;
      exit(1);
    }
    Aeon::SetBaseMask(mMask);
  }
  else		// No base needed when running in cross-sectional mode
    Aeon::SetBaseMask(0);

  // Read diffusion data, anatomical segmentation, and transform to atlas
  // for each time point
  if (InDirList.empty())
    InDirList.push_back(0);

  mDwi.resize(InDirList.size());
  idir = InDirList.begin();

  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++) {
    idwi->ReadData(*idir, DwiFile, GradientFile, BvalueFile,
                          MaskFile, BedpostDir, NumTract, FminPath,
                          BaseXfmFile);
    idir++;
  }

  // In cross-sectional case there is no base template mask,
  // so use the single time point as the base
  if (!mMask)
    mMask = mDwi[0].GetMask();

  // Size of base image
  mNx = mMask->width;
  mNy = mMask->height;
  mNz = mMask->depth;
  mNxy = mNx * mNy;

  // Resolution of DWI space relative to base space (used to determine how big
  // the control point perturbations should be in base space)
  mResolution.resize(3);
  mResolution[0] = mDwi[0].GetDx();
  mResolution[1] = mDwi[0].GetDy();
  mResolution[2] = mDwi[0].GetDz();

  for (vector<Aeon>::const_iterator idwi = mDwi.begin() + 1; idwi < mDwi.end();
                                                             idwi++) {
    const float dx = idwi->GetDx(), dy = idwi->GetDy(), dz = idwi->GetDz();

    if (dx > mResolution[0])	mResolution[0] = dx;
    if (dy > mResolution[1])	mResolution[1] = dy;
    if (dz > mResolution[2])	mResolution[2] = dz;
  }

  mResolution[0] /= mMask->xsize;
  mResolution[1] /= mMask->ysize;
  mResolution[2] /= mMask->zsize;

  cout << "INFO: Resolution of DWI space relative to base space is ("
       << mResolution[0] << ", " << mResolution[1] << ", " << mResolution[2]
       << ")" << endl;

  // Set mask for spline interpolation
  mSpline.SetMask(mMask);

  // Allocate space for saving and reusing voxel coordinates in atlas space
  mAtlasCoords.resize(mNxy*mNz);

  // Read start ROI as atlas-space reference volume
  // TODO: Use more general reference volume if ROI isn't specified
  if (!RoiFile1.empty()) {
    cout << "Loading atlas reference volume from " << RoiFile1 << endl;
    atlasref = MRIread(RoiFile1.c_str());
    if (!atlasref) {
      cout << "ERROR: Could not read " << RoiFile1 << endl;
      exit(1);
    }
  }

  // Read DWI-to-atlas registration
#ifndef NO_CVS_UP_IN_HERE
  if (!NonlinXfmFile.empty()) {
    mAffineReg.ReadXfm(AffineXfmFile.c_str(), mMask, 0);
    mNonlinReg.ReadXfm(NonlinXfmFile.c_str(), atlasref);
  } else {
#endif
    if (!AffineXfmFile.empty()) {
      mAffineReg.ReadXfm(AffineXfmFile.c_str(), mMask, atlasref);
    }
#ifndef NO_CVS_UP_IN_HERE
  }
#endif
    
	/*
vector<float> pt(3);
pt[0] = 75; pt[1] = 55; pt[2] = 51;
cout << "In DWI space: " << pt[0] << " " << pt[1] << " " << pt[2] << endl;
mAffineReg.ApplyXfm(pt, pt.begin());
mNonlinReg.ApplyXfm(pt, pt.begin());
cout << "In atlas space: " << pt[0] << " " << pt[1] << " " << pt[2] << endl;
exit(1);
	*/

  // Free atlas-space reference volume
  MRIfree(&atlasref);

  // Read segmentation map
  for (auto ifile = AsegList.begin(); ifile < AsegList.end(); ifile++) {
    MRI *aseg = 0;

    cout << "Loading segmentation map from " << *ifile << endl;
    aseg = MRIread((*ifile).c_str());
    if (!aseg) {
      cout << "ERROR: Could not read " << *ifile << endl;
      exit(1);
    }

    mAseg.push_back(aseg);
  }

  // Create output directory for current pathway for each time point
  SetOutputDir(OutDir);

  // Set atlas-derived information specific to current pathway
  SetPathway(InitFile, RoiFile1, RoiFile2,
             RoiMeshFile1, RoiMeshFile2, RoiRefFile1, RoiRefFile2,
             XyzPriorFile0, XyzPriorFile1,
             TangPriorFile, CurvPriorFile,
             NeighPriorFile, NeighIdFile,
             LocalPriorFile, LocalIdFile);

  // Set parameters for MCMC algorithm
  SetMcmcParameters(NumBurnIn, NumSample,
                    KeepSampleNth, UpdatePropNth, PropStdFile);
}

Coffin::~Coffin() {
  if (mMask != mDwi[0].GetMask())
    MRIfree(&mMask);

  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    idwi->FreeMask();

  for (vector<MRI *>::iterator iaseg = mAseg.begin(); iaseg < mAseg.end();
                                                      iaseg++)
    MRIfree(&(*iaseg));

  MRIfree(&mRoi1);
  MRIfree(&mRoi2);

  if (mXyzPrior0) {
    MRIfree(&mXyzPrior0);
    MRIfree(&mXyzPrior1);
  }
}

//
// Set output directory for each time point
//
void Coffin::SetOutputDir(const std::string OutDir) {
  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    idwi->SetOutputDir(OutDir);

  mOutDir = mDwi[0].GetOutputDir();
}

//
// Set atlas-derived information specific to a given pathway
//
void Coffin::SetPathway(const std::string InitFile,
                        const std::string RoiFile1, const std::string RoiFile2,
                        const std::string RoiMeshFile1, const std::string RoiMeshFile2,
                        const std::string RoiRefFile1, const std::string RoiRefFile2,
                        const std::string XyzPriorFile0, const std::string XyzPriorFile1,
                        const std::string TangPriorFile, const std::string CurvPriorFile,
                        const std::string NeighPriorFile, const std::string NeighIdFile,
                        const std::string LocalPriorFile, const std::string LocalIdFile) {
  int dirs[45] = { 0,  0,  0,
                   1,  0,  0,
                  -1,  0,  0,
                   0,  1,  0,
                   0, -1,  0,
                   0,  0,  1,
                   0,  0, -1,
                   1,  1,  1,
                  -1,  1,  1,
                   1, -1,  1,
                  -1, -1,  1,
                   1,  1, -1,
                  -1,  1, -1,
                   1, -1, -1,
                  -1, -1, -1 };
  ostringstream infostr;
  vector< vector<unsigned int> > segids, neighids, localids;
  vector< vector<float> > neighpr, localpr;

  // Save input info for logging
  infostr << "Initial control point file: " << InitFile << endl
          << "End ROI 1: " << RoiFile1 << endl;
  if (!RoiMeshFile1.empty())
    infostr << "End ROI 1 mesh: " << RoiMeshFile1 << endl
            << "End ROI 1 reference volume: " << RoiRefFile1 << endl;
  infostr << "End ROI 2: " << RoiFile2 << endl;
  if (!RoiMeshFile2.empty())
    infostr << "End ROI 2 mesh: " << RoiMeshFile2 << endl
            << "End ROI 2 reference volume: " << RoiRefFile2 << endl;
  if (!XyzPriorFile0.empty())
    infostr << "Spatial prior (off path): " << XyzPriorFile0 << endl
            << "Spatial prior (on path): " << XyzPriorFile1 << endl;
  if (!TangPriorFile.empty())
    infostr << "Tangent prior: " << TangPriorFile << endl;
  if (!CurvPriorFile.empty())
    infostr << "Curvature prior: " << CurvPriorFile << endl;
  if (!NeighPriorFile.empty())
    infostr << "Neighbor aseg prior: " << NeighPriorFile << endl
            << "Neighbor aseg label ID list: " << NeighIdFile << endl;
  if (!LocalPriorFile.empty())
    infostr << "Local aseg prior: " << LocalPriorFile << endl
            << "Local aseg label ID list: " << LocalIdFile << endl;
  mInfoPathway = infostr.str();

  // Read start ROI
  if (!RoiFile1.empty()) {
    if (mRoi1)
      MRIfree(&mRoi1);

    cout << "Loading end ROI from " << RoiFile1 << endl;
    mRoi1 = MRIread(RoiFile1.c_str());
    if (!mRoi1) {
      cout << "ERROR: Could not read " << RoiFile1 << endl;
      exit(1);
    }

    if (!RoiMeshFile1.empty()) {
      cout << "ERROR: .label ROIs not supported" << endl;
      exit(1);
    }
  }

  // Read end ROI
  if (!RoiFile2.empty()) {
    if (mRoi2) {
      MRIfree(&mRoi2);
    }

    cout << "Loading end ROI from " << RoiFile2 << endl;
    mRoi2 = MRIread(RoiFile2.c_str());
    if (!mRoi2) {
      cout << "ERROR: Could not read " << RoiFile2 << endl;
      exit(1);
    }

    if (!RoiMeshFile2.empty()) {
      cout << "ERROR: .label ROIs not supported" << endl;
      exit(1);
    }
  }

  // Size of atlas space
  mNxAtlas = mRoi1->width;
  mNyAtlas = mRoi1->height;
  mNzAtlas = mRoi1->depth;

  // Read control point initialization
  ReadControlPoints(InitFile);
  if (!mProposalStdInit.empty() &&
      (mProposalStdInit.size() != mControlPoints.size())) {
    mProposalStdInit.resize(mControlPoints.size(), 1.0);
    cout << "WARN: Resized proposal standard deviations to match new number "
         << "of control points" << endl;
  }

  // Read spatial path priors
  if ((!XyzPriorFile0.empty()) && (!XyzPriorFile1.empty())) {
    if (mXyzPrior0) {
      MRIfree(&mXyzPrior0);
    }

    cout << "Loading spatial path prior from " << XyzPriorFile0 << endl;
    mXyzPrior0 = MRIread(XyzPriorFile0.c_str());
    if (!mXyzPrior0) {
      cout << "ERROR: Could not read " << XyzPriorFile0 << endl;
      exit(1);
    }

    if (mXyzPrior1) {
      MRIfree(&mXyzPrior1);
    }

    cout << "Loading spatial path prior from " << XyzPriorFile1 << endl;
    mXyzPrior1 = MRIread(XyzPriorFile1.c_str());
    if (!mXyzPrior1) {
      cout << "ERROR: Could not read " << XyzPriorFile1 << endl;
      exit(1);
    }
  }
  else {
    mXyzPrior0 = 0;
    mXyzPrior1 = 0;
  }

  mNumArc = 0;

  // Read path tangent prior
  if (!TangPriorFile.empty()) {
    const int nbin = (int) ceil(2 / mTangentBinSize),
              nbin2 = nbin*nbin;
    string prline;
    ifstream prfile;

    mPriorTangent.clear();

    cout << "Loading tangent prior from " << TangPriorFile << endl;

    prfile.open(TangPriorFile, ios::in);
    if (!prfile) {
      cout << "ERROR: Could not open " << TangPriorFile << endl;
      exit(1);
    }

    while (getline(prfile, prline)) {
      float pr;
      vector<float> prior;
      istringstream prstr(prline);

      while (prstr >> pr)
        prior.push_back(pr);

      if (prior.size() != (unsigned int) nbin2) {
        cout << "ERROR: Line length is " << prior.size()
             << ", expected" << nbin2 << endl;
        exit(1);
      }

      mPriorTangent.push_back(prior);
    }

    prfile.close();

    mNumArc = (int) mPriorTangent.size();
  }

  // Read path curvature prior
  if (!CurvPriorFile.empty()) {
    string prline;
    ifstream prfile;

    mPriorCurvature.clear();

    cout << "Loading curvature prior from " << CurvPriorFile << endl;

    prfile.open(CurvPriorFile, ios::in);
    if (!prfile) {
      cout << "ERROR: Could not open " << CurvPriorFile << endl;
      exit(1);
    }

    while (getline(prfile, prline)) {
      float pr;
      vector<float> prior;
      istringstream prstr(prline);

      while (prstr >> pr)
        prior.push_back(pr);

      mPriorCurvature.push_back(prior);
    }

    prfile.close();

    if (mNumArc == 0)
      mNumArc = (int) mPriorCurvature.size();
    else if (mNumArc != (int) mPriorCurvature.size()) {
      cout << "ERROR: Mismatch between the numbers of arc segments in "
           << TangPriorFile
           << " (" << mPriorTangent.size() << "), "
           << CurvPriorFile
           << " (" << mPriorCurvature.size() << ")" << endl;
      exit(1);
    }
  }

  // Read neighbor aseg priors
  if ((!NeighPriorFile.empty()) && (!NeighIdFile.empty())) {
    mPriorNear.clear();
    mIdsNear.clear();
    mDirNear.clear();

    // Directions in which to look for neighboring labels
    if (mPriorSetNear == 6)
      mDirNear.insert(mDirNear.begin(), dirs+3, dirs+21);
    else if (mPriorSetNear == 14)
      mDirNear.insert(mDirNear.begin(), dirs+3, dirs+45);

    for (vector<int>::const_iterator idir = mDirNear.begin();
                                     idir < mDirNear.end(); idir += 3) {
      const int idx = idir[0], idy = idir[1], idz = idir[2];
      std::stringstream prname, idname;
      string prline, idline;
      ifstream prfile, idfile;

      prname << NeighPriorFile << '_'
	     << idx << '_'
	     << idy << '_'
	     << idz << ".txt";
      idname << NeighIdFile << '_'
	     << idx << '_'
	     << idy << '_'
	     << idz << ".txt";

      cout << "Loading nearest neighbor prior from " << prname.str()
           << " with label IDs from " << idname.str() << endl;

      prfile.open(prname.str(), ios::in);
      if (!prfile) {
        cout << "ERROR: Could not open " << prname.str() << endl;
        exit(1);
      }

      idfile.open(idname.str(), ios::in);
      if (!idfile) {
        cout << "ERROR: Could not open " << idname.str() << endl;
        exit(1);
      }

      while (getline(prfile, prline) && getline(idfile, idline)) {
        unsigned int id;
        float pr;
        vector<unsigned int> idlist;
        vector<float> prior;
        istringstream prstr(prline), idstr(idline);

        while (prstr >> pr)
          prior.push_back(pr);

        while (idstr >> id)
          idlist.push_back(id);

        if (prior.size() != idlist.size() + 1) {
          cout << "ERROR: Line length mismatch between "
               << prname.str() << " (" << prline << ") and "
               << idname.str() << " (" << idline << ")" << endl;
          exit(1);
        } 

        mPriorNear.push_back(prior);
        mIdsNear.push_back(idlist);
      }

      prfile.close();
      idfile.close();
    }

    if (mNumArc == 0)
      mNumArc = (int) (mPriorNear.size() / mPriorSetNear);
    else if (mNumArc != (int) (mPriorNear.size() / mPriorSetNear)) {
      cout << "ERROR: Mismatch between the numbers of arc segments in ";
      if (!TangPriorFile.empty())
        cout << TangPriorFile
             << " (" << mPriorTangent.size() << "), ";
      if (!CurvPriorFile.empty())
        cout << CurvPriorFile
             << " (" << mPriorCurvature.size() << "), ";
      cout << NeighPriorFile
           << " (" << mPriorNear.size()  / mPriorSetNear  << ")" << endl;
    }
  }

  // Read local aseg priors
  if ((!LocalPriorFile.empty()) && (!LocalIdFile.empty())) {
    mPriorLocal.clear();
    mIdsLocal.clear();
    mDirLocal.clear();

    // Directions in which to look for neighboring labels
    if (mPriorSetLocal == 1)
      mDirLocal.insert(mDirLocal.begin(), dirs, dirs+3);
    else if (mPriorSetLocal == 7)
      mDirLocal.insert(mDirLocal.begin(), dirs, dirs+21);
    else if (mPriorSetLocal == 15)
      mDirLocal.insert(mDirLocal.begin(), dirs, dirs+45);

    for (vector<int>::const_iterator idir = mDirLocal.begin();
                                     idir < mDirLocal.end(); idir += 3) {
      const int idx = idir[0], idy = idir[1], idz = idir[2];
      std::stringstream prname, idname;
      string prline, idline;
      ifstream prfile, idfile;

      prname << LocalPriorFile << "_" << idx << "_" << idy << '_' << idz << ".txt";
      idname << LocalIdFile << '_' << idx << "_" << idy << '_' << idz << ".txt";

      cout << "Loading local prior from " << prname.str()
           << " with label IDs from " << idname.str() << endl;

      prfile.open(prname.str(), ios::in);
      if (!prfile) {
        cout << "ERROR: Could not open " << prname.str() << endl;
        exit(1);
      }

      idfile.open(idname.str(), ios::in);
      if (!idfile) {
        cout << "ERROR: Could not open " << idname.str() << endl;
        exit(1);
      }

      while (getline(prfile, prline) && getline(idfile, idline)) {
        unsigned int id;
        float pr;
        vector<unsigned int> idlist;
        vector<float> prior;
        istringstream prstr(prline), idstr(idline);

        while (prstr >> pr) {
          prior.push_back(pr);
	}

        while (idstr >> id) {
          idlist.push_back(id);
	}

        if (prior.size() != idlist.size() + 1) {
          cout << "ERROR: Line length mismatch between "
               << prname.str() << " (" << prline << ") and "
               << idname.str() << " (" << idline << ")" << endl;
          exit(1);
        } 

        mPriorLocal.push_back(prior);
        mIdsLocal.push_back(idlist);
      }

      prfile.close();
      idfile.close();
    }

    if (mNumArc == 0) {
      mNumArc = (int) (mPriorLocal.size() / mPriorSetLocal);
    } else if (mNumArc != (int) (mPriorLocal.size() / mPriorSetLocal)) {
      cout << "ERROR: Mismatch between the numbers of arc segments in ";
      if (!TangPriorFile.empty())
        cout << TangPriorFile
             << " (" << mPriorTangent.size() << "), ";
      if (!CurvPriorFile.empty())
        cout << CurvPriorFile
             << " (" << mPriorCurvature.size() << "), ";
      if (!NeighPriorFile.empty())
        cout << NeighPriorFile
             << " (" << mPriorNear.size() / mPriorSetNear << "), ";
      cout << LocalPriorFile
           << " (" << mPriorLocal.size() / mPriorSetLocal << ")" << endl;
      exit(1);
    }
  }
}

//
// Set MCMC parameters
//
void Coffin::SetMcmcParameters(const int NumBurnIn, const int NumSample,
                               const int KeepSampleNth, const int UpdatePropNth,
                               const std::string PropStdFile) {
  ostringstream infostr;

  // Save input info for logging
  infostr << "Number of burn-in samples: " << NumBurnIn << endl
          << "Number of post-burn-in samples: " << NumSample << endl
          << "Keep every: " << KeepSampleNth << "-th sample" << endl
          << "Update proposal every: " << UpdatePropNth << "-th sample" << endl;
  if (!PropStdFile.empty()) {
    infostr << "Initial proposal SD file: " << PropStdFile << endl;
  }
  mInfoMcmc = infostr.str();

  // Set sampling parameters
  mNumBurnIn = NumBurnIn;
  mNumSample = NumSample;
  mKeepSampleNth = KeepSampleNth;
  mUpdatePropNth = UpdatePropNth;
 
  // Read proposal SD initialization
  ReadProposalStds(PropStdFile);
  if (mProposalStdInit.size() != mControlPoints.size()) {
    cout << "ERROR: Dimensions of " << PropStdFile << " must be "
         << mControlPoints.size()/3 << " x 3 to match control points" << endl;
    exit(1);
  }
}

//
// Read initial control points
//
void Coffin::ReadControlPoints(const std::string ControlPointFile) {
  float coord;
  ifstream infile(ControlPointFile, ios::in);

  if (!infile) {
    cout << "ERROR: Could not open " << ControlPointFile << endl;
    exit(1);
  }

  cout << "Loading control points from " << ControlPointFile << endl;
  mControlPoints.clear();
  while (infile >> coord)
    mControlPoints.push_back((int) round(coord));

  if (mControlPoints.size() % 3 != 0) {
    cout << "ERROR: File " << ControlPointFile
         << " must contain triplets of coordinates" << endl;
    exit(1);
  }

  mNumControl = mControlPoints.size()/3;

  // Make sure that initial control points are in mask
  for (vector<int>::iterator icpt = mControlPoints.begin();
                             icpt < mControlPoints.end(); icpt += 3) {
    if (icpt[0] < 0) {
      cout << "WARN: Initial control point "
           << icpt[0] << " " << icpt[1] << " " << icpt[2]
           << " is not in DWI volume - is DWI cropped?" << endl;

      icpt[0] = 0;

      cout << "WARN: Replacing with closest point in volume ("
           << icpt[0] << " " << icpt[1] << " " << icpt[2] << ")" << endl;
    }

    if (icpt[0] >= mMask->width) {
      cout << "WARN: Initial control point "
           << icpt[0] << " " << icpt[1] << " " << icpt[2]
           << " is not in DWI volume - is DWI cropped?" << endl;

      icpt[0] = mMask->width - 1;

      cout << "WARN: Replacing with closest point in volume ("
           << icpt[0] << " " << icpt[1] << " " << icpt[2] << ")" << endl;
    }

    if (icpt[1] < 0) {
      cout << "WARN: Initial control point "
           << icpt[0] << " " << icpt[1] << " " << icpt[2]
           << " is not in DWI volume - is DWI cropped?" << endl;

      icpt[1] = 0;

      cout << "WARN: Replacing with closest point in volume ("
           << icpt[0] << " " << icpt[1] << " " << icpt[2] << ")" << endl;
    }

    if (icpt[1] >= mMask->height) {
      cout << "WARN: Initial control point "
           << icpt[0] << " " << icpt[1] << " " << icpt[2]
           << " is not in DWI volume - is DWI cropped?" << endl;

      icpt[1] = mMask->height - 1;

      cout << "WARN: Replacing with closest point in volume ("
           << icpt[0] << " " << icpt[1] << " " << icpt[2] << ")" << endl;
    }

    if (icpt[2] < 0) {
      cout << "WARN: Initial control point "
           << icpt[0] << " " << icpt[1] << " " << icpt[2]
           << " is not in DWI volume - is DWI cropped?" << endl;

      icpt[2] = 0;

      cout << "WARN: Replacing with closest point in volume ("
           << icpt[0] << " " << icpt[1] << " " << icpt[2] << ")" << endl;
    }

    if (icpt[2] >= mMask->depth) {
      cout << "WARN: Initial control point "
           << icpt[0] << " " << icpt[1] << " " << icpt[2]
           << " is not in DWI volume - is DWI cropped?" << endl;

      icpt[2] = mMask->depth - 1;

      cout << "WARN: Replacing with closest point in volume ("
           << icpt[0] << " " << icpt[1] << " " << icpt[2] << ")" << endl;
    }

    if (! MRIgetVoxVal(mMask, icpt[0], icpt[1], icpt[2], 0)) {
      int dmin = 1000000, ixmin = 0, iymin = 0, izmin = 0;

      cout << "WARN: Initial control point "
           << icpt[0] << " " << icpt[1] << " " << icpt[2]
           << " is not in mask" << endl
           << "WARN: Replacing with closest point in mask (";

      for (int iz = 0; iz < mNz; iz++)
        for (int iy = 0; iy < mNy; iy++)
          for (int ix = 0; ix < mNx; ix++)
            if (MRIgetVoxVal(mMask, ix, iy, iz, 0)) {
              const int dx = icpt[0] - ix,
                        dy = icpt[1] - iy,
                        dz = icpt[2] - iz,
                        dist = dx*dx + dy*dy + dz*dz;

              if (dist < dmin) {
                ixmin = ix;
                iymin = iy;
                izmin = iz;
                dmin = dist;
              }
            }

      icpt[0] = ixmin;
      icpt[1] = iymin;
      icpt[2] = izmin;

      cout << icpt[0] << " " << icpt[1] << " " << icpt[2] << ")" << endl;
    }
  }

  // Make sure that initial start point is in start ROI
  if (mRoi1) {
    vector<int>::iterator icpt = mControlPoints.begin();

    if (!IsInRoi(icpt, mRoi1)) {
      int dmin = 1000000, ixmin = 0, iymin = 0, izmin = 0;
      vector<int> newpoint(3);

      cout << "WARN: Initial start point "
           << icpt[0] << " " << icpt[1] << " " << icpt[2]
           << " is not in start ROI" << endl
           << "WARN: Replacing with closest point in start ROI (";

      for (int iz = 0; iz < mNz; iz++)
        for (int iy = 0; iy < mNy; iy++)
          for (int ix = 0; ix < mNx; ix++) {
            newpoint[0] = ix;
            newpoint[1] = iy;
            newpoint[2] = iz;

            if (IsInRoi(newpoint.begin(), mRoi1)) {
              const int dx = icpt[0] - ix,
                        dy = icpt[1] - iy,
                        dz = icpt[2] - iz,
                        dist = dx*dx + dy*dy + dz*dz;

              if (dist < dmin) {
                ixmin = ix;
                iymin = iy;
                izmin = iz;
                dmin = dist;
              }
            }
          }

      icpt[0] = ixmin;
      icpt[1] = iymin;
      icpt[2] = izmin;

      cout << icpt[0] << " " << icpt[1] << " " << icpt[2] << ")" << endl;
    }
  }

  // Make sure that initial end point is in end ROI
  if (mRoi2) {
    vector<int>::iterator icpt = mControlPoints.end() - 3;

    if (!IsInRoi(icpt, mRoi2)) {
      int dmin = 1000000, ixmin = 0, iymin = 0, izmin = 0;
      vector<int> newpoint(3);

      cout << "WARN: Initial end point "
           << icpt[0] << " " << icpt[1] << " " << icpt[2]
           << " is not in end ROI" << endl
           << "WARN: Replacing with closest point in end ROI (";

      for (int iz = 0; iz < mNz; iz++)
        for (int iy = 0; iy < mNy; iy++)
          for (int ix = 0; ix < mNx; ix++) {
            newpoint[0] = ix;
            newpoint[1] = iy;
            newpoint[2] = iz;

            if (IsInRoi(newpoint.begin(), mRoi2)) {
              const int dx = icpt[0] - ix,
                        dy = icpt[1] - iy,
                        dz = icpt[2] - iz,
                        dist = dx*dx + dy*dy + dz*dz;

              if (dist < dmin) {
                ixmin = ix;
                iymin = iy;
                izmin = iz;
                dmin = dist;
              }
            }
          }

      icpt[0] = ixmin;
      icpt[1] = iymin;
      icpt[2] = izmin;

      cout << icpt[0] << " " << icpt[1] << " " << icpt[2] << ")" << endl;
    }
  }
}

//
// Read initial proposal standard deviations for control point perturbations
//
void Coffin::ReadProposalStds(const std::string PropStdFile) {
  mProposalStdInit.clear();

  if (!PropStdFile.empty()) {
    float val;
    ifstream infile(PropStdFile, ios::in);

    if (!infile) {
      cout << "ERROR: Could not open " << PropStdFile << endl;
      exit(1);
    }

    cout << "Loading initial proposal SD's from " << PropStdFile << endl;
    while (infile >> val)
      mProposalStdInit.push_back(val);
  }
  else {	// Default value (one DWI voxel, scaled to base voxel units)
    vector<float>::iterator istd;

    mProposalStdInit.resize(mControlPoints.size());
    istd = mProposalStdInit.begin();

    copy(mResolution.begin(), mResolution.end(), istd);			//x5.0
    for (istd += 3; istd < mProposalStdInit.end() - 3; istd +=3) {
      copy(mResolution.begin(), mResolution.end(), istd);		//x1.0
    }
    copy(mResolution.begin(), mResolution.end(), istd);			//x5.0
  }
}

//
// Run MCMC (full spline updates)
//
bool Coffin::RunMcmcFull() {
  int iprop, ikeep;
  std::string fname;
  string cmdline;

  // Open log file in first time point's output directory
  fname = mOutDir + "/log.txt";
  mLog.open(fname, ios::out | ios::app);
  if (!mLog) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  // Write input parameters to log file
  mLog << mInfoGeneral << mInfoPathway << mInfoMcmc;

  cout << "Initializing MCMC" << endl;
  mLog << "Initializing MCMC" << endl;
  if (! InitializeMcmc()) {
    mLog.flush();
    mLog.close();
    return false;
  }

  if (mDebug) {
    fname = mOutDir + "/Finit.nii.gz";
    mSpline.WriteVolume(fname.c_str(), true);
  }

  cout << "Running MCMC burn-in jumps" << endl;
  mLog << "Running MCMC burn-in jumps" << endl;
  iprop = 1;
  for (int ijump = mNumBurnIn; ijump > 0; ijump--) {
    if (JumpMcmcFull() || mAcceptF || mAcceptTheta) {	// Accept new path
      UpdatePath();
      UpdateAcceptanceRateFull();

      if (mDebug) {
	std::stringstream tmp;
	tmp << mOutDir << '/'
	    << "Faccept_b"
	    << std::setw(5) << std::setfill('0') << mNumBurnIn-ijump+1
	    << ".nii.gz";
	fname = tmp.str();
        mSpline.WriteVolume(fname.c_str(), true);
      }
    }
    else {						// Reject new path
      UpdateRejectionRateFull();

      if (mDebug) {
	std::stringstream tmp;
	tmp << mOutDir << '/'
	    << "Freject_b"
	    << std::setw(5) << std::setfill('0') << mNumBurnIn-ijump+1
	    << ".nii.gz";
	fname = tmp.str();
        mSpline.WriteVolume(fname.c_str(), true);
      }
    }

    if (iprop == mUpdatePropNth) {
      UpdateProposalStd();
      iprop = 1;
    }
    else
      iprop++;
  }

  mPosteriorOnPathMap = numeric_limits<double>::max();

  cout << "Running MCMC main jumps" << endl;
  mLog << "Running MCMC main jumps" << endl;
  iprop = 1;
  ikeep = 1;
  for (int ijump = mNumSample; ijump > 0; ijump--) {
    if (JumpMcmcFull() || mAcceptF || mAcceptTheta) {	// Accept new path
      SavePathPosterior(true);
      UpdatePath();
      UpdateAcceptanceRateFull();

      if (mDebug) {
	std::stringstream tmp;
	tmp << mOutDir << '/'
	    << "Faccept_"
	    << std::setw(5) << std::setfill('0') << mNumSample-ijump+1
	    << ".nii.gz";
	fname = tmp.str();
        mSpline.WriteVolume(fname.c_str(), true);
      }
    }
    else {						// Reject new path
      SavePathPosterior(false);
      UpdateRejectionRateFull();

      if (mDebug) {
	std::stringstream tmp;
	tmp << mOutDir << '/'
	    << "Freject_"
	    << std::setw(5) << std::setfill('0') << mNumSample-ijump+1
	    << ".nii.gz";
	fname = tmp.str();
        mSpline.WriteVolume(fname.c_str(), true);
      }
    }

    if (iprop == mUpdatePropNth) {
      UpdateProposalStd();
      iprop = 1;
    }
    else
      iprop++;

    if (ikeep == mKeepSampleNth) {
      SavePath();
      ikeep = 1;
    }
    else
      ikeep++;
  }

  // Close log file and copy it to other time points's output directories
  mLog.flush();
  mLog.close();

  for (vector<Aeon>::const_iterator idwi = mDwi.begin() + 1; idwi < mDwi.end();
                                                             idwi++) {
    cmdline = "cp -f " + mDwi[0].GetOutputDir() + "/log.txt " +
              idwi->GetOutputDir();

    if (system(cmdline.c_str()) != 0) {
      cout << "ERROR: Could not save log file in " << idwi->GetOutputDir()
           << endl;
      exit(1);
    }
  }

  return true;
}

//
// Run MCMC (single control point updates)
//
bool Coffin::RunMcmcSingle() {
  int iprop, ikeep;
  std::string fname;
  string cmdline;
  vector<int> cptorder(mNumControl);
  vector<int>::const_iterator icpt;

  // Open log file in first time point's output directory
  fname = mOutDir + "/log.txt";
  mLog.open(fname, ios::out | ios::app);
  if (!mLog) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  // Write input parameters to log file
  mLog << mInfoGeneral << mInfoPathway << mInfoMcmc;

  cout << "Initializing MCMC" << endl;
  mLog << "Initializing MCMC" << endl;
  if (! InitializeMcmc()) {
    mLog.flush();
    mLog.close();
    return false;
  }

  if (mDebug) {
    fname = mOutDir + "/Finit.nii.gz";
    mSpline.WriteVolume(fname.c_str(), true);
  }

  cout << "Running MCMC burn-in jumps" << endl;
  mLog << "Running MCMC burn-in jumps" << endl;
  iprop = 1;
  for (int ijump = mNumBurnIn; ijump > 0; ijump--) {
    // Perturb control points in random order
    for (int k = 0; k < mNumControl; k++) {
      cptorder[k] = k;
    }
    random_shuffle(cptorder.begin(), cptorder.end());

    fill(mRejectControl.begin(), mRejectControl.end(), false);

    for (icpt = cptorder.begin(); icpt != cptorder.end(); icpt++) {
      mRejectSpline = false;
      mRejectF = false;
      mAcceptF = false;
      mRejectTheta = false;
      mAcceptTheta = false;
      mRejectPosterior = false;

      if (JumpMcmcSingle(*icpt) || mAcceptF || mAcceptTheta) {	// Accept point
        UpdatePath();

        if (mDebug) {
	  std::stringstream tmp;
	  tmp << mOutDir << '/'
	      << "Faccept_b"
	      << std::setw(5) << std::setfill('0') << mNumBurnIn-ijump+1
	      << '_'
	      << *icpt
	      << ".nii.gz";
	  fname = tmp.str();
          mSpline.WriteVolume(fname.c_str(), true);
        }
      }
      else {							// Reject point
        mRejectControl[*icpt] = true;

        if (mDebug) {
	  std::stringstream tmp;
	  tmp << mOutDir << '/'
	      << "Freject_b"
	      << std::setw(5) << std::setfill('0') << mNumBurnIn-ijump+1
	      << '_'
	      << *icpt
	      << ".nii.gz";
	  fname = tmp.str();
          mSpline.WriteVolume(fname.c_str(), true);
        }
      }
    }

    UpdateAcceptRejectRateSingle();

    if (iprop == mUpdatePropNth) {
      UpdateProposalStd();
      iprop = 1;
    }
    else
      iprop++;
  }

  mPosteriorOnPathMap = numeric_limits<double>::max();

  cout << "Running MCMC main jumps" << endl;
  mLog << "Running MCMC main jumps" << endl;
  iprop = 1;
  ikeep = 1;
  for (int ijump = mNumSample; ijump > 0; ijump--) {
    // Perturb control points in random order
    for (int k = 0; k < mNumControl; k++)
      cptorder[k] = k;
    random_shuffle(cptorder.begin(), cptorder.end());

    fill(mRejectControl.begin(), mRejectControl.end(), false);

    for (icpt = cptorder.begin(); icpt != cptorder.end(); icpt++) {
      mRejectSpline = false;
      mRejectF = false;
      mAcceptF = false;
      mRejectTheta = false;
      mAcceptTheta = false;
      mRejectPosterior = false;

      if (JumpMcmcSingle(*icpt) || mAcceptF || mAcceptTheta) {	// Accept point
        SavePathPosterior(true);
        UpdatePath();

        if (mDebug) {
	  std::stringstream tmp;
	  tmp << mOutDir << '/'
	      << "Faccept_"
	      << std::setw(5) << std::setfill('0') << mNumSample-ijump+1
	      << '_'
	      << *icpt;
	  fname = tmp.str();
          mSpline.WriteVolume(fname.c_str(), true);
        }
      }
      else {							// Reject point
        SavePathPosterior(false);
        mRejectControl[*icpt] = true;

        if (mDebug) {
	  std::stringstream tmp;
	  tmp << mOutDir << '/'
	      << "Freject_"
	      << std::setw(5) << std::setfill('0') << mNumSample-ijump+1
	      << '_'
	      << *icpt;
	  fname = tmp.str();
          mSpline.WriteVolume(fname.c_str(), true);
        }
      }
    }

    UpdateAcceptRejectRateSingle();

    if (iprop == mUpdatePropNth) {
      UpdateProposalStd();
      iprop = 1;
    }
    else
      iprop++;

    if (ikeep == mKeepSampleNth) {
      SavePath();
      ikeep = 1;
    }
    else
      ikeep++;
  }

  // Close log file and copy it to other time points's output directories
  mLog.flush();
  mLog.close();

  for (vector<Aeon>::const_iterator idwi = mDwi.begin() + 1; idwi < mDwi.end();
                                                             idwi++) {
    cmdline = "cp -f " + mDwi[0].GetOutputDir() + "/log.txt " +
              idwi->GetOutputDir();

    if (system(cmdline.c_str()) != 0) {
      cout << "ERROR: Could not save log file in " << idwi->GetOutputDir()
           << endl;
      exit(1);
    }
  }

  return true;
}

//
// Initialize path and MCMC proposals
//
bool Coffin::InitializeMcmc() {
  bool success = true, doinit = true, firstinit = true;
  int failseg = -1;
  vector<int> atlaspoints;
  vector<int>::iterator iptatlas;

  // Initialize control point proposal distribution
  mProposalStd.resize(mProposalStdInit.size());
  copy(mProposalStdInit.begin(), mProposalStdInit.end(), mProposalStd.begin());

  if (mDebug) {
    mLog << "Proposal STDs: ";
    for (vector<float>::const_iterator pstd = mProposalStd.begin();
                                       pstd < mProposalStd.end(); pstd++)
      mLog << *pstd << " ";
    mLog << endl;
  }

  // Initialize jump acceptance statistics
  mRejectControl.resize(mNumControl);
  fill(mRejectControl.begin(), mRejectControl.end(), false);

  mAcceptCount.resize(mNumControl);
  fill(mAcceptCount.begin(), mAcceptCount.end(), 0);
  mRejectCount.resize(mNumControl);
  fill(mRejectCount.begin(), mRejectCount.end(), 0);

  mAcceptSpan.resize(mControlPoints.size());
  fill(mAcceptSpan.begin(), mAcceptSpan.end(), 0.0);
  mRejectSpan.resize(mControlPoints.size());
  fill(mRejectSpan.begin(), mRejectSpan.end(), 0.0);

  mControlPointsNew.resize(mControlPoints.size());
  mControlPointJumps.resize(mControlPoints.size());

  // Interpolate spline from initial control points
  mSpline.SetControlPoints(mControlPoints);

  if (!mSpline.InterpolateSpline()) {	// Initial path goes off brain mask
    failseg = FindErrorSegment();

    if (!InitializeFixOffMask(failseg)) {
      cout << "ERROR: Path from initial control points is not entirely"
           << " in mask or is irregular" << endl
           << "ERROR: Initialization failed" << endl;
      return false;
    }
  }
  else {
    copy(mControlPoints.begin(), mControlPoints.end(),
         mControlPointsNew.begin());

    // Clear path-related variables for all time points
    for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
      idwi->ClearPath();

    // Propagate initial path to all time points
    for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end();
                                                     idwi++) {
      success = idwi->MapPathFromBase(mSpline);
      if (!success) {
        failseg = idwi->FindErrorSegment(mSpline);
        break;
      }
    }

    if (!success) {			// Initial path goes off brain mask
      if (!InitializeFixOffMask(failseg)) {
        cout << "ERROR: Path from initial control points is not entirely"
             << " in mask or is irregular" << endl
             << "ERROR: Initialization failed" << endl;
        return false;
      }
    }
  }

  while (doinit) {
    // Get initial path points
    mPathPointsNew.clear();
    mPathPointsNew.insert(mPathPointsNew.begin(), mSpline.GetAllPointsBegin(),
                                                  mSpline.GetAllPointsEnd());

    // Compute data-fit terms on initial path for all time points
    mDataPosteriorOnPathNew = 0;
    mDataPosteriorOffPathNew = 0;

    for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end();
                                                     idwi++) {
      success = idwi->ComputePathDataFit();
      if (!success) {
        failseg = idwi->FindErrorSegment(mSpline);
        break;
      }

      mDataPosteriorOnPathNew  += idwi->GetPosteriorOnPathNew();
      mDataPosteriorOffPathNew += idwi->GetPosteriorOffPathNew();
    }

    // Map initial path from diffusion/base space to atlas space
    atlaspoints.resize(mPathPointsNew.size());
    iptatlas = atlaspoints.begin();

    for (vector<int>::const_iterator ipt = mPathPointsNew.begin();
                                     ipt < mPathPointsNew.end(); ipt += 3) {
      vector< vector<int> >::iterator icoord =
        mAtlasCoords.begin() + (ipt[0] + ipt[1]*mNx + ipt[2]*mNxy);

      if (icoord->empty()) {
        // Map point coordinates from diffusion/base space to atlas space
        MapPointToAtlas(iptatlas, ipt);

        // Save transformed coordinates for future use
        icoord->insert(icoord->end(), iptatlas, iptatlas+3);
      }
      else
        // Retrieve previously computed coordinates in atlas space
        copy(icoord->begin(), icoord->end(), iptatlas);

      iptatlas += 3;
    }

    // Compute atlas-derived prior terms on initial path
    mXyzPriorOnPathNew = ComputeXyzPriorOnPath(atlaspoints);

    mAnatomicalPriorNew = ComputeAnatomicalPrior(atlaspoints);

    mShapePriorNew = ComputeShapePrior(atlaspoints);

    // Compute posterior on initial path
    mPosteriorOnPathNew  = mDataPosteriorOnPathNew + mXyzPriorOnPathNew
                         + mAnatomicalPriorNew + mShapePriorNew;
    mPosteriorOffPathNew = mDataPosteriorOffPathNew + mXyzPriorOffPathNew;

    // Save initial path
    UpdatePath();

    if (!success && firstinit) {	// Initial path goes off white matter
      doinit = !InitializeFixOffWhite(failseg);
      firstinit = false;
    }
    else
      doinit = false;
  }

  return true;
}

//
// Try to fix the initial path if it is not entirely in the brain mask
//
bool Coffin::InitializeFixOffMask(int FailSegment) {
  bool success = false;
  int failseg = FailSegment,
      perturbfirst = failseg,
      perturblast  = (failseg == mNumControl-1) ? failseg : failseg+1;
  const float maxdist2 = 4 * ( mResolution[0]*mResolution[0] +
                               mResolution[1]*mResolution[1] +
                               mResolution[2]*mResolution[2] ) / 3.0;
  vector<int> controlorig(mControlPoints);

  // Perturb control points until I get a valid initial path
  for (unsigned int itry = 0; itry < mMaxTryMask; itry++) {
    cout << "INFO: Path from initial control points is not entirely in mask"
         << " or is irregular" << endl
         << "INFO: in segment " << failseg << endl
         << "INFO: Attempting to perturb control points" << endl;

    fill(mRejectControl.begin(), mRejectControl.end(), false);

    for (int icpt = perturbfirst; icpt <= perturblast; icpt++) {
      mRejectSpline = false;

      success = ProposePathSingle(icpt);

      if (success) {
        vector<int>::const_iterator inew = mControlPointsNew.begin() + 3*icpt,
                                    iold = controlorig.begin() + 3*icpt;
        const int dx = inew[0] - iold[0],
                  dy = inew[1] - iold[1],
                  dz = inew[2] - iold[2];

        if (dx*dx + dy*dy + dz*dz > maxdist2)	// Too far from original points
          success = false;
        else {
          copy(mControlPointsNew.begin(), mControlPointsNew.end(),
               mControlPoints.begin());
          break;
        }
      }
      else if (mRejectSpline) {		// Error may now be in another segment
        int failsegnew = FindErrorSegment();

        if (failsegnew > -1) {
          failseg = failsegnew;
          perturbfirst = failseg;
          perturblast  = (failseg == mNumControl-1) ? failseg : failseg+1;
        }
        else
          for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end();
                                                           idwi++) {
            failsegnew = idwi->FindErrorSegment(mSpline);

            if (failsegnew > -1) {
              failseg = failsegnew;
              perturbfirst = failseg;
              perturblast  = (failseg == mNumControl-1) ? failseg : failseg+1;
              break;
            }
          }

        if (failsegnew > -1)
          break;
      }
    }

    if (success) {
      cout << "INFO: Success after " << itry+1 << " perturbation(s)" << endl;
      break;
    }
  }

  return success;
}

//
// Try to fix the initial path if it is not entirely in the white matter
//
bool Coffin::InitializeFixOffWhite(int FailSegment) {
  bool success = true, improved = true;
  int failseg = FailSegment,
      perturbfirst = failseg,
      perturblast  = (failseg == mNumControl-1) ? failseg : failseg+1;
  std::string fname;
  vector<int> controlorig(mControlPoints);

  // Set proposal standard deviations to a conservative value for this
  if (mMaxTryWhite > 0) {
    for (vector<float>::iterator istd = mProposalStd.begin();
	 istd < mProposalStd.end(); istd += 3) {
      for (int k = 0; k < 3; k++) {
        if (istd[k] > mResolution[k]) {
          istd[k] = mResolution[k];
	}
      }
    }
  }

  // Perturb control points to find a valid initial path
  for (unsigned int itry = 0; itry < mMaxTryWhite; itry++) {
    cout << "INFO: Path from initial control points is not entirely in "
         << "white matter" << endl
         << "INFO: in segment " << failseg << endl
         << "INFO: Attempting to perturb control points" << endl;

    fill(mRejectControl.begin(), mRejectControl.end(), false);

    for (int icpt = perturbfirst; icpt <= perturblast; icpt++) {
      mRejectSpline = false;
      mRejectF = false;
      mAcceptF = false;
      mRejectTheta = false;
      mAcceptTheta = false;
      mRejectPosterior = false;

      // If the new path is on the white matter, accept it
      success = (JumpMcmcSingle(icpt) || mAcceptF || mRejectTheta); 

      // If the new path is still off the white matter 
      // but less so than the current one,
      // decide to accept it based on the path prior only
      if (mRejectF || mAcceptTheta) {
        unsigned int noffnew = 0, noff = 0,
                     ntotnew = 0, ntot = 0;

        for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end();
                                                         idwi++) {
          noffnew += idwi->GetNumFZerosNew();
          noff    += idwi->GetNumFZeros();
          ntotnew += idwi->GetPathLengthNew();
          ntot    += idwi->GetPathLength();
        }

        improved = ( (noffnew / (float) ntotnew) < (noff / (float) ntot)
                     && AcceptPath(true) );
      }
      else {
        improved = success;
      }

      if (success || improved) {				// Accept point
        UpdatePath();

        if (mDebug) {
	  std::stringstream tmp;
	  tmp << mOutDir << '/'
	      << "Faccept_f"
	      << std::setw(5) << std::setfill('0') << itry+1
	      << '_'
	      << icpt
	      << ".nii.gz";
	  fname = tmp.str();
          mSpline.WriteVolume(fname.c_str(), true);
        }

        if (success)
          break;
        else
          for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end();
                                                           idwi++) {
            const int failsegnew = idwi->FindErrorSegment(mSpline);

            if (failsegnew > -1 && failseg != failsegnew) {
              failseg = failsegnew;
              perturbfirst = failseg;
              perturblast  = (failseg == mNumControl-1) ? failseg : failseg+1;
              break;
            }
          }
      }
      else {							// Reject point
        mRejectControl[icpt] = true;

        if (mDebug) {
	  std::stringstream tmp;
	  tmp << mOutDir << '/'
	      << "Freject_f"
	      << std::setw(5) << std::setfill('0') << itry+1
	      << '_'
	      << icpt
	      << ".nii.gz";
	  fname = tmp.str();
          mSpline.WriteVolume(fname.c_str(), true);
        }
      }
    }

    if (success) {
      cout << "INFO: Success after " << itry+1 << " perturbation(s)" << endl;
      break;
    }
  }

  if (!success) {
    cout << "INFO: Reverting to original control points" << endl;
    copy(controlorig.begin(), controlorig.end(), mControlPointsNew.begin());
    mSpline.SetControlPoints(mControlPointsNew);
    mSpline.InterpolateSpline();
  }

  // Restore proposal standard deviations
  if (mMaxTryWhite > 0)
    copy(mProposalStdInit.begin(), mProposalStdInit.end(),
         mProposalStd.begin());

  return success;
}

//
// Find spline segment corresponding to the point along the path
// where an error occured
//
int Coffin::FindErrorSegment() {
  const int ierr = (mSpline.GetAllPointsEnd() -
		    mSpline.GetAllPointsBegin())/3 - 1;
  int iseg;
  
  if (ierr < 0) {
    iseg = 0;
  } else {
    iseg = (int) mSpline.PointToSegment(ierr);
  }
  
  if (iseg == mNumControl-1) {
    // Spline interpolation was completed w/o error
    iseg = -1;
  }

  return iseg;
}

//
// Perform one MCMC jump (full spline update)
//
bool Coffin::JumpMcmcFull() {
  fill(mRejectControl.begin(), mRejectControl.end(), false);
  mRejectSpline = false;
  mRejectF = false;
  mAcceptF = false;
  mRejectTheta = false;
  mAcceptTheta = false;
  mRejectPosterior = false;

  if (!ProposePathFull())
    return false;

  ProposeDiffusionParameters();

  if (!AcceptPath())
    return false;

  return true;
}

//
// Perform one MCMC jump (single control point update)
//
bool Coffin::JumpMcmcSingle(int ControlIndex) {
  if (!ProposePathSingle(ControlIndex))
    return false;

  ProposeDiffusionParameters();

  if (!AcceptPath())
    return false;

  return true;
}

//
// Propose path by perturbing all control points
//
bool Coffin::ProposePathFull() {
  vector<bool>::iterator isrej;
  vector<int>::const_iterator coord = mControlPoints.begin();
  vector<int>::const_iterator cpoint = mControlPointsNew.begin();
  vector<int>::iterator newcoord = mControlPointsNew.begin();
  vector<float>::const_iterator pstd = mProposalStd.begin();
  vector<float>::iterator jump = mControlPointJumps.begin();

  // Perturb current control points
  if (mDebug)
    mLog << "Jump norms: ";
  for (isrej = mRejectControl.begin(); isrej < mRejectControl.end(); isrej++) {
    double norm = 0;

    for (int ii = 0; ii < 3; ii++) {
      *jump = round((*pstd) * PDFgaussian());
      *newcoord = *coord + (int) *jump;

      *jump *= *jump;
      norm += *jump;

      jump++;
      pstd++;
      coord++;
      newcoord++;
    }

    if (mDebug)
      mLog << sqrt(norm) << " ";

    if (norm > 0) {
      jump -= 3;
      for (int ii = 0; ii < 3; ii++) {
        *jump /= norm;
        jump++;
      }
    }

    // Check that new control point is in mask
    if (!IsInMask(cpoint)) {
      *isrej = true;

      if (mDebug)
        mLog << "Reject due to control point " 
             << isrej - mRejectControl.begin() << " off mask at "
             << cpoint[0] << " " << cpoint[1] << " " << cpoint[2] << endl;
        LogObjectiveNaN(0);
    }

    cpoint += 3;
  }

  if (mDebug)
    mLog << endl;

  // Check that new end points are in end ROIs
  if (!IsInRoi(mControlPointsNew.begin(), mRoi1)) {
    mRejectControl[0] = true;

    if (mDebug) {
      cpoint = mControlPointsNew.begin();
      mLog << "Reject due to control point 0 off end ROI at "
           << cpoint[0] << " " << cpoint[1] << " " << cpoint[2] << endl;
      LogObjectiveNaN(0);
    }
  }

  if (!IsInRoi(mControlPointsNew.end() - 3, mRoi2)) {
    mRejectControl[mNumControl-1] = true;

    if (mDebug) {
      cpoint = mControlPointsNew.end() - 3;
      mLog << "Reject due to control point " << mNumControl-1
           << " off end ROI at "
           << cpoint[0] << " " << cpoint[1] << " " << cpoint[2] << endl;
      LogObjectiveNaN(0);
    }
  }

  for (isrej = mRejectControl.begin(); isrej < mRejectControl.end(); isrej++)
    if (*isrej)
      return false;

  // Check if new control points create zig-zag in path
  if (IsZigZag(mControlPointsNew, mControlPointsNew.begin(), 
                                  mControlPointsNew.end() - 3)) {
    mRejectSpline = true;

    if (mDebug) {
      mLog << "Reject due to zig-zag in spline" << endl;
      LogObjectiveNaN(0);
    }

    return false;
  }

  // Interpolate spline from proposed control points
  mSpline.SetControlPoints(mControlPointsNew);
  if (!mSpline.InterpolateSpline()) {
    mRejectSpline = true;

    if (mDebug) {
      mLog << "Reject due to spline" << endl;
      LogObjectiveNaN(0);
    }

    return false;
  }

  // Get proposed path points
  mPathPointsNew.clear();
  mPathPointsNew.insert(mPathPointsNew.begin(), mSpline.GetAllPointsBegin(),
                                                mSpline.GetAllPointsEnd());

  // Propagate proposed path to all time points
  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    if (!idwi->MapPathFromBase(mSpline)) {
      mRejectSpline = true;

      if (mDebug) {
        mLog << "Reject due to spline" << endl;
        LogObjectiveNaN(0);
      }

      return false;
    }

  return true;
}

//
// Propose path by perturbing a single control point
//
bool Coffin::ProposePathSingle(int ControlIndex) {
  const int offset = ControlIndex * 3;
  double norm = 0;
  vector<int>::const_iterator coord = mControlPoints.begin() + offset;
  vector<int>::const_iterator cpoint = mControlPointsNew.begin() + offset;
  vector<int>::iterator newcoord = mControlPointsNew.begin() + offset;
  vector<float>::const_iterator pstd = mProposalStd.begin() + offset;
  vector<float>::iterator jump = mControlPointJumps.begin() + offset;

  copy(mControlPoints.begin(), mControlPoints.end(), mControlPointsNew.begin());

  // Perturb current control point
  for (int ii = 0; ii < 3; ii++) {
    *jump = round((*pstd) * PDFgaussian());
    *newcoord = *coord + (int) *jump;

    *jump *= *jump;
    norm += *jump;

    jump++;
    pstd++;
    coord++;
    newcoord++;
  }

  if (norm > 0) {
    jump -= 3;
    for (int ii = 0; ii < 3; ii++) {
      *jump /= norm;
      jump++;
    }
  }

  // Check that new control point is in mask and, if end point, in end ROI
  if (!IsInMask(cpoint) ||
      ((ControlIndex == 0) && !IsInRoi(cpoint, mRoi1)) ||
      ((ControlIndex == mNumControl-1) && !IsInRoi(cpoint, mRoi2))) {
    mRejectControl[ControlIndex] = true;

    if (mDebug) {
      mLog << "Reject due to control point " << ControlIndex 
           << (IsInMask(cpoint) ? " off end ROI at " : " off mask at ")
           << cpoint[0] << " " << cpoint[1] << " " << cpoint[2] << endl;
      LogObjectiveNaN(0);
    }

    return false;
  }

  // Check if new control point creates zig-zag in path
  if (IsZigZag(mControlPointsNew, cpoint, cpoint)) {
    mRejectControl[ControlIndex] = true;

    if (mDebug) {
      mLog << "Reject due to zig-zag at control point " << ControlIndex << endl;
      LogObjectiveNaN(0);
    }

    return false;
  }

  // Interpolate spline from proposed control points
  mSpline.SetControlPoints(mControlPointsNew);
  if (!mSpline.InterpolateSpline()) {
    mRejectSpline = true;

    if (mDebug) {
      mLog << "Reject due to spline" << endl;
      LogObjectiveNaN(0);
    }

    return false;
  }
    
  // Get proposed path points
  mPathPointsNew.clear();
  mPathPointsNew.insert(mPathPointsNew.begin(), mSpline.GetAllPointsBegin(),
                                                mSpline.GetAllPointsEnd());

  // Propagate proposed path to all time points
  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    if (!idwi->MapPathFromBase(mSpline)) {
      mRejectSpline = true;

      if (mDebug) {
        mLog << "Reject due to spline" << endl;
        LogObjectiveNaN(0);
      }

      return false;
    }

  return true;
}

//
// Propose diffusion parameters along the path for all time points
//
void Coffin::ProposeDiffusionParameters() {
  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    idwi->ProposeDiffusionParameters();
}

//
// Determine if proposed path will be accepted
//
bool Coffin::AcceptPath(bool UsePriorOnly) {
  double neglogratio;
  vector<int> atlaspoints;
  vector<int>::iterator iptatlas;

  mDataPosteriorOnPathNew = 0;
  mDataPosteriorOffPathNew = 0;
  mDataPosteriorOnPath = 0;
  mDataPosteriorOffPath = 0;

  // Compute data-fit terms for all time points on proposed and current path
  if (!UsePriorOnly)
    for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end();
                                                     idwi++) {
      if (!idwi->ComputePathDataFit()) {
        mRejectF = idwi->RejectF();
        mAcceptF = idwi->AcceptF();
        mRejectTheta = idwi->RejectTheta();
        mAcceptTheta = idwi->AcceptTheta();

        if (mDebug) {
          const unsigned int itime = idwi - mDwi.begin();
          mLog << idwi->GetLog() << ", time point " << itime << endl;
          LogObjectiveNaN(itime);
        }

        return false;
      }

      mDataPosteriorOnPathNew  += idwi->GetPosteriorOnPathNew();
      mDataPosteriorOffPathNew += idwi->GetPosteriorOffPathNew();
      mDataPosteriorOnPath     += idwi->GetPosteriorOnPath();
      mDataPosteriorOffPath    += idwi->GetPosteriorOffPath();
    }

  // Map proposed path from diffusion/base space to atlas space
  atlaspoints.resize(mPathPointsNew.size());
  iptatlas = atlaspoints.begin();

  for (vector<int>::iterator ipt = mPathPointsNew.begin();
                             ipt < mPathPointsNew.end(); ipt += 3) {
    vector< vector<int> >::iterator icoord =
      mAtlasCoords.begin() + (ipt[0] + ipt[1]*mNx + ipt[2]*mNxy);

    if (icoord->empty()) {
      // Map point coordinates from diffusion/base space to atlas space
      MapPointToAtlas(iptatlas, ipt);

      // Save transformed coordinates for future use
      icoord->insert(icoord->end(), iptatlas, iptatlas+3);
    }
    else
      // Retrieve previously computed coordinates in atlas space
      copy(icoord->begin(), icoord->end(), iptatlas);

    iptatlas += 3;
  }

  // Compute atlas-derived prior terms on proposed path
  mXyzPriorOffPathNew = ComputeXyzPriorOffPath(atlaspoints);
  mXyzPriorOnPathNew  = ComputeXyzPriorOnPath(atlaspoints);

  mAnatomicalPriorNew = ComputeAnatomicalPrior(atlaspoints);

  mShapePriorNew = ComputeShapePrior(atlaspoints);

  // Compute posterior on proposed path
  mPosteriorOnPathNew  = mDataPosteriorOnPathNew + mXyzPriorOnPathNew
                       + mAnatomicalPriorNew + mShapePriorNew;
  mPosteriorOffPathNew = mDataPosteriorOffPathNew + mXyzPriorOffPathNew;

  // Compute posterior on current path
  mPosteriorOnPath  = mDataPosteriorOnPath + mXyzPriorOnPath
                    + mAnatomicalPrior + mShapePrior;
  mPosteriorOffPath = mDataPosteriorOffPath + mXyzPriorOffPath;

  // Compute ratio of proposed/current path posteriors
  neglogratio = mPosteriorOnPathNew - mPosteriorOffPathNew
              + mPosteriorOffPath   - mPosteriorOnPath;

  // Accept or reject proposed path based on ratio of posteriors
  if (drand48() < exp(-neglogratio)) {
    if (mDebug) {
      mLog << "Accept due to posterior (alpha = " << exp(-neglogratio) << ")"
           << endl;
      LogObjective();
    }

    return true;
  }

  mRejectPosterior = true;

  if (mDebug) {
    mLog << "Reject due to posterior (alpha = " << exp(-neglogratio) << ")"
         << endl;
    LogObjective();
  }

  return false;
}
 
//
// Compute off-path spatial prior for a set of points
//
double Coffin::ComputeXyzPriorOffPath(std::vector<int> &PathAtlasPoints) {
  double prior = 0;

  if (!mXyzPrior0)
    return 0;

  for (vector<int>::const_iterator ipt = PathAtlasPoints.begin();
                                   ipt < PathAtlasPoints.end(); ipt += 3)
    prior += (double) MRIgetVoxVal(mXyzPrior0, ipt[0], ipt[1], ipt[2], 0);

  return prior / (PathAtlasPoints.size()/3);
}

//
// Compute on-path spatial prior for a set of points
//
double Coffin::ComputeXyzPriorOnPath(std::vector<int> &PathAtlasPoints) {
  double prior = 0;

  if (!mXyzPrior1)
    return 0;

  for (vector<int>::const_iterator ipt = PathAtlasPoints.begin();
                                   ipt < PathAtlasPoints.end(); ipt += 3)
    prior += (double) MRIgetVoxVal(mXyzPrior1, ipt[0], ipt[1], ipt[2], 0);

  return prior / (PathAtlasPoints.size()/3);
}

//
// Compute prior on path given its tangent vector and curvature
//
double Coffin::ComputeShapePrior(vector<int> &PathAtlasPoints) {
  const int nbin = (int) ceil(2 / mTangentBinSize);
  const double darc = mNumArc / (double) (PathAtlasPoints.size()/3);
  double larc = 0, prior = 0;
  vector<float>::const_iterator id1, id2;
  vector< vector<float> >::const_iterator iprtang = mPriorTangent.begin(),
                                          iprcurv = mPriorCurvature.begin();
  vector<float> pathsmooth(PathAtlasPoints.size()),
                diff1(PathAtlasPoints.size()),
                diff2(PathAtlasPoints.size());

  if (mPriorTangent.empty() && mPriorCurvature.empty())
    return 0;

// FILL THE PATH FOR ALL PRIORS?
//vector<int> pathfilled = CurveFill(PathAtlasPoints);

  // Smooth discrete point coordinates
  CurveSmooth(pathsmooth, PathAtlasPoints);

  // Approximate first derivative by smoothed finite differences
  // of the point coordinates
  CurveFiniteDifferences(diff1, pathsmooth, mDiffStep);

  // Approximate second derivative by smoothed finite differences
  // of the first derivative
  CurveFiniteDifferences(diff2, diff1, mDiffStep);

  id2 = diff2.begin();

  for (id1 = diff1.begin(); id1 < diff1.end(); id1 += 3) {
    int ix, iy;
    float tangx, tangy, curv;
    const float nrm = sqrt(id1[0]*id1[0] + id1[1]*id1[1] + id1[2]*id1[2]);

    if (nrm > 0) {
      // Tangent vector
      if (id1[2] < 0) {		// Use only z>0 octants, direction is irrelevant
        tangx = - id1[0] / nrm;
        tangy = - id1[1] / nrm;
      }
      else {
        tangx = id1[0] / nrm;
        tangy = id1[1] / nrm;
      }

      // Curvature = |r' x r''| / |r'|^3
      curv =
        sqrt( pow(id1[1] * id2[2] - id1[2] * id2[1], 2) +
              pow(id1[2] * id2[0] - id1[0] * id2[2], 2) +
              pow(id1[0] * id2[1] - id1[1] * id2[0], 2) ) / pow(nrm, 3);
    }
    else {
      tangx = 0;
      tangy = 0;
      curv = 0;
    }

    // Find prior given tangent vector
    ix = (int) floor((tangx + 1) / mTangentBinSize);
    iy = (int) floor((tangy + 1) / mTangentBinSize);

    prior += iprtang->at( ((iy<nbin) ? iy : (nbin-1)) * nbin +
                          ((ix<nbin) ? ix : (nbin-1)) );

    // Find prior given curvature
    ix = (int) floor(curv / mCurvatureBinSize);

    prior += iprcurv->at( (ix < (int) iprcurv->size()) ?
                          ix : ((int) iprcurv->size()-1) );

    id2 += 3;

    larc += darc;

    if (larc > 1)  {    // Move to the next segment
      larc -= 1;
      iprtang++;
      iprcurv++;
    }
  }

  return prior / (PathAtlasPoints.size()/3);
}

//
// Compute prior on path given anatomical segmentation labels around path
//
double Coffin::ComputeAnatomicalPrior(vector<int> &PathAtlasPoints) {
  const double darc = mNumArc / (double) (PathAtlasPoints.size()/3);
  double larc = 0, prior = 0;
  vector<float>::iterator iseg0;
  vector<unsigned int>::const_iterator imatch;
  vector< vector<unsigned int> >::const_iterator iidlocal = mIdsLocal.begin(),
                                                 iidnear  = mIdsNear.begin(),
                                                 iid;
  vector< vector<float> >::const_iterator iprlocal = mPriorLocal.begin(), 
                                          iprnear  = mPriorNear.begin(),
                                          ipr;
  vector<float> seg0(mAseg.size());

  if (mPriorLocal.empty() && mPriorNear.empty())
    return 0;

  for (vector<int>::const_iterator ipt = PathAtlasPoints.begin();
                                   ipt < PathAtlasPoints.end(); ipt += 3) {
    const int ix0 = ipt[0], iy0 = ipt[1], iz0 = ipt[2];

    // Find prior given local neighbor labels
    iid = iidlocal;
    ipr = iprlocal;

    for (vector<int>::const_iterator idir = mDirLocal.begin();
                                     idir != mDirLocal.end(); idir += 3) {
      const int ix = ix0 + idir[0],
                iy = iy0 + idir[1],
                iz = iz0 + idir[2];

      for (vector<MRI *>::const_iterator iaseg = mAseg.begin();
                                         iaseg < mAseg.end(); iaseg++) {
        imatch = find(iid->begin(), iid->end(),
                      (unsigned int) MRIgetVoxVal(*iaseg,
                                     ((ix > -1 && ix < mNxAtlas) ? ix : ix0),
                                     ((iy > -1 && iy < mNyAtlas) ? iy : iy0),
                                     ((iz > -1 && iz < mNzAtlas) ? iz : iz0),
                                     0));

        if (imatch < iid->end())
          prior += ipr->at(imatch - iid->begin());
        else
          prior += *(ipr->end() - 1);
      }

      iid += mNumArc;
      ipr += mNumArc;
    }

    // Find prior given nearest neighbor labels
    iid = iidnear;
    ipr = iprnear;

    iseg0 = seg0.begin();
    for (vector<MRI *>::const_iterator iaseg = mAseg.begin();
                                       iaseg < mAseg.end(); iaseg++) {
      *iseg0 = MRIgetVoxVal(*iaseg, ix0, iy0, iz0, 0);
      iseg0++;
    }

    for (vector<int>::const_iterator idir = mDirNear.begin();
                                     idir != mDirNear.end(); idir += 3) {
      int dist = 0, ix = ix0 + idir[0],
                    iy = iy0 + idir[1],
                    iz = iz0 + idir[2];

      iseg0 = seg0.begin();
      for (vector<MRI *>::const_iterator iaseg = mAseg.begin();
                                         iaseg < mAseg.end(); iaseg++) {
        float seg = *iseg0;

        while ((ix > -1) && (ix < mNxAtlas) &&
               (iy > -1) && (iy < mNyAtlas) &&
               (iz > -1) && (iz < mNzAtlas) && (seg == *iseg0)) {
          seg = MRIgetVoxVal(*iaseg, ix, iy, iz, 0);
          dist++;

          ix += idir[0];
          iy += idir[1];
          iz += idir[2];
        }

        imatch = find(iid->begin(), iid->end(), (unsigned int) seg);

        if (imatch < iid->end())
          prior += ipr->at(imatch - iid->begin());
        else
          prior += *(ipr->end() - 1);

        iseg0++;
      }

      iid += mNumArc;
      ipr += mNumArc;
    }

    larc += darc;

    if (larc > 1)  {	// Move to the next segment
      larc -= 1;
      iidlocal++;
      iprlocal++;
      iidnear++;
      iprnear++;
    }
  }

  return prior / (PathAtlasPoints.size()/3);
}

//
// Copy newly accepted path over current path for all time points and base
//
void Coffin::UpdatePath() {
  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    idwi->UpdatePath();

  copy(mControlPointsNew.begin(), mControlPointsNew.end(), 
       mControlPoints.begin());

  mPathPoints.resize(mPathPointsNew.size());
  copy(mPathPointsNew.begin(), mPathPointsNew.end(), mPathPoints.begin());

  mDataPosteriorOnPath  = mDataPosteriorOnPathNew; 
  mDataPosteriorOffPath = mDataPosteriorOffPathNew; 
  mXyzPriorOnPath       = mXyzPriorOnPathNew;
  mXyzPriorOffPath      = mXyzPriorOffPathNew;
  mAnatomicalPrior      = mAnatomicalPriorNew; 
  mShapePrior           = mShapePriorNew; 
  mPosteriorOnPath      = mPosteriorOnPathNew;
  mPosteriorOffPath     = mPosteriorOffPathNew;
}

//
// Update control point acceptance rate (for full spline updates)
//
void Coffin::UpdateAcceptanceRateFull() {
  vector<float>::const_iterator jump = mControlPointJumps.begin();
  
  if (!(mAcceptF || mAcceptTheta))
    // Update length of accepted jumps for all control points
    for (vector<float>::iterator span = mAcceptSpan.begin();
                                 span < mAcceptSpan.end(); span++) {
      *span += *jump;
      jump++;
    }
}

//
// Update control point rejection rate (for full spline updates)
//
void Coffin::UpdateRejectionRateFull() {
  vector<float>::const_iterator jump = mControlPointJumps.begin();

  if (mRejectPosterior || mRejectSpline)
    // Update length of rejected jumps for all control points
    for (vector<float>::iterator span = mRejectSpan.begin();
                                 span < mRejectSpan.end(); span++) {
      *span += *jump;
      jump++;
    }
  else {
    vector<float>::iterator span = mRejectSpan.begin();

    // Update length of rejected jumps for specific rejected control points
    for (vector<bool>::const_iterator isrej = mRejectControl.begin();
                                      isrej < mRejectControl.end(); isrej++)
      if (*isrej)
        for (int ii = 0; ii < 3; ii++) {
          *span += *jump;
          jump++;
          span++;
        }
      else {
        jump += 3;
        span += 3;
      }
  }
}

//
// Update control point acceptance/rejection rates (for single control updates)
//
void Coffin::UpdateAcceptRejectRateSingle() {
  vector<int>::iterator acount = mAcceptCount.begin();
  vector<int>::iterator rcount = mRejectCount.begin();
  vector<float>::const_iterator jump = mControlPointJumps.begin();
  vector<float>::iterator aspan = mAcceptSpan.begin();
  vector<float>::iterator rspan = mRejectSpan.begin();

  for (vector<bool>::const_iterator isrej = mRejectControl.begin();
                                    isrej < mRejectControl.end(); isrej++)
    if (*isrej) {	// Increment rejection count and jump length
      for (int ii = 0; ii < 3; ii++) {
        *rspan += *jump;
        jump++;
        aspan++;
        rspan++;
      }

      (*rcount)++;
      acount++;
      rcount++;
    }
    else {		// Increment acceptance count and jump length
      for (int ii = 0; ii < 3; ii++) {
        *aspan += *jump;
        jump++;
        aspan++;
        rspan++;
      }

      (*acount)++;
      acount++;
      rcount++;
    }

  if (mDebug) {
    mLog << "Normalized jumps: ";
    for (jump = mControlPointJumps.begin(); jump < mControlPointJumps.end();
                                            jump++)
      mLog << sqrt(*jump) << " ";
    mLog << endl;

    acount = mAcceptCount.begin();
    rcount = mRejectCount.begin();

    mLog << "Acceptance counts: ";
    for (vector<bool>::const_iterator isrej = mRejectControl.begin();
                                      isrej < mRejectControl.end(); isrej++) {
      mLog << *acount << "/" << *acount + *rcount << "="
           << *acount / float(*acount + *rcount) << " ";
      acount++;
      rcount++;
    }
    mLog << endl;

    aspan = mAcceptSpan.begin();
    rspan = mRejectSpan.begin();

    mLog << "Acceptance spans: ";
    for (vector<float>::const_iterator pstd = mProposalStd.begin();
                                       pstd < mProposalStd.end(); pstd++) {
      mLog << *aspan << "/" << *aspan + *rspan << "="
           << *aspan / (*aspan + *rspan) << " ";

      aspan++;
      rspan++;
    }
    mLog << endl;
  }
}

//
// Update control point proposal distribution
//
void Coffin::UpdateProposalStd() {
  vector<float>::const_iterator aspan = mAcceptSpan.begin();
  vector<float>::const_iterator rspan = mRejectSpan.begin();

  for (vector<float>::iterator pstd = mProposalStd.begin();
                               pstd < mProposalStd.end(); pstd++) {
    *pstd *= sqrt((*aspan + 1) / (*rspan + 1));

    aspan++;
    rspan++;
  }

  if (mDebug) {
    mLog << "Proposal STDs: ";
    for (vector<float>::const_iterator pstd = mProposalStd.begin();
                                       pstd < mProposalStd.end(); pstd++)
      mLog << *pstd << " ";
    mLog << endl;
  }

  fill(mAcceptCount.begin(), mAcceptCount.end(), 0);
  fill(mRejectCount.begin(), mRejectCount.end(), 0);

  fill(mAcceptSpan.begin(), mAcceptSpan.end(), 0.0);
  fill(mRejectSpan.begin(), mRejectSpan.end(), 0.0);
}

//
// Save data-fit and prior terms of accepted and rejected path
// for all time points
//
void Coffin::SavePathPosterior(bool IsPathAccepted) {
  vector<float> priors(6);

  // Save path data-fit terms for each time point
  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    idwi->SavePathDataFit(IsPathAccepted);

  // Save atlas-derived path prior terms, common among all time points
  if (IsPathAccepted) {			// The newly sampled path was accepted
    priors[0] = (float) mXyzPriorOnPathNew - mXyzPriorOffPathNew;
    priors[1] = (float) mAnatomicalPriorNew;
    priors[2] = (float) mShapePriorNew;
    priors[3] = (float) mXyzPriorOnPath - mXyzPriorOffPath;
    priors[4] = (float) mAnatomicalPrior;
    priors[5] = (float) mShapePrior;
  }
  else {				// The newly sampled path was rejected
    priors[0] = (float) mXyzPriorOnPath - mXyzPriorOffPath;
    priors[1] = (float) mAnatomicalPrior;
    priors[2] = (float) mShapePrior;
    priors[3] = (float) mXyzPriorOnPathNew - mXyzPriorOffPathNew;
    priors[4] = (float) mAnatomicalPriorNew;
    priors[5] = (float) mShapePriorNew;
  }

  Aeon::SavePathPriors(priors);
}

//
// Save current path as an MCMC sample
//
void Coffin::SavePath() {
  // If in longitudinal mode, downsample from base to native space
  if (mDwi[0].GetBaseMask())
    RemoveDuplicatePathPoints();

  // Save current path for each time point
  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    idwi->SavePath();

  // If in longitudinal mode, also save current path in base space
  if (mDwi[0].GetBaseMask())
    Aeon::SaveBasePath(mPathPoints);

  // Keep track of MAP path
  if (mPosteriorOnPath < mPosteriorOnPathMap) {
    Aeon::SetPathMap(mDwi[0].GetNumSample() - 1);
    mPosteriorOnPathMap = mPosteriorOnPath;
  }
}

//
// Remove duplicate points from path
// Duplicate points can arise when mapping a path from base space to
// a lower-resolution native space 
//
void Coffin::RemoveDuplicatePathPoints() {
  unsigned int newsize = 0;
  vector<int>::const_iterator ipt;
  vector<int>::iterator iptnew;
  vector<bool> isdup(mPathPoints.size()/3, false);
  vector<int> newpath;

  // Find duplicate path points
  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    idwi->FindDuplicatePathPoints(isdup);

  // Remove duplicate path points in base space
  for (vector<bool>::const_iterator idup = isdup.begin(); idup < isdup.end();
                                                          idup++)
    if (! *idup)
      newsize += 3;

  newpath.resize(newsize);

  ipt = mPathPoints.begin();
  iptnew = newpath.begin();

  for (vector<bool>::const_iterator idup = isdup.begin(); idup < isdup.end();
                                                          idup++) {
    if (! *idup) {
      copy(ipt, ipt+3, iptnew);
      iptnew += 3;
    }

    ipt += 3;
  }

  mPathPoints.resize(newsize);
  copy(newpath.begin(), newpath.end(), mPathPoints.begin());

  // Remove duplicate path points in the space of each time point
  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++)
    idwi->RemovePathPoints(isdup, newsize);
}

// Check that a point is inside the mask
//
bool Coffin::IsInMask(vector<int>::const_iterator Point) {
  return (Point[0] > -1) && (Point[0] < mMask->width) &&
         (Point[1] > -1) && (Point[1] < mMask->height) &&
         (Point[2] > -1) && (Point[2] < mMask->depth) &&
         (MRIgetVoxVal(mMask, Point[0], Point[1], Point[2], 0) > 0);
}

//
// Check that a point is inside an ROI
//
bool Coffin::IsInRoi(vector<int>::const_iterator Point, MRI *Roi) {
  vector<int> outpoint(3);

  MapPointToAtlas(outpoint.begin(), Point);

  return (outpoint[0] > -1) && (outpoint[0] < Roi->width) &&
         (outpoint[1] > -1) && (outpoint[1] < Roi->height) &&
         (outpoint[2] > -1) && (outpoint[2] < Roi->depth) &&
         (MRIgetVoxVal(Roi, outpoint[0], outpoint[1], outpoint[2], 0) > 0);
}

//
// Check if perturbing a control point (or series of control points)
// would cause a zig-zag in the path:
// A zig-zag is detected as two consecutive acute angles between path segments
//
bool Coffin::IsZigZag(vector<int> &ControlPoints,
                      vector<int>::const_iterator FirstPerturbedPoint,
                      vector<int>::const_iterator LastPerturbedPoint) {
  vector<int>::const_iterator curpoint, checkfirst, checklast;

  // The zig-zag test will also be positive if a perturbed control point
  // overlaps with one of its two neighbors, in all but these two cases
  // (when the overlap is between the first two or last two control points)
  if ( ( FirstPerturbedPoint <= ControlPoints.begin() + 3 &&
         *ControlPoints.begin()       == *(ControlPoints.begin() + 3) &&
         *(ControlPoints.begin() + 1) == *(ControlPoints.begin() + 4) &&
         *(ControlPoints.begin() + 2) == *(ControlPoints.begin() + 5) ) ||
       ( LastPerturbedPoint >= ControlPoints.end() - 6 &&
         *(ControlPoints.end() - 6) == *(ControlPoints.end() - 3) &&
         *(ControlPoints.end() - 5) == *(ControlPoints.end() - 2) &&
         *(ControlPoints.end() - 4) == *(ControlPoints.end() - 1) ) )
    return true;

  // When a single control point is perturbed, there are
  // up to 4 neighboring spots that need to be checked for a potential zig-zag
  // First spot to check:
  curpoint = FirstPerturbedPoint - 3;
  checkfirst = ControlPoints.begin() + 6;
  if (curpoint > checkfirst)
    checkfirst = curpoint;

  // Last spot to check:
  curpoint = LastPerturbedPoint + 6;
  checklast = ControlPoints.end() - 6;
  if (curpoint < checklast)
    checklast = curpoint;

  if (checkfirst > checklast)	// If spline has fewer than 4 control points
    return false;
  else {  
    int dot1, dot2;
    vector<int>::const_iterator x1, x2, x3, x4;
    vector<int> diff12(3), diff32(3), diff43(3);

    // Check at first zig-zag candidate spot
    curpoint = checkfirst;

    x1 = curpoint - 6;
    x2 = curpoint - 3;
    x3 = curpoint;
    x4 = curpoint + 3;

    dot1 = 0; dot2 = 0;

    for (int k = 0; k < 3; k++) {
      diff12[k] = x1[k] - x2[k];
      diff32[k] = x3[k] - x2[k];
      diff43[k] = x4[k] - x3[k];

      dot1 += diff12[k] * diff32[k];
      dot2 += diff43[k] * (-diff32[k]);
    }

    if (dot1 >= 0 && dot2 >= 0)		// If both angles are acute
      return true;

    // Check at remaining zig-zag candidate spots
    for (curpoint = checkfirst + 3; curpoint <= checklast; curpoint += 3) {
      x3 += 3;
      x4 += 3;

      dot1 = 0; dot2 = 0;

      for (int k = 0; k < 3; k++) {
        diff12[k] = -diff32[k];
        diff32[k] = diff43[k];
        diff43[k] = x4[k] - x3[k];

        dot1 += diff12[k] * diff32[k];
        dot2 += diff43[k] * (-diff32[k]);
      }

      if (dot1 >= 0 && dot2 >= 0)	// If both angles are acute
        return true;
    }

    return false;
  }
}

//
// Map point coordinates from diffusion/base space to atlas space
//
void Coffin::MapPointToAtlas(vector<int>::iterator OutPoint,
                             vector<int>::const_iterator InPoint) {
  if (!mAffineReg.IsEmpty()) {
    vector<float> point(InPoint, InPoint+3);

    mAffineReg.ApplyXfm(point, point.begin());
#ifndef NO_CVS_UP_IN_HERE
    if (!mNonlinReg.IsEmpty())
      mNonlinReg.ApplyXfm(point, point.begin());
#endif

    for (int k = 0; k < 3; k++)
      OutPoint[k] = (int) round(point[k]);
  }
  else
    for (int k = 0; k < 3; k++)
      OutPoint[k] = InPoint[k];
}

//
// Write terms of objective function to log file (all terms are defined)
//
void Coffin::LogObjective() {
  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOnPathNew=" << idwi->GetLikelihoodOnPathNew()
         << " mPriorOnPathNew=" << idwi->GetPriorOnPathNew() << endl;

  mLog << "mXyzPriorOnPathNew=" << mXyzPriorOnPathNew << endl
       << "mAnatomicalPriorNew=" << mAnatomicalPriorNew << endl
       << "mShapePriorNew=" << mShapePriorNew << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOffPathNew=" << idwi->GetLikelihoodOffPathNew()
         << " mPriorOffPathNew=" << idwi->GetPriorOffPathNew() << endl;

  mLog << "mXyzPriorOffPathNew=" << mXyzPriorOffPathNew << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOn-OffPathNew="
         << idwi->GetLikelihoodOnPathNew() - idwi->GetLikelihoodOffPathNew()
         << " mPriorOn-OffPathNew="
         << idwi->GetPriorOnPathNew() - idwi->GetPriorOffPathNew() << endl;

  mLog << "mXyzPriorOn-OffPathNew="
       << mXyzPriorOnPathNew - mXyzPriorOffPathNew << endl
       << "mPathLengthNew=" << mPathPointsNew.size()/3 << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOnPath=" << idwi->GetLikelihoodOnPath()
         << " mPriorOnPath=" << idwi->GetPriorOnPath() << endl;

  mLog << "mXyzPriorOnPath=" << mXyzPriorOnPath << endl
       << "mAnatomicalPrior=" << mAnatomicalPrior << endl
       << "mShapePrior=" << mShapePrior << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOffPath=" << idwi->GetLikelihoodOffPath()
         << " mPriorOffPath=" << idwi->GetPriorOffPath() << endl;

  mLog << "mXyzPriorOffPath=" << mXyzPriorOffPath << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOff-OnPath="
         << idwi->GetLikelihoodOffPath() - idwi->GetLikelihoodOnPath()
         << " mPriorOff-OnPath="
         << idwi->GetPriorOffPath() - idwi->GetPriorOnPath() << endl;

  mLog << "mXyzPriorOn-OffPath=" << mXyzPriorOnPath-mXyzPriorOffPath << endl
       << "mPathLength=" << mPathPoints.size()/3 << endl;
}

//
// Write terms of objective function to log file (some terms are not defined
// because the new path was rejected or accepted before some of the terms were
// computed)
//
void Coffin::LogObjectiveNaN(const unsigned int NumValidData) {
  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.begin() + NumValidData; idwi++)
    mLog << "mLikelihoodOnPathNew=" << idwi->GetLikelihoodOnPathNew()
         << " mPriorOnPathNew=" << idwi->GetPriorOnPathNew() << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin() + NumValidData;
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOnPathNew=NaN mPriorOnPathNew=NaN" << endl;

  mLog << "mXyzPriorOnPathNew=NaN" << endl
       << "mAnatomicalPriorNew=NaN" << endl
       << "mShapePriorNew=NaN" << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.begin() + NumValidData; idwi++)
    mLog << "mLikelihoodOffPathNew=" << idwi->GetLikelihoodOffPathNew()
         << " mPriorOffPathNew=" << idwi->GetPriorOffPathNew() << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin() + NumValidData;
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOffPathNew=NaN mPriorOffPathNew=NaN" << endl;

  mLog << "mXyzPriorOffPathNew=NaN" << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.begin() + NumValidData; idwi++)
    mLog << "mLikelihoodOn-OffPathNew="
         << idwi->GetLikelihoodOnPathNew() - idwi->GetLikelihoodOffPathNew()
         << " mPriorOn-OffPathNew="
         << idwi->GetPriorOnPathNew() - idwi->GetPriorOffPathNew() << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin() + NumValidData;
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOn-OffPathNew=NaN mPriorOn-OffPathNew=NaN" << endl;

  mLog << "mXyzPriorOn-OffPathNew=NaN" << endl
       << "mPathLengthNew=" << mPathPointsNew.size()/3 << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.begin() + NumValidData; idwi++)
    mLog << "mLikelihoodOnPath=" << idwi->GetLikelihoodOnPath()
         << " mPriorOnPath=" << idwi->GetPriorOnPath() << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin() + NumValidData;
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOnPath=NaN mPriorOnPath=NaN" << endl;

  mLog << "mXyzPriorOnPath=" << mXyzPriorOnPath << endl
       << "mAnatomicalPrior=" << mAnatomicalPrior << endl
       << "mShapePrior=" << mShapePrior << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.begin() + NumValidData; idwi++)
    mLog << "mLikelihoodOffPath=" << idwi->GetLikelihoodOffPath()
         << " mPriorOffPath=" << idwi->GetPriorOffPath() << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin() + NumValidData;
                                  idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOffPath=NaN mPriorOffPath=NaN" << endl;

  mLog << "mXyzPriorOffPath=" << mXyzPriorOffPath << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin();
                                    idwi < mDwi.begin() + NumValidData; idwi++)
    mLog << "mLikelihoodOff-OnPath="
         << idwi->GetLikelihoodOffPath() - idwi->GetLikelihoodOnPath()
         << " mPriorOff-OnPath="
         << idwi->GetPriorOffPath() - idwi->GetPriorOnPath() << endl;

  for (vector<Aeon>::const_iterator idwi = mDwi.begin() + NumValidData;
                                    idwi < mDwi.end(); idwi++)
    mLog << "mLikelihoodOff-OnPath=NaN mPriorOff-OnPath=NaN" << endl;

  mLog << "mXyzPriorOn-OffPath=" << mXyzPriorOnPath-mXyzPriorOffPath << endl
       << "mPathLength=" << mPathPoints.size()/3 << endl;
}

//
// Check that a point is inside this time point's mask
//
bool Aeon::IsInMask(vector<int>::const_iterator Point) {
  return (Point[0] > -1) && (Point[0] < mMask->width) &&
         (Point[1] > -1) && (Point[1] < mMask->height) &&
         (Point[2] > -1) && (Point[2] < mMask->depth) &&
         (MRIgetVoxVal(mMask, Point[0], Point[1], Point[2], 0) > 0);
}

//
// Compute leengths of path samples
//
void Aeon::ComputePathLengths(vector<int> &PathLengths,
                              vector< vector<int> > &PathSamples) {
  vector<int>::iterator ilen = PathLengths.begin();

  for (vector< vector<int> >::const_iterator ipath = PathSamples.begin();
                                             ipath < PathSamples.end();
                                             ipath++) {
    *ilen = ipath->size() / 3;

    ilen++;
  }
}

//
// Compute spatial histogram from path samples
//
void Aeon::ComputePathHisto(MRI *HistoVol,
                            vector< vector<int> > &PathSamples) {
  for (vector< vector<int> >::const_iterator ipath = PathSamples.begin();
                                             ipath < PathSamples.end();
                                             ipath++)
    for (vector<int>::const_iterator ipt = ipath->begin();
                                     ipt < ipath->end(); ipt += 3) {
      const int ix = ipt[0], iy = ipt[1], iz = ipt[2];
      
      MRIsetVoxVal(HistoVol, ix, iy, iz, 0,
                   MRIgetVoxVal(HistoVol, ix, iy, iz, 0) + 1);
    }
}

//
// Find maximum a posteriori path
//
int Aeon::FindMaxAPosterioriPath(vector< vector<int> > &PathSamples,
                                 vector<int> &PathLengths, MRI *PathHisto) {
  const int lmin = *min_element(PathLengths.begin(), PathLengths.end()),
            lmax = *max_element(PathLengths.begin(), PathLengths.end());
  float lnorm, pathnorm = 0.0, probmax = 0.0;
  vector<int>::const_iterator ilen;
  vector< vector<int> >::const_iterator ipathmap;
  vector<float> lhisto(lmax-lmin+1, 0), lhistofilt(lhisto.size(), 0);
  vector<float>::iterator ihistofilt = lhistofilt.begin();

  // Find histogram of lengths
  for (ilen = PathLengths.begin(); ilen < PathLengths.end(); ilen++)
    lhisto.at(*ilen - lmin)++;

  // Smooth histogram of lengths
  *ihistofilt = (lhisto[0] + lhisto[1]) / 2.0;
  ihistofilt++;

  for (vector<float>::const_iterator ihisto = lhisto.begin() + 1;
                                     ihisto < lhisto.end() - 1; ihisto++) {
    *ihistofilt = (*(ihisto-1) + *ihisto + *(ihisto+1)) / 3.0;
    ihistofilt++;
  }

  *ihistofilt = (*(lhisto.end()-2) + *(lhisto.end()-1)) / 2.0;

  // Find normalization of length posterior to a sum of 1
  lnorm = accumulate(lhistofilt.begin(), lhistofilt.end(), 0.0);

  // Find normalization of path posterior to a sum of 1
  for (int iz = 0; iz < mNz; iz++)
    for (int iy = 0; iy < mNy; iy++)
      for (int ix = 0; ix < mNx; ix++)
        pathnorm += MRIgetVoxVal(PathHisto, ix, iy, iz, 0);

  // Find maximum a posteriori path
  ilen = PathLengths.begin();
  for (vector< vector<int> >::const_iterator ipath = PathSamples.begin();
                                             ipath < PathSamples.end();
                                             ipath++) {
    float prob = 0.0;

    // Path probability
    for (vector<int>::const_iterator ipt = ipath->begin();
                                     ipt < ipath->end(); ipt += 3)
      prob += MRIgetVoxVal(PathHisto, ipt[0], ipt[1], ipt[2], 0);

    prob /= pathnorm;
    prob /= (*ilen);		 // Bin paths by length eventually?

    // Length probability
    prob += lhistofilt.at(*ilen - lmin) / lnorm;

    if (prob > probmax) {
      ipathmap = ipath;
      probmax = prob;
    }

    ilen++;
  }

  // Return index of maximum a posteriori path
  return ipathmap - PathSamples.begin();
}

//
// Write output files for all time points
//
void Coffin::WriteOutputs() {
  string fname;

  for (vector<Aeon>::iterator idwi = mDwi.begin(); idwi < mDwi.end(); idwi++) {
    fname = idwi->GetOutputDir() + "/log.txt";
    mLog.open(fname.c_str(), ios::out | ios::app);
    if (!mLog) {
      cout << "ERROR: Could not open " << fname << " for writing" << endl;
      exit(1);
    }

    cout << "Writing output files to " << idwi->GetOutputDir() << endl;
    mLog << "Writing output files to " << idwi->GetOutputDir() << endl;

    idwi->WriteOutputs();

    mLog.flush();
    mLog.close();
  }
}

