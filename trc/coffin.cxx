/**
 * @file  coffin.cxx
 * @brief Container of tractography data and functions
 *
 * Container of tractography data and functions
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

#include <coffin.h>

using namespace std;

Coffin::Coffin(const char *OutDir, const char *DwiFile,
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
        const bool Debug) :
        mDebug(Debug),
        mNumBurnIn(NumBurnIn), mNumSample(NumSample),
        mKeepSampleNth(KeepSampleNth), mUpdatePropNth(UpdatePropNth),
        mSpline(InitFile, MaskFile) {
  char fname[PATH_MAX];
  MRI *dwi, *phi[NumTract], *theta[NumTract], *f[NumTract],
      *v0[NumTract], *f0[NumTract], *d0,
      *prior0 = NULL, *prior1 = NULL,
      *aseg = NULL, *asegtr = NULL, *pathtr = NULL;
  vector<MRI *> segprior0, segprior1;
  vector< vector<unsigned int> > segids;

  // Make output directory
  cout << "Creating output directory " << OutDir << endl;
  sprintf(fname, "mkdir -p %s", OutDir);
  if (system(fname) != 0) {
    cout << "ERROR: Could not create directory " << OutDir << endl;
    exit(1);
  }

  mOutDir = (char *) calloc(strlen(OutDir)+1, sizeof(char));
  strcpy(mOutDir, OutDir);

  // Read DWI-to-ROI registration matrix
  ReadDwiToRoi(XfmFile);

  // Read control point initialization
  ReadControlPoints(InitFile);

  // Read proposal SD initialization
  ReadProposalStds(PropStdFile);
  if (mProposalStdInit.size() != mControlPoints.size()) {
    cout << "ERROR: Dimensions of " << PropStdFile << " and " << InitFile
         << " do not match" << endl;
    exit(1);
  }

  // Read diffusion-weighted images
  cout << "Loading DWIs from " << DwiFile << endl;
  dwi = MRIread(DwiFile);
  if (!dwi) {
    cout << "ERROR: Could not read " << DwiFile << endl;
    exit(1);
  }

  // Size of diffusion-weighted images
  mNx = (unsigned int) dwi->width;
  mNy = (unsigned int) dwi->height;
  mNz = (unsigned int) dwi->depth;
  mNxy = mNx * mNy;

  // Resolution of diffusion-weighted images
  mDwiVoxelSize.clear();
  mDwiVoxelSize.push_back(dwi->xsize);
  mDwiVoxelSize.push_back(dwi->ysize);
  mDwiVoxelSize.push_back(dwi->zsize);

  // Read mask
  cout << "Loading mask from " << MaskFile << endl;
  mMask = MRIread(MaskFile);
  if (!mMask) {
    cout << "ERROR: Could not read " << MaskFile << endl;
    exit(1);
  }

  mPathSamples = MRIclone(mMask, NULL);

  // Read parameter samples from BEDPOST directory
  cout << "Loading BEDPOST parameter samples from " << BedpostDir << endl;
  for (unsigned int itract = 0; itract < NumTract; itract++) {
    sprintf(fname, "%s/merged_ph%usamples.nii.gz", BedpostDir, itract+1);
    phi[itract] = MRIread(fname);
    if (!phi[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
    sprintf(fname, "%s/merged_th%usamples.nii.gz", BedpostDir, itract+1);
    theta[itract] = MRIread(fname);
    if (!theta[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
    sprintf(fname, "%s/merged_f%usamples.nii.gz", BedpostDir, itract+1);
    f[itract] = MRIread(fname);
    if (!f[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
    sprintf(fname, "%s/dyads%u.nii.gz", BedpostDir, itract+1);
    v0[itract] = MRIread(fname);
    if (!v0[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
    sprintf(fname, "%s/mean_f%usamples.nii.gz", BedpostDir, itract+1);
    f0[itract] = MRIread(fname);
    if (!f0[itract]) {
      cout << "ERROR: Could not read " << fname << endl;
      exit(1);
    }
  }

  sprintf(fname, "%s/mean_dsamples.nii.gz", BedpostDir);
  d0 = MRIread(fname);
  if (!d0) {
    cout << "ERROR: Could not read " << fname << endl;
    exit(1);
  }

  // Read spatial path priors
  if (PriorFile0 && PriorFile1) {
    cout << "Loading spatial path prior from " << PriorFile0 << endl;
    prior0 = MRIread(PriorFile0);
    if (!prior0) {
      cout << "ERROR: Could not read " << PriorFile0 << endl;
      exit(1);
    }

    cout << "Loading spatial path prior from " << PriorFile1 << endl;
    prior1 = MRIread(PriorFile1);
    if (!prior1) {
      cout << "ERROR: Could not read " << PriorFile1 << endl;
      exit(1);
    }
  }

  // Read aseg priors
  if (AsegPriorFile0 && AsegPriorFile1 && AsegIdFile) {
    MRI *tmp;
    vector<string> suff;
    vector<unsigned int> idlist;

    if (FileExists(AsegPriorFile0)) {
      suff.resize(1);
      suff.push_back("");
    }
    else {
      suff.resize(7);
      suff.push_back("_0_0_0");

      if (AsegPriorType == 2) {		// Local neighborhood aseg prior
        suff.push_back("_1_0_0");
        suff.push_back("_-1_0_0");
        suff.push_back("_0_1_0");
        suff.push_back("_0_-1_0");
        suff.push_back("_0_0_1");
        suff.push_back("_0_0_-1");
      }
      else if (AsegPriorType > 2) {	// Nearest-neighbor aseg prior
        suff.push_back("_x");
        suff.push_back("_-x");
        suff.push_back("_y");
        suff.push_back("_-y");
        suff.push_back("_z");
        suff.push_back("_-z");
      }
    }

    for (vector<string>::const_iterator isuff = suff.begin();
                                        isuff < suff.end(); isuff++) {
      float val;
      sprintf(fname, "%s%s.nii.gz", AsegIdFile, isuff->c_str());
      ifstream idfile(fname, ios::in);

      if (!idfile) {
        cout << "ERROR: Could not open " << fname << endl;
        exit(1);
      }

      cout << "Loading label IDs for aseg prior from " << fname << endl;
      idlist.clear();
      while (idfile >> val)
        idlist.push_back(val);
      segids.push_back(idlist);

      sprintf(fname, "%s%s.nii.gz", AsegPriorFile0, isuff->c_str());
      cout << "Loading aseg prior from " << fname << endl;
      tmp = MRIread(fname);
      if (!tmp) {
        cout << "ERROR: Could not read " << fname << endl;
        exit(1);
      }
      if (tmp->nframes != (int) idlist.size()) {
        cout << "ERROR: Expected " << idlist.size() << " frames in " << fname
             << ", found " << tmp->nframes << endl;
        exit(1);
      }
      segprior0.push_back(tmp);

      sprintf(fname, "%s%s.nii.gz", AsegPriorFile1, isuff->c_str());
      cout << "Loading aseg prior from " << fname << endl;
      tmp = MRIread(fname);
      if (!tmp) {
        cout << "ERROR: Could not read " << fname << endl;
        exit(1);
      }
      if (tmp->nframes != (int) idlist.size()) {
        cout << "ERROR: Expected " << idlist.size() << " frames in " << fname
             << ", found " << tmp->nframes << endl;
        exit(1);
      }
      segprior1.push_back(tmp);
    }
  }

  // Read segmentation map training set
  if (AsegTrainFile) {
    cout << "Loading seg. map training set from " << AsegTrainFile << endl;
    asegtr = MRIread(AsegTrainFile);
    if (!asegtr) {
      cout << "ERROR: Could not read " << AsegTrainFile << endl;
      exit(1);
    }
  }

  // Read path training set
  if (PathTrainFile) {
    cout << "Loading path training set from " << PathTrainFile << endl;
    pathtr = MRIread(PathTrainFile);
    if (!pathtr) {
      cout << "ERROR: Could not read " << PathTrainFile << endl;
      exit(1);
    }
    if (pathtr->nframes != asegtr->nframes) {
      cout << "ERROR: Numbers of training samples in " << PathTrainFile 
           << " and " << AsegTrainFile << " do not match." << endl;
      exit(1);
    }
  }

  // Read segmentation map
  if (AsegFile) {
    cout << "Loading segmentation map from " << AsegFile << endl;
    aseg = MRIread(AsegFile);
    if (!aseg) {
      cout << "ERROR: Could not read " << AsegFile << endl;
      exit(1);
    }
  }

  // Initialize voxel-wise diffusion model
  Bite::InitializeStatic(GradientFile, BvalueFile,
                         NumTract, (unsigned int) phi[0]->nframes,
                         AsegPriorType, segids,
                         pathtr ? (unsigned int) pathtr->nframes : 0);

  if (Bite::GetNumDir() != (unsigned int) dwi->nframes) {
    cout << "ERROR: Dimensions of " << BvalueFile << " and " << DwiFile
         << " do not match" << endl;
    exit(1);
  }

  mData.clear();
  for (unsigned int iz = 0; iz < mNz; iz++)
    for (unsigned int iy = 0; iy < mNy; iy++)
      for (unsigned int ix = 0; ix < mNx; ix++)
        if (MRIgetVoxVal(mMask, (int) ix, (int) iy, (int) iz, 0)) {
          Bite data = Bite(dwi, phi, theta, f, v0, f0, d0,
                           prior0, prior1,
                           segprior0, segprior1,
                           asegtr, pathtr, aseg,
                           ix, iy, iz);
          mData.push_back(data);
        }

  mDataMask.clear();
  mNumVox = 0;
  for (unsigned int iz = 0; iz < mNz; iz++)
    for (unsigned int iy = 0; iy < mNy; iy++)
      for (unsigned int ix = 0; ix < mNx; ix++)
        if (MRIgetVoxVal(mMask, (int) ix, (int) iy, (int) iz, 0)) {
          mDataMask.push_back(&mData[mNumVox]);
          mNumVox++;
        }
        else
          mDataMask.push_back(0);

  // Read end ROIs
  if (RoiFile1) {
    cout << "Loading end ROI from " << RoiFile1 << endl;
    mRoi1 = MRIread(RoiFile1);
    if (!mRoi1) {
      cout << "ERROR: Could not read " << RoiFile1 << endl;
      exit(1);
    }
  }
  else if (RoiMeshFile1) {}

  if (RoiFile2) {
    cout << "Loading end ROI from " << RoiFile2 << endl;
    mRoi2 = MRIread(RoiFile2);
    if (!mRoi2) {
      cout << "ERROR: Could not read " << RoiFile2 << endl;
      exit(1);
    }
  }
  else if (RoiMeshFile2) {}

  // Resolution of end ROIs
  mRoiVoxelSize.clear();
  mRoiVoxelSize.push_back(mRoi1->xsize);
  mRoiVoxelSize.push_back(mRoi1->ysize);
  mRoiVoxelSize.push_back(mRoi1->zsize);

//vector<int> point;
//point.push_back(87);
//point.push_back(61);
//point.push_back(57);
//cout << "in mask " << IsInMask(point.begin()) << " in roi " << IsInRoi(point.begin(), mRoi2) << endl;
//exit(1);

  MRIfree(&dwi);
  for (unsigned int itract = 0; itract < NumTract; itract++) {
    MRIfree(&phi[itract]);
    MRIfree(&theta[itract]);
    MRIfree(&f[itract]);
    MRIfree(&v0[itract]);
    MRIfree(&f0[itract]);
  }
  MRIfree(&d0);
  if (prior0) {
    MRIfree(&prior0);
    MRIfree(&prior1);
  }
  for (unsigned int k = 0; k < segprior0.size(); k++) {
    MRIfree(&segprior0[k]);
    MRIfree(&segprior1[k]);
  }
  if (asegtr)
    MRIfree(&asegtr);
  if (pathtr)
    MRIfree(&pathtr);
  if (aseg)
    MRIfree(&aseg);
}

Coffin::~Coffin() {
  MRIfree(&mMask);
  MRIfree(&mRoi1);
  MRIfree(&mRoi2);
  MRIfree(&mPathSamples);
}

//
// Read DWI-to-ROI registration matrix
//
void Coffin::ReadDwiToRoi(const char *XfmFile) {
  mDwiToRoi.clear();

  if (XfmFile) {
    float val;
    ifstream infile(XfmFile, ios::in);

    if (!infile) {
      cout << "ERROR: Could not open " << XfmFile << endl;
      exit(1);
    }

    cout << "Loading DWI-to-ROI registration from " << XfmFile << endl;
    while (infile >> val)
      mDwiToRoi.push_back(val);

    if (mDwiToRoi.size() != 16) {
      cout << "ERROR: File " << XfmFile << " must contain a 4x4 matrix" << endl;
      exit(1);
    }
  }
  else {	// Identity by default
    float id[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    mDwiToRoi.insert(mDwiToRoi.begin(), id, id+16);
  }
}

//
// Read initial control points
//
void Coffin::ReadControlPoints(const char *ControlPointFile) {
  float coord;
  ifstream infile(ControlPointFile, ios::in);

  if (!infile) {
    cout << "ERROR: Could not open " << ControlPointFile << endl;
    exit(1);
  }

  cout << "Loading control points from " << ControlPointFile << endl;
  mControlPoints.clear();
  while (infile >> coord)
    mControlPoints.push_back(round(coord));

  if (mControlPoints.size() % 3 != 0) {
    cout << "ERROR: File " << ControlPointFile
         << " must contain triplets of coordinates" << endl;
    exit(1);
  }

  mNumControl = mControlPoints.size()/3;
}

//
// Read initial proposal standard deviations for control point perturbations
//
void Coffin::ReadProposalStds(const char *PropStdFile) {
  mProposalStdInit.clear();

  if (PropStdFile) {
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
  else {	// Default values
    mProposalStdInit.resize(mControlPoints.size());
    fill(mProposalStdInit.begin(), mProposalStdInit.begin()+3, 1.0);	//5.0
    fill(mProposalStdInit.begin()+3, mProposalStdInit.end()-3, 1.0);	//1.0
    fill(mProposalStdInit.end()-3, mProposalStdInit.end(), 1.0);	//5.0
  }
}

//
// Run MCMC
//
void Coffin::RunMCMC() {
  unsigned int iprop, ikeep;

  cout << "Initializing MCMC" << endl;
  InitializeMCMC();

  if (mDebug) {
    char fname[PATH_MAX];
    sprintf(fname, "%s/Finit.nii.gz", mOutDir);
    mSpline.WriteVolume(fname, true);
  }

  cout << "Running MCMC burn-in jumps" << endl;
  iprop = 1;
  for (unsigned int ijump = mNumBurnIn; ijump > 0; ijump--) {
    if (JumpMCMC() || mAcceptF || mAcceptTheta) {	// Accept new path
      UpdatePath();
      UpdateAcceptanceRate();

      if (mDebug) {
        char fname[PATH_MAX];
        sprintf(fname, "%s/Faccept_b%05d.nii.gz", mOutDir, mNumBurnIn-ijump+1);
        mSpline.WriteVolume(fname, true);
      }
    }
    else {						// Reject new path
      UpdateRejectionRate();

      if (mDebug) {
        char fname[PATH_MAX];
        sprintf(fname, "%s/Freject_b%05d.nii.gz", mOutDir, mNumBurnIn-ijump+1);
        mSpline.WriteVolume(fname, true);
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
  iprop = 1;
  ikeep = 1;
  for (unsigned int ijump = mNumSample; ijump > 0; ijump--) {
    if (JumpMCMC() || mAcceptF || mAcceptTheta) {	// Accept new path
      UpdatePath();
      UpdateAcceptanceRate();

      if (mDebug) {
        char fname[PATH_MAX];
        sprintf(fname, "%s/Faccept_%05d.nii.gz", mOutDir, mNumSample-ijump+1);
        mSpline.WriteVolume(fname, true);
      }
    }
    else {						// Reject new path
      UpdateRejectionRate();

      if (mDebug) {
        char fname[PATH_MAX];
        sprintf(fname, "%s/Freject_%05d.nii.gz", mOutDir, mNumSample-ijump+1);
        mSpline.WriteVolume(fname, true);
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
}

//
// Run MCMC (single point updates)
//
void Coffin::RunMCMC1() {
  unsigned int iprop, ikeep;
  vector<unsigned int> cptorder(mNumControl);
  vector<unsigned int>::const_iterator icpt;

  cout << "Initializing MCMC" << endl;
  InitializeMCMC();

  if (mDebug) {
    char fname[PATH_MAX];
    sprintf(fname, "%s/Finit.nii.gz", mOutDir);
    mSpline.WriteVolume(fname, true);
  }

  cout << "Running MCMC burn-in jumps" << endl;
  iprop = 1;
  for (unsigned int ijump = mNumBurnIn; ijump > 0; ijump--) {
    // Perturb control points in random order
    for (unsigned int k = 0; k < mNumControl; k++)
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

      if (JumpMCMC1(*icpt) || mAcceptF || mAcceptTheta) {	// Accept point
        UpdatePath();

        if (mDebug) {
          char fname[PATH_MAX];
          sprintf(fname, "%s/Faccept_b%05d_%d.nii.gz",
                  mOutDir, mNumBurnIn-ijump+1, *icpt);
          mSpline.WriteVolume(fname, true);
        }
      }
      else {							// Reject point
        mRejectControl[*icpt] = true;

        if (mDebug) {
          char fname[PATH_MAX];
          sprintf(fname, "%s/Freject_b%05d_%d.nii.gz",
                  mOutDir, mNumBurnIn-ijump+1, *icpt);
          mSpline.WriteVolume(fname, true);
        }
      }
    }

    UpdateAcceptRejectRate1();

    if (iprop == mUpdatePropNth) {
      UpdateProposalStd();
      iprop = 1;
    }
    else
      iprop++;
  }

  mPosteriorOnPathMap = numeric_limits<double>::max();

  cout << "Running MCMC main jumps" << endl;
  iprop = 1;
  ikeep = 1;
  for (unsigned int ijump = mNumSample; ijump > 0; ijump--) {
    // Perturb control points in random order
    for (unsigned int k = 0; k < mNumControl; k++)
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

      if (JumpMCMC1(*icpt) || mAcceptF || mAcceptTheta) {	// Accept point
        UpdatePath();

        if (mDebug) {
          char fname[PATH_MAX];
          sprintf(fname, "%s/Faccept_%05d_%d.nii.gz",
                  mOutDir, mNumSample-ijump+1, *icpt);
          mSpline.WriteVolume(fname, true);
        }
      }
      else {							// Reject point
        mRejectControl[*icpt] = true;

        if (mDebug) {
          char fname[PATH_MAX];
          sprintf(fname, "%s/Freject_%05d_%d.nii.gz",
                  mOutDir, mNumSample-ijump+1, *icpt);
          mSpline.WriteVolume(fname, true);
        }
      }
    }

    UpdateAcceptRejectRate1();

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
}

//
// Initialize path and MCMC proposals
//
void Coffin::InitializeMCMC() {
  vector<float>::iterator phi, theta;

  // Interpolate spline from initial control points
  if (!mSpline.InterpolateSpline()) {
    cout << "Path from initial control points is not entirely in mask" << endl;
    exit (1);
  }

  mSpline.ComputeTangent();

  // Get initial path points
  mPathPoints.clear();
  mPathPoints.insert(mPathPoints.begin(),
                     mSpline.GetAllPointsBegin(), mSpline.GetAllPointsEnd());

  // Get initial path tangent angles
  mPathPhi.resize(mPathPoints.size());
  phi = mPathPhi.begin();;
  mPathTheta.resize(mPathPoints.size());
  theta = mPathTheta.begin();
  for (vector<float>::const_iterator ipt = mSpline.GetTangentBegin();
                                     ipt < mSpline.GetTangentEnd(); ipt += 3) {
    *phi = atan2(ipt[1], ipt[0]);
    *theta = acos(ipt[2] / sqrt(ipt[0]*ipt[0] + ipt[1]*ipt[1] + ipt[2]*ipt[2]));
    if (mDebug)
      cout << (*phi)/M_PI*180 << " " << (*theta)/M_PI*180 << endl;
    phi++;
    theta++;
  }

  // Compute energy on initial path
  mLikelihoodOnPath = 0;
  mPriorOnPath = 0;
  phi = mPathPhi.begin();
  theta = mPathTheta.begin();
  for (vector<int>::const_iterator ipt = mPathPoints.begin();
                                   ipt < mPathPoints.end(); ipt += 3) {
    Bite *ivox = mDataMask[ipt[0] + ipt[1]*mNx + ipt[2]*mNxy];

    ivox->ComputeLikelihoodOnPath(*phi, *theta);
    ivox->ComputePriorOnPath();

    mLikelihoodOnPath += ivox->GetLikelihoodOnPath();
    mPriorOnPath += ivox->GetPriorOnPath();

    phi++;
    theta++;
  }

  mPosteriorOnPath = mLikelihoodOnPath + mPriorOnPath;

  // Initialize control point proposal distribution
  mProposalStd.resize(mProposalStdInit.size());
  copy(mProposalStdInit.begin(), mProposalStdInit.end(), mProposalStd.begin());

  mAcceptCount.resize(mControlPoints.size());
  fill(mAcceptCount.begin(), mAcceptCount.end(), 0.0);
  mRejectCount.resize(mControlPoints.size());
  fill(mRejectCount.begin(), mRejectCount.end(), 0.0);

  mRejectControl.resize(mNumControl);

  mControlPointsNew.resize(mControlPoints.size());
  mControlPointJumps.resize(mControlPoints.size());
  mControlPointsMap.resize(mControlPoints.size());

  // Clear saved path samples
  mControlPointSamples.clear();
  mPathLengthSamples.clear();
  mLikelihoodOnPathSamples.clear();
  mPriorOnPathSamples.clear();
  mPosteriorOnPathSamples.clear();
  MRIclear(mPathSamples);
}

//
// Perform a single MCMC jump
//
bool Coffin::JumpMCMC() {
  fill(mRejectControl.begin(), mRejectControl.end(), false);
  mRejectSpline = false;
  mRejectF = false;
  mAcceptF = false;
  mRejectTheta = false;
  mAcceptTheta = false;
  mRejectPosterior = false;

  if (!ProposePath())
    return false;

  ProposeDiffusionParameters();

  if (!AcceptPath())
    return false;

  return true;
}

//
// Perform a single MCMC jump (single point updates)
//
bool Coffin::JumpMCMC1(unsigned int ControlIndex) {
  if (!ProposePath1(ControlIndex))
    return false;

  ProposeDiffusionParameters();

  if (!AcceptPath())
    return false;

  return true;
}

//
// Propose path by perturbing spline control points
//
bool Coffin::ProposePath() {
  vector<bool>::iterator isrej;
  vector<int>::const_iterator coord = mControlPoints.begin();
  vector<int>::const_iterator cpoint = mControlPointsNew.begin();
  vector<int>::iterator newcoord = mControlPointsNew.begin();
  vector<float>::const_iterator pstd = mProposalStd.begin();
  vector<float>::iterator jump = mControlPointJumps.begin();
  vector<float>::iterator phi, theta;

  // Perturb current control points
  if (mDebug)
    cout << "Jump norms ";
  for (isrej = mRejectControl.begin(); isrej < mRejectControl.end(); isrej++) {
    double norm = 0;

    for (unsigned int ii = 0; ii < 3; ii++) {
      *jump = round((*pstd) * PDFgaussian());
      *newcoord = *coord + *jump;

      *jump *= *jump;
      norm += *jump;

      jump++;
      pstd++;
      coord++;
      newcoord++;
    }

    if (mDebug)
      cout << sqrt(norm) << " ";

    if (norm > 0) {
      jump -= 3;
      for (unsigned int ii = 0; ii < 3; ii++) {
        *jump /= norm;
        jump++;
      }
    }

    // Check that new control point is in mask
    if (!IsInMask(cpoint))
      *isrej = true;

    cpoint += 3;
  }

  if (mDebug)
    cout << endl;

  // Check that new end points are in end ROIs
  if (!IsInRoi(mControlPointsNew.begin(), mRoi1))
    mRejectControl[0] = true;

  if (!IsInRoi(mControlPointsNew.end() - 3, mRoi2))
    mRejectControl[mNumControl-1] = true;

  for (isrej = mRejectControl.begin(); isrej < mRejectControl.end(); isrej++)
    if (*isrej)
      return false;

  // Interpolate spline from proposed control points
  mSpline.SetControlPoints(mControlPointsNew);
  if (!mSpline.InterpolateSpline()) {
    mRejectSpline = true;
    return false;
  }
    
  mSpline.ComputeTangent();

  // Get proposed path points
  mPathPointsNew.clear();
  mPathPointsNew.insert(mPathPointsNew.begin(),
                        mSpline.GetAllPointsBegin(), mSpline.GetAllPointsEnd());

  // Get proposed path tangent angles
  mPathPhiNew.resize(mPathPointsNew.size());
  phi = mPathPhiNew.begin();;
  mPathThetaNew.resize(mPathPointsNew.size());
  theta = mPathThetaNew.begin();
  for (vector<float>::const_iterator ipt = mSpline.GetTangentBegin();
                                     ipt < mSpline.GetTangentEnd(); ipt += 3) {
    *phi = atan2(ipt[1], ipt[0]);
    *theta = acos(ipt[2] / sqrt(ipt[0]*ipt[0] + ipt[1]*ipt[1] + ipt[2]*ipt[2]));
    phi++;
    theta++;
  }

  return true;
}

//
// Propose path by perturbing spline control points (single point updates)
//
bool Coffin::ProposePath1(unsigned int ControlIndex) {
  const unsigned int offset = ControlIndex * 3;
  double norm = 0;
  vector<int>::const_iterator coord = mControlPoints.begin() + offset;
  vector<int>::const_iterator cpoint = mControlPointsNew.begin() + offset;
  vector<int>::iterator newcoord = mControlPointsNew.begin() + offset;
  vector<float>::const_iterator pstd = mProposalStd.begin() + offset;
  vector<float>::iterator jump = mControlPointJumps.begin() + offset;
  vector<float>::iterator phi, theta;

  copy(mControlPoints.begin(), mControlPoints.end(), mControlPointsNew.begin());

  // Perturb current control point
  for (unsigned int ii = 0; ii < 3; ii++) {
    *jump = round((*pstd) * PDFgaussian());
    *newcoord = *coord + *jump;

    *jump *= *jump;
    norm += *jump;

    jump++;
    pstd++;
    coord++;
    newcoord++;
  }

  if (norm > 0) {
    jump -= 3;
    for (unsigned int ii = 0; ii < 3; ii++) {
      *jump /= norm;
      jump++;
    }
  }

  // Check that new control point is in mask and, if end point, in end ROI
  if (!IsInMask(cpoint) ||
      ((ControlIndex == 0) && !IsInRoi(cpoint, mRoi1)) ||
      ((ControlIndex == mNumControl-1) && !IsInRoi(cpoint, mRoi2))) {
    mRejectControl[ControlIndex] = true;
    return false;
  }

  // Interpolate spline from proposed control points
  mSpline.SetControlPoints(mControlPointsNew);
  if (!mSpline.InterpolateSpline()) {
    mRejectSpline = true;
    return false;
  }
    
  mSpline.ComputeTangent();

  // Get proposed path points
  mPathPointsNew.clear();
  mPathPointsNew.insert(mPathPointsNew.begin(),
                        mSpline.GetAllPointsBegin(), mSpline.GetAllPointsEnd());

  // Get proposed path tangent angles
  mPathPhiNew.resize(mPathPointsNew.size());
  phi = mPathPhiNew.begin();;
  mPathThetaNew.resize(mPathPointsNew.size());
  theta = mPathThetaNew.begin();
  for (vector<float>::const_iterator ipt = mSpline.GetTangentBegin();
                                     ipt < mSpline.GetTangentEnd(); ipt += 3) {
    *phi = atan2(ipt[1], ipt[0]);
    *theta = acos(ipt[2] / sqrt(ipt[0]*ipt[0] + ipt[1]*ipt[1] + ipt[2]*ipt[2]));
    phi++;
    theta++;
  }

  return true;
}

//
// Propose diffusion parameters by sampling from their marginal posteriors
//
void Coffin::ProposeDiffusionParameters() {
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
// Determine if proposed path will be accepted
//
bool Coffin::AcceptPath() {
  double neglogratio;
  vector<int>::const_iterator ipt;
  vector<float>::const_iterator phi, theta;

  mLikelihoodOnPathNew = 0;
  mPriorOnPathNew = 0;
  mLikelihoodOffPathNew = 0;
  mPriorOffPathNew = 0;

  // Compute posterior on proposed path
  phi = mPathPhiNew.begin();
  theta = mPathThetaNew.begin();
  for (ipt = mPathPointsNew.begin(); ipt < mPathPointsNew.end(); ipt += 3) {
    Bite *ivox = mDataMask[ipt[0] + ipt[1]*mNx + ipt[2]*mNxy];

    ivox->ComputeLikelihoodOffPath();
    ivox->ComputeLikelihoodOnPath(*phi, *theta);
    if (ivox->IsFZero()) {
      mRejectF = true;
      return false;
    }
    if (ivox->IsThetaZero()) {
      mAcceptTheta = true;
      return false;
    }
    ivox->ComputePriorOffPath();
    ivox->ComputePriorOnPath();

    if (mDebug & 0)
      cout << ipt[0] << " " << ipt[1] << " " << ipt[2] << " "
           << ivox->GetLikelihoodOffPath() << " "
           << ivox->GetPriorOffPath() << " "
           << ivox->GetPathPriorOffPath() << " "
           << ivox->GetLikelihoodOnPath() << " "
           << ivox->GetPriorOnPath() << " "
           << ivox->GetPathPriorOnPath() << endl;

    mLikelihoodOnPathNew += ivox->GetLikelihoodOnPath();
    mPriorOnPathNew += ivox->GetPriorOnPath();

    mLikelihoodOffPathNew += ivox->GetLikelihoodOffPath();
    mPriorOffPathNew += ivox->GetPriorOffPath();

    phi++;
    theta++;
  }

  mPosteriorOnPathNew  = mLikelihoodOnPathNew  + mPriorOnPathNew;
  mPosteriorOffPathNew = mLikelihoodOffPathNew + mPriorOffPathNew;

  if (mDebug)
    cout << "mLikelihoodOnPathNew=" << mLikelihoodOnPathNew << " "
         << "mPriorOnPathNew=" << mPriorOnPathNew << endl
         << "mLikelihoodOffPathNew=" << mLikelihoodOffPathNew << " "
         << "mPriorOffPathNew=" << mPriorOffPathNew << endl
         << "mLikelihoodOn-OffPathNew="
         << mLikelihoodOnPathNew-mLikelihoodOffPathNew << " "
         << "mPriorOn-OffPathNew=" << mPriorOnPathNew-mPriorOffPathNew << endl
         << "mPathLengthNew=" << mPathPointsNew.size()/3 << endl;

  mLikelihoodOnPath = 0;
  mPriorOnPath = 0;
  mLikelihoodOffPath = 0;
  mPriorOffPath = 0;

  // Compute posterior on current path
  phi = mPathPhi.begin();
  theta = mPathTheta.begin();
  for (ipt = mPathPoints.begin(); ipt < mPathPoints.end(); ipt += 3) {
    Bite *ivox = mDataMask[ipt[0] + ipt[1]*mNx + ipt[2]*mNxy];

    ivox->ComputeLikelihoodOffPath();
    ivox->ComputeLikelihoodOnPath(*phi, *theta);
    if (ivox->IsFZero()) {
      mAcceptF = true;
      return false;
    }
    if (ivox->IsThetaZero()) {
      mRejectTheta = true;
      return false;
    }
    ivox->ComputePriorOffPath();
    ivox->ComputePriorOnPath();

    mLikelihoodOnPath += ivox->GetLikelihoodOnPath();
    mPriorOnPath += ivox->GetPriorOnPath();

    mLikelihoodOffPath += ivox->GetLikelihoodOffPath();
    mPriorOffPath += ivox->GetPriorOffPath();

    phi++;
    theta++;
  }

  mPosteriorOnPath  = mLikelihoodOnPath  + mPriorOnPath;
  mPosteriorOffPath = mLikelihoodOffPath + mPriorOffPath;

  if (mDebug)
    cout << "mLikelihoodOnPath=" << mLikelihoodOnPath << " "
         << "mPriorOnPath=" << mPriorOnPath << endl
         << "mLikelihoodOffPath=" << mLikelihoodOffPath << " "
         << "mPriorOffPath=" << mPriorOffPath << endl
         << "mLikelihoodOff-OnPath="
         << mLikelihoodOffPath-mLikelihoodOnPath << " "
         << "mPriorOff-OnPath=" << mPriorOffPath-mPriorOnPath << endl
         << "mPathLength=" << mPathPoints.size()/3 << endl;

  neglogratio = (mPosteriorOnPathNew - mPosteriorOffPathNew)/(mPathPointsNew.size()/3) 
              + (mPosteriorOffPath   - mPosteriorOnPath)/(mPathPoints.size()/3);
//neglogratio = mLikelihoodOnPathNew/(mPathPointsNew.size()/3) - mLikelihoodOnPath/(mPathPoints.size()/3);

  if (mDebug)
    cout << "alpha=" << exp(-neglogratio) << endl;

  // Accept or reject based on ratio of posteriors
  if (drand48() < exp(-neglogratio))
    return true;

  mRejectPosterior = true;
  return false;
}

//
// Copy newly accepted path over current path
//
void Coffin::UpdatePath() {
  copy(mControlPointsNew.begin(), mControlPointsNew.end(), 
       mControlPoints.begin());

  mPathPoints.resize(mPathPointsNew.size());
  copy(mPathPointsNew.begin(), mPathPointsNew.end(), mPathPoints.begin());

  mPathPhi.resize(mPathPhiNew.size());
  copy(mPathPhiNew.begin(), mPathPhiNew.end(), mPathPhi.begin());

  mPathTheta.resize(mPathThetaNew.size());
  copy(mPathThetaNew.begin(), mPathThetaNew.end(), mPathTheta.begin());

  mLikelihoodOnPath  = mLikelihoodOnPathNew; 
  mLikelihoodOffPath = mLikelihoodOffPathNew; 
  mPriorOnPath       = mPriorOnPathNew; 
  mPriorOffPath      = mPriorOffPathNew; 
  mPosteriorOnPath   = mPosteriorOnPathNew;
  mPosteriorOffPath  = mPosteriorOffPathNew;
}

//
// Update control point acceptance rate
//
void Coffin::UpdateAcceptanceRate() {
  vector<float>::const_iterator jump = mControlPointJumps.begin();
  
  if (!(mAcceptF || mAcceptTheta))
    // Increment acceptance count for all control points
    for (vector<float>::iterator count = mAcceptCount.begin();
                                 count < mAcceptCount.end(); count++) {
      *count += *jump;
      jump++;
    }
}

//
// Update control point rejection rate
//
void Coffin::UpdateRejectionRate() {
  vector<float>::const_iterator jump = mControlPointJumps.begin();

  if (mDebug)
    if (mRejectPosterior)
      cout << "Reject due to posterior" << endl;
    else if (mRejectSpline)
      cout << "Reject due to spline" << endl;
    else
      cout << "Reject due to control point" << endl;

  if (mRejectPosterior || mRejectSpline)
    // Increment rejection count for all control points
    for (vector<float>::iterator count = mRejectCount.begin();
                                 count < mRejectCount.end(); count++) {
      *count += *jump;
      jump++;
    }
  else {
    vector<float>::iterator count = mRejectCount.begin();

    // Increment rejection count for specific rejected control points
    for (vector<bool>::const_iterator isrej = mRejectControl.begin();
                                      isrej < mRejectControl.end(); isrej++)
      if (*isrej)
        for (unsigned int ii = 0; ii < 3; ii++) {
          *count += *jump;
          jump++;
          count++;
        }
      else {
        jump += 3;
        count += 3;
      }
  }
}

//
// Update control point acceptance and rejection rates (single point update)
//
void Coffin::UpdateAcceptRejectRate1() {
  vector<float>::const_iterator jump = mControlPointJumps.begin();
  vector<float>::iterator acount = mAcceptCount.begin();
  vector<float>::iterator rcount = mRejectCount.begin();

  for (vector<bool>::const_iterator isrej = mRejectControl.begin();
                                    isrej < mRejectControl.end(); isrej++)
    if (*isrej)		// Increment rejection count
      for (unsigned int ii = 0; ii < 3; ii++) {
//        *rcount += *jump;
(*rcount)++;
        jump++;
        acount++;
        rcount++;
      }
    else		// Increment acceptance count
      for (unsigned int ii = 0; ii < 3; ii++) {
//        *acount += *jump;
(*acount)++;
        jump++;
        acount++;
        rcount++;
      }

acount = mAcceptCount.begin();
rcount = mRejectCount.begin();
cout << "Acceptance rates ";
for (vector<bool>::const_iterator isrej = mRejectControl.begin(); isrej < mRejectControl.end(); isrej++) {
cout << *acount << "/" << *acount+*rcount << "=" << *acount/(*acount+*rcount) << " ";
acount+=3; rcount+=3;
}
cout << endl;

  if (mDebug) {
    cout << "Normalized jumps ";

    for (jump = mControlPointJumps.begin(); jump < mControlPointJumps.end();
                                            jump++)
      cout << sqrt(*jump) << " ";

    cout << endl;
  }
}

//
// Update control point proposal distribution
//
void Coffin::UpdateProposalStd() {
  vector<float>::const_iterator acount = mAcceptCount.begin();
  vector<float>::const_iterator rcount = mRejectCount.begin();

  for (vector<float>::iterator pstd = mProposalStd.begin();
                               pstd < mProposalStd.end(); pstd++) {
    *pstd *= sqrt((*acount + 1) / (*rcount + 1));

    acount++;
    rcount++;
  }

  if (mDebug) {
    acount = mAcceptCount.begin();
    rcount = mRejectCount.begin();

    cout << "Acceptance rates ";
    for (vector<float>::const_iterator pstd = mProposalStd.begin();
                                       pstd < mProposalStd.end(); pstd++) {
      cout << *acount << "/" << *acount+*rcount << " ";

      acount++;
      rcount++;
    }
    cout << endl;

    cout << "Updated proposal STDs ";
    for (vector<float>::const_iterator pstd = mProposalStd.begin();
                                       pstd < mProposalStd.end(); pstd++)
      cout << *pstd << " ";
    cout << endl;
  }

  fill(mAcceptCount.begin(), mAcceptCount.end(), 0.0);
  fill(mRejectCount.begin(), mRejectCount.end(), 0.0);
}

//
// Save current path as an MCMC sample
//
void Coffin::SavePath() {
  // Save current control points
  mControlPointSamples.insert(mControlPointSamples.end(),
                              mControlPoints.begin(), mControlPoints.end());

  // Save current path length
  mPathLengthSamples.push_back(mPathPoints.size()/3);

  // Save current path energy
  mLikelihoodOnPathSamples.push_back(mLikelihoodOnPath);
  mPriorOnPathSamples.push_back(mPriorOnPath);
  mPosteriorOnPathSamples.push_back(mPosteriorOnPath);

  // Add current path to path posterior
  for (vector<int>::const_iterator ipt = mPathPoints.begin();
                                   ipt < mPathPoints.end(); ipt += 3) {
    const int ix = (int) ipt[0],
              iy = (int) ipt[1],
              iz = (int) ipt[2];
      
    MRIsetVoxVal(mPathSamples, ix, iy, iz, 0, 
                 MRIgetVoxVal(mPathSamples, ix, iy, iz, 0) + 1);
  }

  // Keep track of MAP path
  if (mPosteriorOnPath < mPosteriorOnPathMap) {
    copy(mControlPoints.begin(), mControlPoints.end(), 
         mControlPointsMap.begin());
    mPathPointsMap.resize(mPathPoints.size());
    copy(mPathPoints.begin(), mPathPoints.end(), mPathPointsMap.begin());
    mPosteriorOnPathMap = mPosteriorOnPath;
  }
}

//
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

  ApplyAffine(outpoint, Point, mDwiToRoi, mRoiVoxelSize, mDwiVoxelSize);

  return (outpoint[0] > -1) && (outpoint[0] < Roi->width) &&
         (outpoint[1] > -1) && (outpoint[1] < Roi->height) &&
         (outpoint[2] > -1) && (outpoint[2] < Roi->depth) &&
         (MRIgetVoxVal(Roi, outpoint[0], outpoint[1], outpoint[2], 0) > 0);
}

//
// Apply an affine transform to a single point
//
void Coffin::ApplyAffine(vector<int> &OutPoint,
                         vector<int>::const_iterator InPoint,
                         vector<float> &In2Out,
                         vector<float> &OutVoxelSize,
                         vector<float> &InVoxelSize) {
  vector<float>::const_iterator in2out = In2Out.begin();
  vector<float> pin, pout;

  pin.resize(4);
  pout.resize(4);

  for (unsigned int i = 0; i < 3; i++)
    pin[i] = InPoint[i] * InVoxelSize[i];
  pin[3] = 1;

  for (unsigned int i = 0; i < 4; i++) {
    float psum = 0;
    for (unsigned int j = 0; j < 4; j++) {
      psum += (*in2out) * pin[j];
      in2out++;
    }
    pout[i] = psum;
  }

  for (unsigned int i = 0; i < 3; i++)
    OutPoint[i] = round(pout[i] / pout[3] / OutVoxelSize[i]);
}

//
// Write output files
//
void Coffin::WriteOutputs() {
  char fname[PATH_MAX];
  MRI *out1 = MRIclone(mMask, NULL);
  MRI *out2 = MRIclone(mMask, NULL);

  cout << "Writing output files to " << mOutDir << endl;

  // Save volume of path samples
  sprintf(fname, "%s/Fsamples_1_1.nii.gz", mOutDir);
  MRIwrite(mPathSamples, fname);

  // Save volumes of end point samples
  for (vector<int>::const_iterator icpt = mControlPointSamples.begin();
                                 icpt < mControlPointSamples.end(); icpt += 3) {
    int ix, iy, iz;

    ix = (int) icpt[0];
    iy = (int) icpt[1];
    iz = (int) icpt[2];
    MRIsetVoxVal(out1, ix, iy, iz, 0, MRIgetVoxVal(out1, ix, iy, iz, 0) + 1);

    icpt += 3 * (mNumControl-1);

    ix = (int) icpt[0];
    iy = (int) icpt[1];
    iz = (int) icpt[2];
    MRIsetVoxVal(out2, ix, iy, iz, 0, MRIgetVoxVal(out2, ix, iy, iz, 0) + 1);
  }

  sprintf(fname, "%s/Lsamples_1_1.nii.gz", mOutDir);
  MRIwrite(out1, fname);

  sprintf(fname, "%s/Tsamples_1_1.nii.gz", mOutDir);
  MRIwrite(out2, fname);

  // Save control point samples
  sprintf(fname, "%s/SAMPLES_1_1.txt", mOutDir);
  ofstream sampfile(fname, ios::out);

  if (!sampfile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  for (vector<int>::const_iterator icpt = mControlPointSamples.begin();
                                   icpt < mControlPointSamples.end(); icpt += 3)
    sampfile << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;

  // Save length of path samples
  sprintf(fname, "%s/LENGTH_1_1.txt", mOutDir);
  ofstream lenfile(fname, ios::out);

  if (!lenfile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  for (vector<int>::const_iterator ilen = mPathLengthSamples.begin();
                                   ilen < mPathLengthSamples.end(); ilen++)
    lenfile << *ilen << endl;

  // Save likelihood and prior of path samples
  sprintf(fname, "%s/FOLLOW_LIKELIHOOD_1_1.txt", mOutDir);
  ofstream likefile(fname, ios::out);

  if (!likefile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  vector<float>::const_iterator ilik = mLikelihoodOnPathSamples.begin();
  for (vector<float>::const_iterator ipri = mPriorOnPathSamples.begin();
                                     ipri < mPriorOnPathSamples.end(); ipri++) {
    likefile << *ilik << " " << *ipri << endl;
    ilik++;
  }

  // Save MAP control point sample
  sprintf(fname, "%s/CONTROLS_1_1.txt", mOutDir);
  ofstream mapfile(fname, ios::out);

  if (!mapfile) {
    cout << "ERROR: Could not open " << fname << " for writing" << endl;
    exit(1);
  }

  for (vector<int>::const_iterator icpt = mControlPointsMap.begin();
                                   icpt < mControlPointsMap.end(); icpt += 3)
    mapfile << icpt[0] << " " << icpt[1] << " " << icpt[2] << endl;

  // Save MAP path sample
  MRIclear(out1);
  for (vector<int>::const_iterator ipt = mPathPointsMap.begin();
                                   ipt < mPathPointsMap.end(); ipt += 3) {
    const int ix = (int) ipt[0],
              iy = (int) ipt[1],
              iz = (int) ipt[2];
    MRIsetVoxVal(out1, ix, iy, iz, 0, MRIgetVoxVal(out1, ix, iy, iz, 0) + 1);
  }

  sprintf(fname, "%s/CONTROLS_1_1_spline.nii.gz", mOutDir);
  MRIwrite(out1, fname);

  MRIfree(&out1);
  MRIfree(&out2);
}

