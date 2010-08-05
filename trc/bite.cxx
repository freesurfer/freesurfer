/**
 * @file  bite.cxx
 * @brief Tractography data and methods for an individual voxel
 *
 * Tractography data and methods for an individual voxel
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

#include <bite.h>

using namespace std;

unsigned int Bite::mNumDir, Bite::mNumB0, Bite::mNumTract, Bite::mNumBedpost,
             Bite::mAsegPriorType, Bite::mNumTrain;
const float Bite::mFminPath = 0.02;
vector<float> Bite::mGradients, Bite::mBvalues;
vector< vector<unsigned int> > Bite::mAsegIds;

Bite::Bite(MRI *Dwi, MRI **Phi, MRI **Theta, MRI **F,
           MRI **V0, MRI **F0, MRI *D0,
           MRI *Prior0, MRI *Prior1,
           vector<MRI *> &AsegPrior0, vector<MRI *> &AsegPrior1,
           MRI *AsegTrain, MRI *PathTrain, MRI *Aseg,
           unsigned int CoordX, unsigned int CoordY, unsigned int CoordZ) :
           mCoordX(CoordX), mCoordY(CoordY), mCoordZ(CoordZ) {
  float fsum, vx, vy, vz;

  mDwi.clear();
  mPhiSamples.clear();
  mThetaSamples.clear();
  mFSamples.clear();
  mPhi.clear();
  mTheta.clear();
  mF.clear();

  // Initialize s0
  mS0 = 0;
  for (unsigned int idir = 0; idir < mNumDir; idir++) {
    mDwi.push_back(MRIgetVoxVal(Dwi, mCoordX, mCoordY, mCoordZ, idir));
    if (mBvalues[idir] == 0)
      mS0 += MRIgetVoxVal(Dwi, mCoordX, mCoordY, mCoordZ, idir);
  }
  mS0 /= mNumB0;

  // Samples of phi, theta, f
  for (unsigned int isamp = 0; isamp < mNumBedpost; isamp++)
    for (unsigned int itract = 0; itract < mNumTract; itract++) {
      mPhiSamples.push_back(MRIgetVoxVal(Phi[itract],
                                         mCoordX, mCoordY, mCoordZ, isamp));
      mThetaSamples.push_back(MRIgetVoxVal(Theta[itract],
                                         mCoordX, mCoordY, mCoordZ, isamp));
      mFSamples.push_back(MRIgetVoxVal(F[itract],
                                         mCoordX, mCoordY, mCoordZ, isamp));
    }

  fsum = 0;
  for (unsigned int itract = 0; itract < mNumTract; itract++) {
    // Initialize phi, theta
    vx = MRIgetVoxVal(V0[itract], mCoordX, mCoordY, mCoordZ, 0),
    vy = MRIgetVoxVal(V0[itract], mCoordX, mCoordY, mCoordZ, 1),
    vz = MRIgetVoxVal(V0[itract], mCoordX, mCoordY, mCoordZ, 2);
    mPhi.push_back(atan2(vy, vx));
    mTheta.push_back(acos(vz / sqrt(vx*vx + vy*vy + vz*vz)));

    // Initialize f
    mF.push_back(MRIgetVoxVal(F0[itract], mCoordX, mCoordY, mCoordZ, 0));
    fsum += MRIgetVoxVal(F0[itract], mCoordX, mCoordY, mCoordZ, 0);
  }

  // Initialize d
  mD = MRIgetVoxVal(D0, mCoordX, mCoordY, mCoordZ, 0);
  //mD = log(mDwi[mNumDir-1] / mS0 / (1-fsum);

  // Spatial path priors
  if (Prior0 && Prior1) {
    mPathPrior0 = MRIgetVoxVal(Prior0, mCoordX, mCoordY, mCoordZ, 0);
    mPathPrior1 = MRIgetVoxVal(Prior1, mCoordX, mCoordY, mCoordZ, 0);
  }
  else {
    mPathPrior0 = 0;
    mPathPrior1 = 0;
  }

  // Aseg priors
  if (!AsegPrior0.empty() && !AsegPrior1.empty()) {
    vector<float> tmp0, tmp1;
    vector<unsigned int>::const_iterator iid;
    vector< vector<unsigned int> >::const_iterator iids;
    vector<MRI *>::const_iterator iprior0 = AsegPrior0.begin(),
                                  iprior1 = AsegPrior1.begin();

    for (iids = mAsegIds.begin(); iids != mAsegIds.end(); iids++) {
      unsigned int k = 0;

      tmp0.clear();
      tmp1.clear();

      for (iid = iids->begin(); iid != iids->end(); iid++) {
        tmp0.push_back(MRIgetVoxVal(*iprior0, mCoordX, mCoordY, mCoordZ, k));
        tmp1.push_back(MRIgetVoxVal(*iprior1, mCoordX, mCoordY, mCoordZ, k));
        k++;
      }

      mAsegPrior0.push_back(tmp0);
      mAsegPrior1.push_back(tmp1);

      iprior0++;
      iprior1++;
    }
  }

  // Training paths and segmentation maps
  if (AsegTrain && PathTrain) {
    for (unsigned int itrain = 0; itrain < mNumTrain; itrain++) {
      mPathTrain.push_back((unsigned int)
                 MRIgetVoxVal(PathTrain, mCoordX, mCoordY, mCoordZ, itrain));

      mAsegTrain.push_back((unsigned int)
                 MRIgetVoxVal(AsegTrain, mCoordX, mCoordY, mCoordZ, itrain));
    }

    if (mAsegPriorType == 2) {
      unsigned int coord;

      coord = (mCoordX < (unsigned int) Aseg->width-1 ? mCoordX+1 : mCoordX);
      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++)
        mAsegTrain.push_back((unsigned int)
                   MRIgetVoxVal(AsegTrain, coord, mCoordY, mCoordZ, itrain));

      coord = (mCoordX > 0 ? mCoordX-1 : mCoordX);
      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++)
        mAsegTrain.push_back((unsigned int)
                   MRIgetVoxVal(AsegTrain, coord, mCoordY, mCoordZ, itrain));

      coord = (mCoordY < (unsigned int) Aseg->height-1 ? mCoordY+1 : mCoordY);
      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++)
        mAsegTrain.push_back((unsigned int)
                   MRIgetVoxVal(AsegTrain, mCoordX, coord, mCoordZ, itrain));

      coord = (mCoordY > 0 ? mCoordY-1 : mCoordY);
      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++)
        mAsegTrain.push_back((unsigned int)
                   MRIgetVoxVal(AsegTrain, mCoordX, coord, mCoordZ, itrain));

      coord = (mCoordZ < (unsigned int) Aseg->depth-1 ? mCoordZ+1 : mCoordZ);
      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++)
        mAsegTrain.push_back((unsigned int)
                   MRIgetVoxVal(AsegTrain, mCoordX, mCoordY, coord, itrain));

      coord = (mCoordZ > 0 ? mCoordZ-1 : mCoordZ);
      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++)
        mAsegTrain.push_back((unsigned int)
                   MRIgetVoxVal(AsegTrain, mCoordX, mCoordY, coord, itrain));
    }
    else if (mAsegPriorType > 2) {
      unsigned int coord;
      float seg, seg0;

      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++) {
        seg0 = MRIgetVoxVal(Aseg, mCoordX, mCoordY, mCoordZ, itrain);
        seg = seg0;
        coord = mCoordX;
        while ((coord < (unsigned int) Aseg->width-1) && (seg == seg0)) {
          coord++;
          seg = MRIgetVoxVal(Aseg, coord, mCoordY, mCoordZ, itrain);
        }
        mAsegTrain.push_back((unsigned int) seg);

        if (mAsegPriorType == 4)
          mAsegDistTrain.push_back(coord - mCoordX);
      }

      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++) {
        seg0 = MRIgetVoxVal(Aseg, mCoordX, mCoordY, mCoordZ, itrain);
        seg = seg0;
        coord = mCoordX;
        while ((coord > 0) && (seg == seg0)) {
          coord--;
          seg = MRIgetVoxVal(Aseg, coord, mCoordY, mCoordZ, itrain);
        }
        mAsegTrain.push_back((unsigned int) seg);

        if (mAsegPriorType == 4)
          mAsegDistTrain.push_back(mCoordX - coord);
      }

      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++) {
        seg0 = MRIgetVoxVal(Aseg, mCoordX, mCoordY, mCoordZ, itrain);
        seg = seg0;
        coord = mCoordY;
        while ((coord < (unsigned int) Aseg->height-1) && (seg == seg0)) {
          coord++;
          seg = MRIgetVoxVal(Aseg, mCoordX, coord, mCoordZ, itrain);
        }
        mAsegTrain.push_back((unsigned int) seg);

        if (mAsegPriorType == 4)
          mAsegDistTrain.push_back(coord - mCoordY);
      }

      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++) {
        seg0 = MRIgetVoxVal(Aseg, mCoordX, mCoordY, mCoordZ, itrain);
        seg = seg0;
        coord = mCoordY;
        while ((coord > 0) && (seg == seg0)) {
          seg = MRIgetVoxVal(Aseg, mCoordX, coord, mCoordZ, itrain);
          coord--;
        }
        mAsegTrain.push_back((unsigned int) seg);

        if (mAsegPriorType == 4)
          mAsegDistTrain.push_back(mCoordY - coord);
      }

      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++) {
        seg0 = MRIgetVoxVal(Aseg, mCoordX, mCoordY, mCoordZ, itrain);
        seg = seg0;
        coord = mCoordZ;
        while ((coord < (unsigned int) Aseg->depth-1) && (seg == seg0)) {
          seg = MRIgetVoxVal(Aseg, mCoordX, mCoordY, coord, itrain);
          coord++;
        }
        mAsegTrain.push_back((unsigned int) seg);

        if (mAsegPriorType == 4)
          mAsegDistTrain.push_back(coord - mCoordZ);
      }

      for (unsigned int itrain = 0; itrain < mNumTrain; itrain++) {
        seg0 = MRIgetVoxVal(Aseg, mCoordX, mCoordY, mCoordZ, itrain);
        seg = seg0;
        coord = mCoordZ;
        while ((coord > 0) && (seg == seg0)) {
          seg = MRIgetVoxVal(Aseg, mCoordX, mCoordY, coord, itrain);
            coord--;
        }
        mAsegTrain.push_back((unsigned int) seg);

        if (mAsegPriorType == 4)
          mAsegDistTrain.push_back(mCoordZ - coord);
      }
    }
   }

  // Segmentation map
  if (Aseg) {
    mAseg.push_back((unsigned int)
          MRIgetVoxVal(Aseg, mCoordX, mCoordY, mCoordZ, 0));

    if (mAsegPriorType == 2) {
      unsigned int coord;

      coord = (mCoordX < (unsigned int) Aseg->width-1 ? mCoordX+1 : mCoordX);
      mAseg.push_back((unsigned int)
            MRIgetVoxVal(Aseg, coord, mCoordY, mCoordZ, 0));

      coord = (mCoordX > 0 ? mCoordX-1 : mCoordX);
      mAseg.push_back((unsigned int)
            MRIgetVoxVal(Aseg, coord, mCoordY, mCoordZ, 0));

      coord = (mCoordY < (unsigned int) Aseg->height-1 ? mCoordY+1 : mCoordY);
      mAseg.push_back((unsigned int)
            MRIgetVoxVal(Aseg, mCoordX, coord, mCoordZ, 0));

      coord = (mCoordY > 0 ? mCoordY-1 : mCoordY);
      mAseg.push_back((unsigned int)
            MRIgetVoxVal(Aseg, mCoordX, coord, mCoordZ, 0));

      coord = (mCoordZ < (unsigned int) Aseg->depth-1 ? mCoordZ+1 : mCoordZ);
      mAseg.push_back((unsigned int)
            MRIgetVoxVal(Aseg, mCoordX, mCoordY, coord, 0));

      coord = (mCoordZ > 0 ? mCoordZ-1 : mCoordZ);
      mAseg.push_back((unsigned int)
            MRIgetVoxVal(Aseg, mCoordX, mCoordY, coord, 0));
    }
    else if (mAsegPriorType > 2) {
      unsigned int coord;
      float seg;
      const float seg0 = MRIgetVoxVal(Aseg, mCoordX, mCoordY, mCoordZ, 0);

      seg = seg0;
      coord = mCoordX;
      while ((coord < (unsigned int) Aseg->width-1) && (seg == seg0)) {
        coord++;
        seg = MRIgetVoxVal(Aseg, coord, mCoordY, mCoordZ, 0);
      }
      mAseg.push_back((unsigned int) seg);

      if (mAsegPriorType == 4)
        mAsegDist.push_back(coord - mCoordX);

      seg = seg0;
      coord = mCoordX;
      while ((coord > 0) && (seg == seg0)) {
        coord--;
        seg = MRIgetVoxVal(Aseg, coord, mCoordY, mCoordZ, 0);
      }
      mAseg.push_back((unsigned int) seg);

      if (mAsegPriorType == 4)
        mAsegDist.push_back(mCoordX - coord);

      seg = seg0;
      coord = mCoordY;
      while ((coord < (unsigned int) Aseg->height-1) && (seg == seg0)) {
        coord++;
        seg = MRIgetVoxVal(Aseg, mCoordX, coord, mCoordZ, 0);
      }
      mAseg.push_back((unsigned int) seg);

      if (mAsegPriorType == 4)
        mAsegDist.push_back(coord - mCoordY);

      seg = seg0;
      coord = mCoordY;
      while ((coord > 0) && (seg == seg0)) {
        seg = MRIgetVoxVal(Aseg, mCoordX, coord, mCoordZ, 0);
        coord--;
      }
      mAseg.push_back((unsigned int) seg);

      if (mAsegPriorType == 4)
        mAsegDist.push_back(mCoordY - coord);

      seg = seg0;
      coord = mCoordZ;
      while ((coord < (unsigned int) Aseg->depth-1) && (seg == seg0)) {
        seg = MRIgetVoxVal(Aseg, mCoordX, mCoordY, coord, 0);
        coord++;
      }
      mAseg.push_back((unsigned int) seg);

      if (mAsegPriorType == 4)
        mAsegDist.push_back(coord - mCoordZ);

      seg = seg0;
      coord = mCoordZ;
      while ((coord > 0) && (seg == seg0)) {
        seg = MRIgetVoxVal(Aseg, mCoordX, mCoordY, coord, 0);
          coord--;
      }
      mAseg.push_back((unsigned int) seg);

      if (mAsegPriorType == 4)
        mAsegDist.push_back(mCoordZ - coord);
    }
  }
}

Bite::~Bite() {
}

//
// Initialize variables that are common for all voxels
//
void Bite::InitializeStatic(const char *GradientFile, const char *BvalueFile,
                            unsigned int NumTract, unsigned int NumBedpost,
                            unsigned int AsegPriorType,
                            const vector< vector<unsigned int> > &AsegIds,
                            unsigned int NumTrain) {
  float val;
  ifstream gfile(GradientFile, ios::in);
  ifstream bfile(BvalueFile, ios::in);

  if (!gfile) {
    cout << "ERROR: Could not open " << GradientFile << endl;
    exit(1);
  }
  if (!bfile) {
    cout << "ERROR: Could not open " << BvalueFile << endl;
    exit(1);
  }

  cout << "Loading b-values from " << BvalueFile << endl;
  mBvalues.clear();
  mNumB0 = 0;
  while (bfile >> val) {
    mBvalues.push_back(val);
    if (val == 0) mNumB0++;
  }
  mNumDir = mBvalues.size();

  cout << "Loading gradients from " << GradientFile << endl;
  mGradients.clear();
  mGradients.resize(3*mNumDir);
  for (unsigned int ii = 0; ii < 3; ii++)
    for (unsigned int idir = 0; idir < mNumDir; idir++) {
      if (!(gfile >> val)) {
        cout << "ERROR: Dimensions of " << BvalueFile << " and "
             << GradientFile << " do not match" << endl;
        exit(1);
      }
      mGradients.at(ii + 3*idir) = val;
    }

  if (gfile >> val) {
    cout << "ERROR: Dimensions of " << BvalueFile << " and " << GradientFile
         << " do not match" << endl;
    exit(1);
  }

  mNumTract = NumTract;

  mNumBedpost = NumBedpost;

  mAsegPriorType = AsegPriorType;

  mAsegIds.clear();
  for (vector< vector<unsigned int> >::const_iterator ipr = AsegIds.begin();
                                                     ipr < AsegIds.end(); ipr++)
    mAsegIds.push_back(*ipr);

  mNumTrain = NumTrain;
}

unsigned int Bite::GetNumTract() { return mNumTract; }

unsigned int Bite::GetNumDir() { return mNumDir; }

unsigned int Bite::GetNumB0() { return mNumB0; }

unsigned int Bite::GetNumBedpost() { return mNumBedpost; }

//
// Draw samples from marginal posteriors of diffusion parameters
//
void Bite::SampleParameters() {
  const unsigned int isamp = (unsigned int) round(drand48() * (mNumBedpost-1))
                                            * mNumTract;
  vector<float>::const_iterator samples;
 
  samples = mPhiSamples.begin() + isamp;
  copy(samples, samples + mNumTract, mPhi.begin());
 
  samples = mThetaSamples.begin() + isamp;
  copy(samples, samples + mNumTract, mTheta.begin());
 
  samples = mFSamples.begin() + isamp;
  copy(samples, samples + mNumTract, mF.begin());
}

//
// Compute likelihood given that voxel is off path
//
void Bite::ComputeLikelihoodOffPath() {
  double like = 0;
  vector<float>::const_iterator ri = mGradients.begin();
  vector<float>::const_iterator bi = mBvalues.begin();
  vector<float>::const_iterator sij = mDwi.begin();

  for (unsigned int idir = mNumDir; idir > 0; idir--) {
    double sbar = 0, fsum = 0;
    const double bidj = (*bi) * mD;
    vector<float>::const_iterator fjl = mF.begin();
    vector<float>::const_iterator phijl = mPhi.begin();
    vector<float>::const_iterator thetajl = mTheta.begin();

    for (unsigned int itract = mNumTract; itract > 0; itract--) {
      const double iprod =
        (ri[0] * cos(*phijl) + ri[1] * sin(*phijl)) * sin(*thetajl)
        + ri[2] * cos(*thetajl);

      sbar += *fjl * exp(-bidj * iprod * iprod);
      fsum += *fjl;

      fjl++;
      phijl++;
      thetajl++;
    }

    sbar += (1-fsum) * exp(-bidj);
    sbar *= mS0;
    like += pow(*sij - sbar, 2);

    ri+=3;
    bi++;
    sij++;
  }

  mLikelihood0 = (float) log(like/2) * mNumDir/2;
}

//
// Compute likelihood given that voxel is on path
//
void Bite::ComputeLikelihoodOnPath(float PathPhi, float PathTheta) {
  double like = 0;
  vector<float>::const_iterator ri = mGradients.begin();
  vector<float>::const_iterator bi = mBvalues.begin();
  vector<float>::const_iterator sij = mDwi.begin();

  // Choose which anisotropic compartment in voxel corresponds to path
  ChoosePathTractAngle(PathPhi, PathTheta);

  // Calculate likelihood by replacing the chosen tract orientation from path
  for (unsigned int idir = mNumDir; idir > 0; idir--) {
    double sbar = 0, fsum = 0;
    const double bidj = (*bi) * mD;
    vector<float>::const_iterator fjl = mF.begin();
    vector<float>::const_iterator phijl = mPhi.begin();
    vector<float>::const_iterator thetajl = mTheta.begin();

    for (unsigned int itract = 0; itract < mNumTract; itract++) {
      double iprod;
      if (itract == mPathTract)
        iprod = (ri[0] * cos(PathPhi) + ri[1] * sin(PathPhi)) * sin(PathTheta)
              + ri[2] * cos(PathTheta);
      else
        iprod = (ri[0] * cos(*phijl) + ri[1] * sin(*phijl)) * sin(*thetajl)
              + ri[2] * cos(*thetajl);

      sbar += *fjl * exp(-bidj * iprod * iprod);
      fsum += *fjl;

      fjl++;
      phijl++;
      thetajl++;
    }

    sbar += (1-fsum) * exp(-bidj);
    sbar *= mS0;
    like += pow(*sij - sbar, 2);

    ri+=3;
    bi++;
    sij++;
  }

  mLikelihood1 = (float) log(like/2) * mNumDir/2;
}

//
// Find tract closest to path orientation
//
void Bite::ChoosePathTractAngle(float PathPhi, float PathTheta) {
  double maxprod = 0;
  vector<float>::const_iterator fjl = mF.begin();
  vector<float>::const_iterator phijl = mPhi.begin();
  vector<float>::const_iterator thetajl = mTheta.begin();

  for (unsigned int itract = 0; itract < mNumTract; itract++) {
    if (*fjl > mFminPath) {
      const double iprod =
        (cos(PathPhi) * cos(*phijl) + sin(PathPhi) * sin(*phijl)) 
        * sin(PathTheta) * sin(*thetajl) + cos(PathTheta) * cos(*thetajl);

      if (iprod > maxprod) {
        mPathTract = itract;
        maxprod = iprod;
      }
    }

    fjl++;
    phijl++;
    thetajl++;
  }

  if (maxprod == 0)
    mPathTract = 0;
}

//
// Find tract that changes the likelihood the least
//
void Bite::ChoosePathTractLike(float PathPhi, float PathTheta) {
  double mindlike = numeric_limits<double>::max();

  for (unsigned int jtract = 0; jtract < mNumTract; jtract++)
    if (mF[jtract] > mFminPath) {
      double dlike, like = 0;
      vector<float>::const_iterator ri = mGradients.begin();
      vector<float>::const_iterator bi = mBvalues.begin();
      vector<float>::const_iterator sij = mDwi.begin();

      // Calculate likelihood by replacing the chosen tract orientation from path
      for (unsigned int idir = mNumDir; idir > 0; idir--) {
        double sbar = 0, fsum = 0;
        const double bidj = (*bi) * mD;
        vector<float>::const_iterator fjl = mF.begin();
        vector<float>::const_iterator phijl = mPhi.begin();
        vector<float>::const_iterator thetajl = mTheta.begin();

        for (unsigned int itract = 0; itract < mNumTract; itract++) {
          double iprod;
          if (itract == jtract)
            iprod = (ri[0] * cos(PathPhi) + ri[1] * sin(PathPhi)) * sin(PathTheta)
                  + ri[2] * cos(PathTheta);
          else
            iprod = (ri[0] * cos(*phijl) + ri[1] * sin(*phijl)) * sin(*thetajl)
                  + ri[2] * cos(*thetajl);

          sbar += *fjl * exp(-bidj * iprod * iprod);
          fsum += *fjl;

          fjl++;
          phijl++;
          thetajl++;
        }

        sbar += (1-fsum) * exp(-bidj);
        sbar *= mS0;
        like += pow(*sij - sbar, 2);

        ri+=3;
        bi++;
        sij++;
      }

      like = log(like/2) * mNumDir/2;
      dlike = fabs(like - (double) mLikelihood0);

      if (dlike < mindlike) {
        mPathTract = jtract;
//        mLikelihood1 = (float) like;
        mindlike = dlike;
      }
    }

  if (mindlike == numeric_limits<double>::max()) {
    mPathTract = 0;
//    mLikelihood1 = 
  }
}

//
// Compute prior given that voxel is off path
//
void Bite::ComputePriorOffPath() {
  vector<float>::const_iterator fjl = mF.begin() + mPathTract;
  vector<float>::const_iterator thetajl = mTheta.begin() + mPathTract;
  vector<unsigned int>::const_iterator iaseg = mAseg.begin();


//cout << (*fjl) << " " << log((*fjl - 1) * log(1 - *fjl)) << " "
//     << log(((double)*fjl - 1) * log(1 - (double)*fjl)) << endl;

if (1) 
  mPrior0 = log((*fjl - 1) * log(1 - *fjl)) - log(fabs(sin(*thetajl)))
          + mPathPrior0;
else  
  mPrior0 = mPathPrior0;

  if (!mAsegPrior0.empty()) {
    vector<unsigned int>::const_iterator iid;
    vector< vector<unsigned int> >::const_iterator iids;
    vector<float>::const_iterator iprior;
    vector< vector<float> >::const_iterator ipriors;

    ipriors = mAsegPrior0.begin();

    for (iids = mAsegIds.begin(); iids != mAsegIds.end(); iids++) {
      bool isinlist = false;

      iprior = ipriors->begin();

      for (iid = iids->begin(); iid != iids->end(); iid++) {
        if (*iaseg == *iid) {
          mPrior0 += *iprior;
          isinlist = true;
          break;
        }

        iprior++;
      }

      if (!isinlist)
        mPrior0 += 0.6931; // Instead should be: mPrior0 -= log((nA+1)/(nA+2));

      iaseg++;
      ipriors++;
    }
  }

  if (!mAsegTrain.empty()) {
    vector<unsigned int>::const_iterator iasegtr = mAsegTrain.begin();

    while (iaseg != mAseg.end()) {
      unsigned int n1 = 0, n2 = 0;
      vector<unsigned int>::const_iterator ipathtr = mPathTrain.begin();

      for (unsigned int itrain = mNumTrain; itrain > 0; itrain--) {
        if (*iaseg == *iasegtr) {
          n2++;
          if (*ipathtr > 0)
            n1++;
        }

        iasegtr++;
        ipathtr++;
      }

      mPrior0 -= log(float(n2-n1+1) / (n2+2));

      iaseg++;
    }
  }
}

//
// Compute prior given that voxel is on path
//
void Bite::ComputePriorOnPath() {
  vector<unsigned int>::const_iterator iaseg = mAseg.begin();

  mPrior1 = mPathPrior1;

  if (!mAsegPrior1.empty()) {
    vector<unsigned int>::const_iterator iid;
    vector< vector<unsigned int> >::const_iterator iids;
    vector<float>::const_iterator iprior;
    vector< vector<float> >::const_iterator ipriors;

    ipriors = mAsegPrior1.begin();

    for (iids = mAsegIds.begin(); iids != mAsegIds.end(); iids++) {
      bool isinlist = false;

      iprior = ipriors->begin();

      for (iid = iids->begin(); iid != iids->end(); iid++) {
        if (*iaseg == *iid) {
          mPrior1 += *iprior;
          isinlist = true;
          break;
        }

        iprior++;
      }

      if (!isinlist)
        mPrior1 += 0.6931; // Instead should be: mPrior1 -= log(1/(nA+2));

      iaseg++;
      ipriors++;
    }
  }

  if (!mAsegTrain.empty()) {
    vector<unsigned int>::const_iterator iasegtr = mAsegTrain.begin();

    while (iaseg != mAseg.end()) {
      unsigned int n1 = 0, n2 = 0;
      vector<unsigned int>::const_iterator ipathtr = mPathTrain.begin();

      for (unsigned int itrain = mNumTrain; itrain > 0; itrain--) {
        if (*iaseg == *iasegtr) {
          n2++;
          if (*ipathtr > 0)
            n1++;
        }

        iasegtr++;
        ipathtr++;
      }

      mPrior1 -= log(float(n1+1) / (n2+2));

      iaseg++;
    }
  }
}

bool Bite::IsFZero() { return (mF[mPathTract] < mFminPath); }

bool Bite::IsThetaZero() { return (mTheta[mPathTract] == 0); }

float Bite::GetLikelihoodOffPath() { return mLikelihood0; }

float Bite::GetLikelihoodOnPath() { return mLikelihood1; }

float Bite::GetPathPriorOffPath() { return mPathPrior0; }

float Bite::GetPathPriorOnPath() { return mPathPrior1; }

float Bite::GetPriorOffPath() { return mPrior0; }

float Bite::GetPriorOnPath() { return mPrior1; }

float Bite::GetPosteriorOffPath() { return mLikelihood0 + mPrior0; }

float Bite::GetPosteriorOnPath() { return mLikelihood1 + mPrior1; }

