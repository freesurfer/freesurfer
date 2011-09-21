/**
 * @file  vial.cxx
 * @brief Holds utilities for probabilistic tractography
 *
 * Holds utilities for probabilistic tractography
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

#include <vial.h>

using namespace std;

//
// Affine registration class
//
AffineReg::AffineReg() {}

AffineReg::~AffineReg() {}

bool AffineReg::IsEmpty() { return mInToOut.empty(); }

//
// Read affine registration matrix from file and set input/output resolutions
//
void AffineReg::ReadXfm(const char *XfmFile, const MRI *InRefVol,
                                             const MRI *OutRefVol) {
  // Read registration matrix from file
  mInToOut.clear();

  if (XfmFile) {
    float val;
    ifstream infile(XfmFile, ios::in);

    if (!infile) {
      cout << "ERROR: Could not open " << XfmFile << endl;
      exit(1);
    }

    cout << "Loading affine registration from " << XfmFile << endl;
    while (infile >> val)
      mInToOut.push_back(val);

    if (mInToOut.size() != 16) {
      cout << "ERROR: File " << XfmFile << " must contain a 4x4 matrix" << endl;
      exit(1);
    }
  }
  else {	// Identity by default
    float id[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    mInToOut.insert(mInToOut.begin(), id, id+16);
  }

  // Get resolution of input images
  mInVoxelSize.resize(3);
  if (InRefVol) {
    mInVoxelSize[0] = InRefVol->xsize;
    mInVoxelSize[1] = InRefVol->ysize;
    mInVoxelSize[2] = InRefVol->zsize;
  }
  else
    fill(mInVoxelSize.begin(), mInVoxelSize.end(), 1.0);

  // Get resolution of output images
  mOutVoxelSize.resize(3);
  if (OutRefVol) {
    mOutVoxelSize[0] = OutRefVol->xsize;
    mOutVoxelSize[1] = OutRefVol->ysize;
    mOutVoxelSize[2] = OutRefVol->zsize;
  }
  else
    fill(mOutVoxelSize.begin(), mOutVoxelSize.end(), 1.0);
}

//
// Apply an affine transform to a single point
//
void AffineReg::ApplyXfm(vector<float> &OutPoint,
                         vector<float>::const_iterator InPoint) {
  vector<float>::const_iterator in2out = mInToOut.begin();
  vector<float> pin, pout;

  pin.resize(4);
  pout.resize(4);

  for (int i = 0; i < 3; i++)
    pin[i] = InPoint[i] * mInVoxelSize[i];
  pin[3] = 1;

  for (int i = 0; i < 4; i++) {
    float psum = 0;
    for (int j = 0; j < 4; j++) {
      psum += (*in2out) * pin[j];
      in2out++;
    }
    pout[i] = psum;
  }

  for (int i = 0; i < 3; i++)
    OutPoint[i] = pout[i] / pout[3] / mOutVoxelSize[i];
}

#ifndef NO_CVS_UP_IN_HERE
	/*
//
// Non-linear registration class
//
NonlinReg::NonlinReg() {
  mMorph = boost::shared_ptr<gmp::VolumeMorph>(new gmp::VolumeMorph);
}

NonlinReg::~NonlinReg() {}

bool NonlinReg::IsEmpty() { return (mMorph->m_template == 0); }

//
// Read a non-linear transform from file
//
void NonlinReg::ReadXfm(const char *XfmFile, MRI *OutRefVol) {
  unsigned int zlibBuffer = 5;

  ifstream xfile(XfmFile, ios::in);	// Just to check if file exists
  if (!xfile) {
    cout << "ERROR: Could not open " << XfmFile << endl;
    exit(1);
  }
  xfile.close();
  
  cout << "Loading non-linear registration from " << XfmFile << endl;
  mMorph->m_template = OutRefVol;

  try {
    mMorph->load(XfmFile, zlibBuffer);
  }
  catch (const char* msg) {
    cout << "Exception caught while loading registration: " << msg << endl;
    exit(1);
  }

  mMorph->m_interpolationType = SAMPLE_NEAREST;
  mMorph->invert();
}

//
// Apply a non-linear transform to a single point
//
void NonlinReg::ApplyXfm(vector<float> &OutPoint,
                         vector<float>::const_iterator InPoint) {
  Coords3d inpt, outpt;

  for (int k = 0; k < 3; k++)
    inpt(k) = InPoint[k];

  outpt = mMorph->image(inpt);

  for (int k = 0; k < 3; k++)
    OutPoint[k] = (float) outpt(k);
}
	*/

//
// Non-linear registration class
//
NonlinReg::NonlinReg() : mMorph(0) {}

NonlinReg::~NonlinReg() {}

bool NonlinReg::IsEmpty() { return (mMorph == 0); }

//
// Read a non-linear transform from file
//
void NonlinReg::ReadXfm(const char *XfmFile, MRI *OutRefVol) {
  ifstream xfile(XfmFile, ios::in);	// Just to check if file exists
  if (!xfile) {
    cout << "ERROR: Could not open " << XfmFile << endl;
    exit(1);
  }
  xfile.close();

  cout << "Loading non-linear registration from " << XfmFile << endl;
  mMorph = GCAMreadAndInvertNonTal(XfmFile);

  if (mMorph == NULL) exit(1);

  mMorph->gca = gcaAllocMax(1, 1, 1, OutRefVol->width, OutRefVol->height,
                                                       OutRefVol->depth, 0, 0);
}

//
// Apply a non-linear transform to a single point
//
void NonlinReg::ApplyXfm(vector<float> &OutPoint,
                         vector<float>::const_iterator InPoint) {
  float inpoint[3];

  copy(InPoint, InPoint+3, inpoint);

  GCAMmorphPlistFromAtlas(1, inpoint, mMorph, &OutPoint[0]);
}

//
// Apply the inverse of a non-linear transform to a single point
//
void NonlinReg::ApplyXfmInv(vector<float> &OutPoint,
                            vector<float>::const_iterator InPoint) {
  float inpoint[3];

  copy(InPoint, InPoint+3, inpoint);

  GCAMmorphPlistToSource(1, inpoint, mMorph, &OutPoint[0]);
}
#endif

