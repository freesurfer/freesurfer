/**
 * @brief Random-forrest classifier for white-matter segmentation
 *
 * Random-forrest classifier for white-matter segmentation
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

#ifndef FORREST_H
#define FORREST_H

#include "mri.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <limits>
#include <math.h>
#include <sstream>
#include <vector>

class Forrest {
public:
  Forrest();
  ~Forrest();
  void ReadTestSubject(const char *TestDir, const char *MaskFile,
                       const char *AsegFile, const char *OrientFile);
  void ReadTrainingSubjects(const char *TrainListFile, const char *MaskFile,
                            const char *AsegFile, const char *OrientFile,
                            std::vector<char *> TractFileList);
  int  GetNx();
  int  GetNy();
  int  GetNz();
  int  GetNumTrain();
  bool IsInMask(int CoordX, int CoordY, int CoordZ);
  std::vector<unsigned int> GetTestAseg(int CoordX, int CoordY, int CoordZ);
  std::vector<float>        GetTestOrient(int CoordX, int CoordY, int CoordZ);
  std::vector<int>          GetTrainXyz(int SampleIndex);
  std::vector<unsigned int> GetTrainAseg(int SampleIndex);
  std::vector<float>        GetTrainOrient(int SampleIndex);
  std::vector<unsigned int> GetTrainTractIds(int SampleIndex);

private:
  static const int mSampleStep;

  int                       mNx, mNy, mNz, mNumTrain, mNumLocal, mNumNear;
  std::vector<int>          mDirLocal, mDirNear, mTrainXyz;
  std::vector<unsigned int> mTrainAsegIdsLocal, mTrainAsegIdsNear;
  std::vector<float>        mTrainAsegDist, mTrainOrient;
  std::vector<std::vector<unsigned int>> mTrainTractIds;
  MRI *                                  mMask, *mAseg, *mOrient;
};

#endif
