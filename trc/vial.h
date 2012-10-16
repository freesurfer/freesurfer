/**
 * @file  vial.h
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

#ifndef VIAL_H
#define VIAL_H

//define NO_CVS_UP_IN_HERE

#ifndef NO_CVS_UP_IN_HERE
// Needed for CVS - these must be included first or else they don't compile
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include <mpi.h>
#include "../fem_elastic/morph.h"
#include "../fem_elastic/surf_utils.h"
#include "../fem_elastic/morph_utils.h"
#include "gcamorph.h"
	/*
#include "../fem_elastic/morph.h"
#include "../fem_elastic/morph_utils.h"
	*/
#endif

#include <vector>
#include <iostream>
#include <fstream>
#include "mri.h"

class AffineReg {
  public:
    AffineReg();
    ~AffineReg();
    bool IsEmpty();
    void ReadXfm(const char *XfmFile, const MRI *InRefVol, 
                                      const MRI *OutRefVol);
    void ApplyXfm(std::vector<float> &OutPoint,
                  std::vector<float>::const_iterator InPoint);

  private:
    std::vector<float> mInToOut,			// [4 x 4]
                       mInVoxelSize, mOutVoxelSize;	// [3]
};

#ifndef NO_CVS_UP_IN_HERE
class NonlinReg {
  public:
    NonlinReg();
    ~NonlinReg();
    bool IsEmpty();
    void ReadXfm(const char *XfmFile, MRI *OutRefVol);
    void ApplyXfm(std::vector<float> &OutPoint,
                  std::vector<float>::const_iterator InPoint);
    void ApplyXfmInv(std::vector<float> &OutPoint,
                     std::vector<float>::const_iterator InPoint);

  private:
//    boost::shared_ptr<gmp::VolumeMorph> mMorph;
    GCAM *mMorph;
};
#endif

#endif

