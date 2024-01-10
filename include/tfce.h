/*
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

/**
 * @brief routines for performing Threshold Free Cluster Enhancement (TFCE)
 *
 * Based on Smith and Nichols, Neuroimage, 2009
 # Example usage
    mri = MRIread(argv[1]);
    tfce.mask = MRIread(argv[2]);
    tfce.surf = MRISread(argv[3]);
    tfce.hmin = FLT_EPSILON;
    tfce.hmax = 4;
    tfce.nh = 50;
    tfce.E = 0.5;
    tfce.H = 2.0;
    tfce.hlistUniform();
    tfcemap = tfce.compute(mri);
    MRIwrite(tfce,argv[4]);
  If surf is NULL, then it assumes a volume topology
  Note: volume topo has not been extensively tested (1/10/24)
 */

#ifndef TFCE_H
#define TFCE_H

#include <vector>
#include "mrisurf.h"

class TFCE {
public:
  MRIS *surf=NULL; //If using surf topo; if null assumes volume topo
  MRI *mask=NULL; // include voxel if mask>0.5
  std::vector<double> hlist; // threshold list
  double hmin=-1, hdelta=0.1, hmax=-1; // limits for threshold
  int nh=50; // number of threshold samples
  int thsign = 0; // 0=abs, +1=pos, -1=neg
  double E=0.5,H=2; // TFCE parameters; these vals are suggested by Smith and Nichols
  int debug=0, vnodebug=0;
  int hlistUniform(void);
  int hlistAuto(MRI *map);
  MRI *compute(MRI *map); // compute the TFCE map
  // Below are useful functions unrelated to TFCE
  std::vector<double> maxstatsim(MRI *temp, int niters);
  int write_vector_double(char *fname,  std::vector<double> vlist);
  std::vector<double> read_vector_double(char *fname);
  MRI *voxcor(MRI *uncormap, std::vector<double> maxstatlist);// voxelwise multiple comp cor
  double getmax(MRI *statmap, MRI *mask);
};

#endif
