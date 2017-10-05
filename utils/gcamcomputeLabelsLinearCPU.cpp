/**
 * @file  gcamcomputeLabelsLinearCPU.cpp
 * @brief Implement GCAMcomputeLabels using the linearised GCA on the CPU
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/26 20:19:38 $
 *    $Revision: 1.1 $
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

#include "error.h"
#include "gca.h"
#include "gcalinearnode.hpp"
#include "gcalinearprior.hpp"
#include "gcamorph.h"
#include "mri.h"

#include "gcamcomputeLabelsLinearCPU.h"

int GCAMcomputeLabelsLinearCPU(MRI *mri, GCA_MORPH *gcam) {
  int nchanged = 0;

  if (gcam->gca == NULL) {
    return (NO_ERROR);
  }

  Freesurfer::GCAlinearNode gcaLN;
  Freesurfer::GCAlinearPrior gcaLP;

  gcaLN.Exhume(gcam->gca);
  gcaLP.Exhume(gcam->gca);

  return nchanged;
}
