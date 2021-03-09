/**
 * @brief functions for finding the best rotation for aligning a MRIS with a target
 *
 */
/*
 * Original Author: Bevin Brett
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

#pragma once

#include "mrisurf.h"

void MRISrigidBodyAlignGlobal_findMinSSE(
    double* mina, double* new_minb, double* new_ming, double* new_sse,  // outputs
    MRI_SURFACE*        mris,
    INTEGRATION_PARMS*  parms,
    float               min_radians,
    float               max_radians,
    double              ext_sse,
    int                 nangles);
