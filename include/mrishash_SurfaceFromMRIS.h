#pragma once

/*
 * Original Author: Bevin R. Brett
 * CVS Revision Info:
 *
 * Copyright © 2019 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "mrishash.h"
#include "mrisurf_SurfaceFromMRIS_generated.h"

MRIS_HASH_TABLE* MHTcreateFaceTable             (SurfaceFromMRIS::Analysis::Surface surface);
MRIS_HASH_TABLE* MHTcreateFaceTable_Resolution  (SurfaceFromMRIS::Analysis::Surface surface, int which, float res);
MRIS_HASH_TABLE* MHTcreateVertexTable           (SurfaceFromMRIS::Analysis::Surface surface, int which);
MRIS_HASH_TABLE* MHTcreateVertexTable_Resolution(SurfaceFromMRIS::Analysis::Surface surface, int which, float res);

