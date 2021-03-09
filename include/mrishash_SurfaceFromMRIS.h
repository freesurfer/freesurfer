#pragma once

/*
 * Original Author: Bevin R. Brett
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

#include "mrishash.h"
#include "mrisurf_SurfaceFromMRIS_generated.h"

namespace Minimal_Surface_MRIS = SurfaceFromMRIS::XYZPosition;

MRIS_HASH_TABLE* MHTcreateFaceTable             (Minimal_Surface_MRIS::Surface surface);
MRIS_HASH_TABLE* MHTcreateFaceTable_Resolution  (Minimal_Surface_MRIS::Surface surface, int which, float res);
MRIS_HASH_TABLE* MHTcreateVertexTable           (Minimal_Surface_MRIS::Surface surface, int which);
MRIS_HASH_TABLE* MHTcreateVertexTable_Resolution(Minimal_Surface_MRIS::Surface surface, int which, float res);

