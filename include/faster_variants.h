/**
 * @brief define the macros that conditionalize faster variants of some algorithms
 *        before those variants completely replace the original variant
 *
 *        keep both variants when the faster variant is 
 *           1) not yet fully tested, or
 *           2) substantially more complex than the original
 *                  because comparing the two results can test the faster one for correctness
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

#define FASTER_gcamSmoothnessEnergy
#define FASTER_MRI_EM_REGISTER
#define FASTER_mriSurfaceRASToVoxel
