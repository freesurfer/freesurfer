/**
 * @file  gcam_mi.h
 * @brief morphed GCA
 *
 */
/*
 * Original Author: Gheorghe Postelnicu, MGH, June 2007
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:23 $
 *    $Revision: 1.2 $
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

#ifndef _h_gcam_mi_grad_h
#define _h_gcam_mi_grad_h


#ifdef __cplusplus
extern "C"
{
#endif

#include "mri.h"
#include "gcamorph.h"

  MRI* GCAMcomputeMIDenseGradient(GCAM* gcam,
                                  MRI* mriMoving);
#ifdef __cplusplus
};
#endif

#endif
