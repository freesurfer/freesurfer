/**
 * @file  gcam_mi.h
 * @brief morphed GCA
 *
 */
/*
 * Original Author: Gheorghe Postelnicu, MGH, June 2007
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2010/02/12 19:33:27 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2010,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
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
