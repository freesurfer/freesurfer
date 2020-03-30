/*
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


/* CVS Version Info
   $Id: PoissonSolver.h,v 1.3 2011/03/02 00:04:23 nicks Exp $
*/

#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include "mri.h"

void multigrid(MRI *u, MRI *f, MRI *hu, int XN, int YN, int ZN);

#endif
