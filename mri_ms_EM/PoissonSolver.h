/* CVS Version Info
   $Id: PoissonSolver.h,v 1.1 2006/02/01 22:00:22 xhan Exp $
*/

#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include "mri.h"

void multigrid(MRI *u, MRI *f, MRI *hu, int XN, int YN, int ZN);

#endif
