
#ifndef FMRIUTILS_INC
#define FMRIUTILS_INC

#include "matrix.h"
#include "mri.h"

#ifdef X
#undef X
#endif

MRI *fMRImatrixMultiply(MRI *inmri, MATRIX *M, MRI *outmri);
MRI *fMRIvariance(MRI *fmri, float DOF, int RmMean, MRI *var);
MRI *fMRIcomputeT(MRI *ces, MATRIX *X, MATRIX *C, MRI *var, MRI *t);
MRI *fMRIsigT(MRI *t, float DOF, MRI *p);
MRI *fMRIsumSquare(MRI *fmri, int Update, MRI *sumsqr);

#endif
