
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
MRI *fMRIcomputeF(MRI *ces, MATRIX *X, MATRIX *C, MRI *var, MRI *F);
MRI *fMRIsigT(MRI *t, float DOF, MRI *p);
MRI *fMRIsigF(MRI *F, float DOF1, float DOF2, MRI *sig);
MRI *fMRIsumSquare(MRI *fmri, int Update, MRI *sumsqr);
MRI *fMRInskip(MRI *inmri, int nskip, MRI *outmri);

#endif
