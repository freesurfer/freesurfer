
#ifndef unwarpGradientNonlinearity_H
#define unwarpGradientNonlinearity_H

#include"mri.h"

MRI *unwarpGradientNonlinearity(MRI *mri, 
				char *unwarp_gradientType, 
				char *unwarp_partialUnwarp, 
				char *unwarp_jacobianCorrection,
				char *unwarp_interpType,
				int unwarp_sincInterpHW);



#ifdef unwarpGradientNonlinearity_SRC
int unwarp_flag = 0;
char unwarp_gradientType[STRLEN];
char unwarp_partialUnwarp[STRLEN];
char unwarp_jacobianCorrection[STRLEN];
char unwarp_interpType[STRLEN];
int unwarp_sincInterpHW = 0;

#else
extern int unwarp_flag;
extern char unwarp_gradientType[STRLEN];
extern char unwarp_partialUnwarp[STRLEN];
extern char unwarp_jacobianCorrection[STRLEN];
extern char unwarp_interpType[STRLEN];
extern int unwarp_sincInterpHW;

#endif /* #ifdef unwarpGradientNonlinearity_SRC */

#endif /* #ifndef unwarpGradientNonlinearity_H */
