// Do not include, intended only to be included into mrisurf_deform.c

#ifdef FS_CUDA

#define MRIS_COMPUTESSE_CUDA
#include "mrisurf_deform_computeSSE.h"
#undef  MRIS_COMPUTESSE_CUDA

#endif /* FS_CUDA */

