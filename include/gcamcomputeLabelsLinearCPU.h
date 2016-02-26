#pragma once

#include "gca.h"
#include "gcamorph.h"
#include "mri.h"


#if defined(__cplusplus)
extern "C" {
#endif


int GCAMcomputeLabelsLinearCPU(MRI *mri, GCA_MORPH *gcam);


#if defined(__cplusplus)
};
#endif
