#ifndef HIPPO_H
#define HIPPO_H

#include "mri.h"

MRI *HIPPOremoveNonHippoLabels(MRI *mri_src, MRI *mri_dst) ;
MRI *HIPPOestimateIntensityImage(MRI *mri_hippo_labels, MRI *mri_aseg, MRI *mri_intensity, MRI *mri_dst);

#endif
