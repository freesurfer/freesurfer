#include "mri.h"
#include "mrisurf.h"

#ifndef TESS_INCLUDED
#define TESS_INCLUDED

MRIS* MRIScreateSurfaceFromVolume(MRI *mri,int label,int connectivity);
MRIS** MRIScreateSurfacesFromVolume(MRI *mri,int number_of_labels, int* labelvalue,int connectivity);

#endif
