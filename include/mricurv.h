#include "mri.h"
#include "mri_topology.h"

#ifndef MRICURV_INCLUDED
#define MRICURV_INCLUDED

float MRIcurvature(MRI *mri,int i,int j,int k,int label,int connectivity);
float Nbhcurvature(Nbh *nbh,int connectivity);


#endif





