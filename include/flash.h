#ifndef FLASH_H
#define FLASH_H

#include "mri.h"

double dFlash_dT1(double T1, double PD, double TR, double fa, double TE) ;
double dFlash_dPD(double T1, double PD, double TR, double fa, double TE) ;
double FLASHforwardModel(double T1, double PD, double TR, double flip_angle, double TE) ;
MRI    *MRIparameterMapsToFlash(MRI *mri_src, MRI *mri_dst, double *TRs, double *TEs, double *FAs, int nflash) ;

#endif
