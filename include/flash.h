#ifndef FLASH_H
#define FLASH_H

#include "mri.h"

double dFlash_dT1(double T1, double PD, double TR, double fa, double TE) ;
double dFlash_dPD(double T1, double PD, double TR, double fa, double TE) ;
double FLASHforwardModel(double T1, double PD, double TR, double flip_angle, double TE) ;
double FLASHforwardModelT2star(double T1, double PD, double T2star, double TR, double flip_angle, double TE) ;
MRI    *MRIparameterMapsToFlash(MRI *mri_src, MRI *mri_dst, double *TRs, double *TEs, double *FAs, int nflash) ;
int    compute_T1_PD(int nvals, float *vals, double *TRs, double *FAs, double *TEs, double *pT1, double *pPD);
int    FlashBuildLookupTables(int nvolumes, double *TRs, double *FAs, double *TEs) ;

#endif
