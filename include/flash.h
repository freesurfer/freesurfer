/**
 * @file  flash.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.5 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


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
