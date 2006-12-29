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
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
