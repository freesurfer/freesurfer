/**
 * @file  mrisutils.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.9 $
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


#include "mri.h"
#include "mrisurf.h"

#ifndef SUTILS_INCLUDED
#define SUTILS_INCLUDED

void MRISsmoothSurface(MRI_SURFACE *mris,int niter,float step);
void MRISsmoothSurface2(MRI_SURFACE *mris,int niter,float step,int avrg);
void MRIScenterCOG(MRI_SURFACE *mris);
void MRIScenterCOG2(MRI_SURFACE *mris,double *xCOG,double *yCOG,double *zCOG);
MRI* MRISpeelVolume(MRIS *mris,MRI *mri_src,MRI *mri_dst,int type,unsigned char val,unsigned long *NbVoxels);
MRIS *MRISmatchSurfaceToLabel(MRIS *mris,MRI *mri_seg,int label,MRI_REGION *mri_region,INTEGRATION_PARMS *integration_parms,int connectivity);
MRIS *MRISloadSurfSubject(char *subj, char *hemi, char *surfid,
                          char *SUBJECTS_DIR);
int MRISfdr2vwth(MRIS *surf, double fdr, int signid,
                 int log10flag, int maskflag, double *vwth);

int MRISfwhm2niters(double fwhm, MRIS *surf);
int MRISfwhm2nitersSubj(double fwhm,char *subject,char *hemi,char *surfname);
double MRISfwhmFromAR1(MRIS *surf, double ar1);
int MRISscale(MRIS *mris, double scale);
int MRISseg2annot(MRIS *mris, MRI *surfseg, COLOR_TABLE *ctab);
MRI *MRISannotIndex2Seg(MRIS *mris);

#endif
