/**
 * @file  mrisutils.h
 * @brief more surface processing utils
 *
 */
/*
 * Original Authors: Segonne and Greve 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/10/03 21:28:35 $
 *    $Revision: 1.21 $
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

#include "mri.h"
#include "mrisurf.h"

#ifndef SUTILS_INCLUDED
#define SUTILS_INCLUDED

void MRISsmoothSurface(MRI_SURFACE *mris,int niter,float step);
void MRISsmoothSurface2(MRI_SURFACE *mris,int niter,float step,int avrg);
void MRIScenterCOG(MRI_SURFACE *mris);
void MRIScenterCOG2(MRI_SURFACE *mris,double *xCOG,double *yCOG,double *zCOG);
MRI* MRISpeelVolume(MRIS *mris,
                    MRI *mri_src,
                    MRI *mri_dst,
                    int type,
                    unsigned char val,
                    unsigned long *NbVoxels);
MRIS *MRISmatchSurfaceToLabel(MRIS *mris,
                              MRI *mri_seg,
                              int label,
                              MRI_REGION *mri_region,
                              INTEGRATION_PARMS *integration_parms,
                              int connectivity);
MRIS *MRISloadSurfSubject(char *subj, char *hemi, char *surfid,
                          char *SUBJECTS_DIR);
int MRISfdr2vwth(MRIS *surf, double fdr, int signid,
                 int log10flag, int maskflag, double *vwth);

int    MRISfwhm2niters(double fwhm, MRIS *surf);
double MRISniters2fwhm(int niters, MRIS *surf);
int MRISfwhm2nitersSubj(double fwhm,char *subject,char *hemi,char *surfname);
double MRISfwhmFromAR1(MRIS *surf, double ar1);
int MRISscale(MRIS *mris, double scale);
int MRISseg2annot(MRIS *mris, MRI *surfseg, COLOR_TABLE *ctab);
MRI *MRISannotIndex2Seg(MRIS *mris);
double MRISvolumeInSurf(MRIS *mris);

LABEL *MRIScortexLabel(MRI_SURFACE *mris, MRI *mri_aseg, int min_vertices);
int MRISripZeros(MRIS *surf, MRI *mri);
int MRISsphericalCoords(MRIS *surf);
int MRISfindPath ( int *vert_vno, int num_vno, int max_path_length,
                   int *path, int *path_length, MRIS *mris );

double *
MRISsampleProfile(MRI_SURFACE *mris, MRI *mri, double *profile, int nsamples, int wm_samples, int outside_samples,
		  double wx,  double wy,  double wz,
		  double l4x, double l4y, double l4z,
		  double px,  double py,  double pz) ;
MRI_SP *MRISmakeTemplate(int nsubjects, char **subjlist, int nhemis, char **hemilist, char *surfregname);
#endif
