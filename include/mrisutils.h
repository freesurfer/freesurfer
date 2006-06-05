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
