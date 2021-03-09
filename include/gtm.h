/**
 * @brief Routines to create and analyze the Geometric Transfer Matrix (GTM)
 *
 */
/*
 * Original Author: Douglas N. Greve
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef GTM_INC
#define GTM_INC

#include "matrix.h"
#include "mri.h"
#include "transform.h"
#include "mrisurf.h"

#ifdef X
#undef X
#endif

typedef struct
{
  char *subject;
  int USF; // used internally
  int OutputUSF; 
  const char *apasfile;
  const char *ctxannotfile;
  int ctxlhbase,ctxrhbase;
  MRIS *lhw, *lhp, *rhw, *rhp;
  int KeepHypo;
  int KeepCC;
  int SubSegWM;
  const char *wmannotfile;
  int wmlhbase,wmrhbase;
  float dmax;
  int nlist,srclist[300],targlist[300];
  MRI *seg;
  LTA *anat2seg;
  int lhmin, lhmax, rhmin, rhmax;
  int *segidlist, nsegs;
  MRI *anat; // header of anatomical to get geom
} GTMSEG;

typedef struct 
{
  char *name; // contrast name
  MATRIX *C;
  MATRIX *gamma;
  MATRIX *gammavar;
  MATRIX *t;
  MATRIX *p;
} GTM_CONTRAST, GTMCON;

typedef struct 
{
  int nrad;
  double Xthresh;
  MRI *res; // residual at voxel
  MRI *rvar; // residual variance across neighborhood
} LGTM;

typedef struct 
{
  MRI *yvol; // source (PET) data
  MRI *anatseg; // segmentation in anatomical space, may be hires (ie, not conformed)
  LTA *anat2pet; // registration from anat to pet space
  LTA *anat2seg; // registration from FS conformed space to seg
  LTA *seg2anat; // inverse registration
  LTA *seg2pet;  // registration from segmentation to source volume
  MRI *mask;     // binary mask in PET space
  double cFWHM, rFWHM, sFWHM; // assumed FWHM of PSF
  double cStd, rStd, sStd; // PSF FWHM converted to standard dev
  int UseMBrad;// = 1 for radial blurring
  MB2D *mbrad; // for radial/motion blurring
  int UseMBtan;// = 1 for tangential blurring
  MB2D *mbtan; // for tangential/motion blurring
  int nframes; // one place to get the number of input frames
  MRI *anatconf;   // header of conformed anatomical

  double automask_fwhm,automask_thresh; // Use this FWHM instead of PSF when computing mask
  int reduce_fov; // Flag: reduce PET FoV to be tight to mask.
  MRI_REGION *automaskRegion; // Keep reduction region in case reinstate at original FoV
  MRI *yvol_full_fov; // Used to keep header of original FoV source data
  double PadThresh; // Used to dilate mask based on automask_fwhm
  int nPad; // Actual amount of dilation

  int rescale,n_scale_refids,scale_refids[100];
  double scale, scale_refval; //scale = scale factor, scalerefval = 1 or 100, etc

  MRI *segpvf; // PVF of each seg (one in each frame)
  MRI *ttpvf; // PVF of each tissue type (one in each frame), used by MG
  MRI *gtmseg; // Segmentation in PET space
  int nmask; // number of voxels in the mask

  int nsegs,*segidlist; // number of segments in segmentation, list of segids
  int *nperseg; // number of voxels per seg measured in pet space
  MATRIX *nvox; // same as nperseg, doh
  MRI *volperseg; // volume of each segment (nperseg * voxsize) measure in anat seg space
  COLOR_TABLE *ctGTMSeg; // color table of segments
  int nReplace , SrcReplace[1000], TrgReplace[1000]; // for replacing segs
  MATRIX *vrf;  // variance reduction factor for each seg
  MATRIX *segrvar; // residual variance in each seg
  MATRIX *ttpct; // percent of the signal in each seg from each tt

  // GLM stuff for GTM
  MATRIX *X,*X0;
  MATRIX *y, *XtX, *iXtX, *Xty, *beta, *res, *yhat,*betavar;
  MATRIX *rvar,*rvargm,*rvarbrain,*rvarUnscaled; // residual variance, all vox and only GM
  MATRIX *som; // spillover matrix
  int dof;
  double XtXcond;
  MRI *ysynth,*ysynthsm; // synthesized vs yhat
  MATRIX *skew,*kurtosis; // for QAe
  int DoVoxFracCor; // Flag to correct for volume fraction effect

  int nContrasts;
  GTMCON *contrasts[100];

  MATRIX *glob_gm, *glob_gmwm, *glob_gmwmcsf; // global means

  MRI *rbv; // RBV computed volume
  MRI *rbvsegmean; // seg mean in RBV, used for QA
  MRI *rbvseg; // may be different than anatseg if rbvsegres used
  LTA *rbvseg2pet; // may be different than seg2pet if rbvsegres used
  LTA *anat2rbv; // conformed anat to RBV for convenience with mri_vol2vol
  LTA *pet2bbpet; // maps from original pet space to boundingbox
  double rbvsegres;
  int mask_rbv_to_brain; // Reduce FoV of RBV to be tight to brain
  MRI *rbvsegmasked; 

  int DoMGPVC; // Muller-Gartner
  MRI *mg; // MG output volume
  double mg_gmthresh; // GM PVF threshold
  int n_mg_refids,mg_refids[100]; // WM reference seg IDs
  MATRIX *mg_reftac; // WM reference TAC
  MRI *gmpvfpsf; // GM PVF smoothed by PSF
  int DoMGXPVC; // Muller-Gartner X 
  double mgx_gmthresh; // GM PVF threshold

  int DoMeltzerPVC; // Meltzer method
  MRI *meltzer; // MG output volume
  double MeltzerBinThresh; // Binarize GM+WM PVF before smoothing, 0 turns off bin
  double MeltzerMaskThresh; // Post-smoothing GM+WM PVF mask threshold 
  int MeltzerNDil; // Dilate binarized mask
  MRI *mzseg; // means of MZ over each seg

  int DoLGTMPVC; // Local GTM
  LGTM *lgtm;

  int DoKMRef; // Kinetic Modeling Reference TAC
  int n_km_refids,km_refids[100]; // KM reference seg IDs
  MATRIX *km_reftac; // KM reference TAC

  int DoKMHB; // Kinetic Modeling HiBinding TAC
  int n_km_hbids,km_hbids[100];// KM HiBinding seg IDs
  MATRIX *km_hbtac; // KM HiBinding TAC

  //Steady-state blood plasma concentration, unit scale, decay correction factor
  double ss_bpc, ss_scale, ss_dcf;
  int DoSteadyState; // flag, set to 1 to do steady state analysis

  char *OutDir, *AuxDir; // output folder, auxiliary output folder
  FILE *logfp;  // log file pointer

  int Optimizing; // set to 1 if currently optimzing
} GTM;


int MRIgtmSeg(GTMSEG *gtmseg);
int GTMSEGprint(GTMSEG *gtmseg, FILE *fp);
int GTMdefaultSegReplacmentList(int *nReplace, int *ReplaceThis, int *WithThat);
int GTMoptSegReplacmentList(int *nReplace, int *ReplaceThis, int *WithThat);
COLOR_TABLE *GTMSEGctab(GTMSEG *gtmseg, COLOR_TABLE *ctSubCort);

GTM *GTMalloc();
int GTMfree(GTM **pGTM);
int GTMmatrixY(GTM *gtm);
int GTMsetNMask(GTM *gtm);
int GTMpsfStd(GTM *gtm);
int GTMsegidlist(GTM *gtm);
int GTMnPad(GTM *gtm);
int GTMbuildX(GTM *gtm);
int GTMsolve(GTM *gtm);
int GTMsegrvar(GTM *gtm);
int GTMsynth(GTM *gtm, int NoiseSeed, int nReps);
int GTMsmoothSynth(GTM *gtm);
MRI *GTMsegSynth(GTM *gtm, int frame, MRI *synth);
int GTMrbvseg(GTM *gtm);
int GTMrbv(GTM *gtm);
int GTMmgRefTAC(GTM *gtm);
int GTMmgpvc(GTM *gtm);
MRI *GTMmgxpvc(GTM *gtm, int Target);
int GTMmeltzerpvc(GTM *gtm);
MATRIX *GTMvol2mat(GTM *gtm, MRI *vol, MATRIX *m);
MRI *GTMmat2vol(GTM *gtm, MATRIX *m, MRI *vol);
int GTMprintReplaceList(FILE *fp, const int nReplace, const int *ReplaceThis, const int *WithThat);
int GTMcheckReplaceList(const int nReplace, const int *ReplaceThis, const int *WithThat);
int GTMloadReplacmentList(const char *fname, int *nReplace, int *ReplaceThis, int *WithThat);
int GTMcheckX(MATRIX *X);
int GTMautoMask(GTM *gtm);
int GTMrvarGM(GTM *gtm);
int GTMttest(GTM *gtm);
int GTMmgRefIds(GTM *gtm);
int GTMprintMGRefTAC(GTM *gtm, FILE *fp);
int GTMwriteMGRefTAC(GTM *gtm, char *filename);
int GTMrescale(GTM *gtm);
int GTMsteadyState(GTM *gtm);
int GTMwriteContrasts(GTM *GTM);
int GTMprintRefIds(GTM *gtm, FILE *fp);
int GTMcheckRefIds(GTM *gtm);
int GTMrefTAC(GTM *gtm);
int VRFStats(GTM *gtm, double *vrfmean, double *vrfmin, double *vrfmax);
int WriteVRFStats(char *fname, GTM *gtm);
int GTMglobalStats(GTM *gtm);
MRI **GTMlocal(GTM *gtm, MRI **pvc);
int GTMttPercent(GTM *gtm);
int GTMsom(GTM *gtm);
int GTMsegid2nthseg(GTM *gtm, int segid);
int GTMwriteText(GTM *gtm, char *OutDir, int DeMean);

#endif
