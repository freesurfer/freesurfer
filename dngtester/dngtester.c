/**
 * @file  dngtester.c
 * @brief dougs super special test code
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2014/03/03 02:59:55 $
 *    $Revision: 1.53 $
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "mrisurf.h"
#include "utils.h"
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "region.h"
#include "machine.h"
#include "fio.h"
#include "mri_identify.h"
#include "mrisurf.h"
#include "fmriutils.h"
#include "gca.h"
#include "gcsa.h"
#include "fsgdf.h"
#include "icosahedron.h"
#include "gca.h"
#include "gcamorph.h"
#include "DICOMRead.h"
#include "fsenv.h"
#include "qdecutils.h"
#include "dti.h"
#include "registerio.h"
#include "timer.h"
#include "evschutils.h"
#include "matrix.h"
#include "matfile.h"
#include "randomfields.h"
#include "mri2.h"
#include "annotation.h"
#include "mrisutils.h"
#include "image.h"
#include "retinotopy.h"
#include "bfileio.h"
#include "cma.h"
#include "transform.h"
#include "mrinorm.h"
#include "resample.h"
#ifdef _OPENMP
#include <omp.h>
#endif
//#include "dbfa.h"
#include "cpputils.h"

double round(double);

// setenv SUBJECTS_DIR /space/greve/1/users/greve/subjects
// /autofs/space/greve_001/users/greve/dev/trunk/dngtester
// ./dngtester register.dat func.mgz func-in-m3z.mgh
// ./dngtester identity.dat ~/subjects/fbirn-anat-101/mri/orig.mgz orig-in-m3z.mgh
// tkregister2 --targ orig-in-m3z.mgh --mov func-in-m3z.mgh --regheader --reg tmp.reg
// Atlas: $SUBJECTS_DIR/avgtst/mri/T1MLV.mgz

MRI *MRIsetSliceNo(MRI *mri, MRI *out);

LTA *TransformRegDat2LTA(MRI *targ, MRI *mov, MATRIX *R);
char *Progname = "dngtester";

char *outvolfile=NULL;
char *subject = NULL;
char *regfile = NULL;
TRANSFORM *Rtransform;  //types : M3D, M3Z, LTA, FSLMAT, DAT, OCT(TA), XFM
FSENV     *fsenv;
char gcamfile[1000];
char origfile[1000];
char *sourcefile;
GCAM      *gcam;
MRI *mri, *mri2, *mri3, *mri4, *mriref, *mask=NULL, *mriarray[100];
MATRIX *V, *W, *m_tmp;
float ipr, bpr, intensity;
MATRIX *R;
int float2int;
int err,n;
struct timeb  mytimer;
int msecFitTime;
double a,b,zthresh,pthresh;
char *SUBJECTS_DIR;
char tmpstr[2000];
MRIS *surf, *surf2, *sphere;
char *hemi;

MRI *fMRIvariance2(MRI *fmri, float DOF, int RmMean, MRI *var);
void printrgb(void);

COLOR_TABLE *CTABaddEntry(COLOR_TABLE *ctold, const char *name);
int MRISmercator(MRIS *surf);
IMAGE *I;
MHT *lhwhite_hash;
MRIS *lhwhite,*white,*sphere;
MRIS *surfs[100], *mris, *mris1, *mris2;
int *XNbrVtxNo, nXNbrs;
double *XNbrDotProd;

MRIS *MRISaverageSurfaces(int nsurfaces, MRIS **surfs);
MATRIX *MatrixLoadFSL(char *fname);
int find_path ( int* vert_vno, int num_vno, int max_path_length,
                int* path, int* path_length, MRIS *mris );
MRI *seg1, *seg2;
int mygd(char *fname);
MRI *MRIaddB(MRI *mri1, MRI *mri2, MRI *mriadd);
MRI *MRInnFill(MRI *src, MRI *seed, MRI *out, MRI *mask);
int mrisComputeDistanceTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
MRI *surf2surf_nnfrM(MRI *SrcSurfVals, MRI_SURFACE **SurfReg, int nsurfs,
		     int ReverseMapFlag, int DoJac, int UseHash);
MRI *MRISmercatorGrid(MRIS *mris, double dtheta, double dphi);
int GCAprint(GCA *gca, FILE *fp);
MRI * GCAmri2(GCA *gca, MRI *mri);
MRI *GCAnlables(GCA *gca, MRI *mri);
MRI *MRIpolyfitBiasField(MRI *vol, int order, MRI *mask, MRI *bias);
int MRIsegCount(MRI *seg, int id, int frame);
float MRIvoxelsInLabelWithPartialVolumeEffects2( const MRI *mri,
						const MRI *mri_vals, 
						const int label, 
						MRI *mri_mixing_coef, 
						 MRI *mri_nbr_labels );
int ProjSurf(MRIS *surf, MRI *d, double f);
MRI_SP  *mrisp, /* *mrisp_aligned,*/ *mrisp_template ;
MRI_SP *MRISmakeTemplateB(int nsubjects, char **subjlist, int nhemis, char **hemilist, char *surfregname);
MRI_SP *MRISmakeTemplateC(int nsurfaces, char **surfregflist, char **annotflist, 
			  int nparams, char **paramflist, int navgs);

MATRIX *MatrixSVD2(MATRIX *mA, VECTOR *v_z, MATRIX *mV)
{
  if(mA->rows < mA->cols){
    MATRIX *Q, *sQ=NULL, *VQ=NULL;
    Q = MatrixTranspose(mA,NULL);
    MatrixSVD2(Q, sQ, VQ); // Q becomes U
    mA = MatrixCopy(VQ,mA);
    v_z = MatrixTranspose(sQ,v_z);
    mV = MatrixCopy(Q,mV);
    MatrixFree(&Q);
    MatrixFree(&sQ);
    MatrixFree(&VQ);
  }

  if(mV==NULL)  mV = MatrixAlloc(mA->rows,mA->rows,MATRIX_REAL);
  if(v_z==NULL) v_z = MatrixAlloc(mA->rows,mA->cols,MATRIX_REAL);
  OpenSvdcmp ( mA, v_z, mV ) ;

  return ( mV ) ;
}

MATRIX *gdfContrastDOSSDeleteMe(FSGD *fsgd, float *wClass, float *wCovar) {
  MATRIX *C;
  float w;
  int c, v, n;

  if(strcasecmp(fsgd->DesignMatMethod,"doss") != 0){
    printf("ERROR: gdfContrastDOSS() cannot be used with %s\n",fsgd->DesignMatMethod);
    //return(NULL);
  }

  /* Contrast matrix*/
  C = MatrixAlloc(1, fsgd->nvariables+fsgd->nclasses, MATRIX_REAL);

  n = 0;
  for (c=0; c < fsgd->nclasses; c++) {
    w = 1;
    if (wClass != NULL) w *= wClass[c];
    C->rptr[1][n+1] = w;
    n++;
  }
  for (v=0; v < fsgd->nvariables+1; v++) {
    w = 1;
    if (wCovar != NULL) w *= wCovar[v];
    C->rptr[1][n+1] = w;
    n++;
  }
  return(C);
}

FSGD *LoadFactorFile(char *fname);
int FSGDFactorsPrint(FSGD *gd, FILE *fp);
MRI *MRIsurf2VolOpt(MRI *tempvol, MRIS **surfs, MRI **overlays, int nsurfs, 
		    MRI *ribbon, MATRIX *R, MRI *volsurf);
MRI *FillSeg(MRI *seg, COLOR_TABLE *ctab, MRI *vol);
MRI *MRIdownSmoothUp(MRI *src, int Fc, int Fr, int Fs, 
		     double cFWHM, double rFWHM, double sFWHM, 
		     MRI *dst);
MRI *MRIfcIntrinsicLI(MRI *lh, MRI *rh, double DenThresh);
MRI *NNGLMPVC(MRI *src, MRI *pvf, MRI *mask, int nrad, MRI *pvc);
int MRISfillInteriorTest(char *subject, int UseNew, FILE *fp);

/*----------------------------------------*/
int main(int argc, char **argv) 
{
  int err,area32p,area32v,superiorfrontal,medialorbitofrontal;
  int rostralanteriorcingulate, rostralmiddlefrontal;
  int index,k,c,r,s,nvox,a1,a2,k1,k2,nmask,nhits,nsup;
  int *nunits;
  char *parcnames[10], *annot1, *annot2;
  COLOR_TABLE *ct ;
  VERTEX *vtx,*vtx1w,*vtx2w,*vtx1s,*vtx2s;
  float dlhw,DotProdThresh;
  int  lhwvtx;
  double sumval;
  int nsegid1, *segidlist1;
  int nsegid2, *segidlist2;
  int nsegs, *segidlist,vtxno1,vtxno2;
  double *area1, *area2, *area12, *dice;
  double f,radius,radius2,DotProd,theta,d2,d3,d3Sqr;
  double fwhm,fwhmSqr,*stats;
  COLOR_TABLE *ctab = NULL;
  int stg, ttg, bsts, mtg, insula;
  MRI_REGION region;
  int msec;
  GCA *gca;
  float inplaneres, betplaneres, intensity;
  int float2int;
  SURFCLUSTERSUM *SurfClustList;
  int nClusters;
  MATRIX *vox2ras, *Qsrc, *C;
  FILE *fp, *fp2;
  char *subjlist[10], *hemilist[2], *surfregflist[100], *annotflist[100], *paramflist[300];
  float *wClass, *wCovar;
  FSGD *fsgd;
  int FLevels[10], FactorList[10];
  int nFactors, nFactorList;

  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  Gdiag = DIAG_WRITE | DIAG_VERBOSE_ON;

  fp = fopen(argv[1],"r"); // slist
  fp2 = fopen(argv[2],"w"); // outfile
  while(fp){
    fscanf(fp,"%s",tmpstr);
    MRISfillInteriorRibbonTest(tmpstr, 1, fp2);
  }
  fclose(fp2);

  exit(0);

  //---------------------------------------------------------

  mri = MRIread(argv[1]); // lh
  mri2 = MRIfisherTransform(mri, NULL, NULL);
  MRIwrite(mri2,"ft.nii");
  exit(0);


  //---------------------------------------------------------
  if(argc != 4){
    printf("add_geometry_to_surf volume insurface outsurface\n");
    exit(1);
  }
  mri = MRIread(argv[1]);
  if(mri == NULL) exit(1);
  mris = MRISread(argv[2]);
  if(mris == NULL) exit(1);
  initVolGeom(&mris->vg);
  getVolGeom(mri, &mris->vg);
  err = MRISwrite(mris,argv[3]);
  if(err) exit(1);

  exit(0);

  mri = MRIread(argv[1]); // seg
  ctab = CTABreadASCII(argv[2]);
  mri2 = FillSeg(mri, ctab, NULL);
  MRIwrite(mri2,"segfill.mgh");

  exit(0);

  mri = MRIread(argv[1]); // template
  mri2 = MRIread(argv[2]); // ribbon
  regio_read_register(argv[3], &subject, &inplaneres, &betplaneres, &intensity,  &R, &float2int);
  surfs[0]    = MRISread(argv[4]); // surf
  mriarray[0] = MRIread(argv[5]); // overlay
  surfs[1]    = MRISread(argv[6]); // surf
  mriarray[1] = MRIread(argv[7]); // overlay
  mri4 = MRIsurf2VolOpt(mri, surfs, mriarray, 2, mri2, R, NULL);
  MRIwrite(mri4,"mri4.nii");

  exit(0);

  // FSGD ---------------------------------------
  fsgd = LoadFactorFile(argv[1]);
  printf("----------------------------------\n");
  FSGDFactorsPrint(fsgd, stdout);
  exit(0);

  for(nFactors = 1; nFactors < 4; nFactors++){
    for(k=0;k<nFactors;k++)  FLevels[k] = k+2;
    for(nFactorList = 1; nFactorList <= nFactors; nFactorList++){
      for(k=0;k<nFactorList;k++)  FactorList[k] = k+1;
      C = ANOVAContrast(FLevels, nFactors, FactorList, nFactorList);
      printf("nF=%d nFL=%d --------------------\n",nFactors,nFactorList);
      printf("c = fast_anovamtx([");
      for(k=0;k<nFactors;k++) printf("%d ",FLevels[k]);
      printf("],[");
      for(k=0;k<nFactorList;k++) printf("%d ",FactorList[k]);
      printf("])\n");
      MatrixPrint(stdout,C);
    }
  }

  exit(0);


  exit(0);
  fsgd = gdfRead(argv[1],0);
  //strcpy(fsgd->DesignMatMethod,"dods");
  wClass = (float *)calloc(10,sizeof(float));
  wClass[0] = 0;
  wClass[1] = 0;
  wCovar = (float *)calloc(10,sizeof(float));
  wCovar[0] = 1;
  wCovar[1] = 0;
  //C = gdfContrastDODS(fsgd, wClass, wCovar);
  C = gdfContrastDOSS(fsgd, wClass, wCovar);
  printf("---------------------\n");
  MatrixPrint(stdout,C);
  printf("---------------------\n");

  exit(0);

  mri  = MRIread(argv[1]);
  mris = MRISread(argv[2]);
  mri3 = MRIread(argv[3]);
  Gdiag_no = 1;
  mri2 = MRISsmoothMRIFastD(mris, mri, 50, mri3,NULL);
  MRIwrite(mri2,"new.mgh");
  mri2 = MRISsmoothMRIFast(mris, mri, 50, mri3,NULL);
  MRIwrite(mri2,"old.mgh");

  exit(0);

  int srclist[20], targlist[20], nlist=11;
  srclist[0]  = 4;  targlist[0] = 24;
  srclist[1]  = 5;  targlist[1] = 24;
  srclist[2]  = 14;  targlist[2] = 24;
  srclist[3]  = 15;  targlist[3] = 24;
  srclist[4]  = 31;  targlist[4] = 24;
  srclist[5]  = 43;  targlist[5] = 24;
  srclist[6]  = 44;  targlist[6] = 24;
  srclist[7]  = 72;  targlist[7] = 24;
  srclist[8]  = 63;  targlist[8] = 24;
  srclist[9]  = 30;  targlist[9] = 24;
  srclist[10]  = 62;  targlist[10] = 24;
  mri2 = MRIreplaceList(mri, srclist, targlist, nlist, NULL);
  MRIwrite(mri2,"newseg.mgh");

  exit(0);

  R = MatrixAlloc(2,1,MATRIX_REAL);
  R->rptr[1][1] = 1;
  R->rptr[1][2] = .5;
  //R->rptr[2][1] = 4;
  //R->rptr[2][2] = 5;
  printf("R ------------------\n");
  MatrixPrint(stdout,R);
  vox2ras = NULL; //MatrixAlloc(2,1,MATRIX_REAL);
  Qsrc = NULL; //MatrixAlloc(2,2,MATRIX_REAL);
  MatrixSVD2(R, vox2ras, Qsrc);
  printf("U ------------------\n");
  MatrixPrint(stdout,R);
  printf("S ------------------\n");
  MatrixPrint(stdout,vox2ras);
  printf("V ------------------\n");
  MatrixPrint(stdout,Qsrc);
  exit(1);

surfregflist[0] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/lh.fssym.i00.reg";
annotflist[0] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/lh.aparc.annot";
paramflist[0] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/lh.inflated.H";
paramflist[1] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/lh.sulc";
paramflist[2] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/lh.smoothwm.H";
surfregflist[1] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/rh.fssym.i00.reg";
annotflist[1] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/rh.aparc.annot";
paramflist[3] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/rh.inflated.H";
paramflist[4] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/rh.sulc";
paramflist[5] = "/cluster/con/9/users/greve/fsaverage_sym-data/073/surf/lrreg/rh.smoothwm.H";
surfregflist[2] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/lh.fssym.i00.reg";
annotflist[2] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/lh.aparc.annot";
paramflist[6] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/lh.inflated.H";
paramflist[7] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/lh.sulc";
paramflist[8] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/lh.smoothwm.H";
surfregflist[3] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/rh.fssym.i00.reg";
annotflist[3] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/rh.aparc.annot";
paramflist[9] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/rh.inflated.H";
paramflist[10] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/rh.sulc";
paramflist[11] = "/cluster/con/9/users/greve/fsaverage_sym-data/080/surf/lrreg/rh.smoothwm.H";
surfregflist[4] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/lh.fssym.i00.reg";
annotflist[4] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/lh.aparc.annot";
paramflist[12] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/lh.inflated.H";
paramflist[13] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/lh.sulc";
paramflist[14] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/lh.smoothwm.H";
surfregflist[5] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/rh.fssym.i00.reg";
annotflist[5] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/rh.aparc.annot";
paramflist[15] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/rh.inflated.H";
paramflist[16] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/rh.sulc";
paramflist[17] = "/cluster/con/9/users/greve/fsaverage_sym-data/123/surf/lrreg/rh.smoothwm.H";
surfregflist[6] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/lh.fssym.i00.reg";
annotflist[6] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/lh.aparc.annot";
paramflist[18] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/lh.inflated.H";
paramflist[19] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/lh.sulc";
paramflist[20] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/lh.smoothwm.H";
surfregflist[7] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/rh.fssym.i00.reg";
annotflist[7] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/rh.aparc.annot";
paramflist[21] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/rh.inflated.H";
paramflist[22] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/rh.sulc";
paramflist[23] = "/cluster/con/9/users/greve/fsaverage_sym-data/144/surf/lrreg/rh.smoothwm.H";
surfregflist[8] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/lh.fssym.i00.reg";
annotflist[8] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/lh.aparc.annot";
paramflist[24] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/lh.inflated.H";
paramflist[25] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/lh.sulc";
paramflist[26] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/lh.smoothwm.H";
surfregflist[9] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/rh.fssym.i00.reg";
annotflist[9] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/rh.aparc.annot";
paramflist[27] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/rh.inflated.H";
paramflist[28] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/rh.sulc";
paramflist[29] = "/cluster/con/9/users/greve/fsaverage_sym-data/102/surf/lrreg/rh.smoothwm.H";
surfregflist[10] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/lh.fssym.i00.reg";
annotflist[10] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/lh.aparc.annot";
paramflist[30] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/lh.inflated.H";
paramflist[31] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/lh.sulc";
paramflist[32] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/lh.smoothwm.H";
surfregflist[11] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/rh.fssym.i00.reg";
annotflist[11] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/rh.aparc.annot";
paramflist[33] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/rh.inflated.H";
paramflist[34] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/rh.sulc";
paramflist[35] = "/cluster/con/9/users/greve/fsaverage_sym-data/140/surf/lrreg/rh.smoothwm.H";
surfregflist[12] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/lh.fssym.i00.reg";
annotflist[12] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/lh.aparc.annot";
paramflist[36] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/lh.inflated.H";
paramflist[37] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/lh.sulc";
paramflist[38] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/lh.smoothwm.H";
surfregflist[13] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/rh.fssym.i00.reg";
annotflist[13] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/rh.aparc.annot";
paramflist[39] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/rh.inflated.H";
paramflist[40] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/rh.sulc";
paramflist[41] = "/cluster/con/9/users/greve/fsaverage_sym-data/084/surf/lrreg/rh.smoothwm.H";
surfregflist[14] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/lh.fssym.i00.reg";
annotflist[14] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/lh.aparc.annot";
paramflist[42] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/lh.inflated.H";
paramflist[43] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/lh.sulc";
paramflist[44] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/lh.smoothwm.H";
surfregflist[15] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/rh.fssym.i00.reg";
annotflist[15] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/rh.aparc.annot";
paramflist[45] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/rh.inflated.H";
paramflist[46] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/rh.sulc";
paramflist[47] = "/cluster/con/9/users/greve/fsaverage_sym-data/145/surf/lrreg/rh.smoothwm.H";
surfregflist[16] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/lh.fssym.i00.reg";
annotflist[16] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/lh.aparc.annot";
paramflist[48] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/lh.inflated.H";
paramflist[49] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/lh.sulc";
paramflist[50] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/lh.smoothwm.H";
surfregflist[17] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/rh.fssym.i00.reg";
annotflist[17] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/rh.aparc.annot";
paramflist[51] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/rh.inflated.H";
paramflist[52] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/rh.sulc";
paramflist[53] = "/cluster/con/9/users/greve/fsaverage_sym-data/092/surf/lrreg/rh.smoothwm.H";
surfregflist[18] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/lh.fssym.i00.reg";
annotflist[18] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/lh.aparc.annot";
paramflist[54] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/lh.inflated.H";
paramflist[55] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/lh.sulc";
paramflist[56] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/lh.smoothwm.H";
surfregflist[19] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/rh.fssym.i00.reg";
annotflist[19] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/rh.aparc.annot";
paramflist[57] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/rh.inflated.H";
paramflist[58] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/rh.sulc";
paramflist[59] = "/cluster/con/9/users/greve/fsaverage_sym-data/004/surf/lrreg/rh.smoothwm.H";
surfregflist[20] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/lh.fssym.i00.reg";
annotflist[20] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/lh.aparc.annot";
paramflist[60] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/lh.inflated.H";
paramflist[61] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/lh.sulc";
paramflist[62] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/lh.smoothwm.H";
surfregflist[21] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/rh.fssym.i00.reg";
annotflist[21] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/rh.aparc.annot";
paramflist[63] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/rh.inflated.H";
paramflist[64] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/rh.sulc";
paramflist[65] = "/cluster/con/9/users/greve/fsaverage_sym-data/008/surf/lrreg/rh.smoothwm.H";
surfregflist[22] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/lh.fssym.i00.reg";
annotflist[22] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/lh.aparc.annot";
paramflist[66] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/lh.inflated.H";
paramflist[67] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/lh.sulc";
paramflist[68] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/lh.smoothwm.H";
surfregflist[23] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/rh.fssym.i00.reg";
annotflist[23] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/rh.aparc.annot";
paramflist[69] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/rh.inflated.H";
paramflist[70] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/rh.sulc";
paramflist[71] = "/cluster/con/9/users/greve/fsaverage_sym-data/017/surf/lrreg/rh.smoothwm.H";
surfregflist[24] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/lh.fssym.i00.reg";
annotflist[24] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/lh.aparc.annot";
paramflist[72] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/lh.inflated.H";
paramflist[73] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/lh.sulc";
paramflist[74] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/lh.smoothwm.H";
surfregflist[25] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/rh.fssym.i00.reg";
annotflist[25] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/rh.aparc.annot";
paramflist[75] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/rh.inflated.H";
paramflist[76] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/rh.sulc";
paramflist[77] = "/cluster/con/9/users/greve/fsaverage_sym-data/021/surf/lrreg/rh.smoothwm.H";
surfregflist[26] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/lh.fssym.i00.reg";
annotflist[26] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/lh.aparc.annot";
paramflist[78] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/lh.inflated.H";
paramflist[79] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/lh.sulc";
paramflist[80] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/lh.smoothwm.H";
surfregflist[27] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/rh.fssym.i00.reg";
annotflist[27] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/rh.aparc.annot";
paramflist[81] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/rh.inflated.H";
paramflist[82] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/rh.sulc";
paramflist[83] = "/cluster/con/9/users/greve/fsaverage_sym-data/032/surf/lrreg/rh.smoothwm.H";
surfregflist[28] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/lh.fssym.i00.reg";
annotflist[28] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/lh.aparc.annot";
paramflist[84] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/lh.inflated.H";
paramflist[85] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/lh.sulc";
paramflist[86] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/lh.smoothwm.H";
surfregflist[29] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/rh.fssym.i00.reg";
annotflist[29] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/rh.aparc.annot";
paramflist[87] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/rh.inflated.H";
paramflist[88] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/rh.sulc";
paramflist[89] = "/cluster/con/9/users/greve/fsaverage_sym-data/039/surf/lrreg/rh.smoothwm.H";
surfregflist[30] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/lh.fssym.i00.reg";
annotflist[30] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/lh.aparc.annot";
paramflist[90] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/lh.inflated.H";
paramflist[91] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/lh.sulc";
paramflist[92] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/lh.smoothwm.H";
surfregflist[31] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/rh.fssym.i00.reg";
annotflist[31] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/rh.aparc.annot";
paramflist[93] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/rh.inflated.H";
paramflist[94] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/rh.sulc";
paramflist[95] = "/cluster/con/9/users/greve/fsaverage_sym-data/040/surf/lrreg/rh.smoothwm.H";
surfregflist[32] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/lh.fssym.i00.reg";
annotflist[32] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/lh.aparc.annot";
paramflist[96] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/lh.inflated.H";
paramflist[97] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/lh.sulc";
paramflist[98] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/lh.smoothwm.H";
surfregflist[33] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/rh.fssym.i00.reg";
annotflist[33] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/rh.aparc.annot";
paramflist[99] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/rh.inflated.H";
paramflist[100] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/rh.sulc";
paramflist[101] = "/cluster/con/9/users/greve/fsaverage_sym-data/045/surf/lrreg/rh.smoothwm.H";
surfregflist[34] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/lh.fssym.i00.reg";
annotflist[34] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/lh.aparc.annot";
paramflist[102] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/lh.inflated.H";
paramflist[103] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/lh.sulc";
paramflist[104] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/lh.smoothwm.H";
surfregflist[35] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/rh.fssym.i00.reg";
annotflist[35] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/rh.aparc.annot";
paramflist[105] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/rh.inflated.H";
paramflist[106] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/rh.sulc";
paramflist[107] = "/cluster/con/9/users/greve/fsaverage_sym-data/049/surf/lrreg/rh.smoothwm.H";
surfregflist[36] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/lh.fssym.i00.reg";
annotflist[36] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/lh.aparc.annot";
paramflist[108] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/lh.inflated.H";
paramflist[109] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/lh.sulc";
paramflist[110] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/lh.smoothwm.H";
surfregflist[37] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/rh.fssym.i00.reg";
annotflist[37] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/rh.aparc.annot";
paramflist[111] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/rh.inflated.H";
paramflist[112] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/rh.sulc";
paramflist[113] = "/cluster/con/9/users/greve/fsaverage_sym-data/067/surf/lrreg/rh.smoothwm.H";
surfregflist[38] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/lh.fssym.i00.reg";
annotflist[38] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/lh.aparc.annot";
paramflist[114] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/lh.inflated.H";
paramflist[115] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/lh.sulc";
paramflist[116] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/lh.smoothwm.H";
surfregflist[39] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/rh.fssym.i00.reg";
annotflist[39] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/rh.aparc.annot";
paramflist[117] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/rh.inflated.H";
paramflist[118] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/rh.sulc";
paramflist[119] = "/cluster/con/9/users/greve/fsaverage_sym-data/074/surf/lrreg/rh.smoothwm.H";
surfregflist[40] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/lh.fssym.i00.reg";
annotflist[40] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/lh.aparc.annot";
paramflist[120] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/lh.inflated.H";
paramflist[121] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/lh.sulc";
paramflist[122] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/lh.smoothwm.H";
surfregflist[41] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/rh.fssym.i00.reg";
annotflist[41] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/rh.aparc.annot";
paramflist[123] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/rh.inflated.H";
paramflist[124] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/rh.sulc";
paramflist[125] = "/cluster/con/9/users/greve/fsaverage_sym-data/091/surf/lrreg/rh.smoothwm.H";
surfregflist[42] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/lh.fssym.i00.reg";
annotflist[42] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/lh.aparc.annot";
paramflist[126] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/lh.inflated.H";
paramflist[127] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/lh.sulc";
paramflist[128] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/lh.smoothwm.H";
surfregflist[43] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/rh.fssym.i00.reg";
annotflist[43] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/rh.aparc.annot";
paramflist[129] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/rh.inflated.H";
paramflist[130] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/rh.sulc";
paramflist[131] = "/cluster/con/9/users/greve/fsaverage_sym-data/093/surf/lrreg/rh.smoothwm.H";
surfregflist[44] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/lh.fssym.i00.reg";
annotflist[44] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/lh.aparc.annot";
paramflist[132] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/lh.inflated.H";
paramflist[133] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/lh.sulc";
paramflist[134] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/lh.smoothwm.H";
surfregflist[45] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/rh.fssym.i00.reg";
annotflist[45] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/rh.aparc.annot";
paramflist[135] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/rh.inflated.H";
paramflist[136] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/rh.sulc";
paramflist[137] = "/cluster/con/9/users/greve/fsaverage_sym-data/095/surf/lrreg/rh.smoothwm.H";
surfregflist[46] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/lh.fssym.i00.reg";
annotflist[46] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/lh.aparc.annot";
paramflist[138] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/lh.inflated.H";
paramflist[139] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/lh.sulc";
paramflist[140] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/lh.smoothwm.H";
surfregflist[47] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/rh.fssym.i00.reg";
annotflist[47] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/rh.aparc.annot";
paramflist[141] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/rh.inflated.H";
paramflist[142] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/rh.sulc";
paramflist[143] = "/cluster/con/9/users/greve/fsaverage_sym-data/097/surf/lrreg/rh.smoothwm.H";
surfregflist[48] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/lh.fssym.i00.reg";
annotflist[48] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/lh.aparc.annot";
paramflist[144] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/lh.inflated.H";
paramflist[145] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/lh.sulc";
paramflist[146] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/lh.smoothwm.H";
surfregflist[49] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/rh.fssym.i00.reg";
annotflist[49] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/rh.aparc.annot";
paramflist[147] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/rh.inflated.H";
paramflist[148] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/rh.sulc";
paramflist[149] = "/cluster/con/9/users/greve/fsaverage_sym-data/099/surf/lrreg/rh.smoothwm.H";
surfregflist[50] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/lh.fssym.i00.reg";
annotflist[50] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/lh.aparc.annot";
paramflist[150] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/lh.inflated.H";
paramflist[151] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/lh.sulc";
paramflist[152] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/lh.smoothwm.H";
surfregflist[51] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/rh.fssym.i00.reg";
annotflist[51] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/rh.aparc.annot";
paramflist[153] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/rh.inflated.H";
paramflist[154] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/rh.sulc";
paramflist[155] = "/cluster/con/9/users/greve/fsaverage_sym-data/103/surf/lrreg/rh.smoothwm.H";
surfregflist[52] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/lh.fssym.i00.reg";
annotflist[52] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/lh.aparc.annot";
paramflist[156] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/lh.inflated.H";
paramflist[157] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/lh.sulc";
paramflist[158] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/lh.smoothwm.H";
surfregflist[53] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/rh.fssym.i00.reg";
annotflist[53] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/rh.aparc.annot";
paramflist[159] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/rh.inflated.H";
paramflist[160] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/rh.sulc";
paramflist[161] = "/cluster/con/9/users/greve/fsaverage_sym-data/106/surf/lrreg/rh.smoothwm.H";
surfregflist[54] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/lh.fssym.i00.reg";
annotflist[54] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/lh.aparc.annot";
paramflist[162] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/lh.inflated.H";
paramflist[163] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/lh.sulc";
paramflist[164] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/lh.smoothwm.H";
surfregflist[55] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/rh.fssym.i00.reg";
annotflist[55] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/rh.aparc.annot";
paramflist[165] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/rh.inflated.H";
paramflist[166] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/rh.sulc";
paramflist[167] = "/cluster/con/9/users/greve/fsaverage_sym-data/108/surf/lrreg/rh.smoothwm.H";
surfregflist[56] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/lh.fssym.i00.reg";
annotflist[56] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/lh.aparc.annot";
paramflist[168] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/lh.inflated.H";
paramflist[169] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/lh.sulc";
paramflist[170] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/lh.smoothwm.H";
surfregflist[57] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/rh.fssym.i00.reg";
annotflist[57] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/rh.aparc.annot";
paramflist[171] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/rh.inflated.H";
paramflist[172] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/rh.sulc";
paramflist[173] = "/cluster/con/9/users/greve/fsaverage_sym-data/111/surf/lrreg/rh.smoothwm.H";
surfregflist[58] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/lh.fssym.i00.reg";
annotflist[58] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/lh.aparc.annot";
paramflist[174] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/lh.inflated.H";
paramflist[175] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/lh.sulc";
paramflist[176] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/lh.smoothwm.H";
surfregflist[59] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/rh.fssym.i00.reg";
annotflist[59] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/rh.aparc.annot";
paramflist[177] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/rh.inflated.H";
paramflist[178] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/rh.sulc";
paramflist[179] = "/cluster/con/9/users/greve/fsaverage_sym-data/114/surf/lrreg/rh.smoothwm.H";
surfregflist[60] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/lh.fssym.i00.reg";
annotflist[60] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/lh.aparc.annot";
paramflist[180] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/lh.inflated.H";
paramflist[181] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/lh.sulc";
paramflist[182] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/lh.smoothwm.H";
surfregflist[61] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/rh.fssym.i00.reg";
annotflist[61] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/rh.aparc.annot";
paramflist[183] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/rh.inflated.H";
paramflist[184] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/rh.sulc";
paramflist[185] = "/cluster/con/9/users/greve/fsaverage_sym-data/124/surf/lrreg/rh.smoothwm.H";
surfregflist[62] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/lh.fssym.i00.reg";
annotflist[62] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/lh.aparc.annot";
paramflist[186] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/lh.inflated.H";
paramflist[187] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/lh.sulc";
paramflist[188] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/lh.smoothwm.H";
surfregflist[63] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/rh.fssym.i00.reg";
annotflist[63] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/rh.aparc.annot";
paramflist[189] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/rh.inflated.H";
paramflist[190] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/rh.sulc";
paramflist[191] = "/cluster/con/9/users/greve/fsaverage_sym-data/128/surf/lrreg/rh.smoothwm.H";
surfregflist[64] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/lh.fssym.i00.reg";
annotflist[64] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/lh.aparc.annot";
paramflist[192] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/lh.inflated.H";
paramflist[193] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/lh.sulc";
paramflist[194] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/lh.smoothwm.H";
surfregflist[65] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/rh.fssym.i00.reg";
annotflist[65] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/rh.aparc.annot";
paramflist[195] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/rh.inflated.H";
paramflist[196] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/rh.sulc";
paramflist[197] = "/cluster/con/9/users/greve/fsaverage_sym-data/129/surf/lrreg/rh.smoothwm.H";
surfregflist[66] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/lh.fssym.i00.reg";
annotflist[66] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/lh.aparc.annot";
paramflist[198] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/lh.inflated.H";
paramflist[199] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/lh.sulc";
paramflist[200] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/lh.smoothwm.H";
surfregflist[67] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/rh.fssym.i00.reg";
annotflist[67] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/rh.aparc.annot";
paramflist[201] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/rh.inflated.H";
paramflist[202] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/rh.sulc";
paramflist[203] = "/cluster/con/9/users/greve/fsaverage_sym-data/130/surf/lrreg/rh.smoothwm.H";
surfregflist[68] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/lh.fssym.i00.reg";
annotflist[68] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/lh.aparc.annot";
paramflist[204] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/lh.inflated.H";
paramflist[205] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/lh.sulc";
paramflist[206] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/lh.smoothwm.H";
surfregflist[69] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/rh.fssym.i00.reg";
annotflist[69] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/rh.aparc.annot";
paramflist[207] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/rh.inflated.H";
paramflist[208] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/rh.sulc";
paramflist[209] = "/cluster/con/9/users/greve/fsaverage_sym-data/131/surf/lrreg/rh.smoothwm.H";
surfregflist[70] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/lh.fssym.i00.reg";
annotflist[70] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/lh.aparc.annot";
paramflist[210] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/lh.inflated.H";
paramflist[211] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/lh.sulc";
paramflist[212] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/lh.smoothwm.H";
surfregflist[71] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/rh.fssym.i00.reg";
annotflist[71] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/rh.aparc.annot";
paramflist[213] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/rh.inflated.H";
paramflist[214] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/rh.sulc";
paramflist[215] = "/cluster/con/9/users/greve/fsaverage_sym-data/133/surf/lrreg/rh.smoothwm.H";
surfregflist[72] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/lh.fssym.i00.reg";
annotflist[72] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/lh.aparc.annot";
paramflist[216] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/lh.inflated.H";
paramflist[217] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/lh.sulc";
paramflist[218] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/lh.smoothwm.H";
surfregflist[73] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/rh.fssym.i00.reg";
annotflist[73] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/rh.aparc.annot";
paramflist[219] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/rh.inflated.H";
paramflist[220] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/rh.sulc";
paramflist[221] = "/cluster/con/9/users/greve/fsaverage_sym-data/136/surf/lrreg/rh.smoothwm.H";
surfregflist[74] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/lh.fssym.i00.reg";
annotflist[74] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/lh.aparc.annot";
paramflist[222] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/lh.inflated.H";
paramflist[223] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/lh.sulc";
paramflist[224] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/lh.smoothwm.H";
surfregflist[75] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/rh.fssym.i00.reg";
annotflist[75] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/rh.aparc.annot";
paramflist[225] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/rh.inflated.H";
paramflist[226] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/rh.sulc";
paramflist[227] = "/cluster/con/9/users/greve/fsaverage_sym-data/138/surf/lrreg/rh.smoothwm.H";
surfregflist[76] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/lh.fssym.i00.reg";
annotflist[76] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/lh.aparc.annot";
paramflist[228] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/lh.inflated.H";
paramflist[229] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/lh.sulc";
paramflist[230] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/lh.smoothwm.H";
surfregflist[77] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/rh.fssym.i00.reg";
annotflist[77] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/rh.aparc.annot";
paramflist[231] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/rh.inflated.H";
paramflist[232] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/rh.sulc";
paramflist[233] = "/cluster/con/9/users/greve/fsaverage_sym-data/141/surf/lrreg/rh.smoothwm.H";
surfregflist[78] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/lh.fssym.i00.reg";
annotflist[78] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/lh.aparc.annot";
paramflist[234] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/lh.inflated.H";
paramflist[235] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/lh.sulc";
paramflist[236] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/lh.smoothwm.H";
surfregflist[79] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/rh.fssym.i00.reg";
annotflist[79] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/rh.aparc.annot";
paramflist[237] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/rh.inflated.H";
paramflist[238] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/rh.sulc";
paramflist[239] = "/cluster/con/9/users/greve/fsaverage_sym-data/149/surf/lrreg/rh.smoothwm.H";


  mrisp_template=MRISmakeTemplateC(80, surfregflist, annotflist, 3, paramflist, 0);
  MRISPwrite(mrisp_template, "lh.junk.new.tiff") ;
  exit(0);


  surfregflist[0] = "/space/tanha/1/users/greve/vrfp-oct11/vrfp-oct11-anat/surf/lh.sphere.reg";
  annotflist[0] = "/space/tanha/1/users/greve/vrfp-oct11/vrfp-oct11-anat/label/lh.aparc.annot";
  paramflist[0] = "/space/tanha/1/users/greve/vrfp-oct11/vrfp-oct11-anat/surf/lh.inflated.H";
  paramflist[1] = "/space/tanha/1/users/greve/vrfp-oct11/vrfp-oct11-anat/surf/lh.sulc";
  paramflist[2] = "/space/tanha/1/users/greve/vrfp-oct11/vrfp-oct11-anat/surf/lh.smoothwm.H";
  mrisp_template=MRISmakeTemplateC(1, surfregflist, annotflist, 3, paramflist, 0);
  MRISPwrite(mrisp_template, "lh.junk.new.tiff") ;

  exit(0);

  subjlist[0] = "fsaverage";
  hemilist[0] = "lh";
  hemilist[1] = "rh";
  mrisp_template=MRISmakeTemplate(1, subjlist, 2, hemilist, "sphere.left_right");
  MRISPwrite(mrisp_template, "lh.fsaverage.left-righ.tiff") ;
  MRISPfree(&mrisp_template) ;

  hemilist[0] = "rh";
  hemilist[1] = "lh";
  mrisp_template=MRISmakeTemplate(1, subjlist, 2, hemilist, "sphere.left_right");
  MRISPwrite(mrisp_template, "rh.fsaverage.left-righ.tiff") ;
  MRISPfree(&mrisp_template) ;

  exit(0);

  hemilist[0] = "lh";
  mrisp_template=MRISmakeTemplateB(1, subjlist, 1, hemilist, "sphere.reg");
  MRISPwrite(mrisp_template, "lh.junk.tiff") ;
  MRISPfree(&mrisp_template) ;
  exit(0);

  hemilist[0] = "rh";
  mrisp_template=MRISmakeTemplate(2, subjlist, 1, hemilist, "sphere.reg");
  MRISPwrite(mrisp_template, "rh.junk.tiff") ;
  MRISPfree(&mrisp_template) ;

  exit(0);

  stats = ComputeBrainVolumeStats(argv[1]);
  for(k=0; k<11; k++) printf("%2d  %g\n",k+1,stats[k]);
  exit(0);

  mri  = MRIread(argv[1]);
  mri2 = MRIread(argv[2]);
  mri3 = MRIfixAsegWithRibbon(mri, mri2,NULL);
  MRIwrite(mri3,"aseg.fixed2.mgz");
  exit(0);

  setenv("SUBJECTS_DIR","/space/freesurfer/subjects/fsfast-tutorial.subjects",1);
  regfile = argv[1];
  regio_read_register(regfile, &subject, &inplaneres, &betplaneres, &intensity,  &R, &float2int);

  sprintf(tmpstr,"%s/%s/surf/lh.white",getenv("SUBJECTS_DIR"),subject);
  printf("%s\n",tmpstr);
  mris = MRISread(tmpstr);

  sprintf(tmpstr,"%s/%s/surf/lh.sphere",getenv("SUBJECTS_DIR"),subject);
  sphere = MRISread(tmpstr);

  mri  = MRIread(argv[2]); // mask
  mri2 = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, 1);
  MRIcopyHeader(mri,mri2);
  mri2 = MRIconst(mri->width, mri->height, mri->depth, 1, 0, mri2);
  mri3 = NULL;

  vox2ras = MRIxfmCRS2XYZtkreg(mri);
  Qsrc = MatrixInverse(vox2ras,NULL);

  nmask = 0;
  nvox = 0;
  fp = fopen("results.dat","w");
  for(c=0; c < mri->width; c++){
    for(r=0; r < mri->height; r++){
      for(s=0; s < mri->depth; s++){
	if(MRIgetVoxVal(mri,c,r,s,0) < 0.5) continue;
	nmask++;
	MRIsetVoxVal(mri2,c,r,s,0, 1.0);
        //mri3 = MRIvol2surfVSM(mri2, R, mris, NULL, SAMPLE_TRILINEAR, NULL, 0, 0, 1, mri3);
        mri3 = MRIvol2surfVSM(mri2, R, mris, NULL, SAMPLE_NEAREST, NULL, 0, 0, 1, mri3);
	nhits = 0;
	for(k=0; k < mris->nvertices; k++){
	  mris->vertices[k].val   = MRIFseq_vox(mri3,k,0,0,0);
	  sphere->vertices[k].val = MRIFseq_vox(mri3,k,0,0,0);
	  if(MRIFseq_vox(mri3,k,0,0,0) > 0.0001) nhits ++;
	}
	if(nhits < 1) {
	  MRIsetVoxVal(mri2,c,r,s,0, 0.0);
	  for(k=0; k < mris->nvertices; k++) MRIsetVoxVal(mri3,k,0,0,0, 0.0); 
	  continue;
	}
	SurfClustList = sclustMapSurfClusters(sphere,10e-10,-1,1,0,&nClusters,NULL);
	if(SurfClustList[0].nmembers < 2){
	  MRIsetVoxVal(mri2,c,r,s,0, 0.0); 
	  for(k=0; k < mris->nvertices; k++) MRIsetVoxVal(mri3,k,0,0,0, 0.0); 
	  continue;
	}
	printf("%5d %2d %2d %2d   %2d %2d\n",nmask,c,r,s,nClusters,SurfClustList[0].nmembers);
	fprintf(fp,"%5d %2d %2d %2d   %2d %2d\n",nmask,c,r,s,nClusters,SurfClustList[0].nmembers);
	if(nClusters > 2){
	  nsup = 0;
	  for(n=0; n<nClusters;n++) if(SurfClustList[n].maxval > .1) nsup++;
	  if(nsup > 1){
	    nvox ++;
	    fprintf(fp,"nvox = %5d   %5d %2d %2d %2d   %2d %2d\n",nvox,nmask,c,r,s,nClusters,SurfClustList[0].nmembers);
	    DumpSurfClusterSum(fp,SurfClustList,nClusters);
	    fflush(fp);
	    printf("nvox = %5d   %5d %2d %2d %2d   %2d %2d\n",nvox,nmask,c,r,s,nClusters,SurfClustList[0].nmembers);
	    DumpSurfClusterSum(stdout,SurfClustList,nClusters);
	    fflush(stdout);
	    //MRIwrite(mri2,"mri2.mgh");
	    //MRIwrite(mri3,"mri3.mgh");
	    //exit(1);
	  }
	}
	// set back to 0
	free(SurfClustList);
	for(k=0; k < mris->nvertices; k++) MRIsetVoxVal(mri3,k,0,0,0, 0.0);
	MRIsetVoxVal(mri2,c,r,s,0, 0.0); 
      }
    }
  }
  printf("nvox = %d\n",nvox);
  fprintf(fp,"nvox = %d\n",nvox);

  exit(0);

  gca = GCAread(argv[1]);
  GCAprint(gca, stdout);

  mri = GCAnlables(gca,NULL);
  MRIwrite(mri,"gca.nlabels.mgh");

  //mri = GCAmri2(gca,NULL);
  //MRIwrite(mri,"gca.mgh");
  exit(1);

  printf("UFSS %s\n",getenv("USE_FAST_SURF_SMOOTHER"));
  mris = MRISread(argv[1]);  

  printf("Init\n");
  for(k=0; k<mris->nvertices; k++){
    mris->vertices[k].dx = drand48();
    mris->vertices[k].dy = drand48();
    mris->vertices[k].dz = drand48();
    mris->vertices[k].tdx = 0;
    mris->vertices[k].tdy = 0;
    mris->vertices[k].tdz = 0;
  }
  mris->vertices[100].ripflag = 1; // Make sure 1 is ripped
  mri = MRIcopyMRIS(NULL, mris, 2, "dz");
  mri = MRIcopyMRIS(mri,  mris, 1, "dy");
  mri = MRIcopyMRIS(mri,  mris, 0, "dx");

  printf("Running\n");
  TimerStart(&mytimer) ;
  MRISaverageGradientsFast(mris, 10);
  msec = TimerStop(&mytimer) ;
  printf("done %6.2f min\n",msec/(1000*60.0));

  printf("dxtx after %g %g\n",mris->vertices[100].dx,mris->vertices[100].tdx);

  mri2 = MRIcopyMRIS(NULL, mris, 2, "tdz");
  mri2 = MRIcopyMRIS(mri2, mris, 1, "tdy");
  mri2 = MRIcopyMRIS(mri2, mris, 0, "tdx");
  MRIwrite(mri2,"new.t.mgh");

  // restore
  MRIScopyMRI(mris, mri, 0, "dx");
  MRIScopyMRI(mris, mri, 1, "dy");
  MRIScopyMRI(mris, mri, 2, "dz");

  printf("dxtx before %g %g\n",mris->vertices[100].dx,mris->vertices[100].tdx);

  printf("Running\n");
  TimerStart(&mytimer) ;
  MRISaverageGradients(mris, 10);
  msec = TimerStop(&mytimer) ;
  printf("done %6.2f min\n",msec/(1000*60.0));

  printf("dxtx after %g %g\n",mris->vertices[100].dx,mris->vertices[100].tdx);

  mri2 = MRIcopyMRIS(mri2, mris, 2, "tdz");
  mri2 = MRIcopyMRIS(mri2, mris, 1, "tdy");
  mri2 = MRIcopyMRIS(mri2, mris, 0, "tdx");
  MRIwrite(mri2,"old.t.mgh");

  exit(0);

  //MRISaverageGradientsFastCheck(10);
  //exit(1);

  //MRISsmoothMRIFastCheck(10);
  //exit(1);

  printf("Reading\n");
  mris = MRISread(argv[1]);  
  mri = MRIread(argv[2]);  
  mask = MRIread(argv[3]);  

  printf("Fast\n");
  mri2 = MRISsmoothMRIFast(mris, mri, 10, mask, NULL);
  //mri2 = MRISsmoothMRIFast(mris, mri, 1, NULL, NULL);
  MRIwrite(mri2,"fast.mgh");

  printf("Slow\n");
  mri2 = MRISsmoothMRI(mris, mri, 10, mask, NULL);
  //mri2 = MRISsmoothMRI(mris, mri, 1, NULL, NULL);
  MRIwrite(mri2,"slow.mgh");

  exit(1);


  mri = MRISmercatorGrid(mris, 10*M_PI/180,10*M_PI/180);
  MRIwrite(mri,argv[2]);

  exit(0);

  region.x = 10;
  region.y = 12;
  region.z = 4;
  region.dx = 40;
  region.dy = 50;
  region.dz = 20;

  seg1 = MRIread(argv[1]);
  seg2 = MRIextractRegionAndPad(seg1, NULL, &region, 4);
  MRIwrite(seg2,"tmp.mgh");

  exit(0);
  subject = argv[1];
  hemi = argv[2];
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
  white = MRISread(tmpstr);
  if(!white) exit(1);
  sprintf(tmpstr,"%s/%s/label/%s.%s.annot",SUBJECTS_DIR,subject,hemi,"aparc");
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(white, tmpstr);
  if(err) exit(1);

  stg = CTABentryNameToAnnotation("superiortemporal", white->ct);
  ttg = CTABentryNameToAnnotation("transversetemporal", white->ct);
  mtg = CTABentryNameToAnnotation("middletemporal", white->ct);
  bsts = CTABentryNameToAnnotation("bankssts", white->ct);
  insula = CTABentryNameToAnnotation("insula", white->ct);

  c = 0;
  for(k=0; k < white->nvertices;k++){
    if(white->vertices[k].annotation == stg ||
       white->vertices[k].annotation == ttg ||
       white->vertices[k].annotation == mtg ||
       white->vertices[k].annotation == insula ||
       white->vertices[k].annotation == bsts){
      c++;
    }
    else{
      white->vertices[k].ripflag = 1;
    }
  }
  printf("n = %d\n",c);
  MRISripFaces(white);

  sprintf(tmpstr,"%s/%s/surf/%s.temporal.patch",SUBJECTS_DIR,subject,hemi);
  MRISwritePatch(white,tmpstr);



  exit(0);

  k = 0;
  while(1){
    printf("%5d\n",k);
    mri = MRIread(argv[1]);
    if(mri == NULL) exit(1);
    k = k + 1;
  }

  mri = MRIread(argv[1]);
  MRIsetVox2RASFromMatrixUnitTest(mri);
  exit(1);


  subject = argv[1];
  hemi = argv[2];
  fwhm = 10;
  fwhmSqr = fwhm*fwhm;

  sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
  white = MRISread(tmpstr);
  if(!white) exit(1);
  printf("White Total Area = %g \n",white->total_area);

  sprintf(tmpstr,"%s/%s/surf/%s.sphere",SUBJECTS_DIR,subject,hemi);
  sphere = MRISread(tmpstr);
  if(!sphere) exit(1);
  printf("Sphere Total Area = %g \n",sphere->total_area);

  f = sqrt(white->total_area/sphere->total_area);
  printf("f = %g \n",f);

  // normalize radius to 1
  vtx1s = &(sphere->vertices[100]);
  for(vtxno1 = 0; vtxno1 < sphere->nvertices; vtxno1++){
    vtx1s = &(sphere->vertices[vtxno1]);
    radius = sqrt(vtx1s->x*vtx1s->x + vtx1s->y*vtx1s->y + vtx1s->z*vtx1s->z);
    vtx1s->x *= f;
    vtx1s->y *= f;
    vtx1s->z *= f;
  }
  MRIScomputeMetricProperties(sphere);
  printf("Sphere Total Area = %g \n",sphere->total_area);

  radius = sqrt(vtx1s->x*vtx1s->x + vtx1s->y*vtx1s->y + vtx1s->z*vtx1s->z);
  printf("Radius %f\n",radius);
  radius2 = radius*radius;

  printf("Alloc\n");
  mri = MRIallocSequence(white->nvertices, 1, 1, MRI_FLOAT, 1);

  printf("loop\n");
  for(vtxno1 = 0; vtxno1 < white->nvertices-1; vtxno1++){
    if(vtxno1%100 ==0) printf("%6d \n",vtxno1);
    vtx1w = &(white->vertices[vtxno1]);
    vtx1s = &(sphere->vertices[vtxno1]);

    for(vtxno2 = vtxno1+1; vtxno2 < white->nvertices; vtxno2++){
      vtx2w = &(white->vertices[vtxno2]);
      vtx2s = &(sphere->vertices[vtxno2]);

      d3Sqr = SQR(vtx1w->x-vtx2w->x)+SQR(vtx1w->y-vtx2w->y)+SQR(vtx1w->z-vtx2w->z);
      if(d3Sqr > fwhmSqr) continue;
      d3 = sqrt(d3Sqr);

      DotProd = (vtx1s->x*vtx2s->x + vtx1s->y*vtx2s->y + vtx1s->z*vtx2s->z)/radius2;
      if(DotProd > +1) DotProd = +1;
      if(DotProd < -1) DotProd = -1;
      theta = acos(DotProd);
      d2 = f * radius * theta;
      if(d2 < fwhm) continue;

      k = MRIgetVoxVal(mri,vtxno1,0,0,0);
      MRIsetVoxVal(mri,vtxno1,0,0,0, k+1);

      //printf("%6d %6d  %4.1f %4.1f\n",vtxno1,vtxno2,d3,d2);

    }
    //if(vtxno1 > 5000 ) break;
  }

  MRIwrite(mri,"count.mgh");

  exit(0);




  k = 1;
  printf("%d %s\n",k,argv[k]);
  seg1 = MRIread(argv[k]);
  if(seg1 == NULL) exit(1);

  k = 2;
  printf("%d %s\n",k,argv[k]);
  seg2 = MRIread(argv[k]);
  if(seg2 == NULL) exit(1);

  mri = MRIadd(seg1,seg2,NULL);
  if(mri == NULL) exit(1);

  for(k=3; k < argc-1; k++){
    printf("%d %s\n",k,argv[k]);
    seg2 = MRIread(argv[k]);
    if(seg2 == NULL) exit(1);
    mri = MRIadd(mri,seg2,mri);
    if(mri == NULL) exit(1);
  }

  f = 1.0/((double)(argc-2));
  printf("Dividing by %lf\n",f);
  MRImultiplyConst(mri, f, mri);
  printf("Saving to %s\n",argv[argc-1]);
  MRIwrite(mri,argv[argc-1]);
  exit(0);


  mri = MRIsegDiff(seg1,seg2,&k);
  printf("DiffFlag %d\n",k);
  MRIwrite(mri,"aseg.diff.mgz");

  mri2 = MRIsegMergeDiff(seg1, mri);
  MRIwrite(mri2,"aseg.new.mgz");

  exit(0);

  seg1 = MRIread(argv[1]);
  seg2 = MRIread(argv[2]);
  dice = MRIsegDice(seg1, seg2, &nsegs, &segidlist);

  ctab = CTABreadASCII("/space/greve/1/users/greve/freesurfer/FreeSurferColorLUT.txt");
  if (ctab == NULL) {
    printf("ERROR: reading ctab\n");
    exit(1);
  }

  for(k=0; k < nsegs; k++){
    if(segidlist[k] >= 0){
      CTABcopyName(ctab,segidlist[k],tmpstr,sizeof(tmpstr));
      printf("%2d %4d %-30s %7.5lf\n",k,segidlist[k],tmpstr,dice[k]);
    }
  }


  return(0);


  //----------------------------------------------------
  subject = argv[1];
  hemi = argv[2];
  annot1 = argv[3];
  annot2 = argv[4];

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) exit(1);
  MRIScomputeMetricProperties(surf);

  surf2 = MRISread(tmpstr);
  if(!surf2) exit(1);

  sprintf(tmpstr,"%s/%s/label/%s.%s.annot",SUBJECTS_DIR,subject,hemi,annot1);
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(surf, tmpstr);
  if(err) exit(1);

  sprintf(tmpstr,"%s/%s/label/%s.%s.annot",SUBJECTS_DIR,subject,hemi,annot2);
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(surf2, tmpstr);
  if(err) exit(1);

  dice = MRISannotDice(surf, surf2, &nsegs, &segidlist);

  for(k=0; k < nsegs; k++){
    if(segidlist[k] >= 0){
      printf("%2d %4d %-25s  %7.5lf\n",k,segidlist[k],
	     surf->ct->entries[segidlist[k]]->name,
	     dice[k]);
    }
  }



  return(0);
  //----------------------------------------------------

  seg1 = MRISannot2seg(surf,1000);
  seg2 = MRISannot2seg(surf2,1000);

  printf("Generating list of segmentation ids\n");
  segidlist1 = MRIsegIdList(seg1, &nsegid1,0);
  printf("Found %d\n",nsegid1);
  fflush(stdout);

  printf("Generating list of segmentation ids\n");
  segidlist2 = MRIsegIdList(seg1, &nsegid2,0);
  printf("Found %d\n",nsegid2);
  fflush(stdout);

  area1 = (double *) calloc(nsegid1,sizeof(double));
  area2 = (double *) calloc(nsegid1,sizeof(double));
  area12 = (double *) calloc(nsegid1,sizeof(double));

  for(c=0; c < seg1->width; c++){
    for(r=0; r < seg1->height; r++){
      for(s=0; s < seg1->depth; s++){
	a1 = MRIgetVoxVal(seg1,c,r,s,0);
	k1 = -1;
	for(k=0; k < nsegid1; k++) {
	  if(a1 == segidlist1[k]){
	    k1 = k; 
	    break;
	  }
	}
	a2 = MRIgetVoxVal(seg2,c,r,s,0);
	k2 = -1;
	for(k=0; k < nsegid2; k++) {
	  if(a2 == segidlist2[k]){
	    k2 = k; 
	    break;
	  }
	}
	area1[k1] += surf->vertices[c].area;
	area2[k2] += surf->vertices[c].area;
	if(a1 == a2) area12[k1] += surf->vertices[c].area;
      }
    }
  }

  for(k=0; k < nsegid1; k++) {
    printf("%2d %4d   %7.1lf %7.1lf %7.1lf   %6.4f\n",
	   k,segidlist1[k],
	   area1[k],area2[k],area12[k],
	   (float)area12[k]/((area1[k]+area2[k])/2.0));
  }


  exit(1);

  sprintf(tmpstr,"%s/%s/label/%s.aparc.annot",SUBJECTS_DIR,subject,hemi);
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(surf, tmpstr);
  if(err) exit(1);


  MRISfbirnAnnot(surf);
  sprintf(tmpstr,"%s.fbirn.annot",hemi);
  MRISwriteAnnotation(surf, tmpstr);

  exit(1);
  //----------------------------------------------

  mri = MRIread(argv[1]);
  if(mri == NULL) exit(1);

  nvox = 0;
  sumval = 0;
  for (c=0; c < mri->width; c++) {
    for (r=0; r < mri->height; r++) {
      for (s=0; s < mri->depth; s++) {
	sumval += MRIgetVoxVal(mri,c,r,s,0);
	nvox ++;
      }
    }
  }
  printf("%lf\n",sumval/nvox);
  exit(0);

  if(argc <= 1){
    printf("dngtester subject hemi\n");
    exit(1);
  }
  subject = argv[1];
  hemi = argv[2];
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/surf/%s.sphere",SUBJECTS_DIR,subject,hemi);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) exit(1);

  sprintf(tmpstr,"%s/%s/surf/%s.occip.patch.flat",SUBJECTS_DIR,subject,hemi);
  MRISreadPatchNoRemove(surf, tmpstr) ;

  RETlogMap(surf, 15, .7, 0, 0);
  mri = MRIcopyMRIS(NULL, surf, 0, "val");
  MRIwrite(mri,"logmap.mgh");

  exit(1);

  //-----------------------------------------------------
  printf("nsurfs %d\n",argc-1);
  for(k=1; k<argc; k++){
    printf("Loading %s\n",argv[k]);
    surf = MRISread(argv[k]);
    if(surf == NULL) exit(1);
    surfs[k-1] = surf;
  }
  surf = MRISaverageSurfaces(argc-1, surfs);
  MRISwrite(surf,"lh.avgsurf");

  exit(1);

  vtx = &(surf->vertices[12282]);
  dlhw = sqrt((vtx->x * vtx->x) + (vtx->y * vtx->y) + (vtx->z * vtx->z));
  DotProdThresh = (100*100)*cos(10/dlhw)*(1.0001);

  XNbrVtxNo   = (int *) calloc(surf->nvertices,sizeof(int));
  XNbrDotProd = (double *) calloc(surf->nvertices,sizeof(double));
  mri = MRIallocSequence(surf->nvertices, 1, 1, MRI_FLOAT, 1);

  //-----------------------------------------------
  nXNbrs = 0;
  MRISextendedNeighbors(surf, 12282, 12282, DotProdThresh,
			XNbrVtxNo, XNbrDotProd, 
			&nXNbrs,surf->nvertices,1);
  printf("Found %d neighbors\n",nXNbrs);

  for(k = 0; k < nXNbrs; k++)
    MRIsetVoxVal(mri,XNbrVtxNo[k],0,0,0,1);
  MRIwrite(mri,"xnbr.mgh");

  //-----------------------------------------------
  err = MRISreadPatch(surf, argv[3]);

  for(lhwvtx = 0; lhwvtx < surf->nvertices; lhwvtx++){
    if(surf->vertices[lhwvtx].ripflag) continue;
    for(k = 0; k < surf->nvertices; k++) surf->vertices[k].val2bak = -1;
    nXNbrs = 0;
    MRISextendedNeighbors(surf, lhwvtx, lhwvtx, 5*5,
			  XNbrVtxNo, XNbrDotProd, 
			  &nXNbrs,surf->nvertices,2);
    if(lhwvtx%1000 == 1) printf("%5d %5d neighbors\n",lhwvtx,nXNbrs);
  }

  for(k = 0; k < nXNbrs; k++)
    MRIsetVoxVal(mri,XNbrVtxNo[k],0,0,0,1);

  MRIwrite(mri,"xnbr2.mgh");

  exit (1);

  //----------------------------------------------------
  //R = MatrixDRand48(20,20,NULL);
  R = MatrixIdentity(100,NULL);
  MatrixPrint(stdout,R);
  I = ImageFromMatrix(R, NULL);
  ImageWrite(I, "./myimage.jpg") ;
  return(1);


  //----------------------------------------------------
  subject = argv[1];
  hemi = argv[2];
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/surf/%s.white",SUBJECTS_DIR,subject,hemi);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) exit(1);

  sprintf(tmpstr,"%s/%s/mri/aseg.mgz",SUBJECTS_DIR,subject);
  mri = MRIread(tmpstr);
  if(!mri) exit(1);

  //#define MIN_NONCORTEX_VERTICES 10
  //lcortex = MRIScortexLabel(surf, mri, MIN_NONCORTEX_VERTICES) ;
  //sprintf(tmpstr,"%s/%s/label/%s.%s.label",SUBJECTS_DIR, subject,hemi,"cortex");
  //printf("writing cortex label to %s...\n", tmpstr) ;
  //LabelWrite(lcortex, tmpstr) ;

  sprintf(tmpstr,"%s/%s/label/%s.aparc.annot",SUBJECTS_DIR,subject,hemi);
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(surf, tmpstr);
  if(err) exit(1);

  parcnames[0] = "caudalmiddlefrontal";
  parcnames[1] = "superiorfrontal";
  parcnames[2] = "rostralmiddlefrontal";
  parcnames[3] = "parsopercularis";
  parcnames[4] = "parstriangularis";
  parcnames[5] = "parsorbitalis";
  parcnames[6] = "lateralorbitofrontal";
  parcnames[7] = "medialorbitofrontal";
  parcnames[8] = "paracentral";
  parcnames[9] = "frontalpole";
  MRISmergeAnnotations(surf, 10, parcnames, "frontal");

  parcnames[0] = "superiortemporal";
  parcnames[1] = "entorhinal";
  parcnames[2] = "temporalpole";
  parcnames[3] = "fusiform";
  parcnames[4] = "inferiortemporal";
  parcnames[5] = "middletemporal";
  parcnames[6] = "parahippocampal";
  parcnames[7] = "bankssts";
  parcnames[8] = "transversetemporal";
  MRISmergeAnnotations(surf, 9, parcnames, "temporal");

  parcnames[0] = "supramarginal";
  parcnames[1] = "inferiorparietal";
  parcnames[2] = "superiorparietal";
  parcnames[3] = "precuneus";
  MRISmergeAnnotations(surf, 4, parcnames, "parietal");

  parcnames[0] = "pericalcarine";
  parcnames[1] = "cuneus";
  parcnames[2] = "lingual";
  parcnames[3] = "lateraloccipital";
  MRISmergeAnnotations(surf, 4, parcnames, "occipital");

  parcnames[0] = "isthmuscingulate";  
  parcnames[1] = "posteriorcingulate";
  parcnames[2] = "caudalanteriorcingulate";
  parcnames[3] = "rostralanteriorcingulate";
  MRISmergeAnnotations(surf, 4, parcnames, "cingulate");

  sprintf(tmpstr,"%s.lobes.annot",hemi);
  MRISwriteAnnotation(surf, tmpstr);

  //------------------------------------------------------------
  ct = CTABaddEntry(surf->ct,"area32p");
  CTABfree(&surf->ct);
  surf->ct = ct;

  ct = CTABaddEntry(surf->ct,"area32v");
  CTABfree(&surf->ct);
  surf->ct = ct;

  CTABfindName(surf->ct, "area32p", &area32p);
  CTABfindName(surf->ct, "area32v", &area32v);
  CTABfindName(surf->ct, "superiorfrontal", &superiorfrontal);
  CTABfindName(surf->ct, "medialorbitofrontal", &medialorbitofrontal);
  CTABfindName(surf->ct, "rostralanteriorcingulate", &rostralanteriorcingulate);
  CTABfindName(surf->ct, "rostralmiddlefrontal", &rostralmiddlefrontal);

  mri = MRISfbirnMask_MOF_RACing(surf);
  MRISdilateConfined(surf, mri, medialorbitofrontal, 12, area32v);

  mri = MRISfbirnMask_SFG_Cing(surf);
  MRISdilateConfined(surf, mri, superiorfrontal, 12, area32p);


  nunits = (int *)calloc(surf->ct->nentries, sizeof(int)) ;
  nunits[rostralanteriorcingulate] = 3;
  nunits[rostralmiddlefrontal] = 3;
  nunits[superiorfrontal] = 5;
  nunits[area32p] = 2;
  //nunits[area32v] = 2;

  MRISdivideAnnotation(surf, nunits) ;

  CTABfindName(surf->ct, "area32p", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area32p_pseudo");

  CTABfindName(surf->ct, "area32p_div2", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area32a_pseudo");

  CTABfindName(surf->ct, "area32v", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area32v_pseudo");

  //CTABfindName(surf->ct, "area32v_div2", &index);
  //sprintf(surf->ct->entries[index]->name, "%s","area32v_pseudo");


  CTABfindName(surf->ct, "rostralanteriorcingulate", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area24d_pseudo");

  CTABfindName(surf->ct, "rostralanteriorcingulate_div2", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area24pg_pseudo");

  CTABfindName(surf->ct, "rostralanteriorcingulate_div3", &index);
  sprintf(surf->ct->entries[index]->name, "%s","area24v_pseudo");

  sprintf(tmpstr,"%s.fbirn0.annot",hemi);
  MRISwriteAnnotation(surf, tmpstr);

  return(0);

  printrgb();
  return(0);

  printf("Reading\n");
  mri = MRIread(argv[1]);
  printf("chunck %d\n",mri->ischunked);

  printf("Done\n");

  exit(0);


  mri2 = MRIvote(mri, mask, NULL);
  MRIwrite(mri2,"vote.mgh");

  return(0);
  exit(0);
}

/*--------------------------------------------------------*/
MRI *fMRIvariance2(MRI *fmri, float DOF, int RmMean, MRI *var) {
  int c, r, s, f;
  double val,sumsqval, sumval;
  int nvox_per_row, nvox_per_slice, bytes_per_vol;
  void *p;
  val = 0;

  if (DOF < 0) DOF = fmri->nframes;

  if (var==NULL) {
    var = MRIallocSequence(fmri->width, fmri->height, fmri->depth,
                           MRI_FLOAT, 1);
    if (var==NULL) {
      printf("ERROR: fMRIvariance: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(fmri,var);
  } else {
    if (var->width  != fmri->width ||
        var->height != fmri->height ||
        var->depth  != fmri->depth) {
      printf("ERROR: fMRIvariance: output dimension mismatch\n");
      return(NULL);
    }
    if (var->type != MRI_FLOAT) {
      printf("ERROR: fMRIvariance: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  nvox_per_row = fmri->width;
  nvox_per_slice = fmri->width * fmri->height;
  bytes_per_vol = fmri->width * fmri->height * fmri->depth * fmri->bytes_per_vox;
  for (c=0; c < fmri->width; c++) {
    for (r=0; r < fmri->height; r++) {
      for (s=0; s < fmri->depth; s++) {
        sumval = 0;
        sumsqval = 0;
        p = fmri->chunk + (c + r*nvox_per_row + s*nvox_per_slice)*fmri->bytes_per_vox;
        for (f=0; f < fmri->nframes; f++) {
          if (fmri->ischunked) {
            switch (fmri->type) {
            case MRI_UCHAR:
              val = (double)(*((char *)p));
              break;
            case MRI_SHORT:
              val = (double)(*((short*)p));
              break;
            case MRI_INT:
              val = (double)(*((int  *)p));
              break;
            case MRI_LONG:
              val = (double)(*((long *)p));
              break;
            case MRI_FLOAT:
              val = (double)(*((float*)p));
              break;
            }
          } else val = MRIgetVoxVal(fmri, c, r, s, f);
          //printf("%d  %lu   %g  %g\n",f,(long int)p,val,MRIgetVoxVal(fmri,c,r,s,f));
          sumsqval += (val*val);
          if (RmMean) sumval += val;
          p += bytes_per_vol;
        }
        MRIFseq_vox(var,c,r,s,0) = sumsqval/DOF;
        if (RmMean)
          MRIFseq_vox(var,c,r,s,0) -= ((sumval/DOF)*(sumval/DOF));
      }
    }
  }

  return(var);
}



/*---------------------------------------------------------------*/
/*
LTA *TransformRegDat2LTA(MRI *targ, MRI *mov, MATRIX *R)
{
  LTA *lta;
  MATRIX *vox2vox; // Targ->Mov
  MATRIX *Ttarg, *Tmov, *invTmov;

  Ttarg = MRIxfmCRS2XYZtkreg(targ);
  Tmov  = MRIxfmCRS2XYZtkreg(mov);
  invTmov = MatrixInverse(Tmov,NULL);

  // vox2vox = invTmov * R * Ttarg
  vox2vox = MatrixMultiply(invTmov,R,NULL);
  MatrixMultiply(vox2vox,Ttarg,vox2vox);

  lta = LTAalloc(1,NULL);
  lta->type = LINEAR_VOX_TO_VOX;
  lta->xforms[0].type = LINEAR_VOX_TO_VOX;
  getVolGeom(targ,&lta->xforms[0].src);
  getVolGeom(mov,&lta->xforms[0].dst);
  lta->xforms[0].m_L = MatrixCopy(vox2vox,NULL);

  MatrixFree(&Ttarg);
  MatrixFree(&Tmov);
  MatrixFree(&invTmov);
  MatrixFree(&vox2vox);

  return(lta);
}
*/

/*---------------------------------------------------------------*/
MRI *MRIsetSliceNo(MRI *mri, MRI *out) {
  int c, r, s, f, n, ncols, nrows, nslices,nframes;
  void   *pmri=NULL, *pout=NULL;
  int sz, szout;

  ncols   = mri->width;
  nrows   = mri->height;
  nslices = mri->depth;
  nframes = mri->nframes;

  if (out==NULL) {
    out = MRIallocSequence(ncols, nrows, nslices, mri->type, nframes);
    if (out==NULL) {
      printf("ERROR: MRIsetSliceNo: could not alloc output\n");
      return(NULL);
    }
    MRIcopyHeader(mri,out); // ordinarily would need to change nframes
  }
  if (out->width != ncols   || out->height != nrows ||
      out->depth != nslices || out->nframes != nframes) {
    printf("ERROR: MRIsetSliceNo: dimension mismatch\n");
    return(NULL);
  }

  // Number of bytes in the mri data types
  sz   = MRIsizeof(mri->type);
  szout = MRIsizeof(out->type);

  n = 0;
  for (f=0; f<nframes; f++) {
    for (s=0; s<nslices; s++) {
      for (r=0; r<nrows; r++) {
        // Pointers to the start of the column
        pmri  = (void *) mri->slices[n][r];
        pout  = (void *) out->slices[n][r];
        for (c=0; c<ncols; c++) {
          MRIdbl2ptr(s, pout, out->type);
          pmri += sz;
          pout  += szout;
        } // cols
      } // rows
      n++;
    } // slices
  } // frames

  return(out);
}

/*-------------------------------------------------------*/
//  if(1){
//  MRISmercator(surf);
//  sprintf(tmpstr,"%s/%s/surf/%s.mercator",SUBJECTS_DIR,subject,hemi);
//  MRISwrite(surf,tmpstr);
//  exit(0);
//  }
int MRISmercator(MRIS *surf)
{
  int vtxno;
  VERTEX *vtx;
  double x,y,z;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  double dthresh;
  dthresh = 5;

  xmin=ymin=zmin =  1e10;
  xmax=ymax=zmax = -1e10;
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    vtx = &(surf->vertices[vtxno]);
    x = vtx->x;
    y = vtx->y;
    z = vtx->z;
    vtx->x = 0;
    vtx->y = (100*(atan2(y,-x)))/M_PI;
    // z stays the same

    x = vtx->x;
    y = vtx->y;
    z = vtx->z;
    if(xmin > x) xmin = x;
    if(xmax < x) xmax = x;
    if(ymin > y) ymin = y;
    if(ymax < y) ymax = y;
    if(zmin > z) zmin = z;
    if(zmax < z) zmax = z;
  }

  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    vtx = &(surf->vertices[vtxno]);
    if(0){
      if(fabs(vtx->x - xmax) < dthresh) vtx->y = -100;
      if(fabs(vtx->x - xmin) < dthresh) vtx->y = -100;
      if(fabs(vtx->z - zmax) < dthresh) vtx->y = -100;
      if(fabs(vtx->z - zmin) < dthresh) vtx->y = -100;
    } else {
      if(fabs(vtx->y - ymax) < dthresh) vtx->x = -100;
      if(fabs(vtx->y - ymin) < dthresh) vtx->x = -100;
      if(fabs(vtx->z - zmax) < dthresh) vtx->x = -100;
      if(fabs(vtx->z - zmin) < dthresh) vtx->x = -100;
    }
  }

  return(0);
}

MRIS *MRISaverageSurfaces(int nsurfaces, MRIS **surfs)
{
  MRIS *surf, *avgsurf;
  int n,k;
  float average_surface_area;

  surf = surfs[0];
  avgsurf = MRISalloc(surf->nvertices, surf->nfaces) ;

  // Make sure xyz is 0, copy faces and vertices
  for(k=0; k < surf->nvertices; k++){
    avgsurf->vertices[k].x = 0.0;
    avgsurf->vertices[k].y = 0.0;
    avgsurf->vertices[k].z = 0.0;
    avgsurf->vertices[k].num = surf->vertices[k].num;
    avgsurf->vertices[k].f = (int*) calloc(surf->vertices[k].num,sizeof(int));
    avgsurf->vertices[k].n = (uchar*) calloc(surf->vertices[k].num,sizeof(uchar));
    for(n=0; n < surf->vertices[k].num; n++){
      avgsurf->vertices[k].f[n] = surf->vertices[k].f[n];
      avgsurf->vertices[k].n[n] = surf->vertices[k].n[n];
    }
    avgsurf->vertices[k].vnum = surf->vertices[k].vnum;
    avgsurf->vertices[k].v = (int*) calloc(surf->vertices[k].vnum,sizeof(int));
    avgsurf->vertices[k].dist = (float*) calloc(surf->vertices[k].vnum,sizeof(float));
    for(n=0; n < surf->vertices[k].vnum; n++)
      avgsurf->vertices[k].v[n] = surf->vertices[k].v[n];
  }
  for(k=0; k < surf->nfaces; k++) {
    for(n=0; n < VERTICES_PER_FACE; n++)
      avgsurf->faces[k].v[n] = surf->faces[k].v[n];
  }

  // Loop thru all surfaces, sume xyz 
  for(n=0; n < nsurfaces; n++){
    printf("%2d  surface area %g\n",n,surf->total_area);
    surf = surfs[0];
    if(surf->nvertices != surfs[0]->nvertices){
      printf("ERROR: MRISaverageSurfaces(): dimension mismatch surface %d\n",n);
      printf(" number of vertices %d vs %d\n",surf->nvertices,surfs[0]->nvertices);
      return(NULL);
    }
    for(k=0; k < surf->nvertices; k++){
      avgsurf->vertices[k].x += surf->vertices[k].x;
      avgsurf->vertices[k].y += surf->vertices[k].y;
      avgsurf->vertices[k].z += surf->vertices[k].z;
    }
    average_surface_area += surf->total_area ;
  }

  average_surface_area /= nsurfaces;
  printf("average surface area %g\n",surf->total_area);

  // Now divide by number of surfaces
  for(k=0; k < surf->nvertices; k++){
    avgsurf->vertices[k].x /= nsurfaces;
    avgsurf->vertices[k].y /= nsurfaces;
    avgsurf->vertices[k].z /= nsurfaces;
  }

  MRIScomputeMetricProperties(avgsurf);
  printf("avg  surface area %g\n",avgsurf->total_area);
  return(avgsurf);
}
/*!
  \fn MATRIX *MatrixLoadFSL(char *fname)
  \brief Loads in an FSL matrix. Not a registration matrix but one
  stored in their design.con or design.mat files used in the GLM 
  estimation. The format is pretty simple, it just has a /Matrix 
  keyword followed by the data.
*/
MATRIX *MatrixLoadFSL(char *fname)
{
  FILE *fp;
  MATRIX *M0,*M;
  char line[2000], *s, tag[2000];
  int hit,Nc,Nc0,c,nrows,r;

  fp = fopen(fname,"r");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",fname);
    return(NULL);
  }

  // Get past the "/Matrix" string
  hit = 0;
  while(1){
    s = fgets(line, 1999, fp);
    if(s == NULL){
      printf("ERROR: %s is not formatted correctly\n",fname);
      fclose(fp);
      return(NULL);
    }
    sscanf(line,"%s",tag);
    if(! strcmp(tag,"/Matrix")){
      hit = 1;
      break;
    }
  }
  if(!hit){
    printf("ERROR: %s is not formatted correctly\n",fname);
    fclose(fp);
    return(NULL);
  }

  // Now read in each line
  M = NULL;
  Nc0 = 0; Nc=0;
  nrows = 0;
  while(1){
    if(nrows > 5000){
      printf("ERROR: %s, nrows exceeds 5000\n",fname);
      fclose(fp);
      MatrixFree(&M0);
      return(NULL);
    }
    // read in line
    s = fgets(line, 1999, fp);
    if(s == NULL) break;
    // Count number of columns
    Nc = gdfCountItemsInString(line);
    if(Nc == 0){
      printf("ERROR: %s is not formatted correctly\n",fname);
      fclose(fp);
      return(NULL);
    }
    // If first pass, alloc matrix, etc
    if(nrows == 0){
      Nc0 = Nc;
      //printf("%d cols\n",Nc);
      M0 = MatrixAlloc(5000,Nc,MATRIX_REAL);
      if(M0==NULL){
	printf("ERROR: could not alloc %d cols\n",Nc);
	fclose(fp);
	return(NULL);
      }
    }
    // Make sure this row is conistent with previous rows
    if(Nc0 != Nc){
      printf("ERROR: %s is not formatted correctly\n",fname);
      fclose(fp);
      MatrixFree(&M0);
      return(NULL);
    }
    // Read in each colum for this row
    for(c=0; c < Nc0; c++){
      s = gdfGetNthItemFromString(line, c);
      sscanf(s,"%f",&M0->rptr[nrows+1][c+1]);
      free(s);
    }
    nrows ++;
  }
  fclose(fp);
  //printf("%d rows\n",nrows);

  if(nrows == 0){
    printf("ERROR: %s is not formatted correctly\n",fname);
    return(NULL);
  }

  // Pack data into a matrix of the correct size
  M = MatrixAlloc(nrows,Nc,MATRIX_REAL);
  for(r=0; r < nrows; r++)
    for(c=0; c < Nc0; c++)
      M->rptr[r+1][c+1] = M0->rptr[r+1][c+1];
  MatrixFree(&M0);

  MatrixPrint(stdout,M);
  
  return(M);
}
/* ---------------------------------------------------------------------- */

int find_path ( int* vert_vno, int num_vno, int max_path_length,
                int* path, int* path_length, MRIS *mris ) {
  int cur_vert_vno;
  int src_vno;
  int dest_vno;
  int vno;
  char* check;
  float* dist;
  int* pred;
  char done;
  VERTEX* v;
  VERTEX* u;
  float closest_dist;
  int closest_vno;
  int neighbor;
  int neighbor_vno;
  float dist_uv;
  int path_vno;
  int num_path = 0;
  int num_checked;
  float vu_x, vu_y, vu_z;
  int flag2d = 0; // for flattend surface?

  dist = (float*) calloc (mris->nvertices, sizeof(float));
  pred = (int*) calloc (mris->nvertices, sizeof(int));
  check = (char*) calloc (mris->nvertices, sizeof(char));
  num_path = 0;
  num_checked = 0;
  (*path_length) = 0;

  for (cur_vert_vno = 0; cur_vert_vno < num_vno-1; cur_vert_vno++) {
    /* clear everything */
    for (vno = 0; vno < mris->nvertices; vno++) {
      dist[vno] = 999999;
      pred[vno] = -1;
      check[vno] = FALSE;
    }

    /* Set src and dest */
    src_vno = vert_vno[cur_vert_vno+1];
    dest_vno = vert_vno[cur_vert_vno];

    /* make sure both are in range. */
    if (src_vno < 0 || src_vno >= mris->nvertices ||
        dest_vno < 0 || dest_vno >= mris->nvertices)
      continue;

    if (src_vno == dest_vno)
      continue;

    /* pull the src vertex in. */
    dist[src_vno] = 0;
    pred[src_vno] = vno;
    check[src_vno] = TRUE;

    done = FALSE;
    while (!done) {

      /* find the vertex with the shortest edge. */
      closest_dist = 999999;
      closest_vno = -1;
      for (vno = 0; vno < mris->nvertices; vno++)
        if (check[vno])
          if (dist[vno] < closest_dist) {
            closest_dist = dist[vno];
            closest_vno = vno;
          }
      v = &(mris->vertices[closest_vno]);
      check[closest_vno] = FALSE;

      /* if this is the dest node, we're done. */
      if (closest_vno == dest_vno) {
        done = TRUE;
      } else {
        /* relax its neighbors. */
        for (neighbor = 0; neighbor < v->vnum; neighbor++) {
          neighbor_vno = v->v[neighbor];
          u = &(mris->vertices[neighbor_vno]);

          /* calc the vector from u to v. */
          vu_x = u->x - v->x;
          vu_y = u->y - v->y;
	  if (flag2d)	    vu_z = 0;
	  else     	    vu_z = u->z - v->z;

          /* recalc the weight. */
	  if (flag2d)
	    dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
			   ((v->y - u->y) * (v->y - u->y)));
	  else
	    dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
			   ((v->y - u->y) * (v->y - u->y)) +
			   ((v->z - u->z) * (v->z - u->z)));

          /* if this is a new shortest path, update the predecessor,
             weight, and add it to the list of ones to check next. */
          if (dist_uv + dist[closest_vno] < dist[neighbor_vno]) {
            pred[neighbor_vno] = closest_vno;
            dist[neighbor_vno] = dist_uv + dist[closest_vno];
            check[neighbor_vno] = TRUE;
          }
        }
      }
      num_checked++;
      if ((num_checked % 100) == 0) {
        printf (".");
        fflush (stdout);
      }
    }

    /* add the predecessors from the dest to the src to the path. */
    path_vno = dest_vno;
    path[(*path_length)++] = dest_vno;
    while (pred[path_vno] != src_vno &&
           (*path_length) < max_path_length ) {
      path[(*path_length)++] = pred[path_vno];
      path_vno = pred[path_vno];
    }
  }
  printf (" done\n");
  fflush (stdout);

  free (dist);
  free (pred);
  free (check);

  return (ERROR_NONE);
}

int mygd(char *fname)
{
  FILE *fp;
  char *match = "DiffusionGradientDirection\0";
  char *check;
  int matchlen, n, m;
  int c, hit;
  float f;
  char buf[4];
  
  fp = fopen(fname,"r");
  
  matchlen = strlen(match);
  check = calloc(sizeof(char),matchlen);

  n = -1;
  while(fp){
    n++;
    c = fgetc(fp);
    if(n < matchlen){
      check[n] = c;
      continue;
    }
    for(m=0; m < matchlen-1; m++) check[m] = check[m+1];
    check[m] = c;

    hit = 1;
    for(m=0; m < matchlen; m++){
      if(check[m] != match[m]) hit = 0;
    }
    if(hit){
      for(m = 0; m < 10; m++){
	fread(buf,sizeof(float),1,fp);
	byteswapbuffloat(buf,1);
	f = (float)buf[0];
	//printf("%s %7.4f %7.4f %7.4f \n",fname,f[0],f[1],f[2]);
	printf("%d %7.4f\n",m,f);
      }
      return(0);
    }
  }
  return(1);
}

MRI *MRIaddB(MRI *mri1, MRI *mri2, MRI *mriadd)
{
  int c,r,s,f;
  double v1,v2,vadd;

  if(mriadd == NULL){
    mriadd = MRIallocSequence(mri1->width, mri1->height,mri1->depth,
			      MRI_FLOAT, mri1->nframes);
    MRIcopyHeader(mri1, mriadd);
  }

  for(c=0; c < mri1->width; c++){
    for(r=0; r < mri1->height; r++){
      for(s=0; s < mri1->depth; s++){
	for(f=0; f < mri1->nframes; f++){
	  v1 = MRIgetVoxVal(mri1,c,r,s,f);
	  v2 = MRIgetVoxVal(mri2,c,r,s,f);
	  vadd = v1+v2;
	  MRIsetVoxVal(mriadd,c,r,s,f,vadd);
	}
      }
    }
  }
  return(mriadd);
}

MRI *MRInnFill(MRI *src, MRI *seed, MRI *out, MRI *mask)
{
  int c,r,s, cseed, rseed, sseed;
  double d2, d2min, x, y, z;
  int nseeds,nthseed,nthseedmin,nmask;
  double *xseed,*yseed,*zseed, val;

  if(out == NULL){
    out = MRIallocSequence(src->width, src->height, src->depth,MRI_FLOAT, 1);
    if (out==NULL){
      printf("ERROR: MRInnFill(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(src,out);
  }

  printf("Looking for seeds\n");
  nseeds = 0;
  nmask = 0;
  for(c=0; c < seed->width; c++){
    for(r=0; r < seed->height; r++){
      for(s=0; s < seed->depth; s++){
	if(mask != NULL && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	nmask++;
	if(MRIgetVoxVal(seed,c,r,s,0) > 0.5) nseeds++;
      }
    }
  }
  printf("Found %d seeds (%d mask), computing xyz\n",nseeds,nmask);
  xseed = (double *)calloc(sizeof(double),nseeds);
  yseed = (double *)calloc(sizeof(double),nseeds);
  zseed = (double *)calloc(sizeof(double),nseeds);
  nseeds = 0;
  for(c=0; c < seed->width; c++){
    for(r=0; r < seed->height; r++){
      for(s=0; s < seed->depth; s++){
	if(mask != NULL && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	if(MRIgetVoxVal(seed,c,r,s,0) > 0.5){
	  xseed[nseeds] = c*seed->xsize;
	  yseed[nseeds] = r*seed->ysize;
	  zseed[nseeds] = s*seed->zsize;
	  nseeds++;
	}
      }
    }
  }

  printf("Filling \n");
  for(c=0; c < seed->width; c++){
    printf("%d ",c);
    fflush(stdout);
    for(r=0; r < seed->height; r++){
      for(s=0; s < seed->depth; s++){
	if(mask != NULL && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	if(MRIgetVoxVal(seed,c,r,s,0) > 0.5) {
	  val = MRIgetVoxVal(src,c,r,s,0);
	  MRIsetVoxVal(out,c,r,s,0,val);
	  continue;
	}
	x = c*seed->xsize;
	y = r*seed->ysize;
	z = s*seed->zsize;
	d2min = 10e10;
	nthseedmin = 0;
	for(nthseed = 0; nthseed < nseeds; nthseed++){
	  d2 = (x-xseed[nthseed])*(x-xseed[nthseed]) + 
	       (y-yseed[nthseed])*(y-yseed[nthseed]) + 
  	       (z-zseed[nthseed])*(z-zseed[nthseed]);
	  if(d2 < d2min){
	    d2min = d2;
	    nthseedmin = nthseed;
	  }
	}
	cseed = round(xseed[nthseedmin]/seed->xsize);
	rseed = round(yseed[nthseedmin]/seed->ysize);
	sseed = round(zseed[nthseedmin]/seed->zsize);
	val = MRIgetVoxVal(src,cseed,rseed,sseed,0);
	MRIsetVoxVal(out,c,r,s,0,val);
      }
    }
  }
  printf("\n");
  printf("Done filling \n");

  return(out);
}


int mrisComputeDistanceTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR  *v_y, *v_delta, *v_n ;
  float   l_dist=1, d0, d0b, dt, delta, nc, scale, norm;
  VERTEX  *v, *vn ;
  int     vno, n, vnum;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_y = VectorAlloc(3, MATRIX_REAL) ;
  v_delta = VectorAlloc(3, MATRIX_REAL) ;
  norm = 1.0f / mris->avg_nbrs ;

  if (mris->patch)
    scale = 1.0f ;
  else
    if (mris->status == MRIS_PARAMETERIZED_SPHERE || mris->status == MRIS_SPHERE)
      scale = sqrt(mris->orig_area / mris->total_area) ;
    else
      scale = mris->neg_area < mris->total_area ?
              sqrt(mris->orig_area / (mris->total_area-mris->neg_area)) :
              sqrt(mris->orig_area / mris->total_area) ;

  printf("scale = %g\n",scale);
  for (vno = 0 ; vno < mris->nvertices ; vno++){
    v = &mris->vertices[vno] ;
    vnum = v->vtotal ;
    if (v->ripflag || vnum <= 0)  continue ;

    V3_CLEAR(v_delta) ;
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;

    printf("%d %d -----------------------\n",vno,vnum);
    for (n = 0 ; n < vnum ; n++) {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag) continue ;
      d0 = v->dist_orig[n]/scale ;
      d0b = sqrt(pow(v->origx-vn->origx,2)+pow(v->origy-vn->origy,2)+
		 pow(v->origz-vn->origz,2));
      printf("%4d %7.4f %7.4f\n",n,d0,d0b);
      continue;

      dt = v->dist[n] ;
      delta = dt - d0 ;

      VECTOR_LOAD(v_y, vn->x - v->x, vn->y - v->y, vn->z - v->z) ;
      if ((V3_LEN_IS_ZERO(v_y))) continue ;
      V3_NORMALIZE(v_y, v_y) ;   /* make it a unit vector */
      V3_SCALAR_MUL(v_y, delta, v_y) ;
      V3_ADD(v_y, v_delta, v_delta) ;
      printf("%d %d %4d %g %g %g   (%g,%g,%g)\n",vno,n,v->v[n],
	     d0,dt,delta,
	     v_y->rptr[1][1],v_y->rptr[2][1],v_y->rptr[3][1]);
      printf("  c=[%g %g %g]; n=[%g %g %g];\n",
	     v->origx,v->origy,v->origz,
	     vn->origx,vn->origy,vn->origz);
      d0b = sqrt(pow(v->origx-vn->origx,2)+pow(v->origy-vn->origy,2)+
		 pow(v->origz-vn->origz,2));


    }

    V3_SCALAR_MUL(v_delta, norm, v_delta) ;

    /* take out normal component */
    nc = V3_DOT(v_n, v_delta) ;
    V3_SCALAR_MUL(v_n, -nc, v_n) ;
    V3_ADD(v_delta, v_n, v_delta) ;

    v->dx += l_dist * V3_X(v_delta) ;
    v->dy += l_dist * V3_Y(v_delta) ;
    v->dz += l_dist * V3_Z(v_delta) ;
    printf("   %d %g %g %g\n",vno,v->dx,v->dy,v->dz);
  exit(1);

  }

  VectorFree(&v_n) ;
  VectorFree(&v_y) ;
  VectorFree(&v_delta) ;
  return(NO_ERROR) ;
}

/*----------------------------------------------------*/
MRI *MRISmercatorGrid(MRIS *mris, double dtheta, double dphi)
{
  double r,x,y,z,d,theta,phi,f,mark;
  MRI *iso;
  int vno;
  VERTEX *vertex;

  f = 0.05;
  iso = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 1);
  r = MRISaverageRadius(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++){
    vertex = &mris->vertices[vno] ;
    x = vertex->x ;
    y = vertex->y ;
    z = vertex->z ;
    theta = atan2(y/r, x/r) ;
    if (theta < 0.0f) theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = r*r - z*z ;
    if (d < 0.0) d = 0 ;
    phi = atan2(sqrt(d), z) ;
    mark = 0;
    if(fabs(fmod(theta,dtheta))<f*dtheta) mark++;
    if(fabs(fmod(phi,dphi))<f*dphi) mark++;
    mark *= SIGN(z);
    MRIsetVoxVal(iso,vno,0,0,0, mark);
  }
  return(iso);
}


MRI *MRISsmoothMRIFast2(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *IncMask,  MRI *Targ)
{
  int nnbrs, nthstep, frame, vno, nthnbr, num, nvox, nbrvno;
  MRI *SrcTmp;
  struct timeb  mytimer;
  int msecTime;
  int *nNbrs, *nNbrs0, *rip, *rip0, nNbrsMax;
  float **pF, **pF0, *tF, *tF0, sumF;
  VERTEX *v, *vn;

  nvox = Src->width * Src->height * Src->depth;
  if (Surf->nvertices != nvox){
    printf("ERROR: MRISsmoothMRIFast(): Surf/Src dimension mismatch\n");
    return(NULL);
  }

  if(Targ == NULL){
    Targ = MRIallocSequence(Src->width, Src->height, Src->depth,
                            MRI_FLOAT, Src->nframes);
    if (Targ==NULL){
      printf("ERROR: MRISsmoothMRIFast(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(Src,Targ);
  }
  if(MRIdimMismatch(Src,Targ,1)){
    printf("ERROR: MRISsmoothFast(): output dimension mismatch\n");
    return(NULL);
  }
  if(Targ->type != MRI_FLOAT){
    printf("ERROR: MRISsmoothFast(): structure passed is not MRI_FLOAT\n");
    return(NULL);
  }

  // Alloc arrays. If there are ripped vertices, then only rip 
  // needs nvertices elements
  nNbrsMax = 12; // Should measure this, but overalloc does not hurt
  pF = (float **) calloc(Surf->nvertices*nNbrsMax,sizeof(float *));
  tF = (float *)  calloc(Surf->nvertices,sizeof(float));
  nNbrs = (int *) calloc(Surf->nvertices,sizeof(int));
  rip = (int *) calloc(Surf->nvertices,sizeof(int));

  pF0 = pF;
  tF0 = tF;
  rip0 = rip;
  nNbrs0 = nNbrs;

  TimerStart(&mytimer) ;
  SrcTmp = MRIcopy(Src,NULL);

  // Loop through frames
  for (frame = 0; frame < Src->nframes; frame ++){

    // Set up pointers for this frame
    pF = pF0;
    rip = rip0;
    nNbrs = nNbrs0;
    for (vno = 0 ; vno < Surf->nvertices ; vno++){
      v = &Surf->vertices[vno] ;
      if(IncMask && MRIgetVoxVal(IncMask,vno,0,0,0) < 0.5){
	// Mask is inclusive, so look for out of mask
	// should exclude rips here too? Original does not.
	rip[vno] = 1;
	MRIFseq_vox(SrcTmp,vno,0,0,frame) = 0;
	MRIFseq_vox(Targ,vno,0,0,frame) = 0;
	continue ;
      }
      rip[vno] = 0;
      *pF++ = (float *)(&(MRIFseq_vox(SrcTmp,vno,0,0,frame)));
      nnbrs = Surf->vertices[vno].vnum;
      num = 1;
      for (nthnbr = 0 ; nthnbr < nnbrs ; nthnbr++) {
	nbrvno = Surf->vertices[vno].v[nthnbr];
	vn = &Surf->vertices[nbrvno] ;
	if(vn->ripflag) continue ;
	if(IncMask && MRIgetVoxVal(IncMask,nbrvno,0,0,0) < 0.5) continue ;
	*pF++ = (float *)(&(MRIFseq_vox(SrcTmp,nbrvno,0,0,frame)));
	num++ ;
      }
      *nNbrs++ = num; // num takes into account all rips/masks
    }

    // Step through the iterations
    for(nthstep = 0; nthstep < nSmoothSteps; nthstep++){
      // Init pointers for this iteration
      pF  = pF0;
      rip = rip0;
      tF  = tF0;
      nNbrs = nNbrs0;
      // Loop through vertices, average nearest neighbors
      for (vno = 0 ; vno < Surf->nvertices ; vno++){
	if(*rip++) continue ;
	sumF = *(*pF++);
	for(nthnbr = 0 ; nthnbr < (*nNbrs)-1 ; nthnbr++) sumF += *(*pF++);
	*tF++ = sumF/(*nNbrs);
	nNbrs++;
      }
      // Load up for the next step
      rip = rip0;
      tF  = tF0;
      for (vno = 0 ; vno < Surf->nvertices ; vno++){
	if(*rip++) continue ;
	MRIsetVoxVal(SrcTmp,vno,0,0,frame, *tF++);
      }
    }

  }/* end loop over frame */

  // Copy to the output
  MRIcopy(SrcTmp,Targ);

  msecTime = TimerStop(&mytimer) ;
  if(Gdiag_no > 0){
    printf("MRISsmoothFast() nsteps = %d, tsec = %g\n",nSmoothSteps,msecTime/1000.0);
    fflush(stdout);
  }

  MRIfree(&SrcTmp);
  free(pF0);
  free(tF0);
  free(rip0);
  free(nNbrs0);

  return(Targ);
}
/*----------------------------------------------------------*/
int MRISsmoothMRIFastCheck2(int nSmoothSteps)
{
  char tmpstr[2000], *UFSM;
  MRIS *mris;
  MRI *src, *mri1, *mri2, *mask;
  int k,nerrs,c,r,s,f, cmax, rmax, smax, fmax;
  float val1, val2, diff,maxdiff;

  // Make sure to turn off override (restored later)
  UFSM = getenv("USE_FAST_SURF_SMOOTHER");
  unsetenv("USE_FAST_SURF_SMOOTHER");

  printf("MRISsmoothMRIFastCheck() nSmoothSteps = %d\n",nSmoothSteps);
  
  sprintf(tmpstr,"%s/subjects/fsaverage/surf/lh.white",getenv("FREESURFER_HOME"));
  printf("Reading surface %s\n",tmpstr);
  mris = MRISread(tmpstr);
  if(mris == NULL){
    printf("ERROR: could not read %s\n",tmpstr);
    return(-1);
  }

  // Use 3 frames
  src = MRIrandn(mris->nvertices, 1, 1, 3, .5, 1, NULL);

  // Create mask
  mask = MRIconst(mris->nvertices, 1, 1, 1, 1.0, NULL);
  for(k=0; k < mris->nvertices-1; k++){
    MRIsetVoxVal(mask,k,0,0,0, 0.0); // turn off mask
    mris->vertices[k+1].ripflag = 1; // rip a few
  }

  printf("Running slow smoother\n");
  mri1 = MRISsmoothMRI(mris, src, nSmoothSteps, mask, NULL);
  printf("Running fast smoother\n");
  mri2 = MRISsmoothMRIFast(mris, src, nSmoothSteps, mask, NULL);

  printf("Checking differences\n");
  nerrs = 0;
  cmax = 0;  rmax = 0;   smax = 0;  fmax = 0;
  maxdiff = 0.0;
  for (c=0; c < src->width; c++) {
    for (r=0; r < src->height; r++) {
      for (s=0; s < src->depth; s++) {
	for (f=0; f < src->nframes; f++) {
	  val1 = MRIgetVoxVal(mri1,c,r,s,f);
	  val2 = MRIgetVoxVal(mri2,c,r,s,f);
	  diff = val1-val2;
	  if(fabs(maxdiff) < fabs(diff)) {
	    maxdiff = diff;
	    cmax = c;
	    rmax = r;
	    smax = s;
	    fmax = f;
	  }
	  if(!FZERO(diff)) nerrs++;
	}
      }
    }
  }
  printf("nerrs = %d, maxdiff %f at %d %d %d %d\n",nerrs, maxdiff,cmax,rmax,smax,fmax);

  setenv("USE_FAST_SURF_SMOOTHER",UFSM,1);

  return(nerrs);
}
/*---------------------------------------------------------------------------------*/
int GCAprint(GCA *gca, FILE *fp)
{
  int n;
  GCA_TISSUE_PARMS *tp;

  fprintf(fp,"node_spacing %f\n",gca->node_spacing);
  fprintf(fp,"prior_spacing %f\n",gca->prior_spacing);
  fprintf(fp,"node_width %d\n",gca->node_width);
  fprintf(fp,"node_height %d\n",gca->node_height);
  fprintf(fp,"node_depth %d\n",gca->node_depth);
  fprintf(fp,"prior_width %d\n",gca->prior_width);
  fprintf(fp,"prior_height %d\n",gca->prior_height);
  fprintf(fp,"prior_depth %d\n",gca->prior_depth);
  fprintf(fp,"ninputs %d\n",gca->ninputs);
  for(n=0; n<gca->ninputs; n++){
    tp = &gca->tissue_parms[n];
    fprintf(fp,"TP %d  %d %f %f %f  %f %f %f\n",n,tp->label,
	    tp->T1_mean,tp->PD_mean,tp->T2_mean,
	    tp->T1_var,tp->PD_var,tp->T2_var);
  }
  //fprintf(fp," %f\n",gca->);
  return(0);
}

MRI * GCAmri2(GCA *gca, MRI *mri)
{
  int       frame, width, height, depth, x, y, z, xp, yp, zp, n, xn, yn, zn ;
  float     val ;
  GC1D      *gc ;
  GCA_PRIOR *gcap ;
  int nlabelsmax;

  width  = gca->node_width ;
  height = gca->node_height ;
  depth  = gca->node_depth ;

  if (!mri)
  {
    mri = MRIallocSequence(gca->node_width,
                           gca->node_height,
                           gca->node_depth,
                           MRI_UCHAR,
                           gca->ninputs) ;
    mri->xsize = gca->xsize*gca->node_spacing;
    mri->ysize = gca->ysize*gca->node_spacing;
    mri->zsize = gca->zsize*gca->node_spacing;
  }
  // in order to create the gca volume,
  // the volume must have the same direction cosines
  GCAcopyDCToMRI(gca, mri);

  nlabelsmax=0;
  for (x = 0 ; x < width ; x++)    {
    for (y = 0 ; y < height ; y++)      {
      for (z = 0 ; z < depth ; z++)        {
	if (!GCAvoxelToPrior(gca, mri, x, y, z, &xp, &yp, &zp)){
	  if (!GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn)){
	    gcap = &gca->priors[xp][yp][zp] ;
	    if (gcap==NULL) continue;
	    if(nlabelsmax < gcap->nlabels) nlabelsmax = gcap->nlabels;
	  }
	}
      }
    }
  }
  printf("nlabelsmax = %d\n",nlabelsmax);

  if (!mri)
  {
    mri = MRIallocSequence(gca->node_width,
                           gca->node_height,
                           gca->node_depth,
                           MRI_UCHAR,
                           gca->ninputs) ;
    mri->xsize = gca->xsize*gca->node_spacing;
    mri->ysize = gca->ysize*gca->node_spacing;
    mri->zsize = gca->zsize*gca->node_spacing;
  }
  // in order to create the gca volume,
  // the volume must have the same direction cosines
  GCAcopyDCToMRI(gca, mri);

  for (frame = 0 ; frame < gca->ninputs ; frame++)  {
    for (x = 0 ; x < width ; x++)    {
      for (y = 0 ; y < height ; y++)      {
        for (z = 0 ; z < depth ; z++)        {
          if (!GCAvoxelToPrior(gca, mri, x, y, z, &xp, &yp, &zp)){
            if (!GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn)){
              gcap = &gca->priors[xp][yp][zp] ;
              if (gcap==NULL) continue;
	      val = 0.0; 
              for (n = 0 ; n < gcap->nlabels ; n++) {
                gc = GCAfindGC(gca, xn, yn, zn, gcap->labels[n]) ;
                if(gc) val += gc->means[frame] * gcap->priors[n] ;
              }
              MRIsetVoxVal(mri, x, y, z, frame, val) ;
            }
	  }
        }
      }
    }
  }
  return(mri) ;
}

MRI *GCAnlables(GCA *gca, MRI *mri)
{
  int       width, height, depth, x, y, z, nthlabel;
  GCA_NODE *gcan;

  width  = gca->node_width ;
  height = gca->node_height ;
  depth  = gca->node_depth ;

  if (!mri)
  {
    mri = MRIallocSequence(gca->node_width,
                           gca->node_height,
                           gca->node_depth,
                           MRI_FLOAT,1);
    mri->xsize = gca->xsize*gca->node_spacing;
    mri->ysize = gca->ysize*gca->node_spacing;
    mri->zsize = gca->zsize*gca->node_spacing;
  }
  // in order to create the gca volume,
  // the volume must have the same direction cosines
  GCAcopyDCToMRI(gca, mri);

  for (x = 0 ; x < width ; x++) {
    for (y = 0 ; y < height ; y++) {
      for (z = 0 ; z < depth ; z++) {
	gcan = &gca->nodes[x][y][z];
	for(nthlabel = 0; nthlabel < gcan->nlabels; nthlabel++){
	  if(gcan->labels[nthlabel] == 17)
	    MRIsetVoxVal(mri,x,y,z,0,gcan->gcs->means[nthlabel]);
	}
      }
    }
  }
  return(mri) ;
}

/*----------------------------------------------------*/
MRI *MRIpolyfitBiasField(MRI *vol, int order, MRI *seg, MRI *bias)
{
  int c,r,s, nseg, nX, *segids, nsegids,z,k,segid,*tmplist,nthseg, npoly=0;
  MATRIX *X, *y, *Xt, *XtX, *iXtX, *Xty, *beta;
  double dc,dr,ds,val;
  double nchalf,nrhalf,nshalf;

  if(!bias) {
    bias = MRIallocSequence(vol->width, vol->height, vol->depth,MRI_FLOAT, 1);
    MRIcopyHeader(vol,bias);
  }
  if(MRIdimMismatch(vol, bias, 0)){
    printf("ERROR: MRIpolyfitBiasField(): vol/bias dim mismatch\n");
    return(NULL);
  }
  if(MRIdimMismatch(vol, seg, 0)){
    printf("ERROR: MRIpolyfitBiasField(): vol/seg dim mismatch\n");
    return(NULL);
  }

  // Count number of voxels in seg so can alloc
  nseg = 0;
  for (c=0; c < seg->width; c++)
    for (r=0; r < seg->height; r++)
      for (s=0; s < seg->depth; s++)
	if(MRIgetVoxVal(seg,c,r,s,0) > 0.5) nseg++;
  printf("MRIpolyfitBiasField(): found %d voxels in seg\n",nseg);

  // Get number of unique list segmentation IDs
  segids = MRIsegmentationList(seg, &nsegids);
  // Check whether there is a segmentation 0
  z = 0;
  for(k=0; k<nsegids;k++) if(segids[k] == 0) z = 1;
  if(z){
    tmplist = (int *) calloc(nsegids-1,sizeof(int));
    nthseg = 0;
    for(k=0; k<nsegids;k++) {
      if(segids[k] == 0) continue;
      tmplist[nthseg] = segids[k];
      nthseg ++;
    }
    free(segids);
    segids = tmplist;
    nsegids = nthseg;
  }
  printf("Found %d non-zero segids\n",nsegids);
  //for(k=0; k<nsegids;k++) printf("%2d %5d\n",k,segids[k]);
  
  // Alloc
  if(order == 2) npoly = 9;
  if(order == 3) npoly = 9+10;
  nX = nsegids + npoly; 
  X = MatrixAlloc(nseg,nX,MATRIX_REAL);
  y = MatrixAlloc(nseg,1,MATRIX_REAL);

  //Scale CRS to make X better conditioned
  nchalf = seg->width/2.0;
  nrhalf = seg->height/2.0;
  nshalf = seg->depth/2.0;

  // Set up the matrices to do the estimation
  nseg = 0;
  for(s=0; s < seg->depth; s++){
    ds = (s - nshalf)/nshalf; 
    for(c=0; c < seg->width; c++){
      dc = (c - nchalf)/nchalf;
      for(r=0; r < seg->height; r++){
	dr = (r - nrhalf)/nrhalf;
	segid = MRIgetVoxVal(seg,c,r,s,0);
	if(segid == 0) continue;
	for(k=0; k<nsegids;k++) if(segid == segids[k]) break;
	// Linear
	X->rptr[nseg+1][1] = dc;
	X->rptr[nseg+1][2] = dr;
	X->rptr[nseg+1][3] = ds;
	// Pure Quadratic
	X->rptr[nseg+1][4] = dc*dc;
	X->rptr[nseg+1][5] = dr*dr;
	X->rptr[nseg+1][6] = ds*ds;
	// Quadratic cross-products
	X->rptr[nseg+1][7] = dc*dr;
	X->rptr[nseg+1][8] = dc*ds;
	X->rptr[nseg+1][9] = dr*ds;
	if(order > 2){
	  // Cubic
	  X->rptr[nseg+1][10] = dc*dc*dc;
	  X->rptr[nseg+1][11] = dr*dr*dr;
	  X->rptr[nseg+1][12] = ds*ds*ds;
	  X->rptr[nseg+1][13] = dc*dc*dr;
	  X->rptr[nseg+1][14] = dc*dc*ds;
	  X->rptr[nseg+1][15] = dc*dr*ds;
	  X->rptr[nseg+1][16] = dc*dr*dr;
	  X->rptr[nseg+1][17] = dc*ds*ds;
	  X->rptr[nseg+1][18] = dr*dr*ds;
	  X->rptr[nseg+1][19] = dr*ds*ds;
	}
	// Constant
	X->rptr[nseg+1][npoly+k+1] = 1.0;
	// Input data
	y->rptr[nseg+1][1] = MRIgetVoxVal(vol,c,r,s,0);
	nseg++;
      }
    }
  }

  // Do the estimation
  Xt   = MatrixTranspose(X,NULL);
  XtX  = MatrixMultiplyD(Xt,X,NULL);
  iXtX = MatrixInverse(XtX,NULL);
  Xty  = MatrixMultiplyD(Xt,y,NULL);
  beta = MatrixMultiplyD(iXtX,Xty,NULL);

  //MatrixWrite(X,"X.mat","X");
  //MatrixWrite(y,"y.mat","y");
  //MatrixWrite(beta,"beta.mat","beta");
  //MatrixWrite(Xt,"Xt.mat","Xt");
  //MatrixWrite(XtX,"XtX.mat","XtX");
  //MatrixWrite(iXtX,"iXtX.mat","iXtX");

  MatrixFree(&y);
  MatrixFree(&X);
  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&Xty);

  // Use the estimation to compute bias field at all voxels
  X = MatrixAlloc(1,nX,MATRIX_REAL);
  for(s=0; s < seg->depth; s++){
    ds = (s - nshalf)/nshalf; //norm makes X better conditioned
    for(c=0; c < seg->width; c++){
      dc = (c - nchalf)/nchalf;
      for(r=0; r < seg->height; r++){
	dr = (r - nrhalf)/nrhalf;
	X->rptr[1][1] = dc;
	X->rptr[1][2] = dr;
	X->rptr[1][3] = ds;
	X->rptr[1][4] = dc*dc;
	X->rptr[1][5] = dr*dr;
	X->rptr[1][6] = ds*ds;
	X->rptr[1][7] = dc*dr;
	X->rptr[1][8] = dc*ds;
	X->rptr[1][9] = dr*ds;
	if(order > 2){
	  // Cubic
	  X->rptr[1][10] = dc*dc*dc;
	  X->rptr[1][11] = dr*dr*dr;
	  X->rptr[1][12] = ds*ds*ds;
	  X->rptr[1][13] = dc*dc*dr;
	  X->rptr[1][14] = dc*dc*ds;
	  X->rptr[1][15] = dc*dr*ds;
	  X->rptr[1][16] = dc*dr*dr;
	  X->rptr[1][17] = dc*ds*ds;
	  X->rptr[1][18] = dr*dr*ds;
	  X->rptr[1][19] = dr*ds*ds;
	}
	// Note: all constant terms excluded
	y = MatrixMultiply(X,beta,y);
	val = MRIgetVoxVal(vol,c,r,s,0) - y->rptr[1][1];
	MRIsetVoxVal(bias,c,r,s,0,val);
      }
    }
  }
  MatrixFree(&y);
  MatrixFree(&X);
  free(segids);
  return(bias);

}



/* ----------------------------------------------------------
   MRIsegCount() - returns the number of times the given
   segmentation id appears in the volume.
   --------------------------------------------------------- */
int MRIsegCount(MRI *seg, int id, int frame)
{
  int nhits, v, c,r,s;
  nhits = 0;
  for (c=0; c < seg->width; c++)
  {
    for (r=0; r < seg->height; r++)
    {
      for (s=0; s < seg->depth; s++)
      {
        v = (int) MRIgetVoxVal(seg,c,r,s,frame);
        if (v == id)
        {
          nhits ++;
        }
      }
    }
  }
  return(nhits);
}

float MRIvoxelsInLabelWithPartialVolumeEffects2( const MRI *mri,
						const MRI *mri_vals, 
						const int label, 
						MRI *mri_mixing_coef, 
						MRI *mri_nbr_labels ) {
  enum { maxlabels = 20000 };
  float volume;
  int     x, y, z;
  MRI     *mri_border ;
  // DNG 6/7/07 : had to use maxlabels instead of MAX_CMA_LABELS here
  // so that segmentations with values > MAX_CMA_LABELS can be
  // accessed. This includes the cortical segmentations as well as
  // white matter segs. Currently, the max seg no is 4181, but this
  // could easily change.
  // NJS 2/17/10 : Indeed, it did change... the Destrieux a2009s atlas has
  // label values up to about 15000.

  if(label >= maxlabels){
    printf("ERROR: MRIvoxelsInLabelWithPartialVolumeEffects()\n");
    printf(" label %d exceeds maximum label number %d\n",label,maxlabels);
    return(-100000);
  }
  const float vox_vol = mri->xsize*mri->ysize*mri->zsize ;

  /* first find border voxels */
  mri_border = MRImarkLabelBorderVoxels(mri, NULL, label, 1, 1);

  if( DIAG_VERBOSE_ON && (Gdiag & DIAG_WRITE) ) {
    MRIwrite(mri_border, "b.mgz");
  }

  volume = 0 ;
  for( x = 0 ; x < mri->width ; x++ ) {
    for( y = 0 ; y < mri->height ; y++ ) {
      for( z = 0 ; z < mri->depth ; z++ ) {
        const int vox_label = MRIgetVoxVal(mri, x, y, z, 0);
        const int border = MRIgetVoxVal(mri_border, x, y, z, 0);
	/* Note that these are all zeroed at the start of MRIcomputeLabelNbhd*/
	int nbr_label_counts[maxlabels], label_counts[maxlabels];
	float label_means[maxlabels];

	// Not in label and not a border
        if( (vox_label != label) && (border == 0) ) continue;

        if( border == 0 ) {
	  volume += vox_vol;
        } else { /* compute partial volume */
          MRIcomputeLabelNbhd( mri, mri_vals, x, y, z,
			       nbr_label_counts, label_means, 1, maxlabels);

          MRIcomputeLabelNbhd( mri, mri_vals, x, y, z,
			       label_counts, label_means, 7, maxlabels );

	  // Compute partial volume based on intensity
          const float val = MRIgetVoxVal( mri_vals, x, y, z, 0 );
          float mean_label = label_means[vox_label];
          int nbr_label = -1 ;
          int max_count = 0 ;
	  float pv, mean_nbr;

	  /*for a label that is a nbr and is on the other side of val from the label mean */
	  int this_label;
          for( this_label = 0 ; this_label < maxlabels;  this_label++ ) {
            if( this_label == vox_label ) continue ;
            if( nbr_label_counts[this_label] == 0 )  continue ; /* not a nbr */

            if( (label_counts[this_label] > max_count) &&
                ((label_means[this_label] - val) *
                 (mean_label - val) < 0) ) {
              max_count = label_counts[this_label] ;
              nbr_label = this_label ;
            }
          }
	  //printf("%3d %3d %3d nbr_label %d, maxc = %d\n",x,y,z,nbr_label,max_count);

          if( vox_label != label && nbr_label != label ) {
            continue; // this struct not in voxel 
	  }

          if( max_count == 0 ) {
            volume += vox_vol ; // couldn't find an appropriate label

            if (mri_nbr_labels) {
	      // find max nbr label anyway for caller
              for( this_label = 0; this_label < maxlabels;  this_label++ ) {

                if( this_label == vox_label ) {
                  continue;
		}

                if( nbr_label_counts[this_label] == 0 ) {
                  continue ; /* not a nbr */
		}
                
                if( label_counts[this_label] > max_count ) {
                  max_count = label_counts[this_label] ;
                  nbr_label = this_label ;
                }
              }

              MRIsetVoxVal(mri_nbr_labels, x, y, z, 0, nbr_label) ;
              if( mri_mixing_coef ) {
                MRIsetVoxVal(mri_mixing_coef, x, y, z, 0, 1.0);
	      }

            }

          } else {
	    // compute partial volume pct 
            mean_nbr = label_means[nbr_label] ;
            pv = (val - mean_nbr) / (mean_label - mean_nbr) ;

            if (pv < 0 || pv > 1) {
              DiagBreak() ;
	    }

            if (pv > 1) {
              pv = 1 ;
	    }

            if (pv < 0) {
              continue ;  // shouldn't happen
	    }

            if( vox_label != label ) {
              pv = 1-pv ;
	    }

            volume += vox_vol * pv ;

            if( mri_mixing_coef ) {
              MRIsetVoxVal(mri_mixing_coef, x, y, z, 0, pv);
	    }

            if (mri_nbr_labels) {
	      // return nbr label to caller
              if (vox_label != label) {
                MRIsetVoxVal(mri_nbr_labels, x, y, z, 0, vox_label);
              } else {
                MRIsetVoxVal(mri_nbr_labels, x, y, z, 0, nbr_label);
	      }
            }
          }
        }
      }
    }
  }
  
  MRIfree(&mri_border) ;

  return(volume) ;
}
int ProjSurf(MRIS *surf, MRI *d, double f)
{
  VERTEX *v;
  int vtx;
  float Tx, Ty, Tz, pd;

  for (vtx = 0; vtx < surf->nvertices; vtx++)
  {
    v = &surf->vertices[vtx] ;
    pd = f*MRIgetVoxVal(d,vtx,0,0,0);
    //if(pd < 0) continue;
    ProjNormDist(&Tx,&Ty,&Tz,surf,vtx,pd);
    v->x = Tx;
    v->y = Ty;
    v->z = Tz;
  }
  return(0);
}



MRI_SP *MRISmakeTemplateB(int nsubjects, char **subjlist, int nhemis, char **hemilist, char *surfregname)
{
  static char *surface_names[] = {"inflated","smoothwm","smoothwm"} ;
  //static char *curvature_names[] =  {"inflated.H","sulc",NULL} ;
  static char *curvature_names[] =  {"inflated.H","sulc","smoothwm.H"} ;
  char tmpstr[2000];
  int images_per_surface = 3;
  int nsurfaces = sizeof(curvature_names) / sizeof(curvature_names[0]);
  int nparam_images = images_per_surface*nsurfaces;
  float scale = 1 ;
  char *annot_name = "aparc", *SUBJECTS_DIR, *hemi, *subject ;
  INTEGRATION_PARMS parms ;
  int which_norm = NORM_MEAN;
  int navgs = 0, nthhemi, sno, nthsubject,err, nbrs=3 ; 
  MRI_SP  *mrisp, /* *mrisp_aligned,*/ *mrisp_template ;
  MRIS *mris;

  /* default template fields*/
  memset(&parms, 0, sizeof(parms)) ;
  parms.nfields=3;
  SetFieldLabel(&parms.fields[0], INFLATED_CURV_CORR_FRAME,0,0.0,0.0,0,which_norm);
  /* only use sulc for rigid registration */
  SetFieldLabel(&parms.fields[1],SULC_CORR_FRAME,1,1.0,0.0,0,which_norm);
  SetFieldLabel(&parms.fields[2],CURVATURE_CORR_FRAME,2,0.0,0.0,0,which_norm);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  mrisp_template = MRISPalloc(scale, nparam_images);
  for(nthsubject = 0; nthsubject < nsubjects; nthsubject++){
    subject = subjlist[nthsubject];

    for(nthhemi = 0; nthhemi < nhemis; nthhemi++){
      hemi = hemilist[nthhemi];
      printf("subject %s hemi %s\n",subject,hemi);
      sprintf(tmpstr, "%s/%s/surf/%s.%s",SUBJECTS_DIR, subject, hemi, surfregname);
      printf("   reading surface %s...\n",tmpstr);
      mris = MRISread(tmpstr);
      if(mris == NULL){
	printf("ERROR: could not load %s\n",tmpstr);
	return(NULL);
      }
      
      err = MRISreadAnnotation(mris, annot_name);
      if(err){
	printf("ERROR: could not load %s\n",annot_name);
	return(NULL);
      }

      MRISripMedialWall(mris) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
      MRISstoreMetricProperties(mris) ;
      
      for (sno = 0; sno < nsurfaces ; sno++){
	if(curvature_names[sno]){
	  /* read in precomputed curvature file */
	  sprintf(tmpstr, "%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,curvature_names[sno]) ;
	  err = MRISreadCurvatureFile(mris, tmpstr);
	  if(err){
	    printf("ERROR: could not load %s\n",tmpstr);
	    return(NULL);
	  }
	  MRISaverageCurvatures(mris, navgs) ;
	  MRISnormalizeCurvature(mris, which_norm) ;
	}
	else {
	  sprintf(tmpstr, "%s/%s/surf/%s.%s",SUBJECTS_DIR,subject, hemi,surface_names[sno]) ;
	  printf("   Loading surface %s\n",tmpstr);
	  err = MRISreadVertexPositions(mris, tmpstr);
	  if(err){
	    printf("ERROR: could not load %s\n",tmpstr);
	    return(NULL);
	  }
	  if (nbrs > 1) MRISsetNeighborhoodSize(mris, nbrs) ;
	  MRIScomputeMetricProperties(mris) ;
	  MRIScomputeSecondFundamentalForm(mris) ;
	  MRISuseMeanCurvature(mris) ;

	  sprintf(tmpstr,"%s.%s.curv",hemi,subject);
	  printf("writing %s\n",tmpstr);
	  MRISwriteCurvature(mris, tmpstr) ;

	  MRISaverageCurvatures(mris, navgs) ;
	  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
	  MRISnormalizeCurvature(mris, which_norm) ;
	}
	printf("  computing parameterization for surface %s...\n",tmpstr);
	mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
	MRISPcombine(mrisp, mrisp_template, sno*3) ;
	MRISPfree(&mrisp) ;
      }
      MRISfree(&mris) ;
    }
  }
  return(mrisp_template);
}

MRI_SP *MRISmakeTemplateC(int nsurfaces, char **surfregflist, char **annotflist, 
			  int nparams, char **paramflist, int navgs)
{
  int which_norm = NORM_MEAN;
  float scale = 1 ;
  INTEGRATION_PARMS parms ;
  int nthparam,nthparamfile,nthsurf, err ;
  MRIS *mris;
  MRI_SP  *mrisp, /* *mrisp_aligned,*/ *mrisp_template ;

  /* default template fields*/
  memset(&parms, 0, sizeof(parms)) ;
  parms.nfields=3;
  SetFieldLabel(&parms.fields[0], INFLATED_CURV_CORR_FRAME,0,0.0,0.0,0,which_norm);
  /* only use sulc for rigid registration */
  SetFieldLabel(&parms.fields[1],SULC_CORR_FRAME,1,1.0,0.0,0,which_norm);
  SetFieldLabel(&parms.fields[2],CURVATURE_CORR_FRAME,2,0.0,0.0,0,which_norm);

  mrisp_template = MRISPalloc(scale, 3*nparams);// 3 = avg,std,dof

  nthparamfile = 0;
  for(nthsurf=0; nthsurf < nsurfaces; nthsurf++){
    printf("surf %2d ------------------------\n",nthsurf);
    printf("   %s\n",surfregflist[nthsurf]);
    mris = MRISread(surfregflist[nthsurf]);
    if(mris == NULL){
      printf("ERROR: could not load %s\n",surfregflist[nthsurf]);
      return(NULL);
    }

    if(annotflist){
      printf("   annot %s\n",annotflist[nthsurf]);
      err = MRISreadAnnotation(mris, annotflist[nthsurf]);
      if(err){
	printf("ERROR: could not load %s\n",annotflist[nthsurf]);
	return(NULL);
      }
      MRISripMedialWall(mris) ;
    }

    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    MRISstoreMetricProperties(mris) ;

    for(nthparam=0; nthparam < nparams; nthparam++){
      printf("   surf %2d param %2d %s\n",nthsurf,nthparam,paramflist[nthparam]);
      err = MRISreadCurvatureFile(mris, paramflist[nthparam]);
      if(err){
	printf("ERROR: could not load %s\n",paramflist[nthparam]);
	return(NULL);
      }
      MRISaverageCurvatures(mris, navgs) ;
      MRISnormalizeCurvature(mris, which_norm) ;

      mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
      MRISPcombine(mrisp, mrisp_template, nthparam*3) ;
      MRISPfree(&mrisp) ;
      nthparamfile++;
    }
    MRISfree(&mris) ;
  }
  return(mrisp_template);
}

MRI *fMRIrepMat(MATRIX *y, MRI *template, MRI *mask, MRI *out)
{
  int c, r, s, f;

  if(out==NULL){
    out = MRIallocSequence(template->width,template->height,template->depth,MRI_FLOAT,y->rows);
    if(out==NULL){
      printf("ERROR: fMRIrepMat(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(template,out);
  }

  for(c=0; c < out->width; c++)  {
    for(r=0; r < out->height; r++)    {
      for(s=0; s < out->depth; s++)   {
	if(mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5){
	  for(f=0; f < out->nframes; f++) MRIFseq_vox(out,c,r,s,f) = 0;
	  continue;
	}
	for(f=0; f < out->nframes; f++)
	  MRIsetVoxVal(out,c,r,s,f,y->rptr[f+1][1]);
      }
    }
  }
  return(out);
}


FSGD *LoadFactorFile(char *fname)
{
  FILE *fp;
  char tag[1000],tmpstr[1000];
  int r,n;
  FSGD_FACTOR *f;
  FSGD_FACTOR_CONTRAST *c;
  FSGD *gd;

  gd = gdfAlloc(1);
  gd->nFactors = 0;
  gd->nFC = 0;

  fp = fopen(fname,"r");
  if (fp==NULL) {
    printf("ERROR: LoadFactorFile: cannot open %s for reading\n",fname);
    return(NULL);
  }

  /*------- begin input loop --------------*/
  while (1) {

    r = fscanf(fp,"%s",tag);
    if(r==EOF) break;

    printf("fsgd tag: %s\n",tag);
    //if(Gdiag_no > 0) printf("fsgd tag: %s\n",tag);

    if(!strcasecmp(tag,"dfactor")) {
      f = (FSGD_FACTOR*)calloc(1,sizeof(FSGD_FACTOR));
      r = fscanf(fp,"%s",f->name);
      if(r==EOF) goto formaterror;
      f->type = FSGD_FACTOR_DISCRETE;
      f->nLevels = gdfCountItemsOnLine(fp);
      if(f->nLevels < 2){
	printf("ERROR: LoadFactorFile: factor %s only has %d levels\n",f->name,f->nLevels);
	goto formaterror;	
      }
      printf("D %s nlevels = %d\n",f->name,f->nLevels);
      for(n=0; n < f->nLevels; n++){
	fscanf(fp,"%s",tmpstr);
	f->Levels[n] = strcpyalloc(tmpstr);
	printf("  %d %s\n",n+1,f->Levels[n]);
      }
      gd->Factors[gd->nFactors] = f;
      gd->nFactors++;
      continue;
    }
    if(!strcasecmp(tag,"cfactor")) {
      f = (FSGD_FACTOR*)calloc(1,sizeof(FSGD_FACTOR));
      r = fscanf(fp,"%s",f->name);
      if(r==EOF) goto formaterror;
      f->type = FSGD_FACTOR_CONTINUOUS;
      f->nFactors = gdfCountItemsOnLine(fp);
      printf("C %s %d\n",f->name,f->nFactors);
      for(n=0; n < f->nFactors; n++){
	fscanf(fp,"%s",tmpstr);
	f->FactorNames[n] = strcpyalloc(tmpstr);
	printf("  %d %s\n",n+1,f->FactorNames[n]);
      }
      gd->Factors[gd->nFactors] = f;
      gd->nFactors++;
      continue;
    }
    if(!strcasecmp(tag,"maineffect")) {
      c = (FSGD_FACTOR_CONTRAST*)calloc(1,sizeof(FSGD_FACTOR_CONTRAST));
      r = fscanf(fp,"%s",tmpstr);
      if(r==EOF) goto formaterror;
      c->FactorNames[0] = strcpyalloc(tmpstr);
      c->type = FSGD_FACTOR_CONTRAST_MAIN;
      c->nFactors = 1;
      printf("MF %s \n",c->FactorNames[0]);
      gd->fc[gd->nFC] = c;
      gd->nFC++;
      continue;
    }
    if(!strcasecmp(tag,"interaction")) {
      c = (FSGD_FACTOR_CONTRAST*)calloc(1,sizeof(FSGD_FACTOR_CONTRAST));
      //r = fscanf(fp,"%s",c->name);
      //if(r==EOF) goto formaterror;
      c->type = FSGD_FACTOR_CONTRAST_INTERACTION;
      c->nFactors = gdfCountItemsOnLine(fp);
      printf("IX %d\n",c->nFactors);
      if(c->nFactors < 2){
	printf("ERROR: LoadFactorFile: interaction only has %d factors\n",c->nFactors);
	goto formaterror;	
      }
      for(n=0; n < c->nFactors; n++){
	fscanf(fp,"%s",tmpstr);
	c->FactorNames[n] = strcpyalloc(tmpstr);
	printf("  %d %s\n",n+1,c->FactorNames[n]);
      }
      gd->fc[gd->nFC] = c;
      gd->nFC++;
      continue;
    }
  }

  return(gd);

 formaterror:
  printf("FSGDF Format Error: file = %s, tag=%s\n",fname,tag);
  return(NULL);
}

int FSGDFactorsPrint(FSGD *gd, FILE *fp)
{
  int n,m;
  FSGD_FACTOR *f;
  FSGD_FACTOR_CONTRAST *c;

  for(n=0; n < gd->nFactors; n++){
    f = gd->Factors[n];
    if(f->type == FSGD_FACTOR_DISCRETE){
      fprintf(fp,"D %s nlevels = %d\n",f->name,f->nLevels);
      for(m=0; m < f->nLevels; m++){
	fprintf(fp,"  %d %s\n",m+1,f->Levels[m]);
      }
    }
    if(f->type == FSGD_FACTOR_CONTINUOUS){
      fprintf(fp,"C %s %d\n",f->name,f->nFactors);
      for(m=0; m < f->nFactors; m++){
	fprintf(fp,"  %d %s\n",m+1,f->FactorNames[m]);
      }
    }
  }

  for(n=0; n < gd->nFC; n++){
    c = gd->fc[n];
    if(c->type == FSGD_FACTOR_CONTRAST_MAIN)
      fprintf(fp,"MainEffect %s",c->FactorNames[0]);
    if(c->type == FSGD_FACTOR_CONTRAST_INTERACTION){
      fprintf(fp,"Interaction %d ",c->nFactors);
      for(m=0; m < c->nFactors; m++){
	fprintf(fp,"  %s",c->FactorNames[m]);
      }
    }
    fprintf(fp,"\n");
  }

  return(0);
}

// Checks:
// 1. Level names not repeated
// 2. Factor names not repeated
// 3. Contrasts are distinct
// 4. Contrast factors exist
// 5. Interactions only have 1 CF
// 6. Interactions with CF are within the model

// Add: tContrast and FContrast

// Spec Interaction Model instead of full DODS
//   -- read dfactors after continuous factor


// dmax, proj

MRI *MRIsurf2VolOpt(MRI *tempvol, MRIS **surfs, MRI **overlays, int nsurfs, 
		    MRI *ribbon, MATRIX *R, MRI *volsurf)
{
  int n,c,r,s,f,nmin, vtxno,vtxnomin=0, nframes, ribval;
  MHT **hash=NULL;
  int UseHash = 1;
  MATRIX *T, *invR, *M, *surfRAS=NULL,*crs;
  VERTEX v;
  float dmin, d, val;

  for(n=0; n<nsurfs; n++){
    //if(surfs[n]->hemisphere == LEFT_HEMISPHERE) continue;
    //if(surfs[n]->hemisphere == RIGHT_HEMISPHERE) continue;
    if(overlays[n]->nframes != overlays[0]->nframes){
      printf("ERROR: MRIsurf2VolOpt(): overlay dim mismatch %d\n",n);
      return(NULL);
    }
  }
  if(MRIdimMismatch(tempvol,ribbon,0)){
    printf("ERROR: MRIsurf2VolOpt(): tempvol/ribbon mismatch\n");
    return(NULL);
  }

  nframes = overlays[0]->nframes;
  if(volsurf == NULL){
    volsurf = MRIallocSequence(tempvol->width, tempvol->height, tempvol->depth,
                              MRI_FLOAT, nframes);
    if (volsurf==NULL){
      printf("ERROR: MRIsurf2VolOpt(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(tempvol,volsurf);
  }

  if(UseHash){
    hash = (MHT **) calloc(sizeof(MHT *),nsurfs);
    for(n=0; n<nsurfs; n++){
      hash[n] = MHTfillVertexTableRes(surfs[n], NULL,CURRENT_VERTICES,16);
    }
  }

  // M converts tempvol CRS to surface RAS
  T = MRIxfmCRS2XYZtkreg(tempvol);
  invR = MatrixInverse(R,NULL);
  M = MatrixMultiply(invR,T,NULL);
  MatrixFree(&T);
  MatrixFree(&invR);

  crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[4][1] = 1;
  for(c=0; c < tempvol->width; c++){
    for(r=0; r < tempvol->height; r++){
      for(s=0; s < tempvol->depth; s++){
	ribval = MRIgetVoxVal(ribbon,c,r,s,0);
	if(ribval != 3 && ribval != 42) {
	  //for(f=0; f < nframes; f++) MRIsetVoxVal(volsurf,c,r,s,f, 0);
	  continue;
	}
	crs->rptr[1][1] = c;
	crs->rptr[2][1] = r;
	crs->rptr[3][1] = s;
	surfRAS = MatrixMultiply(M,crs,surfRAS);
	v.x = surfRAS->rptr[1][1];
	v.y = surfRAS->rptr[2][1];
	v.z = surfRAS->rptr[3][1];
	dmin = 1000;
	nmin = -1;
	for(n=0; n<nsurfs; n++){
	  if(surfs[n]->hemisphere == LEFT_HEMISPHERE  && ribval !=  3) continue;
	  if(surfs[n]->hemisphere == RIGHT_HEMISPHERE && ribval != 42) continue;

	  if(UseHash) vtxno = MHTfindClosestVertexNo(hash[n],surfs[n],&v,&d);
	  else        vtxno = MRISfindClosestVertex(surfs[n],v.x,v.y,v.z,&d);
	  //printf("%d %d %d %d %g %g %g  %5d\n",n,c,r,s,v.x,v.y,v.z,vtxno);
	  if(vtxno < 0){
	    printf("MRIsurf2VolOpt(): No Match: %3d %3d %3d  %6.2f %6.2f %6.2f\n",c,r,s,v.x,v.y,v.z);
	    vtxno = MRISfindClosestVertex(surfs[n],v.x,v.y,v.z,&d);
	    printf("%d %d %d %d %g %g %g  %5d\n",n,c,r,s,v.x,v.y,v.z,vtxno);
	    //continue;
	  }
	  if(d < dmin){
	    dmin = d;
	    nmin = n;
	    vtxnomin = vtxno;
	  }
	} // surfs
	if(nmin == -1){
	  //printf("MRIsurf2VolOpt(): No Match: %3d %3d %3d  %6.2f %6.2f %6.2f\n",c,r,s,v.x,v.y,v.z);
	  continue;
	}
	for(f=0; f < nframes; f++){
	  val = MRIgetVoxVal(overlays[nmin],vtxnomin,0,0,f);
	  MRIsetVoxVal(volsurf,c,r,s,f, val);
	}
      } // slice
    } // row
  } //col

  if(UseHash) for(n=0; n<nsurfs; n++) if(UseHash) MHTfree(&hash[n]);

  MatrixFree(&surfRAS);
  return(volsurf);
}

MRI *FillSeg(MRI *seg, COLOR_TABLE *ctab, MRI *vol)
{
  int c,r,s,f,segid;
  CTE *cte;

  if(vol == NULL){
    vol = MRIallocSequence(seg->width, seg->height, seg->depth,
			   MRI_FLOAT, 1);
    if (vol==NULL){
      printf("ERROR: FillSeg(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(seg,vol);
  }
  if(MRIdimMismatch(seg,vol,0)){
      printf("ERROR: FillSeg(): dim mismatch\n");
    return(NULL);
  }
  
  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	segid = MRIgetVoxVal(seg,c,r,s,0);
	if(segid >= ctab->nentries){
	  printf("ERROR: segid = %d, nentries = %d\n",segid,ctab->nentries);
	  return(NULL);
	}
	cte = ctab->entries[segid];
	if(cte == NULL){
	  printf("ERROR: segid = %d not in ctab,\n",segid);
	  return(NULL);
	}
	for(f=0; f < vol->nframes; f++)
	  MRIsetVoxVal(vol,c,r,s,f, cte->ai);
      }
    }
  }
  return(vol);
}

MRI *MRIfcIntrinsicLI(MRI *lh, MRI *rh, double DenThresh)
{
  MRI *lhn, *rhn, *iLI;
  int c,f, nrois, roi1, roi2, nframes;
  double v,ss, v1,v2,den,num, LL, LR, RL, RR;

  if(lh->nframes != rh->nframes){
    printf("ERROR: MRIfcIntrinsicLI(): frame mismatch\n");
    return(NULL);
  }
  nframes = lh->nframes;

  if(lh->width != rh->width){
    printf("ERROR: MRIfcIntrinsicLI(): roi mismatch\n");
    return(NULL);
  }
  nrois = lh->width;

  iLI = MRIalloc(nrois,nrois,1,MRI_FLOAT);
  if(iLI == NULL){
    printf("ERROR: MRIfcIntrinsicLI(): could not alloc %d\n",nrois);
    return(NULL);
  }

  lhn = MRIcopy(lh,NULL);
  rhn = MRIcopy(rh,NULL);

  for(c=0; c < nrois; c++){
    ss = 0;
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(lh,c,0,0,f);
      ss += (v*v);
    }
    ss = sqrt(ss);
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(lh,c,0,0,f);
      MRIsetVoxVal(lhn,c,0,0,f, v/ss);
    }
    ss = 0;
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(rh,c,0,0,f);
      ss += (v*v);
    }
    ss = sqrt(ss);
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(rh,c,0,0,f);
      MRIsetVoxVal(rhn,c,0,0,f, v/ss);
    }
  }
  //MRIwrite(lhn,"lhn.nii");
  //MRIwrite(rhn,"rhn.nii");


  for(roi1 = 0; roi1 < nrois; roi1++){
    for(roi2 = 0; roi2 < nrois; roi2++){

      LL = 0;
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(lhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(lhn,roi2,0,0,f);
	LL += (v1*v2);
      }
      LR = 0;
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(lhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(rhn,roi2,0,0,f);
	LR += (v1*v2);
      }
      RL = 0; // not the same as LR (matrices are transposes)
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(rhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(lhn,roi2,0,0,f);
	RL += (v1*v2);
      }
      RR = 0;
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(rhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(rhn,roi2,0,0,f);
	RR += (v1*v2);
      }

      num = ((LL-RL)-(RR-LR));
      den = (fabs(LL)+fabs(LR)+fabs(RR)+fabs(RL));
      if(den>DenThresh) v = num/den;
      else              v = 0.0;
      MRIsetVoxVal(iLI,roi1,roi2,0,0, v);

    } // roi2
  } // roi1

  MRIfree(&lhn);
  MRIfree(&rhn);
  return(iLI);
}

MRI *NNGLMPVC(MRI *src, MRI *pvf, MRI *mask, int nrad, MRI *pvc)
{
  MATRIX *X,  *y, *beta=NULL, *Xt=NULL, *XtX=NULL, *Xty=NULL, *iXtX=NULL, *Xsum,*ytmp,*Xtmp;
  int c,r,s, dc,dr,ds, tt, nvmax, nth, nhits, nv, nkeep;
  double v;
  
  if(MRIdimMismatch(src,pvf,0)){
      printf("ERROR: NNGLMPVC(): src-pvf dim mismatch\n");
    return(NULL);
  }
  if(pvc == NULL){
    pvc = MRIallocSequence(pvf->width, pvf->height, pvf->depth,
			   MRI_FLOAT, pvf->nframes);
    if(pvc==NULL){
      printf("ERROR: NNGLMPVC(): could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(src,pvc);
  }
  if(MRIdimMismatch(src,pvc,0)){
      printf("ERROR: NNGLMPVC(): src-pvc dim mismatch\n");
    return(NULL);
  }

  nvmax = (2*nrad+1)*(2*nrad+1)*(2*nrad+1);
  printf("nvmax = %d\n",nvmax);

  nhits = 0;
  for(c=0; c < src->width; c++){
    printf("%2d ",c);
    for(r=0; r < src->height; r++){
      for(s=0; s < src->depth; s++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue; 
	y = MatrixAlloc(nvmax,1,MATRIX_REAL);
	X = MatrixAlloc(nvmax,pvf->nframes,MATRIX_REAL);

	nth = 0;
	y->rptr[nth+1][1] = MRIgetVoxVal(src,c,r,s,0);
	for(tt=0; tt < pvf->nframes; tt++)
	  X->rptr[nth+1][tt+1] = MRIgetVoxVal(pvf,c,r,s,tt);
	for(dc = -nrad; dc <= nrad; dc++){
	  if(c+dc < 0 || c+dc >= src->width) continue; 
	  for(dr = -nrad; dr <= nrad; dr++){
	    if(r+dr < 0 || r+dr >= src->height) continue; 
	    for(ds = -nrad; ds <= nrad; ds++){
	      if(s+ds < 0 || s+ds >= src->depth) continue; 
	      y->rptr[nth+1][1] = MRIgetVoxVal(src,c+dc,r+dr,s+ds,0);
	      for(tt=0; tt < pvf->nframes; tt++){
		v  = MRIgetVoxVal(pvf,c+dc,r+dr,s+ds,tt);
		X->rptr[nth+1][tt+1] = v;
	      }
	      nth++;
	    } // ds
	  } //dr
	} //dc
	nv = nth;
	Xsum = MatrixSum(X,1,NULL);
	nkeep = 0;
	for(tt=0; tt < pvf->nframes; tt++) if(Xsum->rptr[1][tt+1]>.1) nkeep++;
	if(nkeep != pvf->nframes || nv != nvmax){
	  ytmp = MatrixAlloc(nv,1,MATRIX_REAL);
	  Xtmp = MatrixAlloc(nv,nkeep,MATRIX_REAL);
	  for(nth = 0; nth < nv; nth++){
	    ytmp->rptr[nth+1] = y->rptr[nth+1];
	    nkeep = 0;
	    for(tt=0; tt < pvf->nframes; tt++){
	      if(Xsum->rptr[1][tt+1]>.1){
		Xtmp->rptr[nth+1][nkeep+1] = Xtmp->rptr[nth+1][tt+1];
		nkeep++;
	      }
	    }
	  }
	  MatrixFree(&y);	
	  MatrixFree(&X);	
	  y = ytmp;
	  X = Xtmp;
	}
	Xt = MatrixTranspose(X,NULL);
	XtX = MatrixMultiply(Xt,X,NULL);
	iXtX = MatrixInverse(XtX,NULL);
	if(iXtX == NULL) {
	  MatrixFree(&y);
	  MatrixFree(&X);
	  MatrixFree(&Xt);
	  MatrixFree(&XtX);
	  continue;
	}
	Xty = MatrixMultiply(Xt,y,NULL);
	beta = MatrixMultiply(iXtX,Xty,beta);
	for(tt=0; tt < pvf->nframes; tt++)
	  MRIsetVoxVal(pvc, c,r,s,tt, beta->rptr[tt+1][1]);

	if(0){
	  printf("c=%d; r=%d; s=%d;\n",c+1,r+1,s+1);
	  MatrixPrint(stdout,y);
	  MatrixPrint(stdout,X);
	  MatrixPrint(stdout,beta);
	  //exit(1);
	}

	MatrixFree(&y);	
	MatrixFree(&X);	
	MatrixFree(&Xt);
	MatrixFree(&XtX);
	MatrixFree(&Xty);
	MatrixFree(&iXtX);
	nhits++;
      } //s
    } // r
  } // c
  printf("\n");
  printf("nhits = %d\n",nhits);

  MatrixFree(&beta);

  return(pvc);
}

