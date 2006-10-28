// $Id: dti.c,v 1.11 2006/10/28 18:24:02 greve Exp $

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pwd.h>
#include <time.h>
#include "fsenv.h"
#include "utils.h"
#include "version.h"
#include "fio.h"
#include "mri.h"
#include "DICOMRead.h"
#include "dti.h"

/* --------------------------------------------- */
// Return the CVS version of this file.
const char *DTIsrcVersion(void) { 
  return("$Id: dti.c,v 1.11 2006/10/28 18:24:02 greve Exp $");
}


/*-----------------------------------------------------------------
  DTIparamsFromSiemensAscii() - reads in diffusion parameters from the
  Siemens ASCII header stored in fname. fname may be a siemens dicom
  file or an infodump file as produced by mri_probedicom run on a
  siemens dicom file. 
  -----------------------------------------------------------------*/
int DTIparamsFromSiemensAscii(char *fname, float *bValue, 
			      int *nAcq, int *nDir, int *nB0)
{
  char *tag, *pc;

  if(!fio_FileExistsReadable(fname)){
    printf("ERROR: cannot read %s\n",fname);
    return(1);
  }

  tag = "sDiffusion.alBValue[1]";
  pc = SiemensAsciiTag(fname,tag);
  if(pc == NULL){
    printf("ERROR: cannot extract %s from %s\n",tag,fname);
    return(1);
  }
  sscanf(pc, "%f", bValue);
  printf("bValue = %g\n",*bValue);
  free(pc);

  tag = "sWiPMemBlock.alFree[8]";
  pc = SiemensAsciiTag(fname,tag);
  if(pc == NULL){
    printf("ERROR: cannot extract %s from %s\n",tag,fname);
    return(1);
  }
  sscanf(pc, "%d", nB0);
  printf("nB0 = %d\n",*nB0);
  free(pc);

  tag = "sDiffusion.lDiffDirections";
  pc = SiemensAsciiTag(fname,tag);
  if(pc == NULL){
    printf("ERROR: cannot extract %s from %s\n",tag,fname);
    return(1);
  }
  sscanf(pc, "%d", nDir);
  printf("nDir = %d\n",*nDir);
  free(pc);

  *nAcq = 1;
  printf("nAcq = %d (HARDWIRED)\n",*nAcq);

  return(0);
}
/*------------------------------------------------------------*/
int DTIloadGradients(DTI *dti, char *GradFile)
{
  static char tmpstr[2000];
  FILE *fp;
  int c,r,n;
  FSENV *fsenv;

  fsenv = FSENVgetenv();

  if(GradFile) dti->GradFile = strcpyalloc(GradFile);
  if(dti->GradFile == NULL){
    //sprintf(tmpstr,"%s/diffusion/graddir/gradient_mgh_dti%02d.gdt",
    //fsenv->FREESURFER_HOME,dti->nDir);
    sprintf(tmpstr,"%s/gradient_mgh_dti%02d.gdt",
	    getenv("FS_DIFF_GRAD_DIR"),dti->nDir);
    dti->GradFile = strcpyalloc(tmpstr);
    printf("GradFile %s\n",dti->GradFile);
  }
  
  fp = fopen(dti->GradFile,"r");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",dti->GradFile);
    return(1);
  }

  dti->GradDir = MatrixAlloc(dti->nDir+dti->nB0,3,MATRIX_REAL);

  // Set the first nB0 rows to be all 0s (no gradients)
  for(r=1; r < dti->nB0+1; r++){
    dti->GradDir->rptr[r][1] = 0;
    dti->GradDir->rptr[r][2] = 0;
    dti->GradDir->rptr[r][3] = 0;
  }

  for(r=dti->nB0+1; r <= dti->GradDir->rows; r++){
    for(c=1; c <= 3; c++){
      n = fscanf(fp,"%f",&(dti->GradDir->rptr[r][c]));
      if(n != 1){
	printf("ERROR: reading gradients from %s\n",dti->GradFile);
	printf("  row = %d, col = %d\n",r-1,c);
	fclose(fp);
	return(1);
      }
    }
  }
  fclose(fp);

  FSENVfree(&fsenv);
  return(0);
}

/*--------------------------------------------------------*/
DTI *DTIstructFromSiemensAscii(char *fname)
{
  int err;
  DTI *dti;

  dti = (DTI *) calloc(sizeof(DTI),1);

  err = DTIparamsFromSiemensAscii(fname, &dti->bValue, &dti->nAcq, 
				  &dti->nDir, &dti->nB0);
  if(err){
    free(dti);
    dti = NULL;
    return(NULL);
  }
  err = DTIloadGradients(dti, NULL);
  if(err){
    free(dti);
    dti = NULL;
    return(NULL);
  }
  
  DTInormGradDir(dti);
  DTIdesignMatrix(dti);

  return(dti);
}

/*--------------------------------------------------------*/
int DTInormGradDir(DTI *dti)
{
  int r,c;
  double len, maxlen;

  dti->GradDirNorm = MatrixAlloc(dti->GradDir->rows,dti->GradDir->cols,MATRIX_REAL);

  maxlen = 0;
  for(r=1; r <= dti->GradDir->rows; r++){
    len = 0;
    for(c=1; c <= dti->GradDir->cols; c++) 
      len += pow(dti->GradDir->rptr[r][c],2.0);
    len = sqrt(len);
    if(maxlen < len) maxlen = len;
  }

  for(r=1; r <= dti->GradDir->rows; r++){
    len = 0;
    for(c=1; c <= dti->GradDir->cols; c++) 
      dti->GradDirNorm->rptr[r][c] = dti->GradDir->rptr[r][c]/maxlen;
  }

  return(0);
}
/*--------------------------------------------------------*/
int DTIdesignMatrix(DTI *dti)
{
  int r,xr,nthacq;
  MATRIX *g;

  g = dti->GradDirNorm;
  dti->B = MatrixAlloc((g->rows)*dti->nAcq,7,MATRIX_REAL);

  xr = 1;
  for(nthacq=0; nthacq < dti->nAcq; nthacq++){
    for(r=1; r <= g->rows; r ++){

      dti->B->rptr[xr][1] = dti->bValue * pow(g->rptr[r][1],2.0);
      dti->B->rptr[xr][2] = 2 * dti->bValue * g->rptr[r][1]*g->rptr[r][2];
      dti->B->rptr[xr][3] = 2 * dti->bValue * g->rptr[r][1]*g->rptr[r][3];
      
      dti->B->rptr[xr][4] = dti->bValue * pow(g->rptr[r][2],2.0);
      dti->B->rptr[xr][5] = 2 * dti->bValue * g->rptr[r][2]*g->rptr[r][3];
      
      dti->B->rptr[xr][6] = dti->bValue * pow(g->rptr[r][3],2.0);
      
      dti->B->rptr[xr][7] = 1;
      xr++;
    }

  }
  //MatrixWriteTxt("G.dat",dti->GradDirNorm);
  //MatrixWriteTxt("B.dat",dti->B);

  return(0);
}


/*---------------------------------------------------------*/
MRI *DTIbeta2Tensor(MRI *beta, MRI *mask, MRI *tensor)
{
  int c,r,s;
  double m,v;

  if(beta->nframes < 6){
    printf("ERROR: beta must have at least 6 frames\n");
    return(NULL);
  }
  if(tensor == NULL){
    tensor = MRIcloneBySpace(beta, MRI_FLOAT, 9); // 9 = 3x3
    if(!tensor) return(NULL);
  }
  // should check consistency with spatial

  for(c=0; c < beta->width; c++){
    for(r=0; r < beta->height; r++){
      for(s=0; s < beta->depth; s++){
	if(mask){
	  m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m < 0.5) continue;
	}

	// 0 1 2 --> 0 1 2
	// 1 3 4 --> 3 4 5
	// 2 4 5 --> 6 7 8

	// 0 -> 0
	v = MRIgetVoxVal(beta,c,r,s,0);
	MRIsetVoxVal(tensor,c,r,s,0,v);

	// 1 -> 1, 3
	v = MRIgetVoxVal(beta,c,r,s,1);
	MRIsetVoxVal(tensor,c,r,s,1,v);
	MRIsetVoxVal(tensor,c,r,s,3,v);

	// 2 -> 2, 6
	v = MRIgetVoxVal(beta,c,r,s,2);
	MRIsetVoxVal(tensor,c,r,s,2,v);
	MRIsetVoxVal(tensor,c,r,s,6,v);

	// 3 -> 4
	v = MRIgetVoxVal(beta,c,r,s,3);
	MRIsetVoxVal(tensor,c,r,s,4,v);

	// 4 -> 5, 7
	v = MRIgetVoxVal(beta,c,r,s,4);
	MRIsetVoxVal(tensor,c,r,s,5,v);
	MRIsetVoxVal(tensor,c,r,s,7,v);

	// 5 -> 8
	v = MRIgetVoxVal(beta,c,r,s,5);
	MRIsetVoxVal(tensor,c,r,s,8,v);

      }
    }
  }

  return(tensor);
}
/*---------------------------------------------------------*/
int DTItensor2Eig(MRI *tensor, MRI *mask,   MRI **evals, 
		  MRI **evec1, MRI **evec2, MRI **evec3)
{
  int c,r,s,a,b,n;
  double m;
  MATRIX *T,*Evec;
  float eval[3];

  if(tensor->nframes != 9){
    printf("ERROR: tensor must have 9 frames\n");
    return(1);
  }
  if(*evals == NULL){
    *evals = MRIcloneBySpace(tensor, MRI_FLOAT,3);
    if(!*evals) return(1);
  }
  if(*evec1 == NULL){
    *evec1 = MRIcloneBySpace(tensor, MRI_FLOAT,3);
    if(!*evec1) return(1);
  }
  if(*evec2 == NULL){
    *evec2 = MRIcloneBySpace(tensor, MRI_FLOAT,3);
    if(!*evec2) return(1);
  }
  if(*evec3 == NULL){
    *evec3 = MRIcloneBySpace(tensor, MRI_FLOAT,3);
    if(!*evec3) return(1);
  }
  // should check consistency with spatial

  T = MatrixAlloc(3,3,MATRIX_REAL);
  Evec = MatrixAlloc(3, 3, MATRIX_REAL);

  for(c=0; c < tensor->width; c++){
    for(r=0; r < tensor->height; r++){
      for(s=0; s < tensor->depth; s++){
	if(mask){
	  m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m < 0.5) continue;
	}

	// Load up the tensor into a matrix struct
	// 0 1 2
	// 3 4 5
	// 6 7 8
	n = 0;
	for(a=1; a <= 3; a++){
	  for(b=1; b <= 3; b++){
	    T->rptr[a][b] = MRIgetVoxVal(tensor,c,r,s,n);
	    n++;
	  }
	}

	/* Do eigen-decomposition */
	MatrixEigenSystem(T, eval, Evec);
	DTIsortEV(eval, Evec);

	for(a=0; a<3; a++){
	  MRIsetVoxVal(*evals,c,r,s,a,eval[a]);
	  MRIsetVoxVal(*evec1,c,r,s,a,Evec->rptr[a+1][1]);
	  MRIsetVoxVal(*evec2,c,r,s,a,Evec->rptr[a+1][2]);
	  MRIsetVoxVal(*evec3,c,r,s,a,Evec->rptr[a+1][3]);
	}

      } // slice
    } // row
  } // col

  MatrixFree(&T);
  MatrixFree(&Evec);

  return(0);
}

/*---------------------------------------------------------
  DTIsortEV() - sorts the eigenvalues and eigenvectors from
  max to min.
  ---------------------------------------------------------*/
int DTIsortEV(float *EigVals, MATRIX *EigVecs)
{
  int r;
  static MATRIX *EigVecsTmp=NULL;
  static float EigValsTmp[3];

  for(r=0; r<3; r++) EigValsTmp[r] = EigVals[r];
  EigVecsTmp = MatrixCopy(EigVecs,EigVecsTmp);

  if(EigVals[0] > EigVals[1] && EigVals[0] > EigVals[2]){
    // 1st is max
    if(EigVals[1] > EigVals[2]){
      // 1st > 2nd > 3rd -- nothing to do
      return(0);
    } else {
      // 1st > 3rd > 2nd -- swap 2nd and 3rd cols
      for(r=1; r<=3; r++){
	EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][3];
	EigVecs->rptr[r][3] = EigVecsTmp->rptr[r][2];
      }
      EigVals[2-1] = EigValsTmp[3-1];
      EigVals[3-1] = EigValsTmp[2-1];
      return(0);
    }
  }

  if(EigVals[1] > EigVals[0] && EigVals[1] > EigVals[2]){
    // 2nd is max
    if(EigVals[0] > EigVals[2]){
      // 2nd > 1st > 3rd -- swap 1st and 2nd
      for(r=1; r<=3; r++){
	EigVecs->rptr[r][1] = EigVecsTmp->rptr[r][2];
	EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][1];
      }
      EigVals[1-1] = EigValsTmp[2-1];
      EigVals[2-1] = EigValsTmp[1-1];
      return(0);
    } else {
      // 2nd > 3rd > 1st 
      for(r=1; r<=3; r++){
	EigVecs->rptr[r][1] = EigVecsTmp->rptr[r][2];
	EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][3];
	EigVecs->rptr[r][3] = EigVecsTmp->rptr[r][1];
      }
      EigVals[1-1] = EigValsTmp[2-1];
      EigVals[2-1] = EigValsTmp[3-1];
      EigVals[3-1] = EigValsTmp[1-1];
      return(0);
    }
  }

  // 3rd is max if it gets here
  if(EigVals[0] > EigVals[1]){
    // 3rd > 1st > 2nd 
    for(r=1; r<=3; r++){
      EigVecs->rptr[r][1] = EigVecsTmp->rptr[r][3];
      EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][1];
      EigVecs->rptr[r][3] = EigVecsTmp->rptr[r][2];
    }
    EigVals[1-1] = EigValsTmp[3-1];
    EigVals[2-1] = EigValsTmp[1-1];
    EigVals[3-1] = EigValsTmp[2-1];
    return(0);
  } else {
    // 3rd > 2nd > 1st 
    for(r=1; r<=3; r++){
      EigVecs->rptr[r][1] = EigVecsTmp->rptr[r][3];
      EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][2];
      EigVecs->rptr[r][3] = EigVecsTmp->rptr[r][1];
    }
    EigVals[1-1] = EigValsTmp[3-1];
    EigVals[2-1] = EigValsTmp[2-1];
    EigVals[3-1] = EigValsTmp[1-1];
    return(0);
  }

  printf("DTIsortEV(): ERROR: should never get here\n");
  for(r=1; r<=3; r++) printf("%g ",EigValsTmp[r]);
  printf("\n");

  return(1);
}

/*---------------------------------------------------------
  DTIbeta2LowB() - computes exp(-v) where v is the offset
  regression parameter. This should be the average volume
  when bvalue=0.
  ---------------------------------------------------------*/
MRI *DTIbeta2LowB(MRI *beta, MRI *mask, MRI *lowb)
{
  int c,r,s;
  double m,v;

  if(beta->nframes < 7){
    printf("ERROR: beta must have at least 7 frames\n");
    return(NULL);
  }
  if(lowb == NULL){
    lowb = MRIcloneBySpace(beta, MRI_FLOAT, 1);
    if(!lowb) return(NULL);
  }
  // should check consistency with spatial

  for(c=0; c < beta->width; c++){
    for(r=0; r < beta->height; r++){
      for(s=0; s < beta->depth; s++){
	if(mask){
	  m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m < 0.5) continue;
	}
	v = MRIgetVoxVal(beta,c,r,s,6);
	MRIsetVoxVal(lowb,c,r,s,0,exp(-v));
      }
    }
  }

  return(lowb);
}
/*---------------------------------------------------------
  DTItensor2ADC() - computes apparent diffusion coefficient
  as the trace/3.
  ---------------------------------------------------------*/
MRI *DTItensor2ADC(MRI *tensor, MRI *mask, MRI *adc)
{
  int c,r,s;
  double m,v1,v2,v3,vadc;

  if(tensor->nframes != 9){
    printf("ERROR: tensor must have at least 9 frames\n");
    return(NULL);
  }
  if(adc == NULL){
    adc = MRIcloneBySpace(tensor, MRI_FLOAT, 1);
    if(!adc) return(NULL);
  }
  // should check consistency with spatial

  for(c=0; c < tensor->width; c++){
    for(r=0; r < tensor->height; r++){
      for(s=0; s < tensor->depth; s++){
	if(mask){
	  m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m < 0.5) continue;
	}
	v1 = MRIgetVoxVal(tensor,c,r,s,0);
	v2 = MRIgetVoxVal(tensor,c,r,s,4);
	v3 = MRIgetVoxVal(tensor,c,r,s,8);
	vadc = (v1+v2+v3)/3;
	MRIsetVoxVal(adc,c,r,s,0,vadc);
      }
    }
  }

  return(adc);
}
/*------------------------------------------------------------*/
MRI *DTIeigvals2FA(MRI *evals, MRI *mask, MRI *FA)
{
  int c,r,s;
  double m,v1,v2,v3,vmean,vsse,vnorm,v;

  if(evals->nframes != 3){
    printf("ERROR: evals must have 3 frames\n");
    return(NULL);
  }
  if(FA == NULL){
    FA = MRIcloneBySpace(evals, MRI_FLOAT, 1);
    if(!FA) return(NULL);
  }
  // should check consistency with spatial

  for(c=0; c < evals->width; c++){
    for(r=0; r < evals->height; r++){
      for(s=0; s < evals->depth; s++){
	if(mask){
	  m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m < 0.5) continue;
	}
	v1 = MRIgetVoxVal(evals,c,r,s,0);
	v2 = MRIgetVoxVal(evals,c,r,s,1);
	v3 = MRIgetVoxVal(evals,c,r,s,2);
	vmean = (v1+v2+v3)/3.0;
	vsse  = pow(v1-vmean,2.0) + pow(v2-vmean,2.0) + pow(v3-vmean,2.0);
	vnorm = pow(v1,2.0) + pow(v2,2.0) + pow(v3,2.0);
	v = sqrt(1.5*vsse/vnorm); // correct formula?
	MRIsetVoxVal(FA,c,r,s,0,v);
      }
    }
  }

  return(FA);
}
/*------------------------------------------------------------
  DTIeigvals2RA() - relative anisotropy
  ------------------------------------------------------------*/
MRI *DTIeigvals2RA(MRI *evals, MRI *mask, MRI *RA)
{
  int c,r,s;
  double m,v1,v2,v3,vmean,vsse,v;

  if(evals->nframes != 3){
    printf("ERROR: evals must have 3 frames\n");
    return(NULL);
  }
  if(RA == NULL){
    RA = MRIcloneBySpace(evals, MRI_FLOAT, 1);
    if(!RA) return(NULL);
  }
  // should check consistency with spatial

  for(c=0; c < evals->width; c++){
    for(r=0; r < evals->height; r++){
      for(s=0; s < evals->depth; s++){
	if(mask){
	  m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m < 0.5) continue;
	}
	v1 = MRIgetVoxVal(evals,c,r,s,0);
	v2 = MRIgetVoxVal(evals,c,r,s,1);
	v3 = MRIgetVoxVal(evals,c,r,s,2);
	vmean = (v1+v2+v3)/3.0;
	if(vmean != 0){
	  vsse  = pow(v1-vmean,2.0) + pow(v2-vmean,2.0) + pow(v3-vmean,2.0);
	  v = sqrt(vsse/(3.0*vmean));
	}
	else v = 0;
	MRIsetVoxVal(RA,c,r,s,0,v);
      }
    }
  }

  return(RA);
}
/*------------------------------------------------------------
  DTIeigvals2VR() - volume ratio measure of anisotropy. Actually,
  1-VR is used so that it increases with anisotropy.
  ------------------------------------------------------------*/
MRI *DTIeigvals2VR(MRI *evals, MRI *mask, MRI *VR)
{
  int c,r,s;
  double m,v1,v2,v3,vmean,v;

  if(evals->nframes != 3){
    printf("ERROR: evals must have 3 frames\n");
    return(NULL);
  }
  if(VR == NULL){
    VR = MRIcloneBySpace(evals, MRI_FLOAT, 1);
    if(!VR) return(NULL);
  }
  // should check consistency with spatial

  for(c=0; c < evals->width; c++){
    for(r=0; r < evals->height; r++){
      for(s=0; s < evals->depth; s++){
	if(mask){
	  m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m < 0.5) continue;
	}
	v1 = MRIgetVoxVal(evals,c,r,s,0);
	v2 = MRIgetVoxVal(evals,c,r,s,1);
	v3 = MRIgetVoxVal(evals,c,r,s,2);
	vmean = (v1+v2+v3)/3.0;
	if(vmean != 0)	v = 1-(v1*v2*v3)/pow(vmean,3.0);
	else            v = 0.0;
	MRIsetVoxVal(VR,c,r,s,0,v);
      }
    }
  }

  return(VR);
}
/*----------------------------------------------------------------
  DTIfslBValFile() -- saves bvalues in a format that can be
  read in by FSL's dtifit with -b option. They put all the bvalues
  on one line.
  ----------------------------------------------------------------*/
int DTIfslBValFile(DTI *dti, char *bvalfname)
{
  FILE *fp;
  int n;

  fp = fopen(bvalfname,"w");
  if(!fp){
    printf("ERROR: opening %s for writing\n",bvalfname);
    return(1);
  }

  for(n=0; n < dti->nB0; n++)  fprintf(fp,"0 ");
  for(n=0; n < dti->nDir; n++) fprintf(fp,"%f ",dti->bValue);
  fprintf(fp,"\n");
  fclose(fp);

  return(0);
}
/*----------------------------------------------------------------
  DTIfslBVecFile() -- saves directions in a format that can be read in
  by FSL's dtifit with -r option. They put all the gradients on three
  lines (ie, there are 3 rows and nsamples columns).
  ----------------------------------------------------------------*/
int DTIfslBVecFile(DTI *dti, char *bvecfname)
{
  FILE *fp;
  int n,c;

  fp = fopen(bvecfname,"w");
  if(!fp){
    printf("ERROR: opening %s for writing\n",bvecfname);
    return(1);
  }

  for(c=0; c < 3; c++){
    //for(n=0; n < dti->nB0; n++)  fprintf(fp,"0 ");
    for(n=0; n < dti->GradDir->rows; n++) 
      fprintf(fp,"%f ",dti->GradDir->rptr[n+1][c+1]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  return(0);
}
