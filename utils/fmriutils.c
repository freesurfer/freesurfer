/* 
   fmriutils.c 
   $Id: fmriutils.c,v 1.5 2004/03/08 22:05:29 greve Exp $
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "mri.h"
#include "MRIio_old.h"
#include "sig.h"
#include "fmriutils.h"

#ifdef X
#undef X
#endif

/*--------------------------------------------------------*/
MRI *fMRImatrixMultiply(MRI *inmri, MATRIX *M, MRI *outmri)
{
  int c, r, s, fin, fout;
  int nframesout;
  float val;
  int nin, nout;
  int nin0, nout0;
  float *pin=NULL, *pout=NULL;

  if(inmri->nframes != M->cols){
    printf("ERROR: fMRImatrixMultiply: input dimension mismatch\n");
    return(NULL);
  }

  nframesout = M->rows;
  if(outmri==NULL){
    outmri = MRIallocSequence(inmri->width, inmri->height, inmri->depth, 
			      MRI_FLOAT, nframesout);
    if(outmri==NULL){
      printf("ERROR: fMRImatrixMultiply: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(inmri,outmri);
  }
  else{
    if(outmri->width  != inmri->width || 
       outmri->height != inmri->height || 
       outmri->depth  != inmri->depth || 
       outmri->nframes != nframesout){
      printf("ERROR: fMRImatrixMultiply: output dimension mismatch\n");
      return(NULL);
    }
    if(outmri->type != MRI_FLOAT){
      printf("ERROR: fMRImatrixMultiply: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  MRIclear(outmri);
  for(fout=0; fout < outmri->nframes; fout++){
    nout0 = fout*outmri->depth;
    for(fin=0; fin < inmri->nframes; fin++){
      val = M->rptr[fout+1][fin+1];
      nin0 = fin*inmri->depth;
      for(s=0; s < outmri->depth; s++){
	nin  = s + nin0;
	nout = s + nout0;
	for(r=0; r < outmri->height; r++){
	  pin = (float*)inmri->slices[nin][r];
	  pout = (float*)outmri->slices[nout][r];
	  for(c=0; c < outmri->width; c++)
	    (*pout++) += val*(*pin++);
	}
      }
    }
  }

  return(outmri);
}
/*--------------------------------------------------------*/
MRI *fMRIvariance(MRI *fmri, float DOF, int RmMean, MRI *var)
{
  int c, r, s, f;
  float val,sumsqval, sumval;

  if(DOF < 0) DOF = fmri->nframes;

  if(var==NULL){
    var = MRIallocSequence(fmri->width, fmri->height, fmri->depth, 
			      MRI_FLOAT, 1);
    if(var==NULL){
      printf("ERROR: fMRIvariance: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(fmri,var);
  }
  else{
    if(var->width  != fmri->width || 
       var->height != fmri->height || 
       var->depth  != fmri->depth){
      printf("ERROR: fMRIvariance: output dimension mismatch\n");
      return(NULL);
    }
    if(var->type != MRI_FLOAT){
      printf("ERROR: fMRIvariance: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  for(c=0; c < fmri->width; c++){
    for(r=0; r < fmri->height; r++){
      for(s=0; s < fmri->depth; s++){
	sumval = 0;
	sumsqval = 0;
	for(f=0; f < fmri->nframes; f++){
	  val = MRIgetVoxVal(fmri, c, r, s, f);
	  sumsqval += (val*val);
	  if(RmMean) sumval += val;
	}
	MRIFseq_vox(var,c,r,s,0) = sumsqval/DOF;
	if(RmMean)
	  MRIFseq_vox(var,c,r,s,0) -= ((sumval/DOF)*(sumval/DOF));
      }
    }
  }

  return(var);
}
/*--------------------------------------------------------
  fMRIsumSquare() - computes the sum of the squares over the
  frames. If the Update flag is set, then the sum of the 
  squares is added to that already in sumsqr. If sumsqr
  is NULL, it will be allocated.
  --------------------------------------------------------*/
MRI *fMRIsumSquare(MRI *fmri, int Update, MRI *sumsqr)
{
  int c, r, s, f, n;
  float v;
  float *pfmri=NULL, *psumsqr = NULL;

  if(sumsqr==NULL){
    sumsqr = MRIallocSequence(fmri->width, fmri->height, fmri->depth, 
			      MRI_FLOAT, 1);
    if(sumsqr==NULL){
      printf("ERROR: fMRIsumSquare: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(fmri,sumsqr);
  }
  else{
    if(sumsqr->width  != fmri->width || 
       sumsqr->height != fmri->height || 
       sumsqr->depth  != fmri->depth){
      printf("ERROR: fMRIsumsqriance: output dimension mismatch\n");
      return(NULL);
    }
    if(sumsqr->type != MRI_FLOAT){
      printf("ERROR: fMRIsumsqriance: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  if(Update) MRIclear(sumsqr);
  n = 0;
  for(f=0; f < fmri->nframes; f++){
    for(s=0; s < fmri->depth; s++){
      for(r=0; r < fmri->height; r++){
	pfmri   = (float *)   fmri->slices[n][r];
	psumsqr = (float *) sumsqr->slices[s][r];
	for(c=0; c < fmri->width; c++){
	  v = (*pfmri++);
	  (*psumsqr++) += (v*v);
	}
      }
      n++;
    }
  }

  return(sumsqr);
}
/*--------------------------------------------------------------------*/
MRI *fMRIcomputeT(MRI *ces, MATRIX *X, MATRIX *C, MRI *var, MRI *t)
{
  int c, r, s;
  MATRIX *Xt, *XtX, *iXtX, *CiXtX, *Ct, *CiXtXCt;
  float srf, cesval, std;

  if(C->rows != 1){
    printf("ERROR: fMRIcomputeT: contrast matrix has more than 1 row.\n");
    return(NULL);
  }

  if(ces->nframes != 1){
    printf("ERROR: fMRIcomputeT: contrast effect size has "
	   "more than 1 frame.\n");
    return(NULL);
  }

  if(t==NULL){
    t = MRIallocSequence(ces->width, ces->height, ces->depth, 
			      MRI_FLOAT, 1);
    if(t==NULL){
      printf("ERROR: fMRIcomputeT: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(ces,t);
  }
  else{
    if(t->width  != ces->width || 
       t->height != ces->height || 
       t->depth  != ces->depth){
      printf("ERROR: fMRIcomputeT: output dimension mismatch\n");
      return(NULL);
    }
    if(t->type != MRI_FLOAT){
      printf("ERROR: fMRIcomputeT: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  Xt      = MatrixTranspose(X,NULL);
  XtX     = MatrixMultiply(Xt,X,NULL);
  iXtX    = MatrixInverse(XtX,NULL);
  CiXtX   = MatrixMultiply(C,iXtX,NULL);
  Ct      = MatrixTranspose(C,NULL);
  CiXtXCt = MatrixMultiply(CiXtX,Ct,NULL);
  srf = sqrt(CiXtXCt->rptr[1][1]);
  //printf("fMRIcomputeT: srf = %g\n",srf);

  for(c=0; c < ces->width; c++){
    for(r=0; r < ces->height; r++){
      for(s=0; s < ces->depth; s++){
	std = sqrt(MRIgetVoxVal(var, c, r, s, 0));
	if(std == 0)
	  MRIFseq_vox(t,c,r,s,0) = 0;
	else{
	  cesval = MRIgetVoxVal(ces, c, r, s, 0);
	  MRIFseq_vox(t,c,r,s,0) = cesval/(srf*std);
	}
      }
    }
  }

  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&CiXtX);
  MatrixFree(&Ct);
  MatrixFree(&CiXtXCt);

  return(t);
}

/*--------------------------------------------------------*/
MRI *fMRIsigT(MRI *t, float DOF, MRI *sig)
{
  int c, r, s, f;
  float tval, sigtval;

  if(sig==NULL){
    sig = MRIallocSequence(t->width, t->height, t->depth, 
			   MRI_FLOAT, t->nframes);
    if(sig==NULL){
      printf("ERROR: fMRIsigT: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(t,sig);
  }
  else{
    if(t->width   != sig->width  || 
       t->height  != sig->height || 
       t->depth   != sig->depth  ||
       t->nframes != sig->nframes){
      printf("ERROR: fMRIsigT: output dimension mismatch\n");
      return(NULL);
    }
    if(sig->type != MRI_FLOAT){
      printf("ERROR: fMRIsigT: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  for(c=0; c < t->width; c++){
    for(r=0; r < t->height; r++){
      for(s=0; s < t->depth; s++){
	for(f=0; f < t->nframes; f++){
	  tval = MRIFseq_vox(t,c,r,s,f);
          sigtval = sigt(tval, rint(DOF));
	  if(tval < 0) sigtval *= -1;
	  MRIFseq_vox(sig,c,r,s,f) = sigtval;
	}
      }
    }
  }

  return(sig);
}
/*--------------------------------------------------------*/
MRI *fMRInskip(MRI *inmri, int nskip, MRI *outmri)
{
  int c, r, s, fin, fout;
  int nframesout;
  float val;

  if(inmri->nframes <= nskip){
    printf("ERROR: fMRInskip: nskip >= nframes\n");
    return(NULL);
  }

  nframesout = inmri->nframes - nskip;
  if(outmri==NULL){
    outmri = MRIallocSequence(inmri->width, inmri->height, inmri->depth, 
			      inmri->type, nframesout);
    if(outmri==NULL){
      printf("ERROR: fMRInskip: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(inmri,outmri);
  }
  else{
    if(outmri->width  != inmri->width || 
       outmri->height != inmri->height || 
       outmri->depth  != inmri->depth || 
       outmri->nframes != nframesout){
      printf("ERROR: fMRInskip: output dimension mismatch\n");
      return(NULL);
    }
    if(outmri->type != inmri->type){
      printf("ERROR: fMRInskip: structure type mismatch\n");
      return(NULL);
    }
  }

  MRIclear(outmri);
  for(fout=0; fout < outmri->nframes; fout++){
    fin = fout + nskip;
    for(s=0; s < outmri->depth; s++){
      for(r=0; r < outmri->height; r++){
	for(c=0; c < outmri->width; c++){
	  val = MRIgetVoxVal(inmri, c, r, s, fin);
	  MRIsetVoxVal(outmri,c, r, s, fout, val);
	}
      }
    }
  }

  return(outmri);
}
