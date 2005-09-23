// fsglm.c - routines to perform GLM analysis.
// $Id: fsglm.c,v 1.8 2005/09/23 03:27:40 greve Exp $

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <float.h>

#include "utils.h"
#include "fsglm.h"
#include "timer.h"
#include "gsl/gsl_cdf.h"

/* --------------------------------------------- */
// Return the CVS version of this file.
const char *GLMSrcVersion(void) { 
  return("$Id: fsglm.c,v 1.8 2005/09/23 03:27:40 greve Exp $"); 
}

/*-------------------------------------------------------
  GLMalloc - allocs a GLMMAT struct and makes sure
  that everything is either NULL or 0. Does not alloc any
  of the matrices (that's done by GLMfit(), GLMtest(), and
  GLMtransposeC(), and GLMmatrices.
  -------------------------------------------------------*/
GLMMAT *GLMalloc(void)
{
  int n;
  GLMMAT *glm;

  glm = (GLMMAT *) calloc(sizeof(GLMMAT),1);
  glm->y = NULL;
  glm->X = NULL;
  glm->beta = NULL;
  glm->yhat = NULL;
  glm->eres = NULL;
  glm->rvar = 0;
  glm->dof  = 0;
  glm->ill_cond_flag  = 0;

  glm->Xt   = NULL;
  glm->XtX  = NULL;
  glm->iXtX = NULL;
  glm->Xty  = NULL;

  glm->ncontrasts = 0;

  for(n=0; n < GLMMAT_NCONTRASTS_MAX; n++){
    glm->C[n] = NULL;
    glm->Cname[n] = NULL;
    glm->Ccond[n] = -1;

    glm->Mpmf[n] = NULL;
    glm->ypmfflag[n] = 0;
    glm->ypmf[n] = NULL;

    glm->gamma[n] = NULL;
    glm->gCVM[n] = NULL;

    glm->F[n] = 0;
    glm->p[n] = 0;

    glm->Ct[n] = NULL;
    glm->CiXtX[n] = NULL;
    glm->CiXtXCt[n] = NULL;
    glm->igCVM[n] = NULL;

    glm->igCVM[n] = NULL;
    glm->gammat[n] = NULL;
    glm->gtigCVM[n] = NULL;
  }
  return(glm);
}

/*---------------------------------------------------------------------
  GLMfree() - frees all the matrices associcated with the GLM struct,
  and the GLM struct itself.
  ------------------------------------------------------------------*/
int GLMfree(GLMMAT **pglm)
{
  int n;
  GLMMAT *glm;
  glm = *pglm;

  if(glm->y) MatrixFree(&glm->y);
  if(glm->X) MatrixFree(&glm->X);
  if(glm->beta) MatrixFree(&glm->beta);
  if(glm->yhat) MatrixFree(&glm->yhat);
  if(glm->eres) MatrixFree(&glm->eres);

  if(glm->Xt)   MatrixFree(&glm->Xt);
  if(glm->XtX)  MatrixFree(&glm->XtX);
  if(glm->iXtX) MatrixFree(&glm->iXtX);
  if(glm->Xty)  MatrixFree(&glm->Xty);

  for(n=0; n < GLMMAT_NCONTRASTS_MAX; n++){
    if(glm->C[n])           MatrixFree(&glm->C[n]);
    if(glm->Cname[n])       free(&glm->Cname[n]);
    if(glm->Mpmf[n])        MatrixFree(&glm->Mpmf[n]);
    if(glm->ypmf[n])        MatrixFree(&glm->ypmf[n]);
    if(glm->Ct[n])          MatrixFree(&glm->Ct[n]);
    if(glm->CiXtX[n])       MatrixFree(&glm->CiXtX[n]);
    if(glm->CiXtXCt[n])     MatrixFree(&glm->CiXtXCt[n]);
    if(glm->gCVM[n])        MatrixFree(&glm->gCVM[n]);
    if(glm->igCVM[n])       MatrixFree(&glm->igCVM[n]);
    if(glm->gamma[n])       MatrixFree(&glm->gamma[n]);
    if(glm->gammat[n])      MatrixFree(&glm->gammat[n]);
    if(glm->gtigCVM[n])     MatrixFree(&glm->gtigCVM[n]);
  }
  free(*pglm);
  *pglm = NULL;
  return(0);
}

/*-----------------------------------------------------------------
  GLMtransposeC() - given all the C's computes all the Ct's.  Also
  computes condition number of each C as well as it's PMF.  This
  includes the allocation of Ct[n] and PMF. It would be possible to do
  this within GLMtest(), but GLMtest() may be run many times whereas
  Ct only needs to be computed once. 
  ----------------------------------------------------------------*/
int GLMtransposeC(GLMMAT *glm)
{
  int n;
  for(n=0; n < glm->ncontrasts; n++){
    glm->Ct[n]   = MatrixTranspose(glm->C[n],NULL);
    glm->Mpmf[n] = GLMpmfMatrix(glm->C[n],&glm->Ccond[n],NULL);
  }
  return(0);
}

/*---------------------------------------------------------------
  GLMmatrices() - compute all the matrices needed to do the
  estimation and testing, but does not do estimation or testing.
  This could be done inside of GLMfit() and/or GLMtest(). However,
  it may or may not be necessary to compute these matrices more
  than once. Eg, when doing an analysis where the desgin matrix
  is the same at all voxels, then it is not necessary.
  ---------------------------------------------------------------*/
int GLMmatrices(GLMMAT *glm)
{
  int n;
  MATRIX *Mtmp;

  glm->dof = glm->X->rows - glm->X->cols;

  glm->Xt   = MatrixTranspose(glm->X,glm->Xt);
  glm->XtX  = MatrixMultiply(glm->Xt,glm->X,glm->XtX);
  Mtmp = MatrixInverse(glm->XtX,glm->iXtX);
  if(Mtmp == NULL){
    // Ill-conditioned
    MatrixPrint(stdout,glm->X);
    glm->ill_cond_flag = 1;
    return(1);
  }
  glm->ill_cond_flag  = 0;
  glm->iXtX = Mtmp;

  for(n = 0; n < glm->ncontrasts; n++){
    // gamma = C*beta
    // gCVM  = rvar*J*C*inv(X'*X)*C'
    // F     = gamma' * inv(gCVM) * gamma;
    glm->CiXtX[n]    = MatrixMultiply(glm->C[n],glm->iXtX,glm->CiXtX[n]);
    glm->CiXtXCt[n]  = MatrixMultiply(glm->CiXtX[n],glm->Ct[n],glm->CiXtXCt[n]);
  }
  return(0);
}

/*---------------------------------------------------------------
  GLMfit() - fit linear parameters (betas). Also computes yhat,
  eres, and rvar. May want to defer rvar at some point (eg, to
  do spatial filtering on the eres). XtX and iXtX must have been
  computed by GLMmatrices() first.
  ---------------------------------------------------------------*/
int GLMfit(GLMMAT *glm)
{
  int f;

  if(glm->ill_cond_flag) return(0);

  // Compute X'*y
  glm->Xty  = MatrixMultiply(glm->Xt,glm->y,glm->Xty);

  // Now do the actual parameter (beta) estmation
  // beta = inv(X'*X)*X'*y
  glm->beta = MatrixMultiply(glm->iXtX,glm->Xty,glm->beta);

  // Compute yhat, eres, and residual variance
  glm->yhat = MatrixMultiply(glm->X, glm->beta, glm->yhat);
  glm->eres = MatrixSubtract(glm->y, glm->yhat, glm->eres);
  glm->rvar = 0;
  for(f = 1; f <= glm->eres->rows; f++)
    glm->rvar += (glm->eres->rptr[f][1] * glm->eres->rptr[f][1]);
  glm->rvar /= glm->dof;
  if(glm->rvar < FLT_MIN) glm->rvar = FLT_MIN; // not quite 0

  return(0);
}
/*------------------------------------------------------------------------
  GLMtest() - tests all the contrasts for the given GLM. Must have already
  run GLMtransposeC(), GLMmatrices(), and GLMfit().
  ------------------------------------------------------------------------*/
int GLMtest(GLMMAT *glm)
{
  int n;
  double dtmp;
  static MATRIX *F=NULL;

  if(glm->ill_cond_flag){
    // If it's ill cond, just return F=0
    for(n = 0; n < glm->ncontrasts; n++){
      glm->F[n] = 0;
      glm->p[n] = 1;
    }
    return(0);
  }

  for(n = 0; n < glm->ncontrasts; n++){
    // gamma = C*beta
    // gCVM  = rvar*J*C*inv(X'*X)*C'
    // F     = gamma' * inv(gCVM) * gamma;
    // CiXtX and CiXtXCt are now computed by GLMmatrices().
    //glm->CiXtX[n]    = MatrixMultiply(glm->C[n],glm->iXtX,glm->CiXtX[n]);
    //glm->CiXtXCt[n]  = MatrixMultiply(glm->CiXtX[n],glm->Ct[n],glm->CiXtXCt[n]);
    dtmp = glm->rvar*glm->C[n]->rows;
    glm->gamma[n]    = MatrixMultiply(glm->C[n],glm->beta,glm->gamma[n]);
    glm->gammat[n]   = MatrixTranspose(glm->gamma[n],glm->gammat[n]);
    glm->gCVM[n]     = MatrixScalarMul(glm->CiXtXCt[n],dtmp,glm->gCVM[n]);
    glm->igCVM[n]    = MatrixInverse(glm->gCVM[n],glm->igCVM[n]);
    glm->gtigCVM[n]  = MatrixMultiply(glm->gammat[n],glm->igCVM[n],glm->gtigCVM[n]);
    F                = MatrixMultiply(glm->gtigCVM[n],glm->gamma[n],F);
    glm->F[n]        = F->rptr[1][1];
    glm->p[n]        = gsl_cdf_fdist_Q(glm->F[n],glm->C[n]->rows,glm->dof);
    if(glm->ypmfflag[n])
      glm->ypmf[n] = MatrixMultiply(glm->Mpmf[n],glm->beta,glm->ypmf[n]);
  }
  return(0);
}

/*-----------------------------------------------------------
  GLMprofile() - this can be used as both a profile and
  a memory leak tester. Design matrix is nrows-by-ncols
  (which forces y to be nrows-by-1). ncon contrasts are
  tested, where each contrast matrix is 2-by-ncols. Returns
  the number of msec used.
  -----------------------------------------------------------*/
int GLMprofile(int nrows, int ncols, int ncon, int niters)
{
  int n,c,msec;
  GLMMAT *glm;
  struct timeb  then ;

  TimerStart(&then) ;
  for(n=0; n < niters; n++){
    glm = GLMalloc();
    glm->y = MatrixDRand48(nrows, 1, NULL);
    glm->X = MatrixDRand48(nrows, ncols, NULL );
    glm->ncontrasts = ncon;
    for(c=0; c < glm->ncontrasts; c++) {
      glm->C[c] = MatrixDRand48(2, ncols, NULL );
      glm->ypmfflag[c] = 1;
    }
    GLMtransposeC(glm);
    GLMmatrices(glm);
    GLMfit(glm);
    GLMtest(glm);
    GLMfree(&glm);
  }
  msec = TimerStop(&then) ;

  printf("GLMprofile: nrows=%d, ncols=%d, ncon=%d, niters=%d, msec=%d, avgmsec=%g\n",
	 nrows,ncols,ncon,niters,msec,(float)msec/niters);

  return(msec);
}
/*---------------------------------------------------------
  GLMsynth() - synthesizes y, X, and Cs and fits and tests. 
  ---------------------------------------------------------*/
GLMMAT *GLMsynth(void)
{
  static char tmpstr[1000];
  int nrows, ncols, ncon, c;
  GLMMAT *glm;
  
  nrows = 100;
  ncols = 10;
  ncon = 3;

  glm = GLMalloc();
  glm->y = MatrixDRand48(nrows, 1, NULL);
  glm->X = MatrixDRand48(nrows, ncols, NULL );
  glm->ncontrasts = ncon;
  for(c=0; c < ncon; c++){
    glm->C[c] = MatrixDRand48(c+1, ncols, NULL );
    glm->ypmfflag[c] = 1;
    sprintf(tmpstr,"contrast%02d",c);
    glm->Cname[c] = strcpyalloc(tmpstr);
  }
  GLMtransposeC(glm);
  GLMmatrices(glm);
  GLMfit(glm);
  GLMtest(glm);

  return(glm);
}
/*----------------------------------------------------------------------
  GLMresynthTest() - tests GLM by synthesizing y, X, fitting, then
  setting y=yhat (ie, resynth), refitting. rvar should be 0, but less
  than 10e-9 suffices. This is repeated niter times.  If any of the
  niters are over tolerance, then 1 is returned immediately. Otherwise
  0 is returned after all iters.  If prvar is non-null then rvar of
  the failure is passed back. If there is no failure, then the max
  rvar is passed. 
  ---------------------------------------------------------*/
int GLMresynthTest(int niters, double *prvar)
{
  int nrows, ncols, n;
  GLMMAT *glm;
  double rvarmax;
  
  nrows = 100;
  ncols = 10;
  glm = GLMalloc();

  rvarmax = 0;
  for(n=0; n<niters; n++){
    // synthesize a GLM
    glm->y = MatrixDRand48(nrows, 1, glm->y);
    glm->X = MatrixDRand48(nrows, ncols, glm->X);
    // Fit
    GLMmatrices(glm);
    GLMfit(glm);
    // Copy yhat into y
    glm->y = MatrixCopy(glm->yhat,glm->y);
    // Re-Fit
    GLMfit(glm);
    if(glm->rvar > 10e-9){
      // rvar should be 0, but at least less than 10^-9.
      // Report an error if not.
      printf("GLMresynth failure: rvar = %le\n",glm->rvar);
      if(prvar != NULL) *prvar = glm->rvar;
      GLMfree(&glm);
      return(1);
    }
    if(glm->rvar > rvarmax) rvarmax = glm->rvar;
  }
  GLMfree(&glm);
  if(prvar != NULL) *prvar = rvarmax;
  return(0);
}

/*---------------------------------------------------------
  GLMdump() - saves a lot of the stuff from the GLMMAT
  struct into ascii files in the given directory.
  ---------------------------------------------------------*/
int GLMdump(char *dumpdir, GLMMAT *glm)
{
  char fname[1000];
  FILE *fp;
  int c;

  mkdir(dumpdir,(mode_t)-1);

  sprintf(fname,"%s/y.dat",dumpdir);
  MatrixWriteTxt(fname, glm->y);

  sprintf(fname,"%s/X.dat",dumpdir);
  MatrixWriteTxt(fname, glm->X);

  sprintf(fname,"%s/dof.dat",dumpdir);
  fp = fopen(fname,"w");
  fprintf(fp,"%lf",glm->dof);
  fclose(fp);

  sprintf(fname,"%s/ill_cond_flag.dat",dumpdir);
  fp = fopen(fname,"w");
  fprintf(fp,"%d",glm->ill_cond_flag);
  fclose(fp);
  if(glm->ill_cond_flag) return(0);

  sprintf(fname,"%s/beta.dat",dumpdir);
  MatrixWriteTxt(fname, glm->beta);

  sprintf(fname,"%s/yhat.dat",dumpdir);
  MatrixWriteTxt(fname, glm->yhat);

  sprintf(fname,"%s/eres.dat",dumpdir);
  MatrixWriteTxt(fname, glm->eres);

  sprintf(fname,"%s/rvar.dat",dumpdir);
  fp = fopen(fname,"w");
  fprintf(fp,"%lf",glm->rvar);
  fclose(fp);

  sprintf(fname,"%s/ncontrasts.dat",dumpdir);
  fp = fopen(fname,"w");
  fprintf(fp,"%d",glm->ncontrasts);
  fclose(fp);
 
  for(c=0; c < glm->ncontrasts; c++){
    sprintf(fname,"%s/C%03d.dat",dumpdir,c+1);
    MatrixWriteTxt(fname, glm->C[c]);

    sprintf(fname,"%s/Ccond%03d.dat",dumpdir,c+1);
    fp = fopen(fname,"w");
    fprintf(fp,"%f",glm->Ccond[c]);
    fclose(fp);

    sprintf(fname,"%s/Mpmf%03d.dat",dumpdir,c+1);
    MatrixWriteTxt(fname, glm->Mpmf[c]);

    sprintf(fname,"%s/gamma%03d.dat",dumpdir,c+1);
    MatrixWriteTxt(fname, glm->gamma[c]);

    sprintf(fname,"%s/F%03d.dat",dumpdir,c+1);
    fp = fopen(fname,"w");
    fprintf(fp,"%lf",glm->F[c]);
    fclose(fp);

    sprintf(fname,"%s/p%03d.dat",dumpdir,c+1);
    fp = fopen(fname,"w");
    fprintf(fp,"%le",glm->p[c]);
    fclose(fp);

    if(glm->Cname[c] != NULL){
      sprintf(fname,"%s/Cname%03d.dat",dumpdir,c+1);
      fp = fopen(fname,"w");
      fprintf(fp,"%s",glm->Cname[c]);
      fclose(fp);
    }

    if(glm->ypmfflag[c]){
      sprintf(fname,"%s/ypmf%03d.dat",dumpdir,c+1);
      MatrixWriteTxt(fname, glm->ypmf[c]);
    }
  }

  return(0);
}
/*--------------------------------------------------------------------
  GLMpmfMatrix() - compute matrix P used for partial model fit (pmf)
  of a contrast. P = C'*inv(C*C')*C, where C is the contrast
  matrix. Then ypmf = X*P*beta, where ypmf is y projected into the
  contrast space. To do this, C must be well-conditioned. The
  condition number is computed and returned through cond so that the
  calling program can decide whether it is  ill-conditioned.
  ------------------------------------------------------------------*/
MATRIX *GLMpmfMatrix(MATRIX *C, double *cond, MATRIX *P)
{
  MATRIX *Ct=NULL, *CCt=NULL, *iCCt=NULL, *CtiCCt=NULL;

  if(P != NULL){
    if(P->rows != C->cols || P->cols != C->cols){
      printf("ERROR: GLMpmfMatrix: dimension mismatch\n");
      return(NULL);
    }
  }

  Ct     = MatrixTranspose(C,Ct);
  CCt    = MatrixMultiply(C,Ct,CCt);
  *cond  = MatrixConditionNumber(CCt);
  iCCt   = MatrixInverse(CCt,iCCt);
  CtiCCt = MatrixMultiply(Ct,iCCt,CtiCCt);
  P      = MatrixMultiply(CtiCCt,C,P);

  MatrixFree(&Ct);
  MatrixFree(&CCt);
  MatrixFree(&iCCt);
  MatrixFree(&CtiCCt);

  return(P);
}
