// fsglm.c
// $Id: fsglm.c,v 1.1 2005/09/17 23:10:27 greve Exp $

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <float.h>

#include "fsglm.h"
#include "timer.h"
#include "gsl/gsl_cdf.h"

/* --------------------------------------------- */
// Return the CVS version of this file.
const char *GLMSrcVersion(void) { 
  return("$Id: fsglm.c,v 1.1 2005/09/17 23:10:27 greve Exp $"); 
}

/*-------------------------------------------------------
  GLMalloc - allocs a GLMMAT struct and makes sure
  that everything is either NULL or 0. Does not alloc any
  of the matrices (that's done by GLMfit(), GLMtest(), and
  GLMcontrastTranspose()).
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
  glm->DOF  = 0;
  glm->ill_cond_flag  = 0;

  glm->Xt   = NULL;
  glm->XtX  = NULL;
  glm->iXtX = NULL;
  glm->Xty  = NULL;

  glm->ncontrasts = 0;

  for(n=0; n < GLMMAT_NCONTRASTS_MAX; n++){
    glm->C[n] = NULL;
    glm->Ct[n] = NULL;
    glm->CiXtX[n] = NULL;
    glm->CiXtXCt[n] = NULL;
    glm->gCVM[n] = NULL;
    glm->igCVM[n] = NULL;

    glm->igCVM[n] = NULL;
    glm->gamma[n] = NULL;
    glm->gammat[n] = NULL;
    glm->gtigCVM[n] = NULL;
    glm->F[n] = 0;
    glm->p[n] = 0;
  }
  return(glm);
}

/*------------------------------------------------*/
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
    if(glm->Ct[n])          MatrixFree(&glm->Ct[n]);
    if(glm->CiXtX[n])       MatrixFree(&glm->CiXtX[n]);
    if(glm->CiXtX[n])       MatrixFree(&glm->CiXtXCt[n]);
    if(glm->gCVM[n])    MatrixFree(&glm->gCVM[n]);
    if(glm->igCVM[n])   MatrixFree(&glm->igCVM[n]);
    if(glm->gamma[n])       MatrixFree(&glm->gamma[n]);
    if(glm->gammat[n])      MatrixFree(&glm->gammat[n]);
    if(glm->gtigCVM[n]) MatrixFree(&glm->gtigCVM[n]);
  }
  free(*pglm);
  *pglm = NULL;
  return(0);
}

/*-----------------------------------------------------------------
  GLMcontrastTranspose() - given all the C's computes all the Ct's.
  This includes the allocation of Ct[n]. It would be possible
  to do this within GLMtest(), but GLMtest() may be run many times
  where as Ct only needs to be computed once.
  ----------------------------------------------------------------*/
int GLMtransposeC(GLMMAT *glm)
{
  int n;
  for(n=0; n < glm->ncontrasts; n++)
    glm->Ct[n] = MatrixTranspose(glm->C[n],NULL);
  return(0);
}


/*---------------------------------------------------------------
  GLMfit() - fit linear parameters (betas). Also computes yhat,
  eres, and rvar.
  ---------------------------------------------------------------*/
int GLMfit(GLMMAT *glm)
{
  int f;
  MATRIX *Mtmp;

  glm->DOF = glm->X->rows - glm->X->cols;

  glm->Xt   = MatrixTranspose(glm->X,glm->Xt);
  glm->XtX  = MatrixMultiply(glm->Xt,glm->X,glm->XtX);
  Mtmp = MatrixInverse(glm->XtX,glm->iXtX);
  if(Mtmp == NULL){
    // Ill-conditioned
    MatrixPrint(stdout,glm->X);
    exit(1);
    glm->ill_cond_flag = 1;
    return(1);
  }
  glm->ill_cond_flag  = 0;
  glm->iXtX = Mtmp;
  glm->Xty  = MatrixMultiply(glm->Xt,glm->y,glm->Xty);

  // Now do the actual parameter estmation
  glm->beta = MatrixMultiply(glm->iXtX,glm->Xty,glm->beta);

  // Compute residual variance
  glm->yhat = MatrixMultiply(glm->X,glm->beta,glm->yhat);
  glm->eres = MatrixSubtract(glm->y, glm->yhat, glm->eres);
  glm->rvar = 0;
  for(f = 1; f <= glm->eres->rows; f++)
    glm->rvar += (glm->eres->rptr[f][1] * glm->eres->rptr[f][1]);
  glm->rvar /= glm->DOF;
  if(glm->rvar < FLT_MIN) glm->rvar = FLT_MIN; // not quite 0

  return(0);
}

/*------------------------------------------------------------------------
  GLMtest() - tests all the contrasts for the given GLM. Must have already
  run GLMfit() and GLMcontrastTranspose().
  ------------------------------------------------------------------------*/
int GLMtest(GLMMAT *glm)
{
  int n;
  double dtmp;
  static MATRIX *F=NULL;

  if(glm->ill_cond_flag){
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
    dtmp = glm->rvar*glm->C[n]->rows;
    glm->gamma[n]    = MatrixMultiply(glm->C[n],glm->beta,glm->gamma[n]);
    glm->gammat[n]   = MatrixTranspose(glm->gamma[n],glm->gammat[n]);
    glm->CiXtX[n]    = MatrixMultiply(glm->C[n],glm->iXtX,glm->CiXtX[n]);
    glm->CiXtXCt[n]  = MatrixMultiply(glm->CiXtX[n],glm->Ct[n],glm->CiXtXCt[n]);
    glm->gCVM[n]     = MatrixScalarMul(glm->CiXtXCt[n],dtmp,glm->gCVM[n]);
    glm->igCVM[n]    = MatrixInverse(glm->gCVM[n],glm->igCVM[n]);
    glm->gtigCVM[n]  = MatrixMultiply(glm->gammat[n],glm->igCVM[n],glm->gtigCVM[n]);
    F                = MatrixMultiply(glm->gtigCVM[n],glm->gamma[n],F);
    glm->F[n]        = F->rptr[1][1];
    glm->p[n]        = gsl_cdf_fdist_Q(glm->F[n],glm->C[n]->rows,glm->DOF);
  }

  return(0);
}

/*-----------------------------------------------------------
  GLMprofile() - this can be used as both a profile and
  a memory leak tester. Design matrix is nrows-by-ncols
  (which forces y to be nrows-by-1). ncon contrasts are
  tested, where each contrast matrix is 2-by-ncols.
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
    for(c=0; c < glm->ncontrasts; c++) 
      glm->C[c] = MatrixDRand48(2, ncols, NULL );
    GLMtransposeC(glm);
    GLMfit(glm);
    GLMtest(glm);
    GLMfree(&glm);
  }
  msec = TimerStop(&then) ;

  printf("GLMprofile: nrows=%d, ncols=%d, ncon=%d, niters=%d, msec=%d, avgmsec=%g\n",
	 nrows,ncols,niters,ncon,msec,(float)msec/niters);

  return(msec);
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
  fprintf(fp,"%lf",glm->DOF);
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
  }

  return(0);
}

/*---------------------------------------------------------
  GLMsynth() - synthesizes y, X, and Cs and fits and tests. 
  ---------------------------------------------------------*/
GLMMAT *GLMsynth(void)
{
  int nrows, ncols, ncon, c;
  GLMMAT *glm;
  
  nrows = 100;
  ncols = 10;
  ncon = 3;

  glm = GLMalloc();
  glm->y = MatrixDRand48(nrows, 1, NULL);
  glm->X = MatrixDRand48(nrows, ncols, NULL );
  glm->ncontrasts = ncon;
  for(c=0; c < ncon; c++) 
    glm->C[c] = MatrixDRand48(c+1, ncols, NULL );
  GLMtransposeC(glm);
  GLMfit(glm);
  GLMtest(glm);

  return(glm);
}

/*----------------------------------------------------------------------
  GLMresynthTest() - tests GLM by synthesizing y, X, fitting, then
  setting y=yhat (resynth), refitting. rvar should be 0, but less than
  10e-9 suffices. This is repeated niter times.  If any of the niters
  are over tolerance, then 1 is returned immediately. Otherwise 0 is
  returned after all iters.  If prvar is non-null then rvar is passed
  back. 
  ---------------------------------------------------------*/
int GLMresynthTest(int niters, double *prvar)
{
  int nrows, ncols, n;
  GLMMAT *glm;
  
  nrows = 100;
  ncols = 10;
  glm = GLMalloc();

  for(n=0; n<niters; n++){

    // synthesize a GLM
    glm->y = MatrixDRand48(nrows, 1, glm->y);
    glm->X = MatrixDRand48(nrows, ncols, glm->X);
    
    // Fit
    GLMfit(glm);
    
    // Copy yhat into y
    glm->y = MatrixCopy(glm->yhat,glm->y);
    
    // Re-Fit
    GLMfit(glm);
    
    if(glm->rvar > 10e-9){
      printf("GLMresynth failure: rvar = %le\n",glm->rvar);
      if(prvar != NULL) *prvar = glm->rvar;
      GLMfree(&glm);
      return(1);
    }
  }
    
  GLMfree(&glm);

  return(0);
}

