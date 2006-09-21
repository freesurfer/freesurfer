// $Id: dti.c,v 1.4 2006/09/21 04:15:08 greve Exp $

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
  return("$Id: dti.c,v 1.4 2006/09/21 04:15:08 greve Exp $");
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
    sprintf(tmpstr,"%s/diffusion/graddir/gradient_mgh_dti%02d.gdt",
	    fsenv->FREESURFER_HOME,dti->nDir);
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


