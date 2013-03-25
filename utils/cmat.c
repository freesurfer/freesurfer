/**
 * @file  cmat.c
 * @brief utilities for reading/writing a Connectome MATrix structure
 *
 * Reading and writing and utilities for the Connectome Matrix (CMAT)
 * structure.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2013/03/25 17:28:06 $
 *    $Revision: 1.2 $
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

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>

#include "cmat.h"
#include "error.h"


CMAT  *
CMATread(char *fname)
{
  CMAT   *cmat ;
  int    nlabels, i, j, ind1, ind2 ;
  FILE   *fp ;

  fp = fopen(fname, "r") ;

  fscanf(fp, "CMAT - %d\n", &nlabels) ;
  cmat = CMATalloc(nlabels, NULL) ;

  for (i = 0 ; i < nlabels ; i++)
    fscanf(fp, "%d\n", &cmat->labels[i]) ;

  for (i = 0 ; i < cmat->nlabels ; i++)
  {
    for (j = i+1 ; j < cmat->nlabels ; j++)
      fscanf(fp, "%lf", &cmat->weights[i][j]) ;
    fscanf(fp, "\n") ;
  }
  for (i = 0 ; i < cmat->nlabels ; i++)
  {
    for (j = i+1 ; j < cmat->nlabels ; j++)
    {
      fscanf(fp, "%d %d\n", &ind1, &ind2) ;
      cmat->splines[ind1][ind2] = LabelReadFrom(NULL, fp) ;
      if (feof(fp))
	break ;
    }
    if (feof(fp))
      break ;
  }
  fclose(fp) ;
  return(cmat) ;
}

int 
CMATwrite(CMAT *cmat, char *fname)
{
  FILE   *fp ;
  int    i, j ;

  fp = fopen(fname, "w") ;

  fprintf(fp, "CMAT - %d\n", cmat->nlabels) ;
  for (i = 0 ; i < cmat->nlabels ; i++)
    fprintf(fp, "%d\n", cmat->labels[i]) ;

  for (i = 0 ; i < cmat->nlabels ; i++)
  {
    for (j = i+1 ; j < cmat->nlabels ; j++)
      fprintf(fp, "%f", cmat->weights[i][j]) ;
    fprintf(fp, "\n") ;
  }

  for (i = 0 ; i < cmat->nlabels ; i++)
    for (j = i+1 ; j < cmat->nlabels ; j++)
    {
      if (cmat->splines[i][j] == NULL)
	continue ;
      fprintf(fp, "%d %d\n", i, j) ;
      LabelWriteInto(cmat->splines[i][j], fp) ;
    }
  fclose(fp) ;
  return(NO_ERROR) ;
}

CMAT *
CMATalloc(int nlabels, int *labels)
{
  CMAT *cmat ;
  int  i ;
  
  cmat = (CMAT *)calloc(1, sizeof(CMAT)) ;
  if (cmat == NULL)
    ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat", nlabels) ;

  cmat->nlabels = nlabels ;
  cmat->labels = (int *)calloc(nlabels, sizeof(int)) ;
  if (cmat->labels == NULL)
    ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->labels", nlabels) ;
  cmat->splines = (LABEL ***)calloc(nlabels, sizeof(LABEL **)) ;
  if (cmat->splines == NULL)
    ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->splines", nlabels) ;
  cmat->weights = (double **)calloc(nlabels, sizeof(double *)) ;
  if (cmat->weights == NULL)
    ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->weights", nlabels) ;
  for (i = 0 ; i < nlabels ; i++)
  {
    if (labels)
      cmat->labels[i] = labels[i] ;
    cmat->splines[i] = (LABEL **)calloc(nlabels, sizeof(LABEL *)) ;
    cmat->weights[i] = (double *)calloc(nlabels, sizeof(double)) ;
    if (cmat->weights[i] == NULL)
      ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->weights[%d]", nlabels, i) ;
    if (cmat->splines[i] == NULL)
      ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->splines[%d]", nlabels, i) ;
  }

  return(cmat) ;
}


int
CMATfree(CMAT **pcmat)
{
  CMAT *cmat ;
  int  i, j ;

  cmat = *pcmat ;
  *pcmat = NULL ;

  free(cmat->labels) ;
  for (i = 0 ; i < cmat->nlabels ; i++)
  {
    for (j = i+1 ; j < cmat->nlabels ; j++)
    {
      if (cmat->splines[i])
	LabelFree(&cmat->splines[i][j]) ;
    }
    free(cmat->splines[i]) ;
    free(cmat->weights[i]) ;
  }

  free(cmat->splines) ;
  free(cmat->weights) ;
  free(cmat) ;
  return(NO_ERROR) ;
}
