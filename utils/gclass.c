/*
 *       FILE NAME:   gclass.c
 *
 *       DESCRIPTION: Gaussian Classification
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/97
 *
*/

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "matrix.h"
#include "diag.h"
#include "error.h"
#include "const.h"
#include "gclass.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
GCLASSIFY *
GCalloc(int nclasses, int nvars, char *class_names[])
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        cno ;

  gc = (GCLASSIFY *)calloc(1, sizeof(GCLASSIFY)) ;
  if (!gc)
    ErrorReturn(NULL,
              (ERROR_NO_MEMORY,"GCalloc(%d): could not allocate GC",nclasses));

  gc->nclasses = nclasses ;
  gc->nvars = nvars ;
  gc->classes = (GCLASS *)calloc(nclasses, sizeof(GCLASS)) ;
  if (!gc->classes)
    ErrorReturn(NULL, 
                (ERROR_NO_MEMORY, 
                 "GFalloc(%d): could not allocated class table",nclasses));

  gc->log_probabilities = (float *)calloc(nclasses, sizeof(float)) ;
  if (!gc->log_probabilities)
  {
    GCfree(&gc) ;
    ErrorReturn(NULL, 
                (ERROR_NO_MEMORY, 
                 "GFalloc(%d): could not probability table",nclasses));
  }

  for (cno = 0 ; cno < nclasses ; cno++)
  {
    gcl = &gc->classes[cno] ;
    gcl->classno = cno ;
    gcl->m_covariance = MatrixAlloc(nvars,nvars,MATRIX_REAL);
    if (!gcl->m_covariance)
    {
      GCfree(&gc) ;
      ErrorReturn(NULL, 
                  (ERROR_NO_MEMORY, 
                   "GFalloc(%d): could not allocated %d x %d cov. matrix"
                   "for %dth class", nclasses, nvars, nvars, cno)) ;
    }
    gcl->m_u = MatrixAlloc(nvars,1,MATRIX_REAL);
    if (!gcl->m_u)
    {
      GCfree(&gc) ;
      ErrorReturn(NULL, 
                  (ERROR_NO_MEMORY, 
                   "GFalloc(%d): could not allocated %d x 1 mean vector"
                   "for %dth class", nclasses, nvars, cno)) ;
    }
    if (class_names)
      strncpy(gcl->class_name, class_names[cno], 29) ;
    else
      sprintf(gcl->class_name, "class %d", cno) ;
  }
  
  return(gc) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
GCtrain(GCLASSIFY *gc, int class, MATRIX *m_inputs)
{
  GCLASS  *gcl ;
  MATRIX  *m_sigma_inverse, *m_uT, *m_tmp, *m_tmp2 ;
  float   det ;

  gcl = &gc->classes[class] ;
  gcl->nobs = m_inputs->rows ;
  MatrixCovariance(m_inputs, gcl->m_covariance, gcl->m_u) ;

  m_sigma_inverse = MatrixInverse(gcl->m_covariance, NULL) ;
  m_uT = MatrixTranspose(gcl->m_u, NULL) ;
  det = MatrixDeterminant(gcl->m_covariance) ;
  gcl->m_W = MatrixScalarMul(m_sigma_inverse, -0.5f, NULL) ;
  m_tmp = MatrixMultiply(m_sigma_inverse, gcl->m_u, NULL) ;
  gcl->m_wT = MatrixTranspose(m_tmp, NULL) ;
  MatrixFree(&m_tmp) ;
  m_tmp = MatrixMultiply(m_sigma_inverse, gcl->m_u, NULL) ;
  m_tmp2 = MatrixMultiply(m_uT, m_tmp, NULL) ;
  det = MatrixDeterminant(m_sigma_inverse) ;
  gcl->w0 = -0.5*(gc->nvars * log(2*PI) + m_tmp2->rptr[1][1] + log(det)) ;

#if 0
fprintf(stdout, "\nclass %d:\n", class) ;
MatrixPrint(stdout, gcl->m_covariance) ;
MatrixPrint(stdout, m_sigma_inverse) ;
fprintf(stdout, "means: \n") ;
MatrixPrint(stdout, gcl->m_u) ;
fprintf(stdout, "det = %2.3f\n", det) ;
fprintf(stdout, "Wi = \n") ;
MatrixPrint(stdout, gcl->m_W) ;
fprintf(stdout, "w = \n") ;
MatrixPrint(stdout, gcl->m_wT) ;
fprintf(stdout, "m_tmp:\n") ;
MatrixPrint(stdout, m_tmp) ;
fprintf(stdout, "m_tmp2:\n") ;
MatrixPrint(stdout, m_tmp2) ;
fprintf(stdout, "w0: %2.3f\n", gcl->w0) ;
#endif

  MatrixFree(&m_sigma_inverse) ;
  MatrixFree(&m_uT) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&m_tmp2) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
GCfree(GCLASSIFY **pgc)
{
  GCLASSIFY *gc ;
  GCLASS    *gcl ;
  int       cno ;

  gc = *pgc ;
  *pgc = NULL ;

  if (!gc)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "GCfree: NULL pointer")) ;

  if (gc->classes) for (cno = 0 ; cno < gc->nclasses ; cno++)
  {
    gcl = &gc->classes[cno] ;
    if (gcl->m_covariance)
      MatrixFree(&gcl->m_covariance) ;
    if (gcl->m_u)
      MatrixFree(&gcl->m_u) ;
    if (gcl->m_W)
      MatrixFree(&gcl->m_W) ;
    if (gcl->m_wT)
      MatrixFree(&gcl->m_wT) ;
  }

  if (gc->log_probabilities)
    free(gc->log_probabilities) ;
  if (gc->classes)
    free(gc->classes) ;
  free(gc) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
GCclassify(GCLASSIFY *gc, MATRIX *m_x, float *prisk)
{
  int     cno, class = -1 ;
  GCLASS  *gcl ;
  MATRIX  *m_xT, *m_tmp, *m_tmp2, *m_tmp3 ;
  float   log_p, max_p ;

#if 0
fprintf(stdout, "GCclassify(%2.3f)\n", m_x->rptr[1][1]) ;
#endif

/*
   see Duda and Hart page 30
*/
  m_xT = MatrixTranspose(m_x, NULL) ;
  max_p = -100000.0f ;
  m_tmp = m_tmp2 = m_tmp3 = NULL ;
  for (cno = 0 ; cno < gc->nclasses ; cno++)
  {
    gcl = &gc->classes[cno] ;
    m_tmp = MatrixMultiply(gcl->m_W, m_x, m_tmp) ;
    m_tmp2 = MatrixMultiply(m_xT, m_tmp, m_tmp2) ;
    m_tmp3 = MatrixMultiply(gcl->m_wT, m_x, m_tmp3) ;
    log_p = gcl->w0 + m_tmp2->rptr[1][1] + m_tmp3->rptr[1][1] ;
#if 0
fprintf(stdout, "class %d: log(p) = %2.3f + %2.3f + %2.3f = %2.3f\n",
        cno, gcl->w0, m_tmp2->rptr[1][1], m_tmp3->rptr[1][1], log_p) ;
#endif
    gc->log_probabilities[cno] = log_p ;
    if (log_p > max_p)
    {
      max_p = log_p ;
      class = cno ;
    }
  }

  if (prisk)
    *prisk = exp(max_p) ;

  MatrixFree(&m_xT) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&m_tmp2) ;
  MatrixFree(&m_tmp3) ;
  return(class) ;
}

