/*
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diag.h"
#include "error.h"
#include "gclass.h"
#include "macros.h"
#include "matrix.h"
#include "proto.h"

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
GCLASSIFY *GCalloc(int nclasses, int nvars, const char *class_names[])
{
  GCLASSIFY *gc;
  GCLASS *gcl;
  int cno;

  gc = (GCLASSIFY *)calloc(1, sizeof(GCLASSIFY));
  if (!gc) ErrorReturn(NULL, (ERROR_NO_MEMORY, "GCalloc(%d): could not allocate GC", nclasses));

  gc->nclasses = nclasses;
  gc->nvars = nvars;
  gc->classes = (GCLASS *)calloc(nclasses, sizeof(GCLASS));
  if (!gc->classes) ErrorReturn(NULL, (ERROR_NO_MEMORY, "GFalloc(%d): could not allocated class table", nclasses));

  gc->log_probabilities = (float *)calloc(nclasses, sizeof(float));
  if (!gc->log_probabilities) {
    GCfree(&gc);
    ErrorReturn(NULL, (ERROR_NO_MEMORY, "GFalloc(%d): could not probability table", nclasses));
  }

  for (cno = 0; cno < nclasses; cno++) {
    gcl = &gc->classes[cno];
    gcl->classno = cno;
    gcl->m_covariance = MatrixAlloc(nvars, nvars, MATRIX_REAL);
    if (!gcl->m_covariance) {
      GCfree(&gc);
      ErrorReturn(NULL,
                  (ERROR_NO_MEMORY,
                   "GFalloc(%d): could not allocated %d x %d cov. matrix"
                   "for %dth class",
                   nclasses,
                   nvars,
                   nvars,
                   cno));
    }
    gcl->m_u = MatrixAlloc(nvars, 1, MATRIX_REAL);
    if (!gcl->m_u) {
      GCfree(&gc);
      ErrorReturn(NULL,
                  (ERROR_NO_MEMORY,
                   "GFalloc(%d): could not allocated %d x 1 mean vector"
                   "for %dth class",
                   nclasses,
                   nvars,
                   cno));
    }
    gcl->m_W = MatrixAlloc(nvars, nvars, MATRIX_REAL);
    if (!gcl->m_W) {
      GCfree(&gc);
      ErrorReturn(NULL,
                  (ERROR_NO_MEMORY,
                   "GFalloc(%d): could not allocated %d x %d matrix"
                   "for %dth class",
                   nclasses,
                   nvars,
                   nvars,
                   cno));
    }
#if 1
    gcl->m_wT = MatrixAlloc(gcl->m_u->cols, gcl->m_covariance->rows, MATRIX_REAL);
    if (!gcl->m_wT) {
      GCfree(&gc);
      ErrorReturn(NULL,
                  (ERROR_NO_MEMORY,
                   "GFalloc(%d): could not allocated %d x %d matrix"
                   "for %dth class",
                   nclasses,
                   gcl->m_u->cols,
                   gcl->m_covariance->rows,
                   cno));
    }
#endif
    if (class_names)
      strncpy(gcl->class_name, class_names[cno], 29);
    else
      sprintf(gcl->class_name, "class %d", cno);
  }

  return (gc);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
             each col of m_inputs is a vector containing
             observations of one variable. Thus, the # of rows
             is the # of observations, and the # of cols is the #
             of variables.
------------------------------------------------------*/
int GCtrain(GCLASSIFY *gc, int classnum, MATRIX *m_inputs)
{
  GCLASS *gcl;

  gcl = &gc->classes[classnum];
  gcl->nobs = m_inputs->rows;
  if (gcl->nobs == 0) /* no training data */
  {
    gcl->ill_cond = 1;
    return (NO_ERROR);
  }
#if 0
  else if (gcl->nobs <= gc->nvars)  /* not enough training data */
    gcl->m_covariance = MatrixIdentity(gc->nvars, NULL) ;
#endif
  else
    MatrixCovariance(m_inputs, gcl->m_covariance, gcl->m_u);
  return (GCinit(gc, classnum));
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int GCfree(GCLASSIFY **pgc)
{
  GCLASSIFY *gc;
  GCLASS *gcl;
  int cno;

  gc = *pgc;
  *pgc = NULL;

  if (!gc) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GCfree: NULL pointer"));

  if (gc->classes)
    for (cno = 0; cno < gc->nclasses; cno++) {
      gcl = &gc->classes[cno];
      if (gcl->m_covariance) MatrixFree(&gcl->m_covariance);
      if (gcl->m_u) MatrixFree(&gcl->m_u);
      if (gcl->m_W) MatrixFree(&gcl->m_W);
      if (gcl->m_wT) MatrixFree(&gcl->m_wT);
    }

  if (gc->log_probabilities) free(gc->log_probabilities);
  if (gc->classes) free(gc->classes);
  free(gc);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int GCclassify(GCLASSIFY *gc, MATRIX *m_x, MATRIX *m_priors, float *prisk)
{
  int cno, classnum = -1;
  GCLASS *gcl;
  static MATRIX *m_xT = NULL, *m_tmp, *m_tmp2, *m_tmp3;
  float log_p, max_p, sum_p, prior;

  if (m_x->cols != 1 || m_x->rows != gc->nvars)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GCclassify: inappropriately sized m_x"));
  /*
     see Duda and Hart page 30
  */
  if (m_xT && ((m_xT->rows != m_x->cols) || (m_xT->cols != m_x->rows))) {
    if (m_xT) MatrixFree(&m_xT);
    if (m_tmp) MatrixFree(&m_tmp);
    if (m_tmp2) MatrixFree(&m_tmp2);
    if (m_tmp3) MatrixFree(&m_tmp3);
  }
  m_xT = MatrixTranspose(m_x, m_xT);
  max_p = -100000.0f;
  classnum = -1;
  sum_p = 0.0f;
  for (cno = 0; cno < gc->nclasses; cno++) {
    gcl = &gc->classes[cno];

    /* check to see if covariance matrix was ill-conditioned */
    if (FZERO(gcl->w0) || gcl->nobs <= gc->nvars + 1) {
      gc->log_probabilities[cno] = -10000.0f;
      continue;
    }

    m_tmp = MatrixMultiply(gcl->m_W, m_x, m_tmp);
    m_tmp2 = MatrixMultiply(m_xT, m_tmp, m_tmp2);
    m_tmp3 = MatrixMultiply(gcl->m_wT, m_x, m_tmp3);
    if (m_priors)
      prior = m_priors->rptr[cno + 1][1];
    else
      prior = 1.0f;

    log_p = gcl->w0 + m_tmp2->rptr[1][1] + m_tmp3->rptr[1][1];
    if (!FZERO(prior))
      log_p += log(prior);
    else
      log_p -= 100000.0f;

    gc->log_probabilities[cno] = log_p;
    sum_p += exp(log_p);
    if (log_p > max_p) /* tentatively set this as the most probable class */
    {
      max_p = log_p;
      classnum = cno;
    }
  }

  if (prisk && classnum >= 0 && !FZERO(sum_p)) *prisk = exp(max_p) / sum_p;

  return (classnum);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          save a classifier to disk in an already opened file
------------------------------------------------------*/
int GCasciiWriteInto(FILE *fp, GCLASSIFY *gc)
{
  int classno;

  fprintf(fp, "%d %d %d\n", gc->nclasses, gc->nvars, gc->type);
  for (classno = 0; classno < gc->nclasses; classno++) GCasciiWriteClassInto(fp, &gc->classes[classno]);

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           read a classifier from disk in ASCII format from
           an already opened file.
------------------------------------------------------*/
GCLASSIFY *GCasciiReadFrom(FILE *fp, GCLASSIFY *gc)
{
  int classno, nclasses, nvars, type;

  if (fscanf(fp, "%d %d %d\n", &nclasses, &nvars, &type) != 3)
    ErrorReturn(NULL, (ERROR_BADFILE, "GCasciiReadFrom: could not scan parms"));

  if (!gc) {
    gc = GCalloc(nclasses, nvars, NULL);
    if (!gc) ErrorReturn(NULL, (ERROR_BADFILE, "GCasciiReadFrom: GCalloc failed"));
  }
  else if ((gc->nclasses != nclasses) || (gc->nvars != nvars))
    ErrorReturn(NULL, (ERROR_BADPARM, "GCasciiReadFrom: specified classifier is of wrong form"));

  for (classno = 0; classno < gc->nclasses; classno++) GCasciiReadClassFrom(fp, &gc->classes[classno]);

  return (gc);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          write a single class to disk in ascii format
------------------------------------------------------*/
int GCasciiWriteClassInto(FILE *fp, GCLASS *gcl)
{
  fprintf(fp, "%d %f %d\n", gcl->classno, gcl->w0, gcl->nobs);
  fprintf(fp, "%s\n", gcl->class_name);
  MatrixAsciiWriteInto(fp, gcl->m_covariance);
  MatrixAsciiWriteInto(fp, gcl->m_u);
  MatrixAsciiWriteInto(fp, gcl->m_W);
  MatrixAsciiWriteInto(fp, gcl->m_wT);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           read a single class from disk in ascii format
------------------------------------------------------*/
GCLASS *GCasciiReadClassFrom(FILE *fp, GCLASS *gcl)
{
  char class_name[CLASS_NAME_LEN], *cp;

  if (fscanf(fp, "%d %f %d\n", &gcl->classno, &gcl->w0, &gcl->nobs) != 3)
    ErrorReturn(NULL, (ERROR_BADFILE, "GCasciiReadClassFrom: could not scan parms from file"));

  if ((cp = fgetl(class_name, CLASS_NAME_LEN - 1, fp)) == NULL)
    ErrorReturn(NULL, (ERROR_BADFILE, "GCasciiReadClassFrom: could not scan class name from file"));
  strcpy(gcl->class_name, cp);

  /* some of these matrices may already be allocated */
  gcl->m_covariance = MatrixAsciiReadFrom(fp, gcl->m_covariance);
  gcl->m_u = MatrixAsciiReadFrom(fp, gcl->m_u);
  gcl->m_W = MatrixAsciiReadFrom(fp, gcl->m_W);
  gcl->m_wT = MatrixAsciiReadFrom(fp, gcl->m_wT);
  return (gcl);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           training is finished for this class, compute it's
           static components. After this call, the classifier
           should be ready for use.
------------------------------------------------------*/
int GCinit(GCLASSIFY *gc, int classnum)
{
  GCLASS *gcl;
  MATRIX *m_sigma_inverse, *m_uT, *m_tmp, *m_tmp2;
  float det;

  gcl = &gc->classes[classnum];
  if (gcl->nobs <= gc->nvars) /* not enough training data */
    MatrixMakeDiagonal(gcl->m_covariance, gcl->m_covariance);

  det = MatrixDeterminant(gcl->m_covariance);
  if (FZERO(det) || (det < 0.0f)) /* matrix is singular or ill-conditioned */
  {
    gcl->ill_cond = 1;
    return (NO_ERROR);
  }

  m_sigma_inverse = MatrixInverse(gcl->m_covariance, NULL);
  if (!m_sigma_inverse) /* don't really know what to do.... */
  {
    m_sigma_inverse = MatrixIdentity(gc->nvars, NULL);
    gcl->ill_cond = 1;
  }
  m_uT = MatrixTranspose(gcl->m_u, NULL);
  gcl->m_W = MatrixScalarMul(m_sigma_inverse, -0.5f, gcl->m_W);
  m_tmp = MatrixMultiply(m_sigma_inverse, gcl->m_u, NULL);
  gcl->m_wT = MatrixTranspose(m_tmp, gcl->m_wT);
  MatrixFree(&m_tmp);
  m_tmp = MatrixMultiply(m_sigma_inverse, gcl->m_u, NULL);
  m_tmp2 = MatrixMultiply(m_uT, m_tmp, NULL);
  gcl->w0 = -0.5 * (gc->nvars * log(2 * M_PI) + m_tmp2->rptr[1][1] + log(det));

/* log of prior can be added to gcl->w0 */

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

  MatrixFree(&m_sigma_inverse);
  MatrixFree(&m_uT);
  MatrixFree(&m_tmp);
  MatrixFree(&m_tmp2);
  return (NO_ERROR);
}
