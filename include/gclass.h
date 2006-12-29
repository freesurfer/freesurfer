/**
 * @file  gclass.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef GCLASS_H
#define GCLASS_H

#include "matrix.h"

/*
   see Duda and Hart (1st ed.) page 30 for these naming conventions.
*/
#define CLASS_NAME_LEN  30
typedef struct
{
  int     classno ;
  float   w0 ;             /* -.5 uT sigma^-1 u - .5 log |sigma^-1| +log(p) */
  int     nobs ;           /* # of observations in training set */
  char    class_name[CLASS_NAME_LEN] ;
  MATRIX  *m_covariance ;  /* covariance matrix */
  MATRIX  *m_u ;           /* vector of means */
  MATRIX  *m_W ;           /* -0.5 Sigma^-1 */
  MATRIX  *m_wT ;          /* sigma^-1 u */
  int     ill_cond ;       /* is this class ill-conditioned */
}
GAUSSIAN_CLASS, GCLASS ;

typedef struct
{
  int     nclasses ;
  int     nvars ;
  int     type ;
  GCLASS  *classes ;
  float   *log_probabilities ;
}
GAUSSIAN_CLASSIFIER, GCLASSIFY ;

GCLASSIFY *GCalloc(int nclass, int nvars, char *class_names[]) ;
int       GCtrain(GCLASSIFY *gc, int class, MATRIX *m_inputs) ;
int       GCfree(GCLASSIFY **pgc) ;
int       GCclassify(GCLASSIFY *gc, MATRIX *m_input, MATRIX *m_priors,
                     float *prisk) ;
int       GCasciiWriteInto(FILE *fp, GCLASSIFY *gc) ;
GCLASSIFY *GCasciiReadFrom(FILE *fp, GCLASSIFY *gc) ;
int       GCasciiWriteClassInto(FILE *fp, GCLASS *gcl) ;
GCLASS    *GCasciiReadClassFrom(FILE *fp, GCLASS *gcl) ;
int       GCinit(GCLASSIFY *gc, int class) ;

#endif
