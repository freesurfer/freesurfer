/**
 * @brief utilities for Gaussian Classifier Surface Arrays.
 *
 * the heart of the cortical parcellation code - implements an
 * anisotropic nonstationary MRF on the surface.
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
  The GCSA structure allows for computing the posterior probability of
  a given label from potentially multimodal data (ninputs). There might
  be a limit of 3.

  posterior = prior * likelihood, 
  likelihood = exp((x-u)'*inv(S)*(x-u))/(((2pi)^-(k/2))*sqrt(det(S)))
  where x is the intensity vector of the given input, u is the mean vector of
  that label at that location, and S is the covariance matrix of the input
  intensities for that label at that location. In the GCSA class, this
  information is stored in this way:

  Priors: mris_priors is a surface (typically ico7) and *cp_nodes is a
  CP_NODE array with the same number of vertices. Each item in the
  cp_nodes array has an array of CLASS_PRIORS *cp. The CP array is of
  length equal to the number of non-zero labels. For each label, there
  is a prior and an annotation (int label) indicating which label the
  prior refers to.

  Likelihoods: mris_classifiers is a surface (typically
  ico4). gc_nodes is an array of GCSA_NODEs, one for each vertex in
  the classifier surface.  Each node has an array of GCS structures,
  one for each label with non-zero prior; each GCS has a vector of
  means (length = number of inputs) and a covariance matrix
  (ninputs^2). Each node also has an array of annotations (*labels)
  that indicate to what label the corresponding gcs belongs.

  Be prepared: using the annotations is a huge pain.
 */

#ifndef GCSA_H
#define GCSA_H

#include <stdio.h>

#include "mrisurf.h"
#include "mrishash.h"
#include "label.h"
#include "colortab.h"

#define GCSA_MAGIC 0xababcdcd


#define GCSA_INPUT_CURVATURE         0
#define GCSA_INPUT_CURV_FILE         1

/* input flags */
#define GCSA_NORMALIZE               0x00000001
typedef struct
{
  int    type ;
  char   fname[STRLEN] ;
  int    navgs ;
  int    flags ;
}
GCSA_INPUT ;

#define GIBBS_SURFACE_NEIGHBORHOOD   4
#define GIBBS_SURFACE_NEIGHBORS      GIBBS_SURFACE_NEIGHBORHOOD

#define GCSA_MAX_INPUTS    10

typedef struct
{
  VECTOR  *v_means ;
  MATRIX  *m_cov ;
  int     total_training ;
  int     regularized ;
}
GCS, SURFACE_GAUSSIAN_CLASSIFIER ;

typedef struct
{
  float   prior ;
  float   *label_priors[GIBBS_SURFACE_NEIGHBORHOOD] ;
  int    *labels[GIBBS_SURFACE_NEIGHBORHOOD] ;
  short   nlabels[ GIBBS_SURFACE_NEIGHBORHOOD];
  int     total_nbrs[GIBBS_SURFACE_NEIGHBORHOOD] ;
}
CLASS_PRIORS, CP ;

typedef struct
{
  int     nlabels ;
  int     max_labels ;      /* amount allocated */
  int     *labels ;
  CP      *cps ;            /* class priors */
  int     total_training ;  /* total # of times this node was was accessed */
}
CP_NODE ;

typedef struct
{
  int     nlabels ;
  int     max_labels ;   /* amount allocated */
  int     *labels ;
  GCS     *gcs ;       /* classifiers */
  int     total_training ;  /* total # of times this node was was accessed */
}
GCSA_NODE ;

typedef struct
{
  MRI_SURFACE      *mris_priors ;
  MRI_SURFACE      *mris_classifiers ;
  CP_NODE          *cp_nodes;
  GCSA_NODE        *gc_nodes ;
  int              ninputs ;
  int              icno_priors ;
  int              icno_classifiers ;
  MRIS_HASH_TABLE  *mht_priors ;
  MRIS_HASH_TABLE  *mht_classifiers ;
  GCSA_INPUT       inputs[GCSA_MAX_INPUTS] ;
  char             *ptable_fname ;   /* name of color lookup table */
  COLOR_TABLE      *ct ;
}
GAUSSIAN_CLASSIFIER_SURFACE_ARRAY, GCSA ;


GCSA  *GCSAalloc(int ninputs, int icno_priors, int icno_classifiers) ;
int   GCSAfree(GCSA **pgcas) ;
int   GCSAtrainMeans(GCSA *gcsa, MRI_SURFACE *mris) ;
int   GCSAtrainCovariances(GCSA *gcsa, MRI_SURFACE *mris);
int   GCSAnormalizeMeans(GCSA *gcsa) ;
int   GCSAnormalizeCovariances(GCSA *gcsa) ;
int   GCSAwrite(GCSA *gcsa, char *fname) ;
GCSA  *GCSAread(char *fname) ;
int   GCSAlabel(GCSA *gcsa, MRI_SURFACE *mris) ;
int   GCSAdump(GCSA *gcsa, int vno, MRI_SURFACE *mris, FILE *fp) ;
int   GCSAreclassifyUsingGibbsPriors(GCSA *gcsa, MRI_SURFACE *mris) ;
int   GCSAreclassifyLabel(GCSA *gcsa, MRI_SURFACE *mris, LABEL *area) ;
int   GCSAputInputType(GCSA *gcsa, int type, const char *fname, int navgs, int ino,
                       int flags) ;
int   GCSAsetCovariancesToIdentity(GCSA *gcsa) ;
VERTEX *GCSAsourceToPriorVertex     (GCSA *gcsa, VERTEX const *v);
int     GCSAsourceToPriorVertexNo   (GCSA *gcsa, VERTEX const *v);
VERTEX *GCSAsourceToClassifierVertex(GCSA *gcsa, VERTEX const *v);

int dump_gcsan(GCSA_NODE *gcsan, CP_NODE *cpn, FILE *fp, int verbose) ;
int GCSAbuildMostLikelyLabels(GCSA *gcsa, MRI_SURFACE *mris) ;
int GCSArelabelWithAseg(GCSA *gcsa, MRI_SURFACE *mris, MRI *mri_aseg) ;
int GCSAreclassifyMarked(GCSA *gcsa, MRI_SURFACE *mris,int mark, int *exclude_list, int nexcluded) ;
int GCSANclassify(GCSA_NODE *gcsan, CP_NODE *cpn, double *v_inputs, int ninputs, double *pprob, int *exclude_list, int nexcluded);
MRI *GCSApriors2MRI(GCSA *gcsa);
MRI *GCSAlikelihoodMeans2MRI(GCSA *gcsa, int inputno);
int GCSAfill_cpn_holes(GCSA *gcsa, int nitersmax);
int GCSAfill_gcsan_holes(GCSA *gcsa, int nitersmax);

#endif
