/*
 *       FILE NAME:   rbf.h
 *
 *       DESCRIPTION: prototypes and structures for Radial Basis Function
 *                    utilities.
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        5/19/97
 *
*/

#ifndef RBF_H
#define RBF_H

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include "cluster.h"

#define MAX_OUTPUTS  30

typedef struct
{
  int          noutputs ;
  int          ninputs ;
  int          nhidden ;          /* number of RBFs */
  char         *class_names[MAX_OUTPUTS] ;
  CLUSTER_SET  *cs[MAX_OUTPUTS] ; /* one per class */
  float        *min_inputs ;      /* minimum feature values in input space */
  float        *max_inputs ;      /* maximum feature values in input space */
  MATRIX       *m_wij ;           /* weights for linear part of RBF */
  MATRIX       *m_delta_wij ;     /* delta weights for linear part of RBF */
  MATRIX       **m_delta_sigma_inv ;/* for use in momentum */
  VECTOR       **v_delta_means ;   /* for use in momentum */
  VECTOR       **v_z ;            /* zero-mean vectors */
  VECTOR       *v_biases ;        /* biases of linear classifier */
  VECTOR       *v_delta_biases ;  /* delta biases of linear classifier */
  VECTOR       *v_outputs ;       /* outputs of the RBF */
  VECTOR       *v_hidden ;        /* activation of hidden layer */
  void         *parm ;            /* for clustering */
  int          (*get_observation_func)(VECTOR *v_obs, int no, void *parm);
  int          current_class ;    /* for clustering classes separately */
  float        momentum ;         /* current value of momentum */
  float        base_momentum ;
  float        trate ;            /* training rate */
  CLUSTER      **clusters ;       /* pointers to CLUSTER_SET clusters */
  unsigned char *observed ;       /* used for randomizing training order */
  int           nobs ;            /* # of observations in training set */
} RBF ;

RBF   *RBFinit(int ninputs, int noutputs, int max_clusters[], char *names[]) ;
int   RBFtrain(RBF *rbf, int (*get_observation_func)
               (VECTOR *v_obs, int no, void *parm, int same_class,int *class),
               void *parm, float momentum) ;
int   RBFfree(RBF **prbf) ;
int   RBFprint(RBF *rbf, FILE *fp) ;
int   RBFclassify(RBF *rbf, VECTOR *v_obs) ;
int   RBFwrite(RBF *rbf, char *fname) ;
RBF   *RBFread(char *fname) ;
RBF   *RBFcopyWeights(RBF *rbf_src, RBF *rbf_dst) ;
               

#endif
