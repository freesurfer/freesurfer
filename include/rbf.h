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

#define MAX_OUTPUTS  10

typedef struct
{
  int          noutputs ;
  int          ninputs ;
  int          nhidden ;          /* number of RBFs */
  char         *class_names[MAX_OUTPUTS] ;
  CLUSTER_SET  *cs[MAX_OUTPUTS] ; /* one per class */
  MATRIX       *m_wij ;           /* weights for linear part of RBF */
  MATRIX       *m_delta_wij ;     /* delta weights for linear part of RBF */
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
} RBF ;

RBF   *RBFinit(int ninputs, int noutputs, int max_clusters[], char *names[]) ;
int   RBFtrain(RBF *rbf, int (*get_observation_func)
               (VECTOR *v_obs, int no, void *parm, int *class), void *parm) ;
int   RBFfree(RBF **prbf) ;
int   RBFprint(RBF *rbf, FILE *fp) ;
               

#endif
