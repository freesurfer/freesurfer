/**
 * @file  rforest.c
 * @brief types and prototypes for random forest classifier
 *
 * Base on Leo Breiman's random forest classification algorithm
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/03/21 00:56:37 $
 *    $Revision: 1.5 $
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "rforest.h"
#include "error.h"
#include "diag.h"
#include "const.h"
#include "macros.h"
#include "utils.h"


#define MAX_CLASSES 50
static double entropy(int *class_counts, int nclasses, int *Nc) ;

RANDOM_FOREST *
RFalloc(int ntrees, int nfeatures, int nclasses, int max_depth, char **class_names, int nsteps) 
{
  RANDOM_FOREST *rf ;
  TREE          *tree ;
  int           n, c ;

  rf = calloc(1, sizeof(RANDOM_FOREST)) ;
  if (rf == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate RF",
              ntrees,nfeatures,nclasses, max_depth) ;

  rf->ntrees = ntrees ; rf->nfeatures = nfeatures ; rf->nclasses = nclasses ;
  rf->max_depth = max_depth ;
  rf->nsteps = nsteps ;
  rf->class_names = (char **)calloc(rf->nclasses, sizeof(rf->class_names[0])) ;
  if (rf->class_names == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate class names",
              ntrees,nfeatures,nclasses, max_depth) ;
  if (class_names)
  {
    for (c = 0 ; c < rf->nclasses ; c++)
    {
      rf->class_names[c] = (char *)calloc(strlen(class_names[c])+1, sizeof(char)) ;
      if (rf->class_names[c] == NULL)
	ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate class names %d",
		  ntrees,nfeatures,nclasses, max_depth, c) ;
      strcpy(rf->class_names[c], class_names[c]) ;
    }
  }
  

  rf->trees = (TREE *)calloc(rf->ntrees, sizeof(TREE)) ;
  if (rf->trees == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate trees",
              ntrees,nfeatures,nclasses, max_depth) ;
  rf->feature_min = (double *)calloc(rf->nfeatures, sizeof(rf->feature_min[0])) ;
  if (rf->feature_min == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate fmin",
              ntrees,nfeatures,nclasses, max_depth) ;
  rf->feature_max = (double *)calloc(rf->nfeatures, sizeof(rf->feature_max[0])) ;
  if (rf->feature_max == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate fmax",
              ntrees,nfeatures,nclasses, max_depth) ;
  for (n = 0 ; n < rf->ntrees ; n++)
  {
    tree = &rf->trees[n] ;
    tree->root.class_counts = calloc(rf->nclasses, sizeof(tree->root.class_counts[0])) ;
    if (tree->root.class_counts == NULL)
      ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate root %d",
                ntrees,nfeatures,nclasses, max_depth, c) ;
  }
  return(rf) ;
}


static double
compute_info_gain(RF *rf, TREE *tree, NODE *parent, NODE *left, NODE *right, double **training_data, int fno, double thresh) 
{
  double entropy_before, entropy_after, w1, w2 ;
  int    i, tno ;
  NODE   *node ;

  entropy_before = entropy(parent->class_counts, rf->nclasses, tree->root.class_counts) ;
  memset(left->class_counts, 0, rf->nclasses*sizeof(left->class_counts[0])) ;
  memset(right->class_counts, 0, rf->nclasses*sizeof(right->class_counts[0])) ;
  left->total_counts = right->total_counts = 0 ;
  for (tno = 0 ; tno < parent->total_counts ; tno++)
  {
    i = parent->training_set[tno] ;
    if (training_data[i][fno] < thresh)
      node = left ;
    else
      node = right ;
      
    node->class_counts[rf->training_classes[i]]++ ;
    node->training_set[node->total_counts] = i ;
    node->total_counts++ ;
  }
  w1 = (double)left->total_counts / parent->total_counts ;
  w2 = (double)right->total_counts / parent->total_counts ;
  entropy_after = 
    w1 * entropy(left->class_counts, rf->nclasses, tree->root.class_counts) +
    w2 * entropy(right->class_counts, rf->nclasses, tree->root.class_counts) ;
  return(entropy_before - entropy_after) ;
}
static int
adjust_optimal_threshold(RF *rf, TREE *tree, NODE *parent, NODE *left, NODE *right, double **training_data, int fno, double *pbest_thresh) 
{
  double previous_thresh, next_thresh, info_gain, best_info_gain, step, thresh, best_thresh ;

  best_thresh = *pbest_thresh ;
  best_info_gain = compute_info_gain(rf, tree, parent, left, right, training_data, fno,best_thresh) ;
  step = (rf->feature_max[fno] - rf->feature_min[fno]) / (10*rf->nsteps-1);
  for (thresh = best_thresh ; thresh <= rf->feature_max[fno] ; thresh += step)
  {
    info_gain = compute_info_gain(rf, tree, parent, left, right, training_data, fno,thresh) ;
    if (info_gain > best_info_gain)
    {
      best_thresh = thresh ;
      best_info_gain = info_gain ;
    }

    if (info_gain < best_info_gain)
      break ;
  }
  next_thresh = thresh - step ;
  for (thresh = best_thresh ; thresh >= rf->feature_min[fno] ; thresh -= step)
  {
    info_gain = compute_info_gain(rf, tree, parent, left, right, training_data, fno,thresh) ;
    if (info_gain > best_info_gain)
    {
      best_thresh = thresh ;
      best_info_gain = info_gain ;
    }

    if (info_gain < best_info_gain)
      break ;
  }
  previous_thresh = thresh+step ;

  thresh = (next_thresh + previous_thresh) / 2 ;  // maximize margin
  info_gain = compute_info_gain(rf, tree, parent, left, right, training_data, fno,thresh) ;
  if (info_gain < best_info_gain) // don't use it
    compute_info_gain(rf, tree, parent, left, right, training_data, fno,best_thresh) ;
  else
    best_thresh = thresh ;
  *pbest_thresh = best_thresh ;
  
  return(NO_ERROR) ;
}

static double
find_optimal_threshold(RF *rf, TREE *tree, NODE *parent, NODE *left, NODE *right, double **training_data, int ntraining, int fno, double *pinfo_gain, int nsteps)
{
  double step, thresh, best_thresh, info_gain, best_info_gain ;

  step = (rf->feature_max[fno] - rf->feature_min[fno]) / (nsteps-1) ;
  if (FZERO(step))
    return(0.0) ;

  best_info_gain = -1e10; best_thresh = 0 ;
  for (thresh = rf->feature_min[fno] ; thresh < rf->feature_max[fno] ; thresh += step)
  {
    info_gain = compute_info_gain(rf, tree, parent, left, right, training_data, fno,thresh) ;
    if (info_gain > best_info_gain && left->total_counts>0 && right->total_counts>0)
    {
      best_info_gain = info_gain ;
      best_thresh = thresh ;
    }
  }
  *pinfo_gain = best_info_gain ;
  return(best_thresh) ;
}

static int
find_optimal_feature_and_threshold(RANDOM_FOREST *rf, TREE *tree, NODE *parent, NODE *left, NODE *right, int *training_classes, double **training_data, int ntraining)
{
  double  info_gain, best_info_gain, thresh, best_thresh ;
  int     f, best_f, fno, nsteps;

  best_info_gain = -1 ; best_f = -1 ; best_thresh = 0 ;
  for (f = 0 ; f < tree->nfeatures ; f++)
  {
    fno = tree->feature_list[f] ;
    nsteps = rf->nsteps ;
    do
    {
      thresh = find_optimal_threshold(rf, tree, parent, left, right, rf->training_data, rf->ntraining,
				      fno, &info_gain, nsteps);
      if (info_gain < 0)
	DiagBreak() ;
      nsteps *= 5 ;
      if (nsteps > 1000*rf->nsteps)
	break ;  // don't keep trying forever
    } while (info_gain <= 0) ;

    if (info_gain > best_info_gain)
    {
      best_info_gain = info_gain ;
      best_f = fno ;
      best_thresh = thresh ;
    }
  }
  if (best_f < 0)
    return(0) ;
  info_gain = compute_info_gain(rf, tree, parent, left, right, training_data, best_f,best_thresh) ;
  adjust_optimal_threshold(rf, tree, parent, left, right, training_data, best_f,&best_thresh);
  parent->thresh = best_thresh ;
  parent->feature = best_f ;
  return(1) ;
}


static NODE *
rfAllocateNode(int ntraining, int depth, int nclasses) 
{
  NODE *node ;
  node = (NODE *)calloc(1, sizeof(NODE)) ;
  if (node == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFtrain: could not allocate node") ;
  node->class_counts = (int *)calloc(nclasses, sizeof(node->class_counts[0])) ;
  if (node->class_counts == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFtrain: could not allocate node class counts") ;
  node->training_set = (int *)calloc(ntraining, sizeof(node->training_set[0])) ;
  if (node->training_set == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFtrain: could not allocate node training_set") ;
  node->depth = depth ;
  return(node) ;
}

static int
rfFreeNode(NODE **pnode)
{
  NODE *node ;

  node = *pnode ; *pnode = NULL ;
  free(node->class_counts) ;
  free(node->training_set) ;
  free(node) ;
  return(0) ;
}

static int
rfTrainNode(RANDOM_FOREST *rf, TREE *tree, NODE *node, int *training_classes, 
            double **training_data, int ntraining)
{
  if (node->left == NULL)   // not trained yet
  {
    node->left = rfAllocateNode(node->total_counts, node->depth+1, rf->nclasses) ;
    node->right = rfAllocateNode(node->total_counts, node->depth+1, rf->nclasses) ;
    if (find_optimal_feature_and_threshold(rf, tree, node, node->left, node->right, 
					   training_classes, 
					   training_data, rf->ntraining) == 0)
    {
      // couldn't find a threshold to improve separation
      rfFreeNode(&node->left) ; rfFreeNode(&node->right) ;
      return(0) ;
    }
  }
  if (node->depth > tree->depth)
    tree->depth = node->depth ;
  if (node->depth < rf->max_depth)
  {
    if (entropy(node->left->class_counts, rf->nclasses, tree->root.class_counts) > 0)
      rfTrainNode(rf, tree, node->left, training_classes, training_data, ntraining) ;
    if (entropy(node->right->class_counts, rf->nclasses, tree->root.class_counts) > 0)
      rfTrainNode(rf, tree, node->right, training_classes, training_data, ntraining) ;
  }
  return(1) ;
}


static int
rfFindNodeLeaves(TREE *tree, NODE *node) 
{
  if (node->left == NULL)
  {
    if (tree->leaves)
    {
      if (tree->nleaves >= tree->max_leaves &&
	  tree->max_leaves > 0)
	DiagBreak() ;
      tree->leaves[tree->nleaves] = node ;
    }
    tree->nleaves++ ;
  }
  else
  {
    rfFindNodeLeaves(tree, node->left) ;
    rfFindNodeLeaves(tree, node->right) ;
  }
  return(NO_ERROR) ;
}

static int
rfFindLeaves(TREE *tree)
{
  if (tree->leaves)
  {
    free(tree->leaves) ;
    tree->leaves = NULL ;
  }
  tree->nleaves = 0 ;
  rfFindNodeLeaves(tree, &tree->root) ;  // when tree->leaves==NULL just counts #
  tree->max_leaves = tree->nleaves ;
  tree->leaves = (NODE **)calloc(tree->nleaves, sizeof(NODE *)) ;
  if (tree->leaves == NULL)
    ErrorExit(ERROR_NOMEMORY, "couldn't allocate %d leaves\n", tree->nleaves) ;
  tree->nleaves = 0 ;
  rfFindNodeLeaves(tree, &tree->root) ;  // now actually fill them in

  return(NO_ERROR) ;
}
static int
rfTrainTree(RANDOM_FOREST *rf, TREE *tree, int *training_classes, double **training_data, int ntraining)
{
  int    done = 0, n, iter, tno, fno, f, i ;
  double total_entropy, last_f ;

  total_entropy = entropy(tree->root.class_counts, rf->nclasses, tree->root.class_counts) ;
  
  // make sure there is at least one feature with nonzero range
  for (f = 0 ; f < tree->nfeatures ; f++)
  {
    fno = tree->feature_list[f] ;
    last_f = training_data[tree->root.training_set[0]][fno] ;
    for (i = 1 ; i < tree->root.total_counts ; i++)
    {
      tno = tree->root.training_set[i] ;
      if (training_data[tno][fno] != last_f)
	break ;
    }
    if (i < tree->root.total_counts) 
      break ;  // found one feature with some spread
  }
  if (f >= tree->nfeatures)  // all features are identical - can't train
  {
    rfFindLeaves(tree) ;
    return(ERROR_BADPARM) ;
  }
  iter = 0 ;
  while (!done && !FZERO(total_entropy))
  {
    done = rfTrainNode(rf, tree, &tree->root, training_classes, training_data, rf->ntraining);
    rfFindLeaves(tree) ;
    for (total_entropy = 0.0, n = 0 ; n < tree->nleaves ; n++)
      total_entropy += entropy(tree->leaves[n]->class_counts, rf->nclasses, tree->root.class_counts) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      printf("\taverage leaf entropy = %2.4f, nleaves = %d, max depth %d\n",
	     total_entropy/tree->nleaves, tree->nleaves, tree->depth) ;
    if (iter++ > 10)
      break ;
  }
  if (tree->nleaves == 0)  // only if loop above never executed
    rfFindLeaves(tree) ;
  return(NO_ERROR) ;
}
int
RFtrain(RANDOM_FOREST *rf, double feature_fraction, double training_fraction, int *training_classes, double **training_data, int ntraining)
{
  int    n, i, nfeatures_per_tree, *feature_permutation, *training_permutation, f, index, start_no, end_no,
         ntraining_per_tree ;
  TREE   *tree ;

  rf->ntraining = ntraining ;
  rf->training_data = training_data ;
  rf->training_classes = training_classes ;
  rf->feature_fraction = feature_fraction ;
  rf->training_fraction = training_fraction ;

  for (f = 0 ; f < rf->nfeatures ; f++)
  {
    rf->feature_min[f] = 1e20 ;
    rf->feature_max[f] = -1e20 ;
    for (i = 0 ; i < ntraining ; i++)
    {
      if (training_data[i][f] < rf->feature_min[f])
        rf->feature_min[f] = training_data[i][f] ;
      if (training_data[i][f] > rf->feature_max[f])
        rf->feature_max[f] = training_data[i][f] ;
    }
  }

  feature_permutation = calloc(rf->nfeatures, sizeof(int)) ;
  if (feature_permutation == NULL)
    ErrorExit(ERROR_NOMEMORY, "could not allocate feature permutation vector");
  training_permutation = calloc(ntraining, sizeof(int)) ;
  if (training_permutation == NULL)
    ErrorExit(ERROR_NOMEMORY, "could not allocate feature permutation vector");

  nfeatures_per_tree = nint((double)rf->nfeatures * feature_fraction) ;
  ntraining_per_tree = nint((double)rf->ntraining * training_fraction) ;

  compute_permutation(rf->nfeatures, feature_permutation) ;
  compute_permutation(ntraining, training_permutation) ;
  for (n = 0 ; n < rf->ntrees ; n++)  // train each tree
  {
    tree = &rf->trees[n] ;

    // randomize what features this tree will use
    tree->feature_list = (int *)calloc(nfeatures_per_tree, sizeof(tree->feature_list[0]));
    if (tree->feature_list == NULL)
      ErrorExit(ERROR_NOMEMORY, "RFtrain: could not allocate feature list %d (%d)",
                n,nfeatures_per_tree) ;
    tree->nfeatures = nfeatures_per_tree ;
    if (rf->ntrees > 1)
      start_no = nint(n*((double)(rf->nfeatures-nfeatures_per_tree))/(rf->ntrees-1.0)) ;
    else
      start_no = 0 ; // only 1 tree
    end_no = MIN(rf->nfeatures-1, start_no+nfeatures_per_tree-1) ;
    for (i = start_no  ; i <= end_no  ; i++)
      tree->feature_list[i-start_no] = feature_permutation[i] ;
    

    // randomize what training data this tree will use
    tree->root.training_set = (int *)calloc(ntraining, sizeof(tree->root.training_set[0])) ;
    if (tree->root.training_set == NULL)
      ErrorExit(ERROR_NOMEMORY, "RFtrain: could not allocate root training set") ;
    tree->root.total_counts = 0 ;
    if (rf->ntrees > 1)
      start_no = nint(n*((double)(rf->ntraining-ntraining_per_tree))/(rf->ntrees-1.0)) ;
    else
      start_no = 0 ; // only 1 tree
    end_no = MIN(rf->ntraining-1, start_no+ntraining_per_tree-1) ;
    for (i = start_no  ; i <= end_no  ; i++)
    {
      index = training_permutation[i] ;
      tree->root.class_counts[training_classes[index]]++ ;
      tree->root.training_set[tree->root.total_counts] = index ;
      tree->root.total_counts++ ;
    }

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      printf("tree %d: initial entropy = %f\n", n, 
	     entropy(tree->root.class_counts, rf->nclasses, tree->root.class_counts)) ;
    rfTrainTree(rf, tree, training_classes, training_data, rf->ntraining) ;
  }

  free(feature_permutation) ;
  free(training_permutation) ;
  return(NO_ERROR) ;
}
int
RFtrainTree(RANDOM_FOREST *rf, int tno, int *training_classes, double **training_data, int ntraining)
{
  int    i, f ;
  TREE   *tree ;


  rf->training_data = training_data ;
  rf->training_classes = training_classes ;

  for (f = 0 ; f < rf->nfeatures ; f++)
  {
    rf->feature_min[f] = 1e20 ;
    rf->feature_max[f] = -1e20 ;
    for (i = 0 ; i < ntraining ; i++)
    {
      if (training_data[i][f] < rf->feature_min[f])
        rf->feature_min[f] = training_data[i][f] ;
      if (training_data[i][f] > rf->feature_max[f])
        rf->feature_max[f] = training_data[i][f] ;
    }
  }

  tree = &rf->trees[tno] ;

  tree->feature_list = (int *)calloc(rf->nfeatures, sizeof(tree->feature_list[0]));
  if (tree->feature_list == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFtrain: could not allocate feature list %d (%d)",
	      tno,rf->nfeatures) ;
  tree->nfeatures = rf->nfeatures ;
  
  for (i = 0  ; i < rf->nfeatures  ; i++)
    tree->feature_list[i] = i ;
  

  tree->root.training_set = (int *)calloc(ntraining, sizeof(tree->root.training_set[0])) ;
  if (tree->root.training_set == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFtrainTree: could not allocate root training set") ;

  for (i = 0  ; i < ntraining  ; i++)
  {
    tree->root.class_counts[training_classes[i]]++ ;
    tree->root.training_set[tree->root.total_counts] = i ; // should be +ntraining
    tree->root.total_counts++ ;
  }
  rf->ntraining += ntraining ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("tree %d: initial entropy = %f\n", tno, 
	   entropy(tree->root.class_counts, rf->nclasses, tree->root.class_counts)) ;
  rfTrainTree(rf, tree, training_classes, training_data, rf->ntraining) ;

  return(NO_ERROR) ;
}

static double
entropy(int *class_counts, int nclasses, int *Nc)
{
  int    c ;
  double ent, p, total ;

  for (total = 0.0, c = 0 ; c < nclasses ; c++)
    if (Nc[c] > 0)
      total += (double)class_counts[c]/(double)Nc[c] ;
  if (FZERO(total))
    return(0.0) ;
  for (ent = 0.0, c = 0 ; c < nclasses ; c++)
  {
    if (Nc[c] > 0)
      p = ((double)class_counts[c]/(double)Nc[c]) / total ;
    else
      p = 0 ;
    if (p > 0)
      ent -= p * log(p) ;
  }
  return(ent) ;
}

static int
rfWriteNode(RANDOM_FOREST *rf, NODE *node, FILE *fp)
{
  int c ;

  fprintf(fp, "NODE %d %d %d %lf\n", node->left == NULL, node->depth, node->feature, node->thresh) ;
  for (c = 0 ; c < rf->nclasses ; c++)
    fprintf(fp, "%d ", node->class_counts[c]) ;
  fprintf(fp, "\n") ;
#if 0
  for (n = 0 ; n < node->total_counts ; n++)
    fprintf(fp, "%d ", node->training_set[n]) ;
#endif
  fprintf(fp, "\nNODE: END\n") ;
  if (node->left)
  {
    rfWriteNode(rf, node->left, fp) ;
    rfWriteNode(rf, node->right, fp) ;
  }
  return(NO_ERROR) ;
}
static int 
rfWriteTree(RANDOM_FOREST *rf, TREE *tree, FILE *fp)
{
  int f ;

  fprintf(fp, "TREE %d %d %d\n", tree->depth, tree->nleaves, tree->nfeatures) ;
  for (f = 0 ; f < tree->nfeatures ; f++)
    fprintf(fp, "%d %d\n", f, tree->feature_list[f]) ;
  rfWriteNode(rf, &tree->root, fp) ;
  fprintf(fp, "TREE: END\n") ;
  return(NO_ERROR) ;
}
static int
rfReadNode(RANDOM_FOREST *rf, NODE *node, FILE *fp)
{
  int   c, leaf ;
  char  line[MAX_LINE_LEN], *cp ;

  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  sscanf(cp, "NODE %d %d %d %lf\n", &leaf, &node->depth, &node->feature, &node->thresh) ;
  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  cp = strtok(line, " ") ;
  for (node->total_counts = 0, c = 0 ; c < rf->nclasses ; c++)
  {
    sscanf(cp, "%d ", &node->class_counts[c]) ;
    cp = strtok(NULL, " ") ;
    node->total_counts += node->class_counts[c] ;
  }

  node->training_set = (int *)calloc(node->total_counts, sizeof(node->training_set[0])) ;
#if 0
  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  cp = strtok(line, " ") ;
  for (n = 0 ; n < node->total_counts ; n++)
  {
    sscanf(cp, "%d ", &node->training_set[n]) ;
    cp = strtok(NULL, " ") ;
  }
#endif
  cp = fgetl(line, MAX_LINE_LEN, fp) ;  // NODE: END line
  if (leaf == 0)  // read the next pair of nodes
  {
    node->left = rfAllocateNode(rf->ntraining, node->depth+1, rf->nclasses) ;
    node->right = rfAllocateNode(rf->ntraining, node->depth+1, rf->nclasses) ;
    rfReadNode(rf, node->left, fp) ;
    rfReadNode(rf, node->right, fp) ;
  }
  return(NO_ERROR) ;
}
static int 
rfReadTree(RANDOM_FOREST *rf, TREE *tree, FILE *fp)
{
  int   f ;
  char  line[MAX_LINE_LEN], *cp ;

  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  sscanf(cp, "TREE %d %d %d\n", &tree->depth, &tree->nleaves, &tree->nfeatures) ;
  tree->feature_list = (int *)calloc(tree->nfeatures, sizeof(tree->feature_list[0]));
  if (tree->feature_list == NULL)
    ErrorExit(ERROR_NOMEMORY, "rfReadTree: could not allocate feature list  (%d)",
	      tree->nfeatures) ;
  for (f = 0 ; f < tree->nfeatures ; f++)
  {
    cp = fgetl(line, MAX_LINE_LEN, fp) ;
    sscanf(cp, "%*d %d", &tree->feature_list[f]) ;
  }

  rfReadNode(rf, &tree->root, fp) ;
  cp = fgetl(line, MAX_LINE_LEN, fp) ;   // TREE: END line
  rfFindLeaves(tree) ;
  return(NO_ERROR) ;
}
int
RFwrite(RANDOM_FOREST *rf, char *fname) 
{
  FILE   *fp ;
  int    err ;

  fp = fopen(fname, "w") ;
  if (fp == NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "RFwrite(%s): could not open file",fname));
  err = RFwriteInto(rf, fp);
  fclose(fp) ;
  return(err) ;
}
int
RFwriteInto(RANDOM_FOREST *rf, FILE *fp)
{
  int    c, f, n ;

  fprintf(fp, "%d %d %d %d %d %d %lf\n", 
          rf->nfeatures, rf->nclasses, rf->ntrees, rf->max_depth, rf->ntraining,rf->nsteps,rf->feature_fraction) ;
  for (c = 0 ; c < rf->nclasses ; c++)
    fprintf(fp, "%s\n", rf->class_names[c]) ;

  for (f = 0 ; f < rf->nfeatures ; f++)
    fprintf(fp, "%lf %lf\n", rf->feature_min[f], rf->feature_max[f]) ;

  for (n = 0 ; n < rf->ntrees ; n++)
    rfWriteTree(rf, &rf->trees[n], fp) ;

  return(NO_ERROR) ;
}


static NODE *
rfFindLeaf(RF *rf, NODE *node, double *feature)
{
  if (DIAG_VERBOSE_ON)
  {
    if (node->left)  // not a leaf
      printf("\tnode: feature #%d: thresh = %2.1f, input val = %2.1f\n",
	      node->feature, node->thresh, feature[node->feature]) ;
  }
  if (node->left == NULL)
    return(node) ;
  if (feature[node->feature] < node->thresh)
    return(rfFindLeaf(rf, node->left, feature)) ;
  else
    return(rfFindLeaf(rf, node->right, feature)) ;
}

int
RFcomputeOutOfBagCorrect(RANDOM_FOREST *rf, int *training_classes, double **training_data,int ntraining)
{
  int    max_class = -1, n, c, max_count, i, correct, i2 ;
  NODE   *node ;
  TREE   *tree ;
  double class_counts[MAX_CLASSES], total_count;
  

  for (correct = i = 0 ; i < ntraining ; i++)
  {
    memset(class_counts, 0, sizeof(class_counts)) ;

    if (i == Gdiag_no)
      Gdiag |= DIAG_VERBOSE ;
    for (n = 0 ; n < rf->ntrees ; n++)
    {
      tree = &rf->trees[n] ;
      for (i2 = 0 ; i2 < tree->root.total_counts ; i2++)
	if (tree->root.training_set[i2] == i)  
	  break;  // this tree was trained on this data - don't use it for accuracy estimation
      if (i2 < tree->root.total_counts)
	continue ;// this tree was trained on this data - don't use it for accuracy estimation

      node = rfFindLeaf(rf, &rf->trees[n].root, training_data[i]) ;
      for (total_count = c = 0 ; c < rf->nclasses ; c++)
      {
	if (tree->root.class_counts[c] > 0)
	  total_count += (double)node->class_counts[c]/(double)tree->root.class_counts[c] ;
      }
      for (c = 0 ; c < rf->nclasses ; c++)
      {
	if (tree->root.class_counts[c] > 0)
	  class_counts[c] += ((double)node->class_counts[c]/(double)tree->root.class_counts[c])/total_count;

      }
      
      if (DIAG_VERBOSE_ON)
      {
	printf("tree %d: ", n) ;
	for (c = 0 ; c < rf->nclasses ; c++)
	  printf("C%d=%d ", c, tree->root.class_counts[c]);

	for (c = 0 ; c < tree->nleaves ; c++)  // diagnostic
	{
	  int c2 ;
	  
	  if (node == tree->leaves[c])
	  {
	    printf("leaf %d: ", c) ;
	    for (c2 = 0 ; c2 < rf->nclasses ; c2++)
	      printf("%2.2f ", ((double)node->class_counts[c2]/(double)tree->root.class_counts[c2])/total_count) ;
	    printf("\n") ;
	    break ;
	  }
	}
      }
    }

    
    for (total_count = max_count = c = 0 ; c < rf->nclasses ; c++)
    {
      total_count += class_counts[c] ;
      if (class_counts[c] > max_count)
      {
	max_count = class_counts[c] ;
	max_class = c ;
      }
    }

    if (training_classes[i] == max_class)
      correct++ ;

#if 0
    if (DIAG_VERBOSE_ON)
    {
      for (c = 0 ; c < rf->nclasses ; c++)
      {
	if (true_class >= 0)
	{
	  if (c == max_class)
	  {
	    printf("%s: p = %2.2f ", 
		   rf->class_names[c], 100*(double)class_counts[c]/total_count) ;
	    if (c==true_class)
	      printf("CORRECT\n") ;
	    else
	      printf("(true = %s, p = %2.1f)\n",
		     rf->class_names[true_class], 
		     100*(double)class_counts[true_class]/total_count) ;
	  }
	}
	else
	  printf("%s: p = %2.2f %s\n", 
		 rf->class_names[c], 100*(double)class_counts[c]/total_count, 
		 c==max_class?"MAX":"") ;
      }
    }
#endif
    if (i == Gdiag_no)
      Gdiag &= ~DIAG_VERBOSE ;
  }

  return(correct) ;
}

int
RFclassify(RANDOM_FOREST *rf, double *feature, double *p_pval, int true_class)
{
  int    max_class = -1, n, c ;
  NODE   *node ;
  double max_count, class_counts[MAX_CLASSES] ;
  TREE   *tree ;

  memset(class_counts, 0, sizeof(class_counts)) ;
  for (n = 0 ; n < rf->ntrees ; n++)
  {
    tree = &rf->trees[n] ;

    // now prune trees that only had one training class as they just bias labeling
    for (c = 0 ; c < rf->nclasses ; c++)
      if (tree->root.class_counts[c] == tree->root.total_counts)
	break ;
    if (c < rf->nclasses)  // one class has all counts - prune this tree
      continue ;  // don't use this tree in classification
    node = rfFindLeaf(rf, &tree->root, feature) ;
    for (c = 0 ; c < rf->nclasses ; c++)
      class_counts[c] += (double)node->class_counts[c]/tree->root.total_counts ;

    if (DIAG_VERBOSE_ON)
    {
      int t ;

      printf("tree %d: ", n) ;
      for  (c = 0 ; c < rf->nclasses ; c++)
	printf(" C%d-%d ", c, tree->root.class_counts[c]) ;
      for (t = 0 ; t < tree->nleaves ; t++)  // diagnostic
	if (node == tree->leaves[t])
	{
	  printf("leaf %d: ", t) ;
	  for  (c = 0 ; c < rf->nclasses ; c++)
	    printf(" C%d=%.0f ", c, 
		   100.0*(double)node->class_counts[c]/tree->root.total_counts) ;
	  printf("\n") ;
	}
    }
  }

  for (max_count = c = 0 ; c < rf->nclasses ; c++)
  {
    class_counts[c] /= rf->ntrees ;
    if (class_counts[c] > max_count)
    {
      max_count = class_counts[c] ;
      max_class = c ;
      if (p_pval)
	*p_pval = (double)class_counts[c] ;
    }
  }

  if (DIAG_VERBOSE_ON) for (c = 0 ; c < rf->nclasses ; c++)
  {
    if (true_class >= 0)
    {
      if (c == max_class)
      {
        printf("%s: p = %2.2f ", rf->class_names[c], 100.0*class_counts[c]) ;
        if (c==true_class)
          printf("CORRECT\n") ;
        else
          printf("(true = %s, p = %2.1f)\n", rf->class_names[true_class], 100.0*class_counts[true_class]) ;
      }
    }
    else
      printf("%s: p = %2.2f %s\n", 
             rf->class_names[c], 100*(double)class_counts[c], c==max_class?"MAX":"") ;

  }
  return(max_class) ;
}

RANDOM_FOREST *
RFread(char *fname)
{
  FILE   *fp ;
  RF     *rf ;

  fp = fopen(fname, "r") ;
  if (fp == NULL)
    ErrorExit(ERROR_NOFILE, "RFread(%s): could not open file", fname) ;

  rf = RFreadFrom(fp) ;
  fclose(fp) ;
  return(rf) ;
}
RANDOM_FOREST *
RFreadFrom(FILE *fp)
{
  RF     *rf ;
  int    nfeatures, nclasses, ntrees, max_depth, ntraining, nsteps, c, n, f ;
  char   line[MAX_LINE_LEN], *cp, *class_names[MAX_CLASSES] ;
  double feature_fraction ;

  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  sscanf(cp, "%d %d %d %d %d %d %lf", 
         &nfeatures, &nclasses, &ntrees, &max_depth, &ntraining,&nsteps, &feature_fraction) ;
  for (c = 0 ; c < nclasses ; c++)
  {
    cp = fgetl(line, MAX_LINE_LEN, fp) ;
    class_names[c] = (char *)calloc(strlen(cp)+1, sizeof(char)) ;
    strcpy(class_names[c], cp) ;
  }
  rf = RFalloc(ntrees, nfeatures, nclasses, max_depth, class_names, nsteps) ;
  rf->feature_fraction = feature_fraction ;
  rf->ntraining = ntraining ;

  rf->feature_min = (double *)calloc(rf->nfeatures, sizeof(rf->feature_min[0])) ;
  if (rf->feature_min == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate fmin",
              ntrees,nfeatures,nclasses, max_depth) ;
  rf->feature_max = (double *)calloc(rf->nfeatures, sizeof(rf->feature_max[0])) ;
  if (rf->feature_max == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate fmax",
              ntrees,nfeatures,nclasses, max_depth) ;
  for (f = 0 ; f < rf->nfeatures ; f++)
  {
    cp = fgetl(line, MAX_LINE_LEN, fp) ;
    sscanf(cp, "%lf %lf\n", &rf->feature_min[f], &rf->feature_max[f]) ;
  }

  for (n = 0 ; n < rf->ntrees ; n++)
    rfReadTree(rf, &rf->trees[n], fp) ;
  for (c = 0 ; c < nclasses ; c++)
    free(class_names[c]) ;
  return(rf) ;
}

int
RFsetNumberOfClasses(RANDOM_FOREST *rf, int nclasses)
{
  int   n, c ;
  char  **old_class_names ;
  int   *old_class_counts ;
  TREE  *tree ;

  if (nclasses > rf->nclasses)
  {
    old_class_names = rf->class_names ;
    rf->class_names = (char **)calloc(nclasses, sizeof(rf->class_names[0])) ;
    if (rf->class_names == NULL)
      ErrorExit(ERROR_NOMEMORY, "RFsetNumberOfClasses(%d): could not allocate class names",
		nclasses) ;
    
    for (c = 0 ; c < rf->nclasses ; c++)
    {
      if (old_class_names[c] == NULL)
	continue ;
      rf->class_names[c] = (char *)calloc(strlen(old_class_names[c])+1, sizeof(char)) ;
      strcpy(rf->class_names[c], old_class_names[c]) ;
      free(old_class_names[c]) ;
    }
    free(old_class_names) ;

    for (n = 0 ; n < rf->ntrees ; n++)
    {
      tree = &rf->trees[n] ;
      old_class_counts = tree->root.class_counts ;
      tree->root.class_counts = calloc(nclasses, sizeof(tree->root.class_counts[0])) ;
      if (tree->root.class_counts == NULL)
      ErrorExit(ERROR_NOMEMORY, "RFsetNumberOfClasses(%d): could not allocate class names",
		nclasses) ;
      memmove(tree->root.class_counts, old_class_counts, 
	      rf->nclasses*sizeof(tree->root.class_counts[0])) ;
    }
  }
  rf->nclasses = nclasses ;
  return(NO_ERROR) ;
}
