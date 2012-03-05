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
 *    $Date: 2012/03/05 16:03:11 $
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
static double entropy(int *class_counts, int nclasses) ;

RANDOM_FOREST *
RFalloc(int ntrees, int nfeatures, int nclasses, int max_depth, char **class_names) 
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
  rf->nsteps = 100 ;
  rf->class_names = (char **)calloc(rf->nclasses, sizeof(rf->class_names[0])) ;
  if (rf->class_names == NULL)
    ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate class names",
              ntrees,nfeatures,nclasses, max_depth) ;
  for (c = 0 ; c < rf->nclasses ; c++)
  {
    rf->class_names[c] = (char *)calloc(strlen(class_names[c])+1, sizeof(char)) ;
    if (rf->class_names[c] == NULL)
      ErrorExit(ERROR_NOMEMORY, "RFalloc(%d, %d, %d, %d): could not allocate class names %d",
                ntrees,nfeatures,nclasses, max_depth, c) ;
    strcpy(rf->class_names[c], class_names[c]) ;
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

static int
compute_permutation(int num, int *vec) 
{
  int n, index, tmp ;

  for (n = 0 ; n < num ; n++)
    vec[n] = n ;

  for (n = 0 ; n < num ; n++)
  {  
    index = (int)randomNumber(0.0, (double)(num-0.0001)) ;
    tmp = vec[index] ;
    vec[index] = vec[n] ;
    vec[n] = tmp ;
  }
  return(NO_ERROR) ;
}

static double
compute_info_gain(RF *rf, NODE *parent, NODE *left, NODE *right, double **training_data, int ntraining, int fno, double thresh) 
{
  double entropy_before, entropy_after, w1, w2 ;
  int    i ;
  NODE   *node ;

  entropy_before = entropy(parent->class_counts, rf->nclasses) ;
  memset(left->training_set, 0, rf->ntraining*sizeof(left->training_set[0])) ;
  memset(left->class_counts, 0, rf->nclasses*sizeof(left->class_counts[0])) ;
  memset(right->class_counts, 0, rf->nclasses*sizeof(right->class_counts[0])) ;
  memset(right->training_set, 0, rf->ntraining*sizeof(right->training_set[0])) ;
  left->total_counts = right->total_counts = 0 ;
  for (i = 0 ; i < ntraining ; i++)
  {
    if (parent->training_set[i])  
    {
      if (training_data[i][fno] < thresh)
        node = left ;
      else
        node = right ;
      
      node->training_set[i] = 1 ;
      node->class_counts[rf->training_classes[i]]++ ;
      node->total_counts++ ;
    }
  }
  w1 = (double)left->total_counts / parent->total_counts ;
  w2 = (double)right->total_counts / parent->total_counts ;
  entropy_after = 
    w1 * entropy(left->class_counts, rf->nclasses) +
    w2 * entropy(right->class_counts, rf->nclasses) ;
  return(entropy_before - entropy_after) ;
}

static double
find_optimal_threshold(RF *rf, NODE *parent, NODE *left, NODE *right, double **training_data, int ntraining, int fno, int *training_set, double *pinfo_gain)
{
  double step, thresh, best_thresh, info_gain, best_info_gain ;

  step = (rf->feature_max[fno] - rf->feature_min[fno]) / (rf->nsteps-1) ;

  best_info_gain = -1e10; best_thresh = 0 ;
  for (thresh = rf->feature_min[fno] ; thresh < rf->feature_max[fno] ; thresh += step)
  {
    info_gain = compute_info_gain(rf, parent, left, right, training_data, ntraining, fno,thresh) ;
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
  int     f, best_f, fno;

  best_info_gain = -1 ; best_f = -1 ; best_thresh = 0 ;
  for (f = 0 ; f < tree->nfeatures ; f++)
  {
    fno = tree->feature_list[f] ;
    thresh = find_optimal_threshold(rf, parent, left, right, rf->training_data, rf->ntraining,
                                    fno, parent->training_set, &info_gain);
    if (info_gain > best_info_gain)
    {
      best_info_gain = info_gain ;
      best_f = fno ;
      best_thresh = thresh ;
    }
  }
  info_gain = compute_info_gain(rf, parent, left, right, training_data, rf->ntraining, 
                                best_f,best_thresh) ;
  parent->thresh = best_thresh ;
  parent->feature = best_f ;
  return(NO_ERROR) ;
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
    ErrorExit(ERROR_NOMEMORY, "RFtrain: could not allocate node training set") ;
  node->depth = depth ;
  return(node) ;
}

static int
rfTrainNode(RANDOM_FOREST *rf, TREE *tree, NODE *node, int *training_classes, 
            double **training_data, int ntraining)
{
  if (node->left == NULL)   // not trained yet
  {
    node->left = rfAllocateNode(rf->ntraining, node->depth+1, rf->nclasses) ;
    node->right = rfAllocateNode(rf->ntraining, node->depth+1, rf->nclasses) ;
    find_optimal_feature_and_threshold(rf, tree, node, node->left, node->right, 
                                       training_classes, 
                                       training_data, rf->ntraining) ;
  }
  if (node->depth > tree->depth)
    tree->depth = node->depth ;
  if (node->depth < rf->max_depth-1)
  {
    if (entropy(node->left->class_counts, rf->nclasses) > 0)
      rfTrainNode(rf, tree, node->left, training_classes, training_data, ntraining) ;
    if (entropy(node->right->class_counts, rf->nclasses) > 0)
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
      tree->leaves[tree->nleaves] = node ;
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
  tree->nleaves = 0 ;
  rfFindNodeLeaves(tree, &tree->root) ;  // when tree->leaves==NULL just counts #
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
  int    done = 0, n ;
  double total_entropy ;

  tree->ntraining = ntraining ;
  do
  {
    done = rfTrainNode(rf, tree, &tree->root, training_classes, training_data, rf->ntraining);
    rfFindLeaves(tree) ;
    for (total_entropy = 0.0, n = 0 ; n < tree->nleaves ; n++)
      total_entropy += entropy(tree->leaves[n]->class_counts, rf->nclasses) ;
    printf("training complete, leaf entropy = %2.1f, nleaves = %d, max depth %d\n",
	   total_entropy, tree->nleaves, tree->depth) ;
  } while (!done) ;
  return(NO_ERROR) ;
}
int
RFtrain(RANDOM_FOREST *rf, double feature_overlap, int *training_classes, double **training_data, int ntraining_total)
{
  int    n, i, nfeatures, *feature_permutation, *training_permutation, f, ntraining,index ;
  TREE   *tree ;

  rf->ntraining = ntraining_total ;
  rf->training_data = training_data ;
  rf->training_classes = training_classes ;
  rf->feature_overlap = feature_overlap ;

  for (f = 0 ; f < rf->nfeatures ; f++)
  {
    rf->feature_min[f] = 1e20 ;
    rf->feature_max[f] = -1e20 ;
    for (i = 0 ; i < ntraining_total ; i++)
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
  training_permutation = calloc(ntraining_total, sizeof(int)) ;
  if (training_permutation == NULL)
    ErrorExit(ERROR_NOMEMORY, "could not allocate feature permutation vector");
  for (n = 0 ; n < rf->ntrees ; n++)  // train each tree
  {
    tree = &rf->trees[n] ;
    tree->root.training_set = calloc(rf->ntraining, sizeof(tree->root.training_set[0])) ;
    if (tree->root.training_set == NULL)
      ErrorExit(ERROR_NOMEMORY, "RFtrain(%d): could not allocate root %d training set",
                ntraining_total, n) ;
    // randomize what features this tree will use
    nfeatures = MIN(rf->nfeatures, feature_overlap*rf->nfeatures) ;
    compute_permutation(rf->nfeatures, feature_permutation) ;
    tree->feature_list = (int *)calloc(nfeatures, sizeof(tree->feature_list[0]));
    if (tree->feature_list == NULL)
      ErrorExit(ERROR_NOMEMORY, "RFtrain: could not allocate feature list %d (%d)",
                n,nfeatures) ;
    tree->nfeatures = nfeatures ;
    for (i = 0 ; i < nfeatures ; i++)
      tree->feature_list[i] = feature_permutation[i] ;

    // randomize what training data this tree will use
    ntraining = MIN(ntraining_total, 2.0*ntraining_total / rf->ntrees) ;
    ntraining = ntraining_total ; // disable randomness in training set
    compute_permutation(ntraining_total, training_permutation) ;
    tree->root.total_counts = 0 ;
    for (i = 0 ; i < ntraining ; i++)
    {
      index = training_permutation[i] ;
      index = i ;  // disable randomness in training set
      tree->root.class_counts[training_classes[index]]++ ;
      tree->root.training_set[index] = 1 ;
      tree->root.total_counts++ ;
    }

    printf("tree %d: initial entropy = %f\n", n, 
           entropy(tree->root.class_counts, rf->nclasses)) ;
    rfTrainTree(rf, tree, training_classes, training_data, rf->ntraining) ;
  }

  free(feature_permutation) ;
  free(training_permutation) ;
  return(NO_ERROR) ;
}

static double
entropy(int *class_counts, int nclasses)
{
  int    total, c ;
  double ent, p ;

  for (total = c = 0 ; c < nclasses ; c++)
    total += class_counts[c] ;
  if (total == 0)
    return(0.0) ;
  for (ent = 0.0, c = 0 ; c < nclasses ; c++)
  {
    p = (double)class_counts[c] / total ;
    if (p > 0)
      ent -= p * log(p) ;
  }
  return(ent) ;
}

static int
rfWriteNode(RANDOM_FOREST *rf, NODE *node, FILE *fp)
{
  int c, n ;

  fprintf(fp, "NODE %d %d %d %lf\n", node->left == NULL, node->depth, node->feature, node->thresh) ;
  for (c = 0 ; c < rf->nclasses ; c++)
    fprintf(fp, "%d ", node->class_counts[c]) ;
  fprintf(fp, "\n") ;
  for (n = 0 ; n < rf->ntraining ; n++)
    fprintf(fp, "%d ", node->training_set[n]) ;
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
    fprintf(fp, "%d ", tree->feature_list[f]) ;
  fprintf(fp, "\n") ;
  rfWriteNode(rf, &tree->root, fp) ;
  fprintf(fp, "TREE: END\n") ;
  return(NO_ERROR) ;
}
static int
rfReadNode(RANDOM_FOREST *rf, NODE *node, FILE *fp)
{
  int   c, n, leaf ;
  char  line[MAX_LINE_LEN], *cp ;

  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  sscanf(cp, "NODE %d %d %d %lf\n", &leaf, &node->depth, &node->feature, &node->thresh) ;
  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  cp = strtok(line, " ") ;
  for (c = 0 ; c < rf->nclasses ; c++)
  {
    sscanf(cp, "%d ", &node->class_counts[c]) ;
    cp = strtok(NULL, " ") ;
  }

  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  cp = strtok(line, " ") ;
  for (n = 0 ; n < rf->ntraining ; n++)
  {
    sscanf(cp, "%d ", &node->training_set[n]) ;
    cp = strtok(NULL, " ") ;
  }

  cp = fgetl(line, MAX_LINE_LEN, fp) ;
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
  char  line[10*MAX_LINE_LEN], *cp ;

  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  sscanf(cp, "TREE %d %d %d\n", &tree->depth, &tree->nleaves, &tree->nfeatures) ;
  tree->feature_list = (int *)calloc(tree->nfeatures, sizeof(tree->feature_list[0]));
  if (tree->feature_list == NULL)
    ErrorExit(ERROR_NOMEMORY, "rfReadTree: could not allocate feature list  (%d)",
	      tree->nfeatures) ;
  cp = fgetl(line, 10*MAX_LINE_LEN, fp) ;   // feature list
  if (strlen(cp) >= 10*MAX_LINE_LEN-1)
    ErrorExit(ERROR_UNSUPPORTED, "rfReadTree: max line len exceeded - cannot read tree") ;
  cp = strtok(line, " ") ;
  for (f = 0 ; f < tree->nfeatures ; f++)
  {
    sscanf(cp, "%d ", &tree->feature_list[f]) ;
    cp = strtok(NULL, " ") ;
    if ( f == 836)
      DiagBreak() ;
    if (cp == NULL)
      DiagBreak() ;
  }
  tree->root.training_set = (int *)calloc(rf->ntraining, sizeof(tree->root.training_set[0])) ;
  rfReadNode(rf, &tree->root, fp) ;
  rfFindLeaves(tree) ;
  return(NO_ERROR) ;
}
int
RFwrite(RANDOM_FOREST *rf, char *fname) 
{
  FILE   *fp ;
  int    c, f, n ;

  fp = fopen(fname, "w") ;
  if (fp == NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "RFwrite(%s): could not open file",fname));

  fprintf(fp, "%d %d %d %d %d %d %lf\n", 
          rf->nfeatures, rf->nclasses, rf->ntrees, rf->max_depth, rf->ntraining,rf->nsteps,rf->feature_overlap) ;
  for (c = 0 ; c < rf->nclasses ; c++)
    fprintf(fp, "%s\n", rf->class_names[c]) ;

  for (f = 0 ; f < rf->nfeatures ; f++)
    fprintf(fp, "%lf %lf\n", rf->feature_min[f], rf->feature_max[f]) ;

  for (n = 0 ; n < rf->ntrees ; n++)
    rfWriteTree(rf, &rf->trees[n], fp) ;

  fclose(fp) ;
  return(NO_ERROR) ;
}


static NODE *
rfFindLeaf(RF *rf, NODE *node, double *feature)
{
  if (node->left == NULL)
    return(node) ;
  if (feature[node->feature] < node->thresh)
    return(rfFindLeaf(rf, node->left, feature)) ;
  else
    return(rfFindLeaf(rf, node->right, feature)) ;
}

int
RFclassify(RANDOM_FOREST *rf, double *feature, int true_class)
{
  int max_class = -1, n, c, class_counts[MAX_CLASSES], max_count, total_count ;
  NODE *node ;

  memset(class_counts, 0, sizeof(class_counts)) ;
  for (n = 0 ; n < rf->ntrees ; n++)
  {
    node = rfFindLeaf(rf, &rf->trees[n].root, feature) ;
    for (c = 0 ; c < rf->nclasses ; c++)
      class_counts[c] += node->class_counts[c] ;
    if (DIAG_VERBOSE_ON)
      for (c = 0 ; c < rf->trees[n].nleaves ; c++)  // diagnostic
	if (node == rf->trees[n].leaves[c])
	  printf("tree %d: leaf %d\n", n, c) ;
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
  return(max_class) ;
}

RF *
RFread(char *fname)
{
  RF     *rf ;
  int    nfeatures, nclasses, ntrees, max_depth, ntraining, nsteps, c, n, f ;
  char   line[MAX_LINE_LEN], *cp, *class_names[MAX_CLASSES] ;
  FILE   *fp ;
  double feature_overlap ;

  fp = fopen(fname, "r") ;
  if (fp == NULL)
    ErrorExit(ERROR_NOFILE, "RFread(%s): could not open file", fname) ;

  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  sscanf(cp, "%d %d %d %d %d %d %lf", 
         &nfeatures, &nclasses, &ntrees, &max_depth, &ntraining,&nsteps, &feature_overlap) ;
  for (c = 0 ; c < nclasses ; c++)
  {
    cp = fgetl(line, MAX_LINE_LEN, fp) ;
    class_names[c] = (char *)calloc(strlen(cp)+1, sizeof(char)) ;
    strcpy(class_names[c], cp) ;
  }
  rf = RFalloc(ntrees, nfeatures, nclasses, max_depth, class_names) ;
  rf->feature_overlap = feature_overlap ;
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

