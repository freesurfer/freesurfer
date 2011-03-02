/**
 * @file  rbftest.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:35 $
 *    $Revision: 1.7 $
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
#include <ctype.h>
#include <memory.h>

#include "diag.h"
#include "error.h"
#include "macros.h"
#include "utils.h"
#include "proto.h"
#include "cluster.h"
#include "matrix.h"
#include "rbf.h"
#include "version.h"

static int  verbose = 0 ;

char *Progname ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
/*static int readAllData(FILE *fp, CLUSTER_SET *cs) ;*/
static int read_line(VECTOR *v_obs, int obs_no, void *vfp, int same_class,
                     int *pclass) ;
static void classify_all(RBF *rbf, FILE *fp) ;
static int  count_classes(FILE *fp) ;

#define NINPUTS             2
#define NCLUSTERS           2
#define NORMALIZE_CLUSTERS  0
#define RBF_FNAME           "test.rbf"

#define MAX_CLASSES            10
#define NCLASSES               2

static float momentum = 0.0f ;
static int nclusters = 0 ;
static int nclasses = 0 ;

static char *class_names[MAX_CLASSES] = { "gray", "white" , "two", "three",
                                        "four"
                                        } ;
static int max_clusters[MAX_CLASSES] = {
                                         NCLUSTERS, NCLUSTERS, NCLUSTERS
                                       } ;

int
main(int argc, char *argv[]) {
  char         *input_file_name ;
  int          nargs, class ;
  FILE         *fp ;
  /*  CLUSTER_SET  *cs ;*/
  RBF          *rbf ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: rbftest.c,v 1.7 2011/03/02 00:04:35 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    input_file_name = "cluster.dat" ;
  else
    input_file_name = argv[1] ;

  fp = fopen(input_file_name, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not read input from %s",
              Progname, input_file_name) ;

  if (!nclasses)
    nclasses = count_classes(fp) ;

  if (nclusters > 0)
    for (class = 0 ; class < nclasses ; class++)
      max_clusters[class] = nclusters ;
  rbf = RBFinit(NINPUTS, nclasses, max_clusters, class_names) ;
  if (RBFtrain(rbf, read_line, (void *)fp, momentum) != NO_ERROR)
    exit(1) ;
  RBFwrite(rbf, RBF_FNAME) ;
  RBFfree(&rbf) ;
  rbf = RBFread(RBF_FNAME) ;
  classify_all(rbf, fp) ;
  RBFfree(&rbf) ;
  fclose(fp) ;

  exit(0) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  switch (toupper(*option)) {
  case 'C':
    sscanf(argv[2], "%d", &nclasses) ;
    nargs = 1 ;
    break ;
  case 'M':
    sscanf(argv[2], "%f", &momentum) ;
    nargs = 1 ;
    break ;
  case 'N':
    sscanf(argv[2], "%d", &nclusters) ;
    nargs = 1 ;
    break ;
  case 'V':
    verbose = !verbose ;
    break ;
  case '?':
  case 'U':
    fprintf(stderr,
            "usage: %s <classifier file> <input volume> <output volume>\n",
            Progname) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
#if 0
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
readAllData(FILE *fp, CLUSTER_SET *cs) {
  char *cp, line[200] ;
  VECTOR *v_obs ;


  rewind(fp) ;
  v_obs = VectorAlloc(cs->ninputs, MATRIX_REAL) ;

  while ((cp = fgetl(line, 199, fp)) != NULL) {
    sscanf(cp, "%f %f", &VECTOR_ELT(v_obs,1), &VECTOR_ELT(v_obs,2)) ;
    MatrixPrint(stderr, MatrixTranspose(v_obs,NULL)) ;
    CSnewObservation(cs, v_obs) ;
  }

  VectorFree(&v_obs) ;
  return(NO_ERROR) ;
}
#endif
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
read_line(VECTOR *v_obs, int obs_no, void *vfp, int same_class,int *pclass) {
  char *cp, line[200] ;
  FILE *fp ;
  int  n, class ;

  fp = (FILE *)vfp ;

  /* find the obs_no th line */
  rewind(fp) ;
  for (n = 0 ; n <= obs_no ; n++) {
    cp = fgetl(line, 199, fp) ;
    if (!cp)
      return(-1) ;
    sscanf(cp, "%f %f %d", &VECTOR_ELT(v_obs,1), &VECTOR_ELT(v_obs,2), &class);
    /* if same_class is set, only count those samples that are in the class
       specified in *pclass. This lets the caller pick the nth sample from
       a given class.
       */
    if (same_class && (class != *pclass))
      n-- ;
  }


  *pclass = class ;
  return(NO_ERROR) ;
}
static void
classify_all(RBF *rbf, FILE *fp) {
  int class, target_class, obs_no ;
  VECTOR *v_obs, *v_error ;
  float  sse, error, rms ;

  v_obs = VectorAlloc(NINPUTS, MATRIX_REAL) ;
  v_error = VectorAlloc(nclasses, MATRIX_REAL) ;

  fprintf(stderr, "classifying inputs...\n") ;
  obs_no = 0 ;
  sse = 0.0f ;
  while (read_line(v_obs, obs_no++, (void *)fp, 0, &target_class) == NO_ERROR) {
    class = RBFclassify(rbf, v_obs) ;
    fprintf(stderr, "%+2.3f %+2.3f --> class %d (%2.3f)%s\n",
            VECTOR_ELT(v_obs,1), VECTOR_ELT(v_obs,2), class,
            RVECTOR_ELT(rbf->v_outputs, class+1),
            class == target_class ? "" : " ERROR") ;
    error = RBFcomputeErrors(rbf, target_class, v_error) ;
    sse += error ;
    RBFprintActivations(rbf, v_obs, v_error, target_class, stderr) ;
  }
  rms = sqrt(sse / (float)obs_no) ;
  fprintf(stderr, "rms = %2.5f\n", rms) ;
  VectorFree(&v_obs) ;
  VectorFree(&v_error) ;
}

static int
count_classes(FILE *fp) {
  VECTOR *v_obs ;
  int    class, classes[100], obs_no = 0, nclasses = 0 ;

  memset(classes, 0, sizeof(classes)) ;
  v_obs = VectorAlloc(NINPUTS, MATRIX_REAL) ;

  while (read_line(v_obs, obs_no++, (void*)fp, 0, &class) == NO_ERROR) {
    if (!classes[class]) {
      classes[class] = 1 ;
      nclasses++ ;
    }
  }
  VectorFree(&v_obs) ;
  return(nclasses) ;
}

