/**
 * @file  mris_divide_parcellation.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/03/15 19:25:08 $
 *    $Revision: 1.7 $
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



/*!
\author Bruce Fischl

*/

/*
BEGINUSAGE --------------------------------------------------------------

mris_divide_parcellation [options] subject hemi sourceannot [splitfile|areathresh] outannot

ENDUSAGE --------------------------------------------------------------

BEGINHELP --------------------------------------------------------------

This program divides one or more parcellations into divisions
perpendicular to the long axis of the label.  The number of divisions
can be specified in one of two ways, depending upon the nature of the
fourth argument.

First, a splitfile can be specified as the fourth argument. The
splitfile is a text file with two columns. The first column is the
name of the label in the source annotation, and the second column is
the number of units that label should be divided into. The names of
the labels depends upon the source parcellation.  For aparc.annot and
aparc.a2005.annot, the names can be found in
$FREESURFER_HOME/FreeSurferColorLUT.txt.  For aparc.annot, the labels
are between the ranges of 1000-1034.  For aparc.a2005s.annot, the
labels are between the ranges of 1100-1181.  The name for the label is
the name of the segmentation without the 'ctx-lh'. Note that the name
included in the splitfile does not indicate the hemisphere. For
example, 1023 is 'ctx-lh-posteriorcingulate'.  You should put
'posteriorcingulate' in the splitfile. Eg, to divide it into three
segments, the following line should appear in the splitfile:

posteriorcingulate 3

Only labels that should be split need be specified in the splitfile.

The second method is to specify an area threshold (in mm^2) as the
fourth argument, in which case each label is divided until each
subdivision is smaller than the area threshold.

The output label name will be the original name with _divN appended,
where N is the division number. N will go from 2 to the number of
divisions. The first division has the same name as the original label.


EXAMPLE:

cd $SUBJECTS_DIR/mysubj001/label

Create a split file to divide the superior frontal gyrus into 4
segements and the precentral gyrus into 3 segments:

echo superiorfrontal 4 >  splittable.txt
echo precentral      3 >> splittable.txt

Run program:

mris_divide_parcellation mysubj001 rh aparc.annot splittable.txt rh.aparc.split-sfg+pc

This reads in rh.aparc.annot, splits the SFG and PC, and and saves the
result to rh.aparc.split-sfg+pc

View

tksurfer mysubj001 rh inflated -annot aparc.split-sfg+pc

The SFG divisions will have the following names: superiorfrontal,
superiorfrontal_div2, superiorfrontal_div3, superiorfrontal_div4. The
PC divisions will be precentral, precentral_div2, precentral_div3.



ENDHELP ----------------------------------------------------------------

*/

// $Id: mris_divide_parcellation.c,v 1.7 2007/03/15 19:25:08 greve Exp $



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "version.h"
#include "utils.h"
#include "mrisurf.h"
#include "error.h"
#include "diag.h"
#include "annotation.h"
#include "matrix.h"

static int get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mris_divide_parcellation.c,v 1.7 2007/03/15 19:25:08 greve Exp $";
char *Progname = NULL;

static char sdir[STRLEN] = "" ;

static int rgb_scale = 30 ;
int  MRISdivideAnnotation(MRI_SURFACE *mris, int *nunits) ;
int  MRISdivideAnnotationUnit(MRI_SURFACE *mris, int annot, int nunits) ;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int         ac, nargs, *nunits, index ;
  MRI_SURFACE *mris ;
  char        **av, fname[STRLEN], *annot_name, *out_fname, *cp,
  *subject, *hemi ;
  float       area_thresh ;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1)  exit (0);
  argc -= nargs;

  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  Progname = argv[0] ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)  usage_exit() ;
  subject = argv[1] ;
  hemi = argv[2] ;
  annot_name = argv[3] ;
  area_thresh = atof(argv[4]) ;
  out_fname = argv[5] ;

  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  sprintf(fname, "%s/%s/surf/%s.sphere", sdir, subject, hemi) ;
  mris = MRISread(fname) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s",
              Progname, fname) ;
  if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read annotation %s",
              Progname, annot_name) ;
  if (mris->ct == NULL)
    ErrorExit(ERROR_NOFILE, "%s: annotation %s must have embedded color table",
              Progname, annot_name) ;

  nunits = (int *)calloc(mris->ct->nentries, sizeof(int)) ;
  if (!nunits)
    ErrorExit(ERROR_BADPARM, "%s: could not allocate %d nunits table", Progname,mris->ct->nentries) ;

  MRIScomputeMetricProperties(mris) ;
  if (isdigit(*argv[4]))  // area threshold specified
  {
    int    vno ;
    VERTEX *v ;
    float  *area ;
    area = (float *)calloc(mris->ct->nentries, sizeof(float)) ;
    if (!area)
      ErrorExit(ERROR_BADPARM, "%s: could not allocate %d area table", Progname,mris->ct->nentries) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      CTABfindAnnotation(mris->ct, v->annotation, &index) ;
      if (index < 0 || v->ripflag)
        continue ;
      area[index] += v->area ;
    }
    for (index = 0 ; index < mris->ct->nentries ; index++)
      nunits[index] =  (int)(area[index] / area_thresh)+1 ;
    free(area) ;
  } else                 // interpret it as a file with parcellation names and # of units
  {
    char  line[STRLEN], *cp, name[STRLEN] ;
    FILE  *fp ;
    int   num ;

    printf("interpreting 4th command line arg as split file name") ;
    fp = fopen(argv[4], "r") ;
    if (fp == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not open parcellation division file %s", Progname,argv[4]) ;
    while ((cp = fgetl(line, STRLEN-1, fp)) != NULL) {
      if (sscanf(line, "%s %d", name, &num) != 2)
        ErrorExit(ERROR_BADFILE, "%s: could not parse name/num pair from '%s'", Progname, line) ;
      CTABfindName(mris->ct, name, &index) ;
      if (index < 0)
        ErrorExit(ERROR_BADFILE, "%s: could not find name '%s' in color table", Progname, name) ;
      nunits[index] = num ;
    }
  }

  for (index = 0 ; index < mris->ct->nentries ; index++)
    if (nunits[index] > 1)
      printf("dividing %s into %d parts\n", mris->ct->entries[index]->name, nunits[index]) ;

  MRISdivideAnnotation(mris, nunits) ;

  free(nunits) ;
  printf("saving annotation to %s\n", out_fname) ;
  MRISwriteAnnotation(mris, out_fname) ;
  return 0;
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
  if (!stricmp(option, "-help") || !stricmp(option, "help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "scale")) {
    rgb_scale = atoi(argv[2]) ;
    printf("scaling rgb values by %d\n", rgb_scale) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      usage_exit() ;
      break ;
    }
  return(nargs) ;
}
static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  printf("\n");
  printf("mris_divide_parcellation [options] subject hemi sourceannot [splitfile|areathresh] outannot\n");
  printf("\n");
  printf("\n");
  printf("options\n");
  printf( "  -help\n");
  printf( "  -scale <scale>   specify offset scaling for rgb values (default=20)\n");
  printf( "  -l <label name>  only process the label <label name> (not implemented yet)\n");
}

static void
print_help(void) {
  print_usage() ;
printf("\n");
printf("This program divides one or more parcellations into divisions\n");
printf("perpendicular to the long axis of the label.  The number of divisions\n");
printf("can be specified in one of two ways, depending upon the nature of the\n");
printf("fourth argument.\n");
printf("\n");
printf("First, a splitfile can be specified as the fourth argument. The\n");
printf("splitfile is a text file with two columns. The first column is the\n");
printf("name of the label in the source annotation, and the second column is\n");
printf("the number of units that label should be divided into. The names of\n");
printf("the labels depends upon the source parcellation.  For aparc.annot and\n");
printf("aparc.a2005.annot, the names can be found in\n");
printf("$FREESURFER_HOME/FreeSurferColorLUT.txt.  For aparc.annot, the labels\n");
printf("are between the ranges of 1000-1034.  For aparc.a2005s.annot, the\n");
printf("labels are between the ranges of 1100-1181.  The name for the label is\n");
printf("the name of the segmentation without the 'ctx-lh'. Note that the name\n");
printf("included in the splitfile does not indicate the hemisphere. For\n");
printf("example, 1023 is 'ctx-lh-posteriorcingulate'.  You should put\n");
printf("'posteriorcingulate' in the splitfile. Eg, to divide it into three\n");
printf("segments, the following line should appear in the splitfile:\n");
printf("\n");
printf("posteriorcingulate 3\n");
printf("\n");
printf("Only labels that should be split need be specified in the splitfile.\n");
printf("\n");
printf("The second method is to specify an area threshold (in mm^2) as the\n");
printf("fourth argument, in which case each label is divided until each\n");
printf("subdivision is smaller than the area threshold.\n");
printf("\n");
printf("The output label name will be the original name with _divN appended,\n");
printf("where N is the division number. N will go from 2 to the number of\n");
printf("divisions. The first division has the same name as the original label.\n");
printf("\n");
printf("\n");
printf("EXAMPLE:\n");
printf("\n");
printf("cd $SUBJECTS_DIR/mysubj001/label\n");
printf("\n");
printf("Create a split file to divide the superior frontal gyrus into 4\n");
printf("segements and the precentral gyrus into 3 segments:\n");
printf("\n");
printf("echo superiorfrontal 4 >  splittable.txt\n");
printf("echo precentral      3 >> splittable.txt\n");
printf("\n");
printf("Run program:\n");
printf("\n");
printf("mris_divide_parcellation mysubj001 rh aparc.annot splittable.txt rh.aparc.split-sfg+pc\n");
printf("\n");
printf("This reads in rh.aparc.annot, splits the SFG and PC, and and saves the\n");
printf("result to rh.aparc.split-sfg+pc\n");
printf("\n");
printf("View\n");
printf("\n");
printf("tksurfer mysubj001 rh inflated -annot aparc.split-sfg+pc\n");
printf("\n");
printf("The SFG divisions will have the following names: superiorfrontal,\n");
printf("superiorfrontal_div2, superiorfrontal_div3, superiorfrontal_div4. The\n");
printf("PC divisions will be precentral, precentral_div2, precentral_div3.\n");
printf("\n");
printf("\n");
printf("\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

int
MRISdivideAnnotation(MRI_SURFACE *mris, int *nunits) {
  int   *done, vno, index, nadded, i, num, j, annot, new_annot ;
  VERTEX *v ;
  COLOR_TABLE *ct ;

  MRIScomputeMetricProperties(mris) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  done = (int *)calloc(mris->ct->nentries, sizeof(int)) ;
  if (done == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d index table",
              Progname, mris->ct->nentries) ;

  MRISclearMarks(mris) ;
  MRISsetNeighborhoodSize(mris, 2) ;
  for (nadded = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation <= 0)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    CTABfindAnnotation(mris->ct, v->annotation, &index) ;
    if (index <= 0 || done[index])  // don't do unknown (index = 0)
      continue ;
    if (index == Gdiag_no)
      DiagBreak() ;
#if 0
    if (stricmp("postcentral", mris->ct->entries[index]->name))
      continue ;
#endif
    num = MRISdivideAnnotationUnit(mris, v->annotation, nunits[index]) ;
    nadded += num ;
    done[index] = 1+num ;
  }

  printf("allocating new colortable with %d additional units...\n", nadded) ;
  ct = CTABalloc(mris->ct->nentries+nadded) ;
  index = mris->ct->nentries ;
  for (i = 0 ; i < mris->ct->nentries ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    *(ct->entries[i]) = *(mris->ct->entries[i]) ;
    for (j = 1 ; j < done[i] ; j++) {
      int offset, new_index, ri, gi, bi, found ;

      *(ct->entries[index]) = *(mris->ct->entries[i]) ;
      sprintf(ct->entries[index]->name, "%s_div%d", ct->entries[i]->name, j+1) ;
      offset = j ;
      found = 0 ;
      do {
        ri = (ct->entries[i]->ri+rgb_scale*offset) % 256 ;
        gi = (ct->entries[i]->gi+rgb_scale*offset) % 256 ;
        bi = (ct->entries[i]->bi+rgb_scale*offset) % 256 ;
        CTABfindRGBi(ct, ri, gi, bi, &new_index) ;
        if (new_index < 0)  // couldn't find this r,g,b set - can use it for new entry
        {
          ct->entries[index]->ri = ri ;
          ct->entries[index]->gi = gi ;
          ct->entries[index]->bi = bi ;

          ct->entries[index]->rf = (float)ri/255.0f;
          ct->entries[index]->gf = (float)gi/255.0f;
          ct->entries[index]->bf = (float)bi/255.0f;
          found = 1 ;

          CTABannotationAtIndex(ct, i, &annot) ;
          CTABannotationAtIndex(ct, index, &new_annot) ;
          // translate old annotations to new ones
          for (vno = 0 ; vno < mris->nvertices ; vno++) {
            v = &mris->vertices[vno] ;
            if (v->ripflag || v->marked != j || v->annotation != annot)
              continue ;
            v->annotation = new_annot ;
          }
        }
        else {
          offset++ ;
          found = 0 ;
        }
      } while (!found) ;
      index++ ;
    }
  }
  CTABfree(&mris->ct) ;
  mris->ct = ct ;
  return(NO_ERROR) ;
}

/*
  will split up parcellation units based on area and write the index of the
  new unit into the marked field (marked=0->old annot, marked=1-> 1st new subdivision, etc..)
  and will also put the continuous measure of how far along the eigendirection the
  vertex is into v->curv (in [0 1]).

  return the # of additional parcellation units that have been added
*/
int
MRISdivideAnnotationUnit(MRI_SURFACE *mris, int annot, int nunits) {
  int    vno, num, min_vno ;
  VERTEX *v, *vc ;
  float  cx, cy, cz, dist, min_dist, evalues[3], u, w, dx, dy, dz,
  min_dot, max_dot, e1x, e1y, e1z, dot, mx ;
  MATRIX *m_obs, *m_obs_T, *m_cov, *m_eig ;


  if (nunits < 2)
    return(0);

  // compute centroid of annotation
  cx = cy = cz = 0 ;
  for (num = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    cx += v->x ;
    cy += v->y ;
    cz += v->z ;
    num++ ;
  }
  if (num == 0) // unused parcellation
    return(0) ;

  cx /= num ;
  cy /= num ;
  cz /= num ;

  // find vertex in annotation closest to centroid
  min_dist = 100000 ;
  min_vno = -1 ;
  vc = NULL ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    dist = sqrt(SQR(v->x-cx)+SQR(v->y-cy)+SQR(v->z-cz));
    if (dist < min_dist) {
      min_vno = vno;
      min_dist = dist ;
      vc = v ;
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("using v %d as closest (%2.3f mm) from centroid (%2.1f, %2.1f, %2.1f)\n",
           min_vno, min_dist, cx, cy, cz) ;

  // now compute eigensystem around this vertex
  m_obs = MatrixAlloc(2, num, MATRIX_REAL) ;
  for (num = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    dx = v->x - cx ;
    dy = v->y - cy ;
    dz = v->z - cz ;
    u = vc->e1x * dx + vc->e1y * dy + vc->e1z*dz ;
    w = vc->e2x * dx + vc->e2y * dy + vc->e2z*dz ;
    *MATRIX_RELT(m_obs, 1, num+1) = u ;
    *MATRIX_RELT(m_obs, 2, num+1) = w ;
    num++ ;
  }

  m_obs_T = MatrixTranspose(m_obs, NULL) ;
  m_cov = MatrixMultiply(m_obs,m_obs_T, NULL) ;
  m_eig = MatrixEigenSystem(m_cov, evalues, NULL) ;
  e1x = *MATRIX_RELT(m_eig, 1,1) * vc->e1x + *MATRIX_RELT(m_eig, 2,1) * vc->e2x;
  e1y = *MATRIX_RELT(m_eig, 1,1) * vc->e1y + *MATRIX_RELT(m_eig, 2,1) * vc->e2y;
  e1z = *MATRIX_RELT(m_eig, 1,1) * vc->e1z + *MATRIX_RELT(m_eig, 2,1) * vc->e2z;
  if (fabs(e1x) > fabs(e1y) &&  fabs(e1x) > fabs(e1z))
    mx = e1x ;
  else if (fabs(e1y) > fabs(e1z))
    mx = e1y ;
  else
    mx = e1z ;
  //  if (mx < 0)
  if (e1y < 0)  // orient them from posterior to anterior
  {
    e1x *= -1 ;
    e1y *= -1 ;
    e1z *= -1 ;
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    MatrixPrint(stdout, m_eig) ;
  dist = sqrt(e1x*e1x + e1y*e1y + e1z*e1z) ;
  e1x /= dist ;
  e1y /= dist ;
  e1z /= dist ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("principle eigendirection = (%2.3f, %2.3f, %2.3f)\n",e1x,e1y,e1z) ;
  min_dot = 10000 ;
  max_dot = -min_dot ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    dx = v->x - cx ;
    dy = v->y - cy ;
    dz = v->z - cz ;
    dot = dx*e1x + dy*e1y + dz*e1z ;
    if (dot > max_dot)
      max_dot = dot ;
    if (dot < min_dot)
      min_dot = dot ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    if  (vno == Gdiag_no)
      DiagBreak() ;
    dx = v->x - cx ;
    dy = v->y - cy ;
    dz = v->z - cz ;
    dot = dx*e1x + dy*e1y + dz*e1z ;
    v->curv = (dot - min_dot) / (max_dot-min_dot) ;
    v->marked = (int)(nunits * v->curv) ;
    if (v->marked >= nunits)
      v->marked = nunits-1 ;  // one with max_dot will be just too big
  }

  MatrixFree(&m_obs) ;
  MatrixFree(&m_cov) ;
  MatrixFree(&m_obs_T) ;
  return(nunits-1) ;
}

