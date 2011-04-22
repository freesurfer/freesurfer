/**
 * @file  mris_segment.c
 * @brief program for segmenting cortical areas based on connectivity/correlation
 *  profiles.
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/04/22 13:47:30 $
 *    $Revision: 1.3 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "mrisurf.h"
#include "label.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

static char cormat_name[STRLEN] = "cormat.mgz" ;
static char label_name[STRLEN] = "MT.fsaverage5.label" ;
static char prior_name[STRLEN] = "MT.label" ;

static int ico_no = 5 ;
static char sdir[STRLEN] = "" ;

static int gaussian_classifier = 1 ;
static char *hemi_name = "lh" ;

static double cor_thresh = 0.6 ;
static double label_thresh = .3 ;
static int nerodes = 0 ;
static double logodds_thresh = .001 ;
static double logodds_slope = .1 ;

#define MAX_SUBJECTS 500

VECTOR *
VectorFromMRIcol(MRI *mri_cmat, VECTOR *v, int col, int frame)
{
  int row ;

  if (v == NULL)
    v = VectorAlloc(mri_cmat->height, MATRIX_REAL) ;

  for (row = 0 ; row < mri_cmat->height ; row++)
    VECTOR_ELT(v, row+1) = MRIgetVoxVal(mri_cmat, col, row, 0, frame) ;
  return(v) ;
}

MATRIX *
MatrixFromMRI(MRI *mri_cmat, MATRIX *cmat, int frame)
{
  int r, c ;

  if (cmat == NULL)
    cmat = MatrixAlloc(mri_cmat->height, mri_cmat->width, MATRIX_REAL) ;

  for (r = 0 ; r < mri_cmat->height ; r++)
    for (c = 0 ; c < mri_cmat->width ; c++)
      *MATRIX_RELT(cmat, r+1, c+1) = (float)MRIgetVoxVal(mri_cmat, c, r, 0, frame) ;

  return(cmat) ;
}

MRI *
MatrixToMRI(MATRIX *cmat, MRI *mri_cmat, int frame)
{
  int r, c ;

  if (mri_cmat == NULL)
    mri_cmat = MRIalloc(cmat->cols, cmat->rows, 1, MATRIX_REAL) ;


  for (r = 0 ; r < mri_cmat->height ; r++)
    for (c = 0 ; c < mri_cmat->width ; c++)
      MRIsetVoxVal(mri_cmat, c, r, 0, frame, *MATRIX_RELT(cmat, r+1,c +1)) ;

  return(mri_cmat) ;
}

VECTOR *
MRIcmatDotProductFrames(MRI *mri1, int frame1, MRI *mri2, int frame2, VECTOR *v_dot)
{
  int    r1, c1 ;
  double dot ;

  if (v_dot == NULL)
    v_dot = VectorAlloc(mri1->height, MATRIX_REAL) ;

  for (r1 = 0 ; r1 < mri1->height ; r1++)
  {
    for (dot = 0.0, c1 = 0 ; c1 < mri1->width ; c1++)
    {
      dot += 
        MRIgetVoxVal(mri1, c1, r1, 0, 0) *
        MRIgetVoxVal(mri2, r1, c1, 0, 0) ;
        
    }
    VECTOR_ELT(v_dot, r1+1) = dot ;
  }
  return(v_dot) ;
}
int
MRISlogOdds(MRI_SURFACE *mris, LABEL *area, double slope) 
{
  int    vno ;
  VERTEX *v ;
  double p ;

  MRISdistanceTransform(mris, area, DTRANS_MODE_SIGNED) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    p = v->val ;
    if (p < 0)
      p = 0 ;
    p = exp(-p*slope) ;
    v->val = p ;
  }
  return(NO_ERROR) ;
}

int
MRIcmatNormalizeRows(MRI *mri_cmat)
{
  int   row, col, frame ;
  double norm ;

  for (frame = 0 ; frame < mri_cmat->nframes ; frame++)
  {
    for (row = 0 ; row < mri_cmat->height ; row++)
    {
      for (norm = 0.0, col = 0 ; col < mri_cmat->width ; col++)
        norm += MRIgetVoxVal(mri_cmat, col, row, 0, frame) ;
      if (DZERO(norm))
        norm = 1.0 ;
      for (col = 0 ; col < mri_cmat->width ; col++)
        MRIsetVoxVal(mri_cmat, col, row, 0, frame, MRIgetVoxVal(mri_cmat, col, row, 0, frame)/norm) ;
    }
  }
  return(NO_ERROR) ;
}
int
MatrixMaxRowIndex(MATRIX *m, int row)
{
  int    c, max_c ;
  double max_val ;

  max_val = *MATRIX_RELT(m, row, 1) ; max_c = 1 ;
  for (c = 2 ; c <= m->cols ; c++)
    if (*MATRIX_RELT(m, row, c) > max_val)
    {
      max_val = *MATRIX_RELT(m, row, c) ;
      max_c = c ;
    }
  return(max_c) ;
}

VECTOR *
MatrixRowNorm(MATRIX *cmat, VECTOR *v_norm)
{
  double   norm, val ;
  int      r, c ;

  if (v_norm == NULL)
    v_norm = VectorAlloc(cmat->cols, MATRIX_REAL) ;

  for (c = 1 ; c <= cmat->cols ; c++)
  {
    for (norm = 0.0, r = 1 ; r <= cmat->rows ; r++)
    {
      val = *MATRIX_RELT(cmat, r, c) ;
      norm += val*val ;
    }
    VECTOR_ELT(v_norm, c) = sqrt(norm) ;
  }
  return(v_norm) ;
}

/*!
  \fn double VectorLogLikelihood(VECTOR *v, VECTOR *v_mean, *v_var)
  \brief compute average of (v-v_mean)^2 / (2*v_var)
*/
double
VectorLogLikelihood(VECTOR *v, VECTOR *v_mean, VECTOR *v_var)
{
  double ll, val, mean, var, total_ll ;
  int    row ;

  for (total_ll = 0.0, row = 1 ; row <= v->rows ; row++)
  {
    val = VECTOR_ELT(v, row) ;
    mean = VECTOR_ELT(v_mean, row) ;
    var = VECTOR_ELT(v_var, row) ;
    if (DZERO(var))
      var = 1.0 ;
    ll = SQR(val-mean) / (2*var) ;
    if (!finite(ll) || !finite(total_ll))
      DiagBreak() ;
    total_ll += ll ;
  }
  return(-total_ll/v->rows) ;
}

static MRI *
classify_subject(MRI_SURFACE *mris, LABEL *prior_label, MRI *mri_cmat, 
                 LABEL **labels, int nsubjects, double thresh) 
{
  int    start_index, end_index, sno, vno, ind, vno2, nvertices, nevals ;
  MRI    *mri_out ;
  MATRIX *m_train ;
  double dot, max_dot, logp, var, val, val2, ll_in, ll_out ;
  VECTOR *v_test, *v_in_vars, *v_in_means, *v_out_means, *v_out_vars ;
  LABEL  *area ;
  VERTEX *v ;

  mri_out = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT) ;

  if (mris->hemisphere == LEFT_HEMISPHERE)
  {
    start_index = 0 ; end_index = mris->nvertices-1 ;
  }
  else   // processing rh
  {
    start_index = mris->nvertices ; end_index = 2*mris->nvertices-1 ;
  }

  MRIthreshold(mri_cmat, mri_cmat, cor_thresh) ;
  MRIcmatNormalizeRows(mri_cmat) ;

  sno = 0 ;
  LabelMark(labels[sno], mris) ;

  m_train = MatrixFromMRI(mri_cmat, NULL, sno) ;

  // compute variance within training MT and assume homoskedasticity
  v_in_vars = VectorAlloc(m_train->rows, MATRIX_REAL) ;
  v_in_means = VectorAlloc(m_train->rows, MATRIX_REAL) ;
  for (ind = 0 ; ind < labels[sno]->n_points ; ind++)
  {
    vno = labels[sno]->lv[ind].vno ;
    for (vno2 = 0 ; vno2 < mri_cmat->height ; vno2++)
    {
      val = MRIgetVoxVal(mri_cmat, vno, vno2, 0, 0) ;
      VECTOR_ELT(v_in_vars, vno2+1) += val*val ;
      VECTOR_ELT(v_in_means, vno2+1) += val ;
    }
  }
  for (vno2 = 0 ; vno2 < mri_cmat->height ; vno2++)
  {
    val = VECTOR_ELT(v_in_means, vno2+1) ;
    val /= labels[sno]->n_points ;  // mean
    val2 = VECTOR_ELT(v_in_vars, vno2+1) ;
    val2 = (val2 / labels[sno]->n_points - (val*val)) ; // variance
    VECTOR_ELT(v_in_means, vno2+1) = val ;
    VECTOR_ELT(v_in_vars, vno2+1) = val2 ;
  }
  // compute mean and variance vectors in the region close to but outside MT
  v_out_vars = VectorAlloc(m_train->rows, MATRIX_REAL) ;
  v_out_means = VectorAlloc(m_train->rows, MATRIX_REAL) ;
  MRISdistanceTransform(mris, labels[sno], DTRANS_MODE_SIGNED) ;
  for (nevals = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->val > 1 && v->val < 10) // near MT, but not in it
    {
      nevals++ ;
      for (vno2 = 0 ; vno2 < mri_cmat->height ; vno2++)
      {
        val = MRIgetVoxVal(mri_cmat, vno, vno2, 0, 0) ;
        VECTOR_ELT(v_out_vars, vno2+1) += val*val ;
        VECTOR_ELT(v_out_means, vno2+1) += val ;
      }
    }
  }
  for (vno2 = 0 ; vno2 < mri_cmat->height ; vno2++)
  {
    val = VECTOR_ELT(v_out_means, vno2+1) ;
    val /= nevals ;  // mean
    val2 = VECTOR_ELT(v_out_vars, vno2+1) ;
    val2 = (val2 / nevals - (val*val)) ; // variance
    VECTOR_ELT(v_out_means, vno2+1) = val ;
    VECTOR_ELT(v_out_vars, vno2+1) = val2 ;
  }

  printf("classifying vertices...\n") ;
  // build log-odds priors
  LabelMarkWithThreshold(prior_label, mris, label_thresh) ;
  MRIScloseMarked(mris, 1) ;
  if (nerodes > 0)
    MRISerodeMarked(mris, nerodes) ;
  area = LabelFromMarkedSurface(mris) ;
  MRISlogOdds(mris, area, logodds_slope) ;
  LabelFree(&area) ;

  for (nvertices = 0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    mris->vertices[vno].ripflag = 1 ;
    for (ind = 0 ; ind < mri_cmat->width ; ind++)
      if (!DZERO(*MATRIX_RELT(m_train, vno+1, ind+1)) && 
          !DEQUAL(*MATRIX_RELT(m_train, vno+1, ind+1),1))
      {
        mris->vertices[vno].ripflag = 0 ;
        nvertices++ ;
        break ;
      }
  }
  printf("processing %d vertices in the FOV (%2.1f%%)\n",
         nvertices, 100.0*nvertices/(float)mris->nvertices) ;

  for (nevals = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (!(vno % 100))
    {
      printf("%d of %d complete (%d considered)\r", 
             vno-start_index, end_index-start_index, nevals) ;
      fflush(stdout) ;
    }
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < logodds_thresh)
      continue ;
    nevals++ ;
    logp = log(v->val) ;

    // extract correlation pattern for test subject at this vertex
    v_test = VectorFromMRIcol(mri_cmat, NULL, vno, nsubjects-1) ;

    // use gaussian classifier for in and out of area
    if (gaussian_classifier)
    {
      ll_in = VectorLogLikelihood(v_test, v_in_means, v_in_vars) ;
      ll_out = VectorLogLikelihood(v_test, v_out_means, v_out_vars) ;
      v->val2 = ll_in ;
      v->val2bak = ll_out ;
      if (vno == Gdiag_no)
        printf("v %d: ll_in = %2.4f, ll_out = %2.4f, prior = %2.4f\n", ll_in, ll_out, v->val) ;
      ll_in += log(v->val) ;
      if (FEQUAL(v->val, 1))
        ll_out += log(1e-10) ;
      else
        ll_out += log(1-v->val) ;
      v->imag_val = ll_in ;
      v->stat = ll_out ;
      val = exp(ll_in) / (exp(ll_in) + exp(ll_out)) ;
      MRIsetVoxVal(mri_out, vno, 0, 0, 0, val) ;
    }
    else
    {
      // compute dot product with every vertex in training subject
      // and find max
      max_dot = -1e10 ; ind = -1 ;
      for (vno2 = 0 ; vno2 < mris->nvertices ; vno2++)
      {
        dot = MatrixRowDotProduct(m_train, vno2+1, v_test) ;
        if (dot > max_dot)
        {
          max_dot = dot ;
          ind = vno2 ;
        }
      }
      
      // if this index is in the training MT, set it to 1, otherwise 0
      if (mris->vertices[ind].marked)
        MRIsetVoxVal(mri_out, vno, 0, 0, 0, 1) ;
      // force mapping to be single-valued. This is order-dependent!
      MatrixSetRegion(m_train, m_train, ind+1, 1, 1, m_train->cols, 0) ;
    }
    VectorFree(&v_test) ;
  }
  printf("\n") ;
  MatrixFree(&m_train) ;

  return(mri_out) ;
}
int
main(int argc, char *argv[])
{
  char         **av, *out_fname, *subject, fname[STRLEN], *cp ;
  int          ac, nargs, sno ;
  int          msec, minutes, seconds, nsubjects ;
  struct timeb start ;
  MRI_SURFACE  *mris ;
  MRI          *mri_frame, *mri, *mri_out ;
  LABEL        *labels[MAX_SUBJECTS], *prior_label ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mris_segment.c,v 1.3 2011/04/22 13:47:30 fischl Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
  {
    usage_exit(1) ;
  }

  if (strlen(sdir) == 0)
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (cp == NULL)
      ErrorExit(ERROR_UNSUPPORTED, "%s: must define SUBJECTS_DIR in env",
                Progname) ;
    strcpy(sdir, cp) ;
  }
  out_fname = argv[argc-1] ;
  nsubjects = argc-2 ;
  printf("processing %d subjects and writing output to %s\n", 
         nsubjects,out_fname) ;

  mri = NULL ; // get rid of compiler warning
  for (sno = 0 ; sno < nsubjects ; sno++)
  {
    subject = argv[sno+1] ;
    printf("processing subject %s, %d of %d\n", subject, sno+1, nsubjects) ;
    sprintf(fname, "%s/%s/fmri/%s", sdir, subject, cormat_name) ;
    mri_frame = MRIread(fname) ;
    if (mri_frame == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load correlation matrix from %s", 
                Progname,fname) ;
    if (sno == 0)
    {
      mri = MRIallocSequence(mri_frame->width, mri_frame->height, 
                             mri_frame->depth, mri_frame->type, nsubjects) ;
      if (mri == NULL)
        ErrorExit(ERROR_NOMEMORY, 
                  "%s: could not allocate (%d x %d x %d x %d) array",
                  mri_frame->width, mri_frame->height, mri_frame->depth,
                  nsubjects) ;
    }
    else
      if (mri_frame->width != mri->width ||
          mri_frame->height != mri->height ||
          mri_frame->depth != mri->depth)
        ErrorExit(ERROR_BADPARM, 
                  "%s: volume %s has incompatible dimensions (%d x %d x %d)",
                  mri_frame->width, mri_frame->height, mri_frame->depth) ;

    MRIcopyFrame(mri_frame, mri, 0, sno) ;
    MRIfree(&mri_frame) ;

    sprintf(fname, "%s/%s/label/%s.%s", sdir, subject, hemi_name,label_name) ;
    labels[sno] = LabelRead(NULL, fname) ;
    if (labels[sno] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load label from %s", 
                Progname,fname) ;
  }
  sprintf(fname, "%s/fsaverage%d/surf/%s.inflated", sdir, ico_no, hemi_name) ;
  mris = MRISread(fname) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load surface %s", Progname,fname) ;
  sprintf(fname, "%s/fsaverage%d/label/%s.%s", 
          sdir, ico_no, hemi_name, prior_name) ;
  prior_label = LabelRead(NULL, fname) ;
  if (prior_label == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load prior label %s",Progname,fname);
  mri_out = classify_subject(mris, prior_label, mri, labels, nsubjects, cor_thresh);
  printf("writing output to %s\n", out_fname) ;
  MRIwrite(mri_out, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "segmentation took %d minutes and %d seconds.\n",
          minutes,seconds) ;
  exit(0) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "cmat"))
  {
    nargs = 1 ;
    strcpy(cormat_name, argv[2]) ;
    printf("using fmri/%s as name of correlation matrices\n", cormat_name) ;
  }
  else if (!stricmp(option, "prior"))
  {
    strcpy(prior_name, argv[2]) ;
    printf("using label/%s as name of prior label\n", prior_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "sdir"))
  {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "label"))
  {
    strcpy(label_name, argv[2]) ;
    printf("using label/%s as name of task-based training label\n", label_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "lthresh"))
  {
    label_thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("using label thresh = %2.3f\n", label_thresh);
  }
  else if (!stricmp(option, "erode"))
  {
    nerodes = atoi(argv[2]) ;
    nargs = 1 ;
    printf("eroding prior label %d times\n", nerodes);
  }
  else switch (toupper(*option))
  {
  case 'T':
    cor_thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("using correlation threshold %2.3f\n", cor_thresh) ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    printf("debugging vertex %d\n", Gdiag_no) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code)
{
  printf("usage: %s [options] <subject 1> ... <output subject> <output file>\n",
         Progname) ;
  exit(code) ;
}





