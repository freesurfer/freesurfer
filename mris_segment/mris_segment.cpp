/**
 * @brief segments cortical areas based on connectivity/correlation/intensity profiles
 *
 */
/*
 * Original Author: Bruce Fischl
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
#include "mri2.h"
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

const char *Progname ;
static void usage_exit(int code) ;

static char data_name[STRLEN] = "cormat.mgz" ;
static char label_name[STRLEN] = "MT.fsaverage5.label" ;
static char prior_name[STRLEN] = "invivo.MT.logodds.mgz" ;

static int ico_no = 5 ;
static char sdir[STRLEN] = "" ;

#define CLASSIFY_GAUSSIAN      0
#define CLASSIFY_SIMILARITY    1
#define CLASSIFY_LABEL_FUSION  2

static int prior_number_of_vertices = 30 ;
static int nsmooth = 1 ;
static int nclose = 1 ;
static int classifier = CLASSIFY_GAUSSIAN ;
static const char *hemi_name = "lh" ;

static int prior_only = 0 ;
static double cor_thresh = 0.6 ;
static double logodds_thresh = .001 ;
static const char *data_dir = "fmri" ;

#define MAX_SUBJECTS 1500

VECTOR *
VectorFromMRIcol(MRI *mri_cmat, VECTOR *v, int col, int frame)
{
  int row ;

  if (v == NULL)
  {
    v = VectorAlloc(mri_cmat->height, MATRIX_REAL) ;
  }

  for (row = 0 ; row < mri_cmat->height ; row++)
  {
    VECTOR_ELT(v, row+1) = MRIgetVoxVal(mri_cmat, col, row, 0, frame) ;
  }
  return(v) ;
}

MATRIX *
MatrixFromMRI(MRI *mri_cmat, MATRIX *cmat, int frame)
{
  int r, c ;

  if (cmat == NULL)
  {
    cmat = MatrixAlloc(mri_cmat->height, mri_cmat->width, MATRIX_REAL) ;
  }

  for (r = 0 ; r < mri_cmat->height ; r++)
    for (c = 0 ; c < mri_cmat->width ; c++)
    {
      *MATRIX_RELT(cmat, r+1, c+1) = (float)MRIgetVoxVal(mri_cmat, c, r, 0, frame) ;
    }

  return(cmat) ;
}

MRI *
MatrixToMRI(MATRIX *cmat, MRI *mri_cmat, int frame)
{
  int r, c ;

  if (mri_cmat == NULL)
  {
    mri_cmat = MRIalloc(cmat->cols, cmat->rows, 1, MATRIX_REAL) ;
  }


  for (r = 0 ; r < mri_cmat->height ; r++)
    for (c = 0 ; c < mri_cmat->width ; c++)
    {
      MRIsetVoxVal(mri_cmat, c, r, 0, frame, *MATRIX_RELT(cmat, r+1,c +1)) ;
    }

  return(mri_cmat) ;
}

VECTOR *
MRIcmatDotProductFrames(MRI *mri1, int frame1, MRI *mri2, int frame2, VECTOR *v_dot)
{
  int    r1, c1 ;
  double dot ;

  if (v_dot == NULL)
  {
    v_dot = VectorAlloc(mri1->height, MATRIX_REAL) ;
  }

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
MRIcmatNormalizeRows(MRI *mri_cmat)
{
  int   row, col, frame ;
  double norm ;

  for (frame = 0 ; frame < mri_cmat->nframes ; frame++)
  {
    for (row = 0 ; row < mri_cmat->height ; row++)
    {
      for (norm = 0.0, col = 0 ; col < mri_cmat->width ; col++)
      {
        norm += MRIgetVoxVal(mri_cmat, col, row, 0, frame) ;
      }
      if (DZERO(norm))
      {
        norm = 1.0 ;
      }
      for (col = 0 ; col < mri_cmat->width ; col++)
      {
        MRIsetVoxVal(mri_cmat, col, row, 0, frame, MRIgetVoxVal(mri_cmat, col, row, 0, frame)/norm) ;
      }
    }
  }
  return(NO_ERROR) ;
}
int
MatrixMaxRowIndex(MATRIX *m, int row)
{
  int    c, max_c ;
  double max_val ;

  max_val = *MATRIX_RELT(m, row, 1) ;
  max_c = 1 ;
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
  {
    v_norm = VectorAlloc(cmat->cols, MATRIX_REAL) ;
  }

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
    {
      var = 1.0 ;
    }
    ll = SQR(val-mean) / (2*var) ;
    if (!std::isfinite(ll) || !std::isfinite(total_ll))
    {
      DiagBreak() ;
    }
    total_ll += ll ;
  }
  return(-total_ll/v->rows) ;
}

static MRI *
classify_vertices(MRI_SURFACE *mris, MRI *mri_prior, MRI *mri_cmat,
                  LABEL **labels, int nsubjects, double thresh, int prior_only)
{
  int    start_index, end_index, sno, vno, ind, vno2, nvertices, nevals,
         nin, nout, *in_label[MAX_SUBJECTS] ;
  MRI    *mri_out ;
  MATRIX *m_train = NULL, *m_trains[MAX_SUBJECTS] ;
  double dot, max_dot, val=0.0, val2, ll_in, ll_out, prior;
  VECTOR *v_test, *v_in_vars, *v_in_means, *v_out_means, *v_out_vars ;
  VERTEX *v ;

  mri_out = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT) ;

  if (mris->hemisphere == LEFT_HEMISPHERE)
  {
    start_index = 0 ;
    end_index = mris->nvertices-1 ;
  }
  else   // processing rh
  {
    start_index = mris->nvertices ;
    end_index = 2*mris->nvertices-1 ;
  }


  v_in_vars = VectorAlloc(mri_cmat->width, MATRIX_REAL) ;
  v_in_means = VectorAlloc(mri_cmat->width, MATRIX_REAL) ;
  v_out_vars = VectorAlloc(mri_cmat->width, MATRIX_REAL) ;
  v_out_means = VectorAlloc(mri_cmat->width, MATRIX_REAL) ;
  m_train = MatrixAlloc(mri_cmat->width, mri_cmat->width, MATRIX_REAL) ;

  if (prior_only <= 0)
  {
    MRIthreshold(mri_cmat, mri_cmat, cor_thresh) ;
    MRIcmatNormalizeRows(mri_cmat) ;
    for (nin = nout = 0, sno = 0 ; sno < nsubjects-1 ; sno++)
    {
      MRISclearMarks(mris) ;
      LabelMark(labels[sno], mris) ;

      MatrixFromMRI(mri_cmat, m_train, sno) ;

      // compute variance within training MT and assume homoskedasticity
      for (ind = 0 ; ind < labels[sno]->n_points ; ind++)
      {
        vno = labels[sno]->lv[ind].vno ;
        for (vno2 = 0 ; vno2 < mri_cmat->height ; vno2++)
        {
          val = MRIgetVoxVal(mri_cmat, vno, vno2, 0, sno) ;
          VECTOR_ELT(v_in_vars, vno2+1) += val*val ;
          VECTOR_ELT(v_in_means, vno2+1) += val ;
        }
      }
      nin += labels[sno]->n_points ;
      // compute mean and variance vectors in region close to but outside MT
      MRISdistanceTransform(mris, labels[sno], DTRANS_MODE_SIGNED) ;
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->val > 1 && v->val < 10) // near MT, but not in it
        {
          nout++ ;
          for (vno2 = 0 ; vno2 < mri_cmat->height ; vno2++)
          {
            val = MRIgetVoxVal(mri_cmat, vno, vno2, 0, 0) ;
            VECTOR_ELT(v_out_vars, vno2+1) += val*val ;
            VECTOR_ELT(v_out_means, vno2+1) += val ;
          }
        }
      }
    }

    // normalize perimeter/outside distrubution
    for (vno2 = 0 ; vno2 < mri_cmat->height ; vno2++)
    {
      val = VECTOR_ELT(v_out_means, vno2+1) ;
      val /= nout ;  // mean
      val2 = VECTOR_ELT(v_out_vars, vno2+1) ;
      val2 = (val2 / nout - (val*val)) ; // variance
      VECTOR_ELT(v_out_means, vno2+1) = val ;
      VECTOR_ELT(v_out_vars, vno2+1) = val2 ;
    }
    // normalize MT distribution
    for (vno2 = 0 ; vno2 < mri_cmat->height ; vno2++)
    {
      val = VECTOR_ELT(v_in_means, vno2+1) ;
      val /= nin ;  // mean
      val2 = VECTOR_ELT(v_in_vars, vno2+1) ;
      val2 = (val2 / nin - (val*val)) ; // variance
      VECTOR_ELT(v_in_means, vno2+1) = val ;
      VECTOR_ELT(v_in_vars, vno2+1) = val2 ;
    }

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
  }
  else
  {
    nvertices = mris->nvertices ;
  }
  printf("classifying %d vertices in the FOV (%2.1f%%)\n",
         nvertices, 100.0*nvertices/(float)mris->nvertices) ;

  /* build the in_labels[subject][vno] array to list
     which vertices are in which subject's label */
  if (classifier == CLASSIFY_LABEL_FUSION)
  {
    for (sno = 0 ; sno < nsubjects-1 ; sno++)
    {
      in_label[sno] = (int *)calloc(mris->nvertices, sizeof(in_label[sno][0])) ;
      for (ind = 0 ; ind < labels[sno]->n_points ; ind++)
      {
        in_label[sno][labels[sno]->lv[ind].vno] = 1 ;
      }
      m_trains[sno] = MatrixFromMRI(mri_cmat, NULL, sno) ;
    }
  }

  for (nevals = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    if (!(vno % 100))
    {
      printf("%d of %d complete (%d considered)\r",
             vno, end_index-start_index, nevals) ;
      fflush(stdout) ;
    }
    v = &mris->vertices[vno] ;
    prior = MRIgetVoxVal(mri_prior, vno, 0, 0, 0) ;
    if (v->ripflag || prior < logodds_thresh)
    {
      continue ;
    }
    nevals++ ;

    // extract correlation pattern for test subject at this vertex
    v_test = VectorFromMRIcol(mri_cmat, NULL, vno+start_index, nsubjects-1) ;

    // use gaussian classifier for in and out of area
    if (classifier == CLASSIFY_GAUSSIAN)
    {
      ll_in = VectorLogLikelihood(v_test, v_in_means, v_in_vars) ;
      ll_out = VectorLogLikelihood(v_test, v_out_means, v_out_vars) ;
      v->val2 = ll_in ;
      v->val2bak = ll_out ;
      if (vno == Gdiag_no)
        printf("v %d: ll_in = %2.4f, ll_out = %2.4f, prior = %2.4f\n",
               Gdiag_no, ll_in, ll_out, prior) ;
      if (prior_only >= 0)
      {
        if (FZERO(prior))
        {
          ll_in += log(1e-10) ;
        }
        else
        {
          ll_in += log(prior) ;
        }
        if (FEQUAL(prior, 1))
        {
          ll_out += log(1e-10) ;
        }
        else
        {
          ll_out += log(1-prior) ;
        }
      }

      v->imag_val = ll_in ;
      v->stat = ll_out ;
      val = exp(ll_in) / (exp(ll_in) + exp(ll_out)) ;
      val = exp(ll_in) ;
      MRIsetVoxVal(mri_out, vno, 0, 0, 0, val) ;
    }
    else if (classifier == CLASSIFY_LABEL_FUSION)
    {
      int nsubj_at_this_vertex ;
      for (nsubj_at_this_vertex = sno = 0 ; sno < nsubjects-1 ; sno++)
      {
        if (in_label[sno][vno] > 0)
        {
          nsubj_at_this_vertex++ ;
        }
      }
      if (nsubj_at_this_vertex > 0)  // in at least one subject's label
      {
        double corrs[MAX_SUBJECTS], total_corr ;
        // extract correlation pattern for test subject at this vertex
        for (total_corr = 0.0, sno = 0 ; sno < nsubjects-1 ; sno++)
        {
          corrs[sno] = MatrixRowDotProduct(m_trains[sno], vno+1+start_index, v_test) ;
          total_corr += corrs[sno] ;
        }
        if (!DZERO(total_corr))
          for (val = 0.0, sno = 0 ; sno < nsubjects-1 ; sno++)
          {
            corrs[sno] /= total_corr ;
            if (in_label[sno][vno])
            {
              val += corrs[sno] ;
            }
          }
        else
          val = 0.0 ;
        MRIsetVoxVal(mri_out, vno, 0, 0, 0, val) ;
      }
    }
    else // use max similarity
    {
      // compute dot product with every vertex in training subject
      // and find max
      max_dot = -1e10 ;
      ind = -1 ;
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
      {
        MRIsetVoxVal(mri_out, vno, 0, 0, 0, 1) ;
      }
      // force mapping to be single-valued. This is order-dependent!
      MatrixSetRegion(m_train, m_train, ind+1, 1, 1, m_train->cols, 0) ;
    }
    VectorFree(&v_test) ;
  }
  if (classifier == CLASSIFY_LABEL_FUSION)
  {
    for (sno = 0 ; sno < nsubjects-1 ; sno++)
    {
      free(in_label[sno]) ;
      MatrixFree(&m_trains[sno]) ;
    }
  }
  printf("\n") ;
  MatrixFree(&m_train) ;
  VectorFree(&v_in_vars) ;
  VectorFree(&v_in_means) ;
  VectorFree(&v_out_vars) ;
  VectorFree(&v_out_means) ;

  return(mri_out) ;
}

static LABEL *
segment_area(MRI_SURFACE *mris, MRI *mri, LABEL *area, int nvertices)
{
  int          vno, max_vno, vno2, n, m ;
  double       max_p, p ;
  LABEL_VERTEX *lv ;

  area = LabelAlloc(nvertices, NULL, "segmented area") ;

  max_p = MRIgetVoxVal(mri, 0, 0, 0, 0) ;
  max_vno = 0 ;
  for (vno = 1 ; vno < mris->nvertices ; vno++)
  {
    p = MRIgetVoxVal(mri, vno, 0, 0, 0) ;
    if (p > max_p)
    {
      max_p = p ;
      max_vno = vno ;
    }
  }

  area->n_points = 1 ;
  lv = &area->lv[0] ;
  {
    VERTEX const * const v = &mris->vertices[max_vno] ;
    lv->x = v->x ;
    lv->y = v->y ;
    lv->z = v->z ;
    lv->stat = max_p ;
    lv->vno = max_vno ;
  }
  mris->vertices[max_vno].marked = 1 ;

  do
  {
    max_vno = -1 ;
    max_p = -100000 ;
    MRISclearMarks(mris) ;
    LabelMark(area, mris) ;
    for (n = 0 ; n < area->n_points ; n++)
    {
      vno = area->lv[n].vno ;
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      for (m = 0 ; m < vt->vnum ; m++)
      {
        vno2 = vt->v[m] ;
        if (mris->vertices[vno2].marked == 1)
        {
          continue ;  // already in the label
        }
        p = MRIgetVoxVal(mri, vno2, 0, 0, 0) ;
        if (p > max_p || max_vno < 0)
        {
          max_p = p ;
          max_vno = vno2 ;
        }
      }
    }
    lv = &area->lv[area->n_points++] ;
    { 
      VERTEX const * const v = &mris->vertices[max_vno] ;
      lv->vno = max_vno ;
      lv->x = v->x ;
      lv->y = v->y ;
      lv->z = v->z ;
      lv->stat = max_p ;
    }
    mris->vertices[max_vno].marked = 1 ;
  }
  while (area->n_points < nvertices) ;
  return(area) ;
}

int
main(int argc, char *argv[])
{
  char         **av, *out_fname, *subject, fname[STRLEN], *cp ;
  int          ac, nargs, sno ;
  int          msec, minutes, seconds, nsubjects ;
  Timer start ;
  MRI_SURFACE  *mris ;
  MRI          *mri_frame, *mri, *mri_out, *mri_prior ;
  LABEL        *labels[MAX_SUBJECTS], *area ;

  nargs = handleVersionOption(argc, argv, "mris_segment");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

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

  int req = snprintf(fname, STRLEN, "%s/fsaverage%d/surf/%s.inflated", sdir, ico_no, hemi_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mris = MRISread(fname) ;
  if (mris == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not load surface %s", Progname,fname) ;
  }
  mri = NULL ; // get rid of compiler warning
  for (sno = 0 ; sno < nsubjects ; sno++)
  {
    subject = argv[sno+1] ;
    printf("processing subject %s, %d of %d\n", subject, sno+1, nsubjects) ;
    int req = snprintf(fname, STRLEN, "%s/%s/%s/%s", sdir, subject, data_dir, data_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (prior_only <= 0)
    {
      mri_frame = MRIread(fname) ;
    }
    else
    {
      mri_frame = MRIalloc(2*mris->nvertices, 2*mris->nvertices, 1, MRI_FLOAT) ;
      MRIsetValues(mri_frame, .5) ;  // 0 and 1 are used as indicator values
    }
    if (mri_frame == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load correlation matrix from %s",
                Progname,fname) ;
    if (mri_frame->height == 1 && mri_frame->nframes == mri_frame->width)
    {
      MRI *mri_tmp  ;
      mri_tmp = mri_reshape(mri_frame,mri_frame->width,mri_frame->nframes,1,1);
      MRIfree(&mri_frame) ;
      mri_frame = mri_tmp ;
    }
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
    else if (mri_frame->width != mri->width ||
             mri_frame->height != mri->height ||
             mri_frame->depth != mri->depth)
      ErrorExit(ERROR_BADPARM,
                "%s: volume %s has incompatible dimensions (%d x %d x %d)",
                mri_frame->width, mri_frame->height, mri_frame->depth) ;

    MRIcopyFrame(mri_frame, mri, 0, sno) ;
    MRIfree(&mri_frame) ;

    req = snprintf(fname, STRLEN, "%s/%s/label/%s.%s", sdir, subject, hemi_name,label_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    labels[sno] = LabelRead(NULL, fname) ;
    if (labels[sno] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load label from %s",
                Progname,fname) ;
  }
  req = snprintf(fname, STRLEN, "%s/fsaverage%d/label/%s.%s",
		 sdir, ico_no, hemi_name, prior_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  mri_prior = MRIread(fname) ;
  if (mri_prior == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not load prior from %s",Progname,fname);
  }
  mri_out = classify_vertices(mris, mri_prior, mri, labels, nsubjects, cor_thresh, prior_only);
  printf("writing output to %s\n", out_fname) ;
  MRIwrite(mri_out, out_fname) ;
  MRISsmoothMRI(mris, mri_out, nsmooth, NULL, mri_out) ;
  area = segment_area(mris, mri_out, NULL, prior_number_of_vertices) ;
  LabelDilate(area, mris, nclose, CURRENT_VERTICES) ;
  LabelErode(area, mris, nclose) ;
  FileNameRemoveExtension(out_fname, out_fname) ;
  LabelWrite(area, out_fname) ;
  msec = start.milliseconds() ;
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
  if (!stricmp(option, "cmat") || !stricmp(option,"input") || !stricmp(option,"data"))
  {
    nargs = 1 ;
    strcpy(data_name, argv[2]) ;
    printf("using fmri/%s as name of correlation matrices\n", data_name) ;
  }
  else if (!stricmp(option, "smooth"))
  {
    nargs = 1 ;
    nsmooth = atoi(argv[2]) ;
    printf("smoothing posterior %d times\n", nsmooth) ;
  }
  else if (!stricmp(option, "prior"))
  {
    strcpy(prior_name, argv[2]) ;
    printf("using label/%s as name of prior label\n", prior_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "prioronly") || !stricmp(option, "prior_only"))
  {
    prior_only = 1 ;
    printf("only using priors in generating segmentation\n") ;
  }
  else if (!stricmp(option, "noprior") || !stricmp(option, "no_prior"))
  {
    prior_only = -1 ;
    printf("not using  priors in generating segmentation\n") ;
  }
  else if (!stricmp(option, "rh") || !stricmp(option, "lh"))
  {
    hemi_name = option ;
    printf("processing %s\n", hemi_name) ;
  }
  else if (!stricmp(option, "dir") || !stricmp(option, "input_dir"))
  {
    data_dir = argv[2] ;
    printf("reading data from subdirectory %s\n", data_dir) ;
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
    logodds_thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("using log odds thresh = %2.3f\n", logodds_thresh);
  }
  else switch (toupper(*option))
    {
    case 'G':
      classifier = CLASSIFY_GAUSSIAN ;
      printf("using Gaussian classifier\n") ;
      break ;
    case 'L':
      classifier = CLASSIFY_LABEL_FUSION ;
      printf("using label fusion classifier\n") ;
      break ;
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





