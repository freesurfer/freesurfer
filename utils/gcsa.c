#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mrisurf.h"
#include "mrishash.h"
#include "error.h"
#include "diag.h"
#include "utils.h"
#include "gcsa.h"
#include "proto.h"
#include "fio.h"
#include "transform.h"
#include "macros.h"
#include "utils.h"
#include "cma.h"
#include "icosahedron.h"


static int dump_gcsan(GCSA_NODE *gcsan, FILE *fp) ;
static int    GCSAupdateNodeMeans(GCSA_NODE *gcsa_node, 
                                  int label, double *v_inputs, int ninputs) ;
static int    GCSAupdateNodeCovariance(GCSA_NODE *gcsan, int label, 
                                       double *v_inputs, int ninputs) ;
static VERTEX *GCSAsourceToAverageVertex(GCSA *gcsa, VERTEX *v) ;
static int    GCSANclassify(GCSA_NODE *gcsa_node, double *v_inputs, 
                            int ninputs, double *pprob) ;

GCSA  *
GCSAalloc(int ninputs, int icno)
{
  char      fname[STRLEN], *cp ;
  int       nvertices, i, n ;
  GCSA      *gcsa ;
  GCSA_NODE *gcsan ;
  GCS       *gcs ;
  double    max_len ;

  gcsa = (GCSA *)calloc(1, sizeof(GCSA)) ;
  if (!gcsa)
    ErrorExit(ERROR_NOMEMORY, "GCSAalloc(%d, %d): could not allocate gcsa",
              ninputs, icno) ;

  gcsa->icno = icno ;
  cp = getenv("MRI_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM,"%s: MRI_DIR not defined in environment",Progname);

  sprintf(fname, "%s/lib/bem/ic%d.tri", cp, icno) ;
  gcsa->mris = ICOread(fname) ; 
  if (!gcsa->mris)
    ErrorExit(ERROR_NOMEMORY, "GCSAalloc(%d, %d): could not read ico %s",
              ninputs, icno, fname) ;
  nvertices = gcsa->mris->nvertices ;
  MRISprojectOntoSphere(gcsa->mris, gcsa->mris, DEFAULT_RADIUS) ;
  MRIScomputeVertexSpacingStats(gcsa->mris, NULL, NULL, &max_len, NULL,NULL);
  gcsa->mht = MHTfillVertexTableRes(gcsa->mris, NULL, CURRENT_VERTICES,
                                    2*max_len);
  
  gcsa->ninputs = ninputs ;
  gcsa->nodes = (GCSA_NODE *)calloc(nvertices, sizeof(GCSA_NODE)) ;
  if (!gcsa->nodes)
    ErrorExit(ERROR_NOMEMORY, "GCSAalloc(%d, %d): could not allocate nodes",
              ninputs, nvertices) ;

#define DEFAULT_GCS  1

  for (i = 0 ; i < nvertices ; i++)
  {
    if (i == Gdiag_no)
      DiagBreak() ;
    gcsan = &gcsa->nodes[i] ;
    gcsan->max_labels = DEFAULT_GCS ;
    gcsan->gcs = (GCS *)calloc(DEFAULT_GCS, sizeof(GCS)) ;
    if (!gcsan->gcs)
      ErrorExit(ERROR_NOMEMORY, "GCSAalloc(%d, %d): could not allocate gcs %d",
                ninputs, nvertices, i) ;
    gcsan->labels = (int *)calloc(DEFAULT_GCS, sizeof(int)) ;
    if (!gcsan->labels)
      ErrorExit(ERROR_NOMEMORY, 
                "GCSAalloc(%d, %d): could not allocate labels %d",
                ninputs, nvertices, i) ;
    for (n = 0 ; n < DEFAULT_GCS ; n++)
    {
      gcs = &gcsan->gcs[n] ;
      gcs->v_means = VectorAlloc(ninputs, MATRIX_REAL) ;
      gcs->m_cov = MatrixAlloc(ninputs, ninputs, MATRIX_REAL) ;
    }
  }

  return(gcsa) ;
}

int
GCSAfree(GCSA **pgcsa)
{
  int       nvertices, i, n ;
  GCSA      *gcsa ;
  GCSA_NODE *gcsan ;
  GCS       *gcs ;

  gcsa = *pgcsa ; *pgcsa = NULL ;

  nvertices = gcsa->mris->nvertices ;
  MRISfree(&gcsa->mris) ;
  MHTfree(&gcsa->mht) ;


  for (i = 0 ; i < nvertices ; i++)
  {
    if (i == Gdiag_no)
      DiagBreak() ;
    gcsan = &gcsa->nodes[i] ;
    
    for (n = 0 ; n < gcsan->nlabels ; n++)
    {
      gcs = &gcsan->gcs[n] ;
      VectorFree(&gcs->v_means) ; MatrixFree(&gcs->m_cov) ;
    }
    free(gcsan->gcs) ;
  }

  free(gcsa->nodes) ; free(gcsa) ;
  return(NO_ERROR) ;
}

/*
  v->x,y,z     should have canonical coordinates in them and
  v->origx,y,z should have original coordinates.
  second fundamental form should have been computed on
*/
int
GCSAtrainMeans(GCSA *gcsa, MRI_SURFACE *mris)
{
  int        vno, vno_avg ;
  VERTEX     *v, *v_avg ;
  GCSA_NODE  *gcsa_node ;
  double     v_inputs[100] ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    
    v_avg = GCSAsourceToAverageVertex(gcsa, v) ;
    vno_avg = v_avg - gcsa->mris->vertices ;
    if (vno_avg == Gdiag_no)
      DiagBreak() ;
    gcsa_node = &gcsa->nodes[vno_avg] ;
    v_inputs[0] = v->val ;
    if (gcsa->ninputs > 1)
      v_inputs[1] = v->val2 ;
    if (gcsa->ninputs > 2)
      v_inputs[2] = v->imag_val ;
    if (vno == Gdiag_no && v->annotation == 1336341)
      DiagBreak() ;
    if (v->annotation == 0)  /* not labeled */
      continue ;
    GCSAupdateNodeMeans(gcsa_node, v->annotation, v_inputs, gcsa->ninputs) ;
  }
  return(NO_ERROR) ;
}

int
GCSAtrainCovariances(GCSA *gcsa, MRI_SURFACE *mris)
{
  int        vno, vno_avg ;
  VERTEX     *v, *v_avg ;
  GCSA_NODE  *gcsa_node ;
  double     v_inputs[100] ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    
    v_avg = GCSAsourceToAverageVertex(gcsa, v) ;
    vno_avg = v_avg - gcsa->mris->vertices ;
    if (vno_avg == Gdiag_no)
      DiagBreak() ;
    gcsa_node = &gcsa->nodes[vno_avg] ;
    v_inputs[0] = v->val ;
    if (gcsa->ninputs > 1)
      v_inputs[1] = v->val2 ;
    if (gcsa->ninputs > 2)
      v_inputs[2] = v->imag_val ;
    if (vno == Gdiag_no && v->annotation == 1336341)
      DiagBreak() ;
    if (v->annotation == 0)  /* not labeled */
      continue ;
    GCSAupdateNodeCovariance(gcsa_node, v->annotation, v_inputs,gcsa->ninputs);
  }
  return(NO_ERROR) ;
}
static VERTEX *
GCSAsourceToAverageVertex(GCSA *gcsa, VERTEX *v)
{
  VERTEX *vdst ;

  vdst = MHTfindClosestVertex(gcsa->mht, gcsa->mris, v) ;
  return(vdst) ;
}


static int
GCSAupdateNodeMeans(GCSA_NODE *gcsan, int label, double *v_inputs, int ninputs)
{
  int   n, i ;
  GCS   *gcs ;

  if (label == 0)
    DiagBreak() ;

  for (n = 0 ; n < gcsan->nlabels ; n++)
  {
    if (gcsan->labels[n] == label)
      break ;
  }
  if (n >= gcsan->nlabels)  /* have to allocate a new classifier */
  {
    if (n >= gcsan->max_labels)
    {
      int   old_max_labels, i ;
      int   *old_labels ;
      GCS   *old_gcs ;

      if (gcsan->max_labels >= 6)
        DiagBreak() ;
      old_max_labels = gcsan->max_labels ; 
      gcsan->max_labels += 2 ;
      old_labels = gcsan->labels ;
      old_gcs = gcsan->gcs ;

      /* allocate new ones */
      gcsan->gcs = (GCS *)calloc(gcsan->max_labels, sizeof(GCS)) ;
      if (!gcsan->gcs)
        ErrorExit(ERROR_NOMEMORY, 
                  "GCANupdateNodeMeans: couldn't expand gcs to %d",
                  gcsan->max_labels) ;
      gcsan->labels = (int *)calloc(gcsan->max_labels, sizeof(int)) ;
      if (!gcsan->labels)
        ErrorExit(ERROR_NOMEMORY, 
                  "GCANupdateNodeMeans: couldn't expand labels to %d",
                  gcsan->max_labels) ;
      /* copy the old ones over */
      for (i = 0 ; i < old_max_labels ; i++)
      {
        gcsan->gcs[i].v_means = old_gcs[i].v_means ;
        gcsan->gcs[i].m_cov = old_gcs[i].m_cov ;
        gcsan->gcs[i].prior = old_gcs[i].prior ;

        /* copy over Gibbs stuff XXX */
      }

      /* allocate new ones */
      for (i = old_max_labels ; i < gcsan->max_labels ; i++)
      {
        gcs = &gcsan->gcs[i] ;
        gcs->v_means = VectorAlloc(ninputs, MATRIX_REAL) ;
        gcs->m_cov = MatrixAlloc(ninputs, ninputs, MATRIX_REAL) ;
      }

      memmove(gcsan->labels, old_labels, old_max_labels*sizeof(int)) ;

      /* free the old ones */
      free(old_gcs) ; free(old_labels) ;
    }
    gcsan->nlabels++ ;
  }

  gcs = &gcsan->gcs[n] ;

  gcsan->labels[n] = label ;
  gcs->prior += 1.0f ;
  gcsan->total_training++ ;

  /* these will be updated when training is complete */
  for (i = 0 ; i < ninputs ; i++)
    VECTOR_ELT(gcs->v_means, i+1) += v_inputs[i] ; 
  

  return(NO_ERROR) ;
}

int
GCSAnormalizeMeans(GCSA *gcsa)
{
  int        vno, n, total_gcs ;
  GCSA_NODE  *gcsan ;
  GCS        *gcs ;

  for (total_gcs = vno = 0 ; vno < gcsa->mris->nvertices ; vno++)
  {
    gcsan = &gcsa->nodes[vno] ;
    total_gcs += gcsan->nlabels ;
    for (n = 0 ; n < gcsan->nlabels ; n++)
    {
      gcs = &gcsan->gcs[n] ;
      MatrixScalarMul(gcs->v_means, 1.0/gcs->prior, gcs->v_means) ;

      /* the priors will get updated later */
    }
  }

  printf("%d classifier means computed, %2.1f/vertex\n",
         total_gcs, (float)total_gcs/gcsa->mris->nvertices) ;
  return(NO_ERROR) ;
}

static int
GCSAupdateNodeCovariance(GCSA_NODE *gcsan, int label, double *v_inputs, 
                         int ninputs)
{
  int   n, i, j ;
  GCS   *gcs ;
  double cov, mean1, mean2 ;

  for (n = 0 ; n < gcsan->nlabels ; n++)
  {
    if (gcsan->labels[n] == label)
      break ;
  }
  if (n >= gcsan->nlabels) 
  {
    dump_gcsan(gcsan, stderr) ;
    ErrorExit(ERROR_BADPARM, "GCSAupdateNodeCovariance(%d): unknown label\n",
              label) ;
  }
  gcs = &gcsan->gcs[n] ;

  /* these will be updated when training is complete */
  for (i = 0 ; i < ninputs ; i++)
  {
    mean1 = VECTOR_ELT(gcs->v_means, i+1) ;
    for (j = 0 ; j <= i ; j++)
    {
      mean2 = VECTOR_ELT(gcs->v_means, j+1) ;
      cov = (v_inputs[i]-mean1) * (v_inputs[j]-mean2) ;
      *MATRIX_RELT(gcs->m_cov, i+1, j+1) += cov ;
      *MATRIX_RELT(gcs->m_cov, j+1, i+1) += cov ;
    }
  }
  

  return(NO_ERROR) ;
}

int
GCSAnormalizeCovariances(GCSA *gcsa)
{
  int        vno, n, cno ;
  GCSA_NODE  *gcsan ;
  GCS        *gcs ;
  MATRIX     *m_inv ;

  for (vno = 0 ; vno < gcsa->mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    gcsan = &gcsa->nodes[vno] ;
    for (n = 0 ; n < gcsan->nlabels ; n++)
    {
      gcs = &gcsan->gcs[n] ;
      if (gcs->prior > 1)
        MatrixScalarMul(gcs->m_cov, 1.0/(gcs->prior-1), gcs->m_cov) ;
      else  /* only 1 observation - can't estimate covariances */
        MatrixIdentity(gcsa->ninputs, gcs->m_cov) ;

      /* check to see if covariance is singular, and if so regularize it */
      m_inv = MatrixInverse(gcs->m_cov, NULL) ;
      cno = MatrixConditionNumber(gcs->m_cov) ;
      if (m_inv == NULL || (cno > 100) || (cno <= 0))
      {
        m_inv = MatrixIdentity(gcsa->ninputs, NULL) ;
        MatrixScalarMul(m_inv, 0.1, m_inv) ;
        MatrixAdd(m_inv, gcs->m_cov, gcs->m_cov) ;
      }
      MatrixFree(&m_inv) ;
      gcs->prior /= (float)gcsan->total_training ;
    }
  }

  return(NO_ERROR) ;
}
int
GCSAwrite(GCSA *gcsa, char *fname)
{
  FILE       *fp ;
  int        vno, n ;
  GCSA_NODE  *gcsan ;
  GCS        *gcs ;

  fp = fopen(fname, "wb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "GCSAwrite(%s): could not open file", fname)) ;

  fwriteInt(GCSA_MAGIC, fp) ;
  fwriteInt(gcsa->ninputs, fp) ;
  fwriteInt(gcsa->icno, fp) ;
  for (vno = 0 ; vno < gcsa->mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    gcsan = &gcsa->nodes[vno] ;
    fwriteInt(gcsan->nlabels, fp) ;
    fwriteInt(gcsan->total_training, fp) ;
    
    for (n = 0 ; n < gcsan->nlabels ; n++)
    {
      gcs = &gcsan->gcs[n] ;
      fwriteInt(gcsan->labels[n], fp) ;
      fwriteFloat(gcs->prior, fp) ;
      MatrixAsciiWriteInto(fp, gcs->v_means) ;
      MatrixAsciiWriteInto(fp, gcs->m_cov) ;
    }
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}

GCSA *
GCSAread(char *fname)
{
  FILE       *fp ;
  int        vno, n, ninputs, icno, magic ;
  GCSA_NODE  *gcsan ;
  GCS        *gcs ;
  GCSA       *gcsa ;

  
  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL, 
                (ERROR_NOFILE, "GCSAread(%s): could not open file", fname)) ;

  magic = freadInt(fp) ;
  if (magic != GCSA_MAGIC)
    ErrorReturn(NULL,
                (ERROR_BADFILE, 
                 "GCSAread(%s): file magic #%x != GCSA_MAGIC (%x)", 
                 fname, magic, GCSA_MAGIC)) ;
  ninputs = freadInt(fp) ;
  icno = freadInt(fp) ;

  gcsa = GCSAalloc(ninputs,icno) ;

  for (vno = 0 ; vno < gcsa->mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    gcsan = &gcsa->nodes[vno] ;
    for (n = 0 ; n < gcsan->nlabels ; n++)
    {
      MatrixFree(&gcsan->gcs[n].v_means) ;
      MatrixFree(&gcsan->gcs[n].m_cov) ;
    }
    gcsan->nlabels = freadInt(fp) ;
    free(gcsan->labels) ; free(gcsan->gcs) ;
    gcsan->gcs = (GCS *)calloc(gcsan->nlabels, sizeof(GCS)) ;
    if (!gcsan->gcs)
      ErrorExit(ERROR_NOMEMORY, "GCSAread(%s): could not allocate gcs %d",
                fname, vno) ;
    gcsan->labels = (int *)calloc(gcsan->nlabels, sizeof(int)) ;
    if (!gcsan->labels)
      ErrorExit(ERROR_NOMEMORY, 
                "GCSAread(%s): could not allocate labels %d", fname, vno) ;
    
    gcsan->total_training = freadInt(fp) ;
    
    for (n = 0 ; n < gcsan->nlabels ; n++)
    {
      gcs = &gcsan->gcs[n] ;
      gcsan->labels[n] = freadInt(fp) ;
      gcs->prior = freadFloat(fp) ;
      gcs->v_means = MatrixAsciiReadFrom(fp, NULL) ;
      gcs->m_cov = MatrixAsciiReadFrom(fp, NULL) ;
    }
  }

  fclose(fp) ;
  return(gcsa) ;
}

int
GCSAlabel(GCSA *gcsa, MRI_SURFACE *mris)
{
  int        vno, vno_avg, label ;
  VERTEX     *v, *v_avg ;
  GCSA_NODE  *gcsa_node ;
  double     v_inputs[100], p ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    
    v_avg = GCSAsourceToAverageVertex(gcsa, v) ;
    vno_avg = v_avg - gcsa->mris->vertices ;
    if (vno_avg == Gdiag_no)
      DiagBreak() ;
    gcsa_node = &gcsa->nodes[vno_avg] ;
    v_inputs[0] = v->val ;
    if (gcsa->ninputs > 1)
      v_inputs[1] = v->val2 ;
    if (gcsa->ninputs > 2)
      v_inputs[2] = v->imag_val ;
    label = GCSANclassify(gcsa_node, v_inputs, gcsa->ninputs, &p) ;
    v->annotation = label ;
  }
  return(NO_ERROR) ;
}
static int
GCSANclassify(GCSA_NODE *gcsan, double *v_inputs, int ninputs, 
              double *pprob)
{
  int    n, best_label ;
  double p, ptotal, max_p ;
  GCS    *gcs ;
  static MATRIX *m_cov_inv = NULL ;
  static VECTOR *v_means_T = NULL, *v_tmp = NULL ;

  if (v_means_T && ninputs != v_means_T->cols)
  {
    VectorFree(&v_means_T) ; MatrixFree(&m_cov_inv) ; MatrixFree(&v_tmp) ;
  }

  ptotal = 0.0 ; max_p = -10000 ; best_label = 0 ;
  for (n = 0 ; n < gcsan->nlabels ; n++)
  {
    gcs = &gcsan->gcs[n] ;
    v_means_T = MatrixTranspose(gcs->v_means, v_means_T) ;
    m_cov_inv = MatrixInverse(gcs->m_cov, m_cov_inv) ;
    if (!m_cov_inv)
    {
      MATRIX *m_tmp ;
      fprintf(stderr, "Singular matrix in GCSAclassify:\n") ;
      MatrixPrint(stderr, gcs->m_cov) ;
#if 0
      ErrorExit(ERROR_BADPARM, "") ;
#else
      m_tmp = MatrixIdentity(ninputs, NULL) ;
      MatrixScalarMul(m_tmp, 0.1, m_tmp) ;
      MatrixAdd(m_tmp, gcs->m_cov, m_tmp) ;
      m_cov_inv = MatrixInverse(m_tmp, NULL) ;
      if (!m_cov_inv)
      {
        ErrorExit(ERROR_BADPARM, "GCSANclassify: could not regularize matrix") ;
      }
#endif
    }
    v_tmp = MatrixMultiply(m_cov_inv, gcs->v_means, v_tmp) ;
    p = VectorDot(v_means_T, v_tmp) ;
    p = gcs->prior * exp(-0.5 * p) ; ptotal += p ;
    if (p > max_p)
    {
      max_p = p ; best_label = gcsan->labels[n] ;
    }
  }
  if (pprob)
    *pprob = max_p / ptotal ;
  return(best_label) ;
}

static int
dump_gcsan(GCSA_NODE *gcsan, FILE *fp)
{
  int   n ;
  GCS   *gcs ;

  fprintf(fp, "GCSAN with %d labels (%d training examples)\n", 
          gcsan->nlabels, gcsan->total_training) ;
  for (n = 0 ; n < gcsan->nlabels ; n++)
  {
    fprintf(fp, "\t%d: label %d\n", n, gcsan->labels[n]) ;
    gcs = &gcsan->gcs[n] ;
    fprintf(fp, "\tprior %2.1f\n", gcs->prior) ;
    fprintf(fp, "mean:\n") ; MatrixPrint(fp, gcs->v_means) ;
    fprintf(fp, "covariance:\n") ; MatrixPrint(fp, gcs->m_cov) ;
  }
  return(NO_ERROR) ;
}
