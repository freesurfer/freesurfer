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

#define BIG_AND_NEGATIVE            -10000000.0

static int edge_to_index(VERTEX *v, VERTEX *vn) ;
static int dump_gcsan(GCSA_NODE *gcsan, FILE *fp) ;
static int    MRIScomputeVertexPermutation(MRI_SURFACE *mris, int *indices) ;
static int    GCSAupdateNodeMeans(GCSA_NODE *gcsan, 
                                  int label, double *v_inputs, int ninputs) ;
static int    GCSAupdateNodeGibbsPriors(GCSA_NODE *gcsan, int label, 
                                        MRI_SURFACE *mris, int vno) ;
static int    GCSAupdateNodeCovariance(GCSA_NODE *gcsan, int label, 
                                       double *v_inputs, int ninputs) ;
static VERTEX *GCSAsourceToAverageVertex(GCSA *gcsa, VERTEX *v) ;
static int    GCSANclassify(GCSA_NODE *gcsa_node, double *v_inputs, 
                            int ninputs, double *pprob) ;
static double gcsaNbhdGibbsLogLikelihood(GCSA *gcsa, MRI_SURFACE *mris, 
                                         double *v_inputs, int vno, 
                                         double gibbs_coef,
                                         int label) ;
static double gcsaVertexGibbsLogLikelihood(GCSA *gcsa, MRI_SURFACE *mris, 
                                           double *v_inputs, int vno, 
                                           double gibbs_coef) ;

GCSA  *
GCSAalloc(int ninputs, int icno)
{
  char      fname[STRLEN], *cp ;
  int       nvertices, i, n ;
  GCSA      *gcsa ;
  GCSA_NODE *gcsan ;
  GCS       *gcs ;
  double    max_len ;

  read_annotation_table() ;  /* for debugging */

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
    GCSAupdateNodeGibbsPriors(gcsa_node, v->annotation, mris, vno) ;
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
  int        vno, n, cno, i, j ;
  GCSA_NODE  *gcsan ;
  GCS        *gcs ;
  MATRIX     *m_inv ;
  double     mean, min_var ;

#define MIN_VAR 0.01

  for (vno = 0 ; vno < gcsa->mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    gcsan = &gcsa->nodes[vno] ;
    for (n = 0 ; n < gcsan->nlabels ; n++)
    {
      gcs = &gcsan->gcs[n] ;
      if (gcs->prior >= 2)
        MatrixScalarMul(gcs->m_cov, 1.0/(gcs->prior-1), gcs->m_cov) ;

      if (gcs->prior < 5)    /* not that many samples - regularize */
      {
        for (i = 1 ; i <= gcsa->ninputs ; i++)
        {
          mean = VECTOR_ELT(gcs->v_means, i) ;
          min_var = (.1*mean * .1*mean) ;
          if (gcs->prior > 1)
            min_var *= (gcs->prior-1) ;
          if (min_var < MIN_VAR)
            min_var = MIN_VAR ;
          if (*MATRIX_RELT(gcs->m_cov, i, i) < min_var)
            *MATRIX_RELT(gcs->m_cov, i, i) = min_var ;
        }
      }
      else
      {
        for (i = 1 ; i <= gcsa->ninputs ; i++)
        {
          mean = VECTOR_ELT(gcs->v_means, i) ;
          min_var = mean*.1 ;
          if (min_var < MIN_VAR)
            min_var = MIN_VAR ;
          if (FZERO(*MATRIX_RELT(gcs->m_cov, i, i)))
            *MATRIX_RELT(gcs->m_cov, i, i) = min_var ;
        }
      }

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
      m_inv = MatrixIdentity(gcsa->ninputs, NULL) ;
      if (!m_inv)
      {
        fprintf(stderr, "vno %d, n = %d, m_cov is singular!\n",
                vno, n) ;
      }
      else
        MatrixFree(&m_inv) ;
      for (i = 0 ; i < GIBBS_SURFACE_NEIGHBORS ; i++)
      {
        for (j = 0 ; j < gcs->nlabels[i] ; j++)
          gcs->label_priors[i][j] /= (float)gcs->total_nbrs[i] ;
      }
      gcs->prior /= (float)gcsan->total_training ;
    }
  }

  return(NO_ERROR) ;
}
int
GCSAwrite(GCSA *gcsa, char *fname)
{
  FILE       *fp ;
  int        vno, n, i, j ;
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
      
      for (i = 0 ; i < GIBBS_SURFACE_NEIGHBORS ; i++)
      {
        /*        long where ;*/

        /*        where = ftell(fp) ;*/
        fwriteInt(gcs->total_nbrs[i], fp) ;
        fwriteInt(gcs->nlabels[i], fp) ;
        /*        where = ftell(fp) ;*/
        for (j = 0 ; j < gcs->nlabels[i] ; j++)
        {
          fwriteInt(gcs->labels[i][j], fp) ;
          fwriteFloat(gcs->label_priors[i][j], fp) ;
        }
      }
    }
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}

GCSA *
GCSAread(char *fname)
{
  FILE       *fp ;
  int        vno, n, ninputs, icno, magic, i, j ;
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
      for (i = 0 ; i < GIBBS_SURFACE_NEIGHBORS ; i++)
      {
        /*        long  where ;*/
        /*        where = ftell(fp) ;*/
        gcs->total_nbrs[i] = freadInt(fp) ;
        gcs->nlabels[i] = freadInt(fp) ;
        /*        where = ftell(fp) ;*/
        gcs->labels[i] = (int *)calloc(gcs->nlabels[i], sizeof(int)) ;
        if (!gcs->labels[i])
          ErrorExit(ERROR_NOMEMORY, 
                    "GCSAread: couldn't allocate labels to %d",
                    gcs->nlabels[i]+1) ;

        gcs->label_priors[i] = (float *)calloc(gcs->nlabels[i], sizeof(float));
        if (!gcs->label_priors[i])
          ErrorExit(ERROR_NOMEMORY, "GCSAread: couldn't allocate gcs to %d" ,
                    gcs->nlabels[i]+1) ;
        for (j = 0 ; j < gcs->nlabels[i] ; j++)
        {
          gcs->labels[i][j] = freadInt(fp) ;
          gcs->label_priors[i][j] = freadFloat(fp) ;
        }
      }
    }
  }

  fclose(fp) ;
  return(gcsa) ;
}

static int Gvno = -1 ;
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
    Gvno = vno_avg ;
    label = GCSANclassify(gcsa_node, v_inputs, gcsa->ninputs, &p) ;
    v->annotation = label ;
  }
  return(NO_ERROR) ;
}
static int
GCSANclassify(GCSA_NODE *gcsan, double *v_inputs, int ninputs, 
              double *pprob)
{
  int    n, best_label, i ;
  double p, ptotal, max_p ;
  GCS    *gcs ;
  static MATRIX *m_cov_inv = NULL ;
  static VECTOR *v_x_T = NULL, *v_tmp = NULL, *v_x = NULL ;

  if (v_x_T && ninputs != v_x_T->cols)
  {
    VectorFree(&v_x_T) ; MatrixFree(&m_cov_inv) ; MatrixFree(&v_tmp) ;
    VectorFree(&v_x) ;
  }

  ptotal = 0.0 ; max_p = -10000 ; best_label = 0 ;
  for (n = 0 ; n < gcsan->nlabels ; n++)
  {
    gcs = &gcsan->gcs[n] ;
    v_x = VectorCopy(gcs->v_means, v_x) ;
    for (i = 0 ; i < ninputs ; i++)
      VECTOR_ELT(v_x, i+1) -= v_inputs[i] ;
    v_x_T = MatrixTranspose(v_x, v_x_T) ;
    m_cov_inv = MatrixInverse(gcs->m_cov, m_cov_inv) ;
    if (!m_cov_inv)
    {
      MATRIX *m_tmp ;
      fprintf(stderr, "Singular matrix in GCSAclassify,vno=%d,n=%d:\n",Gvno,n);
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
        ErrorExit(ERROR_BADPARM, "GCSANclassify: could not regularize matrix");
      }
#endif
    }
    v_tmp = MatrixMultiply(m_cov_inv, v_x, v_tmp) ;
    p = VectorDot(v_x_T, v_tmp) ;
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
  int   n, index ;
  GCS   *gcs ;
  char  *name ;

  fprintf(fp, "GCSAN with %d labels (%d training examples)\n", 
          gcsan->nlabels, gcsan->total_training) ;
  for (n = 0 ; n < gcsan->nlabels ; n++)
  {
    name = annotation_to_name(gcsan->labels[n],&index) ;
    fprintf(fp, "  %d: label %s (%d, %d)\n", n, name, index,gcsan->labels[n]) ;
    gcs = &gcsan->gcs[n] ;
    fprintf(fp, "\tprior %2.1f\n", gcs->prior) ;
    if (gcs->v_means->rows > 1)
    {
      fprintf(fp, "\tmean:\n") ;       MatrixPrint(fp, gcs->v_means) ;
      fprintf(fp, "\tcovariance:\n") ; MatrixPrint(fp, gcs->m_cov) ;
    }
    else
    {
      fprintf(fp, "\tmean %2.2f +- %2.1f\n",
              VECTOR_ELT(gcs->v_means,1), sqrt(*MATRIX_RELT(gcs->m_cov,1,1))) ;
    }
  }
  return(NO_ERROR) ;
}
int
GCSAdump(GCSA *gcsa, int vno, MRI_SURFACE *mris, FILE *fp)
{
  int        vno_avg ;
  VERTEX     *vavg, *v ;
  GCSA_NODE  *gcsan ;

  v = &mris->vertices[vno] ; v->tx = v->x ; v->ty = v->y ; v->tz = v->z ;
  v->x = v->cx ; v->y = v->cy ; v->z = v->cz ;
  vavg = GCSAsourceToAverageVertex(gcsa, v) ;
  vno_avg = vavg - gcsa->mris->vertices ;
  gcsan = &gcsa->nodes[vno_avg] ;
  fprintf(fp, "v %d --> vavg %d\n", vno, vno_avg) ;
  dump_gcsan(gcsan, fp) ;
  v->x = v->tx ; v->y = v->ty ; v->z = v->tz ;
  return(NO_ERROR) ;
}

typedef struct
{
  int    index ;
  int    r, g, b ;
  int    annotation ;
  char   name[100] ;
} ATABLE_ELT ;

static ATABLE_ELT *atable ;
static int num_entries = 0 ;

int
read_annotation_table(void)
{
  FILE  *fp ;
  char  *cp, fname[STRLEN], line[STRLEN] ;
  int   i ;

  if (num_entries)
    return(NO_ERROR) ;   /* already read */

  cp = getenv("MRI_DIR") ;
  if (!cp)
    cp = "." ;

  sprintf(fname, "%s/christophe_parc.txt", cp) ;
  fp = fopen(fname, "r") ;
  if (!fp)
  {
    fprintf(stderr, "could not open translation file %s\n", fname) ;
    return(ERROR_NO_FILE) ;
  }

  num_entries = 0 ;
  do
  {
    cp = fgetl(line, 199, fp) ;
    if (!cp)
      break ;
    num_entries++ ;
  } while (cp && !feof(fp)) ;

  rewind(fp) ;

  atable = (ATABLE_ELT *)calloc(num_entries, sizeof(ATABLE_ELT)) ;
  for (i = 0 ; i < num_entries ; i++)
  {
    cp = fgetl(line, 199, fp) ;
    if (!cp)
      break ;
    sscanf(cp, "%d %s %d %d %d %*d",
           &atable[i].index,
           atable[i].name,
           &atable[i].r,
           &atable[i].g,
           &atable[i].b) ;
    atable[i].annotation = atable[i].r+(atable[i].g << 8)+(atable[i].b << 16);
  }
  return(NO_ERROR) ;
}
char *
annotation_to_name(int annotation, int *pindex)
{
  int   i ;

  if (num_entries <= 0)
  {
    static char name[100] ;

    if (pindex)
      *pindex = -1 ;
    sprintf(name, "%d", annotation) ;
    return(name) ;
  }

  for (i = 0 ; i < num_entries ; i++)
  {
    if (atable[i].annotation == annotation)
    {
      if (pindex)
        *pindex = atable[i].index ;
      return(atable[i].name) ;
    }
  }
  if (pindex)
    *pindex = -1 ;
  return("NOT_FOUND") ;
}

static int
GCSAupdateNodeGibbsPriors(GCSA_NODE *gcsan,int label, MRI_SURFACE *mris, 
                          int vno)
{
  int     n, i, m, nbr_label ;
  GCS    *gcs ;
  VERTEX *v, *vn ;

  if (label == 0)
    DiagBreak() ;

  for (n = 0 ; n < gcsan->nlabels ; n++)
  {
    if (gcsan->labels[n] == label)
      break ;
  }
  if (n >= gcsan->nlabels)  /* have to allocate a new classifier */
    ErrorExit(ERROR_BADPARM, 
              "GCSAupdateNodeGibbsPriors(%d): could not find label %d",
              vno, label) ;

  gcs = &gcsan->gcs[n] ;
  v = &mris->vertices[vno] ;
  for (m = 0 ; m < v->vnum ; m++)
  {
    vn = &mris->vertices[v->v[m]] ;
    i = edge_to_index(v, vn) ;
    gcs->total_nbrs[i]++ ;
    nbr_label = vn->annotation ;
    
    /* see if this label already exists */
    for (n = 0 ; n < gcs->nlabels[i] ; n++)
      if (gcs->labels[i][n] == nbr_label)
        break ;
    if (n >= gcs->nlabels[i])   /* not there - reallocate stuff */
    {
      int   *old_labels ;
      float *old_label_priors ;

      old_labels = gcs->labels[i] ;
      old_label_priors = gcs->label_priors[i] ;

      /* allocate new ones */
      gcs->label_priors[i] = (float *)calloc(gcs->nlabels[i]+1, sizeof(float));
      if (!gcs->label_priors[i])
        ErrorExit(ERROR_NOMEMORY, "GCAupdateNodeGibbsPriors: "
                  "couldn't expand gcs to %d" ,gcs->nlabels[i]+1) ;
      gcs->labels[i] = (int *)calloc(gcs->nlabels[i]+1, sizeof(int)) ;
      if (!gcs->labels[i])
        ErrorExit(ERROR_NOMEMORY, 
                  "GCSAupdateNodeGibbsPriors: couldn't expand labels to %d",
                  gcs->nlabels[i]+1) ;

      if (gcs->nlabels[i] > 0)  /* copy the old ones over */
      {
        memmove(gcs->label_priors[i], old_label_priors, 
                gcs->nlabels[i]*sizeof(float)) ;
        memmove(gcs->labels[i], old_labels, gcs->nlabels[i]*sizeof(int)) ;

        /* free the old ones */
        free(old_label_priors) ; free(old_labels) ;
      }
      gcs->labels[i][gcs->nlabels[i]++] = nbr_label ;
    }
    gcs->label_priors[i][n] += 1.0f ;
  }

  return(NO_ERROR) ;
}
static int 
edge_to_index(VERTEX *v, VERTEX *vn)
{
  float  dx, dy, dz, dot1, dot2, index, x_dist, y_dist, z_dist ;

  /* first find out which of the principal directions the edge connecting v 
     and vn is most closely parallel to */
  dx = vn->origx - v->origx ; dy = vn->origy - v->origy ; 
  dz = vn->origz - v->origz ;

  dot1 = dx*v->e1x + dy*v->e1y + dz*v->e1z ;
  dot2 = dx*v->e2x + dy*v->e2y + dz*v->e2z ;

  x_dist = fabs(1 - dx) ; y_dist = fabs(1 - dy) ; z_dist = fabs(1 - dz) ;
  if (fabs(dot1) > fabs(dot2))   /* 1st principal direction - 0 or 1 */
    index = 0 ;
  else   /* second principal direction - 1 or 2 */
    index = 2 ;

  return(index) ;
  if (x_dist > y_dist && x_dist > z_dist)   /* use sign of x */
    return(dx > 0 ? index : index+1) ;
  else if (y_dist > z_dist)
    return(dy > 0 ? index : index+1) ;
  else
    return(dz > 0 ? index : index+1) ;
}

#define MIN_CHANGED 0

int gcsa_write_iterations = 0 ;
char *gcsa_write_fname = NULL ;

int
GCSAreclassifyUsingGibbsPriors(GCSA *gcsa, MRI_SURFACE *mris)
{
  int       *indices ;
  int        n, vno, i, nchanged, label, best_label, old_label, vno_avg, niter;
  double     ll, max_ll ;
  VERTEX     *v, *v_avg, *vn ;
  GCSA_NODE  *gcsan ;
  double     v_inputs[100] ;

  indices = (int *)calloc(mris->nvertices, sizeof(int)) ;

  niter = 0 ;
  if (gcsa_write_iterations != 0)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s%03d.annot", gcsa_write_fname, niter) ;
    printf("writing snapshot to %s...\n", fname) ;
    MRISwriteAnnotation(mris, fname) ;
  }

  /* mark all vertices, so they will all be considered the first time through*/
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].marked = 1 ;
  do
  {
    nchanged = 0 ;
    MRIScomputeVertexPermutation(mris, indices) ;
    for (i = 0 ; i < mris->nvertices ; i++)
    {
      vno = indices[i] ;
      v = &mris->vertices[vno] ;

      if (vno == Gdiag_no)
        DiagBreak() ;

      v_avg = GCSAsourceToAverageVertex(gcsa, v) ;
      vno_avg = v_avg - gcsa->mris->vertices ;
      if (vno_avg == Gdiag_no)
        DiagBreak() ;

      v_inputs[0] = v->val ;
      if (gcsa->ninputs > 1)
        v_inputs[1] = v->val2 ;
      if (gcsa->ninputs > 2)
        v_inputs[2] = v->imag_val ;
      
      gcsan = &gcsa->nodes[vno_avg] ;
      if (gcsan->nlabels <= 1)
        continue ;

      best_label = old_label = v->annotation ; 
      max_ll =gcsaNbhdGibbsLogLikelihood(gcsa,mris,v_inputs,vno,1.0,old_label);
      for (n = 0 ; n < gcsan->nlabels ; n++)
      {
        label = gcsan->labels[n] ;
        ll = gcsaNbhdGibbsLogLikelihood(gcsa, mris, v_inputs, vno, 1.0,label) ;
        if (ll > max_ll)
        {
          max_ll = ll ;
          best_label = label ;
        }
      }
      if (best_label != old_label)
      {
        v->marked = 1 ;
        nchanged++ ; v->annotation = best_label ;
      }
      else
        v->marked = 0 ;
    }
    printf("%03d: %d changed...\n", niter, nchanged) ;
    niter++ ;
    if (gcsa_write_iterations && (niter % gcsa_write_iterations) == 0)
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s%03d.annot", gcsa_write_fname, niter) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRISwriteAnnotation(mris, fname) ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ; 
      if (v->marked != 1)
        continue ;
      for (n = 0 ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked == 1)
          continue ;
        vn->marked = 2 ;
      }
    }
  } while (nchanged > MIN_CHANGED) ;
  free(indices) ;
  return(NO_ERROR) ;
}
int
MRIScomputeVertexPermutation(MRI_SURFACE *mris, int *indices)
{
  int i, index, tmp ;


  for (i = 0 ; i < mris->nvertices ; i++)
  {
    indices[i] = i ;
  }
  for (i = 0 ; i < mris->nvertices ; i++)
  {
    index = (int)randomNumber(0.0, (double)(mris->nvertices-0.0001)) ;
    tmp = indices[index] ; 
    indices[index] = indices[i] ; indices[i] = tmp ;
  }
  return(NO_ERROR) ;
}

static double
gcsaNbhdGibbsLogLikelihood(GCSA *gcsa, MRI_SURFACE *mris, double *v_inputs, 
                             int vno, double gibbs_coef, int label)
{
  double   total_ll, ll ;
  int      n, old_annotation ;
  VERTEX   *v ;

  v = &mris->vertices[vno] ;
  old_annotation = v->annotation ;
  v->annotation = label ;

  total_ll = gcsaVertexGibbsLogLikelihood(gcsa, mris, v_inputs, vno,
                                          gibbs_coef) ;

  for (n = 0 ; n < v->vnum ; n++)
  {
    ll = gcsaVertexGibbsLogLikelihood(gcsa, mris, v_inputs, v->v[n],
                                      gibbs_coef) ;
    total_ll += ll ;
  }

  v->annotation = old_annotation ;
  return(total_ll) ;
}
static double
gcsaVertexGibbsLogLikelihood(GCSA *gcsa, MRI_SURFACE *mris, double *v_inputs,
                            int vno, double gibbs_coef)
{
  double    ll, nbr_prior ;
  int       nbr_label, label, i,j, n, vno_avg ;
  GCSA_NODE *gcsan ;
  GCS       *gcs ;
  VERTEX    *v, *v_avg ;
  static MATRIX *m_cov_inv = NULL ;
  static VECTOR *v_x_T = NULL, *v_tmp = NULL, *v_x = NULL ;

  if (v_x_T && gcsa->ninputs != v_x_T->cols)
  {
    VectorFree(&v_x_T) ; MatrixFree(&m_cov_inv) ; MatrixFree(&v_tmp) ;
    VectorFree(&v_x) ;
  }


  v = &mris->vertices[vno] ;
  v_avg = GCSAsourceToAverageVertex(gcsa, v) ;
  vno_avg = v_avg - gcsa->mris->vertices ;
  if (vno_avg == Gdiag_no)
    DiagBreak() ;
  gcsan = &gcsa->nodes[vno_avg] ;
  label = v->annotation ;

  for (n = 0 ; n < gcsan->nlabels ; n++)
  {
    if (gcsan->labels[n] == label)
      break ;
  }
  if (n >= gcsan->nlabels)  /* never occured here */
    return(BIG_AND_NEGATIVE) ;

  gcs = &gcsan->gcs[n] ;
  
  /* compute Mahalanobis distance */
  v_x = VectorCopy(gcs->v_means, v_x) ;
  for (i = 0 ; i < gcsa->ninputs ; i++)
    VECTOR_ELT(v_x, i+1) -= v_inputs[i] ;
  v_x_T = MatrixTranspose(v_x, v_x_T) ;
  m_cov_inv = MatrixInverse(gcs->m_cov, m_cov_inv) ;
  if (!m_cov_inv)
    ErrorExit(ERROR_BADPARM,
              "GCSAvertexLogLikelihood: could not invert matrix");

  v_tmp = MatrixMultiply(m_cov_inv, v_x, v_tmp) ;
  ll = -0.5 *VectorDot(v_x_T, v_tmp) ;

  nbr_prior = 0.0 ;
  for (n = 0 ; n < v->vnum ; n++)
  {
    nbr_label = mris->vertices[v->v[n]].annotation ;
    i = edge_to_index(v, &mris->vertices[v->v[n]]) ;
    for (j = 0 ; j < gcs->nlabels[i] ; j++)
    {
      if (nbr_label == gcs->labels[i][j])
        break ;
    }
    if (j < gcs->nlabels[i])
    {
      if (!FZERO(gcs->label_priors[i][j]))
        nbr_prior += log(gcs->label_priors[i][j]) ;
      else
        nbr_prior += BIG_AND_NEGATIVE ;
    }
    else   /* never occurred - make it unlikely */
    {
      nbr_prior += BIG_AND_NEGATIVE ;
    }
  }
  ll += (gibbs_coef * nbr_prior + log(gcs->prior)) ;

  return(ll) ;
}
