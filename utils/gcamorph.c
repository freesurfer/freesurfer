#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "error.h"
#include "fio.h"
#include "diag.h"
#include "gca.h"
#include "gcamorph.h"
#include "transform.h"
#include "proto.h"
#include "mrimorph.h"
#include "mrinorm.h"
#include "matrix.h"
#include "cma.h"
#include "mri.h"
#define GCAM_VERSION   1.0
#define MIN_STD (2.0)
#define MIN_VAR (MIN_STD*MIN_STD)

#define GCAM_X_GRAD 0
#define GCAM_Y_GRAD 1
#define GCAM_Z_GRAD 2

#ifndef FSIGN
#define FSIGN(f)  (((f) < 0) ? -1 : 1)
#endif

static MRI *gcamWriteMRI(GCA_MORPH *gcam, MRI *mri, int which) ;
static int gcamClearMomentum(GCA_MORPH *gcam) ;
static int gcamRemoveNegativeNodes(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms) ;
static int gcamMLElabelAtLocation(GCA_MORPH *gcam, int x, int y, int z, float *vals) ;
static double gcamFindOptimalTimeStep(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms, MRI *mri) ;
static int  gcamAreaTermAtNode(GCA_MORPH *gcam, MRI *mri, double l_area, 
                                   int i, int j, int k, double *pdx, double *pdy, 
                                   double *pdz) ;
static int  gcamJacobianTermAtNode(GCA_MORPH *gcam, MRI *mri, double l_jacobian, 
                                   int i, int j, int k, double *pdx, double *pdy, 
                                   double *pdz) ;
static int   finitep(float f) ;
static double gcamComputeSSE(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms) ;
static double gcamComputeRMS(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms) ;
static int write_snapshot(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms, int iter) ;
static int  log_integration_parms(FILE *fp, GCA_MORPH_PARMS *parms) ;
static int gcamLimitGradientMagnitude(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms, MRI *mri) ;
static int gcamComputeGradient(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, 
                               GCA_MORPH_PARMS *parms) ;
static int gcamLogLikelihoodTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth,
                              double l_log_likelihood) ;
static int gcamLikelihoodTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth,
                              double l_likelihood, GCA_MORPH_PARMS *parms) ;
static double gcamMapEnergy(GCA_MORPH *gcam, MRI *mri) ;
static double gcamLabelEnergy(GCA_MORPH *gcam, MRI *mri, double label_dist) ;
static int gcamMapTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, double l_map) ;
static int gcamLabelTerm(GCA_MORPH *gcam, MRI *mri, double l_label, double label_dist) ;
static int gcamDistanceTerm(GCA_MORPH *gcam, MRI *mri, double l_distance) ;
static int gcamClearGradient(GCA_MORPH *gcam) ;
static int gcamComputeMetricProperties(GCA_MORPH *gcam) ;
static double gcamLogLikelihoodEnergy(GCA_MORPH *gcam, MRI *mri) ;

static int gcamDistanceTerm(GCA_MORPH *gcam, MRI *mri, double l_distance) ;
static double gcamDistanceEnergy(GCA_MORPH *gcam, MRI *mri) ;

static int gcamSmoothnessTerm(GCA_MORPH *gcam, MRI *mri, double l_smoothness) ;
static double gcamSmoothnessEnergy(GCA_MORPH *gcam, MRI *mri) ;

static int gcamAreaTerm(GCA_MORPH *gcam, MRI *mri, double l_jacobian) ;
static double gcamAreaEnergy(GCA_MORPH *gcam, MRI *mri) ;

static int gcamJacobianTerm(GCA_MORPH *gcam, MRI *mri, double l_jacobian, double ratio_thresh) ;
static double gcamJacobianEnergy(GCA_MORPH *gcam, MRI *mri) ;
static int gcamApplyGradient(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms) ;
static int gcamSmoothGradient(GCA_MORPH *gcam, int navgs) ;
static int gcamUndoGradient(GCA_MORPH *gcam) ;
static int check_gcam(GCAM *gcam) ;
#define DEFAULT_PYRAMID_LEVELS 3
#define MAX_PYRAMID_LEVELS     6
#define MAX_EXP          200
#define GCAMN_SUB(mns1, mns2, v)     \
    V3_LOAD(v, mns1->x - mns2->x, mns1->y - mns2->y, mns1->z - mns2->z)

int
GCAMwrite(GCA_MORPH *gcam, char *fname)
{
  FILE            *fp ;
  int             x, y, z ;
  GCA_MORPH_NODE  *gcamn ;

  fp = fopen(fname, "wb") ;
  if (!fp)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GCAMwrite(%s): could not open file",
                                fname)) ;

  fwriteFloat(GCAM_VERSION, fp) ;
  fwriteInt(gcam->width, fp) ;
  fwriteInt(gcam->height, fp) ;
  fwriteInt(gcam->depth, fp) ;
  fwriteInt(gcam->spacing, fp) ;
  fwriteFloat(gcam->exp_k, fp) ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        fwriteFloat(gcamn->origx, fp) ;
        fwriteFloat(gcamn->origy, fp) ;
        fwriteFloat(gcamn->origz, fp) ;

        fwriteFloat(gcamn->x, fp) ;
        fwriteFloat(gcamn->y, fp) ;
        fwriteFloat(gcamn->z, fp) ;

        fwriteInt(gcamn->xn, fp) ;
        fwriteInt(gcamn->yn, fp) ;
        fwriteInt(gcamn->zn, fp) ;
      }
    }
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}

int
GCAMregister(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  char   fname[STRLEN] ;
  int    level, i, level_steps, navgs, l2, relabel, orig_relabel ;
  MRI    *mri_smooth = NULL, *mri_kernel ;
  double base_sigma, pct_change, rms, last_rms, orig_dt,l_smooth ;
  
  navgs = parms->navgs ; orig_dt = parms->dt ; l_smooth = parms->l_smoothness;
  if (FZERO(parms->exp_k))
    parms->exp_k = EXP_K ;
  if (parms->levels < 0)
    parms->levels = DEFAULT_PYRAMID_LEVELS ;
  else if (parms->levels >= MAX_PYRAMID_LEVELS)
    parms->levels = MAX_PYRAMID_LEVELS ;

  gcam->exp_k = parms->exp_k ;
  parms->mri = mri ;
  if (Gdiag & DIAG_WRITE)
  {
    sprintf(fname, "%s.log", parms->base_name) ;
    if (parms->start_t == 0)
      parms->log_fp = fopen(fname, "w") ;
    else
      parms->log_fp = fopen(fname, "a") ;
  }
  else
    parms->log_fp = NULL ;

  base_sigma = parms->sigma ;
  if (parms->relabel_avgs >= parms->navgs)
    GCAMcomputeLabels(mri, gcam) ;
  else
    GCAMcomputeMaxPriorLabels(gcam) ;
  orig_relabel = relabel = parms->relabel ;
  for (level = parms->levels-1 ; level >= 0 ; level--)
  {
    if (parms->reset_avgs == parms->navgs)
    {
      printf("resetting metric properties...\n") ;
      GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, ORIGINAL_POSITIONS) ;
      gcamComputeMetricProperties(gcam) ;
      GCAMstoreMetricProperties(gcam) ;
    }
    parms->relabel = (parms->relabel_avgs >= parms->navgs) ;
    parms->sigma = base_sigma ;
#if 0
    parms->l_smoothness = l_smooth / (sqrt(parms->navgs)+1) ;
#else
    parms->l_smoothness = l_smooth / ((parms->navgs)+1) ;
#endif
    printf("setting smoothness coefficient to %2.2f\n", parms->l_smoothness);
    for (l2 = 0 ; l2 < parms->levels ; l2++)
    {
      if (parms->sigma < 0.4)
	break ;
      if (FZERO(parms->sigma))
	mri_smooth = MRIcopy(mri, mri_smooth) ;
      else
      {
	printf("blurring input image with Gaussian with sigma=%2.3f...\n", parms->sigma) ;
	mri_kernel = MRIgaussian1d(parms->sigma, 100) ;
      	mri_smooth = MRIconvolveGaussian(mri, NULL, mri_kernel) ;
	MRIfree(&mri_kernel) ;
      }
      i = 0 ;
      do
      {
	if (((level != (parms->levels-1)) || (i > 0)) && parms->relabel)
	{
      	  GCAMcomputeLabels(mri, gcam) ;
	  if (parms->write_iterations != 0)
	  {
	    char fname[STRLEN] ;
	    MRI  *mri_gca ;
	    mri_gca = MRIclone(mri, NULL) ;
	    GCAMbuildMostLikelyVolume(gcam, mri_gca) ;
	    sprintf(fname, "%s_target%d", parms->base_name, level) ;
	    MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
	    sprintf(fname, "%s_target%d.mgh", parms->base_name, level) ;
	    printf("writing target volume to %s...\n", fname) ;
	    MRIwrite(mri_gca, fname) ;
	    MRIfree(&mri_gca) ;
	  }
	}
	last_rms = gcamComputeRMS(gcam, mri, parms) ;
	level_steps = parms->start_t ;
	GCAMregisterLevel(gcam, mri, mri_smooth, parms) ;
      	rms = gcamComputeRMS(gcam, mri, parms) ;
	level_steps = parms->start_t - level_steps ;   /* # of steps taken in GCAMregisterLevel */
	if (level_steps == 0)
	  level_steps = 1 ;
      	pct_change = 100.0*(last_rms-rms)/(level_steps*last_rms) ;
#if 0
	printf("iter %d: last rms %2.3f, rms %2.3f, pct_change %2.3f/iter\n",
	       i+1, last_rms, rms, pct_change) ;
#endif
	if (i++ >= 0)  /* don't bother iterating */
      	  break ;
      } while (pct_change > parms->tol) ;
      parms->sigma /= 4 ; 
    }
    parms->navgs /= 4 ;
  }
  MRIfree(&mri_smooth) ;
  
  parms->sigma = base_sigma ;
  if (parms->log_fp)
  {
    fclose(parms->log_fp) ;
    parms->log_fp = NULL ;
  }

  parms->relabel = orig_relabel ;
  parms->navgs = navgs ; parms->dt = orig_dt ;
  return(NO_ERROR) ;
}

GCA_MORPH *
GCAMread(char *fname)
{
  GCA_MORPH       *gcam ;
  FILE            *fp ;
  int             x, y, z, width, height, depth ;
  GCA_MORPH_NODE  *gcamn ;
  float           version ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_BADPARM, "GCAMread(%s): could not open file",
                                fname)) ;

  version = freadFloat(fp) ;
  if (version != GCAM_VERSION)
  {
    fclose(fp) ;
    ErrorReturn(NULL, 
                (ERROR_BADFILE, "GCAMread(%s): invalid version # %2.3f\n", version)) ;
  }
  width = freadInt(fp) ; height = freadInt(fp) ; depth = freadInt(fp) ;
  gcam = GCAMalloc(width, height, depth) ;
  
  gcam->spacing = freadInt(fp) ;
  gcam->exp_k = freadFloat(fp) ;

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        gcamn->origx = freadFloat(fp) ;
        gcamn->origy = freadFloat(fp) ;
        gcamn->origz = freadFloat(fp) ;

        gcamn->x = freadFloat(fp) ;
        gcamn->y = freadFloat(fp) ;
        gcamn->z = freadFloat(fp) ;

        gcamn->xn = freadInt(fp) ;
        gcamn->yn = freadInt(fp) ;
        gcamn->zn = freadInt(fp) ;
      }
    }
  }

  fclose(fp) ;

  return(gcam) ;
}


GCA_MORPH *
GCAMalloc(int width, int height, int depth)
{
  GCA_MORPH  *gcam ;
  int        x, y ;

  gcam = (GCA_MORPH *)calloc(1, sizeof(GCA_MORPH)) ;
  if (!gcam)
    ErrorExit(ERROR_NOMEMORY, "GCAMalloc: could not allocate GCA_MORPH struct") ;

  gcam->width =  width  ; gcam->height = height ; gcam->depth =  depth  ;

  gcam->nodes = (GCA_MORPH_NODE ***)calloc(width, sizeof(GCA_MORPH_NODE **)) ;
  if (!gcam->nodes)
    ErrorExit(ERROR_NOMEMORY, "GCAMalloc: could not allocate nodes") ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    gcam->nodes[x] = (GCA_MORPH_NODE **)calloc(gcam->height, sizeof(GCA_MORPH_NODE *)) ;
    if (!gcam->nodes[x])
      ErrorExit(ERROR_NOMEMORY, "GCAMalloc: could not allocate %dth **",x) ;

    for (y = 0 ; y < gcam->height ; y++)
    {
      gcam->nodes[x][y] = (GCA_MORPH_NODE *)calloc(gcam->depth, sizeof(GCA_MORPH_NODE)) ;
      if (!gcam->nodes[x][y])
        ErrorExit(ERROR_NOMEMORY,"GCAMalloc: could not allocate %d,%dth *",x,y);
    }
  }
  return(gcam) ;
}


int
GCAMinit(GCA_MORPH *gcam, MRI *mri, GCA *gca, TRANSFORM *transform, int relabel)
{
  GCA_MORPH_NODE  *gcamn ;
  GC1D            *gc ;
  GCA_PRIOR       *gcap ;
  int             x, y, z, width, height, depth, n, max_n, max_label, label ;
  float           max_p, ox, oy, oz ;

  width = gcam->width ; height = gcam->height ; depth = gcam->depth ;
  TransformInvert(transform, mri) ;
  gcam->gca = gca ; gcam->spacing = gca->prior_spacing ; 
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        gcap = &gca->priors[x][y][z] ;
        max_p = 0 ;  max_n = -1 ; max_label = 0 ;

        for (n = 0 ; n < gcap->nlabels ; n++)
        {
          label = gcap->labels[n] ;
          if (label == Gdiag_no)
            DiagBreak() ;
          if (gcap->priors[n] >= max_p)
          {
            max_n = n ;
            max_p = gcap->priors[n] ;
            max_label = gcap->labels[n] ;
          }
        }
        gcamn->xn = x ; gcamn->yn = y ; gcamn->zn = z ;
        if (!GCApriorToSourceVoxelFloat(gca, mri, transform, x, y, z,
					&ox, &oy, &oz))
	{
	  gcamn->x = gcamn->origx = ox ; 
	  gcamn->y = gcamn->origy = oy ; 
	  gcamn->z = gcamn->origz = oz ;
#if 1
	  gcamn->label = max_label ;
	  gcamn->n = max_n ;
	  gcamn->prior = max_p ;
	  gc = GCAfindPriorGC(gca, x, y, z, max_label) ;
	  gcamn->gc = gc ;
	  gcamn->log_p = 0 ;
#endif
	  if (x == Gx && y == Gy && z == Gz)
	  {
	    printf("node(%d,%d,%d) --> MRI (%2.1f, %2.1f, %2.1f)\n",
		   x, y, z, ox, oy, oz) ;
	    DiagBreak() ;
#if 0
	    Gvx = nint(ox) ; Gvy = nint(oy) ; Gvz = nint(oz) ;  /* for writing out image views */
#endif
	  }
	}// !GCA
      }
    }
  }

  if (relabel)
    GCAMcomputeLabels(mri, gcam) ;
  else
    GCAMcomputeMaxPriorLabels(gcam) ;
  
  gcamComputeMetricProperties(gcam) ;
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
      for (z = 0 ; z < depth ; z++)
        gcam->nodes[x][y][z].orig_area = gcam->nodes[x][y][z].area ;

  gcamComputeMetricProperties(gcam) ;
  return(NO_ERROR) ;
}

int
GCAMfree(GCA_MORPH **pgcam)
{
  GCA_MORPH *gcam ;
  int       x, y ;

  gcam = *pgcam ; *pgcam = NULL ;

  GCAMfreeInverse(gcam) ;
  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
      free(gcam->nodes[x][y]) ;
    free(gcam->nodes[x]) ;
  }
  free(gcam->nodes) ;
  free(gcam) ;
  return(NO_ERROR) ;
}


/*
	Note:  d [ I(r)' C I(r), r] = delI * C * I(r)
  where delI(r) = 3xn, C=nxn, and I(r)=nx1
*/
static int
gcamLogLikelihoodTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, double l_log_likelihood)
{
  int             x, y, z, n /*,label*/ ;
  Real            dx, dy, dz, norm;
	float           vals[MAX_GCA_INPUTS] ;
  GCA_MORPH_NODE  *gcamn ;
  MATRIX          *m_delI, *m_inv_cov ;
  VECTOR          *v_means, *v_grad ;

  if (FZERO(l_log_likelihood))
    return(NO_ERROR) ;

  m_delI = MatrixAlloc(3, gcam->gca->ninputs, MATRIX_REAL) ;
  m_inv_cov = MatrixAlloc(gcam->gca->ninputs, gcam->gca->ninputs, MATRIX_REAL) ;
  v_means = VectorAlloc(gcam->gca->ninputs, 1) ;
  v_grad = VectorAlloc(3, MATRIX_REAL) ;
  
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
	if (fabs(gcamn->x-Gvx)<1 && fabs(gcamn->y-Gvy)<1  && fabs(gcamn->z-Gvz)<1)
	  DiagBreak() ;
	if (gcamn->status & GCAM_IGNORE_LIKELIHOOD)
	  continue ;
	
	load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->gca->ninputs) ;
	if (!gcamn->gc)
	{
	  MatrixClear(v_means) ; MatrixIdentity(gcam->gca->ninputs, m_inv_cov) ;
	  MatrixScalarMul(m_inv_cov, 1.0/(MIN_VAR), m_inv_cov) ; /* variance=4 is min */
	}
	else
	{
#if 0
	  if (parms->relabel)
	  {
	    label = 
	      GCAcomputeMAPlabelAtLocation(gcam->gca, x,y,z,vals,&n,&gcamn->log_p);
	    if (label == gcamn->label)  /* already correct label - don't move anywhere */
	      continue ;
	  }
#endif
	  load_mean_vector(gcamn->gc, v_means, gcam->gca->ninputs) ;
	  load_inverse_covariance_matrix(gcamn->gc, m_inv_cov, gcam->gca->ninputs) ;
	}
	
	for (n = 0 ; n < gcam->gca->ninputs ; n++)
	{
	  MRIsampleVolumeGradientFrame(mri_smooth, gcamn->x, gcamn->y, gcamn->z, &dx, &dy, &dz, n) ;
	  norm = sqrt(dx*dx+dy*dy+dz*dz) ;
	  if (!FZERO(norm))  /* don't worry about magnitude of gradient */
	  { dx /= norm ; dy /= norm ; dz /= norm ; }
	  *MATRIX_RELT(m_delI, 1, n+1) = dx ;
	  *MATRIX_RELT(m_delI, 2, n+1) = dy ;
	  *MATRIX_RELT(m_delI, 3, n+1) = dz ;
	  VECTOR_ELT(v_means, n+1) -= vals[n] ;
#define MAX_ERROR 100
	  if (fabs(VECTOR_ELT(v_means, n+1)) > MAX_ERROR)
	    VECTOR_ELT(v_means, n+1) = MAX_ERROR * FSIGN(VECTOR_ELT(v_means, n+1)) ;
	}
	
	MatrixMultiply(m_inv_cov, v_means, v_means) ;
	MatrixMultiply(m_delI, v_means, v_grad) ;
	
        gcamn->dx += l_log_likelihood*V3_X(v_grad) ;
        gcamn->dy += l_log_likelihood*V3_Y(v_grad) ;
        gcamn->dz += l_log_likelihood*V3_Z(v_grad) ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_like: node(%d,%d,%d): dI=(%2.1f,%2.1f,%2.1f), grad=(%2.2f,%2.2f,%2.2f), "
                 "node %2.2f+-%2.2f, MRI=%2.1f\n",
                 x, y, z, dx, dy, dz, gcamn->dx, gcamn->dy, gcamn->dz,
                 gcamn->gc ? gcamn->gc->means[0] : 0.0, 
		 gcamn->gc ? sqrt(covariance_determinant(gcamn->gc, gcam->gca->ninputs)) : 0.0, vals[0]) ;
      }
  
  MatrixFree(&m_delI) ; MatrixFree(&m_inv_cov) ; VectorFree(&v_means) ; VectorFree(&v_grad) ;
  return(NO_ERROR) ;
}

/*
*/
static int
gcamLikelihoodTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, double l_likelihood, GCA_MORPH_PARMS *parms)
{
  int             x, y, z, len, xi, yi, zi, x0, y0, z0, half_len ;
  Real            dx, dy, dz, xp1, yp1, zp1, xm1, ym1, zm1 ;
	float           val, mean, var ;
  GCA_MORPH_NODE  *gcamn ;
	MRI             *mri_nbhd, *mri_kernel ;

  if (FZERO(l_likelihood))
    return(NO_ERROR) ;

  len = (int)nint(4.0f * parms->sigma)+1 ;
  if (ISEVEN(len))   /* ensure it's odd */
    len++ ;
	half_len = (len-1)/2 ;

	mri = MRISeqchangeType(mri, MRI_FLOAT, 0, 0, 1) ;
	mri_nbhd = MRIallocSequence(len,len,len, MRI_FLOAT, gcam->gca->ninputs) ;
	mri_kernel = MRIgaussian1d(parms->sigma, 0) ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
				if (fabs(gcamn->x-Gvx)<1 && fabs(gcamn->y-Gvy)<1  && fabs(gcamn->z-Gvz)<1)
					DiagBreak() ;
				if (gcamn->status & GCAM_IGNORE_LIKELIHOOD || (gcamn->gc == NULL))
					continue ;

				x0 = nint(gcamn->x) ; y0 = nint(gcamn->y) ; z0 = nint(gcamn->z) ; 
				if ((x0+half_len >= mri->width) ||
						(y0+half_len >= mri->height) ||
						(z0+half_len >= mri->depth))
					continue ;
				MRIextractInto(mri, mri_nbhd, x0-half_len, y0-half_len, z0-half_len,len, len, len, 0, 0, 0) ;
				if (gcam->gca->ninputs == 1)
				{
					var = gcamn->gc->covars[0] ; mean = gcamn->gc->means[0] ;

					for (xi = 0 ;  xi < mri_nbhd->width ; xi++)
					{
						for (yi = 0 ;  yi < mri_nbhd->height ; yi++)
						{
							for (zi = 0 ;  zi < mri_nbhd->depth ; zi++)
							{
								val = MRIFvox(mri_nbhd, xi, yi, zi) ;
								val = (val - mean) * (val-mean) / var ;
								MRIFvox(mri_nbhd, xi, yi, zi) = sqrt(val) ;
							}
						}
					}
				}
				else  /* haven't written the multi-input case yet */
				{
					ErrorExit(ERROR_UNSUPPORTED, "vector-based likelihood not written yet") ;
				}
				MRIconvolveGaussian(mri_nbhd, mri_nbhd, mri_kernel) ;
				MRIsampleVolumeType(mri_nbhd, half_len+1, half_len, half_len, &xp1, SAMPLE_NEAREST) ;
				MRIsampleVolumeType(mri_nbhd, half_len-1, half_len, half_len, &xm1, SAMPLE_NEAREST) ;
				MRIsampleVolumeType(mri_nbhd, half_len, half_len+1, half_len, &yp1, SAMPLE_NEAREST) ;
				MRIsampleVolumeType(mri_nbhd, half_len, half_len-1, half_len, &ym1, SAMPLE_NEAREST) ;
				MRIsampleVolumeType(mri_nbhd, half_len, half_len, half_len+1, &zp1, SAMPLE_NEAREST) ;
				MRIsampleVolumeType(mri_nbhd, half_len, half_len, half_len-1, &zm1, SAMPLE_NEAREST) ;

				dx = -(xp1 - xm1)/2 ; dy = -(yp1 - ym1)/2 ; dz = -(zp1 - zm1)/2 ;
        gcamn->dx += l_likelihood*dx ;
        gcamn->dy += l_likelihood*dy ;
        gcamn->dz += l_likelihood*dz ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_like: node(%d,%d,%d): dp=(%2.2f,%2.2f,%2.2f), node %2.2f+-%2.2f\n",
                 x, y, z, gcamn->dx, gcamn->dy, gcamn->dz,
                 gcamn->gc ? gcamn->gc->means[0] : 0.0, 
								 gcamn->gc ? sqrt(covariance_determinant(gcamn->gc, gcam->gca->ninputs)) : 0.0) ;
      }

	MRIfree(&mri_kernel) ; MRIfree(&mri_nbhd) ; MRIfree(&mri) ;
  return(NO_ERROR) ;
}

static float last_sse[128][128][128];

static double
gcamLogLikelihoodEnergy(GCA_MORPH *gcam, MRI *mri)
{
  double          sse = 0.0, error ;
  int             x, y, z, max_x, max_y, max_z ;
  Real            max_increase = 0, increase ;
  GCA_MORPH_NODE  *gcamn ;
  float            vals[MAX_GCA_INPUTS] ;

  max_x = max_y = max_z = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
	if (gcamn->status & GCAM_IGNORE_LIKELIHOOD)
	  continue ;
	load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->gca->ninputs) ;
	if (gcamn->gc)
	  error = GCAmahDist(gcamn->gc, vals, gcam->gca->ninputs) + log(covariance_determinant(gcamn->gc, gcam->gca->ninputs)) ;
	else
	{
	  int n ;
	  for (n = 0, error = 0.0 ; n < gcam->gca->ninputs ; n++)
	    error += (vals[n]*vals[n]/MIN_VAR) ;
	}
	
        if (x == Gx && y == Gy && z == Gz)
          printf("E_like: node(%d,%d,%d) -> (%2.1f,%2.1f,%2.1f), target=%2.1f+-%2.1f, val=%2.1f\n",
                 x, y, z, gcamn->x, gcamn->y, gcamn->z,
                 gcamn->gc ? gcamn->gc->means[0] : 0.0, 
		 gcamn->gc ? sqrt(covariance_determinant(gcamn->gc, gcam->gca->ninputs)) : 0.0, vals[0]) ;
	
        if (last_sse[x][y][z] < (.9*error) && !FZERO(last_sse[x][y][z]))
        {
          DiagBreak() ;
          increase = error - last_sse[x][y][z] ;
          if (increase > max_increase)
          {
            max_increase = increase ; max_x = x ; max_y = y ; max_z = z ;
          }
        }
        last_sse[x][y][z] = (error) ;
        sse += (error) ;
        if (!finite(sse))
          DiagBreak() ;
      }
	if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
		printf("max increase %2.2f at (%d, %d, %d)\n",
		       max_increase, max_x, max_y, max_z) ;
  return(sse) ;
}


static int 
gcamDistanceTerm(GCA_MORPH *gcam, MRI *mri, double l_distance)
{
  double          norm, dx, dy, dz, error, d0, d, xdelta, ydelta, zdelta ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  if (FZERO(l_distance))
    return(NO_ERROR) ;
  width = gcam->width ; height = gcam->height ; depth = gcam->depth ; 
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        num = 0 ;
        xdelta = ydelta = zdelta = 0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ; xn = MAX(0,xn) ; xn = MIN(width-1,xn) ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ; yn = MAX(0,yn) ; yn = MIN(height-1,yn) ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              if (!xk && !yk && !zk)
                continue ;
              zn = z+zk ; zn = MAX(0,zn) ; zn = MIN(depth-1,zn) ;
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;
#if 0
              if (gcamn_nbr->label != gcamn->label)
                continue ;
#endif
              dx = gcamn->origx - gcamn_nbr->origx ;
              dy = gcamn->origy - gcamn_nbr->origy ;
              dz = gcamn->origz - gcamn_nbr->origz ;
              d0 = sqrt(dx*dx + dy*dy + dz*dz) ;
              dx = gcamn->x - gcamn_nbr->x ;
              dy = gcamn->y - gcamn_nbr->y ;
              dz = gcamn->z - gcamn_nbr->z ;
              d = sqrt(dx*dx + dy*dy + dz*dz) ;
              norm = d0 ; 
              if (FZERO(norm))
                norm = 1 ;
              dx /= norm ; dy /= norm ; dz /= norm ; 
              error = d0-d ;
              xdelta += error*dx ; 
              ydelta += error*dy; 
              zdelta += error*dz ; 
              num++ ;
            }
          }
        }
        num = 1 ;
        if (num > 0)
        {
          xdelta /= num ; ydelta /= num ; zdelta /= num ; 
        }
        xdelta *= l_distance ; ydelta *= l_distance ; zdelta *= l_distance ; 
        if (x == Gx && y == Gy && z == Gz)
          printf("l_dist: node(%d,%d,%d) distance term (%2.2f,%2.2f,%2.2f)\n",
                 x, y, z, xdelta, ydelta, zdelta) ;
        gcamn->dx += xdelta ; gcamn->dy += ydelta ; gcamn->dz += zdelta ;
      }

  return(NO_ERROR) ;
}

static double
gcamDistanceEnergy(GCA_MORPH *gcam, MRI *mri)
{
  double          sse = 0.0, dx, dy, dz, error, node_sse, d0, d ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  width = gcam->width ; height = gcam->height ; depth = gcam->depth ; 
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        num = 0 ; node_sse = 0.0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ; xn = MAX(0,xn) ; xn = MIN(width-1,xn) ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ; yn = MAX(0,yn) ; yn = MIN(height-1,yn) ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              if (!xk && !yk && !zk)
                continue ;
              zn = z+zk ; zn = MAX(0,zn) ; zn = MIN(depth-1,zn) ;
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;
#if 0
              if (gcamn_nbr->label != gcamn->label)
                continue ;
#endif
              dx = gcamn->origx - gcamn_nbr->origx ;
              dy = gcamn->origy - gcamn_nbr->origy ;
              dz = gcamn->origz - gcamn_nbr->origz ;
              d0 = sqrt(dx*dx + dy*dy + dz*dz) ;
              dx = gcamn->x - gcamn_nbr->x ;
              dy = gcamn->y - gcamn_nbr->y ;
              dz = gcamn->z - gcamn_nbr->z ;
              d = sqrt(dx*dx + dy*dy + dz*dz) ;
              error = d0-d ;
              num++ ;
              node_sse += error*error ;
            }
          }
        }
        num = 1 ;
        if (num > 0)
          sse += node_sse/num ;
        if (x == Gx && y == Gy && z == Gz)
          printf("E_dist: node(%d,%d,%d) distance sse %2.3f\n",
                 x, y, z, node_sse/num) ;
      }

  return(sse) ;
}

#define AREA_NEIGHBORS 8
static int
gcamJacobianTerm(GCA_MORPH *gcam, MRI *mri, double l_jacobian, double ratio_thresh)
{
  int            i, j, k, num, xi, yi, zi, xk, yk, zk ;
  double         dx, dy, dz, norm, orig_area, ratio ;
	GCA_MORPH_NODE *gcamn ;

  
  for (num = i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
	gcamn = &gcam->nodes[i][j][k] ;
	orig_area = gcamn->orig_area ;
	if (FZERO(orig_area))
	  continue ;
	
	ratio = gcamn->area / orig_area ;
	if (ratio < ratio_thresh)
	{
	  num++ ;
	  for (xk = -1 ; xk <= 1 ; xk++)
	  {
	    xi = i+xk ; if (xi < 0) xi = 0 ; if (xi >= gcam->width) xi = gcam->width-1 ;
	    for (yk = -1 ; yk <= 1 ; yk++)
	    {
	      yi = i+yk ; if (yi < 0) yi = 0 ; if (yi >= gcam->height) yi = gcam->height-1 ;
	      for (zk = -1 ; zk <= 1 ; zk++)
	      {
		zi = i+zk ; if (zi < 0) zi = 0 ; if (zi >= gcam->depth) zi = gcam->depth-1 ;
		gcamn = &gcam->nodes[xi][yi][zi] ;
		/*								gcamn->dx = gcamn->dy = gcamn->dz = 0 ;*/
	      }
	    }
	  }
	}
      }
    }
  }
  
  if (DIAG_VERBOSE_ON)
    printf("%d nodes compressed more than %2.2f\n", num, ratio_thresh) ;
  
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
	gcamn = &gcam->nodes[i][j][k] ;
	if (FZERO(gcamn->orig_area))
	  continue ;
        gcamJacobianTermAtNode(gcam, mri, l_jacobian, i, j, k, &dx, &dy, &dz) ;
#if 1
	norm = sqrt(dx*dx + dy*dy + dz*dz) ;
	if (norm > 1)  /* don't let it get too big, otherwise it's the only thing that happens */
	{
	  dx /= norm ; dy /= norm ; dz /= norm ;
	}
#endif
        gcam->nodes[i][j][k].dx += dx ;
        gcam->nodes[i][j][k].dy += dy ;
        gcam->nodes[i][j][k].dz += dz ;
      }
    }
  }
  return(NO_ERROR) ;
}

static int
gcamAreaTerm(GCA_MORPH *gcam, MRI *mri, double l_area)
{
  int    i, j, k ;
  double dx, dy, dz ;
  
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamAreaTermAtNode(gcam, mri, l_area, i, j, k, &dx, &dy, &dz) ;
        gcam->nodes[i][j][k].dx += dx ;
        gcam->nodes[i][j][k].dy += dy ;
        gcam->nodes[i][j][k].dz += dz ;
      }
    }
  }
  return(NO_ERROR) ;
}

static int
gcamJacobianTermAtNode(GCA_MORPH *gcam, MRI *mri, double l_jacobian, 
                       int i, int j, int k, double *pdx, double *pdy, double *pdz)
{
  GCA_MORPH_NODE *gcamn, *gcamni, *gcamnj, *gcamnk ;
  float          delta, ratio ;
  int            n, width, height, depth, num, invert ;
  static VECTOR  *v_i = NULL, *v_j, *v_k, *v_j_x_k, *v_i_x_j,*v_k_x_i,*v_grad,
                 *v_tmp ;
  double         exponent, orig_area ;

  if (FZERO(l_jacobian))
  {
    *pdx = *pdy = *pdz = 0.0 ;
    return(NO_ERROR) ;
  }

  width = gcam->width ; height = gcam->height ; depth = gcam->depth ; 
  if (!v_i)   /* initialize */
  {
    v_i = VectorAlloc(3, MATRIX_REAL) ;
    v_j = VectorAlloc(3, MATRIX_REAL) ;
    v_k = VectorAlloc(3, MATRIX_REAL) ;
    v_grad = VectorAlloc(3, MATRIX_REAL) ;
    v_j_x_k = VectorAlloc(3, MATRIX_REAL) ;
    v_i_x_j = VectorAlloc(3, MATRIX_REAL) ;
    v_k_x_i = VectorAlloc(3, MATRIX_REAL) ;
    v_tmp = VectorAlloc(3, MATRIX_REAL) ;
  }
  else
  { V3_CLEAR(v_grad) ; }

  for (num = n = 0 ; n < AREA_NEIGHBORS ; n++)
  {
    /* assign gcamn pointers to appropriate nodes */
    invert = 1 ;
    switch (n)
    {
    default:
    case 0:    /* first do central node */
      if ((i+1 >= width) || (j+1 >= height) || (k+1 >= depth))
        continue ;
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i+1][j][k] ;
      gcamnj = &gcam->nodes[i][j+1][k] ;
      gcamnk = &gcam->nodes[i][j][k+1] ;
      break ;
    case 1:       /*  i-1 */
      if ((i == 0) || (j+1 >= height) || (k+1 >= depth))
        continue ;
      gcamn = &gcam->nodes[i-1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i-1][j+1][k] ;
      gcamnk = &gcam->nodes[i-1][j][k+1] ;
      break ;
    case 2:       /* j-1 */
      if ((i+1 >= width) || (j == 0) || (k+1 >= depth))
        continue ;
      gcamn = &gcam->nodes[i][j-1][k] ;
      gcamni = &gcam->nodes[i+1][j-1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j-1][k+1] ;
      break ;
    case 3:      /* k-1 */
      if ((i+1 >= width) || (j+1 >= height) || (k == 0))
        continue ;
      gcamn = &gcam->nodes[i][j][k-1] ;
      gcamni = &gcam->nodes[i+1][j][k-1] ;
      gcamnj = &gcam->nodes[i][j+1][k-1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    case 4:
      if ((i == 0) || (j == 0) || (k == 0))
        continue ;
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i-1][j][k] ;
      gcamnj = &gcam->nodes[i][j-1][k] ;
      gcamnk = &gcam->nodes[i][j][k-1] ;
      break ;
    case 5:       /*  i+1 */
      if ((i+1 >= width) || (j == 0) || (k == 0))
        continue ;
      invert = -1 ;
      gcamn = &gcam->nodes[i+1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i+1][j-1][k] ;
      gcamnk = &gcam->nodes[i+1][j][k-1] ;
      break ;
    case 6:       /* j+1 */
      if ((i == 0) || (j+1 >= height) || (k == 0))
        continue ;
      invert = -1 ;
      gcamn = &gcam->nodes[i][j+1][k] ;
      gcamni = &gcam->nodes[i-1][j+1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j+1][k-1] ;
      break ;
    case 7:      /* k+1 */
      if ((i == 0) || (j == 0) || (k+1 >= depth))
        continue ;
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k+1] ;
      gcamni = &gcam->nodes[i-1][j][k+1] ;
      gcamnj = &gcam->nodes[i][j-1][k+1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    }
    orig_area = gcamn->orig_area ;
    if (FZERO(orig_area))
      orig_area = 0.1 ;

    num++ ;

    /* compute cross products and area delta */
    GCAMN_SUB(gcamni, gcamn, v_i) ; GCAMN_SUB(gcamnj, gcamn, v_j) ; 
    GCAMN_SUB(gcamnk, gcamn, v_k) ;

#if 0
    delta = invert * (orig_area - gcamn->area) * gcam->exp_k ;
#endif

    ratio = gcamn->area / orig_area ;
    if (ratio < 0.1 && ratio > 0)
      DiagBreak() ;
    if (gcamn->area < 0)
      DiagBreak() ;

    exponent = gcam->exp_k*ratio ;
    if (exponent > MAX_EXP)
      exponent = MAX_EXP ;

    /* don't use -k, since we are moving in the negative gradient direction */
    delta = (invert * gcam->exp_k / orig_area) * (1.0 / (1.0+exp(exponent))) ;

    if (fabs(delta) > 10000)
      DiagBreak() ;
    
    /* compute cross-products and add the appropriate 
       (i.e. scaled by area difference) cross-products to the gradient */
    switch (n)
    {
    default:
    case 4:
    case 0:    /* first do central node */
      V3_CROSS_PRODUCT(v_j, v_k, v_j_x_k) ;
      V3_CROSS_PRODUCT(v_k, v_i, v_k_x_i) ;
      V3_CROSS_PRODUCT(v_i, v_j, v_i_x_j) ;
      V3_ADD(v_i_x_j, v_j_x_k, v_tmp) ;
      V3_ADD(v_k_x_i, v_tmp, v_tmp) ;
      V3_SCALAR_MUL(v_tmp, -delta, v_tmp) ;
      break ;
    case 5:       /*  i+1 */
    case 1:       /*  i-1 */
      V3_CROSS_PRODUCT(v_j, v_k, v_j_x_k) ;
      V3_SCALAR_MUL(v_j_x_k, delta, v_tmp) ;
      break ;
    case 6:      /* j+1 */
    case 2:      /* j-1 */
      V3_CROSS_PRODUCT(v_k, v_i, v_k_x_i) ;
      V3_SCALAR_MUL(v_k_x_i, delta, v_tmp) ;
      break ;
    case 7:      /* k+1 */
    case 3:      /* k-1 */
      V3_CROSS_PRODUCT(v_i, v_j, v_i_x_j) ;
      V3_SCALAR_MUL(v_i_x_j, delta, v_tmp) ;
      break ;
    }
    V3_ADD(v_tmp, v_grad, v_grad) ;
  }

  *pdx = l_jacobian*V3_X(v_grad) ;
  *pdy = l_jacobian*V3_Y(v_grad) ;
  *pdz = l_jacobian*V3_Z(v_grad) ;

  if (i == Gx && j == Gy && k == Gz)
  {
    gcamn = &gcam->nodes[i][j][k] ;
    printf("l_jaco: node(%d,%d,%d): area=%2.2f, orig_area=%2.2f, grad=(%2.3f,%2.3f,%2.3f)\n", 
           i, j, k, gcamn->area,gcamn->orig_area, *pdx, *pdy, *pdz) ;
    if (fabs(*pdx) > 0.02 || fabs(*pdy) > 0.02 || fabs(*pdz) > 0.02)
      DiagBreak() ;
  }
  return(NO_ERROR) ;
}
static int
gcamAreaTermAtNode(GCA_MORPH *gcam, MRI *mri, double l_area, 
                       int i, int j, int k, double *pdx, double *pdy, double *pdz)
{
  GCA_MORPH_NODE *gcamn, *gcamni, *gcamnj, *gcamnk ;
  float          delta ;
  int            n, width, height, depth, num, invert ;
  static VECTOR  *v_i = NULL, *v_j, *v_k, *v_j_x_k, *v_i_x_j,*v_k_x_i,*v_grad,
                 *v_tmp ;
  double         orig_area ;

  if (FZERO(l_area))
  {
    *pdx = *pdy = *pdz = 0.0 ;
    return(NO_ERROR) ;
  }

  width = gcam->width ; height = gcam->height ; depth = gcam->depth ; 
  if (!v_i)   /* initialize */
  {
    v_i = VectorAlloc(3, MATRIX_REAL) ;
    v_j = VectorAlloc(3, MATRIX_REAL) ;
    v_k = VectorAlloc(3, MATRIX_REAL) ;
    v_grad = VectorAlloc(3, MATRIX_REAL) ;
    v_j_x_k = VectorAlloc(3, MATRIX_REAL) ;
    v_i_x_j = VectorAlloc(3, MATRIX_REAL) ;
    v_k_x_i = VectorAlloc(3, MATRIX_REAL) ;
    v_tmp = VectorAlloc(3, MATRIX_REAL) ;
  }
  else
  { V3_CLEAR(v_grad) ; }

  for (num = n = 0 ; n < AREA_NEIGHBORS ; n++)
  {
    /* assign gcamn pointers to appropriate nodes */
    invert = 1 ;
    switch (n)
    {
    default:
    case 0:    /* first do central node */
      if ((i+1 >= width) || (j+1 >= height) || (k+1 >= depth))
        continue ;
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i+1][j][k] ;
      gcamnj = &gcam->nodes[i][j+1][k] ;
      gcamnk = &gcam->nodes[i][j][k+1] ;
      break ;
    case 1:       /*  i-1 */
      if ((i == 0) || (j+1 >= height) || (k+1 >= depth))
        continue ;
      gcamn = &gcam->nodes[i-1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i-1][j+1][k] ;
      gcamnk = &gcam->nodes[i-1][j][k+1] ;
      break ;
    case 2:       /* j-1 */
      if ((i+1 >= width) || (j == 0) || (k+1 >= depth))
        continue ;
      gcamn = &gcam->nodes[i][j-1][k] ;
      gcamni = &gcam->nodes[i+1][j-1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j-1][k+1] ;
      break ;
    case 3:      /* k-1 */
      if ((i+1 >= width) || (j+1 >= height) || (k == 0))
        continue ;
      gcamn = &gcam->nodes[i][j][k-1] ;
      gcamni = &gcam->nodes[i+1][j][k-1] ;
      gcamnj = &gcam->nodes[i][j+1][k-1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    case 4:
      if ((i == 0) || (j == 0) || (k == 0))
        continue ;
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k] ;
      gcamni = &gcam->nodes[i-1][j][k] ;
      gcamnj = &gcam->nodes[i][j-1][k] ;
      gcamnk = &gcam->nodes[i][j][k-1] ;
      break ;
    case 5:       /*  i+1 */
      if ((i+1 >= width) || (j == 0) || (k == 0))
        continue ;
      invert = -1 ;
      gcamn = &gcam->nodes[i+1][j][k] ;
      gcamni = &gcam->nodes[i][j][k] ;
      gcamnj = &gcam->nodes[i+1][j-1][k] ;
      gcamnk = &gcam->nodes[i+1][j][k-1] ;
      break ;
    case 6:       /* j+1 */
      if ((i == 0) || (j+1 >= height) || (k == 0))
        continue ;
      invert = -1 ;
      gcamn = &gcam->nodes[i][j+1][k] ;
      gcamni = &gcam->nodes[i-1][j+1][k] ;
      gcamnj = &gcam->nodes[i][j][k] ;
      gcamnk = &gcam->nodes[i][j+1][k-1] ;
      break ;
    case 7:      /* k+1 */
      if ((i == 0) || (j == 0) || (k+1 >= depth))
        continue ;
      invert = -1 ;
      gcamn = &gcam->nodes[i][j][k+1] ;
      gcamni = &gcam->nodes[i-1][j][k+1] ;
      gcamnj = &gcam->nodes[i][j-1][k+1] ;
      gcamnk = &gcam->nodes[i][j][k] ;
      break ;
    }
    orig_area = gcamn->orig_area ;
    if (FZERO(orig_area))
      orig_area = 0.1 ;

    num++ ;

    /* compute cross products and area delta */
    GCAMN_SUB(gcamni, gcamn, v_i) ; GCAMN_SUB(gcamnj, gcamn, v_j) ; 
    GCAMN_SUB(gcamnk, gcamn, v_k) ;

    delta = invert * (orig_area - gcamn->area) ;

    if (fabs(delta) > 10000)
      DiagBreak() ;
    
    /* compute cross-products and add the appropriate 
       (i.e. scaled by area difference) cross-products to the gradient */
    switch (n)
    {
    default:
    case 4:
    case 0:    /* first do central node */
      V3_CROSS_PRODUCT(v_j, v_k, v_j_x_k) ;
      V3_CROSS_PRODUCT(v_k, v_i, v_k_x_i) ;
      V3_CROSS_PRODUCT(v_i, v_j, v_i_x_j) ;
      V3_ADD(v_i_x_j, v_j_x_k, v_tmp) ;
      V3_ADD(v_k_x_i, v_tmp, v_tmp) ;
      V3_SCALAR_MUL(v_tmp, -delta, v_tmp) ;
      break ;
    case 5:       /*  i+1 */
    case 1:       /*  i-1 */
      V3_CROSS_PRODUCT(v_j, v_k, v_j_x_k) ;
      V3_SCALAR_MUL(v_j_x_k, delta, v_tmp) ;
      break ;
    case 6:      /* j+1 */
    case 2:      /* j-1 */
      V3_CROSS_PRODUCT(v_k, v_i, v_k_x_i) ;
      V3_SCALAR_MUL(v_k_x_i, delta, v_tmp) ;
      break ;
    case 7:      /* k+1 */
    case 3:      /* k-1 */
      V3_CROSS_PRODUCT(v_i, v_j, v_i_x_j) ;
      V3_SCALAR_MUL(v_i_x_j, delta, v_tmp) ;
      break ;
    }
    V3_ADD(v_tmp, v_grad, v_grad) ;
  }

  *pdx = l_area*V3_X(v_grad) ;
  *pdy = l_area*V3_Y(v_grad) ;
  *pdz = l_area*V3_Z(v_grad) ;

  if (i == Gx && j == Gy && k == Gz)
  {
    gcamn = &gcam->nodes[i][j][k] ;
    printf("l_area: node(%d,%d,%d): area=%2.2f, orig_area=%2.2f, grad=(%2.3f,%2.3f,%2.3f)\n", 
           i, j, k, gcamn->area,gcamn->orig_area, *pdx, *pdy, *pdz) ;
    if (fabs(*pdx) > 0.02 || fabs(*pdy) > 0.02 || fabs(*pdz) > 0.02)
      DiagBreak() ;
  }
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
gcamComputeMetricProperties(GCA_MORPH *gcam)
{
  double         area ;
  int            i, j, k, width, height, depth, num ;
  GCA_MORPH_NODE *gcamn, *gcamni, *gcamnj, *gcamnk ;
  VECTOR         *v_i, *v_j, *v_k ;

  v_i = VectorAlloc(3, MATRIX_REAL) ;
  v_j = VectorAlloc(3, MATRIX_REAL) ;
  v_k = VectorAlloc(3, MATRIX_REAL) ;
  width = gcam->width ; height = gcam->height ; depth = gcam->depth ; 
  gcam->neg = 0 ;
  for (i = 0 ; i < width ; i++)
  {
    for (j = 0 ; j < height ; j++)
    {
      for (k = 0 ; k < depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        num = 0 ;
        if ((i < width-1) && (j < height-1) && (k < depth-1))
        {
          num++ ;
          gcamni = &gcam->nodes[i+1][j][k] ;
          gcamnj = &gcam->nodes[i][j+1][k] ;
          gcamnk = &gcam->nodes[i][j][k+1] ;
          GCAMN_SUB(gcamni, gcamn, v_i) ; 
          GCAMN_SUB(gcamnj, gcamn, v_j) ; 
          GCAMN_SUB(gcamnk, gcamn, v_k) ;
          area = VectorTripleProduct(v_j, v_k, v_i) ;
        }
        else
          area = 0 ;
        if ((i > 0) && (j > 0) && (k > 0))  /* left-hand coordinate system */
        {
          num++ ;
          gcamni = &gcam->nodes[i-1][j][k] ;
          gcamnj = &gcam->nodes[i][j-1][k] ;
          gcamnk = &gcam->nodes[i][j][k-1] ;

          /* invert v_i so that coordinate system is right-handed */
          GCAMN_SUB(gcamn, gcamni, v_i) ; 
          GCAMN_SUB(gcamnj, gcamn, v_j) ; 
          GCAMN_SUB(gcamnk, gcamn, v_k) ;
          area += VectorTripleProduct(v_j, v_k, v_i) ;
        }
        if (num > 0)
          gcamn->area = area / (float)num ;
        else
          gcamn->area = 0 ;
        if ((area <= 0) && !FZERO(gcamn->orig_area))
          gcam->neg++ ;
      }
    }
  }

  VectorFree(&v_i) ; VectorFree(&v_j) ; VectorFree(&v_k) ;
  return(NO_ERROR) ;
}

static double 
gcamJacobianEnergy(GCA_MORPH *gcam, MRI *mri)
{
  double          sse = 0.0, delta, ratio, exponent ;
  int             i, j, k, width, height, depth ;
  GCA_MORPH_NODE *gcamn ;

  width = gcam->width ; height = gcam->height ; depth = gcam->depth ; 
  for (sse = 0.0, i = 0 ; i < width ; i++)
  {
    for (j = 0 ; j < height ; j++)
    {
      for (k = 0 ; k < depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
#if 0
        if ((gcamn->area <= 0) && !FZERO(gcamn->orig_area))
          gcam->neg++ ;
#endif
        /* scale up the area coefficient if the area of the current node is
           close to 0 or already negative */
        if (!FZERO(gcamn->orig_area))
          ratio = gcamn->area / gcamn->orig_area ;
        else
        {
          ratio = 1 ;
#if 0
          fprintf(stderr, "orig area = 0 at (%d, %d, %d)!!!\n", i, j, k) ;
#endif
        }
        exponent = -gcam->exp_k*ratio ;
        if (exponent > MAX_EXP)
          delta = 0.0 ;
        else
          delta = log(1+exp(exponent)) /*   / gcam->exp_k */ ;

        sse += delta * mri->thick ;
        if (!finitep(delta) || !finitep(sse))
          DiagBreak() ;
        if (i == Gx && j == Gy && k == Gz)
          printf("E_jaco: node(%d,%d,%d): area=%2.2f, error=%2.3f\n", i, j, k, gcamn->area,delta);
        if (!FZERO(delta))
          DiagBreak() ;
      }
    }
  }

  return(sse) ;
}

static double 
gcamAreaEnergy(GCA_MORPH *gcam, MRI *mri)
{
  double          sse = 0.0, error ;
  int             i, j, k, width, height, depth ;
  GCA_MORPH_NODE *gcamn ;

  width = gcam->width ; height = gcam->height ; depth = gcam->depth ; 
  for (sse = 0.0, i = 0 ; i < width ; i++)
  {
    for (j = 0 ; j < height ; j++)
    {
      for (k = 0 ; k < depth ; k++)
      {
        gcamn = &gcam->nodes[i][j][k] ;
        error = gcamn->area - gcamn->orig_area ;

        sse += (error*error) ;
        if (!finitep(error) || !finitep(sse))
          DiagBreak() ;
        if (i == Gx && j == Gy && k == Gz)
          printf("E_area: node(%d,%d,%d): area=%2.2f, error=%2.3f\n", i, j, k, gcamn->area,error);
      }
    }
  }

  return(sse) ;
}


MRI *
GCAMmorphFromAtlas(MRI *mri_in, GCA_MORPH *gcam, MRI *mri_morphed)
{
  int        width, height, depth, x, y, z,
             xm1, ym1, zm1, xp1, yp1, zp1 ;
  float      xd, yd, zd, dx, dy, dz, thick ;
  Real       weight, orig_val, val ;
  MRI        *mri_weights, *mri_ctrl, *mri_s_morphed ;

  width = mri_in->width ; height = mri_in->height ; depth = mri_in->depth ; 
  thick = mri_in->thick ;
#define SCALE_FACTOR 25.0f  /* so we can use byte images for float values */

  mri_weights = MRIalloc(width, height, depth, MRI_UCHAR) ;
  mri_s_morphed = MRIalloc(width, height, depth, MRI_SHORT) ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        /* compute voxel coordinates of this morph point */
        GCAMsampleMorph(gcam, (float)x*thick, (float)y*thick, (float)z*thick, 
                        &xd, &yd, &zd) ;
        xd /= thick ; yd /= thick ; zd /= thick ;  /* voxel coords */
#if 1
        MRIsampleVolumeType(mri_in, x, y, z, &orig_val, SAMPLE_NEAREST);
#else
        orig_val = MRIvox(mri_in, x, y, z) ;
#endif
        if (orig_val > 40)
          DiagBreak() ;

        /* now use trilinear interpolation */
        xm1 = (int)floor(xd) ; ym1 = (int)floor(yd) ; zm1 = (int)floor(zd) ; 
        xp1 = xm1 + 1 ; yp1 = ym1 + 1 ; zp1 = zm1 + 1 ;

        /* make sure they are within bounds */
        xm1 = mri_in->xi[xm1] ; ym1 = mri_in->yi[ym1] ; zm1 = mri_in->zi[zm1] ;
        xp1 = mri_in->xi[xp1] ; yp1 = mri_in->yi[yp1] ; zp1 = mri_in->zi[zp1] ;

        dx = xd - xm1 ; dy = yd - ym1 ; dz = zd - zm1 ;

        /* now compute contribution to each of 8 nearest voxels */
        weight = (1-dx) * (1-dy) * (1-dz) ;
        MRISvox(mri_s_morphed,xm1,ym1,zm1) += weight * orig_val ;
        MRIvox(mri_weights,xm1,ym1,zm1) += nint(weight * SCALE_FACTOR) ;

        weight = (dx) * (1-dy) * (1-dz) ;
        MRISvox(mri_s_morphed,xp1,ym1,zm1) += weight * orig_val ;
        MRIvox(mri_weights,xp1,ym1,zm1) += nint(weight * SCALE_FACTOR) ;
        
        weight = (1-dx) * (dy) * (1-dz) ;
        MRISvox(mri_s_morphed,xm1,yp1,zm1) += weight * orig_val ;
        MRIvox(mri_weights,xm1,yp1,zm1) += nint(weight * SCALE_FACTOR) ;
        
        weight = (1-dx) * (1-dy) * (dz) ;
        MRISvox(mri_s_morphed,xm1,ym1,zp1) += weight * orig_val ;
        MRIvox(mri_weights,xm1,ym1,zp1) += nint(weight * SCALE_FACTOR) ;
        
        weight = (dx) * (dy) * (1-dz) ;
        MRISvox(mri_s_morphed,xp1,yp1,zm1) += weight * orig_val ;
        MRIvox(mri_weights,xp1,yp1,zm1) += nint(weight * SCALE_FACTOR) ;
        
        weight = (dx) * (1-dy) * (dz) ;
        MRISvox(mri_s_morphed,xp1,ym1,zp1) += weight * orig_val ;
        MRIvox(mri_weights,xp1,ym1,zp1) += nint(weight * SCALE_FACTOR) ;
        
        weight = (1-dx) * (dy) * (dz) ;
        MRISvox(mri_s_morphed,xm1,yp1,zp1) += weight * orig_val ;
        MRIvox(mri_weights,xm1,yp1,zp1) += nint(weight * SCALE_FACTOR) ;
        
        weight = (dx) * (dy) * (dz) ;
        MRISvox(mri_s_morphed,xp1,yp1,zp1) += weight * orig_val ;
        MRIvox(mri_weights,xp1,yp1,zp1) += nint(weight * SCALE_FACTOR) ;
      }
    }
  }

  /* now normalize weights and values */
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        weight = (float)MRIvox(mri_weights,x,y,z) / SCALE_FACTOR ;
        if (!FZERO(weight))
        {
          val = (Real)MRISvox(mri_s_morphed,x,y,z) / weight ;
          if (val > 255.0)
            val = 255.0 ;
          MRISvox(mri_s_morphed,x,y,z) = (short)nint(val) ;
        }
      }
    }
  }

  /* copy from short image to BUFTYPE one */
  if (!mri_morphed)
    mri_morphed = MRIclone(mri_in, NULL) ;
  else
    MRIclear(mri_morphed) ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
				switch (mri_morphed->type)
				{
				case MRI_UCHAR:
					MRIvox(mri_morphed,x,y,z) = (BUFTYPE)MRISvox(mri_s_morphed,x,y,z) ;
					break ;
				case MRI_SHORT:
					MRISvox(mri_morphed,x,y,z) = MRISvox(mri_s_morphed,x,y,z) ;
					break ;
				case MRI_FLOAT:
					MRIFvox(mri_morphed,x,y,z) = (float)MRISvox(mri_s_morphed,x,y,z) ;
					break ;
				default:
          ErrorReturn(NULL, 
                      (ERROR_UNSUPPORTED, 
                       "GCAMmorphFromAtlas: unsupported volume type %d",
                       mri_morphed->type)) ;
          break ;
				}
      }
    }
  }

  MRIfree(&mri_s_morphed) ;

  /* run soap bubble to fill in remaining holes */
#if 0
  mri_ctrl = MRIclone(mri_in, NULL) ;
#else
  mri_ctrl = mri_weights ;
#endif
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        weight = (float)MRIvox(mri_weights,x,y,z) / SCALE_FACTOR ;
        if (weight > .1)
          MRIvox(mri_ctrl, x, y, z) = 1 ;
        else
          MRIvox(mri_ctrl, x, y, z) = 0 ;
      }
    }
  }

#if 0
  MRIfree(&mri_weights) ;
#endif
  MRIbuildVoronoiDiagram(mri_morphed, mri_ctrl, mri_morphed) ;
  MRIsoapBubble(mri_morphed, mri_ctrl, mri_morphed, 5) ;

  MRIfree(&mri_ctrl) ;
  return(mri_morphed) ;
}

MRI *
GCAMmorphToAtlas(MRI *mri_src, GCA_MORPH *gcam, MRI *mri_morphed)
{
  int        width, height, depth, x, y, z ;
  float      xd, yd, zd ;
  Real       val ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ; 
  if (!mri_morphed)
    mri_morphed = MRIclone(mri_src, NULL) ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        GCAMsampleMorph(gcam, (float)x*mri_src->thick, 
                        (float)y*mri_src->thick, (float)z*mri_src->thick, 
                        &xd, &yd, &zd) ;
        xd /= mri_src->thick ; yd /= mri_src->thick ; zd /= mri_src->thick ; 
        if (xd > -1 && yd > -1 && zd > 0 &&
            xd < width && yd < height && zd < depth)
          MRIsampleVolume(mri_src, xd, yd, zd, &val);
        else
          val = 0.0 ;
        switch (mri_morphed->type)
        {
        case MRI_UCHAR:
          MRIvox(mri_morphed,x,y,z) = val ;
          break ;
				case MRI_SHORT:
          MRISvox(mri_morphed,x,y,z) = val ;
          break ;
        case MRI_FLOAT:
          MRIFvox(mri_morphed,x,y,z) = val ;
          break ;
        default:
          ErrorReturn(NULL, 
                      (ERROR_UNSUPPORTED, 
                       "GCAMmorphToAtlas: unsupported volume type %d",
                       mri_morphed->type)) ;
          break ;
        }
      }
    }
  }
  return(mri_morphed) ;
}
static int
log_integration_parms(FILE *fp, GCA_MORPH_PARMS *parms)
{
  char  *cp, host_name[STRLEN] ;

  cp = getenv("HOST") ;
  if (cp)
    strcpy(host_name, cp) ;
  else
    strcpy(host_name, "unknown") ;

  if (!FZERO(parms->l_map))
    fprintf(fp,"l_map=%2.2f ", parms->l_map) ;
  if (!FZERO(parms->l_jacobian))
    fprintf(fp,"l_jacobian=%2.2f ", parms->l_jacobian) ;
  if (!FZERO(parms->l_label))
    fprintf(fp,"l_label=%2.2f ", parms->l_label) ;
  if (!FZERO(parms->l_distance))
    fprintf(fp,"l_dist=%2.2f ", parms->l_distance) ;
  if (!FZERO(parms->l_log_likelihood))
    fprintf(fp,"l_likelihood=%2.2f ", parms->l_log_likelihood) ;
  if (!FZERO(parms->l_smoothness))
    fprintf(fp,"l_smoothness=%2.2f ", parms->l_smoothness) ;
  if (!FZERO(parms->l_area))
    fprintf(fp,"l_area=%2.2f ", parms->l_area) ;
  
  fprintf(fp, 
	  "\ntol=%2.2e, dt=%2.2e, exp_k=%2.1f, momentum=%2.2f, levels=%d, niter=%d, "
	  "lbl_dist=%2.2f, avgs=%d, sigma=%2.1f,type=%d, relabel=%d, neg=%s\n", 
	  parms->tol, parms->dt, parms->exp_k,parms->momentum, parms->levels, 
	  parms->niterations, parms->label_dist,parms->navgs,parms->sigma,parms->integration_type,
	  parms->relabel, parms->noneg ? "no" : "yes") ;
  return(NO_ERROR) ;
}

#if 0
static int
gcamCropNegativeAreaNodeGradient(GCA_MORPH *gcam, float crop)
{
  int            i, j, k /*, num, xi, yi, zi, xk, yk, zk*/ ;
/*  double         dx, dy, dz, norm, orig_area, ratio ;*/
	GCA_MORPH_NODE *gcamn ;
  
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
	gcamn = &gcam->nodes[i][j][k] ;
	if (FZERO(gcamn->orig_area) || gcamn->area > 0.0)
	  continue ;
	gcamn->dx *= crop ;
	gcamn->dy *= crop ;
	gcamn->dz *= crop ;
#if 0
	for (xk = -1 ; xk <= 1 ; xk++)
	{
	  xi = i+xk ; if (xi < 0) xi = 0 ; if (xi >= gcam->width) xi = gcam->width-1 ;
	  for (yk = -1 ; yk <= 1 ; yk++)
	  {
	    yi = i+yk ; if (yi < 0) yi = 0 ; if (yi >= gcam->height) yi = gcam->height-1 ;
	    for (zk = -1 ; zk <= 1 ; zk++)
	    {
	      zi = i+zk ; if (zi < 0) zi = 0 ; if (zi >= gcam->depth) zi = gcam->depth-1 ;
	      gcamn = &gcam->nodes[xi][yi][zi] ;
	      gcamn->dx = gcamn->dy = gcamn->dz = 0 ;
	    }
	  }
	}
#endif
      }
    }
  }
  
  return(NO_ERROR) ;
}	
#endif
int
GCAMregisterLevel(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, GCA_MORPH_PARMS *parms)
{
  int             n, nsmall, done = 0, which = GCAM_INTEGRATE_OPTIMAL, max_small, increasing, good_step;
  double          rms, last_rms, pct_change, orig_dt, min_dt, orig_j, tol, last_pct_change ;
  GCA_MORPH_PARMS jacobian_parms ;

  max_small = parms->nsmall ;
  gcamClearMomentum(gcam) ;
  jacobian_parms = *parms ; 
  jacobian_parms.l_likelihood = jacobian_parms.l_label = jacobian_parms.l_distance = jacobian_parms.l_smoothness = 0.0 ;
  orig_dt = parms->dt ;
  nsmall = 0 ; pct_change = 0.0 ;
  if (parms->integration_type == GCAM_INTEGRATE_FIXED)
    increasing = 0/*1*/ ;  /* will be set to 0 if new step is smaller than last */
  else
    increasing = 0 ;  /* don't use it for optimal time-step stuff */
  if (parms->log_fp) 
  {
    fprintf(parms->log_fp, "GCAMregisterLevel: using voxel thickness %2.1f, navgs=%d, relabel=%d\n",
            mri->xsize, parms->navgs, parms->relabel) ;
    log_integration_parms(parms->log_fp, parms) ;
    fflush(parms->log_fp) ;
  }
  if (Gdiag & DIAG_SHOW)
  {
#if 0
    printf("GCAMregisterLevel: using voxel thickness %2.1f, navgs=%d\n",
           mri->xsize, parms->navgs) ;
#endif
    log_integration_parms(stdout, parms) ;
  }

  if (parms->write_iterations && (Gdiag & DIAG_WRITE) && parms->start_t == 0)
    write_snapshot(gcam, mri, parms, 0) ;

  last_rms = gcamComputeRMS(gcam, mri, parms) ;
  if (parms->log_fp)
  {
    fprintf(parms->log_fp, "%03d: dt=%2.3f, rms=%2.3f, neg=%d\n",
            0, 0.0f, last_rms, gcam->neg) ;
    fflush(parms->log_fp) ;
  }
  if (Gdiag & DIAG_SHOW)
    printf("%03d: dt=%2.3f, rms=%2.3f, neg=%d\n",
           0, 0.0f, last_rms, gcam->neg) ;
	orig_j = parms->l_jacobian ;
	tol = parms->tol ; good_step = 0 ;
  for (n = parms->start_t ; n < parms->start_t+parms->niterations ; n++)
  {
		/*		parms->l_jacobian = 1 ;*/
    gcamComputeGradient(gcam, mri, mri_smooth, parms) ;
    parms->l_jacobian = orig_j ;
    switch (parms->integration_type)
    {
    case GCAM_INTEGRATE_OPTIMAL:
      min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
      parms->dt = min_dt ;
      break ;
    case GCAM_INTEGRATE_FIXED:
      min_dt = parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
      break ;
    case GCAM_INTEGRATE_BOTH:
      if (which == GCAM_INTEGRATE_OPTIMAL)
      {
	check_gcam(gcam) ;
	min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
	check_gcam(gcam) ;
	parms->dt = min_dt ;
	max_small = parms->nsmall ;
	tol = parms->tol ;
      }
      else  /* take some momentum steps */
      {
	max_small = 2*parms->nsmall ;
	tol = parms->tol/2 ;
	min_dt = parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
      }
      break ;
    default:
      min_dt = parms->dt ;
      ErrorExit(ERROR_BADPARM, "GCAMregisterLevel: unknown integration type %d",
		parms->integration_type) ;
    }

    gcamApplyGradient(gcam, parms) ;
    gcamComputeMetricProperties(gcam) ;
    if (gcam->neg > 0 && parms->noneg == True)
    {
      int i = 0 ;
      while (gcam->neg > 0)
      {
#if 0
	gcamCropNegativeAreaNodeGradient(gcam, 0.9) ;
#else
	parms->dt *= 0.9 ;
	if (Gdiag & DIAG_SHOW)
	  printf("%3.3d: %d negative nodes - reducing timestep to %2.6f...\n", 
		 i+1, gcam->neg, parms->dt) ;
	
#endif
	gcamUndoGradient(gcam) ;
	gcamComputeMetricProperties(gcam) ;
	gcamApplyGradient(gcam, parms) ;
	gcamComputeMetricProperties(gcam) ;
	if (++i > 50)
	{
	  if (gcam->neg > 0)  /* couldn't find a  step without folds */
	  {
	    gcamUndoGradient(gcam) ;
	    gcamComputeMetricProperties(gcam) ;
	  }
	  break ;
	}
      }
    }
    min_dt = parms->dt ;
    parms->dt = orig_dt ;
    
    gcamRemoveNegativeNodes(gcam, mri, parms) ;
    
    if (parms->write_iterations > 0 && !((n+1) % parms->write_iterations) && (Gdiag & DIAG_WRITE))
      write_snapshot(gcam, mri, parms, n+1) ;
    rms = gcamComputeRMS(gcam, mri, parms) ;
		last_pct_change = pct_change ;
    if (FZERO(last_rms))
      pct_change = 0.0 ;
    else
      pct_change = 100.0*(last_rms-rms)/last_rms ;
    if (pct_change < last_pct_change)
    {
      if (increasing)
	printf("pct change decreased\n") ;
      increasing = 0 ;
    }
    if (pct_change < 0)
    {
      printf("rms increased - undoing step...\n") ;
      done = 1 ;
      increasing = 0 ;
      gcamUndoGradient(gcam) ;
      gcamComputeMetricProperties(gcam) ;
      rms = gcamComputeRMS(gcam, mri, parms) ;
    }
    
    if (parms->log_fp)
    {
      fprintf(parms->log_fp, "%03d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d\n",
              n+1, min_dt, rms, pct_change, gcam->neg) ;
      fflush(parms->log_fp) ;
    }

    if (Gdiag & DIAG_SHOW)
      printf("%03d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d\n",
             n+1, min_dt, rms, pct_change, gcam->neg) ;


    if (pct_change < tol && !increasing)
    {
			if (max_small > 1)
	      printf("pct change < tol %2.3f, nsmall = %d of %d\n",
  	           tol, nsmall+1, max_small) ;
      if ((++nsmall >= max_small) || (pct_change < 0))
      {
	if (parms->integration_type == GCAM_INTEGRATE_BOTH)
	{
	  if (!good_step)
	    done++ ;
	  if (done >= 2)  /* couldn't take a step with either technique */
	  {
	    n++ ;
	    break ;
	  }
	  else
	  {
	    if (which == GCAM_INTEGRATE_FIXED)
	    {
	      increasing = 0 ;
	      which = GCAM_INTEGRATE_OPTIMAL ;
	    }
	    else
	    {
	      increasing = /*1 */0;
	      pct_change = 0.0 ;
	      which = GCAM_INTEGRATE_FIXED ;
	    }
	    printf("switching integration type to %s (done=%d)\n",
		   which == GCAM_INTEGRATE_FIXED ? "fixed" : "optimal", done) ;
	    good_step = nsmall = 0 ; gcamClearMomentum(gcam) ;
	  }
	}
	else
	{
	  n++ ;
	  break ;
	}
      }
    }
    else if (pct_change >= tol)    /* took at least one good step */
    {
      good_step = 1 ;
      done = 0 ;    /* for integration type == BOTH, apply both types again */
      nsmall = 0 ;  /* start counting small steps again */
    }
#if 0
    if (parms->relabel)
    {
      GCAMcomputeLabels(mri, gcam) ;
      last_rms = gcamComputeRMS(gcam, mri, parms) ;
    }
    else
#endif
      last_rms = rms ;
    /*		parms->label_dist -= 0.5 ;*/
    if (parms->label_dist < 0)
      parms->label_dist = 0 ;
  }

  parms->start_t = n ;
  parms->dt = orig_dt ;
  return(NO_ERROR) ;
}

static int
mriFillRegion(MRI *mri, int x, int y, int z, int fill_val, int whalf)
{
  int   xi, xk, yi, yk, zi, zk ;

  for (xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = mri->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = mri->zi[z+zk] ;
        MRIvox(mri, xi, yi, zi) = fill_val ;
      }
    }
  }
  return(NO_ERROR) ;
}
static int
write_snapshot(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms, int iter)
{
  char           fname[STRLEN], base_name[STRLEN] ;
  MRI            *mri_morphed, *mri_samples ;
  GCA_MORPH_NODE *gcamn ;
  int            x, y, z, xv, yv, zv ;
	static         int write_samples = -1, write_labels = -1 ;

  mri_morphed = GCAMmorphToAtlas(parms->mri, gcam, NULL) ;
  sprintf(base_name, "%s_%3.3d", parms->base_name, iter) ;
  sprintf(fname, "%s.mgh", base_name) ;
  printf("writing snapshot to %s\n", fname) ;
  MRIwrite(mri_morphed, fname) ;
  MRIwriteImageViews(mri_morphed, base_name, IMAGE_SIZE) ;
  MRIfree(&mri_morphed) ;

  mri_samples = MRIclone(parms->mri, NULL) ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        xv = mri_samples->xi[nint(gcamn->x)] ; 
        yv = mri_samples->yi[nint(gcamn->y)] ; 
        zv = mri_samples->zi[nint(gcamn->z)] ; 
        MRIvox(mri_samples, xv, yv, zv) = gcamn->label ;
      }
  if (write_samples < 0)
    write_samples = getenv("GCAM_WRITE_SAMPLES") != NULL ;
  if (write_labels < 0)
    write_labels = getenv("GCAM_WRITE_LABELS") != NULL ;
  
  if (write_samples > 0)
  {
    sprintf(fname, "%s_fsamples_%3.3d.mgh", parms->base_name, iter) ;
    printf("writing samples to %s....\n", fname) ;
    MRIwrite(mri_samples, fname) ;
  }
  
  if (write_labels > 0)
  {
    MRIclear(mri_samples) ;
    for (x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
	for (z = 0 ; z < gcam->depth ; z++)
	{
	  gcamn = &gcam->nodes[x][y][z] ;
#if 0
	  GCApriorToVoxel(gcam->gca, mri_samples, x, y, z, &xv, &yv, &zv) ;
#else
	  xv = mri_samples->xi[nint(gcamn->x)] ; 
	  yv = mri_samples->yi[nint(gcamn->y)] ; 
	  zv = mri_samples->zi[nint(gcamn->z)] ; 
#endif
	  if (gcamn->status & GCAM_IGNORE_LIKELIHOOD)
	    mriFillRegion(mri_samples, xv, yv, zv, gcamn->label, 1) ;
	}
    sprintf(fname, "%s_labels_%3.3d.mgh", parms->base_name, iter) ;
    printf("writing label samples to %s....\n", fname) ;
    MRIwrite(mri_samples, fname) ;
  }
  MRIfree(&mri_samples) ;
  
  return(NO_ERROR) ;
}

int
GCAMsampleMorph(GCA_MORPH *gcam, float x, float y, float z, 
                float *pxd, float *pyd, float *pzd)
{
  int            xm, xp, ym, yp, zm, zp, width, height, depth ;
  float          xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */

  /* x, y, z are in MRI voxel coords */
  x /= gcam->spacing ; y /= gcam->spacing ; z /= gcam->spacing ;
  width = gcam->width ; height = gcam->height ; depth = gcam->depth ;

  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;

  *pxd =
    xpd * ypd * zpd * gcam->nodes[xm][ym][zm].x +
    xpd * ypd * zmd * gcam->nodes[xm][ym][zp].x +
    xpd * ymd * zpd * gcam->nodes[xm][yp][zm].x +
    xpd * ymd * zmd * gcam->nodes[xm][yp][zp].x +
    xmd * ypd * zpd * gcam->nodes[xp][ym][zm].x +
    xmd * ypd * zmd * gcam->nodes[xp][ym][zp].x +
    xmd * ymd * zpd * gcam->nodes[xp][yp][zm].x +
    xmd * ymd * zmd * gcam->nodes[xp][yp][zp].x ;
  *pyd =
    xpd * ypd * zpd * gcam->nodes[xm][ym][zm].y +
    xpd * ypd * zmd * gcam->nodes[xm][ym][zp].y +
    xpd * ymd * zpd * gcam->nodes[xm][yp][zm].y +
    xpd * ymd * zmd * gcam->nodes[xm][yp][zp].y +
    xmd * ypd * zpd * gcam->nodes[xp][ym][zm].y +
    xmd * ypd * zmd * gcam->nodes[xp][ym][zp].y +
    xmd * ymd * zpd * gcam->nodes[xp][yp][zm].y +
    xmd * ymd * zmd * gcam->nodes[xp][yp][zp].y ;
  *pzd =
    xpd * ypd * zpd * gcam->nodes[xm][ym][zm].z +
    xpd * ypd * zmd * gcam->nodes[xm][ym][zp].z +
    xpd * ymd * zpd * gcam->nodes[xm][yp][zm].z +
    xpd * ymd * zmd * gcam->nodes[xm][yp][zp].z +
    xmd * ypd * zpd * gcam->nodes[xp][ym][zm].z +
    xmd * ypd * zmd * gcam->nodes[xp][ym][zp].z +
    xmd * ymd * zpd * gcam->nodes[xp][yp][zm].z +
    xmd * ymd * zmd * gcam->nodes[xp][yp][zp].z ;
  return(NO_ERROR) ;
}

static double
gcamComputeRMS(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  float   nvoxels ;
  double  rms, sse ;

check_gcam(gcam) ;
  sse = gcamComputeSSE(gcam, mri, parms) ;
check_gcam(gcam) ;
  nvoxels = gcam->width*gcam->height*gcam->depth ;
  rms = sqrt(sse/nvoxels) ;
  return(rms) ;
}

static double
gcamComputeSSE(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  double sse, l_sse, s_sse, j_sse, d_sse, a_sse, nvox, label_sse, map_sse ;

  nvox = gcam->width*gcam->height*gcam->width ;
  label_sse = map_sse = a_sse = sse = l_sse = s_sse = j_sse = d_sse = 0.0;

	check_gcam(gcam) ;
  gcamComputeMetricProperties(gcam) ;
	check_gcam(gcam) ;
	if (!FZERO(parms->l_log_likelihood) || !FZERO(parms->l_likelihood))
		l_sse = MAX(parms->l_log_likelihood, parms->l_likelihood) * gcamLogLikelihoodEnergy(gcam, mri) ;
	check_gcam(gcam) ;
	if (!FZERO(parms->l_label))
		label_sse = parms->l_label * gcamLabelEnergy(gcam, mri, parms->label_dist) ;
	check_gcam(gcam) ;
	if (!FZERO(parms->l_map))
		map_sse = parms->l_map * gcamMapEnergy(gcam, mri) ;
  if (!FZERO(parms->l_distance))
    d_sse = parms->l_distance * gcamDistanceEnergy(gcam, mri) ;
  if (!FZERO(parms->l_jacobian))
    j_sse = parms->l_jacobian * gcamJacobianEnergy(gcam, mri) ;
  if (!FZERO(parms->l_area))
    a_sse = parms->l_area * gcamAreaEnergy(gcam, mri) ;
  if (!FZERO(parms->l_smoothness))
    s_sse = parms->l_smoothness * gcamSmoothnessEnergy(gcam, mri) ;

  if (Gdiag & DIAG_SHOW)
  {
    if (!FZERO(parms->l_map))
      printf("map_sse = %2.2f ", map_sse/nvox) ;
    if (!FZERO(parms->l_log_likelihood))
      printf("l_sse = %2.2f ", l_sse/nvox) ;
    if (!FZERO(parms->l_distance))
      printf("d_sse = %2.2f ", d_sse/nvox) ;
    if (!FZERO(parms->l_jacobian))
      printf("j_sse = %2.2f ", j_sse/nvox) ;
    if (!FZERO(parms->l_smoothness))
      printf("s_sse = %2.2f ", s_sse/nvox) ;
    if (!FZERO(parms->l_area))
      printf("a_sse = %2.2f ", a_sse/nvox) ;
    printf("\n") ;
  }
  sse = l_sse + s_sse + j_sse + d_sse + a_sse + label_sse + map_sse ;
  return(sse) ;
}

static int
gcamComputeGradient(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, GCA_MORPH_PARMS *parms)
{
  static int i = 0 ;
  
  gcamClearGradient(gcam) ;
  gcamComputeMetricProperties(gcam) ;
  gcamMapTerm(gcam, mri, mri_smooth, parms->l_map)  ;
  gcamLabelTerm(gcam, mri, parms->l_label, parms->label_dist)  ;
  gcamLikelihoodTerm(gcam, mri, mri_smooth, parms->l_likelihood, parms)  ;
  gcamLogLikelihoodTerm(gcam, mri, mri_smooth, parms->l_log_likelihood)  ;
  gcamDistanceTerm(gcam, mri, parms->l_distance)  ;
  gcamAreaTerm(gcam, mri, parms->l_area)  ;
  gcamSmoothnessTerm(gcam, mri, parms->l_smoothness)  ;
  if (parms->write_iterations > 0 && (Gdiag & DIAG_WRITE) && getenv("GCAM_YGRAD") != NULL)
  {
    char fname[STRLEN] ;
    MRI *mri_grad, *mri_tmp ;
    
    mri_tmp = gcamWriteMRI(gcam, NULL, GCAM_Y_GRAD) ;
    mri_grad = MRIupsample2(mri_tmp, NULL) ;
    MRIfree(&mri_tmp) ;
    sprintf(fname, "%s_ygrad_before_%3.3d.mgh", parms->base_name, i) ;
    printf("writing smoothed y gradient to %s...\n", fname) ;
    MRIwrite(mri_grad, fname) ;
    MRIfree(&mri_grad) ;
    
  }
  gcamSmoothGradient(gcam, parms->navgs) ;
  if (parms->write_iterations > 0 && (Gdiag & DIAG_WRITE) && getenv("GCAM_YGRAD") != NULL)
  {
    char fname[STRLEN] ;
    MRI *mri_grad, *mri_tmp ;
    
    mri_tmp = gcamWriteMRI(gcam, NULL, GCAM_Y_GRAD) ;
    mri_grad = MRIupsample2(mri_tmp, NULL) ;
    MRIfree(&mri_tmp) ;
    sprintf(fname, "%s_ygrad_after_%3.3d.mgh", parms->base_name, i++) ;
    printf("writing smoothed y gradient to %s...\n", fname) ;
    MRIwrite(mri_grad, fname) ;
    MRIfree(&mri_grad) ;
    
  }
  gcamJacobianTerm(gcam, mri, parms->l_jacobian,parms->ratio_thresh)  ;

  gcamLimitGradientMagnitude(gcam, parms, mri) ;
  
  return(NO_ERROR) ;
}

static int
gcamLimitGradientMagnitude(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms, MRI *mri)
{
  int            x, y, z, xmax, ymax, zmax ;
  double         norm, max_norm /*, scale*/ ;
  GCA_MORPH_NODE *gcamn ;
  float          dt ;

  dt = parms->dt ;
  max_norm = 0.0 ; xmax = ymax = zmax = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
	if (FZERO(gcamn->orig_area))
	  gcamn->dx = gcamn->dy = gcamn->dz = 0 ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        norm = dt*sqrt(gcamn->dx*gcamn->dx+gcamn->dy*gcamn->dy+gcamn->dz*gcamn->dz) ;
#if 0
        if (norm > 3*parms->max_grad)
        {
          scale = 3*parms->max_grad / norm ;
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          gcamn->dx *= scale ;
          gcamn->dy *= scale ;
          gcamn->dz *= scale ;
          norm = dt*sqrt(gcamn->dx*gcamn->dx+gcamn->dy*gcamn->dy+gcamn->dz*gcamn->dz) ;
        }
#endif

        if (norm > max_norm)
        {
          max_norm = norm ;
          xmax = x ; ymax = y ; zmax = z ;
        }
      }
  {
    float vals[MAX_GCA_INPUTS] ;
    int   r ;
    gcamn = &gcam->nodes[xmax][ymax][zmax] ;
    load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->gca->ninputs) ;
    printf("max gradient %2.3f mm @ (%d, %d, %d), A=%2.1f, R=%2.3f, ", 
	   max_norm, xmax, ymax, zmax,
	   gcam->nodes[xmax][ymax][zmax].area,
	   gcam->nodes[xmax][ymax][zmax].area/gcam->nodes[xmax][ymax][zmax].orig_area) ;
    printf("vals(means) = ") ;
    for (r = 0 ; r < gcam->gca->ninputs ; r++)
      printf("%2.1f (%2.1f)  ", vals[r], gcamn->gc ? gcamn->gc->means[r] :-0.0);
    printf("\n") ;
}

#if 0
  if (max_norm > parms->max_grad)
  {
    scale = parms->max_grad / max_norm ;

    printf("scaling by %2.3f based on max gradient %2.3f mm @ (%d, %d, %d)\n",
           scale, max_norm, xmax, ymax, zmax) ;
    for (x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          gcamn->dx *= scale ;
          gcamn->dy *= scale ;
          gcamn->dz *= scale ;
        }
  }
#endif
  return(NO_ERROR) ;
}

static int
gcamSmoothnessTerm(GCA_MORPH *gcam, MRI *mri, double l_smoothness)
{
  double          vx, vy, vz, vnx, vny, vnz, dx, dy, dz ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  if (DZERO(l_smoothness))
    return(NO_ERROR) ;
  width = gcam->width ; height = gcam->height ; depth = gcam->depth ; 
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        vx = gcamn->x - gcamn->origx ;
        vy = gcamn->y - gcamn->origy ;
        vz = gcamn->z - gcamn->origz ;
        dx = dy = dz = 0.0f ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_smoo: node(%d,%d,%d): V=(%2.2f,%2.2f,%2.2f)\n",
                 x, y, z, vx, vy, vz) ;
        num = 0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ; xn = MAX(0,xn) ; xn = MIN(width-1,xn) ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ; yn = MAX(0,yn) ; yn = MIN(height-1,yn) ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              if (!zk && !yk && !xk)
                continue ;
              zn = z+zk ; zn = MAX(0,zn) ; zn = MIN(depth-1,zn) ;
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;
#if 0
              if (gcamn_nbr->label != gcamn->label)
                continue ;
#endif
              vnx = gcamn_nbr->x - gcamn_nbr->origx ;
              vny = gcamn_nbr->y - gcamn_nbr->origy ;
              vnz = gcamn_nbr->z - gcamn_nbr->origz ;
              dx += (vnx-vx) ; 
              dy += (vny-vy) ; 
              dz += (vnz-vz) ;
              if ((x == Gx && y == Gy && z == Gz) && 
                  (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
                printf("\tnode(%d,%d,%d): V=(%2.2f,%2.2f,%2.2f), DX=(%2.2f,%2.2f,%2.2f)\n",
                       xn, yn, zn, vnx, vny, vnz, vnx-vx, vny-vy, vnz-vz) ;
              num++ ;
            }
          }
        }
				/*        num = 1 ;*/
        if (num)
        {
          dx = dx * l_smoothness / num ;
          dy = dy * l_smoothness / num ;
          dz = dz * l_smoothness / num ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("l_smoo: node(%d,%d,%d): DX=(%2.2f,%2.2f,%2.2f)\n",
                 x, y, z, dx, dy, dz) ;
        gcamn->dx += dx ; gcamn->dy += dy ; gcamn->dz += dz ;
      }
  return(NO_ERROR) ;
}
static double
gcamSmoothnessEnergy(GCA_MORPH *gcam, MRI *mri)
{
  double          sse = 0.0, vx, vy, vz, vnx, vny, vnz, error, node_sse ;
  int             x, y, z, xk, yk, zk, xn, yn, zn, width, height, depth, num ;
  GCA_MORPH_NODE  *gcamn, *gcamn_nbr ;

  width = gcam->width ; height = gcam->height ; depth = gcam->depth ; 
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        vx = gcamn->x - gcamn->origx ;
        vy = gcamn->y - gcamn->origy ;
        vz = gcamn->z - gcamn->origz ;
        num = 0 ; node_sse = 0.0 ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xn = x+xk ; xn = MAX(0,xn) ; xn = MIN(width-1,xn) ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yn = y+yk ; yn = MAX(0,yn) ; yn = MIN(height-1,yn) ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              if (!xk && !yk && !zk)
                continue ;
              zn = z+zk ; zn = MAX(0,zn) ; zn = MIN(depth-1,zn) ;
              gcamn_nbr = &gcam->nodes[xn][yn][zn] ;
#if 0
              if (gcamn_nbr->label != gcamn->label)
                continue ;
#endif
              vnx = gcamn_nbr->x - gcamn_nbr->origx ;
              vny = gcamn_nbr->y - gcamn_nbr->origy ;
              vnz = gcamn_nbr->z - gcamn_nbr->origz ;
              error = SQR(vnx-vx)+SQR(vny-vy)+SQR(vnz-vz) ;
              num++ ;
              node_sse += error ;
            }
          }
        }
				/*        num = 1 ;*/
        if (num > 0)
          sse += node_sse/num ;
        if (x == Gx && y == Gy && z == Gz)
          printf("E_smoo: node(%d,%d,%d) smoothness sse %2.3f (%d nbrs)\n",
                 x, y, z, node_sse/num, num) ;
      }
  return(sse) ;
}

static int
gcamClearGradient(GCA_MORPH *gcam)
{
  int   x, y, z ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
        gcam->nodes[x][y][z].dx = gcam->nodes[x][y][z].dy = gcam->nodes[x][y][z].dz = 0.0;
  return(NO_ERROR) ;
}
static int
gcamApplyGradient(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms)
{
  int            x, y, z ;
  float          dx, dy, dz, dt, momentum ;
  GCA_MORPH_NODE *gcamn ;

  dt = parms->dt ; momentum = parms->momentum ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        dx = gcamn->dx*dt + gcamn->odx*momentum ; 
        dy = gcamn->dy*dt + gcamn->ody*momentum ; 
        dz = gcamn->dz*dt + gcamn->odz*momentum ; 
        gcamn->odx = dx ; gcamn->ody = dy ; gcamn->odz = dz ; 

        if (x == Gx && y == Gy && z == Gz)
          printf("GRAD: node(%d,%d,%d): moving by (%2.3f, %2.3f, %2.3f) from "
                 "(%2.1f,%2.1f,%2.1f) to ",
                 x, y, z, dx, dy, dz, gcamn->x, gcamn->y, gcamn->z) ;
        gcamn->x += dx; gcamn->y += dy; gcamn->z += dz;
        if (x == Gx && y == Gy && z == Gz)
          printf("(%2.1f,%2.1f,%2.1f)\n", gcamn->x, gcamn->y, gcamn->z) ;
      }
  return(NO_ERROR) ;
}
static int
gcamClearMomentum(GCA_MORPH *gcam)
{
  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        gcamn->odx = gcamn->ody = gcamn->odz = 0;

      }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
finitep(float f)
{
  return(1) ;
  if (!finite(f))
    return(0) ;
  if (fabs(f) > 1e5)
    return(0) ;
  return(1) ;
}
int
GCAMcomputeLabels(MRI *mri, GCA_MORPH *gcam)
{
  int            x, y, z, width, height, depth, label, n, nchanged = 0 ;
  float          vals[MAX_GCA_INPUTS] ;
  GCA_MORPH_NODE *gcamn ;
  GCA_PRIOR      *gcap ;
  GC1D           *gc ;

  width = gcam->width  ; height = gcam->height ; depth = gcam->depth ; 
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
      for (z = 0 ; z < depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
				load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->gca->ninputs) ;
        label = 
          GCAcomputeMAPlabelAtLocation(gcam->gca, x,y,z,vals,&n,&gcamn->log_p);
        gcap = &gcam->gca->priors[x][y][z] ;

        if (n >= 0)
        {
          if (label != gcamn->label)
            nchanged++ ;
          gcamn->label = label ;
          gcamn->n = n ;
          gcamn->prior = gcap->priors[n] ;
          gcamn->gc = gc = GCAfindPriorGC(gcam->gca, x, y, z, label) ;

          if (x == Gx && y == Gy && z == Gz)
            printf("RELABEL: node(%d, %d, %d): label %s (%d), mean %2.1f+-%2.1f, prior %2.1f, MRI=%2.0f\n",
                   x, y, z, cma_label_to_name(label), label,
		   gcamn->gc ? gcamn->gc->means[0] : 0.0, 
		   gcamn->gc ? sqrt(covariance_determinant(gcamn->gc, gcam->gca->ninputs)) : 0.0, 
                   gcamn->prior,vals[0]) ;
        }
        else  /* out of FOV probably */
        {
          gcamn->label = label ;
          gcamn->n = 0 ;
          gcamn->prior = 1.0 ;
					gcamn->gc = NULL ;
          if (x == Gx && y == Gy && z == Gz)
            printf("RELABEL: node(%d, %d, %d): label %s (%d), mean %2.1f+-%2.1f, prior %2.1f\n",
                   x, y, z, cma_label_to_name(label), label,
                   0.0, 0.0, gcamn->prior) ;
        }
      }

  printf("label assignment complete, %d changed (%2.2f%%)\n",
         nchanged, 100.0*(float)nchanged/(width*height*depth)) ;
  return(nchanged) ;
}
int
GCAMcomputeMaxPriorLabels(GCA_MORPH *gcam)
{
  int            x, y, z, width, height, depth, label, n, nchanged = 0,max_n ;
  GCA_MORPH_NODE *gcamn ;
  GC1D           *gc ;
  GCA_PRIOR      *gcap ;
  double         max_prior ;

  width = gcam->width  ; height = gcam->height ; depth = gcam->depth ; 
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
      for (z = 0 ; z < depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        gcap = &gcam->gca->priors[x][y][z] ;
        for (max_prior = -1, max_n = -1, n = 0 ; n < gcap->nlabels ; n++)
          if  (gcap->priors[n] >= max_prior)
          {
            max_prior = gcap->priors[n] ;
            max_n = n ;
          }
        n = max_n ; 
        if (n >= 0)
        {
          label = gcap->labels[n] ;
          if (label != gcamn->label)
            nchanged++ ;
          gcamn->label = label ;
          gcamn->n = n ;
          gcamn->prior = gcap->priors[n] ;
          gcamn->gc = gc = GCAfindPriorGC(gcam->gca, x, y, z, label) ;
        }
        else  /* out of FOV probably */
        {
          gcamn->label = label = 0 ;
          gcamn->n = 0 ;
          gcamn->prior = 1.0 ;
					gcamn->gc = NULL ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("MPRIOR: node(%d, %d, %d): label %s (%d), mean %2.1f+-%2.1f, prior %2.1f\n",
                 x, y, z, cma_label_to_name(label), label,
                 gcamn->gc ? gcamn->gc->means[0] : 0.0, 
		 gcamn->gc ? sqrt(covariance_determinant(gcamn->gc, gcam->gca->ninputs)) : 0.0, 
                 gcamn->prior) ;
      }

  printf("label assignment complete, %d changed (%2.2f%%)\n",
         nchanged, 100.0*(float)nchanged/(width*height*depth)) ;
  return(NO_ERROR) ;
}
MRI *
GCAMbuildMostLikelyVolume(GCA_MORPH *gcam, MRI *mri)
{
  int            x,  y, z, xn, yn, zn, width, depth, height, n ;
  GCA_MORPH_NODE *gcamn ;
	float          val ;

  width = mri->width ; depth = mri->depth ; height = mri->height ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if (!GCAvoxelToPrior(gcam->gca, mri, x, y, z, &xn, &yn, &zn))
	{
	  gcamn = &gcam->nodes[xn][yn][zn] ;
	  for (n = 0 ; n < gcam->gca->ninputs ; n++)
	  {
	    if (gcamn->gc)
	      val = gcamn->gc->means[n] ;
	    else
	      val = 0 ;
	    
	    switch (mri->type)
	    {
	    default:
	      ErrorReturn(NULL,
			  (ERROR_UNSUPPORTED, 
			   "GCAMbuildMostLikelyVolume: unsupported image type %d", mri->type)) ;
	      break ;
	    case MRI_SHORT:
	      MRISseq_vox(mri, x, y, z, n) = nint(val) ;
	      break ;
	    case MRI_UCHAR:
	      MRIseq_vox(mri, x, y, z, n) = nint(val) ;
	      break ;
	    case MRI_FLOAT:
	      MRIFseq_vox(mri, x, y, z, n) = val ;
	      break ;
	    }
	  }
        }// !GCA
      }
    }
  }

  return(mri) ;
}
MRI *
GCAMbuildVolume(GCA_MORPH *gcam, MRI *mri)
{
  int            x,  y, z, xn, yn, zn, width, depth, height, n ;
  GCA_MORPH_NODE *gcamn ;
	float          val ;

  width = mri->width ; depth = mri->depth ; height = mri->height ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if (!GCAvoxelToPrior(gcam->gca, mri, x, y, z, &xn, &yn, &zn))
	{
	  gcamn = &gcam->nodes[xn][yn][zn] ;
	  for (n = 0 ; n < gcam->gca->ninputs ; n++)
	  {
	    if (gcamn->gc)
	      val = gcamn->gc->means[n] ;
	    else
	      val = 0 ;
	    
	    switch (mri->type)
	    {
	    default:
	      ErrorReturn(NULL,
			  (ERROR_UNSUPPORTED, 
			   "GCAMbuildMostLikelyVolume: unsupported image type %d", mri->type)) ;
	      break ;
	    case MRI_SHORT:
	      MRISseq_vox(mri, x, y, z, n) = nint(val) ;
	      break ;
	    case MRI_UCHAR:
	      MRIseq_vox(mri, x, y, z, n) = nint(val) ;
	      break ;
	    case MRI_FLOAT:
	      MRIFseq_vox(mri, x, y, z, n) = val ;
	      break ;
	    }
	  }
        } // !GCA
      }
    }
  }

  return(mri) ;
}
int
GCAMinvert(GCA_MORPH *gcam, MRI *mri)
{
  int            x, y, z, width, height, depth, xv, yv, zv, num ;
  MRI            *mri_ctrl ;
  GCA_MORPH_NODE *gcamn ;

  if (gcam->mri_xind)   /* already inverted */
    return(NO_ERROR) ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;

  gcam->mri_xind = MRIalloc(width, height, depth, MRI_SHORT) ;
  gcam->mri_yind = MRIalloc(width, height, depth, MRI_SHORT) ;
  gcam->mri_zind = MRIalloc(width, height, depth, MRI_SHORT) ;
  mri_ctrl = MRIalloc(width, height, depth, MRI_UCHAR) ;
  if (!gcam->mri_xind || !gcam->mri_yind || !gcam->mri_zind || !mri_ctrl)
    ErrorExit(ERROR_NOMEMORY, "GCAMinvert: could not allocated %dx%dx%d index volumes",
              width, height, depth) ;

  for (z = 0 ; z < gcam->depth ; z++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (x = 0 ; x < gcam->width ; x++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        xv = nint(gcamn->x) ; yv = nint(gcamn->y) ; zv = nint(gcamn->z) ;
        if (xv < 0)
          xv = 0 ;
        if (yv < 0)
          yv = 0 ;
        if (zv < 0)
          zv = 0 ;
        if (xv >= width)
          xv = width-1 ;
        if (yv >= height)
          yv = height-1 ;
        if (zv >= depth)
          zv = depth-1 ;
        MRISvox(gcam->mri_xind, xv, yv, zv) += x ;
        MRISvox(gcam->mri_yind, xv, yv, zv) += y ;
        MRISvox(gcam->mri_zind, xv, yv, zv) += z ;
        MRIvox(mri_ctrl, xv, yv, zv) += 1 ;
      }
    }
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        num = MRIvox(mri_ctrl, x, y, z) ;
        if (num == 0)
          continue ;   /* nothing there */
        MRISvox(gcam->mri_xind, x, y, z) = 
          nint((float)MRISvox(gcam->mri_xind, x, y, z)/(float)num) ;
        MRISvox(gcam->mri_yind, x, y, z) = 
          nint((float)MRISvox(gcam->mri_yind, x, y, z)/(float)num) ;
        MRISvox(gcam->mri_zind, x, y, z) = 
          nint((float)MRISvox(gcam->mri_zind, x, y, z)/(float)num) ;
        MRIvox(mri_ctrl, x, y, z) = CONTROL_MARKED ;
      }
    }
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("performing soap bubble of x indices...\n") ;
  MRIbuildVoronoiDiagram(gcam->mri_xind, mri_ctrl, gcam->mri_xind) ;
  MRIsoapBubble(gcam->mri_xind, mri_ctrl, gcam->mri_xind, 5) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("performing soap bubble of y indices...\n") ;
  MRIbuildVoronoiDiagram(gcam->mri_yind, mri_ctrl, gcam->mri_yind) ;
  MRIsoapBubble(gcam->mri_yind, mri_ctrl, gcam->mri_yind, 5) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("performing soap bubble of z indices...\n") ;
  MRIbuildVoronoiDiagram(gcam->mri_zind, mri_ctrl, gcam->mri_zind) ;
  MRIsoapBubble(gcam->mri_zind, mri_ctrl, gcam->mri_zind, 5) ;
  MRIfree(&mri_ctrl) ;
  return(NO_ERROR) ;
}
int
GCAMfreeInverse(GCA_MORPH *gcam)
{
  if (gcam->mri_xind)
    MRIfree(&gcam->mri_xind) ;
  if (gcam->mri_yind)
    MRIfree(&gcam->mri_yind) ;
  if (gcam->mri_zind)
    MRIfree(&gcam->mri_zind) ;
  return(NO_ERROR) ;
}

static int
gcamUndoGradient(GCA_MORPH *gcam)
{
  int            x, y, z ;
  float          dx, dy, dz ;
  GCA_MORPH_NODE *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        dx = gcamn->odx ; dy  = gcamn->ody ; dz = gcamn->odz ; 
        gcamn->odx = gcamn->ody = gcamn->odz = 0 ;   /* turn off momentum */

        if (x == Gx && y == Gy && z == Gz)
          printf("UNGRAD: node(%d,%d,%d): moving by (%2.3f, %2.3f, %2.3f) from "
                 "(%2.1f,%2.1f,%2.1f) to ",
                 x, y, z, dx, dy, dz, gcamn->x, gcamn->y, gcamn->z) ;
        gcamn->x -= dx; gcamn->y -= dy; gcamn->z -= dz;
        if (x == Gx && y == Gy && z == Gz)
          printf("(%2.1f,%2.1f,%2.1f)\n", gcamn->x, gcamn->y, gcamn->z) ;
      }
  return(NO_ERROR) ;
}

static int
gcamReadMRI(GCA_MORPH *gcam, MRI *mri, int which)
{
  int            x, y, z ;
  float          d ;
  GCA_MORPH_NODE *gcamn ;
  
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
	gcamn = &gcam->nodes[x][y][z] ;
	d = MRIFvox(mri, x, y, z)  ;
	switch (which)
	{
	default:
	  ErrorExit(ERROR_BADPARM, "gcamReadMRI(%d): unknown which parameter", which) ;
	  break ;
	case GCAM_X_GRAD: gcamn->dx = d ; break ;
	case GCAM_Y_GRAD: gcamn->dy = d; break ;
	case GCAM_Z_GRAD: gcamn->dz = d; break ;
	}
      }
  return(NO_ERROR) ;
}
static MRI *
gcamWriteMRI(GCA_MORPH *gcam, MRI *mri, int which)
{
  int            x, y, z ;
  float          d ;
  GCA_MORPH_NODE *gcamn ;
  
  if (!mri)
    mri = MRIalloc(gcam->width, gcam->height, gcam->depth, MRI_FLOAT) ;
  if (!mri)
    ErrorExit(ERROR_NOMEMORY, "gcamWrite: could not allocate %dx%dx%d MRI\n",
	      gcam->width, gcam->height, gcam->depth) ;
  
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
	gcamn = &gcam->nodes[x][y][z] ;
	switch (which)
	{
	default:
	  ErrorExit(ERROR_BADPARM, "gcamWriteMRI(%d): unknown which parameter", which) ;
	  d = 0 ;
	  break ;
	case GCAM_X_GRAD: d = gcamn->dx ; break ;
	case GCAM_Y_GRAD: d = gcamn->dy ; break ;
	case GCAM_Z_GRAD: d = gcamn->dz ; break ;
	}
	MRIFvox(mri, x, y, z) = d ;
      }
  return(mri) ;
}

static int
gcamSmoothGradient(GCA_MORPH *gcam, int navgs)
{
#if 1
  MRI   *mri_tmp = NULL, *mri_kernel ;
  int   i ;
  
  if (navgs <= 0)     /* no smoothing */
    return(NO_ERROR) ;
  mri_kernel = MRIgaussian1d(sqrt((float)navgs)*M_PI/2.0, 0) ;
  if ((Gx >= 0 && Gy >= 0 && Gz >= 0) && (Gdiag & DIAG_SHOW))
    printf("before smoothing %d times: grad = (%2.3f, %2.3f, %2.3f)\n",
	   navgs, 
	   gcam->nodes[Gx][Gy][Gz].dx,
	   gcam->nodes[Gx][Gy][Gz].dy,
	   gcam->nodes[Gx][Gy][Gz].dz) ;
  
  for (i = GCAM_X_GRAD ; i <= GCAM_Z_GRAD ; i++)
  {
    mri_tmp = gcamWriteMRI(gcam, mri_tmp, i) ;
    MRIconvolveGaussian(mri_tmp, mri_tmp, mri_kernel) ;
    gcamReadMRI(gcam, mri_tmp, i) ;
  }
  
  MRIfree(&mri_tmp) ; MRIfree(&mri_kernel) ;
  
#else
  static double *dx, *dy, *dz = NULL ;
  int            x, y, z, xk, yk, zk, xi, yi, zi, index, i, num ;
  GCA_MORPH_NODE *gcamn, *gcamn_nbr ;
	
  if (navgs <= 0)
    return(NO_ERROR) ;
  
  if (!dz)
  {
    dx = (double *)calloc(gcam->width*gcam->height*gcam->depth, sizeof(double)) ;
    dy = (double *)calloc(gcam->width*gcam->height*gcam->depth, sizeof(double)) ;
    dz = (double *)calloc(gcam->width*gcam->height*gcam->depth, sizeof(double)) ;
    if (!dx || !dy || !dz)
      ErrorExit(ERROR_NOMEMORY, "gcamSmoothGradient: could not allocate grad buffers") ;
  }
  
  if ((Gx >= 0 && Gy >= 0 && Gz >= 0) && (Gdiag & DIAG_SHOW))
    printf("before smoothing %d times: grad = (%2.3f, %2.3f, %2.3f)\n",
	   navgs, 
	   gcam->nodes[Gx][Gy][Gz].dx,
	   gcam->nodes[Gx][Gy][Gz].dy,
	   gcam->nodes[Gx][Gy][Gz].dz) ;
  for (i = 0 ; i < navgs ; i++)
  {
    for (index = x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
	for (z = 0 ; z < gcam->depth ; z++, index++)
	{
	  if (x == Gx && y == Gy && z == Gz)
	    DiagBreak() ;
	  gcamn = &gcam->nodes[x][y][z] ;
	  dx[index] = dy[index] = dz[index] = 0 ;
#if 0
	  if (x == Gx && y == Gy && z == Gz)
	    printf("GRAD: node(%d,%d,%d): moving by (%2.3f, %2.3f, %2.3f) from "
		   "(%2.1f,%2.1f,%2.1f) to ",
		   x, y, z, dx, dy, dz, gcamn->x, gcamn->y, gcamn->z) ;
	  if (x == Gx && y == Gy && z == Gz)
	    printf("(%2.1f,%2.1f,%2.1f)\n", gcamn->x, gcamn->y, gcamn->z) ;
#endif
	  for (num = 0, xk = -1 ; xk <= 1 ; xk++)
	  {
	    xi = x+xk ; if (xi < 0) xi = 0 ; if (xi >= gcam->width) xi = gcam->width-1 ;
	    for (yk = -1 ; yk <= 1 ; yk++)
	    {
	      yi = y+yk ; if (yi < 0) yi = 0 ; if (yi >= gcam->height) yi = gcam->height-1 ;
	      for (zk = -1 ; zk <= 1 ; zk++)
	      {
		zi = z+zk ; if (zi < 0) zi = 0 ; if (zi >= gcam->depth) zi = gcam->depth-1 ;
		gcamn_nbr = &gcam->nodes[xi][yi][zi] ;
#if 1
		if (gcamn_nbr->label != gcamn->label)
		  continue ;
#endif
		num++ ;
		dx[index] += gcamn_nbr->dx ;
		dy[index] += gcamn_nbr->dy ;
		dz[index] += gcamn_nbr->dz ;
	      }
	    }
	  }
	  dx[index] /= (double)num ;
	  dy[index] /= (double)num ;
	  dz[index] /= (double)num ;
	}
    for (index = x = 0 ; x < gcam->width ; x++)
      for (y = 0 ; y < gcam->height ; y++)
	for (z = 0 ; z < gcam->depth ; z++, index++)
	{
	  gcamn = &gcam->nodes[x][y][z] ;
	  gcamn->dx = dx[index] ; gcamn->dy = dy[index] ; gcamn->dz = dz[index] ;
	  
	}
  }
  
#endif

  if ((Gx >= 0 && Gy >= 0 && Gz >= 0) && (Gdiag & DIAG_SHOW))
    printf("after smoothing %d times: grad = (%2.3f, %2.3f, %2.3f)\n",
	   navgs, 
	   gcam->nodes[Gx][Gy][Gz].dx,
	   gcam->nodes[Gx][Gy][Gz].dy,
	   gcam->nodes[Gx][Gy][Gz].dz) ;
  return(NO_ERROR) ;
}

int
GCAMsetStatus(GCA_MORPH *gcam, int status)
{

  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;
	
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	gcamn = &gcam->nodes[x][y][z] ;
	gcamn->status = status ;
      }
  return(NO_ERROR) ;
}
int
GCAMsetLabelStatus(GCA_MORPH *gcam, int label, int status)
{

  int            x, y, z ;
  GCA_MORPH_NODE *gcamn ;
	
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	gcamn = &gcam->nodes[x][y][z] ;
	if (gcamn->label == label)
	  gcamn->status = status ;
      }
  return(NO_ERROR) ;
}
#define MAX_SAMPLES 100
static double
gcamFindOptimalTimeStep(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms, MRI *mri)
{
  MATRIX   *mX, *m_xTx, *m_xTx_inv, *m_xTy, *mP, *m_xT ;
  double   min_dt, rms, min_rms, dt_in[MAX_SAMPLES], rms_out[MAX_SAMPLES],orig_dt,
    a, b, c, max_dt, start_dt ;
  VECTOR   *vY ;
  int      N, i, Gxs, Gys, Gzs ;
  long     diag ;
  
  gcamClearMomentum(gcam) ;
  Gxs = Gx ; Gys = Gy ; Gzs = Gz ; Gx = Gy = Gz = -1 ; diag = Gdiag ; Gdiag = 0 ;
  
  min_rms = gcamComputeRMS(gcam, mri, parms) ; min_dt = 0 ;
  
  /* first find right order of magnitude for time step */
  orig_dt = parms->dt ;
  
  /* define a pretty broad search range initially */
  start_dt = (sqrt(parms->navgs)+1)*orig_dt / (10*16.0*16.0) ;
  max_dt =   (sqrt(parms->navgs)+1)*orig_dt * (16.0*16.0) ;
  
  for (parms->dt = start_dt ; parms->dt <= max_dt ; parms->dt*=4)
  {
    gcamApplyGradient(gcam, parms) ;
    rms = gcamComputeRMS(gcam, mri, parms) ;
    gcamUndoGradient(gcam) ; 
    if (gcam->neg && parms->noneg == True)
      break ;
    if (rms < min_rms)
    {
      min_rms = rms ; min_dt = parms->dt ;
    }
  }
  
  dt_in[0] = min_dt ; rms_out[0] = min_rms ;
  parms->dt = dt_in[1] = dt_in[0] - dt_in[0]*.2 ; 
  gcamApplyGradient(gcam, parms) ; 
  rms = rms_out[1] = gcamComputeRMS(gcam, mri, parms) ;
  gcamUndoGradient(gcam) ; 
  if (rms < min_rms && (gcam->neg == 0 || parms->noneg == False))
  {
    min_rms = rms ; min_dt = parms->dt ;
  }
  
  parms->dt = dt_in[2] = dt_in[0] - dt_in[0]*.4 ; 
  gcamApplyGradient(gcam, parms) ; 
  rms = rms_out[2] = gcamComputeRMS(gcam, mri, parms) ;
  if (rms < min_rms && (gcam->neg == 0 || parms->noneg == False))
  {
    min_rms = rms ; min_dt = parms->dt ;
  }
  gcamUndoGradient(gcam) ;
  
  parms->dt = dt_in[3] = dt_in[0] + dt_in[0]*.2 ; 
  gcamApplyGradient(gcam, parms) ; 
  rms = rms_out[3] = gcamComputeRMS(gcam, mri, parms) ;
  gcamUndoGradient(gcam) ; 
  if (rms < min_rms && (gcam->neg == 0 || parms->noneg == False))
  {
    min_rms = rms ; min_dt = parms->dt ;
  }
  

  parms->dt = dt_in[4] = dt_in[0] + dt_in[0]*.4 ; 
  gcamApplyGradient(gcam, parms) ; 
  rms = rms_out[4] = gcamComputeRMS(gcam, mri, parms) ;
  gcamUndoGradient(gcam) ; 
  if (rms < min_rms && (gcam->neg == 0 || parms->noneg == False))
  {
    min_rms = rms ; min_dt = parms->dt ;
  }
  
  
  /* now compute location of minimum of best quadratic fit */
  N = 5 ;  /* min_dt +- .1*min_dt and .2*min_dt */
  mX = MatrixAlloc(N, 3, MATRIX_REAL) ;
  vY = VectorAlloc(N, MATRIX_REAL) ;   

  for (i = 1 ; i <= N ; i++)
  {
    *MATRIX_RELT(mX, i, 1) = dt_in[i-1] * dt_in[i-1] ;
    *MATRIX_RELT(mX, i, 2) = 2*dt_in[i-1] ;
    *MATRIX_RELT(mX, i, 3) = 1.0f ;
    
    VECTOR_ELT(vY, i) = rms_out[i-1] ;
  }
  
  m_xT = MatrixTranspose(mX, NULL) ;
  m_xTx = MatrixMultiply(m_xT, mX, NULL) ;
  m_xTx_inv = MatrixInverse(m_xTx, NULL) ;
  if (m_xTx_inv)
  {
    m_xTy = MatrixMultiply(m_xT, vY, NULL) ;
    mP = MatrixMultiply(m_xTx_inv, m_xTy, NULL) ;
    a = RVECTOR_ELT(mP, 1) ; b = RVECTOR_ELT(mP, 2) ; c = RVECTOR_ELT(mP, 3);
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout,
              "(a,b,c) = (%2.3f, %2.3f, %2.3f), predicted min at %2.3f\n",
              a, b, c, -b/a) ;
    if (!finite(a))
      DiagBreak() ;
    MatrixFree(&mP) ; 
    MatrixFree(&m_xTx_inv) ;
    MatrixFree(&m_xTy) ;
    if (finite(a) && !FZERO(a))
    {
      parms->dt = -b/a ;
      gcamApplyGradient(gcam, parms) ; 
      rms = gcamComputeRMS(gcam, mri, parms) ;
      gcamUndoGradient(gcam) ; 
      if (rms < min_rms && (gcam->neg == 0 || parms->noneg == False))
      {
	min_rms = rms ; min_dt = parms->dt ;
      }
    }
  }
  MatrixFree(&m_xT) ; MatrixFree(&m_xTx) ; MatrixFree(&mX) ; VectorFree(&vY) ;
  
  gcamComputeMetricProperties(gcam) ;
  parms->dt = orig_dt ;
  Gx = Gxs ; Gy = Gys ; Gz = Gzs ; Gdiag = diag ; 
  
  return(min_dt) ;
}

static double
gcamLabelEnergy(GCA_MORPH *gcam, MRI *mri, double label_dist)
{
  int             x, y, z, label, num, wm_label ;
  float           error, sse, vals[MAX_GCA_INPUTS], yi, yk, max_dist ;
  GCA_MORPH_NODE  *gcamn, *gcamn_inf, *gcamn_sup ;
	GC1D            *wm_gc ;
	double          log_p, max_log_p ;

	sse = 0.0 ;
  for (num = 0, x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
	if ((y == gcam->height-1) || (y == 0))
	  continue ;
	
	/* only process nodes which are hippocampus superior to something else,
	   or white matter inferior to hippocampus 
	*/
        gcamn = &gcam->nodes[x][y][z] ;
	if (y == 0 || y == gcam->height-1 || gcamn->y == 0 || gcamn->y == mri->height-1)
	  continue ;
	
        gcamn = &gcam->nodes[x][y][z] ;
				if ((gcamn->status & GCAM_LABEL_NODE) == 0)
					continue ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
	if ((y == gcam->height-1) || (y == 0))
	  continue ;
	if ((fabs(gcamn->x-Gvx)<=gcam->spacing) &&
	    (fabs(gcamn->y-Gvy)<=gcam->spacing) &&
	    (fabs(gcamn->z-Gvz)<=gcam->spacing))
	  DiagBreak() ;
	
	/* only process nodes which are hippocampus superior to something else,
	   or white matter inferior to hippocampus 
	*/
	if (y == 0 || y == gcam->height-1 || gcamn->y == 0 || gcamn->y == mri->height-1)
	  continue ;
	gcamn_inf = &gcam->nodes[x][y+1][z] ;
	gcamn_sup = &gcam->nodes[x][y-1][z] ;
	if (IS_HIPPO(gcamn->label))
	  load_vals(mri, gcamn->x, gcamn->y+1, gcamn->z, vals, gcam->gca->ninputs) ;
	else
	  load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->gca->ninputs) ;
	
	label = gcamMLElabelAtLocation(gcam, x, y, z, vals) ;
	if (IS_WM(label))  /* already has wm immediately inferior */
	  continue ;
	
	if ((gcamn->label == Left_Hippocampus) || 
	    (gcamn->label == Left_Cerebral_White_Matter))
	  wm_label = Left_Cerebral_White_Matter ;
	else
	  wm_label = Right_Cerebral_White_Matter ;
	wm_gc = GCAfindPriorGC(gcam->gca, x, y, z, wm_label) ;
	if (wm_gc == NULL)
	  continue ;
	
	max_dist = 0 ; max_log_p = -1e20;
	for (yk = -label_dist ; yk <= label_dist ; yk+=0.1) /* sample superiorly */
	{
	  yi = gcamn->y-yk ;   /* sample superiorly */
	  yi = gcamn->y+yk ;   /* sample inferiorly */
	  if (yi >= mri->height)
	    break ;
	  load_vals(mri, gcamn->x, yi, gcamn->z, vals, gcam->gca->ninputs) ;
	  log_p = GCAcomputeConditionalLogDensity(wm_gc, vals, gcam->gca->ninputs, wm_label);
	  if (log_p > max_log_p)
	  {
	    max_log_p = log_p ;
	    max_dist = yk ;
	  }
	}
	
	error = max_dist ;
	sse += (error*error) ;
	if (IS_HIPPO(gcamn_sup->label) && (gcamn_sup->status &GCAM_LABEL_NODE))
	  sse += (error*error) ;
	if (IS_WM(gcamn_inf->label) && (gcamn_inf->status &GCAM_LABEL_NODE))
	  sse += (error*error) ;
	
        if (x == Gx && y == Gy && z == Gz)
          printf("E_label: node(%d,%d,%d): %2.1f\n", x, y, z, error)  ;

	if (!FZERO(error))
	  num++ ;
      }
  
  return(sse) ;
}

static int
gcamLabelTerm(GCA_MORPH *gcam, MRI *mri, double l_label, double label_dist)
{
  int             x, y, z, label, wm_label, num = 0 ;
  Real            dy;
  GCA_MORPH_NODE  *gcamn, *gcamn_inf, *gcamn_sup ;
	float           vals[MAX_GCA_INPUTS], yk, max_dist, yi ;
	double          log_p, max_log_p ;
	GC1D            *wm_gc ;

  if (FZERO(l_label))
    return(NO_ERROR) ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
				
        gcamn = &gcam->nodes[x][y][z] ;
				gcamn->status = GCAM_USE_LIKELIHOOD ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if (x == Gx && y == (Gy-1) && z == Gz)
          DiagBreak() ;
        if (x == Gx && y == (Gy+1) && z == Gz)
          DiagBreak() ;

	if ((y == gcam->height-1) || (y == 0))
	  continue ;
	if ((fabs(gcamn->x-Gvx)<=gcam->spacing) &&
	    (fabs(gcamn->y-Gvy)<=gcam->spacing) &&
	    (fabs(gcamn->z-Gvz)<=gcam->spacing))
	  DiagBreak() ;
	
	/* only process nodes which are hippocampus superior to something else,
	   or white matter inferior to hippocampus 
	*/
	if (y == 0 || y == gcam->height-1 || gcamn->y == 0 || gcamn->y == mri->height-1)
	  continue ;
	if (!IS_HIPPO(gcamn->label) && !IS_WM(gcamn->label))
	  continue ;
	if (!IS_WM(gcamn->label))   /* only do white matter for now */
	  continue ;
	gcamn_inf = &gcam->nodes[x][y+1][z] ;
	gcamn_sup = &gcam->nodes[x][y-1][z] ;
	if (
	    ((IS_HIPPO(gcamn->label) && IS_WM(gcamn_inf->label)) ||
	     (IS_WM(gcamn->label) && IS_HIPPO(gcamn_sup->label))) == 0)
	  continue ;  /* only hippo above wm, or wm below hippo */
	if (IS_HIPPO(gcamn->label))
	  load_vals(mri, gcamn->x, gcamn->y+1, gcamn->z, vals, gcam->gca->ninputs) ;
	else
	  load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->gca->ninputs) ;
	
	label = gcamMLElabelAtLocation(gcam, x, y, z, vals) ;
	if (IS_WM(label))  /* already has wm immediately inferior */
	  continue ;
	
	if ((IS_HIPPO(gcamn->label) && gcamn->label == Left_Hippocampus) ||
	    (IS_WM(gcamn->label) && gcamn->label == Left_Cerebral_White_Matter))
	  wm_label = Left_Cerebral_White_Matter ;
	else
	  wm_label = Right_Cerebral_White_Matter ;
	wm_gc = GCAfindPriorGC(gcam->gca, x, y, z, wm_label) ;
	if (wm_gc == NULL)
	  continue ;
	
	dy = 0 ; 
	max_dist = 0 ; max_log_p = -1e20;
	for (yk = -label_dist ; yk <= label_dist ; yk += 0.1)
	{
	  yi = gcamn->y-yk ;   /* sample superiorly */
	  yi = gcamn->y+yk ;   /* sample inferiorly */
	  if ((yi >= (mri->height-1)) || (yi <= 0))
	    break ;
#if 0
	  gcamn_sup = &gcam->nodes[x][y-nint((float)yk/gcam->spacing)][z] ;
	  gcamn_inf = &gcam->nodes[x][y+nint((float)yk/gcam->spacing)][z] ;
	  if ((gcamn_sup != gcamn) && (!IS_HIPPO(gcamn_sup->label)))
	    break;
#endif
	  load_vals(mri, gcamn->x, yi, gcamn->z, vals, gcam->gca->ninputs) ;
	  log_p = GCAcomputeConditionalLogDensity(wm_gc, vals, gcam->gca->ninputs, wm_label);
	  if (log_p > max_log_p)
	  {
	    max_log_p = log_p ;
	    max_dist = yk ;
	  }
	  
#if 0
	  label = gcamMLElabelAtLocation(gcam, x, y, z, vals) ;
	  if (IS_WM(label))
	    dy -= 1 /*yk*/ ; /* move upwards proportional to distance to go */
#endif
	}
	
	dy = max_dist ;
	if (!FZERO(max_dist))
	{	
	  gcamn->status = (GCAM_IGNORE_LIKELIHOOD | GCAM_LABEL_NODE) ; 
	  if (fabs(max_dist) >= 1) 
	    num++ ; 
	  if (IS_WM(gcamn_inf->label) && ((gcamn_inf->status & GCAM_LABEL_NODE)==0))
	  {
	    gcamn_inf->status = (GCAM_IGNORE_LIKELIHOOD & GCAM_LABEL_NODE) ;
	    gcamn_inf->dy += (l_label)*dy ;
	    if (fabs(max_dist) >= 1) 
	      num++ ;
	    if (x == Gx && (y+1) == Gy && z == Gz)
	      printf("l_label: node(%d,%d,%d): dy = %2.2f\n", x, y+1, z, gcamn_inf->dy)  ;
	  }
	  if (IS_HIPPO(gcamn_sup->label) && ((gcamn_sup->status & GCAM_LABEL_NODE)==0))
	  {
	    gcamn_sup->status = (GCAM_IGNORE_LIKELIHOOD & GCAM_LABEL_NODE) ;
	    gcamn_sup->dy += (l_label)*dy ;
	    if (fabs(max_dist) >= 1) 
	      num++ ;
	    if (x == Gx && (y-1) == Gy && z == Gz)
	      printf("l_label: node(%d,%d,%d): dy = %2.2f\n", x, y-1, z, gcamn_sup->dy)  ;
	  }
	}
        gcamn->dy += (l_label)*dy ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_label: node(%d,%d,%d): dy = %2.2f\n", x, y, z, gcamn->dy)  ;
			}
  
#if 0
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
				
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
	if ((y == gcam->height-1) || (y == 0))
	  continue ;
	if ((fabs(gcamn->x-Gvx)<=gcam->spacing) &&
	    (fabs(gcamn->y-Gvy)<=gcam->spacing) &&
	    (fabs(gcamn->z-Gvz)<=gcam->spacing))
	  DiagBreak() ;
	
	/* 
	   set wm nodes  inferor to label nodes to also be label nodes so
	   that temporal white matter moves together 
	*/
	if (y == 0 || y == gcam->height-1 || gcamn->y == 0 || gcamn->y == mri->height-1)
	  continue ;
	gcamn_sup = &gcam->nodes[x][y-1][z] ;
      }
#endif
  
  
  if (Gdiag & DIAG_SHOW)
    printf("%d nodes for which label term applies\n", num) ;
  return(NO_ERROR) ;
}
static int
gcamMapTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, double l_map)
{
  int             x, y, z, n, xn, yn, zn, i ;
  double          node_prob, prob, dx, dy, dz, norm ;
  float           vals[MAX_GCA_INPUTS] ;
  GCA_MORPH_NODE  *gcamn ;
  GCA_PRIOR       *gcap ;
  GCA_NODE        *gcan ;
  GC1D            *gc ;
  MATRIX          *m_delI, *m_inv_cov ;
  VECTOR          *v_means, *v_grad ;
  
  if (FZERO(l_map))
    return(0) ;
  
  m_delI = MatrixAlloc(3, gcam->gca->ninputs, MATRIX_REAL) ;
  m_inv_cov = MatrixAlloc(gcam->gca->ninputs, gcam->gca->ninputs, MATRIX_REAL) ;
  v_means = VectorAlloc(gcam->gca->ninputs, 1) ;
  v_grad = VectorAlloc(3, MATRIX_REAL) ;
  
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
	if (gcamn->status & GCAM_IGNORE_LIKELIHOOD)
	  continue ;
	/////////////////////
	if (!GCApriorToNode(gcam->gca, x, y, z, &xn, &yn, &zn))
	{
	  gcan = &gcam->gca->nodes[xn][yn][zn] ;
	  gcap = &gcam->gca->priors[x][y][z] ;
	  load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->gca->ninputs) ;
	  
	  for (n = 0 ; n < gcam->gca->ninputs ; n++)
	  {
	    MRIsampleVolumeGradientFrame(mri_smooth, gcamn->x, gcamn->y, gcamn->z, &dx, &dy, &dz, n) ;
	    norm = sqrt(dx*dx+dy*dy+dz*dz) ;
	    if (!FZERO(norm))  /* don't worry about magnitude of gradient */
	    { dx /= norm ; dy /= norm ; dz /= norm ; }
	    *MATRIX_RELT(m_delI, 1, n+1) = dx ;
	    *MATRIX_RELT(m_delI, 2, n+1) = dy ;
	    *MATRIX_RELT(m_delI, 3, n+1) = dz ;
	  }
	  
	  if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
	    printf("l_map: node(%d,%d,%d), delI=(%2.2f, %2.2f, %2.2f)\n", x, y, z, dx, dy, dz) ;
	  dx = dy = dz = 0.0f ;
	  for (node_prob = 0.0, n = 0 ; n < gcap->nlabels ; n++)
	  {
	    gc = GCAfindGC(gcam->gca, xn, yn, zn, gcap->labels[n]) ;
	    if (!gc)
	      continue ;
	    
	    load_mean_vector(gc, v_means, gcam->gca->ninputs) ;
	    load_inverse_covariance_matrix(gc, m_inv_cov, gcam->gca->ninputs) ;
	    prob = GCAcomputeConditionalDensity(gc, vals, gcam->gca->ninputs, gcap->labels[n]) ;
	    
	    for (i = 0 ; i < gcam->gca->ninputs ; i++)
	      VECTOR_ELT(v_means, i+1) -= vals[i] ;
	    MatrixMultiply(m_inv_cov, v_means, v_means) ;
	    MatrixMultiply(m_delI, v_means, v_grad) ;
	    
	    if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
	      printf("l_map: node(%d,%d,%d), label %s: p=%2.3f (%2.3f), D=(%2.1f,%2.1f,%2.1f)\n",
		     x, y, z, cma_label_to_name(gcap->labels[n]), prob, gcap->priors[n], 
		     prob*V3_X(v_grad), prob*V3_Y(v_grad), prob*V3_Z(v_grad)) ;
	    dx += prob*V3_X(v_grad) ;
	    dy += prob*V3_Y(v_grad) ;
	    dz += prob*V3_Z(v_grad) ;
	    node_prob += prob ;
	  }
	  
	  if (DZERO(node_prob))
	    node_prob = 1e-6 ;
	  
	  dx *= (1.0/node_prob) ;
	  dy *= (1.0/node_prob) ;
	  dz *= (1.0/node_prob) ;
	  if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
	    printf("l_map: node(%d,%d,%d), 1/p=%2.3f, grad=(%2.2f,%2.2f,%2.2f)\n", 
		   x, y, z, 1.0/node_prob, dx, dy, dz) ;
	  gcamn->dx += l_map * dx ; gcamn->dy += l_map * dy ; gcamn->dz += l_map * dz ;
	} //!GCA
      }
  
  MatrixFree(&m_delI) ; MatrixFree(&m_inv_cov) ; VectorFree(&v_means) ; VectorFree(&v_grad) ;
  return(NO_ERROR) ;
}
static double
gcamMapEnergy(GCA_MORPH *gcam, MRI *mri)
{
  int             x, y, z, n, xn, yn, zn ;
  double          sse, node_prob, prob ;
  GCA_MORPH_NODE  *gcamn ;
  GCA_PRIOR       *gcap ;
  GCA_NODE        *gcan ;
  GC1D            *gc ;
  float            vals[MAX_GCA_INPUTS] ;
  
  sse = 0.0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
	if (gcamn->status & GCAM_IGNORE_LIKELIHOOD)
	  continue ;
	if (!GCApriorToNode(gcam->gca, x, y, z, &xn, &yn, &zn))
	{
	  gcan = &gcam->gca->nodes[xn][yn][zn] ;
	  gcap = &gcam->gca->priors[x][y][z] ;
	  load_vals(mri, gcamn->x, gcamn->y, gcamn->z, vals, gcam->gca->ninputs) ;
	  for (node_prob = 0.0, n = 0 ; n < gcap->nlabels ; n++)
	  {
	    gc = GCAfindGC(gcam->gca, xn, yn, zn, gcap->labels[n]) ;
	    if (!gc)
	      continue ;
	    
	    prob = GCAcomputeConditionalDensity(gc, vals, gcam->gca->ninputs, gcap->labels[n]) ;
	    if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
	      printf("E_map: node(%d,%d,%d), label %s: p=%2.3f (%2.3f)\n", x, y, z, 
		     cma_label_to_name(gcap->labels[n]), prob, gcap->priors[n]) ;
	    node_prob += prob ;
	  }
	  
	  if (x == Gx && y == Gy && z == Gz && (Gdiag & DIAG_SHOW))
	    printf("E_map: node(%d,%d,%d), -log(p)=%2.3f\n", x, y, z, -log(node_prob)) ;
	  if (FZERO(node_prob))
	    node_prob = 1e-6 ;  /* something big?? */
	  sse -= log(node_prob) ;
	} // !GCA
      }

  return(sse) ;
}

static int
gcamMLElabelAtLocation(GCA_MORPH *gcam, int x, int y, int z, float *vals)
{
  int   max_n ;
  float log_p ;
  
  return(GCAcomputeMLElabelAtLocation(gcam->gca, x, y, z, vals, &max_n, &log_p)) ;
}

int
GCAMstoreMetricProperties(GCA_MORPH *gcam)
{
  int             x, y, z ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
        gcam->nodes[x][y][z].orig_area = gcam->nodes[x][y][z].area ;
	return(NO_ERROR) ;
}

int
GCAMcomputeOriginalProperties(GCA_MORPH *gcam)
{
  GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, SAVED_POSITIONS) ;
  GCAMcopyNodePositions(gcam, ORIGINAL_POSITIONS, CURRENT_POSITIONS) ;
  gcamComputeMetricProperties(gcam) ;
  GCAMstoreMetricProperties(gcam) ;
  GCAMcopyNodePositions(gcam, SAVED_POSITIONS, CURRENT_POSITIONS) ;
  gcamComputeMetricProperties(gcam) ;
  return(NO_ERROR) ;
}

int
GCAMcopyNodePositions(GCA_MORPH *gcam, int from, int to)
{
  int             x, y, z ;
  GCA_MORPH_NODE  *gcamn ;

  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
	switch (from)
	{
	default:
	  ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GCAMcopyNodePositions: unsupported from %d",from)) ;
	case ORIGINAL_POSITIONS:
	  switch (to)
	  {
	  default:
	    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GCAMcopyNodePositions: unsupported to %d",to)) ;
	  case SAVED_POSITIONS:
	    gcamn->xs = gcamn->origx ; gcamn->ys = gcamn->origy ; gcamn->zs = gcamn->origz ;
	    break ;
	  case CURRENT_POSITIONS:
	    gcamn->x = gcamn->origx ; gcamn->y = gcamn->origy ; gcamn->z = gcamn->origz ;
	    break ;
	    
	  }
	  break ;
	  
	case SAVED_POSITIONS:
	  switch (to)
	  {
	  default:
	    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GCAMcopyNodePositions: unsupported to %d",to)) ;
	  case ORIGINAL_POSITIONS:
	    gcamn->origx = gcamn->xs ; gcamn->origy = gcamn->ys ; gcamn->origz = gcamn->zs ;
	    break ;
	  case CURRENT_POSITIONS:
	    gcamn->x = gcamn->xs ; gcamn->y = gcamn->ys ; gcamn->z = gcamn->zs ;
	    break ;
	  }
	  break ;
	  
	case CURRENT_POSITIONS:
	  switch (to)
	  {
	  default:
	    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GCAMcopyNodePositions: unsupported to %d",to)) ;
	  case ORIGINAL_POSITIONS:
	    gcamn->origx = gcamn->x ; gcamn->origy = gcamn->y ; gcamn->origz = gcamn->z ;
	    break ;
	  case SAVED_POSITIONS:
	    gcamn->xs = gcamn->x ; gcamn->ys = gcamn->y ; gcamn->zs = gcamn->z ;
	    break ;
	  }
	  
	  break ;
	  
	}
      }
  return(NO_ERROR) ;
}


static int
gcamRemoveNegativeNodes(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  GCA_MORPH_PARMS saved_parms = *parms ;
  double          min_dt, orig_dt = parms->dt, rms, last_rms, pct_change ;
  int             old_neg, new_neg, i ;
  
  if (gcam->neg <= 0)
    return(NO_ERROR) ;
  
  parms->noneg = 0 ;
  parms->l_distance = parms->l_log_likelihood = parms->l_area = parms->l_smoothness = parms->l_label = 0.0 ;
  parms->navgs = 0 ; parms->l_area = 0.0 ; 
  parms->l_jacobian = 1 ; parms->dt = 0.1 ;
  
  last_rms = rms = gcamComputeRMS(gcam, mri, parms) ;
  printf("starting rms=%2.3f, neg=%d, removing folds in lattice....\n", rms, gcam->neg) ;
  new_neg = gcam->neg ; i = 0 ;
  do
  {
    old_neg = new_neg ;
    gcamComputeGradient(gcam, mri, mri, parms) ;
    min_dt = gcamFindOptimalTimeStep(gcam, parms, mri) ;
    parms->dt = min_dt ;
    gcamApplyGradient(gcam, parms) ;
    last_rms = rms ;
    rms = gcamComputeRMS(gcam, mri, parms) ;
    parms->dt = orig_dt ;
    new_neg = gcam->neg ;
    pct_change = 100.0*(last_rms-rms)/(last_rms) ;
    printf("iter %d, dt=%2.6f: new neg %d, old_neg %d, delta %d, rms=%2.3f\n", ++i, min_dt, new_neg, old_neg, old_neg-new_neg, rms) ;
  } while ((new_neg > 0) && (pct_change > parms->tol) && (i < 5)) ;
  
  *parms = *(&saved_parms) ;
  return(NO_ERROR) ;
}
static int
check_gcam(GCAM *gcam)
{
  int  x, y, z ;
  GCA_MORPH_NODE *gcamn ;
  
  for (z = 0 ; z < gcam->depth ; z++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (x = 0 ; x < gcam->width ; x++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
	if (gcamn->gc)
	  gcamn->gc->means[0] = gcamn->gc->means[0]*1;
      }
    }
  }
  return(0) ;
}
