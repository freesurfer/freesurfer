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

#define GCAM_VERSION   1.0

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
static int gcamLimitGradientMagnitude(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms) ;
static int gcamComputeGradient(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, 
                               GCA_MORPH_PARMS *parms) ;
static int gcamLikelihoodTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth,
                              double l_likelihood) ;
static int gcamClearGradient(GCA_MORPH *gcam) ;
static int gcamComputeMetricProperties(GCA_MORPH *gcam) ;
static double gcamLikelihoodEnergy(GCA_MORPH *gcam, MRI *mri) ;

static int gcamDistanceTerm(GCA_MORPH *gcam, MRI *mri, double l_distance) ;
static double gcamDistanceEnergy(GCA_MORPH *gcam, MRI *mri) ;

static int gcamSmoothnessTerm(GCA_MORPH *gcam, MRI *mri, double l_smoothness) ;
static double gcamSmoothnessEnergy(GCA_MORPH *gcam, MRI *mri) ;

static int gcamAreaTerm(GCA_MORPH *gcam, MRI *mri, double l_jacobian) ;
static double gcamAreaEnergy(GCA_MORPH *gcam, MRI *mri) ;

static int gcamJacobianTerm(GCA_MORPH *gcam, MRI *mri, double l_jacobian) ;
static double gcamJacobianEnergy(GCA_MORPH *gcam, MRI *mri) ;
static int gcamApplyGradient(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms) ;
static int gcamUndoGradient(GCA_MORPH *gcam) ;
#if 0
static int gcamSmoothGradient(GCA_MORPH *gcam, double sigma) ;
#endif

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
  int    level, i, level_steps ;
  MRI    *mri_smooth = NULL, *mri_kernel ;
  double sigma, pct_change, rms, last_rms ;

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

  sigma = parms->sigma ;
  GCAMcomputeMaxPriorLabels(gcam) ;
  for (level = parms->levels-1 ; level >= 0 ; level--)
  {
    if (FZERO(sigma))
      mri_smooth = MRIcopy(mri, mri_smooth) ;
    else
    {
      printf("blurring input image with Gaussian with sigma=%2.3f...\n", sigma) ;
      mri_kernel = MRIgaussian1d(sigma, 100) ;
      mri_smooth = MRIconvolveGaussian(mri, NULL, mri_kernel) ;
      MRIfree(&mri_kernel) ;
    }
    i = 0 ;
    do
    {
      if ((level != (parms->levels-1)) || (i > 0))
        GCAMcomputeLabels(mri, gcam) ;
      last_rms = gcamComputeRMS(gcam, mri, parms) ;
      level_steps = parms->start_t ;
      GCAMregisterLevel(gcam, mri, mri_smooth, parms) ;
      rms = gcamComputeRMS(gcam, mri, parms) ;
      level_steps = parms->start_t - level_steps ;   /* # of steps taken in GCAMregisterLevel */
      if (level_steps == 0)
        level_steps = 1 ;
      pct_change = 100.0*(last_rms-rms)/(level_steps*last_rms) ;
      printf("iter %d: last rms %2.3f, rms %2.3f, pct_change %2.3f/iter\n",
             i+1, last_rms, rms, pct_change) ;
      if (i++ > 4)
        break ;
    } while (pct_change > parms->tol) ;
    sigma /= 2 ;
  }
  MRIfree(&mri_smooth) ;

  if (parms->log_fp)
  {
    fclose(parms->log_fp) ;
    parms->log_fp = NULL ;
  }
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
GCAMinit(GCA_MORPH *gcam, MRI *mri, GCA *gca, TRANSFORM *transform)
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
          if (gcap->priors[n] > max_p)
          {
            max_n = n ;
            max_p = gcap->priors[n] ;
            max_label = gcap->labels[n] ;
          }
        }
        gcamn->xn = x ; gcamn->yn = y ; gcamn->zn = z ;
        GCApriorToSourceVoxelFloat(gca, mri, transform, x, y, z,
                             &ox, &oy, &oz) ;

        gcamn->x = gcamn->origx = ox ; 
        gcamn->y = gcamn->origy = oy ; 
        gcamn->z = gcamn->origz = oz ;
#if 1
        gcamn->label = max_label ;
        gcamn->n = max_n ;
        gcamn->prior = max_p ;
        gc = GCAfindPriorGC(gca, x, y, z, max_label) ;
        if (gc)
        {
          gcamn->std = gc->var ;
          gcamn->mean = gc->mean ;
        }
        else   /* probably out of field of view */
        {
          gcamn->std = 1 ;
          gcamn->mean = 0 ;
        }
        gcamn->log_p = 0 ;
#endif
        if (x == Gx && y == Gy && z == Gz)
        {
          printf("node(%d,%d,%d) --> MRI (%2.1f, %2.1f, %2.1f)\n",
                 x, y, z, ox, oy, oz) ;
          DiagBreak() ;
        }
      }
    }
  }

#if 0
  GCAMcomputeLabels(mri, gcam) ;
#else
  GCAMcomputeMaxPriorLabels(gcam) ;
#endif
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


static int
gcamLikelihoodTerm(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, double l_likelihood)
{
  int             x, y, z ;
  float           error ;
  Real            dx, dy, dz, val, norm;
  GCA_MORPH_NODE  *gcamn ;

  if (FZERO(l_likelihood))
    return(NO_ERROR) ;


  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        MRIsampleVolume(mri, gcamn->x, gcamn->y, gcamn->z, &val) ;
        error = (gcamn->mean - val) / gcamn->std ;
#define MAX_ERROR 5
        if (fabs(error) > MAX_ERROR)
          error = MAX_ERROR * (error / fabs(error)) ;
        MRIsampleVolumeGradient(mri_smooth, gcamn->x,gcamn->y,gcamn->z,&dx, &dy, &dz) ;
        norm = sqrt(dx*dx+dy*dy+dz*dz) ;
        if (!FZERO(norm))  /* don't worry about magnitude of gradient */
        { dx /= norm ; dy /= norm ; dz /= norm ; }
        gcamn->dx += (l_likelihood)*error*dx ;
        gcamn->dy += (l_likelihood)*error*dy ;
        gcamn->dz += (l_likelihood)*error*dz ;
        if (x == Gx && y == Gy && z == Gz)
          printf("l_like: node(%d,%d,%d): dI=(%2.1f,%2.1f,%2.1f), grad=(%2.1f,%2.1f,%2.1f), "
                 "node %2.2f+-%2.2f, MRI=%2.1f\n",
                 x, y, z, dx, dy, dz, gcamn->dx, gcamn->dy, gcamn->dz,
                 gcamn->mean, gcamn->std, val) ;
      }

  return(NO_ERROR) ;
}

static float last_sse[128][128][128];

static double
gcamLikelihoodEnergy(GCA_MORPH *gcam, MRI *mri)
{
  double          sse = 0.0 ;
  int             x, y, z, max_x, max_y, max_z ;
  Real            val, error, max_increase = 0, increase ;
  GCA_MORPH_NODE  *gcamn ;

  max_x = max_y = max_z = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        MRIsampleVolume(mri, gcamn->x, gcamn->y, gcamn->z, &val) ;
        error = (gcamn->mean - val) / gcamn->std ;
        if (x == Gx && y == Gy && z == Gz)
          printf("E_like: node(%d,%d,%d) -> (%2.1f,%2.1f,%2.1f), target=%2.1f+-%2.1f, val=%2.1f\n",
                 x, y, z, gcamn->x, gcamn->y, gcamn->z, gcamn->mean, gcamn->std, val) ;
        if (last_sse[x][y][z] < (.9*error*error) && !FZERO(last_sse[x][y][z]))
        {
          DiagBreak() ;
          increase = error*error - last_sse[x][y][z] ;
          if (increase > max_increase)
          {
            max_increase = increase ; max_x = x ; max_y = y ; max_z = z ;
          }
        }
        last_sse[x][y][z] = (error*error) ;
        sse += (error*error) ;
        if (!finite(sse))
          DiagBreak() ;
      }
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
              if (gcamn_nbr->label != gcamn->label)
                continue ;
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
              if (gcamn_nbr->label != gcamn->label)
                continue ;
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
gcamJacobianTerm(GCA_MORPH *gcam, MRI *mri, double l_jacobian)
{
  int    i, j, k ;
  double dx, dy, dz ;
  
  for (i = 0 ; i < gcam->width ; i++)
  {
    for (j = 0 ; j < gcam->height ; j++)
    {
      for (k = 0 ; k < gcam->depth ; k++)
      {
        gcamJacobianTermAtNode(gcam, mri, l_jacobian, 
                               i, j, k, &dx, &dy, &dz) ;
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

#if 0
static int
gcamSmoothGradient(GCA_MORPH *gcam, double sigma)
{
  return(NO_ERROR) ;
}
#endif

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
#if 0
        MRIsampleVolume(mri_in, xd, yd, zd, &orig_val);
#else
        orig_val = MRIvox(mri_in, x, y, z) ;
#endif
        if (orig_val > 40)
          DiagBreak() ;

        /* now use trilinear interpolation */
        xm1 = (int)floor(xd) ; ym1 = (int)floor(yd) ; zm1 = (int)floor(zd) ; 
        xp1 = xm1 + 1 ; yp1 = ym1 + 1 ; zp1 = zm1 + 1 ;

        /* make sure they are within bounds */
        xm1 = mri_in->xi[xm1] ; ym1 = mri_in->xi[ym1] ; zm1 = mri_in->xi[zm1] ;
        xp1 = mri_in->xi[xp1] ; yp1 = mri_in->xi[yp1] ; zp1 = mri_in->xi[zp1] ;

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
        MRIvox(mri_morphed,x,y,z) = (BUFTYPE)MRISvox(mri_s_morphed,x,y,z) ;
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

  fprintf(fp,"l_area=%2.3f, l_distance=%2.3f, l_likelihood=%2.3f, l_smoothness=%2.3f, l_jacobian=%2.3f\n", 
          parms->l_area,
          parms->l_distance, parms->l_likelihood, parms->l_smoothness, parms->l_jacobian) ;
  fprintf(fp, "tol=%2.2e, dt=%2.2e, exp_k=%2.1f, momentum=%2.2f, levels=%d, niter=%d, host=%s\n", 
          parms->tol, parms->dt, parms->exp_k,parms->momentum, parms->levels, parms->niterations,
          host_name) ;
  return(NO_ERROR) ;
}

int
GCAMregisterLevel(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, GCA_MORPH_PARMS *parms)
{
  int     n, nsmall, niter ;
  double  rms, last_rms, pct_change, orig_dt ;

  orig_dt = parms->dt ;
  nsmall = 0 ;
  if (parms->log_fp) 
  {
    fprintf(parms->log_fp, "GCAMregisterLevel: using voxel thickness %2.1f\n",
            mri->xsize) ;
    log_integration_parms(parms->log_fp, parms) ;
    fflush(parms->log_fp) ;
  }
  if (Gdiag & DIAG_SHOW)
  {
    printf("GCAMregisterLevel: using voxel thickness %2.1f\n",
           mri->xsize) ;
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
  for (n = parms->start_t ; n < parms->start_t+parms->niterations ; n++)
  {
    gcamClearGradient(gcam) ;
    gcamComputeGradient(gcam, mri, mri_smooth, parms) ;
    niter = 0 ;
    do
    {
      gcamLimitGradientMagnitude(gcam, parms) ;
      gcamApplyGradient(gcam, parms) ;
      gcamComputeMetricProperties(gcam) ;
      if (gcam->neg > 0)
      {
        parms->dt *= .5 ;
        gcamUndoGradient(gcam) ;
      }
      else
        parms->dt = orig_dt ;
    } while ((gcam->neg > 0) && (niter++ < 5));
    if (parms->write_iterations > 0 &&
        !((n+1) % parms->write_iterations) && (Gdiag & DIAG_WRITE))
      write_snapshot(gcam, mri, parms, n+1) ;
    rms = gcamComputeRMS(gcam, mri, parms) ;
    if (FZERO(last_rms))
      pct_change = 0.0 ;
    else
      pct_change = 100.0*(last_rms-rms)/last_rms ;
    if (parms->log_fp)
    {
      fprintf(parms->log_fp, "%03d: dt=%2.3f, rms=%2.3f (%2.3f%%), neg=%d\n",
              n+1, parms->dt, rms, pct_change, gcam->neg) ;
      fflush(parms->log_fp) ;
    }

    if (Gdiag & DIAG_SHOW)
      printf("%03d: dt=%2.3f, rms=%2.3f (%2.3f%%), neg=%d\n",
             n+1, parms->dt, rms, pct_change, gcam->neg) ;
#define MAX_SMALL 5
    if (pct_change < parms->tol)
    {
      printf("pct change < tol %2.3f, nsmall = %d of %d\n",
             parms->tol, nsmall+1, MAX_SMALL) ;
      if (++nsmall >= MAX_SMALL)
      {
        n++ ;
        break ;
      }
    }
    else
      nsmall = 0 ;
#if 0
    GCAMcomputeLabels(mri, gcam) ;
    last_rms = gcamComputeRMS(gcam, mri, parms) ;
#else
    last_rms = rms ;
#endif
  }

  parms->start_t = n ;
  parms->dt = orig_dt ;
  return(NO_ERROR) ;
}

static int
write_snapshot(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms, int iter)
{
  char           fname[STRLEN], base_name[STRLEN] ;
  MRI            *mri_morphed, *mri_samples ;
  GCA_MORPH_NODE *gcamn ;
  int            x, y, z, xv, yv, zv ;

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
  sprintf(fname, "%s_fsamples_%3.3d.mgh", parms->base_name, iter) ;
  printf("writing samples to %s....\n", fname) ;
  MRIwrite(mri_samples, fname) ;
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

  sse = gcamComputeSSE(gcam, mri, parms) ;
  nvoxels = gcam->width*gcam->height*gcam->depth ;
  rms = sqrt(sse/nvoxels) ;
  return(rms) ;
}

static double
gcamComputeSSE(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms)
{
  double sse, l_sse, s_sse, j_sse, d_sse, a_sse, nvox ;

  nvox = gcam->width*gcam->height*gcam->width ;
  a_sse = sse = l_sse = s_sse = j_sse = d_sse = 0.0;

  gcamComputeMetricProperties(gcam) ;
  l_sse = parms->l_likelihood * gcamLikelihoodEnergy(gcam, mri) ;
  if (!FZERO(parms->l_distance))
    d_sse = parms->l_distance * gcamDistanceEnergy(gcam, mri) ;
  if (!FZERO(parms->l_jacobian))
    j_sse = parms->l_jacobian * gcamJacobianEnergy(gcam, mri) ;
  if (!FZERO(parms->l_area))
    a_sse = parms->l_area * gcamAreaEnergy(gcam, mri) ;
  if (!FZERO(parms->l_smoothness))
    s_sse = parms->l_smoothness * gcamSmoothnessEnergy(gcam, mri) ;

  printf("l_sse = %2.2f, d_sse = %2.2f, j_sse = %2.2f, s_sse = %2.2f, a_sse=%2.2f\n",
         l_sse/nvox, d_sse/nvox, j_sse/nvox, s_sse/nvox, a_sse/nvox) ;
  sse = l_sse + s_sse + j_sse + d_sse + a_sse ;
  return(sse) ;
}

static int
gcamComputeGradient(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, GCA_MORPH_PARMS *parms)
{
  gcamLikelihoodTerm(gcam, mri, mri_smooth, parms->l_likelihood)  ;
  gcamDistanceTerm(gcam, mri, parms->l_distance)  ;
  gcamJacobianTerm(gcam, mri, parms->l_jacobian)  ;
  gcamAreaTerm(gcam, mri, parms->l_area)  ;
  gcamSmoothnessTerm(gcam, mri, parms->l_smoothness)  ;
  
  return(NO_ERROR) ;
}

static int
gcamLimitGradientMagnitude(GCA_MORPH *gcam, GCA_MORPH_PARMS *parms)
{
  int            x, y, z, xmax, ymax, zmax ;
  double         norm, max_norm, scale ;
  GCA_MORPH_NODE *gcamn ;
  float          dt ;

  dt = parms->dt ;
  max_norm = 0.0 ; xmax = ymax = zmax = 0 ;
  for (x = 0 ; x < gcam->width ; x++)
    for (y = 0 ; y < gcam->height ; y++)
      for (z = 0 ; z < gcam->depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        norm = dt*sqrt(gcamn->dx*gcamn->dx+gcamn->dy*gcamn->dy+gcamn->dz*gcamn->dz) ;
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

        if (norm > max_norm)
        {
          max_norm = norm ;
          xmax = x ; ymax = y ; zmax = z ;
        }
      }

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
        num = 1 ;
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
        num = 1 ;
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
          printf("GRAD: node(%d,%d,%d): moving by (%2.2f, %2.2f, %2.2f) from "
                 "(%2.1f,%2.1f,%2.1f) to ",
                 x, y, z, dx, dy, dz, gcamn->x, gcamn->y, gcamn->z) ;
        gcamn->x += dx; gcamn->y += dy; gcamn->z += dz;
        if (x == Gx && y == Gy && z == Gz)
          printf("(%2.1f,%2.1f,%2.1f)\n", gcamn->x, gcamn->y, gcamn->z) ;
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
  Real           val ;
  GCA_MORPH_NODE *gcamn ;
  GCA_PRIOR      *gcap ;
  GC1D           *gc ;

  return(GCAMcomputeMaxPriorLabels(gcam)) ;
  width = gcam->width  ; height = gcam->height ; depth = gcam->depth ; 
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
      for (z = 0 ; z < depth ; z++)
      {
        gcamn = &gcam->nodes[x][y][z] ;
        MRIsampleVolume(mri, gcamn->x, gcamn->y, gcamn->z, &val) ;
        label = 
          GCAcomputeMAPlabelAtLocation(gcam->gca, x,y,z,val,&n,&gcamn->log_p);
        gcap = &gcam->gca->priors[x][y][z] ;
        if (n >= 0)
        {
          if (label != gcamn->label)
            nchanged++ ;
          gcamn->label = label ;
          gcamn->n = n ;
          gcamn->prior = gcap->priors[n] ;
          gc = GCAfindPriorGC(gcam->gca, x, y, z, label) ;
          if (gc)
          {
            gcamn->std = sqrt(gc->var) ;
            gcamn->mean = gc->mean ;
          }
          else   /* probably out of field of view */
          {
            gcamn->std = 1 ;
            gcamn->mean = 0 ;
          }

          if (x == Gx && y == Gy && z == Gz)
            printf("RELABEL: node(%d, %d, %d): label %s (%d), mean %2.1f+-%2.1f, prior %2.1f, MRI=%2.0f\n",
                   x, y, z, cma_label_to_name(label), label,
                   gcamn->mean, gcamn->std, gcamn->prior,val) ;
        }
        else  /* out of FOV probably */
        {
          gcamn->label = label ;
          gcamn->n = 0 ;
          gcamn->prior = 1.0 ;
          gcamn->std = 1.0 ;
          gcamn->mean = 0.0 ;
          if (x == Gx && y == Gy && z == Gz)
            printf("RELABEL: node(%d, %d, %d): label %s (%d), mean %2.1f+-%2.1f, prior %2.1f\n",
                   x, y, z, cma_label_to_name(label), label,
                   gcamn->mean, gcamn->std, gcamn->prior) ;
        }
      }

  printf("label assignment complete, %d changed (%2.2f%%)\n",
         nchanged, 100.0*(float)nchanged/(width*height*depth)) ;
  return(NO_ERROR) ;
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
          gc = GCAfindPriorGC(gcam->gca, x, y, z, label) ;
          if (gc)
          {
            gcamn->std = sqrt(gc->var) ;
            gcamn->mean = gc->mean ;
          }
          else   /* probably out of field of view */
          {
            gcamn->std = 1 ;
            gcamn->mean = 0 ;
          }
        }
        else  /* out of FOV probably */
        {
          gcamn->label = label = 0 ;
          gcamn->n = 0 ;
          gcamn->prior = 1.0 ;
          gcamn->std = 1.0 ;
          gcamn->mean = 0.0 ;
        }
        if (x == Gx && y == Gy && z == Gz)
          printf("RELABEL: node(%d, %d, %d): label %s (%d), mean %2.1f+-%2.1f, prior %2.1f\n",
                 x, y, z, cma_label_to_name(label), label,
                 gcamn->mean, gcamn->std, gcamn->prior) ;
      }

  printf("label assignment complete, %d changed (%2.2f%%)\n",
         nchanged, 100.0*(float)nchanged/(width*height*depth)) ;
  return(NO_ERROR) ;
}
MRI *
GCAMbuildMostLikelyVolume(GCA_MORPH *gcam, MRI *mri)
{
  int            x,  y, z, xn, yn, zn, width, depth, height ;
  GCA_MORPH_NODE *gcamn ;

  width = mri->width ; depth = mri->depth ; height = mri->height ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        GCAvoxelToPrior(gcam->gca, mri, x, y, z, &xn, &yn, &zn) ;
        gcamn = &gcam->nodes[xn][yn][zn] ;
        switch (mri->type)
        {
        default:
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED, 
                       "GCAbuildMostLikelyVolume: unsupported image type %d", mri->type)) ;
          break ;
        case MRI_UCHAR:
          MRIvox(mri, x, y, z) = nint(gcamn->mean) ;
          break ;
        case MRI_FLOAT:
          MRIFvox(mri, x, y, z) = gcamn->mean ;
          break ;
        }
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
          printf("UNGRAD: node(%d,%d,%d): moving by (%2.2f, %2.2f, %2.2f) from "
                 "(%2.1f,%2.1f,%2.1f) to ",
                 x, y, z, dx, dy, dz, gcamn->x, gcamn->y, gcamn->z) ;
        gcamn->x -= dx; gcamn->y -= dy; gcamn->z -= dz;
        if (x == Gx && y == Gy && z == Gz)
          printf("(%2.1f,%2.1f,%2.1f)\n", gcamn->x, gcamn->y, gcamn->z) ;
      }
  return(NO_ERROR) ;
}

