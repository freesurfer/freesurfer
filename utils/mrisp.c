#include <stdio.h>
#include <math.h>

#include "diag.h"
#include "error.h"
#include "utils.h"
#include "macros.h"
#include "mrisurf.h"
#include "proto.h"

/*---------------------------- STRUCTURES -------------------------*/

/*---------------------------- CONSTANTS -------------------------*/

#define BIG            100000.0

/* something of a hack... */
#define UNFILLED_ELT   0
#define FILLED_ELT     1

#define DEBUG_VNO -1
#define DEBUG_U  -1  /* -10 25*/
#define DEBUG_V  -1  /* 2 */

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
MRI_SP *
MRIStoParameterization(MRI_SURFACE *mris, MRI_SP *mrisp, float scale,int fno)
{
  float     a, b, c, phi, theta, x, y, z, total, uf, vf, d, du, dv, total_d,
            **distances, sigma, two_sigma_sq ;
  int       vno, u, v, uk, vk, n, unfilled, u0, v0, u1, v1, **filled, voff ;
  VERTEX    *vertex ;
float *dp ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "computing parameterization...") ;

  if (!mrisp)
    mrisp = MRISPalloc(scale, 1) ;

  a = b = c = mris->radius ;

  filled = (int **)calloc(U_DIM(mrisp), sizeof(int *)) ;
  distances = (float **)calloc(U_DIM(mrisp), sizeof(float *)) ;
  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    filled[u] = (int *)calloc(V_DIM(mrisp), sizeof(int)) ;
    distances[u] = (float *)calloc(V_DIM(mrisp), sizeof(float)) ;
    
    for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
      filled[u][v] = UNFILLED_ELT ;
  }

  sigma = scale / 4.0f ;
  two_sigma_sq = 2*sigma*sigma ;

  /* first calculate total distances to a point in parameter space */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(y/b, x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = c*c - z*z ;
    if (d < 0.0)
      d = 0 ;
    phi = atan2(sqrt(d), z) ;
    if (phi < RADIANS(1))
      DiagBreak() ;
    if (vno == DEBUG_VNO)
      DiagBreak() ;
    vertex->phi = phi ; vertex->theta = theta ;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u0 = floor(uf) ; u1 = ceil(uf+0.00001f) ;
    v0 = floor(vf) ; v1 = ceil(vf+0.00001f) ;
    du = uf - (float)u0 ; dv = vf - (float)v0 ;
    voff = 0 ;
    if (u0 < 0)  /* enforce spherical topology  */
    {
      voff = V_DIM(mrisp)/2 ;
      u0 = -u0 ;
    }
    if (u0 >= U_DIM(mrisp))
    {
      voff = V_DIM(mrisp)/2 ;
      u0 = U_DIM(mrisp)-(u0-U_DIM(mrisp)+1) ;
    }
    if (u1 < 0)  /* enforce spherical topology  */
      u1 = -u1 ;
    if (u1 >= U_DIM(mrisp))
      u1 = U_DIM(mrisp)-(u1-U_DIM(mrisp)+1) ;
    if (v0 < 0)  /* enforce spherical topology  */
      v0 += V_DIM(mrisp) ;
    if (v0 >= V_DIM(mrisp))
      v0 -= V_DIM(mrisp) ;
    if (v1 < 0)  /* enforce spherical topology  */
      v1 += V_DIM(mrisp) ;
    if (v1 >= V_DIM(mrisp))
      v1 -= V_DIM(mrisp) ;

    if ((u0 == DEBUG_U || u1 == DEBUG_U) && (v0 == DEBUG_V || v1 == DEBUG_V))
      DiagBreak() ;

    d = du*du+dv*dv ;  /* 1-distance to vertex 0,0 */
    d = exp(-d/two_sigma_sq) ;
    if (u0 >= U_DIM(mrisp) || v0 >= V_DIM(mrisp))
      DiagBreak() ;
    filled[u0][v0] = vno ;
    distances[u0][v0] += d ;
    if ((u0 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;

    d = 1.0 - dv ;         /* distance to v1 */
    d = du*du+d*d ;  /* distance to vertex 0,1 */
    d = exp(-d/two_sigma_sq) ;
    distances[u0][v1] += d ;         /* keep track of total distance */
    filled[u0][v1] = vno ;
    if ((u0 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;

    d = 1.0 - du ;         /* distance to u1 */
    d = d*d+dv*dv ;  /* distance to vertex 1,0 */
    d = exp(-d/two_sigma_sq) ;

    distances[u1][v0] += d ;         /* keep track of total distance */
    filled[u1][v0] = vno ;
    if ((u1 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;

    du = 1.0 - du ; dv = 1.0 - dv ;
    d = du*du+dv*dv ;  /* 1-distance to vertex 1,1 */
    d = exp(-d/two_sigma_sq) ;
    distances[u1][v1] += d ;         /* keep track of total distance */
    filled[u1][v1] = vno ;
    if ((u1 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;
  }

  if (DEBUG_U >= 0)
    fprintf(stderr,"\ndistance[%d][%d] = %2.3f\n\n",
            DEBUG_U, DEBUG_V, distances[DEBUG_U][DEBUG_V]);

dp = distances[37] ;

  /* now add in curvatures proportional to their distance from the point */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
if (dp != distances[37])
  DiagBreak() ;
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(y/b, x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = c*c-z*z ; if (d < 0.0) d = 0.0 ;
    phi = atan2(sqrt(d), z) ;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u0 = floor(uf) ; u1 = ceil(uf+0.00001f) ;
    v0 = floor(vf) ; v1 = ceil(vf+0.00001f) ;
    du = uf - (float)u0 ; dv = vf - (float)v0 ;
if (dp != distances[37])
  DiagBreak() ;
    if ((u0 == DEBUG_U || u1 == DEBUG_U) && (v0 == DEBUG_V || v1 == DEBUG_V))
      DiagBreak() ;
    if (u0 < 0)  /* enforce spherical topology  */
      u0 = -u0 ;
    if (u0 >= U_DIM(mrisp))
      u0 = U_DIM(mrisp)-(u0-U_DIM(mrisp)+1) ;
    if (u1 < 0)  /* enforce spherical topology  */
      u1 = -u1 ;
    if (u1 >= U_DIM(mrisp))
      u1 = U_DIM(mrisp)-(u1-U_DIM(mrisp)+1) ;
    if (v0 < 0)  /* enforce spherical topology  */
      v0 += V_DIM(mrisp) ;
    if (v0 >= V_DIM(mrisp))
      v0 -= V_DIM(mrisp) ;
    if (v1 < 0)  /* enforce spherical topology  */
      v1 += V_DIM(mrisp) ;
    if (v1 >= V_DIM(mrisp))
      v1 -= V_DIM(mrisp) ;

if (dp != distances[37])
  DiagBreak() ;
    /* 0,0 */
    d = du*du+dv*dv ;   /* distance to 0,0 */
    d = exp(-d/two_sigma_sq) ;
    total_d = distances[u0][v0] ;
    d /= total_d ;   /* weight by distance */
if (dp != distances[37])
  DiagBreak() ;
    *IMAGEFseq_pix(mrisp->Ip, u0,v0,fno) += d*vertex->curv ;  
if (dp != distances[37])
  DiagBreak() ;
    if ((u0 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

if (dp != distances[37])
  DiagBreak() ;
    /* 1,0 */
    d = 1.0 - du ;          /* distance to u1 */
    d = d*d+dv*dv ;   /* distance to u1,v0 */
    d = exp(-d/two_sigma_sq) ;
    total_d = distances[u1][v0] ;
    d = d / total_d ;  /* weight by distance */
    *IMAGEFseq_pix(mrisp->Ip, u1, v0, fno) += d*vertex->curv ;  
    if ((u1 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;
if (dp != distances[37])
  DiagBreak() ;

    /* 0,1 */
    d = 1.0 - dv ; 
    d = du*du+d*d ;   /* distance to u0,v1 */
    d = exp(-d/two_sigma_sq) ;
    total_d = distances[u0][v1] ;
    d = d / total_d ;    /* weight by distance */
    *IMAGEFseq_pix(mrisp->Ip, u0, v1, fno) += d*vertex->curv ;
    if ((u0 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;
if (dp != distances[37])
  DiagBreak() ;

    /* 1, 1 */
    du = 1.0 - du ; dv = 1.0 - dv ;
    d = (du*du+dv*dv) ;   /* distance to 1,1 */
    d = exp(-d/two_sigma_sq) ;
    total_d = distances[u1][v1] ;
    d = d / total_d ;   /* weight by distance */
    *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) += d*vertex->curv ;  
    if ((u1 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;
if (dp != distances[37])
  DiagBreak() ;

  }

  if (DEBUG_U >= 0)
    fprintf(stderr,"curv[%d][%d] = %2.3f\n\n", DEBUG_U, DEBUG_V, 
            *IMAGEFseq_pix(mrisp->Ip, DEBUG_U, DEBUG_V, fno));

#if 0
  /* normalize by total distances */
  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    for (v = 0 ; v <= V_MAX_INDEX ; v++)
    {
      if (u == DEBUG_U && v == DEBUG_V)
        DiagBreak() ;
      if (filled[u][v] >= 0)
      {
        *IMAGEFpix(mrisp->Ip, u,v) /= distances[u][v] ;
      }
    }
  }
#endif

  /* fill in values which were unmapped */
  do
  {
    unfilled = 0 ;
    for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
    {
      for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
      {
        if (filled[u][v] == UNFILLED_ELT)
        {
          for (total = 0.0f, n = 0, uk = -1 ; uk <= 1 ; uk++)
          {
            u1 = u + uk ;
            if (u1 < 0)  /* enforce spherical topology  */
              u1 = -u1 ;
            else if (u1 >= U_DIM(mrisp))
              u1 = U_DIM(mrisp)-(u1-U_DIM(mrisp)+1) ;
            for (vk = -1 ; vk <= 1 ; vk++)
            {
              v1 = v + vk ;
              if (v1 < 0)  /* enforce spherical topology  */
                v1 += V_DIM(mrisp) ;
              else if (v1 >= V_DIM(mrisp))
                v1 -= V_DIM(mrisp) ;
              if (filled[u1][v1] != UNFILLED_ELT)
              {
                total += *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) ;
                n++ ;
              }
            }
          }
          if (n)
          {
            filled[u][v] = FILLED_ELT ;
            *IMAGEFseq_pix(mrisp->Ip, u, v, fno) = total / (float)n ;
          }
          else
            unfilled++ ;
        }
      }
    }
  } while (unfilled > 0) ;

  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    free(filled[u]) ;
    free(distances[u]) ;
  }
  free(filled) ; free(distances) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "done.\n") ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    ImageWrite(mrisp->Ip, "sphere.hipl") ;

  /* keep track of the number of degrees of freedom in the last frame */
  *IMAGEFseq_pix(mrisp->Ip, 0, 0, mrisp->Ip->num_frame-1) += 1.0f ;
  return(mrisp) ;
}
#else
MRI_SP *
MRIStoParameterization(MRI_SURFACE *mris, MRI_SP *mrisp, float scale,int fno)
{
  float     a, b, c, phi, theta, x, y, z, uf, vf, d, total_d,
            **distances, sigma, two_sigma_sq ;
  int       vno, u, v, unfilled, **filled ;
  VERTEX    *vertex ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "computing parameterization...") ;

  if (!mrisp)
    mrisp = MRISPalloc(scale, 1) ;
  else
    ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, fno) ;

  a = b = c = mris->radius ;

  filled = (int **)calloc(U_DIM(mrisp), sizeof(int *)) ;
  distances = (float **)calloc(U_DIM(mrisp), sizeof(float *)) ;
  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    filled[u] = (int *)calloc(V_DIM(mrisp), sizeof(int)) ;
    distances[u] = (float *)calloc(V_DIM(mrisp), sizeof(float)) ;
    
    for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
      filled[u][v] = UNFILLED_ELT ;
  }

  sigma = scale / 4.0f ;
  two_sigma_sq = 2*sigma*sigma ;

  /* first calculate total distances to a point in parameter space */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(y/b, x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = c*c - z*z ;
    if (d < 0.0)
      d = 0 ;
    phi = atan2(sqrt(d), z) ;
    if (phi < RADIANS(1))
      DiagBreak() ;
    if (vno == DEBUG_VNO)
      DiagBreak() ;
    vertex->phi = phi ; vertex->theta = theta ;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u = nint(uf) ;
    v = nint(vf) ;
    if (u < 0)  /* enforce spherical topology  */
      u = -u ;
    if (u >= U_DIM(mrisp))
      u = U_DIM(mrisp) - (u-U_DIM(mrisp)+1) ;
    if (v < 0)  /* enforce spherical topology  */
      v += V_DIM(mrisp) ;
    if (v >= V_DIM(mrisp))
      v -= V_DIM(mrisp) ;

    if (u == 0 && v == 56)
      DiagBreak() ;
    if ((u == DEBUG_U) && (v == DEBUG_V))
      DiagBreak() ;

    filled[u][v] = vno ;
    distances[u][v] += 1 ;         /* keep track of total # of nodes */
    if ((u == DEBUG_U) && (v == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;
  }

  if (DEBUG_U >= 0)
    fprintf(stderr,"\ndistance[%d][%d] = %2.3f\n\n",
            DEBUG_U, DEBUG_V, distances[DEBUG_U][DEBUG_V]);


  /* now add in curvatures proportional to their distance from the point */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(y/b, x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = c*c-z*z ; if (d < 0.0) d = 0.0 ;
    phi = atan2(sqrt(d), z) ;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u = nint(uf) ;
    v = nint(vf) ;
    if (u < 0)  /* enforce spherical topology  */
      u = -u ;
    if (u >= U_DIM(mrisp))
      u = U_DIM(mrisp) - (u-U_DIM(mrisp)+1) ;
    if (v < 0)  /* enforce spherical topology  */
      v += V_DIM(mrisp) ;
    if (v >= V_DIM(mrisp))
      v -= V_DIM(mrisp) ;

    if (u == 0 && v == 56)
      DiagBreak() ;

    /* 0,0 */
    total_d = distances[u][v] ;
if ((total_d > 10000.0) || (vertex->curv > 1000.0))
  DiagBreak() ;
    if (total_d > 0.0)
      *IMAGEFseq_pix(mrisp->Ip, u,v,fno) += vertex->curv/total_d ;  
    if ((u == DEBUG_U) && (v == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

  }

  if (DEBUG_U >= 0)
    fprintf(stderr,"curv[%d][%d] = %2.3f\n\n", DEBUG_U, DEBUG_V, 
            *IMAGEFseq_pix(mrisp->Ip, DEBUG_U, DEBUG_V, fno));

#if 1
  /* fill in values which were unmapped by sampling back onto surface */
  for (unfilled = u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    double min_d, radius = mris->radius, xd, yd, zd ;
    int    min_v = -1 ;

    for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
    {
      if (filled[u][v] == UNFILLED_ELT)
      {
        if (u == 0 && v == 56)
          DiagBreak() ;
        unfilled++ ;
        phi   = (double)u*PHI_MAX / PHI_DIM(mrisp) ;
        theta = (double)v*THETA_MAX / THETA_DIM(mrisp) ;
        x = radius * sin(phi) * cos(theta) ;
        y = radius * sin(phi) * sin(theta) ;
        z = radius * cos(phi) ;
        for (min_d = 1000.0f, vno = 0 ; vno < mris->nvertices ; vno++)
        {
          vertex = &mris->vertices[vno] ;
          if (vertex->ripflag)
            continue ;
          xd = vertex->x - x ; yd = vertex->y - y ; zd = vertex->z - z ; 
          d = sqrt(xd*xd+yd*yd+zd*zd) ;
          if (d < min_d)
          {
            min_d = d ;
            min_v = vno ;
          }
        }
        *IMAGEFseq_pix(mrisp->Ip, u, v, fno) = mris->vertices[min_v].curv ;
      }
    }
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "%d holes in parameterization filled\n", unfilled) ;
#else
  /* fill in values which were unmapped using soap bubble */
  for (npasses = 0 ; npasses < 10 ; npasses++)
  {
    unfilled = 0 ; ;
    for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
    {
      for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
      {
        if (filled[u][v] == UNFILLED_ELT)
        {
          unfilled++ ;
          for (total = 0.0f, n = 0, uk = -1 ; uk <= 1 ; uk++)
          {
            u1 = u + uk ;
            if (u1 < 0)  /* enforce spherical topology  */
              u1 = -u1 ;
            else if (u1 >= U_DIM(mrisp))
              u1 = (u1-U_DIM(mrisp)+1)
            for (vk = -1 ; vk <= 1 ; vk++, n++)
            {
              v1 = v + vk ;
              if (v1 < 0)  /* enforce spherical topology  */
                v1 += V_DIM(mrisp) ;
              else if (v1 >= V_DIM(mrisp))
                v1 -= V_DIM(mrisp) ;
              
              total += *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) ;
            }
          }
          total /= (float)n ;
          *IMAGEFseq_pix(mrisp->Ip, u, v, fno) = total ;
        }
      }
    } 
  }
#endif

  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    free(filled[u]) ;
    free(distances[u]) ;
  }
  free(filled) ; free(distances) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "done.\n") ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    ImageWrite(mrisp->Ip, "sphere.hipl") ;

  return(mrisp) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISfromParameterization(MRI_SP *mrisp, MRI_SURFACE *mris, int fno)
{
  float     a, b, c, phi, theta, x, y, z, uf, vf, du, dv, curv, d ;
  int       vno, u0, v0, u1, v1 ;
  VERTEX    *vertex ;

  if (!mris)
    mris = MRISclone(mrisp->mris) ;

  a = b = c = mris->radius ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(vertex->y/b, vertex->x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = c*c-z*z ; if (d < 0.0) d = 0.0 ;
    phi = atan2(sqrt(d), z) ;
    if (phi < RADIANS(1))
      DiagBreak() ;
    if (vno == 60935)
      DiagBreak() ;

    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u0 = floor(uf) ; u1 = ceil(uf) ;
    v0 = floor(vf) ; v1 = ceil(vf) ;
    du = uf - (float)u0 ; dv = vf - (float)v0 ;

#if 0
    if (vno == 48092)
    {
      DiagBreak() ; 
      fprintf(stderr, "vno %d (%2.1f, %2.1f, %2.1f), u,v = (%2.1f, %2.1f)\n",
             vno, x, y, z, uf, vf) ;
    }
#endif

    /* enforce spherical topology  */
    if (u0 < 0)  /* enforce spherical topology  */
      u0 = -u0 ;
    if (u0 >= U_DIM(mrisp))
      u0 = U_DIM(mrisp)-(u0-U_DIM(mrisp)+1) ;
    if (u1 < 0)  /* enforce spherical topology  */
      u1 = -u1 ;
    if (u1 >= U_DIM(mrisp))
      u1 = U_DIM(mrisp)-(u1-U_DIM(mrisp)+1) ;
    if (v0 < 0) 
      v0 += V_DIM(mrisp) ;
    if (v0 >= V_DIM(mrisp))
      v0 -= V_DIM(mrisp) ;
    if (v1 < 0) 
      v1 += V_DIM(mrisp) ;
    if (v1 >= V_DIM(mrisp))
      v1 -= V_DIM(mrisp) ;

    if (((u0 == DEBUG_U) || (u1 == DEBUG_U)) &&
        ((v0 == DEBUG_V) || (v1 == DEBUG_V)))
      DiagBreak() ;
    curv = *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) ;
if (curv  > 100000.0)
  DiagBreak() ;
if (*IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) > 100000.0)
  DiagBreak() ;
    curv = 
      du        * dv        * *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) ;
if ((curv > 10000.0) || (vertex->curv > 1000.0))
  DiagBreak() ;
    curv =
      (1.0f-du) * dv        * *IMAGEFseq_pix(mrisp->Ip, u0, v1, fno) ;
if ((curv > 10000.0) || (vertex->curv > 1000.0))
  DiagBreak() ;
    curv = 
      (1.0f-du) * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, fno) ;
if ((curv > 10000.0) || (vertex->curv > 1000.0))
  DiagBreak() ;
    curv =
      du        * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, fno) ;
if ((curv > 10000.0) || (vertex->curv > 1000.0))
  DiagBreak() ;
    
    /* do bilinear interpolation */
    curv = 
      du        * dv        * *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) +
      (1.0f-du) * dv        * *IMAGEFseq_pix(mrisp->Ip, u0, v1, fno) +
      (1.0f-du) * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, fno) +
      du        * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, fno) ;

    vertex->curv = curv ;
if ((curv > 10000.0) || (vertex->curv > 1000.0))
  DiagBreak() ;
  }

  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISgradientToParameterization(MRI_SURFACE *mris, MRI_SP *mrisp, float scale)
{
  float     a, b, c, phi, theta, x, y, z, uf, vf, d, total_d,
            **distances, sigma, two_sigma_sq ;
  int       vno, u, v, unfilled, **filled ;
  VERTEX    *vertex ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "computing parameterization...") ;

  if (!mrisp)
    mrisp = MRISPalloc(scale, 3) ;
  else
    ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, -1) ;

  a = b = c = mris->radius ;

  filled = (int **)calloc(U_DIM(mrisp), sizeof(int *)) ;
  distances = (float **)calloc(U_DIM(mrisp), sizeof(float *)) ;
  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    filled[u] = (int *)calloc(V_DIM(mrisp), sizeof(int)) ;
    distances[u] = (float *)calloc(V_DIM(mrisp), sizeof(float)) ;
    
    for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
      filled[u][v] = UNFILLED_ELT ;
  }

  sigma = scale / 4.0f ;
  two_sigma_sq = 2*sigma*sigma ;

  /* first calculate total distances to a point in parameter space */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(y/b, x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = c*c - z*z ;
    if (d < 0.0)
      d = 0 ;
    phi = atan2(sqrt(d), z) ;
    if (phi < RADIANS(1))
      DiagBreak() ;
    if (vno == DEBUG_VNO)
      DiagBreak() ;
    vertex->phi = phi ; vertex->theta = theta ;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u = nint(uf) ;
    v = nint(vf) ;
    if (u < 0)  /* enforce spherical topology  */
      u = -u ;
    if (u >= U_DIM(mrisp))
      u = U_DIM(mrisp) - (u-U_DIM(mrisp)+1) ;
    if (v < 0)  /* enforce spherical topology  */
      v += V_DIM(mrisp) ;
    if (v >= V_DIM(mrisp))
      v -= V_DIM(mrisp) ;

    if (u == 0 && v == 56)
      DiagBreak() ;
    if ((u == DEBUG_U) && (v == DEBUG_V))
      DiagBreak() ;

    filled[u][v] = vno ;
    distances[u][v] += 1 ;         /* keep track of total # of nodes */
  }

  if (DEBUG_U >= 0)
    fprintf(stderr,"\ndistance[%d][%d] = %2.3f\n\n",
            DEBUG_U, DEBUG_V, distances[DEBUG_U][DEBUG_V]);


  /* now add in gradients proportional to their distance from the point */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(y/b, x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = c*c-z*z ; if (d < 0.0) d = 0.0 ;
    phi = atan2(sqrt(d), z) ;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u = nint(uf) ;
    v = nint(vf) ;
    if (u < 0)  /* enforce spherical topology  */
      u = -u ;
    if (u >= U_DIM(mrisp))
      u = U_DIM(mrisp) - (u-U_DIM(mrisp)+1) ;
    if (v < 0)  /* enforce spherical topology  */
      v += V_DIM(mrisp) ;
    if (v >= V_DIM(mrisp))
      v -= V_DIM(mrisp) ;

    if (u == 0 && v == 56)
      DiagBreak() ;

    /* 0,0 */
    total_d = distances[u][v] ;
    if (total_d > 0.0)
    {
      *IMAGEFseq_pix(mrisp->Ip, u,v,0) += vertex->dx/total_d ;  
      *IMAGEFseq_pix(mrisp->Ip, u,v,1) += vertex->dy/total_d ;  
      *IMAGEFseq_pix(mrisp->Ip, u,v,2) += vertex->dz/total_d ;  
    }
    if ((u == DEBUG_U) && (v == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

  }

  /* fill in values which were unmapped by sampling back onto surface */
  for (unfilled = u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    double min_d, radius = mris->radius, xd, yd, zd ;
    int    min_v = -1 ;

    for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
    {
      if (filled[u][v] == UNFILLED_ELT)
      {
        if (u == 0 && v == 56)
          DiagBreak() ;
        unfilled++ ;
        phi   = (double)u*PHI_MAX / PHI_DIM(mrisp) ;
        theta = (double)v*THETA_MAX / THETA_DIM(mrisp) ;
        x = radius * sin(phi) * cos(theta) ;
        y = radius * sin(phi) * sin(theta) ;
        z = radius * cos(phi) ;
        for (min_d = 1000.0f, vno = 0 ; vno < mris->nvertices ; vno++)
        {
          vertex = &mris->vertices[vno] ;
          if (vertex->ripflag)
            continue ;
          xd = vertex->x - x ; yd = vertex->y - y ; zd = vertex->z - z ; 
          d = sqrt(xd*xd+yd*yd+zd*zd) ;
          if (d < min_d)
          {
            min_d = d ;
            min_v = vno ;
          }
        }
        vertex = &mris->vertices[min_v] ;
        *IMAGEFseq_pix(mrisp->Ip, u, v, 0) = vertex->dx ;
        *IMAGEFseq_pix(mrisp->Ip, u, v, 1) = vertex->dy ;
        *IMAGEFseq_pix(mrisp->Ip, u, v, 2) = vertex->dz ;
      }
    }
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "%d holes in parameterization filled\n", unfilled) ;

  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    free(filled[u]) ;
    free(distances[u]) ;
  }
  free(filled) ; free(distances) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "done.\n") ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    ImageWrite(mrisp->Ip, "sphere.hipl") ;

  return(mrisp) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SURFACE *
MRISgradientFromParameterization(MRI_SP *mrisp, MRI_SURFACE *mris)
{
  float     a, b, c, phi, theta, x, y, z, uf, vf, du, dv, d ;
  int       vno, u0, v0, u1, v1 ;
  VERTEX    *vertex ;

  if (!mris)
    mris = MRISclone(mrisp->mris) ;

  a = b = c = mris->radius ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(vertex->y/b, vertex->x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = c*c-z*z ; if (d < 0.0) d = 0.0 ;
    phi = atan2(sqrt(d), z) ;
    if (phi < RADIANS(1))
      DiagBreak() ;
    if (vno == 60935)
      DiagBreak() ;

    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u0 = floor(uf) ; u1 = ceil(uf) ;
    v0 = floor(vf) ; v1 = ceil(vf) ;
    du = uf - (float)u0 ; dv = vf - (float)v0 ;

#if 0
    if (vno == 48092)
    {
      DiagBreak() ; 
      fprintf(stderr, "vno %d (%2.1f, %2.1f, %2.1f), u,v = (%2.1f, %2.1f)\n",
             vno, x, y, z, uf, vf) ;
    }
#endif

    /* enforce spherical topology  */
    if (u0 < 0)  /* enforce spherical topology  */
      u0 = -u0 ;
    if (u0 >= U_DIM(mrisp))
      u0 = U_DIM(mrisp)-(u0-U_DIM(mrisp)+1) ;
    if (u1 < 0)  /* enforce spherical topology  */
      u1 = -u1 ;
    if (u1 >= U_DIM(mrisp))
      u1 = U_DIM(mrisp)-(u1-U_DIM(mrisp)+1) ;
    if (v0 < 0) 
      v0 += V_DIM(mrisp) ;
    if (v0 >= V_DIM(mrisp))
      v0 -= V_DIM(mrisp) ;
    if (v1 < 0) 
      v1 += V_DIM(mrisp) ;
    if (v1 >= V_DIM(mrisp))
      v1 -= V_DIM(mrisp) ;

    if (((u0 == DEBUG_U) || (u1 == DEBUG_U)) &&
        ((v0 == DEBUG_V) || (v1 == DEBUG_V)))
      DiagBreak() ;
    
    /* do bilinear interpolation */
    vertex->dx = 
      du        * dv        * *IMAGEFseq_pix(mrisp->Ip, u1, v1, 0) +
      (1.0f-du) * dv        * *IMAGEFseq_pix(mrisp->Ip, u0, v1, 0) +
      (1.0f-du) * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, 0) +
      du        * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, 0) ;


    vertex->dy = 
      du        * dv        * *IMAGEFseq_pix(mrisp->Ip, u1, v1, 1) +
      (1.0f-du) * dv        * *IMAGEFseq_pix(mrisp->Ip, u0, v1, 1) +
      (1.0f-du) * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, 1) +
      du        * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, 1) ;

    vertex->dz = 
      du        * dv        * *IMAGEFseq_pix(mrisp->Ip, u1, v1, 2) +
      (1.0f-du) * dv        * *IMAGEFseq_pix(mrisp->Ip, u0, v1, 2) +
      (1.0f-du) * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, 2) +
      du        * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, 2) ;

  }

  return(mris) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRISPfunctionVal(MRI_SURFACE_PARAMETERIZATION *mrisp, MRI_SURFACE *mris,
                 float x, float y, float z, int fno)
{
  double    r, r2, r4, r6, val, d, x2, y2, z2, g, h, f, dx, dy, dz ;
  float     phi, theta, uf, vf, du, dv ;
  int       u0, v0, u1, v1, u0_voff, u1_voff, u1_v0, u1_v1, u0_v0, u0_v1 ;

  r = sqrt(x*x+y*y+z*z) ;
  if (!FEQUAL(r, mris->radius))  /* project it onto sphere */
  {
    r = mris->radius ;
    r2 = r*r ; r4 = r2*r2 ; r6 = r2*r4 ;
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    f = x2/r6+y2/r6+z2/r6;
    g = 2*(x2/r4+y2/r4+z2/r4);
    h = x2/r2+y2/r2+z2/r2-1;
    d = (-g+(float)sqrt((double)(g*g-4*f*h)))/(2*f);
    dx = d*x/r2;
    dy = d*y/r2;
    dz = d*z/r2;
    x = x+dx ;
    y = y+dy;
    z = z+dz;
  }

  theta = atan2(y/r, x/r) ;
  if (theta < 0.0f)
    theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
  d = r*r-z*z ; if (d < 0.0) d = 0.0 ;
  phi = atan2(sqrt(d), z) ;
  if (phi < RADIANS(1))
    DiagBreak() ;

  uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
  vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
  u0 = floor(uf) ; u1 = ceil(uf) ;
  v0 = floor(vf) ; v1 = ceil(vf) ;
  du = uf - (float)u0 ; dv = vf - (float)v0 ;


  /* enforce spherical topology  */
  if (u0 < 0)  /* enforce spherical topology  */
  {
    u0_voff = V_DIM(mrisp)/2 ;
    u0 = -u0 ;
  }
  else if (u0 >= U_DIM(mrisp))
  {
    u0_voff = V_DIM(mrisp)/2 ;
    u0 = U_DIM(mrisp)-(u0-U_DIM(mrisp)+1) ;
  }
  else
    u0_voff = 0 ;
  if (u1 < 0)  /* enforce spherical topology  */
  {
    u1_voff = V_DIM(mrisp)/2 ;
    u1 = -u1 ;
  }
  else if (u1 >= U_DIM(mrisp))
  {
    u1_voff = V_DIM(mrisp)/2 ;
    u1 = U_DIM(mrisp)-(u1-U_DIM(mrisp)+1) ;
  }
  else
    u1_voff = 0 ;
  if (v0 < 0) 
    v0 += V_DIM(mrisp) ;
  if (v0 >= V_DIM(mrisp))
    v0 -= V_DIM(mrisp) ;
  if (v1 < 0) 
    v1 += V_DIM(mrisp) ;
  if (v1 >= V_DIM(mrisp))
    v1 -= V_DIM(mrisp) ;
  
  u0_v1 = v1 + u0_voff ;
  while (u0_v1 >= V_DIM(mrisp))
    u0_v1 -= V_DIM(mrisp) ;
  u0_v0 = v0 + u0_voff ;
  while (u0_v0 >= V_DIM(mrisp))
    u0_v0 -= V_DIM(mrisp) ;
  u1_v1 = v1 + u1_voff ;
  while (u1_v1 >= V_DIM(mrisp))
    u1_v1 -= V_DIM(mrisp) ;
  u1_v0 = v0 + u1_voff ;
  while (u1_v0 >= V_DIM(mrisp))
    u1_v0 -= V_DIM(mrisp) ;

  /* do bilinear interpolation */
  val = 
    du        * dv        * *IMAGEFseq_pix(mrisp->Ip, u1, u1_v1, fno) +
    (1.0f-du) * dv        * *IMAGEFseq_pix(mrisp->Ip, u0, u0_v1, fno) +
    (1.0f-du) * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u0, u0_v0, fno) +
    du        * (1.0f-dv) * *IMAGEFseq_pix(mrisp->Ip, u1, u1_v0, fno) ;
  return(val) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPclone(MRI_SP *mrisp_src)
{
  MRI_SP   *mrisp_dst ;

  mrisp_dst = (MRI_SP *)calloc(1, sizeof(MRI_SP)) ;
  mrisp_dst->mris = mrisp_src->mris ;
  mrisp_dst->Ip = ImageCopy(mrisp_src->Ip, NULL) ;

  return(mrisp_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_LEN  4

MRI_SP *
MRISPblur(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma, int fno)
{
  int    u, v, cart_klen, klen, khalf, uk, vk, u1, v1, no_sphere, voff, f0, f1;
  double k, total, ktotal, sigma_sq_inv, udiff, vdiff, sin_sq_u, phi ;
  IMAGE  *Ip_src, *Ip_dst ;

  no_sphere = getenv("NO_SPHERE") != NULL ;
  if (no_sphere)
    fprintf(stderr, "disabling spherical geometry\n") ;

  if (!mrisp_dst)
    mrisp_dst = MRISPclone(mrisp_src) ;
  mrisp_dst->sigma = sigma ;

  /* determine the size of the kernel */
  cart_klen = (int)nint(6.0f * sigma)+1 ;
  if (ISEVEN(cart_klen))   /* ensure it's odd */
    cart_klen++ ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "blurring surface, sigma = %2.3f, cartesian klen = %d\n",
            sigma, cart_klen) ;

  if (FZERO(sigma))
    sigma_sq_inv = BIG ;
  else
    sigma_sq_inv = 1.0f / (sigma*sigma) ;

  Ip_src = mrisp_src->Ip ;
  Ip_dst = mrisp_dst->Ip ;
  if (fno < 0)
  { f0 = 0 ; f1 = Ip_src->num_frame-1 ; }
  else
  { f0 = f1 = fno ; }
  for (fno = f0 ; fno <= f1 ; fno++)   /* for each frame */
  {
    for (u = 0 ; u < U_DIM(mrisp_src) ; u++)
    {
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stderr, "\r%3.3d of %d     ", u, U_DIM(mrisp_src)-1) ;
      phi = (double)u*PHI_MAX / PHI_DIM(mrisp_src) ;
      sin_sq_u = sin(phi) ; sin_sq_u *= sin_sq_u ;
      if (!FZERO(sin_sq_u))
      {
        k = cart_klen * cart_klen ;
        klen = sqrt(k + k/sin_sq_u) ;
        if (klen > MAX_LEN*cart_klen)
          klen = MAX_LEN*cart_klen ;
      }
      else
        klen = MAX_LEN*cart_klen ;  /* arbitrary max length */
      if (no_sphere)
        sin_sq_u = 1.0f, klen = cart_klen ;
      if (klen >= U_DIM(mrisp_src))
        klen = U_DIM(mrisp_src)-1 ;
      if (klen >= V_DIM(mrisp_src))
        klen = V_DIM(mrisp_src)-1 ;
      khalf = klen/2 ;
      for (v = 0 ; v < V_DIM(mrisp_src) ; v++)
      {
        /*      theta = (double)v*THETA_MAX / THETA_DIM(mrisp_src) ;*/
        if (u == DEBUG_U && v == DEBUG_V)
          DiagBreak() ;
        
        total = ktotal = 0.0 ;
        for (uk = -khalf ; uk <= khalf ; uk++)
        {
          udiff = (double)(uk*uk) ;  /* distance squared in u */
          
          u1 = u + uk ;
          if (u1 < 0)  /* enforce spherical topology  */
          {
            voff = V_DIM(mrisp_src)/2 ;
            u1 = -u1 ;
          }
          else if (u1 >= U_DIM(mrisp_src))
          {
            u1 = U_DIM(mrisp_src) - (u1-U_DIM(mrisp_src)+1) ;
            voff = V_DIM(mrisp_src)/2 ;
          }
          else
          voff = 0 ;
          
#if 0       
          phi = (double)u1*PHI_MAX / PHI_DIM(mrisp_src) ;
          sin_sq_u = sin(phi) ; sin_sq_u *= sin_sq_u ;
          if (no_sphere)
            sin_sq_u = 1.0f ;
#endif
          
          for (vk = -khalf ; vk <= khalf ; vk++)
          {
            vdiff = (double)(vk*vk) ;
            k = exp(-(udiff+sin_sq_u*vdiff)*sigma_sq_inv) ;
            v1 = v + vk + voff ;
            while (v1 < 0)  /* enforce spherical topology */
              v1 += V_DIM(mrisp_src) ;
            while (v1 >= V_DIM(mrisp_src))
              v1 -= V_DIM(mrisp_src) ;
            ktotal += k ;
            total += k**IMAGEFseq_pix(Ip_src, u1, v1, fno) ;
          }
        }
        if (u == DEBUG_U && v == DEBUG_V)
          DiagBreak() ;
        total /= ktotal ;   /* normalize weights to 1 */
        *IMAGEFseq_pix(mrisp_dst->Ip, u, v, fno) = total ;
      }
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "done.\n") ;

  return(mrisp_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPalloc(float scale, int nfuncs)
{
  MRI_SP   *mrisp ;
  int      u_dim, v_dim ;

  if (FZERO(scale))
    scale = 1.0f ;

  u_dim = nint(scale*64) ; /* for fft */
  v_dim = 2*u_dim ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "allocating %d by %d parameterization\n",u_dim,v_dim) ;
  mrisp = (MRI_SP *)calloc(1, sizeof(MRI_SP)) ;
  if (!mrisp)
    ErrorExit(ERROR_NOMEMORY, "MRISPalloc(%d, %d): allocation failed",
              u_dim, v_dim) ;
  mrisp->Ip = ImageAlloc(v_dim, u_dim, PFFLOAT, nfuncs) ;
  if (!mrisp->Ip)
    ErrorExit(ERROR_NOMEMORY, "MRISPalloc(%d, %d): allocation failed",
              u_dim, v_dim) ;

  return(mrisp) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISPfree(MRI_SP **pmrisp)
{
  MRI_SP  *mrisp ;

  mrisp = *pmrisp ;
  *pmrisp = NULL ;
  ImageFree(&mrisp->Ip) ;
  free(mrisp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPalign(MRI_SP *mrisp_orig, MRI_SP *mrisp_src, MRI_SP *mrisp_tmp, 
           MRI_SP *mrisp_dst)
{
  IMAGE  *Icorr = NULL ;
  int    u_peak, v_peak ;
  float  peak_val, delta_theta, delta_phi ;

  Icorr = ImageCorrelate(mrisp_src->Ip, mrisp_tmp->Ip, 0, NULL) ;
  ImageFindPeak(Icorr, &v_peak, &u_peak, &peak_val) ;
  u_peak -= U_DIM(mrisp_src)/2 ; v_peak -= V_DIM(mrisp_src)/2 ;
  delta_phi = u_peak * PHI_MAX / PHI_DIM(mrisp_src) ;
  delta_theta = v_peak * THETA_MAX / THETA_DIM(mrisp_src) ;
  fprintf(stderr, "peak found at (%d, %d) (%2.2f, %2.2f)\n", 
          u_peak, v_peak, DEGREES(delta_phi), DEGREES(delta_theta)) ;

  mrisp_dst = MRISPtranslate(mrisp_orig, NULL, u_peak, v_peak) ;
  MRISrotate(mrisp_orig->mris, mrisp_orig->mris, delta_phi, delta_theta) ;

#if 0 
  ImageWrite(mrisp_src->Ip, "Isrc.hipl") ;
  ImageWrite(mrisp_tmp->Ip, "Itemplate.hipl") ;
  ImageWrite(mrisp_dst->Ip, "Iorig.hipl") ;
#endif
  ImageWrite(Icorr, "Icorr.hipl") ;
  ImageWrite(mrisp_dst->Ip, "Ialign.hipl") ;

  return(mrisp_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPtranslate(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, int du, int dv)
{
  int    u, v, udst, vdst, voff ;

  if (!mrisp_dst)
    mrisp_dst = MRISPclone(mrisp_src) ;

  /* now form the destination as a translated copy of the source */
  for (u = 0 ; u <= U_MAX_INDEX(mrisp_src) ; u++)
  {
    udst = u + du ;
    if (udst < 0)  /* enforce spherical topology  */
    {
      voff = V_DIM(mrisp_src)/2 ;
      udst = -udst ;
    }
    else if (udst >= U_DIM(mrisp_src))
    {
      voff = V_DIM(mrisp_src)/2 ;
      udst = U_DIM(mrisp_src) - (udst-U_DIM(mrisp_src)+1) ;
    }
    else
      voff = 0 ;
    for (v = 0 ; v <= V_MAX_INDEX(mrisp_src) ; v++)
    {
      vdst = v + dv + voff ;
      while (vdst < 0)  /* enforce spherical topology  */
        vdst += V_DIM(mrisp_src) ;
      while (vdst >= V_DIM(mrisp_src))
        vdst -= V_DIM(mrisp_src) ;
      *IMAGEFpix(mrisp_dst->Ip, udst, vdst) =
        *IMAGEFpix(mrisp_src->Ip, u, v) ;
    }
  }
  return(mrisp_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Combine two surface parameterizations assuming that they
          are arranged in mean/variance pairs with an image of
          degrees of freedom at the end.

          The 1st surface in this case represents a single surface,
          while the second is an average of some (unknown) number.
------------------------------------------------------*/
MRI_SP *
MRISPcombine(MRI_SP *mrisp, MRI_SP *mrisp_template, int fno)
{
  int     u, v, nframes ;
  float   dof, mean, vart, var, val ;
  IMAGE   *Ip, *Ipt ;

  if (!mrisp_template)
    mrisp_template = MRISPclone(mrisp) ;


  Ip = mrisp->Ip ;
  Ipt = mrisp_template->Ip ;
  if ((Ip->rows != Ipt->rows) || (Ip->cols != Ipt->cols))
    ErrorReturn(NULL, 
                 (ERROR_BADPARM, 
                  "MRISPcombine: size mismatch (%d x %d) vs (%d x %d)",
                  Ip->rows, Ip->cols, Ipt->rows, Ipt->cols)) ;

  nframes = mrisp_template->Ip->num_frame ;
  dof = *IMAGEFseq_pix(Ipt, 0, 0, fno+2) ;
  *IMAGEFseq_pix(Ipt, 0, 0, fno+2) += 1.0f ;
  for (u = 0 ; u <= U_MAX_INDEX(mrisp) ; u++)
  {
    for (v = 0 ; v <= V_MAX_INDEX(mrisp) ; v++)
    {
      val = *IMAGEFpix(Ip, u, v) ;
      mean = *IMAGEFseq_pix(Ipt, u, v, fno) ;
      mean = (mean*dof + val) / (dof+1) ;

      vart = *IMAGEFseq_pix(Ipt, u, v, fno+1) ;
      var = (val - mean) ; var *= var ;
      var = (vart*dof + var) / (dof+1) ;
      *IMAGEFseq_pix(Ipt, u, v, fno+1) = var ;

      *IMAGEFseq_pix(Ipt, u, v, fno) = mean ;
    }
  }
  
  return(mrisp_template) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRISPwrite(MRI_SP *mrisp, char *fname)
{
  ImageWrite(mrisp->Ip, fname) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPread(char *fname)
{
  MRI_SP  *mrisp ;

  mrisp = (MRI_SP *)calloc(1, sizeof(MRI_SP)) ;
  if (!mrisp)
    ErrorExit(ERROR_NOMEMORY, "MRISPread(%s): allocation failed",fname) ;
  mrisp->Ip = ImageRead(fname) ;
  if (!mrisp)
    ErrorReturn(NULL, (ERROR_NOFILE,
                       "MRISPread(%s): could not open file",fname)) ;
  
  return(mrisp) ;
}


