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

#define DEBUG_VNO 79221
#define DEBUG_U  -10  /* -10 25*/
#define DEBUG_V  0  /* 2 */

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/

MRI_SP *
MRIStoParameterization(MRI_SURFACE *mris, MRI_SP *mrisp, float scale)
{
  float     a, b, c, phi, theta, x, y, z, total, uf, vf, d, du, dv, total_d,
            **distances, sigma, two_sigma_sq ;
  int       vno, u, v, uk, vk, n, unfilled, u0, v0, u1, v1, **filled ;
  VERTEX    *vertex ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "computing parameterization...") ;

  if (!mrisp)
    mrisp = MRISPalloc(mris, scale) ;

#if 0
  if (FZERO(mris->a))
  {
    a = DEFAULT_A ;
    b = DEFAULT_B ;
    c = DEFAULT_C ;
  }
  else
  {
    a = mris->a ;
    b = mris->b ;
    c = mris->c ;
  }
#else
  vertex = &mris->vertices[0] ;
  x = vertex->x ; y = vertex->y ; z = vertex->z ;
  a = b = c = sqrt(x*x + y*y + z*z) ;
#endif

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
    phi = atan2(sqrt(c*c-z*z), z) ;
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
    if (u0 < 0)  /* enforce spherical topology  */
      u0 += U_DIM(mrisp) ;
    if (u0 >= U_DIM(mrisp))
      u0 -= U_DIM(mrisp) ;
    if (u1 < 0)  /* enforce spherical topology  */
      u1 += U_DIM(mrisp) ;
    if (u1 >= U_DIM(mrisp))
      u1 -= U_DIM(mrisp) ;
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

#if 0
    d = sqrt(du*du+dv*dv) ;  /* 1-distance to vertex 0,0 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+dv*dv ;  /* 1-distance to vertex 0,0 */
    d = exp(-d/two_sigma_sq) ;
#endif
    if (u0 >= U_DIM(mrisp) || v0 >= V_DIM(mrisp))
      DiagBreak() ;
    filled[u0][v0] = vno ;
    distances[u0][v0] += d ;
    if ((u0 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;

    d = 1.0 - dv ;         /* distance to v1 */
#if 0
    d = sqrt(du*du+d*d) ;  /* distance to vertex 0,1 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+d*d ;  /* distance to vertex 0,1 */
    d = exp(-d/two_sigma_sq) ;
#endif
    distances[u0][v1] += d ;         /* keep track of total distance */
    filled[u0][v1] = vno ;
    if ((u0 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;

    d = 1.0 - du ;         /* distance to u1 */
#if 0
    d = sqrt(d*d+dv*dv) ;  /* distance to vertex 1,0 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = d*d+dv*dv ;  /* distance to vertex 1,0 */
    d = exp(-d/two_sigma_sq) ;
#endif

    distances[u1][v0] += d ;         /* keep track of total distance */
    filled[u1][v0] = vno ;
    if ((u1 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n", vno, x, y, z, d, vertex->curv) ;

    du = 1.0 - du ; dv = 1.0 - dv ;
#if 0
    d = sqrt(du*du+dv*dv) ;  /* 1-distance to vertex 1,1 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+dv*dv ;  /* 1-distance to vertex 1,1 */
    d = exp(-d/two_sigma_sq) ;
#endif
    distances[u1][v1] += d ;         /* keep track of total distance */
    filled[u1][v1] = vno ;
    if ((u1 == DEBUG_U) && (v1 == DEBUG_V))
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
    phi = atan2(sqrt(c*c-z*z), z) ;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX ;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX ;
    u0 = floor(uf) ; u1 = ceil(uf+0.00001f) ;
    v0 = floor(vf) ; v1 = ceil(vf+0.00001f) ;
    du = uf - (float)u0 ; dv = vf - (float)v0 ;
    if ((u0 == DEBUG_U || u1 == DEBUG_U) && (v0 == DEBUG_V || v1 == DEBUG_V))
      DiagBreak() ;
    if (u0 < 0)  /* enforce spherical topology  */
      u0 += U_DIM(mrisp) ;
    if (u0 >= U_DIM(mrisp))
      u0 -= U_DIM(mrisp) ;
    if (u1 < 0)  /* enforce spherical topology  */
      u1 += U_DIM(mrisp) ;
    if (u1 >= U_DIM(mrisp))
      u1 -= U_DIM(mrisp) ;
    if (v0 < 0)  /* enforce spherical topology  */
      v0 += V_DIM(mrisp) ;
    if (v0 >= V_DIM(mrisp))
      v0 -= V_DIM(mrisp) ;
    if (v1 < 0)  /* enforce spherical topology  */
      v1 += V_DIM(mrisp) ;
    if (v1 >= V_DIM(mrisp))
      v1 -= V_DIM(mrisp) ;

    /* 0,0 */
#if 0
    d = sqrt(du*du+dv*dv) ;   /* distance to 0,0 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+dv*dv ;   /* distance to 0,0 */
    d = exp(-d/two_sigma_sq) ;
#endif
    total_d = distances[u0][v0] ;
    d /= total_d ;
    *IMAGEFpix(mrisp->Ip, u0,v0) += d*vertex->curv ;  /* weight by distance */
    if ((u0 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

    /* 1,0 */
    d = 1.0 - du ;          /* distance to u1 */
#if 0
    d = sqrt(d*d+dv*dv) ;   /* distance to u1,v0 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = d*d+dv*dv ;   /* distance to u1,v0 */
    d = exp(-d/two_sigma_sq) ;
#endif
    total_d = distances[u1][v0] ;
    d = d / total_d ;
    *IMAGEFpix(mrisp->Ip, u1, v0) += d*vertex->curv ;  /* weight by distance */
    if ((u1 == DEBUG_U) && (v0 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

    /* 0,1 */
    d = 1.0 - dv ; 
#if 0
    d = sqrt(du*du+d*d) ;   /* distance to u0,v1 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = du*du+d*d ;   /* distance to u0,v1 */
    d = exp(-d/two_sigma_sq) ;
#endif
    total_d = distances[u0][v1] ;
    d = d / total_d ;
    *IMAGEFpix(mrisp->Ip, u0, v1) += d*vertex->curv ;  /* weight by distance */
    if ((u0 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

    /* 1, 1 */
    du = 1.0 - du ; dv = 1.0 - dv ;
#if 0
    d = sqrt(du*du+dv*dv) ;   /* distance to 1,1 */
    if (FZERO(d))
      d = BIG ;
    else
      d = 1.0f / d ;
#else
    d = (du*du+dv*dv) ;   /* distance to 1,1 */
    d = exp(-d/two_sigma_sq) ;
#endif
    total_d = distances[u1][v1] ;
    d = d / total_d ;
    *IMAGEFpix(mrisp->Ip, u1, v1) += d*vertex->curv ;  /* weight by distance */
    if ((u1 == DEBUG_U) && (v1 == DEBUG_V))
      fprintf(stderr, "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n", vno, x, y, z, vertex->curv, d) ;

  }

  if (DEBUG_U >= 0)
    fprintf(stderr,"curv[%d][%d] = %2.3f\n\n", DEBUG_U, DEBUG_V, 
            *IMAGEFpix(mrisp->Ip, DEBUG_U, DEBUG_V));

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
              u1 += U_DIM(mrisp) ;
            else if (u1 >= U_DIM(mrisp))
              u1 -= U_DIM(mrisp) ;
            for (vk = -1 ; vk <= 1 ; vk++)
            {
              v1 = v + vk ;
              if (v1 < 0)  /* enforce spherical topology  */
                v1 += V_DIM(mrisp) ;
              else if (v1 >= V_DIM(mrisp))
                v1 -= V_DIM(mrisp) ;
              if (filled[u1][v1] != UNFILLED_ELT)
              {
                total += *IMAGEFpix(mrisp->Ip, u1, v1) ;
                n++ ;
              }
            }
          }
          if (n)
          {
            filled[u][v] = FILLED_ELT ;
            *IMAGEFpix(mrisp->Ip, u, v) = total / (float)n ;
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

  if (Gdiag & DIAG_SHOW)
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
MRISfromParameterization(MRI_SP *mrisp, MRI_SURFACE *mris)
{
  float     a, b, c, phi, theta, x, y, z, uf, vf, du, dv, curv ;
  int       vno, u0, v0, u1, v1 ;
  VERTEX    *vertex ;

  if (!mris)
    mris = MRISclone(mrisp->mris) ;

#if 0
  if (FZERO(mris->a))
  {
    a = DEFAULT_A ;
    b = DEFAULT_B ;
    c = DEFAULT_C ;
  }
  else
  {
    a = mris->a ;
    b = mris->b ;
    c = mris->c ;
  }
#else
  vertex = &mris->vertices[0] ;
  x = vertex->x ; y = vertex->y ; z = vertex->z ;
  a = b = c = sqrt(x*x + y*y + z*z) ;
#endif

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    x = vertex->x ; y = vertex->y ; z = vertex->z ;
    theta = atan2(vertex->y/b, vertex->x/a) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    phi = atan2(sqrt(c*c-z*z), z) ;
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
    if (u0 < 0)
      u0 += U_DIM(mrisp) ;
    if (u0 >= U_DIM(mrisp))
      u0 -= U_DIM(mrisp) ;
    if (u1 >= U_DIM(mrisp))
      u1 -= U_DIM(mrisp) ;
    if (v0 < 0) 
      v0 += V_DIM(mrisp) ;
    if (v0 >= V_DIM(mrisp))
      v0 -= V_DIM(mrisp) ;
    if (v1 < 0) 
      v1 += V_DIM(mrisp) ;
    if (v1 >= V_DIM(mrisp))
      v1 -= V_DIM(mrisp) ;
    
    
    /* do bilinear interpolation */
    curv = 
      du        * dv        * *IMAGEFpix(mrisp->Ip, u1, v1) +
      (1.0f-du) * dv        * *IMAGEFpix(mrisp->Ip, u0, v1) +
      (1.0f-du) * (1.0f-dv) * *IMAGEFpix(mrisp->Ip, u0, v0) +
      du        * (1.0f-dv) * *IMAGEFpix(mrisp->Ip, u1, v0) ;

    vertex->curv = curv ;
  }

  return(mris) ;
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
  mrisp_dst->Ip = ImageClone(mrisp_src->Ip) ;

  return(mrisp_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_LEN  4

MRI_SP *
MRISPblur(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma)
{
  int         u, v, cart_klen, klen, khalf, uk, vk, u1, v1, no_sphere ;
  double      k, a, b, c, total, ktotal, sigma_sq_inv, udiff, vdiff, sin_sq_u, 
              phi ;
  MRI_SURFACE *mris ;
  VERTEX      *vertex ;
  float       x, y, z ;

  no_sphere = getenv("NO_SPHERE") != NULL ;
  if (no_sphere)
    fprintf(stderr, "disabling spherical geometry\n") ;

#if 0
  IMAGE       *Iblur ;

  Iblur = ImageGaussian1d(sigma, 0) ;
  if (!mrisp_dst)
    mrisp_dst = MRISPclone(mrisp_src) ;
  ImageCircularConvolveGaussian(mrisp_src->Ip, Iblur, mrisp_dst->Ip, 0) ;

#else
  if (!mrisp_dst)
    mrisp_dst = MRISPclone(mrisp_src) ;

  /* determine the size of the kernel */
  cart_klen = (int)nint(6.0f * sigma)+1 ;
  if (ISEVEN(cart_klen))   /* ensure it's odd */
    cart_klen++ ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "blurring surface, sigma = %2.3f, cartesian klen = %d\n",
            sigma, cart_klen) ;

  mris = mrisp_src->mris ;
#if 0
  if (FZERO(mris->a))
  {
    a = DEFAULT_A ;
    b = DEFAULT_B ;
    c = DEFAULT_C ;
  }
  else
  {
    a = mris->a ;
    b = mris->b ;
    c = mris->c ;
  }
#else
  vertex = &mris->vertices[0] ;
  x = vertex->x ; y = vertex->y ; z = vertex->z ;
  a = b = c = sqrt(x*x + y*y + z*z) ;
#endif

  if (FZERO(sigma))
    sigma_sq_inv = BIG ;
  else
    sigma_sq_inv = 1.0f / (sigma*sigma) ;
  
  for (u = 0 ; u < U_DIM(mrisp_src) ; u++)
  {
    fprintf(stderr, "\r%3.3d of %d     ", u, U_DIM(mrisp_src)) ;
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
          u1 += U_DIM(mrisp_src) ;
        else if (u1 >= U_DIM(mrisp_src))
          u1 -= U_DIM(mrisp_src) ;

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
          v1 = v + vk ;
          if (v1 < 0)  /* enforce spherical topology */
            v1 += V_DIM(mrisp_src) ;
          else if (v1 >= V_DIM(mrisp_src))
            v1 -= V_DIM(mrisp_src) ;
          ktotal += k ;
          total += k**IMAGEFpix(mrisp_src->Ip, u1, v1) ;
        }
      }
      if (u == DEBUG_U && v == DEBUG_V)
        DiagBreak() ;
      total /= ktotal ;   /* normalize weights to 1 */
      *IMAGEFpix(mrisp_dst->Ip, u, v) = total ;
    }
  }
#endif

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;

  return(mrisp_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *
MRISPalloc(MRI_SURFACE *mris, float scale)
{
  MRI_SP   *mrisp ;
  int      u_dim, v_dim ;

  u_dim = nint(sqrt((float)mris->nvertices/(2*scale))) / 3 ;
  v_dim = 2*u_dim ;
  u_dim = 64 ; v_dim = 128 ; /* for fft */
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "allocating %d by %d parameterization\n",u_dim,v_dim) ;
  mrisp = (MRI_SP *)calloc(1, sizeof(MRI_SP)) ;
  if (!mrisp)
    ErrorExit(ERROR_NOMEMORY, "MRISPalloc(%d, %d): allocation failed",
              u_dim, v_dim) ;
  mrisp->mris = mris ;
  mrisp->Ip = ImageAlloc(v_dim, u_dim, PFFLOAT, 1) ;
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
  int    u, v, udst, vdst ;

  if (!mrisp_dst)
    mrisp_dst = MRISPclone(mrisp_src) ;

  /* now form the destination as a translated copy of the source */
  for (u = 0 ; u <= U_MAX_INDEX(mrisp_src) ; u++)
  {
    udst = u + du ;
    if (udst < 0)  /* enforce spherical topology  */
      udst += U_DIM(mrisp_src) ;
    if (udst >= U_DIM(mrisp_src))
      udst -= U_DIM(mrisp_src) ;
    for (v = 0 ; v <= V_MAX_INDEX(mrisp_src) ; v++)
    {
      vdst = v + dv ;
      if (vdst < 0)  /* enforce spherical topology  */
        vdst += V_DIM(mrisp_src) ;
      if (vdst >= V_DIM(mrisp_src))
        vdst -= V_DIM(mrisp_src) ;
      *IMAGEFpix(mrisp_dst->Ip, udst, vdst) =
        *IMAGEFpix(mrisp_src->Ip, u, v) ;
    }
  }
  return(mrisp_dst) ;
}
