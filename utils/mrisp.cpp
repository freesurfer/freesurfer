/**
 * @brief utilities for spherical parameterization of an MRIS
 *
 * mrisp = MRIS parameterization contains utilities for writing various
 * fields over the surface (e.g. curvature, coordinate functions) into a
 * spherical (longitude/colatitude) parameterization in the form of a 2D
 * image. Also for Gaussian blurring and other utilities.
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <math.h>
#include <stdio.h>

#include "diag.h"
#include "error.h"
#include "macros.h"
#include "mrisurf.h"
#include "proto.h"
#include "utils.h"
#include "mrishash.h"
#include "mri_identify.h"
#include "mrisurf_sphere_interp.h"
#include "romp_support.h"
#include "mrisp.h"

/*---------------------------- STRUCTURES -------------------------*/

/*---------------------------- CONSTANTS -------------------------*/

#define BIG 100000.0

/* something of a hack... */
#define UNFILLED_ELT -2
#define FILLING_ELT -1
#define FILLED_ELT 1

#define DEFAULT_UDIM 256

#define DEBUG_VNO -1
 int DEBUG_U = -1;
 int DEBUG_V = -1;

static int spherical_coordinate(double x, double y, double z, double *pphi, double *ptheta);


// even though all spheres should be centered at the origin, it's possible they aren't, and
// this needs to be accounted for when creating parameterizations
static MRIS* makeCenteredSphere(MRIS *surf)
{
  mrisComputeSurfaceDimensions(surf);
  float distance = std::sqrt((surf->xctr * surf->xctr) + (surf->yctr * surf->yctr) + (surf->zctr * surf->zctr));

  if (distance > 0.1) {
    fs::warning() << "Sphere is not centered (distance = " << distance << "). Moving center to origin - this might be unstable.";
    surf = MRISclone(surf);
    MRISmoveOrigin(surf, 0, 0, 0);
  }

  return surf;
}


// if 'centered' was copied from the original, copy the new curv values and free the centered surface
static void resetCenteredSphere(MRIS *orig, MRIS *centered)
{
  if (orig != centered) {
    for (int vno = 0 ; vno < orig->nvertices ; vno++) orig->vertices[vno].curv = centered->vertices[vno].curv;
    MRISfree(&centered);
  }
}


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *MRIStoParameterization(MRIS *mris, MRI_SP *mrisp, float scale, int fno)
{
  MRIS *original = mris;
  mris = makeCenteredSphere(mris);

  float a, b, c, phi, theta, x, y, z, uf, vf, d, total_d, **distances, *fp;
  int vno, u, v, unfilled, **filled, npasses, nfilled;
  VERTEX *vertex;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "computing parameterization...");

  if (!mrisp)
    mrisp = MRISPalloc(scale, 1);
  else
    ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, fno);

  a = b = c = MRISaverageRadius(mris);

  filled = (int **)calloc(U_DIM(mrisp), sizeof(int *));
  distances = (float **)calloc(U_DIM(mrisp), sizeof(float *));
  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    filled[u] = (int *)calloc(V_DIM(mrisp), sizeof(int));
    distances[u] = (float *)calloc(V_DIM(mrisp), sizeof(float));

    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) filled[u][v] = UNFILLED_ELT;
  }

  fp = IMAGEFseq_pix(mrisp->Ip, DEBUG_U, DEBUG_V, fno);
  /* first calculate total distances to a point in parameter space */
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    if (vno == Gdiag_no) DiagBreak();
    theta = atan2(y / b, x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (vno == DEBUG_VNO) DiagBreak();
    vertex->phi = phi;
    vertex->theta = theta;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u = nint(uf);
    v = nint(vf);
    if (u < 0) /* enforce spherical topology  */
      u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) /* enforce spherical topology  */
      v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    if (u == 0 && v == 56) DiagBreak();
    if ((((u == DEBUG_U) && (v == DEBUG_V)) || (vno == Gdiag_no)) && (fno==0))
    {
      printf("v %d --> [%d, %d] (%2.1f, %2.1f)\n", vno, u, v, theta, phi);
      DiagBreak();
    }

    filled[u][v] = vno;
    distances[u][v] += 1; /* keep track of total # of nodes */
    if ((u == DEBUG_U) && (v == DEBUG_V))
      fprintf(stderr,
              "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n",
              vno,
              x,
              y,
              z,
              d,
              vertex->curv);
  }

  if (DEBUG_U >= 0) fprintf(stderr, "\ndistance[%d][%d] = %2.3f\n\n", DEBUG_U, DEBUG_V, distances[DEBUG_U][DEBUG_V]);

  /* now add in curvatures proportional to their distance from the point */
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(y / b, x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u = nint(uf);
    v = nint(vf);
    if (u < 0) /* enforce spherical topology  */
      u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) /* enforce spherical topology  */
      v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    if (u == 0 && v == 56) DiagBreak();

    /* 0,0 */
    total_d = distances[u][v];
    if ((total_d > 10000.0) || (vertex->curv > 1000.0)) DiagBreak();
    if (total_d > 0.0) 
      *IMAGEFseq_pix(mrisp->Ip, u, v, fno) += vertex->curv / total_d;
    if (devFinite(*IMAGEFseq_pix(mrisp->Ip, u, v, fno)) == 0)
      DiagBreak() ;
    if ((u == DEBUG_U) && (v == DEBUG_V))
      fprintf(stderr,
              "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n",
              vno,
              x,
              y,
              z,
              vertex->curv,
              d);
  }

  if (DEBUG_U >= 0)
    fprintf(stderr, "curv[%d][%d] = %2.3f\n\n", DEBUG_U, DEBUG_V, *IMAGEFseq_pix(mrisp->Ip, DEBUG_U, DEBUG_V, fno));

  /* fill in values which were unmapped using soap bubble */
  nfilled = npasses = 0;
  do {
    IMAGE *Ip, *Itmp;
    int u1, v1, uk, vk, n;
    float total;

    Ip = mrisp->Ip;
    Itmp = ImageClone(Ip);
    ImageCopyFrames(Ip, Itmp, 0, Ip->num_frame, 0);
    unfilled = 0;
    ;
    for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
      for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
        if ((u == DEBUG_U) && (v == DEBUG_V)) 
	  DiagBreak();
	if (devFinite(*IMAGEFseq_pix(mrisp->Ip, u, v, fno)) == 0)
	  DiagBreak() ;
        if (filled[u][v] == UNFILLED_ELT) {
          for (total = 0.0f, n = 0, uk = -1; uk <= 1; uk++) {
            u1 = u + uk;
            if (u1 < 0) /* enforce spherical topology  */
              u1 = -u1;
            else if (u1 >= U_DIM(mrisp))
              u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
            for (vk = -1; vk <= 1; vk++) {
              v1 = v + vk;
              if (v1 < 0) /* enforce spherical topology  */
                v1 += V_DIM(mrisp);
              else if (v1 >= V_DIM(mrisp))
                v1 -= V_DIM(mrisp);

              if (filled[u1][v1] >= 0) {
		if (devFinite(*IMAGEFseq_pix(mrisp->Ip, u1, v1, fno)) == 0)
		  DiagBreak() ;
                total += *IMAGEFseq_pix(Ip, u1, v1, fno);
                n++;
              }
            }
          }
          if (n > 0) {
            total /= (float)n;
	    if (devFinite(*IMAGEFseq_pix(Itmp, u, v, fno)) == 0)
	      DiagBreak() ;
            *IMAGEFseq_pix(Itmp, u, v, fno) = total;
            filled[u][v] = FILLING_ELT;
            nfilled++;
          }
          else
            unfilled++;
        }
        else
          *IMAGEFseq_pix(Itmp, u, v, fno) = *IMAGEFseq_pix(Ip, u, v, fno);
      }
    }
    for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
      for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
        if (filled[u][v] == FILLING_ELT) filled[u][v] = FILLED_ELT;
      }
    }
    mrisp->Ip = Itmp;
    ImageFree(&Ip);
    if (npasses++ > 1000)
      ErrorExit(ERROR_BADFILE,
                "MRISPtoParameterization: could not fill "
                "parameterization");
  } while (unfilled > 0);
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) fprintf(stderr, "filling %d elements took %d passes\n", nfilled, npasses);

  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    free(filled[u]);
    free(distances[u]);
  }
  free(filled);
  free(distances);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "done.\n");

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) ImageWrite(mrisp->Ip, "sphere.hipl");

  resetCenteredSphere(original, mris);

  return (mrisp);
}


MRI_SP *MRIStoParameterizationBarycentric(MRIS *mris, MRI_SP *mrisp, float scale, int frameno)
{
  MRIS *original = mris;
  mris = makeCenteredSphere(mris);

  if (!mrisp) {
    mrisp = MRISPalloc(scale, 1);
  } else {
    ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, frameno);
  }

  // allocate fill-marker and distance arrays
  int **filled = (int **)calloc(U_DIM(mrisp), sizeof(int *));
  float **distances = (float **)calloc(U_DIM(mrisp), sizeof(float *));
  for (int u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    filled[u] = (int *)calloc(V_DIM(mrisp), sizeof(int));
    distances[u] = (float *)calloc(V_DIM(mrisp), sizeof(float));
    // mark all as unfilled
    for (int v = 0; v <= V_MAX_INDEX(mrisp); v++) filled[u][v] = UNFILLED_ELT;
  }

  float radius = MRISaverageRadius(mris);

  // first calculate total distances to a point in parameter space
  for (int vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *vertex = &mris->vertices[vno];
    float x = vertex->x;
    float y = vertex->y;
    float z = vertex->z;

    // translate xyz to spherical coordinates
    float theta = atan2(y / radius, x / radius);
    if (theta < 0.0f) theta = 2 * M_PI + theta;  // make it 0 -> 2PI
    float d = radius * radius - z * z;
    if (d < 0.0) d = 0;
    float phi = atan2(sqrt(d), z);

    // cache in vertex for next loop
    vertex->phi = phi;
    vertex->theta = theta;

    // translate to image coordinates
    float uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    float vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    int u = nint(uf);
    int v = nint(vf);

    // enforce spherical topology
    if (u < 0) u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    // keep track of total # of nodes
    filled[u][v] = vno;
    distances[u][v] += 1;
  }

  // now add in curvatures proportional to their distance from the point
  for (int vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *vertex = &mris->vertices[vno];
    float phi = vertex->phi;
    float theta = vertex->theta;

    float uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    float vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    int u = nint(uf);
    int v = nint(vf);

    if (u < 0) u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    float total_dist = distances[u][v];
    if (total_dist > 0.0) *IMAGEFseq_pix(mrisp->Ip, u, v, frameno) += vertex->curv / total_dist;
  }

  // do backwards sampling to fill in missing pixels
  SphericalInterpolator interpolator = SphericalInterpolator(mris);
  for (int u = 0; u <= U_MAX_INDEX(mrisp); u++)  {
    for (int v = 0; v <= V_MAX_INDEX(mrisp); v++) {
      if (filled[u][v] == UNFILLED_ELT) {
        double phi = u * PHI_MAX / PHI_DIM(mrisp);
        double theta = v * THETA_MAX / THETA_DIM(mrisp);
        *IMAGEFseq_pix(mrisp->Ip, u, v, frameno) = interpolator.interp(phi, theta);
      }
    }
  }

  // free the fill-marker and distance arrays
  for (int u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    free(filled[u]);
    free(distances[u]);
  }
  free(filled);
  free(distances);

  resetCenteredSphere(original, mris);

  return mrisp;
}


MRIS *MRISfromParameterizationBarycentric(MRI_SP *mrisp, MRIS *mris, int fno)
{
  float a, b, c, phi, theta, x, y, z, uf, vf, du, dv, curv, d;
  int vno, u0, v0, u1, v1;
  VERTEX *vertex;

  if (!mris) mris = MRISclone(mrisp->mris);

  MRIS *original = mris;
  mris = makeCenteredSphere(mris);

  a = b = c = MRISaverageRadius(mris);

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(vertex->y / b, vertex->x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (phi > M_PI) DiagBreak();
    if (vno == 60935) DiagBreak();

    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u0 = floor(uf);
    u1 = ceil(uf);
    v0 = floor(vf);
    v1 = ceil(vf);
    du = uf - (float)u0;
    dv = vf - (float)v0;

    /* enforce spherical topology  */
    if (u0 < 0) /* enforce spherical topology  */
      u0 = -u0;
    if (u0 >= U_DIM(mrisp)) u0 = U_DIM(mrisp) - (u0 - U_DIM(mrisp) + 1);
    if (u1 < 0) /* enforce spherical topology  */
      u1 = -u1;
    if (u1 >= U_DIM(mrisp)) u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
    if (v0 < 0) v0 += V_DIM(mrisp);
    if (v0 >= V_DIM(mrisp)) v0 -= V_DIM(mrisp);
    if (v1 < 0) v1 += V_DIM(mrisp);
    if (v1 >= V_DIM(mrisp)) v1 -= V_DIM(mrisp);

    if (((u0 == DEBUG_U) || (u1 == DEBUG_U)) && ((v0 == DEBUG_V) || (v1 == DEBUG_V))) DiagBreak();
    /* do bilinear interpolation */
    curv = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) +
           (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, fno) +
           (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, fno) +
           du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, fno);

    vertex->curv = curv;
  }

  resetCenteredSphere(original, mris);
  return original;
}

int
MRIScoordsFromParameterizationBarycentric(MRIS *mris, MRI_SP *mrisp, int which_vertices)
{
  float    *curvs ;
  int      vno ;

  if (which_vertices != WHITE_VERTICES)
    ErrorExit(ERROR_UNSUPPORTED, "MRIScoordsToParameterizationBarycentric: unsupported which_vertices %d", which_vertices) ;

  curvs = (float *)calloc(mris->nvertices, sizeof(float));
  if (curvs == NULL)
    ErrorExit(ERROR_NOMEMORY, "MRIScoordsToParameterizationBarycentric: could not allocate %d-len curvature vector", mris->nvertices) ;

  MRISextractCurvatureVector(mris, curvs) ;
  MRISfromParameterizationBarycentric(mrisp, mris, 0) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].whitex = mris->vertices[vno].curv ;
  MRISfromParameterization(mrisp, mris, 1) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].whitey = mris->vertices[vno].curv ;
  MRISfromParameterization(mrisp, mris, 2) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].whitez = mris->vertices[vno].curv ;
  MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
  MRISimportCurvatureVector(mris, curvs) ;
  free(curvs) ;
  return(NO_ERROR) ;
}


static double UF = 254.8, VF = 409.5;

MRI_SP *MRIScoordsToParameterization(MRIS *mris, MRI_SP *mrisp, float scale, int which_vertices)
{
  float a, b, c, phi, theta, x, y, z, uf, vf, d, total_d, **distances;
  int vno, u, v, unfilled, **filled, npasses, nfilled;
  VERTEX *vertex;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "computing parameterization...");

  if (!mrisp)
    mrisp = MRISPalloc(scale, 3);
  else {
    ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, 0);
    ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, 1);
    ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, 2);
  }

  mrisp->radius = a = b = c = MRISaverageRadius(mris);

  filled = (int **)calloc(U_DIM(mrisp), sizeof(int *));
  distances = (float **)calloc(U_DIM(mrisp), sizeof(float *));
  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    filled[u] = (int *)calloc(V_DIM(mrisp), sizeof(int));
    distances[u] = (float *)calloc(V_DIM(mrisp), sizeof(float));

    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) filled[u][v] = UNFILLED_ELT;
  }

  /* first calculate total distances to a point in parameter space */
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    if (vno == Gdiag_no) DiagBreak();
    theta = atan2(vertex->y / b, vertex->x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (phi > M_PI) DiagBreak();
    if (vno == DEBUG_VNO) DiagBreak();
    vertex->phi = phi;
    vertex->theta = theta;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    if (fabs(uf - UF) < 1 && fabs(vf - VF) < 1) DiagBreak();
    u = nint(uf);
    v = nint(vf);
    if (u < 0) /* enforce spherical topology  */
      u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) /* enforce spherical topology  */
      v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    if (vno == Gdiag_no) {
      DEBUG_U = u;
      DEBUG_V = v;
      DiagBreak();
    }
    if ((u == DEBUG_U) && (v == DEBUG_V)) 
    {
      printf("v %d --> [%d, %d] (%2.1f, %2.1f)\n", vno, u, v, theta, phi);
      DiagBreak();
    }

    filled[u][v] = vno;
    distances[u][v] += 1; /* keep track of total # of nodes */
    if ((((u == DEBUG_U) && (v == DEBUG_V)) || (Gdiag_no == vno)) && DIAG_VERBOSE_ON)
      fprintf(stderr,
              "v = %d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "x = (%2.3f,%2.3f,%2.3f)\n",
              vno,
              x,
              y,
              z,
              d,
              vertex->origx,
              vertex->origy,
              vertex->origz);
  }

  if (DEBUG_U >= 0) fprintf(stderr, "\ndistance[%d][%d] = %2.3f\n\n", DEBUG_U, DEBUG_V, distances[DEBUG_U][DEBUG_V]);

  /* now add in coords proportional to their distance from the point */
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(y / b, x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    if (phi > M_PI) DiagBreak();
    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u = nint(uf);
    v = nint(vf);
    if (u < 0) /* enforce spherical topology  */
      u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) /* enforce spherical topology  */
      v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    if (u == 0 && v == 56) DiagBreak();

    /* 0,0 */
    total_d = distances[u][v];
    if ((total_d > 10000.0) || (vertex->curv > 1000.0)) DiagBreak();
    switch (which_vertices) {
      case ORIGINAL_VERTICES:
        x = vertex->origx;
        y = vertex->origy;
        z = vertex->origz;
        break;
      case CURRENT_VERTICES:
        x = vertex->x;
        y = vertex->y;
        z = vertex->z;
        break;
      case PIAL_VERTICES:
        x = vertex->pialx;
        y = vertex->pialy;
        z = vertex->pialz;
        break;
      case WHITE_VERTICES:
        x = vertex->whitex;
        y = vertex->whitey;
        z = vertex->whitez;
        break;
    default:
      ErrorExit(ERROR_UNSUPPORTED, "MRIScoordsToParameterization: unsupported vertex set %d", which_vertices) ;
    }
    if (total_d > 0.0) {
      *IMAGEFseq_pix(mrisp->Ip, u, v, 0) += x / total_d;
      *IMAGEFseq_pix(mrisp->Ip, u, v, 1) += y / total_d;
      *IMAGEFseq_pix(mrisp->Ip, u, v, 2) += z / total_d;
    }
    if ((u == DEBUG_U) && (v == DEBUG_V)) fprintf(stderr, "v = %6.6d, x = (%2.3f,%2.3f,%2.3f)\n", vno, x, y, z);
  }

  if ((DEBUG_U >= 0 || (vno == Gdiag_no)))
    fprintf(stderr,
            "x[%d][%d] = (%2.3f,%2.3f,%2.3f)\n\n",
            DEBUG_U,
            DEBUG_V,
            *IMAGEFseq_pix(mrisp->Ip, DEBUG_U, DEBUG_V, 0),
            *IMAGEFseq_pix(mrisp->Ip, DEBUG_U, DEBUG_V, 1),
            *IMAGEFseq_pix(mrisp->Ip, DEBUG_U, DEBUG_V, 2));

  /* fill in values which were unmapped using soap bubble */
  nfilled = npasses = 0;
  do {
    IMAGE *Ip, *Itmp;
    int u1, v1, uk, vk, n;
    float totalx, totaly, totalz;

    Ip = mrisp->Ip;
    Itmp = ImageClone(Ip);
    ImageCopyFrames(Ip, Itmp, 0, Ip->num_frame, 0);
    unfilled = 0;
    ;
    for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
      for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
        if ((u == DEBUG_U) && (v == DEBUG_V)) DiagBreak();
        if (filled[u][v] == UNFILLED_ELT) {
          totalx = totaly = totalz = 0.0f;
          for (n = 0, uk = -1; uk <= 1; uk++) {
            u1 = u + uk;
            if (u1 < 0) /* enforce spherical topology  */
              u1 = -u1;
            else if (u1 >= U_DIM(mrisp))
              u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
            for (vk = -1; vk <= 1; vk++) {
              v1 = v + vk;
              if (v1 < 0) /* enforce spherical topology  */
                v1 += V_DIM(mrisp);
              else if (v1 >= V_DIM(mrisp))
                v1 -= V_DIM(mrisp);

              if (filled[u1][v1] >= 0) {
                float x, y, z;

                x = *IMAGEFseq_pix(Ip, u1, v1, 0);
                y = *IMAGEFseq_pix(Ip, u1, v1, 1);
                z = *IMAGEFseq_pix(Ip, u1, v1, 2);
                if (u == 255 && v == 410) DiagBreak();
                totalx += x;
                totaly += y;
                totalz += z;
                n++;
              }
            }
          }
          if (n > 0) {
            if (u == 1 && v == 3) DiagBreak();
            totalx /= (float)n;
            totaly /= (float)n;
            totalz /= (float)n;
            *IMAGEFseq_pix(Itmp, u, v, 0) = totalx;
            *IMAGEFseq_pix(Itmp, u, v, 1) = totaly;
            *IMAGEFseq_pix(Itmp, u, v, 2) = totalz;
            filled[u][v] = FILLING_ELT;
            nfilled++;
          }
          else
            unfilled++;
        }
        else {
          *IMAGEFseq_pix(Itmp, u, v, 0) = *IMAGEFseq_pix(Ip, u, v, 0);
          *IMAGEFseq_pix(Itmp, u, v, 1) = *IMAGEFseq_pix(Ip, u, v, 1);
          *IMAGEFseq_pix(Itmp, u, v, 2) = *IMAGEFseq_pix(Ip, u, v, 2);
        }
      }
    }
    for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
      for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
        if (filled[u][v] == FILLING_ELT) filled[u][v] = FILLED_ELT;
      }
    }
    mrisp->Ip = Itmp;
    ImageFree(&Ip);
    if (npasses++ > 1000)
      ErrorExit(ERROR_BADFILE,
                "MRISPtoParameterization: could not fill "
                "parameterization");
  } while (unfilled > 0);
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) fprintf(stderr, "filling %d elements took %d passes\n", nfilled, npasses);

  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    free(filled[u]);
    free(distances[u]);
  }
  free(filled);
  free(distances);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "done.\n");

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) ImageWrite(mrisp->Ip, "sphere.hipl");

  return (mrisp);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIS *MRISfromParameterization(MRI_SP *mrisp, MRIS *mris, int fno)
{
  float a, b, c, phi, theta, x, y, z, uf, vf, du, dv, curv, d;
  int vno, u0, v0, u1, v1;
  VERTEX *vertex;

  if (!mris) mris = MRISclone(mrisp->mris);

  MRIS *original = mris;
  mris = makeCenteredSphere(mris);

  a = b = c = MRISaverageRadius(mris);

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(vertex->y / b, vertex->x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (phi > M_PI) DiagBreak();
    if (vno == 60935) DiagBreak();

    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u0 = floor(uf);
    u1 = ceil(uf);
    v0 = floor(vf);
    v1 = ceil(vf);
    du = uf - (float)u0;
    dv = vf - (float)v0;


    /* enforce spherical topology  */
    if (u0 < 0) /* enforce spherical topology  */
      u0 = -u0;
    if (u0 >= U_DIM(mrisp)) u0 = U_DIM(mrisp) - (u0 - U_DIM(mrisp) + 1);
    if (u1 < 0) /* enforce spherical topology  */
      u1 = -u1;
    if (u1 >= U_DIM(mrisp)) u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
    if (v0 < 0) v0 += V_DIM(mrisp);
    if (v0 >= V_DIM(mrisp)) v0 -= V_DIM(mrisp);
    if (v1 < 0) v1 += V_DIM(mrisp);
    if (v1 >= V_DIM(mrisp)) v1 -= V_DIM(mrisp);

    if (((u0 == DEBUG_U) || (u1 == DEBUG_U)) && ((v0 == DEBUG_V) || (v1 == DEBUG_V))) DiagBreak();
    /* do bilinear interpolation */
    curv = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) +
           (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, fno) +
           (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, fno) +
           du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, fno);

    vertex->curv = curv;
  }

  resetCenteredSphere(original, mris);
  return original;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIS *MRIScoordsFromParameterization(MRI_SP *mrisp, MRIS *mrisInit, int which_vertices)
{
  float rad, phi, theta, uf, vf, du, dv, origx, origy, origz, d, val1, val2, val3, val4;
  int vno, u0, v0, u1, v1;

  MRIS * const mris = (!mrisInit) ? MRISclone(mrisp->mris) : mrisInit;

  if (which_vertices == ORIGINAL_VERTICES)
    cheapAssert(mris->origxyz_status == mrisInit->origxyz_status);
  
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const vertex = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();
    float x = vertex->x;
    float y = vertex->y;
    float z = vertex->z;
    // DNG 7/19/18: changed to use a vertex specific radius rather than the mean radius
    rad = sqrt(x*x + y*y + z*z);
    theta = atan2(vertex->y / rad, vertex->x / rad);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = rad * rad - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (phi > M_PI) {
      DiagBreak();
    }
    if (vno == 60935) DiagBreak();

    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u0 = floor(uf);
    u1 = ceil(uf);
    v0 = floor(vf);
    v1 = ceil(vf);
    du = uf - (float)u0;
    dv = vf - (float)v0;


    /* enforce spherical topology  */
    if (u0 < 0) /* enforce spherical topology  */
      u0 = -u0;
    if (u0 >= U_DIM(mrisp)) u0 = U_DIM(mrisp) - (u0 - U_DIM(mrisp) + 1);
    if (u1 < 0) /* enforce spherical topology  */
      u1 = -u1;
    if (u1 >= U_DIM(mrisp)) u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
    if (v0 < 0) v0 += V_DIM(mrisp);
    if (v0 >= V_DIM(mrisp)) v0 -= V_DIM(mrisp);
    if (v1 < 0) v1 += V_DIM(mrisp);
    if (v1 >= V_DIM(mrisp)) v1 -= V_DIM(mrisp);

    if (((u0 == DEBUG_U) || (u1 == DEBUG_U)) && ((v0 == DEBUG_V) || (v1 == DEBUG_V))) {
      printf(" u,v (%d, %d) --> v %d\n", u0, v0, vno);
      DiagBreak();
    }
    if ((u0 == Gx && v0 == Gy) || (u0 == Gy && v0 == Gx)) DiagBreak();

    /* do bilinear interpolation */

    val1 = *IMAGEFseq_pix(mrisp->Ip, u1, v1, 0);
    val2 = *IMAGEFseq_pix(mrisp->Ip, u0, v1, 0);
    val3 = *IMAGEFseq_pix(mrisp->Ip, u0, v0, 0);
    val4 = *IMAGEFseq_pix(mrisp->Ip, u1, v0, 0);

    origx = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, 0) + (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, 0) +
            (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, 0) +
            du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, 0);
    val1 = *IMAGEFseq_pix(mrisp->Ip, u1, v1, 1);
    val2 = *IMAGEFseq_pix(mrisp->Ip, u0, v1, 1);
    val3 = *IMAGEFseq_pix(mrisp->Ip, u0, v0, 1);
    val4 = *IMAGEFseq_pix(mrisp->Ip, u1, v0, 1);
    origy = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, 1) + (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, 1) +
            (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, 1) +
            du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, 1);
    val1 = *IMAGEFseq_pix(mrisp->Ip, u1, v1, 2);
    val2 = *IMAGEFseq_pix(mrisp->Ip, u0, v1, 2);
    val3 = *IMAGEFseq_pix(mrisp->Ip, u0, v0, 2);
    val4 = *IMAGEFseq_pix(mrisp->Ip, u1, v0, 2);
    origz = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, 2) + (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, 2) +
            (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, 2) +
            du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, 2);

    switch (which_vertices) {
      case WHITE_VERTICES:
	vertex->whitex = origx ;
	vertex->whitey = origy ;
	vertex->whitez = origz ;
        break;
      case CURRENT_VERTICES:
	vertex->x = origx ;
	vertex->y = origy ;
	vertex->z = origz ;
        break;
      case ORIGINAL_VERTICES:
        MRISsetOriginalXYZ(mris, vno, origx, origy, origz);
        break;
      case TMP_VERTICES:
        vertex->tx = origx;
        vertex->ty = origy;
        vertex->tz = origz;
        break;
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRIScoordsFromParameterization: which = %d not supported", which_vertices);
    }
  }

  return (mris);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIS *MRISnormalizeFromParameterization(MRI_SP *mrisp, MRIS *mris, int fno)
{
  float a, b, c, phi, theta, x, y, z, uf, vf, du, dv, curv, d, var;
  int vno, u0, v0, u1, v1;
  VERTEX *vertex;

  if (!mris) mris = MRISclone(mrisp->mris);

  MRIS *original = mris;
  mris = makeCenteredSphere(mris);

  a = b = c = MRISaverageRadius(mris);

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(vertex->y / b, vertex->x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (vno == 60935) DiagBreak();

    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u0 = floor(uf);
    u1 = ceil(uf);
    v0 = floor(vf);
    v1 = ceil(vf);
    du = uf - (float)u0;
    dv = vf - (float)v0;

    /* enforce spherical topology  */
    if (u0 < 0) /* enforce spherical topology  */
      u0 = -u0;
    if (u0 >= U_DIM(mrisp)) u0 = U_DIM(mrisp) - (u0 - U_DIM(mrisp) + 1);
    if (u1 < 0) /* enforce spherical topology  */
      u1 = -u1;
    if (u1 >= U_DIM(mrisp)) u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
    if (v0 < 0) v0 += V_DIM(mrisp);
    if (v0 >= V_DIM(mrisp)) v0 -= V_DIM(mrisp);
    if (v1 < 0) v1 += V_DIM(mrisp);
    if (v1 >= V_DIM(mrisp)) v1 -= V_DIM(mrisp);

    if (((u0 == DEBUG_U) || (u1 == DEBUG_U)) && ((v0 == DEBUG_V) || (v1 == DEBUG_V))) DiagBreak();
    /* do bilinear interpolation */
    curv = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) +
           (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, fno) +
           (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, fno) +
           du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, fno);

    var = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno + 1) +
          (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, fno + 1) +
          (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, fno + 1) +
          du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, fno + 1);

    if (FZERO(var)) var = FSMALL;
    vertex->curv = curv / (sqrt(var));
  }

  resetCenteredSphere(original, mris);
  return original;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *MRISgradientToParameterization(MRIS *mris, MRI_SP *mrisp, float scale)
{
  MRIS *original = mris;
  mris = makeCenteredSphere(mris);

  float a, b, c, phi, theta, x, y, z, uf, vf, d, total_d, **distances, sigma, two_sigma_sq;
  int vno, u, v, unfilled, **filled;
  VERTEX *vertex;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "computing parameterization...");

  if (!mrisp)
    mrisp = MRISPalloc(scale, 3);
  else
    ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, -1);

  a = b = c = MRISaverageRadius(mris);

  filled = (int **)calloc(U_DIM(mrisp), sizeof(int *));
  distances = (float **)calloc(U_DIM(mrisp), sizeof(float *));
  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    filled[u] = (int *)calloc(V_DIM(mrisp), sizeof(int));
    distances[u] = (float *)calloc(V_DIM(mrisp), sizeof(float));

    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) filled[u][v] = UNFILLED_ELT;
  }

  sigma = scale / 4.0f;
  two_sigma_sq = 2 * sigma * sigma;

  /* first calculate total distances to a point in parameter space */
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(y / b, x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (vno == DEBUG_VNO) DiagBreak();
    vertex->phi = phi;
    vertex->theta = theta;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u = nint(uf);
    v = nint(vf);
    if (u < 0) /* enforce spherical topology  */
      u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) /* enforce spherical topology  */
      v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    if (u == 0 && v == 56) DiagBreak();
    if ((u == DEBUG_U) && (v == DEBUG_V)) DiagBreak();

    filled[u][v] = vno;
    distances[u][v] += 1; /* keep track of total # of nodes */
  }

  if (DEBUG_U >= 0) fprintf(stderr, "\ndistance[%d][%d] = %2.3f\n\n", DEBUG_U, DEBUG_V, distances[DEBUG_U][DEBUG_V]);

  /* now add in gradients proportional to their distance from the point */
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(y / b, x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u = nint(uf);
    v = nint(vf);
    if (u < 0) /* enforce spherical topology  */
      u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) /* enforce spherical topology  */
      v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    if (u == 0 && v == 56) DiagBreak();

    /* 0,0 */
    total_d = distances[u][v];
    if (total_d > 0.0) {
      *IMAGEFseq_pix(mrisp->Ip, u, v, 0) += vertex->dx / total_d;
      *IMAGEFseq_pix(mrisp->Ip, u, v, 1) += vertex->dy / total_d;
      *IMAGEFseq_pix(mrisp->Ip, u, v, 2) += vertex->dz / total_d;
    }
    if ((u == DEBUG_U) && (v == DEBUG_V))
      fprintf(stderr,
              "v = %6.6d (%2.1f, %2.1f, %2.1f), curv = %2.3f, "
              "proportion = %2.3f\n",
              vno,
              x,
              y,
              z,
              vertex->curv,
              d);
  }

  /* fill in values which were unmapped by sampling back onto surface */
  for (unfilled = u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    double min_d, radius = mris->radius, xd, yd, zd;
    int min_v = -1;

    if (FZERO(radius)) radius = MRISaverageRadius(mris);

    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
      if (filled[u][v] == UNFILLED_ELT) {
        if (u == 0 && v == 56) DiagBreak();
        unfilled++;
        phi = (double)u * PHI_MAX / PHI_DIM(mrisp);
        theta = (double)v * THETA_MAX / THETA_DIM(mrisp);
        x = radius * sin(phi) * cos(theta);
        y = radius * sin(phi) * sin(theta);
        z = radius * cos(phi);
        for (min_d = 1000.0f, vno = 0; vno < mris->nvertices; vno++) {
          vertex = &mris->vertices[vno];
          if (vertex->ripflag) continue;
          xd = vertex->x - x;
          yd = vertex->y - y;
          zd = vertex->z - z;
          d = sqrt(xd * xd + yd * yd + zd * zd);
          if (d < min_d) {
            min_d = d;
            min_v = vno;
          }
        }
        vertex = &mris->vertices[min_v];
        *IMAGEFseq_pix(mrisp->Ip, u, v, 0) = vertex->dx;
        *IMAGEFseq_pix(mrisp->Ip, u, v, 1) = vertex->dy;
        *IMAGEFseq_pix(mrisp->Ip, u, v, 2) = vertex->dz;
      }
    }
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "%d holes in parameterization filled\n", unfilled);

  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    free(filled[u]);
    free(distances[u]);
  }
  free(filled);
  free(distances);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "done.\n");

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) ImageWrite(mrisp->Ip, "sphere.hipl");

  resetCenteredSphere(original, mris);

  return (mrisp);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIS *MRISgradientFromParameterization(MRI_SP *mrisp, MRIS *mris)
{
  float a, b, c, phi, theta, x, y, z, uf, vf, du, dv, d;
  int vno, u0, v0, u1, v1;
  VERTEX *vertex;

  if (!mris) mris = MRISclone(mrisp->mris);

  MRIS *original = mris;
  mris = makeCenteredSphere(mris);

  a = b = c = MRISaverageRadius(mris);

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(vertex->y / b, vertex->x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (vno == 60935) DiagBreak();

    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u0 = floor(uf);
    u1 = ceil(uf);
    v0 = floor(vf);
    v1 = ceil(vf);
    du = uf - (float)u0;
    dv = vf - (float)v0;


    /* enforce spherical topology  */
    if (u0 < 0) /* enforce spherical topology  */
      u0 = -u0;
    if (u0 >= U_DIM(mrisp)) u0 = U_DIM(mrisp) - (u0 - U_DIM(mrisp) + 1);
    if (u1 < 0) /* enforce spherical topology  */
      u1 = -u1;
    if (u1 >= U_DIM(mrisp)) u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
    if (v0 < 0) v0 += V_DIM(mrisp);
    if (v0 >= V_DIM(mrisp)) v0 -= V_DIM(mrisp);
    if (v1 < 0) v1 += V_DIM(mrisp);
    if (v1 >= V_DIM(mrisp)) v1 -= V_DIM(mrisp);

    if (((u0 == DEBUG_U) || (u1 == DEBUG_U)) && ((v0 == DEBUG_V) || (v1 == DEBUG_V))) DiagBreak();

    /* do bilinear interpolation */
    vertex->dx = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, 0) +
                 (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, 0) +
                 (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, 0) +
                 du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, 0);

    vertex->dy = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, 1) +
                 (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, 1) +
                 (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, 1) +
                 du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, 1);

    vertex->dz = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, 2) +
                 (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, 2) +
                 (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, 2) +
                 du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, 2);
  }

  resetCenteredSphere(original, mris);
  return original;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
void MRISPfunctionVal_radiusR(                                                      // returns the value that would be stored in resultsForEachFno[0] for fnoLo
                              MRI_SURFACE_PARAMETERIZATION *mrisp,  
                              MRISPfunctionValResultForAlpha* resultsForEachAlpha,  // must be numAlphas elements
                              float r, float x, float y, float z, 
                              int fnoLo, bool getNextAlso,                          // always fills in resultsForEachAlpha.curr for fno, optionally fills in .next for fno+1
                              float const * alphas, float numAlphas,                // rotate x,y,z by these alphas (radians) and get the values
                              bool trace)                                           // note: this rotation is around the z axis, hence z does not change
{
  float phi;
  { float d = r * r - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2f(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
  }
  
  float uf = PHI_DIM(mrisp) * phi / PHI_MAX;
  int   u0 = floor(uf);
  int   u1 = ceil(uf);

  float du = uf - (float)u0;

  /* enforce spherical topology  */
  int u0_voff, u1_voff;

  if (u0 < 0) /* enforce spherical topology  */
  {
    u0_voff = V_DIM(mrisp) / 2;
    u0 = -u0;
  }
  else if (u0 >= U_DIM(mrisp)) {
    u0_voff = V_DIM(mrisp) / 2;
    u0 = U_DIM(mrisp) - (u0 - U_DIM(mrisp) + 1);
  }
  else
    u0_voff = 0;
  if (u1 < 0) /* enforce spherical topology  */
  {
    u1_voff = V_DIM(mrisp) / 2;
    u1 = -u1;
  }
  else if (u1 >= U_DIM(mrisp)) {
    u1_voff = V_DIM(mrisp) / 2;
    u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
  }
  else
    u1_voff = 0;

  // This is the rotation around the z axis
  //
  float const baseTheta = fastApproxAtan2f(y, x);

  int alphaIndex;
  for (alphaIndex = 0; alphaIndex < numAlphas; alphaIndex++) {
    float theta = baseTheta - alphas[alphaIndex];
    while (theta <  0.0f  ) theta += 2*M_PI; /* make it 0 --> 2*PI */
    while (theta >= 2*M_PI) theta -= 2*M_PI; /* make it 0 --> 2*PI */

    float vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    int   v0 = floor(vf);
    int   v1 = ceil(vf);

    float dv = vf - (float)v0;

    if (v0 < 0            ) v0 += V_DIM(mrisp);
    if (v0 >= V_DIM(mrisp)) v0 -= V_DIM(mrisp);
    if (v1 < 0            ) v1 += V_DIM(mrisp);
    if (v1 >= V_DIM(mrisp)) v1 -= V_DIM(mrisp);

    int u0_v1 = v1 + u0_voff; while (u0_v1 >= V_DIM(mrisp)) u0_v1 -= V_DIM(mrisp);
    int u0_v0 = v0 + u0_voff; while (u0_v0 >= V_DIM(mrisp)) u0_v0 -= V_DIM(mrisp);
    int u1_v1 = v1 + u1_voff; while (u1_v1 >= V_DIM(mrisp)) u1_v1 -= V_DIM(mrisp);
    int u1_v0 = v0 + u1_voff; while (u1_v0 >= V_DIM(mrisp)) u1_v0 -= V_DIM(mrisp);

    if (0) {
      static long statsCount, statsLimit = 1, statsV_DIMSum;

      statsV_DIMSum += V_DIM(mrisp);
      statsCount++;
      if (statsCount == statsLimit) {
        if (statsLimit < 10000000) statsLimit *= 2; else statsLimit += 10000000;
        fprintf(stdout, "%s:%d statsCount:%g V_DIM avg:%g\n", __FILE__, __LINE__,
          (float)statsCount, (double)statsV_DIMSum/(double)statsCount);
      }
    }

    int i;
    for (i = 0; i < (getNextAlso?2:1); i++) {
      int fno = fnoLo + i;
      /* do bilinear interpolation */
      double val = 
                    du  *         dv  * *IMAGEFseq_pix(mrisp->Ip, u1, u1_v1, fno) +
            (1.0f - du) *         dv  * *IMAGEFseq_pix(mrisp->Ip, u0, u0_v1, fno) +
            (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, u0_v0, fno) +
                    du  * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, u1_v0, fno);
      if (trace) {
        fprintf(stdout, "%s:%d x:%g y:%g z:%g returns %g for fno:%d phi:%g alphaIndex:%d theta:%g\n", 
            __FILE__, __LINE__, 
            x,y,z,val,fno,phi,alphaIndex,theta);
      }
      *(i ? &resultsForEachAlpha[alphaIndex].next : &resultsForEachAlpha[alphaIndex].curr) = val;
    }
  }
}


double MRISPfunctionValTraceable(MRI_SURFACE_PARAMETERIZATION *mrisp, float desired_radius, float x, float y, float z, int fno, bool trace)
{
  double r = sqrt(x * x + y * y + z * z);

  if (!FEQUAL(r, desired_radius)) /* project it onto sphere */
  {
    r = desired_radius;
    double r2 = r * r;
    double r4 = r2 * r2;
    double r6 = r2 * r4;
    double x2 = x * x;
    double y2 = y * y;
    double z2 = z * z;
    double f  = x2 / r6 + y2 / r6 + z2 / r6;
    double g  = 2 * (x2 / r4 + y2 / r4 + z2 / r4);
    double h  = x2 / r2 + y2 / r2 + z2 / r2 - 1;
    double d  = (-g + (float)sqrt((double)(g * g - 4 * f * h))) / (2 * f);
    double dx = d * x / r2;
    double dy = d * y / r2;
    double dz = d * z / r2;
    x = x + dx;
    y = y + dy;
    z = z + dz;
  }

  float zero = 0.0f;
  MRISPfunctionValResultForAlpha result;
  MRISPfunctionVal_radiusR(mrisp, &result, r, x, y, z, fno, false, &zero, 1, trace);
  return result.curr;
}

double MRISPfunctionVal(MRI_SURFACE_PARAMETERIZATION *mrisp, float desired_radius, float x, float y, float z, int fno) {
    return MRISPfunctionValTraceable(mrisp, desired_radius, x, y, z, fno, false);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *MRISPclone(MRI_SP *mrisp_src)
{
  MRI_SP *mrisp_dst;

  mrisp_dst = (MRI_SP *)calloc(1, sizeof(MRI_SP));
  mrisp_dst->mris = mrisp_src->mris;
  mrisp_dst->scale = mrisp_src->scale ;
  mrisp_dst->Ip = ImageCopy(mrisp_src->Ip, NULL);

  return (mrisp_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Convolve with a kernel with standard deviation = sigma mm.
------------------------------------------------------*/
#define MAX_LEN 4
#define MAX_KLEN 50

MRI_SP *MRISPconvolveGaussian(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma, float radius, int fno)
{
  int u, v, cart_klen, klen, khalf, uk, vk, u1, v1, voff, f0, f1;
  double d, k, total, ktotal, sigma_sq_inv, theta, phi, theta1, phi1, sin_phi, cos_phi, sin_phi1, cos_phi1;
  float x0, y0, z0, x1, y1, z1, circumference = 0.0f, angle, max_len = 0.0f, min_len = 10000.0f;
  IMAGE *Ip_src, *Ip_dst;
  VECTOR *vec1, *vec2;

  vec1 = VectorAlloc(3, MATRIX_REAL);
  vec2 = VectorAlloc(3, MATRIX_REAL);

  if (!mrisp_dst) mrisp_dst = MRISPclone(mrisp_src);
  mrisp_dst->sigma = sigma;

  /* determine the size of the kernel */
  cart_klen = (int)nint(6.0f * sigma) + 1;
  if (ISEVEN(cart_klen)) /* ensure it's odd */
    cart_klen++;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "blurring surface, sigma = %2.3f, cartesian klen = %d\n", sigma, cart_klen);

  if (FZERO(sigma))
    sigma_sq_inv = BIG;
  else
    sigma_sq_inv = 1.0f / (sigma * sigma);

  Ip_src = mrisp_src->Ip;
  Ip_dst = mrisp_dst->Ip;
  if (fno < 0) {
    f0 = 0;
    f1 = Ip_src->num_frame - 1;
  }
  else {
    f0 = f1 = fno;
  }

  for (fno = f0; fno <= f1; fno++) /* for each frame */
  {
    for (u = 0; u < U_DIM(mrisp_src); u++) {
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "\r%3.3d of %d     ", u, U_DIM(mrisp_src) - 1);
      phi = (double)u * PHI_MAX / PHI_DIM(mrisp_src);
      sin_phi = sin(phi);
      cos_phi = cos(phi);

      for (v = 0; v < V_DIM(mrisp_src); v++) {
        theta = (double)v * THETA_MAX / THETA_DIM(mrisp_src);
        x0 = radius * sin_phi * cos(theta);
        y0 = radius * sin_phi * sin(theta);
        z0 = radius * cos_phi;
        VECTOR_LOAD(vec1, x0, y0, z0); /* radius vector */
        if (FZERO(circumference))      /* only calculate once */
          circumference = M_PI * 2.0 * V3_LEN(vec1);
        if (u == DEBUG_U && v == DEBUG_V) DiagBreak();

        /* compute the distance between adjacent spherical matrix
           elements at this point on the surface (probably easier
           to do with the parameterization, but I'll do it in
           Cartesian space for now.
           */
        u1 = u + 1;
        if (u1 >= U_DIM(mrisp_src)) u1 = U_DIM(mrisp_src) - (u1 - U_DIM(mrisp_src) + 2);
        v1 = v + 1;
        if (v1 >= V_DIM(mrisp_src)) v1 = V_DIM(mrisp_src) - (v1 - V_DIM(mrisp_src) + 2);

        phi1 = (double)u1 * PHI_MAX / PHI_DIM(mrisp_src);
        theta1 = (double)v1 * THETA_MAX / THETA_DIM(mrisp_src);
        x1 = radius * sin(phi1) * cos(theta1);
        y1 = radius * sin(phi1) * sin(theta1);
        z1 = radius * cos(phi1);
        VECTOR_LOAD(vec2, x1, y1, z1); /* radius vector */
        angle = fabs(Vector3Angle(vec1, vec2));
        d = circumference * angle / (2.0 * M_PI); /* geodesic distance */
        if (d > max_len) max_len = d;
        if (d < min_len) min_len = d;

        /* d is now the distance between adjacent cells - compute kernel size*/
        klen = nint(6.0f * sigma / d) + 1;
        if (klen > MAX_KLEN) klen = MAX_KLEN;

        if (ISEVEN(klen)) klen++;
        if (klen >= U_DIM(mrisp_src)) klen = U_DIM(mrisp_src) - 1;
        if (klen >= V_DIM(mrisp_src)) klen = V_DIM(mrisp_src) - 1;
        khalf = klen / 2;

        total = ktotal = 0.0;
        for (uk = -khalf; uk <= khalf; uk++) {
          u1 = u + uk;
          if (u1 < 0) /* enforce spherical topology  */
          {
            voff = V_DIM(mrisp_src) / 2;
            u1 = -u1;
          }
          else if (u1 >= U_DIM(mrisp_src)) {
            u1 = U_DIM(mrisp_src) - (u1 - U_DIM(mrisp_src) + 1);
            voff = V_DIM(mrisp_src) / 2;
          }
          else
            voff = 0;

          phi1 = (double)u1 * PHI_MAX / PHI_DIM(mrisp_src);
          sin_phi1 = sin(phi1);
          cos_phi1 = cos(phi1);

          for (vk = -khalf; vk <= khalf; vk++) {
            theta1 = (double)v * THETA_MAX / THETA_DIM(mrisp_src);
            x1 = radius * sin_phi1 * cos(theta1);
            y1 = radius * sin_phi1 * sin(theta1);
            z1 = radius * cos_phi1;
            VECTOR_LOAD(vec2, x1, y1, z1); /* radius vector */
            angle = fabs(Vector3Angle(vec1, vec2));
            d = circumference * angle / (2.0 * M_PI);
            k = exp(-d * d * sigma_sq_inv);
            v1 = v + vk + voff;
            while (v1 < 0) /* enforce spherical topology */
              v1 += V_DIM(mrisp_src);
            while (v1 >= V_DIM(mrisp_src)) v1 -= V_DIM(mrisp_src);
            ktotal += k;
            total += k * *IMAGEFseq_pix(Ip_src, u1, v1, fno);
          }
        }
        if (u == DEBUG_U && v == DEBUG_V) DiagBreak();
        total /= ktotal; /* normalize weights to 1 */
        *IMAGEFseq_pix(mrisp_dst->Ip, u, v, fno) = total;
      }
    }
  }

  if (Gdiag & DIAG_SHOW) fprintf(stderr, "min_len = %2.3f mm, max_len = %2.3f mm\n", min_len, max_len);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "done.\n");

  VectorFree(&vec1);
  VectorFree(&vec2);
  return (mrisp_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static MRI_SP *MRISPblur_new(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma, int fno);
static MRI_SP *MRISPblur_old(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma, int fno);

MRI_SP *MRISPblur(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma, int fno) {
    static bool once, do_old;
    if (!once) { once = true;
        do_old = getenv("FREESUREFER_MRISPblur_old");
    }
    return 
        (do_old ? MRISPblur_old : MRISPblur_new)(mrisp_src, mrisp_dst, sigma, fno);
}

static MRI_SP *MRISPblur_new(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma, int fno) {
  int fnoLo_init, fnoHi_init;
  int no_sphere_init;
  int cart_klen_init;
  double sigma_sq_inv_init;
  IMAGE *Ip_src_init;
  {
    no_sphere_init = getenv("NO_SPHERE") != NULL;
    if (no_sphere_init) fprintf(stderr, "disabling spherical geometry\n");

    if (!mrisp_dst) mrisp_dst = MRISPclone(mrisp_src);
    mrisp_dst->sigma = sigma;

    /* determine the size of the kernel */
    cart_klen_init = (int)nint(6.0f * sigma) + 1;
    if (ISEVEN(cart_klen_init)) /* ensure it's odd */
      cart_klen_init++;

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "blurring surface, sigma = %2.3f, cartesian klen = %d\n", sigma, cart_klen_init);

    if (FZERO(sigma))
      sigma_sq_inv_init = BIG;
    else
      sigma_sq_inv_init = 1.0f / (sigma * sigma);

    Ip_src_init = mrisp_src->Ip ;

    if (fno < 0) {
      fnoLo_init = 0;
      fnoHi_init =  Ip_src_init->num_frame;
    } else {
      fnoLo_init = fno;
      fnoHi_init = fno + 1;
    }
  }

  int const          fnoLo        = fnoLo_init;
  int const          fnoHi        = fnoHi_init;
  int const          cart_klen    = cart_klen_init;
  int const          no_sphere    = no_sphere_init;
  double const       sigma_sq_inv = sigma_sq_inv_init;
  const IMAGE *const Ip_src       = Ip_src_init;

  static int  uMax                    = 0;
  static int  uToKHalfCache_PHI_DIM   = 0;
  static int  uToKHalfCache_cart_klen = 0;
  static int* uToKHalfCache           = NULL;
  
  static int     kHalfHi              = 0;
  static double* uukvkToExpResult     = NULL;
  
  if (uMax                    != U_DIM  (mrisp_src) 
  ||  uToKHalfCache_PHI_DIM   != PHI_DIM(mrisp_src)
  ||  uToKHalfCache_cart_klen != cart_klen
  ) {
    uMax                    = U_DIM  (mrisp_src);
    uToKHalfCache_PHI_DIM   = PHI_DIM(mrisp_src);
    uToKHalfCache_cart_klen = cart_klen;
    
    uToKHalfCache = (int*)realloc(uToKHalfCache, uMax*sizeof(int));
    
    int u;
    for (u = 0; u < U_DIM(mrisp_src); u++) {
      double const phi      = (double)u * PHI_MAX / PHI_DIM(mrisp_src);
      double const sin_sq_u = no_sphere ? 1.0 : squared(sin(phi));

      int klen;

      if (no_sphere) {
        klen = cart_klen;

      } else if (!FZERO(sin_sq_u)) {
        int k = cart_klen * cart_klen;
        klen = sqrt(k + k / sin_sq_u);
        if (klen > MAX_LEN * cart_klen) klen = MAX_LEN * cart_klen;

      } else {
        klen = MAX_LEN * cart_klen; /* arbitrary max length */          // happens when u is 0
      }

      if (klen >= U_DIM(mrisp_src)) klen = U_DIM(mrisp_src) - 1;
      if (klen >= V_DIM(mrisp_src)) klen = V_DIM(mrisp_src) - 1;

      uToKHalfCache[u] = klen / 2;
      
      if (u > 0 && uToKHalfCache[u] > uToKHalfCache[0]) {
        fprintf(stdout, "%s:%d uMax[0] not max\n", __FILE__, __LINE__);
        exit(1);
      }
    }
    
    kHalfHi          = uToKHalfCache[0] + 1;
    uukvkToExpResult = (double*)realloc(uukvkToExpResult, uMax*kHalfHi*kHalfHi*sizeof(double));
    
    for (u = 0; u < U_DIM(mrisp_src); u++) {
      double* ukvkToExpResult = uukvkToExpResult + u*kHalfHi*kHalfHi;

      double const phi      = (double)u * PHI_MAX / PHI_DIM(mrisp_src);
      double const sin_sq_u = no_sphere ? 1.0 : squared(sin(phi));
       
      int uk,vk;
      for (uk = 0; uk < kHalfHi; uk++) 
      for (vk = 0; vk < kHalfHi; vk++) {
        double const udiff = (double)(uk * uk); /* distance squared in u */
        double const vdiff = (double)(vk * vk); /* distance squared in v */
        double neg_exp_input = (udiff + sin_sq_u * vdiff) * sigma_sq_inv;
        double k             = exp(-neg_exp_input);
        ukvkToExpResult[uk*kHalfHi + vk] = k;
      }
    }
  }
  
  int u;

  ROMP_PF_begin
#if HAVE_OPENMP  
  #pragma omp parallel for if_ROMP(assume_reproducible) collapse(2)
#endif
  for (fno = fnoLo; fno < fnoHi; fno++) /* for each frame */
  {
    for (u = 0; u < U_DIM(mrisp_src); u++) {
      ROMP_PFLB_begin
      double const * const ukvkToExpResult = uukvkToExpResult + u*kHalfHi*kHalfHi;
      
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "\r%3.3d of %d     ", u, U_DIM(mrisp_src) - 1);
      
      double const phi      = (double)u * PHI_MAX / PHI_DIM(mrisp_src);
      double const sin_sq_u = no_sphere ? 1.0 : squared(sin(phi));
      
      // khalf is independent of fno and of v.
      //
      int khalf_init;
      {
        int klen;

        if (no_sphere) {
          klen = cart_klen;

        } else if (!FZERO(sin_sq_u)) {
          int k = cart_klen * cart_klen;
          klen = sqrt(k + k / sin_sq_u);
          if (klen > MAX_LEN * cart_klen) klen = MAX_LEN * cart_klen;
        
        } else {
          klen = MAX_LEN * cart_klen; /* arbitrary max length */
        }

        if (klen >= U_DIM(mrisp_src)) klen = U_DIM(mrisp_src) - 1;
        if (klen >= V_DIM(mrisp_src)) klen = V_DIM(mrisp_src) - 1;
        
        khalf_init = klen / 2;
      }
      int const khalf = khalf_init;

      if (khalf != uToKHalfCache[u]) {
        fprintf(stdout, "%s:%d khalf != uToKHalfCache[u]\n",__FILE__, __LINE__);
        exit(1);
      }
      
      int v;
      for (v = 0; v < V_DIM(mrisp_src); v++) {
        /*      theta = (double)v*THETA_MAX / THETA_DIM(mrisp_src) ;*/
        // if (u == DEBUG_U && v == DEBUG_V) DiagBreak();

        double total = 0.0, ktotal = 0.0;

        int uk;
        for (uk = -khalf; uk <= khalf; uk++) {
          int const absUk = (uk < 0) ? -uk : uk; 
#ifdef CHECK_TABLE
          double const udiff = (double)(uk * uk); /* distance squared in u */
#endif

          int voff;

          int u1 = u + uk;
          if (u1 < 0) /* enforce spherical topology  */
          {
            voff = V_DIM(mrisp_src) / 2;
            u1 = -u1;
          }
          else if (u1 >= U_DIM(mrisp_src)) {
            u1 = U_DIM(mrisp_src) - (u1 - U_DIM(mrisp_src) + 1);
            voff = V_DIM(mrisp_src) / 2;
          } else {
            voff = 0;
          }
          
          int vk;
          for (vk = -khalf; vk <= khalf; vk++) {
            int const absVk = (vk < 0) ? -vk : vk; 
            double kFromTable = ukvkToExpResult[absUk*kHalfHi + absVk];

#ifndef CHECK_TABLE
            double k = kFromTable;
#else
            double vdiff = (double)(vk * vk);
            
            double neg_exp_input = (udiff + sin_sq_u * vdiff) * sigma_sq_inv;
                //
                // udiff is v invariant
                // 
                
            double k = exp(-neg_exp_input);

            if (k != kFromTable) {
              fprintf(stdout, "%s:%d k:%f != kFromTable:%f uk:%d vk:%d\n", __FILE__, __LINE__, k, kFromTable, uk, vk);
              exit(1);
            }
#endif
            
            int v1 = v + vk + voff;
            while (v1 < 0) /* enforce spherical topology */
              v1 += V_DIM(mrisp_src);
            while (v1 >= V_DIM(mrisp_src)) v1 -= V_DIM(mrisp_src);
            
            ktotal += k;
            total  += k * *IMAGEFseq_pix(Ip_src, u1, v1, fno);
          }
        }
        if (u == DEBUG_U && v == DEBUG_V) DiagBreak();
        total /= ktotal; /* normalize weights to 1 */
        *IMAGEFseq_pix(mrisp_dst->Ip, u, v, fno) = total;
      }
      ROMP_PFLB_end
    }
  }
  ROMP_PF_end

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "done.\n");

  return (mrisp_dst);
}

static MRI_SP *MRISPblur_old(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma, int fno)
{
  int f0, f1;
  int no_sphere_init;
  int cart_klen_init;
  double sigma_sq_inv_init;
  IMAGE *Ip_src_init;
  {
    int no_sphere, cart_klen;
    double sigma_sq_inv;
    IMAGE *Ip_src, *Ip_dst;

    no_sphere = getenv("NO_SPHERE") != NULL;
    if (no_sphere) fprintf(stderr, "disabling spherical geometry\n");

    if (!mrisp_dst) mrisp_dst = MRISPclone(mrisp_src);
    mrisp_dst->sigma = sigma;

    /* determine the size of the kernel */
    cart_klen = (int)nint(6.0f * sigma) + 1;
    if (ISEVEN(cart_klen)) /* ensure it's odd */
      cart_klen++;

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "blurring surface, sigma = %2.3f, cartesian klen = %d\n", sigma, cart_klen);

    if (FZERO(sigma))
      sigma_sq_inv = BIG;
    else
      sigma_sq_inv = 1.0f / (sigma * sigma);

    Ip_src = mrisp_src->Ip;
    Ip_dst = mrisp_dst->Ip;
    if (fno < 0) {
      f0 = 0;
      f1 = Ip_src->num_frame - 1;
    }
    else {
      f0 = f1 = fno;
    }

    cart_klen_init = cart_klen;
    Ip_src_init = Ip_src;
    no_sphere_init = no_sphere;
    sigma_sq_inv_init = sigma_sq_inv;
  }
  const int cart_klen = cart_klen_init;
  const int no_sphere = no_sphere_init;
  const double sigma_sq_inv = sigma_sq_inv_init;
  const IMAGE *const Ip_src = Ip_src_init;

  int u;

  ROMP_PF_begin
#if HAVE_OPENMP  
  #pragma omp parallel for if_ROMP(assume_reproducible) collapse(2)
#endif
  for (fno = f0; fno <= f1; fno++) /* for each frame */
  {
    for (u = 0; u < U_DIM(mrisp_src); u++) {
      ROMP_PFLB_begin
      int k, klen, khalf;
      double phi, sin_sq_u;

      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "\r%3.3d of %d     ", u, U_DIM(mrisp_src) - 1);
      phi = (double)u * PHI_MAX / PHI_DIM(mrisp_src);
      sin_sq_u = sin(phi);
      sin_sq_u *= sin_sq_u;
      if (!FZERO(sin_sq_u)) {
        k = cart_klen * cart_klen;
        klen = sqrt(k + k / sin_sq_u);
        if (klen > MAX_LEN * cart_klen) klen = MAX_LEN * cart_klen;
      }
      else
        klen = MAX_LEN * cart_klen; /* arbitrary max length */
      if (no_sphere) sin_sq_u = 1.0f, klen = cart_klen;
      if (klen >= U_DIM(mrisp_src)) klen = U_DIM(mrisp_src) - 1;
      if (klen >= V_DIM(mrisp_src)) klen = V_DIM(mrisp_src) - 1;
      khalf = klen / 2;
      int v;
      for (v = 0; v < V_DIM(mrisp_src); v++) {
        /*      theta = (double)v*THETA_MAX / THETA_DIM(mrisp_src) ;*/
        if (u == DEBUG_U && v == DEBUG_V) DiagBreak();

        double total, ktotal;
        int uk;
        total = ktotal = 0.0;
        for (uk = -khalf; uk <= khalf; uk++) {
          int voff, u1;
          double udiff;
          udiff = (double)(uk * uk); /* distance squared in u */

          u1 = u + uk;
          if (u1 < 0) /* enforce spherical topology  */
          {
            voff = V_DIM(mrisp_src) / 2;
            u1 = -u1;
          }
          else if (u1 >= U_DIM(mrisp_src)) {
            u1 = U_DIM(mrisp_src) - (u1 - U_DIM(mrisp_src) + 1);
            voff = V_DIM(mrisp_src) / 2;
          }
          else
            voff = 0;


          int vk;
          for (vk = -khalf; vk <= khalf; vk++) {
            int v1;
            double vdiff, k;
            vdiff = (double)(vk * vk);
            k = exp(-(udiff + sin_sq_u * vdiff) * sigma_sq_inv);
            v1 = v + vk + voff;
            while (v1 < 0) /* enforce spherical topology */
              v1 += V_DIM(mrisp_src);
            while (v1 >= V_DIM(mrisp_src)) v1 -= V_DIM(mrisp_src);
            ktotal += k;
            total += k * *IMAGEFseq_pix(Ip_src, u1, v1, fno);
          }
        }
        if (u == DEBUG_U && v == DEBUG_V) DiagBreak();
        total /= ktotal; /* normalize weights to 1 */
        *IMAGEFseq_pix(mrisp_dst->Ip, u, v, fno) = total;
      }
      ROMP_PFLB_end
    }
  }
  ROMP_PF_end

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "done.\n");

  return (mrisp_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *MRISPalloc(float scale, int nfuncs)
{
  MRI_SP *mrisp;
  int u_dim, v_dim;

  if (FZERO(scale)) scale = 1.0f;

  u_dim = nint(scale * DEFAULT_UDIM); /* for fft */
  v_dim = 2 * u_dim;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "allocating %d by %d parameterization\n", u_dim, v_dim);
  mrisp = (MRI_SP *)calloc(1, sizeof(MRI_SP));
  if (!mrisp) ErrorExit(ERROR_NOMEMORY, "MRISPalloc(%d, %d): allocation failed", u_dim, v_dim);
  mrisp->Ip = ImageAlloc(v_dim, u_dim, PFFLOAT, nfuncs);
  if (!mrisp->Ip) ErrorExit(ERROR_NOMEMORY, "MRISPalloc(%d, %d): allocation failed", u_dim, v_dim);

  mrisp->scale = scale ;
  return (mrisp);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRISPfree(MRI_SP **pmrisp)
{
  MRI_SP *mrisp;

  mrisp = *pmrisp;
  *pmrisp = NULL;
  ImageFree(&mrisp->Ip);
  free(mrisp);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *MRISPalign(MRI_SP *mrisp_orig, MRI_SP *mrisp_src, MRI_SP *mrisp_tmp, MRI_SP *mrisp_dst)
{
  ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRISPalign: not implemented."));
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *MRISPtranslate(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, int du, int dv)
{
  int u, v, udst, vdst, voff;

  if (!mrisp_dst) mrisp_dst = MRISPclone(mrisp_src);

  /* now form the destination as a translated copy of the source */
  for (u = 0; u <= U_MAX_INDEX(mrisp_src); u++) {
    udst = u + du;
    if (udst < 0) /* enforce spherical topology  */
    {
      voff = V_DIM(mrisp_src) / 2;
      udst = -udst;
    }
    else if (udst >= U_DIM(mrisp_src)) {
      voff = V_DIM(mrisp_src) / 2;
      udst = U_DIM(mrisp_src) - (udst - U_DIM(mrisp_src) + 1);
    }
    else
      voff = 0;
    for (v = 0; v <= V_MAX_INDEX(mrisp_src); v++) {
      vdst = v + dv + voff;
      while (vdst < 0) /* enforce spherical topology  */
        vdst += V_DIM(mrisp_src);
      while (vdst >= V_DIM(mrisp_src)) vdst -= V_DIM(mrisp_src);
      *IMAGEFpix(mrisp_dst->Ip, udst, vdst) = *IMAGEFpix(mrisp_src->Ip, u, v);
    }
  }
  return (mrisp_dst);
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
MRI_SP *MRISPcombine(MRI_SP *mrisp, MRI_SP *mrisp_template, int fno)
{
  int u, v, nframes;
  float dof, mean, vart, var, val;
  IMAGE *Ip, *Ipt;

  if (!mrisp_template) mrisp_template = MRISPclone(mrisp);

  Ip = mrisp->Ip;
  Ipt = mrisp_template->Ip;
  if ((Ip->rows != Ipt->rows) || (Ip->cols != Ipt->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MRISPcombine: size mismatch (%d x %d) vs (%d x %d)",
                 Ip->rows,
                 Ip->cols,
                 Ipt->rows,
                 Ipt->cols));

  nframes = mrisp_template->Ip->num_frame;
  dof = *IMAGEFseq_pix(Ipt, 0, 0, fno + 2);
  *IMAGEFseq_pix(Ipt, 0, 0, fno + 2) += 1.0f;
  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
      val = *IMAGEFpix(Ip, u, v);
      mean = *IMAGEFseq_pix(Ipt, u, v, fno);
      mean = (mean * dof + val) / (dof + 1);

      vart = *IMAGEFseq_pix(Ipt, u, v, fno + 1);
      var = (val - mean);
      var *= var;
      var = (vart * dof + var) / (dof + 1);
      *IMAGEFseq_pix(Ipt, u, v, fno + 1) = var;

      *IMAGEFseq_pix(Ipt, u, v, fno) = mean;
    }
  }

  return (mrisp_template);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRISPwrite(MRI_SP *mrisp, const char *fname)
{
  char ext[STRLEN] ;
  
  if (!strcmp(FileNameExtension(fname, ext), "mgz"))
  {
    MRI *mri ;
    int r, c, f ;

    printf("writing MRISP to mgh file\n") ;
    mri = MRIalloc(mrisp->Ip->cols, mrisp->Ip->rows, mrisp->Ip->num_frame, MRI_FLOAT) ;
    for (f = 0 ; f < mrisp->Ip->num_frame ; f++)
      for (r = 0 ; r < mrisp->Ip->rows ; r++)
	for (c = 0 ; c < mrisp->Ip->cols ; c++)
	  MRIsetVoxVal(mri, c, r, f, 0, *IMAGEFseq_pix(mrisp->Ip, c, r, f)) ;
    MRIwrite(mri, fname) ;
    MRIfree(&mri) ;
  }
  else
    ImageWrite(mrisp->Ip, fname);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *MRISPread(char *fname)
{
  MRI_SP *mrisp;
  int    type ;

  mrisp = (MRI_SP *)calloc(1, sizeof(MRI_SP));
  if (!mrisp) ErrorExit(ERROR_NOMEMORY, "MRISPread(%s): allocation failed", fname);

  type = mri_identify(fname) ;
  if (type == MRI_MGH_FILE || type == NII_FILE)
  {
    MRI *mri ;
    int r, c, f ;

    printf("reading MRISP from mgh file\n") ;
    mri = MRIread(fname) ;
    if (mri == NULL)
      ErrorReturn(NULL, (ERROR_NOFILE, "MRISPread(%s): could not read MGH file",fname));
    mrisp->Ip = ImageAlloc(mri->height, mri->width, PFFLOAT, mri->depth);
    mrisp->scale = (float)mri->width/DEFAULT_UDIM ;
    for (f = 0 ; f < mrisp->Ip->num_frame ; f++)
      for (r = 0 ; r < mrisp->Ip->rows ; r++)
	for (c = 0 ; c < mrisp->Ip->cols ; c++)
	{
	  float val ;
	  if (c == Gx && r == Gy)
	    DiagBreak() ;
	  val = MRIgetVoxVal(mri, c, r, f, 0) ;
	  if (devFinite(val) == 0)
	    val = 0 ;
	  if (val > 0)
	    DiagBreak() ;
	  *IMAGEFseq_pix(mrisp->Ip, c, r, f) = val ;
	}
    MRIfree(&mri) ;
  }
  else
    mrisp->Ip = ImageRead(fname);
  if (!mrisp->Ip) ErrorReturn(NULL, (ERROR_NOFILE, "MRISPread(%s): could not open file", fname));

  return (mrisp);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double MRISParea(MRI_SP *mrisp)
{
  double total_area = 0.0, phi0, theta0, phi1, theta1, radius, l1, l2, x0, y0, z0, x1, y1, z1, x2, y2, z2;
  float val, max_len = 0.0;
  int u, v;

  radius = 1.0;
  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
      val = *IMAGEFpix(mrisp->Ip, u, v);
      if (FZERO(val)) continue;

      phi0 = (double)u * PHI_MAX / PHI_DIM(mrisp);
      theta0 = (double)v * THETA_MAX / THETA_DIM(mrisp);
      phi1 = (double)(u + 1) * PHI_MAX / PHI_DIM(mrisp);
      theta1 = (double)(v + 1) * THETA_MAX / THETA_DIM(mrisp);
      x0 = radius * sin(phi0) * cos(theta0);
      y0 = radius * sin(phi0) * sin(theta0);
      z0 = radius * cos(phi0);
      x1 = radius * sin(phi0) * cos(theta1);
      y1 = radius * sin(phi0) * sin(theta1);
      z1 = radius * cos(phi0);
      x2 = radius * sin(phi1) * cos(theta1);
      y2 = radius * sin(phi1) * sin(theta1);
      z2 = radius * cos(phi1);
      l1 = sqrt(SQR(x1 - x0) + SQR(y1 - y0) + SQR(z1 - z0));
      l2 = sqrt(SQR(x2 - x1) + SQR(y2 - y1) + SQR(z2 - z1));
      if (max_len < l1) max_len = l1;
      if (max_len < l2) max_len = l2;
      total_area += l1 * l2;
    }
  }
  if (Gdiag & DIAG_SHOW) fprintf(stderr, "max parameterization len = %2.3f\n", max_len);
  return (total_area);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *MRISPorLabel(MRI_SP *mrisp, MRIS *mris, LABEL *area)
{
  int n, vno, u, v;
  float r, x, y, z, phi, theta, d, rsq;

  if (!mrisp) mrisp = MRISPalloc(1, 1);

  r = mris->radius;
  rsq = r * r;
  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    x = area->lv[n].x;
    y = area->lv[n].y;
    z = area->lv[n].z;
    d = rsq - z * z;
    if (d < 0.0) d = 0;
    phi = atan2(sqrt(d), z);
    theta = atan2(y / r, x / r);
    u = nint(PHI_DIM(mrisp) * phi / PHI_MAX);
    v = nint(THETA_DIM(mrisp) * theta / THETA_MAX);
    *IMAGEFpix(mrisp->Ip, u, v) = 1.0f;
  }
  return (mrisp);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SP *MRISPandLabel(MRI_SP *mrisp, MRIS *mris, LABEL *area)
{
  int n, vno, u, v;
  float r, x, y, z, phi, theta, d, rsq;

  if (!mrisp) {
    mrisp = MRISPalloc(1, 1);
    ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 1, -1);
  }

  r = mris->radius;
  rsq = r * r;
  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    x = area->lv[n].x;
    y = area->lv[n].y;
    z = area->lv[n].z;
    d = rsq - z * z;
    if (d < 0.0) d = 0;
    phi = atan2(sqrt(d), z);
    theta = atan2(y / r, x / r);
    u = nint(PHI_DIM(mrisp) * phi / PHI_MAX);
    v = nint(THETA_DIM(mrisp) * theta / THETA_MAX);
    *IMAGEFpix(mrisp->Ip, u, v) += 1.0f;
  }
  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
      if (*IMAGEFpix(mrisp->Ip, u, v) > 1.5) /* at least two on */
        *IMAGEFpix(mrisp->Ip, u, v) = 1.0f;
      else
        *IMAGEFpix(mrisp->Ip, u, v) = 0.0f;
    }
  }
  return (mrisp);
}

int MRISPcoordinate(MRI_SP *mrisp, float x, float y, float z, int *pu, int *pv)
{
  double phi, theta, uf, vf;
  int u, v;

  spherical_coordinate(x, y, z, &phi, &theta);
  uf = PHI_DIM(mrisp) * phi / PHI_MAX;
  vf = THETA_DIM(mrisp) * theta / THETA_MAX;
  u = nint(uf);
  v = nint(vf);
  if (u < 0) /* enforce spherical topology  */
    u = -u;
  if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
  if (v < 0) /* enforce spherical topology  */
    v += V_DIM(mrisp);
  if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);
  *pu = u;
  *pv = v;
  return (NO_ERROR);
}
static int spherical_coordinate(double x, double y, double z, double *pphi, double *ptheta)
{
  double r, d;

  r = sqrt(x * x + y * y + z * z);
  d = r * r - z * z;
  if (d < 0.0) d = 0.0;

  *pphi = atan2(sqrt(d), z);
  *ptheta = atan2(y / r, x / r);
  if (*ptheta < 0.0f) *ptheta = 2 * M_PI + *ptheta; /* make it 0 --> 2*PI */
  return (NO_ERROR);
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
MRI_SP *MRISPaccumulate(MRI_SP *mrisp, MRI_SP *mrisp_template, int fno)
{
  int u, v, nframes;
  IMAGE *Ip, *Ipt;

  if (!mrisp_template) mrisp_template = MRISPclone(mrisp);

  Ip = mrisp->Ip;
  Ipt = mrisp_template->Ip;
  if ((Ip->rows != Ipt->rows) || (Ip->cols != Ipt->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MRISPcombine: size mismatch (%d x %d) vs (%d x %d)",
                 Ip->rows,
                 Ip->cols,
                 Ipt->rows,
                 Ipt->cols));

  nframes = mrisp_template->Ip->num_frame;
  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
      *IMAGEFseq_pix(Ipt, u, v, fno) += *IMAGEFseq_pix(Ip, u, v, fno);
      if (u == DEBUG_U && v == DEBUG_V)
        fprintf(stderr, "x[%d][%d][%d] = %2.3f\n", DEBUG_U, DEBUG_V, fno, *IMAGEFseq_pix(Ipt, DEBUG_U, DEBUG_V, fno));
    }
  }

  return (mrisp_template);
}

int MRISPfunctionVectorVals(MRI_SURFACE_PARAMETERIZATION *mrisp,
                            MRIS *mris,
                            float x,
                            float y,
                            float z,
                            int *frames,
                            int nframes,
                            double *vals)
{
  double r, r2, r4, r6, val, d, x2, y2, z2, g, h, f, dx, dy, dz;
  float phi, theta, uf, vf, du, dv;
  int n, fno, u0, v0, u1, v1, u0_voff, u1_voff, u1_v0, u1_v1, u0_v0, u0_v1;

  r = sqrt(x * x + y * y + z * z);
  if (!FEQUAL(r, mris->radius)) /* project it onto sphere */
  {
    r = mris->radius;
    r2 = r * r;
    r4 = r2 * r2;
    r6 = r2 * r4;
    x2 = x * x;
    y2 = y * y;
    z2 = z * z;
    f = x2 / r6 + y2 / r6 + z2 / r6;
    g = 2 * (x2 / r4 + y2 / r4 + z2 / r4);
    h = x2 / r2 + y2 / r2 + z2 / r2 - 1;
    d = (-g + (float)sqrt((double)(g * g - 4 * f * h))) / (2 * f);
    dx = d * x / r2;
    dy = d * y / r2;
    dz = d * z / r2;
    x = x + dx;
    y = y + dy;
    z = z + dz;
  }

  theta = atan2(y / r, x / r);
  if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
  d = r * r - z * z;
  if (d < 0.0) d = 0.0;
  phi = atan2(sqrt(d), z);
  if (phi < RADIANS(1)) DiagBreak();

  uf = PHI_DIM(mrisp) * phi / PHI_MAX;
  vf = THETA_DIM(mrisp) * theta / THETA_MAX;
  u0 = floor(uf);
  u1 = ceil(uf);
  v0 = floor(vf);
  v1 = ceil(vf);
  du = uf - (float)u0;
  dv = vf - (float)v0;

  /* enforce spherical topology  */
  if (u0 < 0) /* enforce spherical topology  */
  {
    u0_voff = V_DIM(mrisp) / 2;
    u0 = -u0;
  }
  else if (u0 >= U_DIM(mrisp)) {
    u0_voff = V_DIM(mrisp) / 2;
    u0 = U_DIM(mrisp) - (u0 - U_DIM(mrisp) + 1);
  }
  else
    u0_voff = 0;
  if (u1 < 0) /* enforce spherical topology  */
  {
    u1_voff = V_DIM(mrisp) / 2;
    u1 = -u1;
  }
  else if (u1 >= U_DIM(mrisp)) {
    u1_voff = V_DIM(mrisp) / 2;
    u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
  }
  else
    u1_voff = 0;
  if (v0 < 0) v0 += V_DIM(mrisp);
  if (v0 >= V_DIM(mrisp)) v0 -= V_DIM(mrisp);
  if (v1 < 0) v1 += V_DIM(mrisp);
  if (v1 >= V_DIM(mrisp)) v1 -= V_DIM(mrisp);

  u0_v1 = v1 + u0_voff;
  while (u0_v1 >= V_DIM(mrisp)) u0_v1 -= V_DIM(mrisp);
  u0_v0 = v0 + u0_voff;
  while (u0_v0 >= V_DIM(mrisp)) u0_v0 -= V_DIM(mrisp);
  u1_v1 = v1 + u1_voff;
  while (u1_v1 >= V_DIM(mrisp)) u1_v1 -= V_DIM(mrisp);
  u1_v0 = v0 + u1_voff;
  while (u1_v0 >= V_DIM(mrisp)) u1_v0 -= V_DIM(mrisp);

  /* do bilinear interpolation */
  for (n = 0; n < nframes; n++) {
    fno = frames[n];
    val = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, u1_v1, fno) +
          (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, u0_v1, fno) +
          (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, u0_v0, fno) +
          du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, u1_v0, fno);
    vals[n] = val;
  }

  return (NO_ERROR);
}

MRIS *MRISfromParameterizations(MRI_SP *mrisp, MRIS *mris, int *frames, int *indices, int nframes)
{
  float a, b, c, phi, theta, x, y, z, uf, vf, du, dv, curv, d;
  int n, fno, vno, u0, v0, u1, v1;
  VALS_VP *vp;
  VERTEX *vertex;

  if (!mris) mris = MRISclone(mrisp->mris);

  MRIS *original = mris;
  mris = makeCenteredSphere(mris);

  a = b = c = MRISaverageRadius(mris);

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(vertex->y / b, vertex->x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (vno == 60935) DiagBreak();

    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u0 = floor(uf);
    u1 = ceil(uf);
    v0 = floor(vf);
    v1 = ceil(vf);
    du = uf - (float)u0;
    dv = vf - (float)v0;

    /* enforce spherical topology  */
    if (u0 < 0) /* enforce spherical topology  */
      u0 = -u0;
    if (u0 >= U_DIM(mrisp)) u0 = U_DIM(mrisp) - (u0 - U_DIM(mrisp) + 1);
    if (u1 < 0) /* enforce spherical topology  */
      u1 = -u1;
    if (u1 >= U_DIM(mrisp)) u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
    if (v0 < 0) v0 += V_DIM(mrisp);
    if (v0 >= V_DIM(mrisp)) v0 -= V_DIM(mrisp);
    if (v1 < 0) v1 += V_DIM(mrisp);
    if (v1 >= V_DIM(mrisp)) v1 -= V_DIM(mrisp);

    if (((u0 == DEBUG_U) || (u1 == DEBUG_U)) && ((v0 == DEBUG_V) || (v1 == DEBUG_V))) DiagBreak();
    vp = (VALS_VP *)vertex->vp;
    /* do bilinear interpolation */
    for (n = 0; n < nframes; n++) {
      fno = frames[n];
      curv = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) +
             (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, fno) +
             (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, fno) +
             du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, fno);
      vp->vals[indices[n]] = curv;
    }
  }

  resetCenteredSphere(original, mris);
  return original;
}

MRI_SP *MRIStoParameterizations(MRIS *mris, MRI_SP *mrisp, float scale, int *frames, int *indices, int nframes)
{
  MRIS *original = mris;
  mris = makeCenteredSphere(mris);

  float a, b, c, phi, theta, x, y, z, uf, vf, d, total_d, **distances, *total;
  int m, vno, u, v, unfilled, **filled, npasses, nfilled;
  VERTEX *vertex;
  VALS_VP *vp;

  total = (float *)malloc(nframes * sizeof(float));

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "computing parameterization...");

  if (!mrisp) {
    ErrorExit(0, "MRIStoParameterizations: NULL mrisp\n");
  }
  else
    for (m = 0; m < nframes; m++) ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, frames[m]);

  a = b = c = MRISaverageRadius(mris);

  filled = (int **)calloc(U_DIM(mrisp), sizeof(int *));
  distances = (float **)calloc(U_DIM(mrisp), sizeof(float *));
  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    filled[u] = (int *)calloc(V_DIM(mrisp), sizeof(int));
    distances[u] = (float *)calloc(V_DIM(mrisp), sizeof(float));

    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) filled[u][v] = UNFILLED_ELT;
  }

  //  fp = IMAGEFseq_pix(mrisp->Ip, DEBUG_U, DEBUG_V,fno) ;
  /* first calculate total distances to a point in parameter space */
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    if (vno == Gdiag_no || vno == 1126) DiagBreak();
    theta = atan2(y / b, x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0;
    phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (vno == DEBUG_VNO) DiagBreak();
    vertex->phi = phi;
    vertex->theta = theta;
    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u = nint(uf);
    v = nint(vf);
    if (u < 0) /* enforce spherical topology  */
      u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) /* enforce spherical topology  */
      v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    if (u == 0 && v == 56) DiagBreak();
    if ((u == DEBUG_U) && (v == DEBUG_V)) DiagBreak();

    filled[u][v] = vno;
    distances[u][v] += 1; /* keep track of total # of nodes */
    if ((u == DEBUG_U) && (v == DEBUG_V))
      fprintf(stderr,
              "v = %6.6d (%2.1f, %2.1f, %2.1f), d = %2.3f, "
              "curv = %2.3f\n",
              vno,
              x,
              y,
              z,
              d,
              vertex->curv);
  }

  if (DEBUG_U >= 0) fprintf(stderr, "\ndistance[%d][%d] = %2.3f\n\n", DEBUG_U, DEBUG_V, distances[DEBUG_U][DEBUG_V]);

  /* now add in curvatures proportional to their distance from the point */
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    theta = atan2(y / b, x / a);
    if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    d = c * c - z * z;
    if (d < 0.0) d = 0.0;
    phi = atan2(sqrt(d), z);
    uf = PHI_DIM(mrisp) * phi / PHI_MAX;
    vf = THETA_DIM(mrisp) * theta / THETA_MAX;
    u = nint(uf);
    v = nint(vf);
    if (u < 0) /* enforce spherical topology  */
      u = -u;
    if (u >= U_DIM(mrisp)) u = U_DIM(mrisp) - (u - U_DIM(mrisp) + 1);
    if (v < 0) /* enforce spherical topology  */
      v += V_DIM(mrisp);
    if (v >= V_DIM(mrisp)) v -= V_DIM(mrisp);

    if (u == 0 && v == 56) DiagBreak();

    /* 0,0 */
    total_d = distances[u][v];
    vp = (VALS_VP *)vertex->vp;
    if (total_d > 0.0) {
      for (m = 0; m < nframes; m++) {
        *IMAGEFseq_pix(mrisp->Ip, u, v, frames[m]) += vp->vals[indices[m]] / total_d;
      }
    }
  }

  /* fill in values which were unmapped using soap bubble */
  nfilled = npasses = 0;
  do {
    IMAGE *Ip, *Itmp;
    int u1, v1, uk, vk, n;

    Ip = mrisp->Ip;
    Itmp = ImageClone(Ip);
    ImageCopyFrames(Ip, Itmp, 0, Ip->num_frame, 0);
    unfilled = 0;
    ;
    for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
      for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
        if ((u == DEBUG_U) && (v == DEBUG_V)) DiagBreak();
        if (filled[u][v] == UNFILLED_ELT) {
          memset(total, 0, nframes * sizeof(float));
          for (n = 0, uk = -1; uk <= 1; uk++) {
            u1 = u + uk;
            if (u1 < 0) /* enforce spherical topology  */
              u1 = -u1;
            else if (u1 >= U_DIM(mrisp))
              u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
            for (vk = -1; vk <= 1; vk++) {
              v1 = v + vk;
              if (v1 < 0) /* enforce spherical topology  */
                v1 += V_DIM(mrisp);
              else if (v1 >= V_DIM(mrisp))
                v1 -= V_DIM(mrisp);

              if (filled[u1][v1] >= 0) {
                for (m = 0; m < nframes; m++) total[m] += *IMAGEFseq_pix(Ip, u1, v1, frames[m]);
                n++;
              }
            }
          }
          if (n > 0) {
            for (m = 0; m < nframes; m++) {
              total[m] /= (float)n;
              *IMAGEFseq_pix(Itmp, u, v, frames[m]) = total[m];
            }
            filled[u][v] = FILLING_ELT;
            nfilled++;
          }
          else
            unfilled++;
        }
        else {
          for (m = 0; m < nframes; m++) *IMAGEFseq_pix(Itmp, u, v, frames[m]) = *IMAGEFseq_pix(Ip, u, v, frames[m]);
        }
      }
    }
    for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
      for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
        if (filled[u][v] == FILLING_ELT) filled[u][v] = FILLED_ELT;
      }
    }
    mrisp->Ip = Itmp;
    ImageFree(&Ip);
    if (npasses++ > 1000)
      ErrorExit(ERROR_BADFILE,
                "MRISPtoParameterization: could not fill "
                "parameterization");
  } while (unfilled > 0);
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) fprintf(stderr, "filling %d elements took %d passes\n", nfilled, npasses);

  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    free(filled[u]);
    free(distances[u]);
  }
  free(filled);
  free(distances);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "done.\n");

  free(total);

  resetCenteredSphere(original, mris);

  return (mrisp);
}

MRI_SP *MRISPblurFrames(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma, int *frames, int nframes)
{
  int n, u, v, cart_klen, klen, khalf, uk, vk, u1, v1, no_sphere, voff;
  double k, *total, ktotal, sigma_sq_inv, udiff, vdiff, sin_sq_u, phi;
  IMAGE *Ip_src, *Ip_dst;

  total = (double *)malloc(nframes * sizeof(double));

  no_sphere = getenv("NO_SPHERE") != NULL;
  if (no_sphere) fprintf(stderr, "disabling spherical geometry\n");

  if (!mrisp_dst) mrisp_dst = MRISPclone(mrisp_src);
  mrisp_dst->sigma = sigma;

  /* determine the size of the kernel */
  cart_klen = (int)nint(6.0f * sigma) + 1;
  if (ISEVEN(cart_klen)) /* ensure it's odd */
    cart_klen++;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "blurring surface, sigma = %2.3f, cartesian klen = %d\n", sigma, cart_klen);

  if (FZERO(sigma))
    sigma_sq_inv = BIG;
  else
    sigma_sq_inv = 1.0f / (sigma * sigma);

  Ip_src = mrisp_src->Ip;
  Ip_dst = mrisp_dst->Ip;

  for (u = 0; u < U_DIM(mrisp_src); u++) {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "\r%3.3d of %d     ", u, U_DIM(mrisp_src) - 1);
    phi = (double)u * PHI_MAX / PHI_DIM(mrisp_src);
    sin_sq_u = sin(phi);
    sin_sq_u *= sin_sq_u;
    if (!FZERO(sin_sq_u)) {
      k = cart_klen * cart_klen;
      klen = sqrt(k + k / sin_sq_u);
      if (klen > MAX_LEN * cart_klen) klen = MAX_LEN * cart_klen;
    }
    else
      klen = MAX_LEN * cart_klen; /* arbitrary max length */
    if (no_sphere) sin_sq_u = 1.0f, klen = cart_klen;
    if (klen >= U_DIM(mrisp_src)) klen = U_DIM(mrisp_src) - 1;
    if (klen >= V_DIM(mrisp_src)) klen = V_DIM(mrisp_src) - 1;
    khalf = klen / 2;
    for (v = 0; v < V_DIM(mrisp_src); v++) {
      /*      theta = (double)v*THETA_MAX / THETA_DIM(mrisp_src) ;*/
      if (u == DEBUG_U && v == DEBUG_V) DiagBreak();

      memset(total, 0, nframes * sizeof(double));
      ktotal = 0.0;
      for (uk = -khalf; uk <= khalf; uk++) {
        udiff = (double)(uk * uk); /* distance squared in u */

        u1 = u + uk;
        if (u1 < 0) /* enforce spherical topology  */
        {
          voff = V_DIM(mrisp_src) / 2;
          u1 = -u1;
        }
        else if (u1 >= U_DIM(mrisp_src)) {
          u1 = U_DIM(mrisp_src) - (u1 - U_DIM(mrisp_src) + 1);
          voff = V_DIM(mrisp_src) / 2;
        }
        else
          voff = 0;

        for (vk = -khalf; vk <= khalf; vk++) {
          vdiff = (double)(vk * vk);
          k = exp(-(udiff + sin_sq_u * vdiff) * sigma_sq_inv);
          v1 = v + vk + voff;
          while (v1 < 0) /* enforce spherical topology */
            v1 += V_DIM(mrisp_src);
          while (v1 >= V_DIM(mrisp_src)) v1 -= V_DIM(mrisp_src);
          ktotal += k;
          for (n = 0; n < nframes; n++) total[n] += k * *IMAGEFseq_pix(Ip_src, u1, v1, frames[n]);
        }
      }
      if (u == DEBUG_U && v == DEBUG_V) DiagBreak();
      for (n = 0; n < nframes; n++) {
        total[n] /= ktotal; /* normalize weights to 1 */
        *IMAGEFseq_pix(mrisp_dst->Ip, u, v, frames[n]) = total[n];
      }
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "done.\n");

  free(total);

  return (mrisp_dst);
}

int MRISPsetFrameVal(MRI_SP *mrisp, int frame, float val)
{
  int u, v;

  for (u = 0; u <= U_MAX_INDEX(mrisp); u++) {
    for (v = 0; v <= V_MAX_INDEX(mrisp); v++) {
      *IMAGEFseq_pix(mrisp->Ip, u, v, frame) = val;
    }
  }
  return (NO_ERROR);
}
float MRISPsample(MRI_SP *mrisp, float x, float y, float z, int fno)
{
  float retval, theta, phi, radius, d, du, dv, uf, vf;
  int u0, u1, v0, v1;

  radius = mrisp->radius;
  theta = atan2(y / mrisp->radius, x / mrisp->radius);
  if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
  d = radius * radius - z * z;
  if (d < 0.0) d = 0.0;
  phi = atan2(sqrt(d), z);
  if (phi < RADIANS(1)) DiagBreak();
  if (phi > M_PI) DiagBreak();

  uf = PHI_DIM(mrisp) * phi / PHI_MAX;
  vf = THETA_DIM(mrisp) * theta / THETA_MAX;
  u0 = floor(uf);
  u1 = ceil(uf);
  v0 = floor(vf);
  v1 = ceil(vf);
  du = uf - (float)u0;
  dv = vf - (float)v0;

  /* enforce spherical topology  */
  if (u0 < 0) /* enforce spherical topology  */
    u0 = -u0;
  if (u0 >= U_DIM(mrisp)) u0 = U_DIM(mrisp) - (u0 - U_DIM(mrisp) + 1);
  if (u1 < 0) /* enforce spherical topology  */
    u1 = -u1;
  if (u1 >= U_DIM(mrisp)) u1 = U_DIM(mrisp) - (u1 - U_DIM(mrisp) + 1);
  if (v0 < 0) v0 += V_DIM(mrisp);
  if (v0 >= V_DIM(mrisp)) v0 -= V_DIM(mrisp);
  if (v1 < 0) v1 += V_DIM(mrisp);
  if (v1 >= V_DIM(mrisp)) v1 -= V_DIM(mrisp);

  retval = du * dv * *IMAGEFseq_pix(mrisp->Ip, u1, v1, fno) +
           (1.0f - du) * dv * *IMAGEFseq_pix(mrisp->Ip, u0, v1, fno) +
           (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, fno) +
           du * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, fno);

  return (retval);
}


/*
  Constructs a projector that associates a surface and a parameterization.
  The vertex->uv mapping is computed and cached here.
*/
SphericalProjector::SphericalProjector(MRIS *surf, MRI_SP *param) :
  mris(surf),
  mrisp(param)
{
  // cache input surface
  original = mris;
  mris = makeCenteredSphere(mris);
  interpolator = new SphericalInterpolator(mris);

  // image dimensions
  udim = U_DIM(mrisp);
  vdim = V_DIM(mrisp);

  // image maps
  nearest = ImageArray<int>(udim, vdim, -1);
  hits = ImageArray<int>(udim, vdim, 0);
  distance = ImageArray<float>(udim, vdim, 5);

  // vertex maps
  vertex_u = std::vector<int>(mris->nvertices);
  vertex_v = std::vector<int>(mris->nvertices);
  vertex_uf = std::vector<float>(mris->nvertices);
  vertex_vf = std::vector<float>(mris->nvertices);

  // first calculate total hits to a point in parameter space
  float radius = MRISaverageRadius(mris);

  for (int vno = 0; vno < mris->nvertices; vno++) {

    VERTEX *vertex = &mris->vertices[vno];
    float x = vertex->x;
    float y = vertex->y;
    float z = vertex->z;

    // translate xyz to spherical coordinates
    float theta = atan2(y / radius, x / radius);
    if (theta < 0.0f) theta = 2 * M_PI + theta;  // make it 0 -> 2PI
    float d = radius * radius - z * z;
    if (d < 0.0) d = 0;
    float phi = atan2(sqrt(d), z);

    // translate to image coordinates
    float uf = udim * phi / PHI_MAX;
    float vf = vdim * theta / THETA_MAX;
    int u = nint(uf);
    int v = nint(vf);

    // get distance to coordinate
    float du = uf - u;
    float dv = vf - v;
    float dist = std::sqrt(du * du + dv * dv);

    // enforce spherical topology
    if (u < 0) u = -u;
    if (u >= udim) u = udim - (u - udim + 1);
    if (v < 0) v += vdim;
    if (v >= vdim) v -= vdim;

    // cache vertex uv values
    vertex_uf[vno] = uf;
    vertex_vf[vno] = vf;
    vertex_u[vno] = u;
    vertex_v[vno] = v;

    // keep track of total # of vertices
    hits.item(u, v) += 1;

    // check if it's the closest point so far
    if (dist < distance.item(u, v)) {
      nearest.item(u, v) = vno;
      distance.item(u, v) = dist;
    }
  }
}


SphericalProjector::~SphericalProjector()
{
  resetCenteredSphere(original, mris);
  delete interpolator;
}


/*
  Projects an overlay array (of size nvertices) into the parameterization at the given frame index.
*/
void SphericalProjector::parameterizeOverlay(const float *overlay, int frameno, InterpMethod interp)
{
  // sample overlay values
  ImageClearArea(mrisp->Ip, -1, -1, -1, -1, 0, frameno);

  // barycentric sampling
  if (interp == Barycentric) {
    for (int vno = 0; vno < mris->nvertices; vno++) {
      int u = vertex_u[vno];
      int v = vertex_v[vno];
      *IMAGEFseq_pix(mrisp->Ip, u, v, frameno) += overlay[vno] / hits.item(u, v);
    }
  }
  // nearest neighbor sampling
  else if (interp == Nearest) {
    for (int u = 0; u < udim; u++)  {
      for (int v = 0; v < vdim; v++) {
        *IMAGEFseq_pix(mrisp->Ip, u, v, frameno) = overlay[nearest.item(u, v)];
      }
    }
  }

  // do backwards sampling to fill in missing pixels
  interpolator->setOverlay(overlay);
  interpolator->nearestneighbor = (interp == Nearest);

  for (int u = 0; u < udim; u++)  {
    for (int v = 0; v < vdim; v++) {
      if (hits.item(u, v) == 0) {
        double phi = u * PHI_MAX / udim;
        double theta = v * THETA_MAX / vdim;
        *IMAGEFseq_pix(mrisp->Ip, u, v, frameno) = interpolator->interp(phi, theta);
      }
    }
  }
}


/*
  Samples a parameterization (at the given frame index) into an overlay (of size nvertices).
*/
void SphericalProjector::sampleParameterization(float *overlay, int frameno, InterpMethod interp)
{
  // barycentric sampling
  if (interp == Barycentric) {

    for (int vno = 0; vno < mris->nvertices; vno++) {
      float uf = vertex_uf[vno];
      float vf = vertex_vf[vno];
      int u0 = floor(uf);
      int u1 = ceil(uf);
      int v0 = floor(vf);
      int v1 = ceil(vf);
      float du = uf - float(u0);
      float dv = vf - float(v0);

      // enforce spherical topology
      if (u0 < 0) u0 = -u0;
      if (u0 >= udim) u0 = udim - (u0 - udim + 1);
      if (u1 < 0) u1 = -u1;
      if (u1 >= udim) u1 = udim - (u1 - udim + 1);
      if (v0 < 0) v0 += vdim;
      if (v0 >= vdim) v0 -= vdim;
      if (v1 < 0) v1 += vdim;
      if (v1 >= vdim) v1 -= vdim;

      // interpolate
      overlay[vno] =  du  *         dv  * *IMAGEFseq_pix(mrisp->Ip, u1, v1, frameno) +
              (1.0f - du) *         dv  * *IMAGEFseq_pix(mrisp->Ip, u0, v1, frameno) +
              (1.0f - du) * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u0, v0, frameno) +
                      du  * (1.0f - dv) * *IMAGEFseq_pix(mrisp->Ip, u1, v0, frameno);
    }

  }
  // nearest neighbor sampling
  else if (interp == Nearest) {

    for (int vno = 0; vno < mris->nvertices; vno++) {
      int u = vertex_u[vno];
      int v = vertex_v[vno];
      overlay[vno] = *IMAGEFseq_pix(mrisp->Ip, u, v, frameno);
    }

  }
}
