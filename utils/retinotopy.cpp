/**
 * @brief Utilities for retinotopy analysis.
 *
 */
/*
 * Original Author: Doug Greve (and Marty and Anders, for now)
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
#include <stdlib.h>

double round(double x);
#include "fsglm.h"
#include "mri.h"
#include "mrisurf.h"
#include "surfcluster.h"

/*------------------------------------------------------------
  \fn void RETcompute_angles(MRIS *mris)
  Upon input:
    val  = real * log10(p) for eccen
    val2 = imag * log10(p) for eccen
    valbak  = real * log10(p) for polar
    val2bak = imag * log10(p) for polar
  Output:
    val  = angle for eccen (0 to 2*pi)
    val2 = log10(p) for eccen
    valbak  = angle for polar (-pi to +pi)
    val2bak = log10(p) for polar
  ------------------------------------------------------------*/
void RETcompute_angles(MRIS *mris, double EccenRotAngleRad, double PolarRotAngleRad)
{
  int k;
  float val, valbak, val2, val2bak;

  for (k = 0; k < mris->nvertices; k++) {
    if (!mris->vertices[k].ripflag) {
      // Eccen: 0-2pi
      val = atan2(mris->vertices[k].val2, mris->vertices[k].val);
      val += EccenRotAngleRad;
      if (val < 0.0) val += (2 * M_PI);
      val2 = sqrt(SQR(mris->vertices[k].val2) + SQR(mris->vertices[k].val));
      // Polar
      valbak = atan2(mris->vertices[k].val2bak, mris->vertices[k].valbak);
      valbak += PolarRotAngleRad;
      val2bak = sqrt(SQR(mris->vertices[k].val2bak) + SQR(mris->vertices[k].valbak));
    }
    else {
      val = 0.0;
      val2 = 0.0;
      valbak = 0.0;
      val2bak = 0.0;
    }
    mris->vertices[k].val = val;
    mris->vertices[k].val2 = val2;
    mris->vertices[k].valbak = valbak;
    mris->vertices[k].val2bak = val2bak;
  }
}
/*------------------------------------------------------------*/
float RETcircsubtract(float a, float b)
{
  float h = a - b;
  if (h < -M_PI)
    h = h + 2 * M_PI;
  else if (h > M_PI)
    h = h - 2 * M_PI;
  return h;
}
/*------------------------------------------------------------*/
void RETcompute_fieldsign(MRIS *mris)
{
  int k, m, n;
  float dv1, dv2, dx, dy, dv1dx, dv1dy, dv2dx, dv2dy;
  float m11, m12, m13, m22, m23, z1, z2, z3, z1b, z2b, z3b, denom;

  MRISremoveTriangleLinks(mris);
  printf("surfer: compute_fieldsign()\n");
  for (k = 0; k < mris->nvertices; k++) {
    if (!mris->vertices[k].ripflag) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
      VERTEX                * const v  = &mris->vertices         [k];
      dv1dx = dv1dy = dv2dx = dv2dy = 0;
      m11 = m12 = m13 = m22 = m23 = z1 = z2 = z3 = z1b = z2b = z3b = 0;
      n = 0;
      for (m = 0; m < vt->vnum; m++) {
        if (mris->vertices[vt->v[m]].ripflag) continue;
        dv1 = RETcircsubtract(v->val,    mris->vertices[vt->v[m]].val);
        dv2 = RETcircsubtract(v->valbak, mris->vertices[vt->v[m]].valbak);
        if (dv1 == 0 || dv2 == 0) continue;
        dx = v->x - mris->vertices[vt->v[m]].x;
        dy = v->y - mris->vertices[vt->v[m]].y;
        m11 += dx * dx;
        m12 += dx * dy;
        m13 += dx;
        m22 += dy * dy;
        m23 += dy;
        z1 += dx * dv1;
        z2 += dy * dv1;
        z3 += dv1;
        z1b += dx * dv2;
        z2b += dy * dv2;
        z3b += dv2;
        n++;
      }  // loop over neighbors
      dv1dx = (m22 * z1 - m23 * m23 * z1 - m12 * z2 + m13 * m23 * z2 - m13 * m22 * z3 + m12 * m23 * z3);
      dv2dx = (m22 * z1b - m23 * m23 * z1b - m12 * z2b + m13 * m23 * z2b - m13 * m22 * z3b + m12 * m23 * z3b);
      dv1dy = (-m12 * z1 + m13 * m23 * z1 + m11 * z2 - m13 * m13 * z2 + m12 * m13 * z3 - m11 * m23 * z3);
      dv2dy = (-m12 * z1b + m13 * m23 * z1b + m11 * z2b - m13 * m13 * z2b + m12 * m13 * z3b - m11 * m23 * z3b);
      denom = -m12 * m12 + m11 * m22 - m13 * m13 * m22 + 2 * m12 * m13 * m23 - m11 * m23 * m23;
      if (denom != 0)
        v->fieldsign = (dv1dx * dv2dy - dv2dx * dv1dy) / (denom * denom);
      else
        v->fieldsign = 0;

      v->fieldsign = ((v->fieldsign < 0) ? -1 : (v->fieldsign > 0) ? 2 : 0);

      v->fsmask = sqrt(v->val2 * v->val2bak); /* geom mean of r,th power */
    }                                         // if
  }                                           // loop over vertices
}

/*--------------------------------------------------------------------------*/
int RETlogMap(MRIS *surf, double k, double a, double xc0, double yc0)
{
  double n, r, theta, xw, yw, r2, theta2, xc2, yc2, xc, yc;
  float dmin;
  int vno;

  xc0 = surf->vertices[16930].x;
  yc0 = surf->vertices[16930].y;

  for (n = 1; n < 5; n += .3) {
    r = exp(n);
    for (theta = -M_PI / 2.0; theta < M_PI / 2.0; theta += .2 / r) {
      xw = r * cos(theta);
      yw = r * sin(theta);
      r2 = sqrt((xw + a) * (xw + a) + yw * yw);
      theta2 = atan2(yw, xw + a);
      xc2 = k * log(r2 * cos(theta2));
      yc2 = k * log(r2 * sin(theta2));
      xc = (xc2 - xc0);
      yc = (yc2 - yc0);
      vno = MRISfindClosestVertex(surf, xc, yc, 0, &dmin, CURRENT_VERTICES);
      if (dmin > 5) continue;
      surf->vertices[vno].val = 1;
      surf->vertices[vno].val2 = dmin;
      printf("%2f %6.1f %6.4f   %6.2f %6.2f\n", n, r, theta, xc, yc);
    }
  }
  for (n = 1; n < 10; n += .1) {
    r = exp(n);
    for (theta = -M_PI / 2.0; theta < M_PI / 2.0; theta += .6) {
      xw = r * cos(theta);
      yw = r * sin(theta);
      r2 = sqrt((xw + a) * (xw + a) + yw * yw);
      theta2 = atan2(yw, xw + a);
      xc2 = k * log(r2 * cos(theta2));
      yc2 = k * log(r2 * sin(theta2));
      xc = (xc2 - xc0);
      yc = (yc2 - yc0);
      vno = MRISfindClosestVertex(surf, xc, yc, 0, &dmin, CURRENT_VERTICES);
      if (dmin > 5) continue;
      surf->vertices[vno].val = -1;
      surf->vertices[vno].val2 = dmin;
      printf("%2f %6.1f %6.4f   %6.2f %6.2f\n", n, r, theta, xc, yc);
    }
  }
  return (0);
}

/*--------------------------------------------------------------------------*/
int RETinvLogMapFunc(double xc, double yc, double xc0, double yc0, double a, double k, double *r, double *theta)
{
  double r2, theta2, xw, yw;
  r2 = exp((xc - xc0) / k);
  theta2 = (yc - yc0) / k;
  xw = r2 * cos(theta2) - a;
  yw = r2 * sin(theta2);
  *r = sqrt(xw * xw + yw * yw);
  *theta = atan2(yw, xw);
  return (0);
}
/*--------------------------------------------------------------------------*/
int RETreverseSign(MRI *mri)
{
  int c, r, s;
  double v;

  for (c = 0; c < mri->width; c++) {
    for (r = 0; r < mri->height; r++) {
      for (s = 0; s < mri->depth; s++) {
        v = MRIgetVoxVal(mri, c, r, s, 0);
        if (v > 0)
          v = -1;
        else if (v < 0)
          v = +2;
        MRIsetVoxVal(mri, c, r, s, 0, v);
      }
    }
  }
  return (0);
}
