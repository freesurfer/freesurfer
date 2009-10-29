/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *   {ch|rudolph}@nmr.mgh.harvard.edu                                      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
// $Id: c_vertex.cpp,v 1.2 2009/10/29 15:30:49 rudolph Exp $

#include "c_vertex.h"
#include "dijkstra.h"

bool
vertex_valLTE(
    VERTEX*     pvertex,
    void*       pv_void
) {
  return((pvertex->val >= 0) &&  (pvertex->val <= *((float*) pv_void)));
}

bool
vertex_alwaysTrue(
    VERTEX*     pvertex,
    void*       pv_void
) {
  return true;
}

bool
vertex_ripFlagIsTrue(
    VERTEX*     pvertex,
    void*       pv_void
) {
  return(pvertex->ripflag == TRUE);
}

void
vertex_ripFlagMark(
    VERTEX*     pvertex,
    void*       pv_mark
) {
  char* pch_mark = (char*) pv_mark;
  pvertex->ripflag = *pch_mark;
}

void
vertex_signumFunctional(
    VERTEX*     apvertex,
    void*       apv_functional
) {
  s_signumFunctional*  ps_signum  =
    (s_signumFunctional*)apv_functional;

  if (apvertex->curv <  0) ps_signum->negCount++;
  if (apvertex->curv >  0) ps_signum->posCount++;
  if (apvertex->curv == 0) ps_signum->zeroCount++;
}

void
vertex_rawCurveSum(
    VERTEX*     apvertex,
    void*       apv_functional
) {
  s_rawCurve*   ps_rawC =
    (s_rawCurve*)apv_functional;

  if (apvertex->curv <  0) ps_rawC->f_negCount+=apvertex->curv;
  if (apvertex->curv >  0) ps_rawC->f_posCount+=apvertex->curv;
}

void
vertex_rawCurveMinMax(
    VERTEX*     apvertex,
    void*       apv_functional
) {
  //
  // PRECONDITIONS
  //  o Minimum values are assumed to be less than zero,
  //    maxmimum values greater than zero.
  //  o A call to surface_annotation_do() is required to
  //   correctly set up the vertex indices.
  //

  s_minMax*   ps_minMax =
    (s_minMax*)apv_functional;

  if (apvertex->curv <  ps_minMax->f_min) {
    ps_minMax->f_min  = apvertex->curv;
    ps_minMax->minVertex = apvertex->annotation;
  }
  if (apvertex->curv >  ps_minMax->f_max) {
    ps_minMax->f_max  = apvertex->curv;
    ps_minMax->maxVertex = apvertex->annotation;
  }
}

void
vertex_curveAreaSum(
    VERTEX*     apvertex,
    void*       apv_functional
) {
  s_integratedCurve*  ps_totalIC =
    (s_integratedCurve*)apv_functional;

  float f_product = 0.0;

  switch (ps_totalIC->e_curvature) {
  case e_gaussian:
    f_product = fabs(apvertex->K) * apvertex->area;
    break;
  case e_mean:
    f_product = fabs(apvertex->H) * apvertex->area;
    break;
  }

  ps_totalIC->f_sum += f_product;
}

bool
vertex_isDijkVirgin(
    VERTEX*     pvertex,
    void*       pv_void
) {
  return(pvertex->marked == DIJK_VIRGIN);
}

void
vertex_annotationMark(
    VERTEX*     pvertex,
    void*       pv_mark
) {
  int* mark = (int*) pv_mark;
  pvertex->annotation = *mark;
}

/* eof */
