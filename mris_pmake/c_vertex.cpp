/**
 * @brief The vertex related object API.
 *
 * Vertex type functions include conditional checks as well as functional
 * modifications on vertex fields.
 */
/*
 * Original Author: Rudolph Pienaar / Christian Haselgrove
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
