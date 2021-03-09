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


#ifndef __C_VERTEX_H__
#define __C_VERTEX_H__

#include "general.h"
#include "env.h"

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"

#include <string>
using namespace std;

/// This structure houses information pertaining to the signum count
/// of curvature values for a particular surface.
typedef struct _signumFunctional {
  int  posCount;
  int  negCount;
  int  zeroCount;
}
s_signumFunctional;

/// This structure houses information pertaining to the raw count
/// of curvature values for a particular surface.
typedef struct _rawCurve {
  float f_posCount;
  float f_negCount;
}
s_rawCurve;

/// This structure houses information pertaining to the raw count
/// of curvature minimum and maximum values for a particular surface.
typedef struct _minMax {
  float f_min;
  int  minVertex;
  float f_max;
  int  maxVertex;
}
s_minMax;

/// This structure houses information pertaining to the integrated curvature
/// for a particular surface.
typedef struct _integratedCurve {
  e_CURVATURE  e_curvature;
  float  f_sum;
}
s_integratedCurve;


/// \fn bool vertex_valLTE(VERTEX* pvertex, void* pv_void)
/// \brief A function (typically used as a callback) that determines if the
///  value parameter of a passed vertex is less than or
///  equal to a float unpack of pv_void*. Implied also is that
///  the value is greater than zero.
/// \param  pvertex The vertex to process
/// \param pv_void A pointer to a void structure. This is unpacked
///   (if necessary) by the function itself.
/// \return    (bool) TRUE if vertex value parameter is LTE
///    pv_void and greater than zero.
bool
vertex_valLTE(
  VERTEX*   pvertex,
  void*   pv_void
);

/// \fn bool vertex_isDijkVirgin(VERTEX* pvertex, void* pv_void)
/// \brief A function (typically used as a callback) that determines if the passed vertex is marked as DIJK_VIRGIN.
/// \param  pvertex The vertex to process
/// \param pv_void A pointer to a void structure. This is unpacked
///   (if necessary) by the function itself.
/// \return    (bool) TRUE if vertex is marked with DIJK_VIRGIN, else FALSE
bool
vertex_isDijkVirgin(
  VERTEX*   pvertex,
  void*   pv_void
);

/// \fn bool vertex_ripFlagIsTrue(VERTEX* pvertex, void* pv_void)
/// \brief A function (typically used as a callback) that determines if the "ripflag" for the given vertex index is TRUE or FALSE.
/// \param  pvertex The vertex to process
/// \param pv_void A pointer to a void structure. This is unpacked
///   (if necessary) by the function itself.
/// \return    (bool) TRUE if ripflag is TRUE, else FALSE
bool
vertex_ripFlagIsTrue(
  VERTEX*   pvertex,
  void*   pv_void
);

/// \fn bool vertex_alwaysTrue(VERTEX* pvertex, void* pv_void)
/// \brief A function (typically used as a callback) that always returns true
///  on its vertex.
/// \param  pvertex The vertex to process
/// \param pv_void A pointer to a void structure. This is unpacked
///   (if necessary) by the function itself.
/// \return    (bool) TRUE (always)
bool
vertex_alwaysTrue(
  VERTEX*   pvertex,
  void*   pv_void
);

/// \fn void vertex_signumFunctional(VERTEX* pvertex, void* apv_functional)
/// \brief A function (typically used as a callback) that counts the signum of the curvature value
/// \param  pvertex The vertex to process
/// \param apv_functional A pointer to void* housing the arguments that
///    the function internals require
/// \return    (void)
void
vertex_signumFunctional(
  VERTEX*   apvertex,
  void*   apv_functional
);

/// \fn void vertex_curvAreaSum(VERTEX* pvertex, void* apv_functional)
/// \brief A function (typically used as a callback) that determines the
///  product of a particular curvature and area at a vertex.
/// \param  pvertex The vertex to process
/// \param apv_functional A pointer to void* housing the arguments that
///    the function internals require
/// \return    (void)
void
vertex_curveAreaSum(
  VERTEX*   apvertex,
  void*   apv_functional
);

/// \fn void vertex_rawCurveSum(VERTEX* pvertex, void* apv_functional)
/// \brief A function (typically used as a callback) that determines the
///  sum of positive and negative curvature values at a vertex.
/// \param  pvertex The vertex to process
/// \param apv_functional A pointer to void* housing the arguments that
///    the function internals require
/// \return    (void)
void
vertex_rawCurveSum(
  VERTEX*   apvertex,
  void*   apv_functional
);

/// \fn void vertex_rawCurveMinMax(VERTEX* pvertex, void* apv_functional)
/// \brief A function (typically used as a callback) that determines the
///  minimum and maximum curvature values on a surface by examining
///  each vertex on the surface.
/// \param  pvertex The vertex to process
/// \param apv_functional A pointer to void* housing the arguments that
///    the function internals require
/// \return    (void)
void
vertex_rawCurveMinMax(
  VERTEX*   apvertex,
  void*   apv_functional
);

/// \fn void vertex_ripFlagMark(VERTEX* pvertex, void* pv_mark)
/// \brief A function (typically used as a callback) that marks the "ripflag" for the given vertex index as 'ch_mark'
/// \param  pvertex The vertex to process
/// \param  ch_mark The mark to write to the 'ripflag' field
/// \return    (void)
void
vertex_ripFlagMark(
  VERTEX*   pvertex,
  void*   pv_mark
);

/// \fn void vertex_annotationMark(VERTEX* pvertex, void* pv_mark)
/// \brief A function (typically used as a callback) that marks the
/// "annotation" for the given vertex index as (int) *pv_mark
/// \param  pvertex The vertex to process
/// \param  ch_mark The mark to write to the 'annotation' field
/// \return    (void)
void
vertex_annotationMark(
  VERTEX*   pvertex,
  void*   pv_mark
);

#endif //__C_VERTEX_H__


