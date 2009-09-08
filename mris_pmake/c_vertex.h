/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *    Center for Morphometric Analysis       *
 *    Massachusetts General Hospital        *
 * Building 149, 13th St.         *
 *  Charlestown, MA 02129         *
 *     {ch|rudolph}@nmr.mgh.harvard.edu      *
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
/// \file c_vertex.h
///
/// \brief Brief description
/// The vertex related object API.
///
/// \b DESCRIPTION
/// Vertex type functions include conditional checks as well as functional
/// modifications on vertex fields.
///
/// \b HISTORY
/// 15 March 2005 - Initial consolidation from several other sources.
/// $Id: c_vertex.h,v 1.1 2009/09/08 22:39:27 nicks Exp $
///
///

#ifndef __C_VERTEX_H__
#define __C_VERTEX_H__

#include "general.h"
#include "env.h"

#ifdef __cplusplus
extern  "C" {
#endif

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"

#ifdef __cplusplus
}
#endif

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


