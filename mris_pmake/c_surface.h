/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *    Center for Morphometric Analysis                                     *
 *    Massachusetts General Hospital                                       *
 *    Building 149, 13th St.                                               *
 *    Charlestown, MA 02129                                                *
 *     {ch|rudolph}@nmr.mgh.harvard.edu                                    *
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
/// \file c_surface.h
///
/// \brief Brief description
/// The surface related object API.
///
/// \b DESCRIPTION
/// Surface type functions setting whole surface function pointers, and
/// function definitions.
///
/// \b HISTORY
/// 15 March 2005 - Initial consolidation from several other sources.
/// $Id: c_surface.h,v 1.2 2009/10/29 15:30:49 rudolph Exp $
///
///

#ifndef __C_SURFACE_H__
#define __C_SURFACE_H__

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

void
surface_vertexFunction_do(
  s_env&   st_env,
  bool   (*vertex_satisfyTestCondition)
  (VERTEX* pvertex,
   void*  pv_extra),
  void*   apv_conditional,
  void   (*vertex_function)
  (VERTEX* pvertex,
   void*  pv_extra),
  void*   apv_functional
);

/// \fn void surface_vertexPatternCopy(s_env& st_env, MRIS* apMS_source, MRIS* apMS_target, bool (*vertex_satisfyTestCondition) (VERTEX* pvertex, void* pv_extra), void* apvsource_extra, void (*vertex_modify) (VERTEX* pvertex, void* pv_extra), void*apvtarget_extra)
/// \fn void surface_vertexPatternCopy(s_env& st_env, MRIS* apMS_source, MRIS* apMS_target, bool (*vertex_satisfyTestCondition) (VERTEX* pvertex, void* pv_extra), void* apvsource_extra, void (*vertex_modify) (VERTEX* pvertex, void* pv_extra), void*apvtarget_extra)
/*/// \fn label_coreLoad */
/// \brief Copies a vertex "pattern" from one surface to another. The pattern
///  is described by a TRUE return on each vertex of the source
///  surface subject to (*vertex_satisfyTestCondition)(). On
///  the target surface corresponding vertices are processed with
///  (*vertex_modify)().
/// \param st_env   the main program environment structure
/// \param apMS_source  source surface containing vertex pattern
/// \param apMS_target  target surface to contain vertex pattern
/// \param (*vertex_satisfyTestCondition) per-vertex conditional test
/// \param apvsource_extra  a pointer to void that embodies the
///     argument list to
///     (*vertex_satisfyTestCondition)()
/// \param (*vertex_modify)  per-vertex modifier
/// \param apvtarget_extra  a pointer to void that embodies the
///     argument list to
///     (*vertex_modify)()
/// \return    (void)
void
surface_vertexPatternCopy(
    s_env&      st_env,
    MRIS*       apMS_source,
    MRIS*       apMS_target,
    bool        (*vertex_satisfyTestCondition)
    (VERTEX*    pvertex,
        void*   pv_extra),
    void*       apvsource_extra,
    void        (*vertex_modify)
    (VERTEX*    pvertex,
        void*   pv_extra),
    void*       apvtarget_extra
);

/// \fn void surface_signumFunction_do(s_env& st_env)
/// \brief The sign of curvature values on the internal active surface
///  is counted and dumped to stdout.
/// \param  st_env  the problem environment
/// \return (void)
void
surface_signumFunction_do(
  s_env&  st_env
);

/// \fn void surface_correlationFunction_do(s_env& st_env)
/// \brief A correlation between the working and auxillary
///  surfaces is determined. The auxillary surface
///  is assumed to the be "reference" surface.
/// \param  st_env  the problem environment
/// \return (void)
void
surface_correlationFunction_do(
  s_env&  st_env
);

/// \fn void surface_rawCurveSum_do(s_env& st_env)
/// \brief Calculates the positive and negative summation of
///  signed curvature values for the active surface.
/// \param  st_env  the problem environment
/// \return (void)
void
surface_rawCurveSum_do(
  s_env&  st_env
);

/// \fn void surface_rawCurveMinMax_do(s_env& st_env)
/// \brief Determines the minimum and maximum curvature values
///  for a given surface.
/// \param  st_env  the problem environment
/// \return (void)
void
surface_rawCurveMinMax_do(
  s_env&  st_env
);

/// \fn void surface_averageIntegratedCurvArea_do(s_env& st_env, e_CURVATURE ae_curvature)
/// \brief Calculates the average integrated curvature measure for a given curvature
///  type. The calculation is merely 1/N*(curvType*curvArea) where curvArea
///  is the discrete "limit" area measure of a vertex, taken over all the
///  vertices in a surface.
/// \param  st_env  the problem environment
/// \param ae_curvature the curvature type to integrate (Gaussian, mean, etc.)
/// \return (void)
void
surface_averageIntegratedCurveArea_do(
  s_env&  st_env,
  e_CURVATURE ae_curvature = e_gaussian
);

/// \fn void surface_workingToAux_ripTrueCopy(s_env& st_env)
/// \brief A pattern of TRUE rips across the vertices on the internal working
///  surface is copied to the internal auxillary surface.
/// \param  st_env  the problem environment
/// \return (void)
void
surface_workingToAux_ripTrueCopy(
  s_env&  st_env
);

/// \fn void surface_annotate(s_env& st_env)
/// \brief This moves in vertex order over a surface and sets each vertex's
///  'annotation' field to its cardinal loop order.
/// \param  st_env  the problem environment
/// \return (void)
void
surface_annotation_do(
  s_env&   st_env
);

/// \fn void surface_ripMark(s_env& st_env)
/// \brief This traces an already solved dijkstra path on a surface and sets
///  each member vertex's ripflag to TRUE
/// \param  st_env  the problem environment
/// \return (void)
void
surface_ripMark(
  s_env&   st_env
);

/// \fn void surface_ripClear(s_env& st_env, bool b_wholeSurfaceForce)
/// \brief Clears the 'ripflag' member of vertices across the surface.
///  By default, only vertices that are part of a calculated
///  Dijkstra path are cleared. By setting b_wholeSurfaceForce to
///  'true', \b each node is cleared.
/// \param  st_env  the problem environment
/// \param b_wholeSurfaceForce if 'true', examine each node and set 'ripflag' to 'FALSE'
/// \return (void)
void
surface_ripClear(
  s_env&   st_env,
  bool   b_wholeSurfaceForce = false
);


#endif //__C_SURFACE_H__


