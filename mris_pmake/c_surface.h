/**
 * @brief The surface related object API.
 *
 * Surface type functions setting whole surface function pointers, and
 * function definitions.
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


#ifndef __C_SURFACE_H__
#define __C_SURFACE_H__

#include "general.h"
#include "env.h"

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"

#include <string>
using namespace std;

void
surface_vertexFunction_do(
  s_env&        st_env,
  bool          (*vertex_satisfyTestCondition)
                (VERTEX*        pvertex,
                void*           pv_extra),
  void*         apv_conditional,
  void          (*vertex_function)
                (VERTEX*        pvertex,
                void*           pv_extra),
  void*         apv_functional
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
surface_primaryToSecondary_ripTrueCopy(
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
/// \return the total cost of traveling along the marked path
float
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
  s_env&        st_env,
  bool          b_wholeSurfaceForce = false
);


#endif //__C_SURFACE_H__


