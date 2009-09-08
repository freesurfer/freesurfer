/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *    Center for Morphometric Analysis       *
 *    Massachusetts General Hospital        *
 *    Building 149, 13th St.         *
 *    Charlestown, MA 02129         *
 *    {ch|rudolph}@nmr.mgh.harvard.edu      *
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
/// \file env.h
///
/// \brief Brief description
/// The environment object API.
///
/// \b DESCRIPTION
/// The environment acts as the main access / interface between the core
/// components of the program and the external system. Several structures
/// and objects, taken together, constitute this environment and define
/// amongst others, operational parameters, cost functions, weight
/// structures, etc.
///
/// \b HISTORY
/// 08 March 2005 - Initial consolidation from several other sources.
/// $Id: env.h,v 1.1 2009/09/08 22:39:27 nicks Exp $
///
///

#ifndef __ENV_H__
#define __ENV_H__

#include "general.h"
#include "c_SMessage.h"
#include "scanopt.h"

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



/// Weights and values for the main cost function polynomial.
typedef struct _weights {
  float wd;   // distance
  float wc;   // curve
  float wh;   // sulcal height
  float wdc;   // distance X curve
  float wdh;   // distance X height
  float wch;   // curve X height
  float wdch;   // distance X curve X height
  float wdir;   // direction vector
}
s_weights;

void
s_weights_scan(
  s_weights&  st_costWeight,
  C_scanopt&  cso_options
);

/// \fn void s_weights_print(s_weights& asw)
/// \brief Print the current internal weight structure to stdout.
/// \param  asw The weight structure to print
/// \return    (void)
void s_weights_print( s_weights& asw);

/// \fn void s_weights_setAll(s_weights& asw, float af)
/// \brief Sets all the weights in the current internal weight structure to the passed value.
/// \param  asw The weight structure to access
/// \param af Value assigned to each member
/// \return    (void)
void  s_weights_setAll( s_weights& asw,
                        float  af);
/// Del weights and values.
///
/// These are used when zero-crossings
/// of curvature and sulcal height occur. By incorporating these
/// transition weights, we can force the path search algorithm
/// to choose to strongly avoid areas of zero-crossing.
typedef struct _Dweights {
  float Dwd;   // distance
  float Dwc;   // curve
  float Dwh;   // sulcal height
  float Dwdc;   // distance X curve
  float Dwdh;   // distance X height
  float Dwch;   // curve X height
  float Dwdch;   // distance X curve X height
  float Dwdir;   // direction vector
}
s_Dweights;

void
s_Dweights_scan(
  s_Dweights&  st_DcostWeight,
  C_scanopt&  cso_options
);

/// \fn void s_Dweights_print(s_Dweights& asw)
/// \brief Print the current internal Dweight structure to stdout.
/// \param  asw The weight structure to print
/// \return    (void)
void s_Dweights_print( s_Dweights& asw);

/// \fn void s_Dweights_setAll(s_Dweights& asw, float af)
/// \brief Sets all the delta weights in the current internal Dweight structure to the passed value.
/// \param  asw The Dweight structure to access
/// \param af Value assigned to each member
/// \return    (void)
void  s_Dweights_setAll( s_Dweights& asw,
                         float  af);
/// s_iterInfo contains information pertaining to a particular iteration
/// of the main program loop. It is accessed during each call to the
/// cost function, and is populated with per-call information.
///
typedef struct _iterInfo {
  int  iter;
  float  f_distance;
  float  f_curvature;
  float  f_sulcalHeight;
  float  f_dir;
}
s_iterInfo;

/// The main environment structure. This structure records important
/// variables, mostly interpreted from the process <b>options</b> file.
typedef struct _env s_env;

typedef enum {
  e_default, e_unity, e_euclid, e_distance
} e_COSTFUNCTION;

typedef enum {
  e_workingCurvature, e_workingSulcal, e_auxillary
} e_SURFACE;

typedef enum {
  e_user, e_sys, e_result
} e_LOG;

typedef struct _env {
  s_weights* pSTw;   // weight structure
  s_Dweights* pSTDw;   // Del weight structure

  int  timeoutSec;  // listen timeout
  int  port;   // port on which to listen
  // for async control

  int  lw;   // left width (for std format)
  int  rw;   // right width (for std format)

  bool  b_syslogPrepend; // prepend syslog style
  C_SMessage* pcsm_syslog;  // log file for "sys" events
  C_SMessage* pcsm_userlog;  // log file for "user" events
  C_SMessage* pcsm_resultlog;  // log file for "result" event

  bool  b_labelFile_save; // flag: save label file
  bool  b_patchFile_save; // flag: save patch file
  bool  b_transitionPenalties; // flag: apply transition
  // penalties

  string  str_workingDir;  // directory containing input
  // and output files
  string  str_patchFileName; // file to contain path "patch"
  string  str_labelFileName; // file to contain path "label"
  string  str_labelFileNameOS; // file to contain path "patch"
  // projected back onto
  // the "Original Surface"
  string  str_optionsFileName; // file containing meta options
  string  str_costFileName; // file containing final costs

  int  startVertex;  // index of the start vertex
  int  endVertex;  // index of the end vertex
  float  plyDepth;  // ply depth around a core path

  e_SURFACE esf_active;  // the current active surface
  string*  pstr_activeName; // names of each surface
  int  totalNumSurfaces; // total number of surfaces
  MRIS*  pMS_active;  // a pointer (used by ply
  // functions) to specifiy
  // a particular surface to
  // process
  MRIS*  pMS_curvature;  // (inflated) curvature surface
  MRIS*  pMS_sulcal;  // (inflated) sulcal height
  // surface
  MRIS*  pMS_auxSurface;  // auxillary (optional)
  // surface

  bool  b_surfacesKeepInSync; // flag: behavioural /
  // conditional. Evaluate
  // if wanting to merge
  // information from one
  // surface to another.

  bool  b_surfacesClear; // flag: behavioural /
  // conditional. Should be
  // true for most cases.
  // If false, prevents the
  // clearing of rips
  // after a path search.
  // Useful when keeping
  // a pattern for cases
  // when correlation needs
  // to be calculated.

  bool  b_costHistoryPreserve; // flag: preserve cost history
  // on surface between
  // successive calls to
  // dijkstra function.
  // Set to TRUE if finding
  // ply distances from
  // existing path
  int  totalNumFunctions; // total number of cost
  // functions
  e_COSTFUNCTION ecf_current;  // the current cost function
  string*  pstr_functionName; // names of each cost function
  float  (*costFunc_do)  // a cost function to operate
  // on this environment
  (
    s_env&  st_env,
    s_iterInfo* pst_iterInfo,
    int   vno_c,
    int   j,
    bool  b_relNextReference
  );
}
s_env;

void
s_env_scan(
  s_env&   st_env,
  C_scanopt&  cso_options
);

/// \fn void s_env_nullify( s_env& st_env);
/// \brief Initialise (nullify, i.e. primitive constructor) a passed environment.
/// \param  st_env The environment structure
/// \return    (void)
void s_env_nullify( s_env& st_env);

bool s_env_b_surfacesKeepInSync_get(
  s_env&   ast_env);

void s_env_b_surfacesKeepInSync_set(
  s_env&   ast_env,
  bool   b_val);

bool s_env_b_surfacesClear_get(
  s_env&   ast_env);

void s_env_b_surfacesClear_set(
  s_env&   ast_env,
  bool   b_val);

bool s_env_surfaceFile_set(
  s_env&   st_env,
  string   astr_fileName
);

bool s_env_surfaceCurvature_set(
  s_env&   st_env,
  string   astr_fileName
);

bool s_env_surfaceSulcal_set(
  s_env&   st_env,
  string   astr_fileName
);

bool s_env_auxSurfaceFile_set(
  s_env&   st_env,
  string   astr_fileName
);

bool s_env_auxSurfaceCurvature_set(
  s_env&   st_env,
  string   astr_fileName
);

void  s_env_log_file_changeTo(
  s_env&   ast_env,
  e_LOG   ae_log,
  string   astr_newName
);

void s_env_activeSurfaceList(
  s_env&   ast_env
);

void s_env_activeSurfaceSetIndex(
  s_env*   apst_env,
  int   aindex
);


void s_env_costFctList(
  s_env&   ast_env
);

void s_env_costFctSetIndex(
  s_env*   apst_env,
  int   aindex
);

void  s_env_costFctSet(
  s_env*   pst_env,
  float (*cost_fct)
  (
    s_env&  st_env,
    s_iterInfo* pst_iterInfo,
    int   vno_c,
    int   j,
    bool  b_relNextReference
  ),
  e_COSTFUNCTION  aecf_new = e_default
);

/// \fn float costFunc_defaultDetermine(s_env& st_env, s_iterInfo* pst_iterInfo, int vno_c, int j, bool b_relNextReference = true);
/// \brief Calculate the cost in moving from one vertex to another.
/// \param  st_env The environment structure
/// \param pst_iterInfo A structure that records vertex transition information (cost, distance, curvature, etc.)
/// \param  vno The current vertex index
/// \param  j The neighbour vertex index - either relative or absolute
/// \param b_relNextReference If true, implies that 'j' is a reference relative to vno, otherwise j is an absolute reference
/// \return    (float) cost of the transition
float   costFunc_defaultDetermine(
  s_env&  st_env,
  s_iterInfo* pst_iterInfo,
  int   vno_c,
  int   j,
  bool  b_relNextReference = true
);

/// \fn float costFunc_unityReturn(s_env& st_env, s_iterInfo* pst_iterInfo, int vno_c, int j, bool b_relNextReference = true);
/// \brief Always returns a '1' as the transition cost between a node and its neighbour.
/// \param  st_env The environment structure
/// \param pst_iterInfo A structure that records vertex transition information (cost, distance, curvature, etc.)
/// \param  vno The current vertex index
/// \param  j The neighbour vertex index - either relative or absolute
/// \param b_relNextReference If true, implies that 'j' is a reference relative to vno, otherwise j is an absolute reference
/// \return    (float) cost of the transition - always 1
float   costFunc_unityReturn(
  s_env&  st_env,
  s_iterInfo* pst_iterInfo,
  int   vno_c,
  int   j,
  bool  b_relNextReference = true
);

/// \fn float costFunc_EuclideanReturn(s_env& st_env, s_iterInfo* pst_iterInfo, int vno_c, int j, bool b_relNextReference = true);
/// \brief Returns the Euclidean distance between a vertex and its neighbour.
/// \param  st_env The environment structure
/// \param pst_iterInfo A structure that records vertex transition information (cost, distance, curvature, etc.)
/// \param  vno The current vertex index
/// \param  j The neighbour vertex index - either relative or absolute
/// \param b_relNextReference If true, implies that 'j' is a reference relative to vno, otherwise j is an absolute reference
/// \return    (float) cost of the transition - always 1
float   costFunc_EuclideanReturn(
  s_env&  st_env,
  s_iterInfo* pst_iterInfo,
  int   vno_c,
  int   j,
  bool  b_relNextReference = true
);

/// \fn float costFunc_distanceReturn(s_env& st_env, s_iterInfo* pst_iterInfo, int vno_c, int j, bool b_relNextReference = true);
/// \brief Returns the distance (as encoded in the MRIS itself) between a vertex and its neighbour. For a given vertex, the
///  distance between itself and neighbouring vertices is recorded in the MRIS data structure itself.
/// \param  st_env The environment structure
/// \param pst_iterInfo A structure that records vertex transition information (cost, distance, curvature, etc.)
/// \param  vno The current vertex index
/// \param  j The neighbour vertex index - either relative or absolute
/// \param b_relNextReference If true, implies that 'j' is a reference relative to vno, otherwise j is an absolute reference
/// \return    (float) cost of the transition - always 1
float   costFunc_distanceReturn(
  s_env&  st_env,
  s_iterInfo* pst_iterInfo,
  int   vno_c,
  int   j,
  bool  b_relNextReference = true
);


#endif //__ENV_H__


