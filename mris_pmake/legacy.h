/**
 * @brief The environment object API.
 *
 * This file contains old legacy code, consolidated here during testing
 * to facilitate its eventual removal from the mainline code.
 *
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

#ifndef __LEGACY_H__
#define __LEGACY_H__


#include "env.h"
#include "scanopt.h"

// Forward declaration of environment structure
typedef struct _env s_env;

/// Weights and values for the main cost function polynomial.
typedef struct _weights {
  float wd;         // distance
  float wc;         // curve
  float wh;         // sulcal height
  float wdc;        // distance X curve
  float wdh;        // distance X height
  float wch;        // curve X height
  float wdch;       // distance X curve X height
  float wdir;       // direction vector
}
s_weights;

void
s_weights_scan(
  s_weights&  st_costWeight,
  C_scanopt&  cso_options
);

void s_weights_copy(
    s_weights&		sw_target,
    s_weights&		sw_source
);

void
s_weights_setAll(
    s_weights&  asw,
    float       af
);

void
s_Dweights_setAll(
    s_Dweights&     asw,
    float           af
);

void s_Dweights_copy(
    s_Dweights&		sw_target,
    s_Dweights&		sw_source
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
  float Dwd;            // distance
  float Dwc;            // curve
  float Dwh;            // sulcal height
  float Dwdc;           // distance X curve
  float Dwdh;           // distance X height
  float Dwch;           // curve X height
  float Dwdch;          // distance X curve X height
  float Dwdir;          // direction vector
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
//typedef struct _iterInfo {
//  int       iter;
//  float     f_distance;
//  float     f_curvature;
//  float     f_sulcalHeight;
//  float     f_dir;
//}
//s_iterInfo;

//typedef enum {
//  e_default, e_unity, e_euclid, e_distance
//} e_COSTFUNCTION;

typedef struct _iterInfo	s_iterInfo;

void s_env_costFctList(
    s_env&      ast_env
);

int s_env_costFctSetIndex(
    s_env*      apst_env,
    int         aindex
);

void  s_env_costFctSet(
    s_env*      pst_env,
    float       (*cost_fct) (
        s_env&          st_env,
        s_iterInfo*     pst_iterInfo,
        int             vno_c,
        int             j,
        bool            b_relNextReference
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


/// \fn void asynchEvent_processDWGHT(s_env& st_env, string str_comms)
/// \brief Process socket-based access to the Dweight structure
/// \param ast_env   Core simulation environment
/// \param astr_comms  Control string parsed from socket
/// \return  true|false  Status: success or fail.
bool
asynchEvent_processDWGHT(
  s_env&    ast_env,
  string    astr_comms
);

/// \fn void asynchEvent_processWGHT(s_env& st_env, string str_comms)
/// \brief Process socket-based access to the core weight structure
/// \param ast_env   Core simulation environment
/// \param astr_comms  Control string parsed from socket
/// \return  true|false  Status: success or fail.
bool
asynchEvent_processWGHT(
  s_env&    ast_env,
  string    astr_comms
);

#endif //__LEGACY_H__


