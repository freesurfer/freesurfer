#pragma once
/**
 * @file  mrisurf_aaa.h
 * @brief MRI_SURFACE underlying type declarations and constants.
 *        Separated from mrisurf.h to break some circular includes.
 *
 * constants and structure declarations and some primitive types 
 * for manipulation and i/o of surfaces derived from MRI volumes.
 */
/*
 * Original Author: Bruce Fischl    extracted by Bevin Brett
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2017/02/16 19:42:54 $
 *    $Revision: 1.391 $
 *
 * Copyright © 2011,2019 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "const.h"
#include "matrix.h"
#include "dmatrix.h"


#define MAX_SURFACES 20
#define TALAIRACH_COORDS     0
#define SPHERICAL_COORDS     1
#define ELLIPSOID_COORDS     2

#define VERTICES_PER_FACE    3
#define ANGLES_PER_TRIANGLE  3

#define INFLATED_NAME        "inflated"
#define SMOOTH_NAME          "smoothwm"
#define SPHERE_NAME          "sphere"
#define ORIG_NAME            "orig"
#define WHITE_MATTER_NAME    "white"
#define GRAY_MATTER_NAME     "gray"
#define LAYERIV_NAME         "graymid"
#define GRAYMID_NAME         LAYERIV_NAME
#define MAX_CMDS 1000

#define NEW_VERSION_MAGIC_NUMBER  16777215 // was in mrisurf.c

#define WHICH_FACE_SPLIT(vno0, vno1) (1*nint(sqrt(1.9*vno0) + sqrt(3.5*vno1)))
    //
    // This is used in a repeatable arbitrary true false selector based on the resulting int being EVEN or ODD


typedef struct _mht                 MRIS_HASH_TABLE, MHT ;
typedef struct LABEL_VERTEX         LABEL_VERTEX,    LV  ;
typedef struct LABEL                LABEL;
typedef struct face_type_           face_type, FACE;
typedef struct face_topology_type_  FACE_TOPOLOGY;
typedef struct vertex_type_         vertex_type, VERTEX;
typedef struct VERTEX_TOPOLOGY      VERTEX_TOPOLOGY;



/*
  the vertices in the face structure are arranged in
  counter-clockwise fashion when viewed from the outside.
*/
typedef FixedSizeArray<int,   VERTICES_PER_FACE>   vertices_per_face_t;
typedef FixedSizeArray<float, ANGLES_PER_TRIANGLE> angles_per_triangle_t;

static void copyAnglesPerTriangle(angles_per_triangle_t& dst, angles_per_triangle_t const & src) {
  dst = src;
}
static int cmpAnglesPerTriangle(angles_per_triangle_t const & dst, angles_per_triangle_t const & src) {
  return memcmp(&dst,&src,sizeof(src));
}

typedef struct MRIS_XYZ {
  float x,y,z;
} MRIS_XYZ;


typedef struct {
    unsigned long hash;
} MRIS_HASH;


// MRIS supplies a rich world, but in a format that causes lots of memory traffic
//
typedef struct MRIS MRIS,MRI_SURFACE;       // Prefer using MRIS


// MRIS_MP is a much more efficient supplier of MetricProperties data than MRIS.
// It is implemented in mrisurf_mp.h
//
// It is optimized to cope with the XYZ changing as the shape is mutated to optimize some SSE.  
// Its representation keeps the data in a format that fills cache lines with immediately
// needed information.
//
typedef struct MRIS_MP MRIS_MP;


// The SSE calculation uses some large subsystems, such as MHT, that are coded
// using the MRIS.  Ideally we would use C++, a class derivation hierachy, and 
// virtual functions or C++ templates to implement these functions on top of both 
// MRIS and MRIS_MP
//
// The following is basically a base class with virtual functions.
// It is implemented in mrisurf_MRISBase.h
//
typedef struct MRISBase {
    MRIS_MP*    mris_mp;            // takes precidence over mris
    MRIS*       mris;
} MRISBase;

typedef struct MRISBaseConst {
    MRIS_MP const*    mris_mp;      // takes precidence over mris
    MRIS const*       mris;
} MRISBaseConst;


static MRISBase      MRISBaseCtr     (      MRIS_MP* mris_mp, MRIS       * mris) { MRISBase      base; base.mris_mp = mris_mp; base.mris = mris; return base; }
static MRISBaseConst MRISBaseConstCtr(const MRIS_MP* mris_mp, MRIS const * mris) { MRISBaseConst base; base.mris_mp = mris_mp; base.mris = mris; return base; }
static MRISBaseConst MRISBaseToConst (      MRISBase src                       ) { return MRISBaseConstCtr(src.mris_mp, src.mris); }



//  UnitizeNormalFace is a global variable used in mrisNormalFace() to allow the
//  output norm to be unitized or not. That function computed the norm
//  using a cross product but then did not normalize the result (cross
//  product is not unit length even if inputs are unit). UnitizeNormalFace
//  allows unitization to be turned on and off for testing. Note: skull
//  stripping uses this function, so may want to UnitizeNormalFace=0 when
//  testing effects on surface placement so that the stream is the same up
//  until surface placement.
extern int UnitizeNormalFace;

//  This variable can be used to turn on the hires options
//  in MRIScomputeBorderValues_new()
extern int BorderValsHiRes;

//  This is used to record the actual value and difference
//  into v->valbak and v->val2bak when running mrisRmsValError()
//  for debugging or evaluation purposes.
extern int RmsValErrorRecord;


