#pragma once
/**
 * @brief MRI_SURFACE underlying type declarations and constants.
 *        Separated from mrisurf.h to break some circular includes.
 *
 * constants and structure declarations and some primitive types 
 * for manipulation and i/o of surfaces derived from MRI volumes.
 */
/*
 * Original Author: Bruce Fischl    extracted by Bevin Brett
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

#include "const.h"
#include "matrix.h"
#include "dmatrix.h"

#define MAX_NEIGHBORS (1024)
#define MRIS_MAX_NEIGHBORHOOD_LINKS 50  // bound on nlinks

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

typedef FixedSizeArray<char   *,MAX_CMDS> MRIS_cmdlines_t;
typedef FixedSizeArray<char    ,STRLEN>   MRIS_subject_name_t;
typedef FixedSizeArray<char    ,STRLEN>   MRIS_fname_t;
    //
    // C arrays can not be returned as function results, but there can


enum MRIS_Status_DistanceFormula {
  MRIS_Status_DistanceFormula_0,    // see utils/mrisComputeVertexDistancesWkr_extracted.h
  MRIS_Status_DistanceFormula_1
};

enum MRIS_Status {
#define MRIS_Status_ELTS \
  ELT(MRIS_SURFACE              ,0) SEP \
  ELT(MRIS_PATCH                ,0) SEP \
  ELT(MRIS_PLANE                ,0) SEP \
  ELT(MRIS_ELLIPSOID            ,0) SEP \
  ELT(MRIS_SPHERE               ,1) SEP \
  ELT(MRIS_PARAMETERIZED_SPHERE ,1) SEP \
  ELT(MRIS_RIGID_BODY           ,0) SEP \
  ELT(MRIS_SPHERICAL_PATCH      ,0) SEP \
  ELT(MRIS_UNORIENTED_SPHERE    ,0) SEP \
  ELT(MRIS_PIAL_SURFACE         ,0)     \
  // end of macro
#define SEP ,
#define ELT(E,D) E
  MRIS_Status_ELTS,
  MRIS_Status__end, 
  MRIS_CUT = MRIS_PATCH
#undef ELT
#undef SEP
};

const char* MRIS_Status_text(MRIS_Status s1);
MRIS_Status_DistanceFormula MRIS_Status_distanceFormula(MRIS_Status s1);
bool areCompatible(MRIS_Status s1, MRIS_Status s2);
void checkOrigXYZCompatibleWkr(MRIS_Status s1, MRIS_Status s2, const char* file, int line);
#define checkOrigXYZCompatible(S1,S2) checkOrigXYZCompatibleWkr((S1),(S2),__FILE__,__LINE__);


struct _mht;
typedef struct MRIS_HASH_TABLE MHT;

typedef struct LABEL_VERTEX         LABEL_VERTEX,    LV  ;
typedef struct LABEL                LABEL;
// See mrisurf_FACE_VERTEX_MRIS_generated.h for FACE, VERTEX, and MRI_SURFACE/MRIS
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


typedef struct _area_label
{
  char     name[STRLEN] ;     /* name of region */
  float    cx ;               /* centroid x */
  float    cy ;               /* centroid y */
  float    cz ;               /* centroid z */
  int      label ;            /* an identifier (used as an index) */
}
MRIS_AREA_LABEL ;

typedef struct FaceNormCacheEntry {
    // inputs
        // may have to capture them if the inputs change
    // flag saying the calculation has been deferred
    // may be better to store these separately...
    // value
        float nx,ny,nz,orig_area;
} FaceNormCacheEntry;

typedef struct FaceNormDeferredEntry {
    char deferred;
} FaceNormDeferredEntry;

// the face norm elements have moved into the FaceNormalCacheEntry
/*  ELTT(float,nx) SEP    \
    ELTT(float,ny) SEP    \
    ELTT(float,nz) SEP    \
    ELTT(float,orig_area) SEP    \
*/

typedef struct edge_type_
{
  // topology
  int edgeno; // this edge no
  int vtxno[4]; // vertex numbers of 2 ends + 2 opposites
  int faceno[2]; // two adjacent faces
  unsigned char corner[4][2]; // corner[nthvtx][faceno]
  // metrics
  double len; // length of the edge
  double dot; // dot product of the adjacent face normals
  double angle; // angle (deg) of the adjacent face normals
  double J; // Angle Cost of this edge
  double u[3]; // unit vector pointing from v0 to v1
  double area[2], maxarea; // area of the two faces and max of the two
  DMATRIX *gradU; // 1x3 grad of unit verctor wrt vertex 0
  DMATRIX *gradDot[4]; // 3x3 grad of dot product wrt 4 vertices
} MRI_EDGE;

// Corner of a triangluar face 
typedef struct corner_type_
{
  // topology
  int cornerno; // this corner number
  int vtxno[3]; // vertex 0 is the source
  int faceno; // face this corner belongs to
  int edgeno[2]; // edge numbers
  int edgedir[2]; // direction of the edge relative to the corner
  // metrics
  double dot; // dot product of the adjacent segments
  double angle; // angle (deg) of the corner
  double J; // cost, eg, (dot-0.5)^2 where 0.5 = cos(60) = equilateral tri
  DMATRIX *gradDot[3]; // 3 1x3 grad of dot wrt each vertex
} MRI_CORNER;


#include "colortab.h" // 'COLOR_TABLE'

#define VERTEX_SULCAL  0x00000001L

typedef struct
{
  int nvertices;
  unsigned int *vertex_indices;
}
STRIP;

#include "transform.h" // TRANSFORM, LTA


typedef int*                    pSeveralInt;
typedef uchar*                  pSeveralUchar;
typedef float*                  pSeveralFloat;
typedef float const*            pSeveralConstFloat;
typedef void*                   p_void;
typedef void**                  p_p_void;

typedef DMATRIX*                PDMATRIX;
typedef MATRIX*                 PMATRIX;
typedef LTA*                    PLTA;
typedef COLOR_TABLE*            PCOLOR_TABLE;

typedef MRI*                    PMRI;
typedef MRI_EDGE*               pSeveralMRI_EDGE;
typedef MRI_CORNER*             pSeveralMRI_CORNER;
typedef MRIS_AREA_LABEL*        PMRIS_AREA_LABEL;
typedef VERTEX *                PVERTEX;
typedef FixedSizeArray<PDMATRIX,3> A3PDMATRIX;

typedef STRIP*                  pSeveralSTRIP;
typedef VERTEX *                pSeveralVERTEX;
typedef VERTEX_TOPOLOGY *       pSeveralVERTEX_TOPOLOGY;
typedef FACE *                  pSeveralFACE;
typedef FaceNormCacheEntry*     pSeveralFaceNormCacheEntry;
typedef FaceNormDeferredEntry*  pSeveralFaceNormDeferredEntry;

 

// MRIS supplies a rich world, but in a format that causes lots of memory traffic
//
// It is defined in mrisurf_FACE_VERTEX_MRIS_generated.h
//
typedef struct MRIS MRIS,MRI_SURFACE;       // Prefer using MRIS


// MRIS_MP is a much more efficient supplier of MetricProperties data than MRIS,
// but can get some of its properties from an underlying Surface so it doesn't have to 
// implement or copy ones that are rarely used.
//
// It is defined in mrisurf_MRIS_MPPropertiesInVectors.h with some additional macros in mrisurf_mp.h
//
// It is optimized to cope with the XYZ changing as the shape is mutated to optimize some SSE.  
// Its representation keeps the data in a format that fills cache lines with immediately
// needed information.
//
typedef struct MRIS_MP MRIS_MP;


// MRISPV is a much more efficient supplier of MetricProperties data than MRIS
// and hold all the data.  It is not yet fully implemented.
//
// It is defined in mrisurf_MRIS_PropertiesInVectors.h
//
// It is optimized to cope with the XYZ changing as the shape is mutated to optimize some SSE.  
// Its representation keeps the data in a format that fills cache lines with immediately
// needed information.
//
typedef struct MRISPV MRISPV;


// The SSE calculation uses some large subsystems, such as MHT, that are coded
// using the MRIS.  Ideally we would use C++, a class derivation hierachy, and 
// virtual functions or C++ templates to implement these functions on top of both 
// MRIS and MRIS_MP
//
// The following precedes the change to C++,
// should be elikminated asap.
// It is a crude implementation of a base class with virtual functions.
//
// It is further defined in mrisurf_MRISBase.h
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


// Class used to control ComputeBorderValues()
class  CBV_OPTIONS {
 public:
  MRIS *cbvsurf;
  MRI *LocalMaxFound; // Keep track of which vertices had a local max found
  // Alternative border low threshold - allows specification of a different
  // border_low threshold in selected regions (eg, for high myelin)
  double AltBorderLowFactor=0;
  double AltBorderLow=0;
  char *AltBorderLowLabelFile=NULL;
  LABEL *AltBorderLowLabel=NULL;
  MRI *AltBorderLowMask=NULL;
  // functions
  int Alloc(void);
  int ReadAltBorderLowLabel(void);
};
extern CBV_OPTIONS CBVO;

//  These variables can be used to turn on peak-finding options
//  in MRIScomputeBorderValues_new()
extern int CBVfindFirstPeakD1;
extern int CBVfindFirstPeakD2;

//  This is used to record the actual value and difference
//  into v->valbak and v->val2bak when running mrisRmsValError()
//  for debugging or evaluation purposes.
extern int RmsValErrorRecord;


