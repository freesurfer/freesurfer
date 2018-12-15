#pragma once
/**
 * @file  mrisurf.h
 * @brief MRI_SURFACE utilities.
 *
 * Utilities, constants and structure definitions for manipulation
 * and i/o of surfaces derived from MRI volumes.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2017/02/16 19:42:54 $
 *    $Revision: 1.391 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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



#include "minc_volume_io.h"
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


typedef struct _area_label
{
  char     name[STRLEN] ;     /* name of region */
  float    cx ;               /* centroid x */
  float    cy ;               /* centroid y */
  float    cz ;               /* centroid z */
  int      label ;            /* an identifier (used as an index) */
}
MRIS_AREA_LABEL ;

struct MRIS;

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

/*
  the vertices in the face structure are arranged in
  counter-clockwise fashion when viewed from the outside.
*/
typedef int   vertices_per_face_t[VERTICES_PER_FACE];
typedef float angles_per_triangle_t[ANGLES_PER_TRIANGLE];

// the face norm elements have moved into the FaceNormalCacheEntry
/*  ELTT(float,nx) SEP    \
    ELTT(float,ny) SEP    \
    ELTT(float,nz) SEP    \
    ELTT(float,orig_area) SEP    \
*/

typedef struct edge_type_
{
  int edgeno; // this face no
  int vtxno[4]; // vertex numbers of 2 ends + 2 opposites
  int faceno[2]; // two adjacent faces
  unsigned char corner[4][2]; // corner[nthvtx][faceno]
  double len; // length of the edge
  double dot; // dot product of the adjacent face normals
  double angle; // angle (deg) of the adjacent face normals
  double J; // Angle Cost of this edge
  DMATRIX *gradDot[4]; // 3x3 grad of dot product wrt 4 vertices
} MRI_EDGE;


#if defined(COMPILING_MRISURF_TOPOLOGY) || defined(COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED)
#define CONST_EXCEPT_MRISURF_TOPOLOGY 
#else
#define CONST_EXCEPT_MRISURF_TOPOLOGY const
#endif

#if defined(COMPILING_MRISURF_METRIC_PROPERTIES) || defined(COMPILING_MRISURF_METRIC_PROPERTIES_FRIEND)
#define CONST_EXCEPT_MRISURF_METRIC_PROPERTIES 
#else
#define CONST_EXCEPT_MRISURF_METRIC_PROPERTIES const
#endif
    //
    // Used to find and control where various fields are written
    
typedef struct face_type_
{
#define LIST_OF_FACE_ELTS_1    \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY vertices_per_face_t,v) SEP               /* vertex numbers of this face */    \
  ELTT(float,area) SEP    \
  ELTT(angles_per_triangle_t,angle) SEP    \
  ELTT(angles_per_triangle_t,orig_angle) SEP    \
  ELTT(char,ripflag) SEP                        /* ripped face */    \
  ELTT(char,oripflag) SEP                       /* stored version */    \
  ELTT(int,marked) SEP                          /* marked face */    \
  ELTP(DMATRIX,norm) SEP  /* 3x1 normal vector */ \
  ELTP(DMATRIX,gradNorm[3]) SEP  /* 3x3 Gradient of the normal wrt each of the 3 vertices*/ 
    // end of macro

#if 0
  float logshear,shearx,sheary;  /* compute_shear */
#endif

// Why does mrishash need these?  Where else are they used?
#if 0
#define LIST_OF_FACE_ELTS_2    \
  ELTT(float,cx) SEP    \
  ELTT(float,cy) SEP    \
  ELTT(float,cz) SEP         /* coordinates of centroid */   \
    // end of macro
#define LIST_OF_FACE_ELTS \
    LIST_OF_FACE_ELTS_1 SEP \
    LIST_OF_FACE_ELTS_2 \
    // end of macro
#else
#define LIST_OF_FACE_ELTS \
    LIST_OF_FACE_ELTS_1
#endif

#define ELTT(T,N) T N;
#define ELTP(TARGET,NAME) TARGET *NAME ;
#define SEP
LIST_OF_FACE_ELTS
#undef SEP
#undef ELTT
#undef ELTP

}
face_type, FACE ;

#ifndef uchar
#define uchar  unsigned char
#endif

#include "colortab.h" // 'COLOR_TABLE'


#define LIST_OF_VERTEX_TOPOLOGY_ELTS \
  /* put the pointers before the ints, before the shorts, before uchars, to reduce size  */ \
  /* the whole fits in much less than one cache line, so further ordering is no use      */ \
  ELTP(int,f) SEP                                               /* array[v->num] the fno's of the neighboring faces         */ \
  ELTP(uchar,n) SEP           	                                /* array[v->num] the face.v[*] index for this vertex        */ \
  ELTP(int,e) SEP                                               /* edge state for neighboring vertices                      */ \
  ELTP(CONST_EXCEPT_MRISURF_TOPOLOGY int,v) SEP                 /* array[v->vtotal or more] of vno, head sorted by hops     */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY short,vnum)                /* number of 1-hop neighbots    should use [p]VERTEXvnum(i) */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY short,v2num) SEP           /* number of 1, or 2-hop neighbors                          */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY short,v3num) SEP           /* number of 1,2,or 3-hop neighbors                         */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY short,vtotal) SEP          /* total # of neighbors. copy of vnum.nsizeCur              */ \
  ELTX(CONST_EXCEPT_MRISURF_TOPOLOGY short,nsizeMaxClock) SEP   /* copy of mris->nsizeMaxClock when v#num                   */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY uchar,nsizeMax) SEP        /* the max nsize that was used to fill in vnum etc          */ \
  ELTT(CONST_EXCEPT_MRISURF_TOPOLOGY uchar,nsizeCur) SEP        /* index of the current v#num in vtotal                     */ \
  ELTT(uchar,num) SEP                                           /* number of neighboring faces                              */ \
  // end of macro

//static short* pVERTEXvnum(VERTEX_TOPOLOGY * v, int i);
//static short VERTEXvnum(VERTEX_TOPOLOGY const * v, int i);

// The above elements historically were in the VERTEX
// and can still be there by
//  having VERTEX_TOPOLOGY be a typedef of VERTEX
//  having the mris->vertices and the mris->vertices_topology be the same pointer
// and this is what the code was doing until the separation was completed.
//
#define SEPARATE_VERTEX_TOPOLOGY
#ifndef SEPARATE_VERTEX_TOPOLOGY

#define LIST_OF_VERTEX_TOPOLOGY_ELTS_IN_VERTEX LIST_OF_VERTEX_TOPOLOGY_ELTS SEP

#else

typedef struct VERTEX_TOPOLOGY {
    // The topology of the vertex describes its neighbors
    // but not its position nor any properties derived from its position.
    //
    // This data is not changed as the vertices are moved during distortion of the polyhedra.
    //

#define SEP
#define ELTX(TYPE,NAME) TYPE NAME ;
#define ELTT(TYPE,NAME) TYPE NAME ;
#define ELTP(TARGET,NAME) TARGET *NAME ;
  LIST_OF_VERTEX_TOPOLOGY_ELTS
#undef ELTP
#undef ELTT
#undef ELTX
#undef SEP
} VERTEX_TOPOLOGY;

#define LIST_OF_VERTEX_TOPOLOGY_ELTS_IN_VERTEX

#endif


typedef struct MRIS_XYZ {
  float x,y,z;
} MRIS_XYZ;


typedef struct vertex_type_
{
// The LIST_OF_VERTEX_ELTS macro used here enables the the mris_hash
// and other algorithms to process all the elements without having to explicitly name them there and here
//
// The order is important because
//      . it affects the hash (until a better hash algorithm is implemented)
//      . it affects the size (each item must be aligned on its appropriate boundary for its size)
//      . it affects the number of cache lines that must be read and written
//
// By putting the ripflag at the end, reading it will cause its whole cache line to be read, perhaps the 
// first few elements of the next vertex earlier, and the ripflag test will probably be correctly predicted
// and so the cpu won't wait.
//
#define LIST_OF_VERTEX_ELTS_1    \
  LIST_OF_VERTEX_TOPOLOGY_ELTS_IN_VERTEX \
  /* managed by MRISfreeDists[_orig] and MRISmakeDists[_orig] */ \
  ELTX(float* const,dist)      SEP                                              /* distance to neighboring vertices */          \
  ELTX(float* const,dist_orig) SEP                                              /* original distance to neighboring vertices */ \
  ELTX(int,dist_orig_capacity) SEP \
  \
  ELTT(const float,origx)                                       SEP             /* original coordinates */                      \
  ELTT(const float,origy)                                       SEP             /* use MRISsetOriginalXYZ() to set */           \
  ELTT(const float,origz)                                       SEP                                                             \
  \
  ELTT(/*CONST_EXCEPT_MRISURF_METRIC_PROPERTIES*/ float,x)          SEP             /* current coordinates */                       \
  ELTT(/*CONST_EXCEPT_MRISURF_METRIC_PROPERTIES*/ float,y)          SEP             /* use MRISsetXYZ() to set */                   \
  ELTT(/*CONST_EXCEPT_MRISURF_METRIC_PROPERTIES*/ float,z)          SEP                                                             \
  \
  ELTT(float,nx) SEP    \
  ELTT(float,ny) SEP    \
  ELTT(float,nz) SEP        /* curr normal */    \
  ELTT(float,pnx) SEP    \
  ELTT(float,pny) SEP    \
  ELTT(float,pnz) SEP     /* pial normal */    \
  /* the above is the first cache line */ \
  ELTT(float,wnx) SEP    \
  ELTT(float,wny) SEP    \
  ELTT(float,wnz) SEP     /* white normal */    \
  ELTT(float,onx) SEP    \
  ELTT(float,ony) SEP    \
  ELTT(float,onz) SEP     /* original normal */    \
  ELTT(float,dx) SEP    \
  ELTT(float,dy) SEP    \
  ELTT(float,dz) SEP     /* current change in position */    \
  ELTT(float,odx) SEP    \
  ELTT(float,ody) SEP    \
  ELTT(float,odz) SEP  /* last change of position (for momentum) */    \
  ELTT(float,tdx) SEP    \
  ELTT(float,tdy) SEP    \
  ELTT(float,tdz) SEP  /* temporary storage for averaging gradient */    \
  ELTT(float,curv) SEP            /* curr curvature */    \
  ELTT(float,curvbak) SEP    \
  ELTT(float,val) SEP             /* scalar data value (file: rh.val, sig2-rh.w) */    \
  ELTT(float,imag_val) SEP       /* imaginary part of complex data value */    \
  ELTT(float,cx) SEP    \
  ELTT(float,cy) SEP    \
  ELTT(float,cz) SEP     /* coordinates in canonical coordinate system */    \
  ELTT(float,tx) SEP    \
  ELTT(float,ty) SEP    \
  ELTT(float,tz) SEP     /* tmp coordinate storage */    \
  ELTT(float,tx2) SEP    \
  ELTT(float,ty2) SEP    \
  ELTT(float,tz2) SEP  /* tmp coordinate storage */    \
  ELTT(float,targx) SEP    \
  ELTT(float,targy) SEP    \
  ELTT(float,targz) SEP   /* target coordinates */    \
  ELTT(float,pialx) SEP    \
  ELTT(float,pialy) SEP    \
  ELTT(float,pialz) SEP   /* pial surface coordinates */    \
  ELTT(float,whitex) SEP    \
  ELTT(float,whitey) SEP    \
  ELTT(float,whitez) SEP/* white surface coordinates */    \
  ELTT(float,l4x) SEP    \
  ELTT(float,l4y) SEP    \
  ELTT(float,l4z) SEP   /* layerIV surface coordinates */    \
  ELTT(float,infx) SEP    \
  ELTT(float,infy) SEP    \
  ELTT(float,infz) SEP /* inflated coordinates */    \
  ELTT(float,fx) SEP    \
  ELTT(float,fy) SEP    \
  ELTT(float,fz) SEP      /* flattened coordinates */    \
  ELTT(int,px) SEP    \
  ELTT(int,qx) SEP    \
  ELTT(int,py) SEP    \
  ELTT(int,qy) SEP    \
  ELTT(int,pz) SEP    \
  ELTT(int,qz) SEP /* rational coordinates for exact calculations */    \
  ELTT(float,e1x) SEP    \
  ELTT(float,e1y) SEP    \
  ELTT(float,e1z) SEP  /* 1st basis vector for the local tangent plane */    \
  ELTT(float,e2x) SEP    \
  ELTT(float,e2y) SEP    \
  ELTT(float,e2z) SEP  /* 2nd basis vector for the local tangent plane */    \
  ELTT(float,pe1x) SEP    \
  ELTT(float,pe1y) SEP    \
  ELTT(float,pe1z) SEP  /* 1st basis vector for the local tangent plane */    \
  ELTT(float,pe2x) SEP    \
  ELTT(float,pe2y) SEP    \
  ELTT(float,pe2z) SEP  /* 2nd basis vector for the local tangent plane */    \
  // end of macro

#if 0
#define LIST_OF_VERTEX_ELTS_2 \
  float dipx;    \
  float ipy;    \
  float ipz;  /* dipole position */    \
  float dipnx;    \
  float ipny;    \
  float ipnz; /* dipole orientation */    \
  // end of macro
#else
#define LIST_OF_VERTEX_ELTS_2
#endif

#define LIST_OF_VERTEX_ELTS_3   \
  ELTT(float,nc) SEP              /* curr length normal comp */    \
  ELTT(float,val2) SEP            /* complex comp data value (file: sig3-rh.w) */    \
  ELTT(float,valbak) SEP          /* scalar data stack */    \
  ELTT(float,val2bak) SEP         /* complex comp data stack */    \
  ELTT(float,stat) SEP            /* statistic */    \
    \
  ELTT(int,undefval) SEP          /* [previously dist=0] */    \
  ELTT(int,old_undefval) SEP      /* for smooth_val_sparse */    \
  ELTT(int,fixedval) SEP          /* [previously val=0] */    \
    \
  ELTT(float,fieldsign) SEP       /* fieldsign--final: -1,0,1 (file: rh.fs) */    \
  ELTT(float,fsmask) SEP          /* significance mask (file: rh.fm) */    \
  ELTT(float,d) SEP              /* for distance calculations */    \
  // end of macro
  
#if 0
#define LIST_OF_VERTEX_ELTS_4 \
  ELTP(float,tri_area) SEP      /* array of triangle areas - num long */    \
  ELTP(float,orig_tri_area) SEP /* array of original triangle areas - num long */    \
  ELTT(float,dist) SEP            /* dist from sampled point [defunct: or 1-cos(a)] */    \
  ELTT(float,ox) SEP    \
  ELTT(float,y) SEP    \
  ELTT(float,z) SEP        /* last position (for undoing time steps) */    \
  ELTT(float,mx) SEP    \
  ELTT(float,y) SEP    \
  ELTT(float,z) SEP        /* last movement */    \
  ELTT(float,onc) SEP             /* last length normal comp */    \
  ELTT(float,oval) SEP            /* last scalar data (for smooth_val) */    \
  ELTP(float,fnx) SEP           /* face normal - x component */    \
  ELTP(float,fny) SEP           /* face normal - y component */    \
  ELTP(float,fnz) SEP           /* face normal - z component */    \
  ELTT(float,bnx) SEP    \
  ELTT(float,ny) SEP    \
  ELTT(float,bnx) SEP    \
  ELTT(float,bny) SEP /* boundary normal */    \
  ELTP(float,tri_angle) SEP     /* angles of each triangle this vertex belongs to */    \
  ELTP(float,orig_tri_angle) SEP/* original values of above */    \
  ELTT(float,stress) SEP          /* explosion */    \
  ELTT(float,logshear) SEP    \
  ELTT(float,hearx) SEP    \
  ELTT(float,heary) SEP    \
  ELTT(float,shearx) SEP    \
  ELTT(float,sheary) SEP  /* for shear term */    \
  ELTT(float,ftmp) SEP           /* temporary floating pt. storage */    \
  ELTT(float,logarat) SEP    \
  ELTT(float,logarat) SEP    \
  ELTT(float,qrtarat) SEP /* for area term */    \
  ELTT(float,smx) SEP    \
  ELTT(float,my) SEP    \
  ELTT(float,mz) SEP    \
  ELTT(float,smx) SEP    \
  ELTT(float,smy) SEP    \
  ELTT(float,smz) SEP/* smoothed curr,last move */    \
  // end of macro
#else
#define LIST_OF_VERTEX_ELTS_4
#endif

#define LIST_OF_VERTEX_ELTS_5   \
  ELTT(int,annotation) SEP      /* area label (defunct--now from label file name!) */    \
  ELTT(char,oripflag) SEP    \
  ELTT(char,origripflag) SEP     /* cuts flags */    \
  // end of macro

#if 0
#define LIST_OF_VERTEX_ELTS_6
  float coords[3];
#else
#define LIST_OF_VERTEX_ELTS_6

#endif

#define LIST_OF_VERTEX_ELTS_7    \
  ELTP(void,vp) SEP                     /* to store user's information */    \
  ELTT(float,theta) SEP    \
  ELTT(float,phi) SEP               /* parameterization */    \
  ELTT(float,area) SEP    \
  ELTT(float,origarea) SEP    \
  ELTT(float,group_avg_area) SEP    \
  ELTT(float,K) SEP                 /* Gaussian curvature */    \
  ELTT(float,H) SEP                 /* mean curvature */    \
  ELTT(float,k1) SEP    \
  ELTT(float,k2) SEP                    /* the principal curvatures */    \
  ELTT(float,mean) SEP    \
  ELTT(float,mean_imag) SEP         /* imaginary part of complex statistic */    \
  ELTT(float,std_error) SEP    \
  ELTT(unsigned int,flags) SEP    \
  ELTT(int,fno) SEP                 /* face that this vertex is in */    \
  ELTT(int,cropped) SEP \
  ELTT(short,marked) SEP            /* for a variety of uses */    \
  ELTT(short,marked2) SEP    \
  ELTT(short,marked3) SEP    \
  ELTT(char,neg) SEP                /* 1 if the normal vector is inverted */    \
  ELTT(char,border) SEP             /* flag */    \
  ELTT(char,ripflag)                /* vertex no longer exists - placed last to load the next vertex into cache */ \
  // end of macro
  
#define LIST_OF_VERTEX_ELTS \
  LIST_OF_VERTEX_ELTS_1   \
  LIST_OF_VERTEX_ELTS_2   \
  LIST_OF_VERTEX_ELTS_3   \
  LIST_OF_VERTEX_ELTS_4   \
  LIST_OF_VERTEX_ELTS_5   \
  LIST_OF_VERTEX_ELTS_6   \
  LIST_OF_VERTEX_ELTS_7   \
  // end of macro

#define SEP
#define ELTX(TYPE,NAME) TYPE NAME ;
#define ELTT(TYPE,NAME) TYPE NAME ;
#define ELTP(TARGET,NAME) TARGET *NAME ;
  LIST_OF_VERTEX_ELTS
#undef ELTP
#undef ELTT
#undef ELTX
#undef SEP

#if defined(__cplusplus)
    // C++ requires const members be initialized
    vertex_type_() : dist(nullptr), dist_orig(nullptr), origx(0), origy(0), origz(0), x(0), y(0), z(0) {}
#endif

}
vertex_type, VERTEX ;


#ifndef SEPARATE_VERTEX_TOPOLOGY
typedef vertex_type VERTEX_TOPOLOGY; 
#endif

static CONST_EXCEPT_MRISURF_TOPOLOGY short* pVERTEXvnum(VERTEX_TOPOLOGY CONST_EXCEPT_MRISURF_TOPOLOGY * v, int i) {
  switch (i) {
  case 1: return &v->vnum;
  case 2: return &v->v2num;
  case 3: return &v->v3num;
  default: cheapAssert(false); return NULL;
  }    
}
static short VERTEXvnum(VERTEX_TOPOLOGY const * v, int i) {
  switch (i) {
  case 1: return v->vnum;
  case 2: return v->v2num;
  case 3: return v->v3num;
  default: cheapAssert(false); return 0;
  }    
}


#define VERTEX_SULCAL  0x00000001L

typedef struct
{
  int nvertices;
  unsigned int *vertex_indices;
}
STRIP;

#include "transform.h" // TRANSFORM, LTA

typedef char   *MRIS_cmdlines_t[MAX_CMDS] ;
typedef char    MRIS_subject_name_t[STRLEN] ;
typedef char    MRIS_fname_t[STRLEN] ;

typedef struct MRIS
{
// The LIST_OF_MRIS_ELTS macro used here enables the the mris_hash
// and other algorithms to process all the elements without having to explicitly name them there and here
//
#define LIST_OF_MRIS_ELTS_1     \
    \
  ELTT(const int,nverticesFrozen) SEP           /* # of vertices on surface is frozen */                                                    \
  ELTT(const int,nvertices) SEP                 /* # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al */         \
  ELTT(const int,nfaces) SEP                    /* # of faces on surface,    change by calling MRISreallocVerticesAndFaces et al */         \
  ELTT(const bool,faceAttachmentDeferred) SEP   /* defer connecting faces to vertices, for performance reasons                   */         \
  ELTT(int,nedges) SEP                          /* # of edges on surface*/    \
  ELTT(int,nstrips) SEP    \
  ELTP(VERTEX_TOPOLOGY,vertices_topology) SEP    \
  ELTT(const int,tempsAssigned) SEP             /* State of various temp fields that can be borrowed if not already in use   */    \
  ELTP(VERTEX,vertices) SEP    \
  ELTP(FACE,faces) SEP    \
  ELTP(MRI_EDGE,edges) SEP    \
  ELTP(FaceNormCacheEntry,faceNormCacheEntries) SEP \
  ELTP(FaceNormDeferredEntry,faceNormDeferredEntries) SEP \
  ELTP(STRIP,strips) SEP    \
  ELTT(float,xctr) SEP    \
  ELTT(float,yctr) SEP    \
  ELTT(float,zctr) SEP    \
  ELTT(float,xlo) SEP    \
  ELTT(float,ylo) SEP    \
  ELTT(float,zlo) SEP    \
  ELTT(float,xhi) SEP    \
  ELTT(float,yhi) SEP    \
  ELTT(float,zhi) SEP    \
  ELTT(float,x0) SEP             /* center of spherical expansion */    \
  ELTT(float,y0) SEP    \
  ELTT(float,z0) SEP    \
  ELTP(VERTEX,v_temporal_pole) SEP    \
  ELTP(VERTEX,v_frontal_pole) SEP    \
  ELTP(VERTEX,v_occipital_pole) SEP    \
  ELTT(float,max_curv) SEP    \
  ELTT(float,min_curv) SEP    \
  ELTT(float,total_area) SEP    \
  ELTT(double,avg_vertex_area) SEP    \
  ELTT(const double,avg_vertex_dist) SEP  /* set by MRIScomputeAvgInterVertexDist */ \
  ELTT(double,std_vertex_dist) SEP    \
  ELTT(float,orig_area) SEP    \
  ELTT(float,neg_area) SEP    \
  ELTT(float,neg_orig_area) SEP   /* amount of original surface in folds */    \
  ELTT(int,zeros) SEP    \
  ELTT(int,hemisphere) SEP      /* which hemisphere */    \
  ELTT(int,initialized) SEP \
  // end of macro

#if 0

#define LIST_OF_MRIS_ELTS_2 \
  ELTT(General_transform,transform) SEP   /* the next two are from this struct (MNI transform) */    \
  ELTP(Transform,linear_transform) SEP    \
  ELTP(Transform,inverse_linear_transform) SEP \
  // end of macro

#else

#define LIST_OF_MRIS_ELTS_2 \
  // end of macro

#endif

#define LIST_OF_MRIS_ELTS_3     \
  ELTP(LTA,lta) SEP    \
  ELTP(MATRIX,SRASToTalSRAS_) SEP    \
  ELTP(MATRIX,TalSRASToSRAS_) SEP    \
  ELTT(int,free_transform) SEP    \
  ELTT(double,radius) SEP           /* radius (if status==MRIS_SPHERE) */    \
  ELTT(float,a) SEP    \
  ELTT(float,b) SEP    \
  ELTT(float,c) SEP                 /* ellipsoid parameters */    \
  ELTT(MRIS_fname_t,fname) SEP      /* file it was originally loaded from */    \
  ELTT(float,Hmin) SEP              /* min mean curvature */    \
  ELTT(float,Hmax) SEP              /* max mean curvature */    \
  ELTT(float,Kmin) SEP              /* min Gaussian curvature */    \
  ELTT(float,Kmax) SEP              /* max Gaussian curvature */    \
  ELTT(double,Ktotal) SEP           /* total Gaussian curvature */    \
  ELTT(int,status) SEP              /* type of surface (e.g. sphere, plane) */    \
  ELTT(int,patch) SEP               /* if a patch of the surface */    \
  ELTT(int,nlabels) SEP    \
  ELTP(MRIS_AREA_LABEL,labels) SEP  /* nlabels of these (may be null) */    \
  \
  ELTT(char,nsize) SEP              /* size of neighborhoods, or -1 */    \
  ELTT(uchar,vtotalsMightBeTooBig) SEP /* MRISsampleDistances sets this */ \
  ELTX(short,nsizeMaxClock) SEP     /* changed whenever an edge is added or removed, which invalidates the vertex v#num values */ \
  ELTT(char,max_nsize) SEP          /* max the neighborhood size has been set to (typically 3) */    \
  ELTT(char,dist_nsize) SEP         /* max MRISsetNeighborhoodSizeAndDist has computed distances out to */ \
  ELTT(char,dist_alloced_flags) SEP /* two flags, set when any dist(1) or dist_orig(2) allocated */ \
  \
  ELTT(float,avg_nbrs) SEP          /* mean # of vertex neighbors */    \
  ELTP(void,vp) SEP                 /* for misc. use */    \
  ELTT(float,alpha) SEP             /* rotation around z-axis */    \
  ELTT(float,beta) SEP             /* rotation around y-axis */    \
  ELTT(float,gamma) SEP            /* rotation around x-axis */    \
  ELTT(float,da) SEP    \
  ELTT(float,db) SEP    \
  ELTT(float,dg) SEP                /* old deltas */    \
  ELTT(int,type) SEP                /* what type of surface was this initially*/    \
  ELTT(const int,max_vertices) SEP  /* may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces */    \
  ELTT(const int,max_faces) SEP     /* may be bigger than nfaces, set by calling MRISreallocVerticesAndFaces */    \
  ELTT(MRIS_subject_name_t,subject_name) SEP /* name of the subject */    \
  ELTT(float,canon_area) SEP    \
  ELTT(int,noscale) SEP          /* don't scale by surface area if true */    \
  ELTP(float,dx2) SEP             /* an extra set of gradient (not always alloced) */    \
  ELTP(float,dy2) SEP    \
  ELTP(float,dz2) SEP    \
  ELTP(COLOR_TABLE,ct) SEP    \
  ELTT(int,useRealRAS) SEP        /* if 0, vertex position is a conformed volume RAS with c_(r,a,s)=0       */    \
                                  /* if 1, vertex position is a real RAS (volume stored RAS)                */    \
                                  /* The default is 0.                                                      */    \
  ELTT(VOL_GEOM,vg) SEP           /* volume info from which this surface is created. valid iff vg.valid = 1 */    \
  ELTX(MRIS_cmdlines_t, cmdlines) SEP    \
  ELTT(int,ncmds) SEP    \
  ELTT(float,group_avg_surface_area) SEP    /* average of total surface area for group */       \
  ELTT(int,group_avg_vtxarea_loaded) SEP    /* average vertex area for group at each vertex */  \
  ELTT(int,triangle_links_removed) SEP      /* for quad surfaces                         */     \
  ELTP(void,user_parms) SEP                 /* for whatever the user wants to hang here  */     \
  ELTP(MATRIX,m_sras2vox) SEP               /* for converting surface ras to voxel       */     \
  ELTP(MRI,mri_sras2vox) SEP                /* volume that the above matrix is for       */     \
  ELTP(void,mht) SEP \
  ELTP(void,temps)  \
  // end of macro
  
#define LIST_OF_MRIS_ELTS       \
    LIST_OF_MRIS_ELTS_1         \
    LIST_OF_MRIS_ELTS_2         \
    LIST_OF_MRIS_ELTS_3         \
    // end of macro

#define SEP ;
#define ELTP(TARGET, MBR)   TARGET *MBR     // pointers 
#define ELTT(TYPE,   MBR)   TYPE    MBR     // other members that should     be included in the hash
#define ELTX(TYPE,   MBR)   TYPE    MBR     // other members that should NOT be included in the hash
LIST_OF_MRIS_ELTS ;
#undef ELTX
#undef ELTT
#undef ELTP
#undef SEP

}
MRI_SURFACE, MRIS ;

typedef const MRIS MRIS_const;
    // Ideally the MRIS and all the things it points to would be unchangeable via this object but C can't express this concept easily.

void MRISctr(MRIS *mris, int max_vertices, int max_faces, int nvertices, int nfaces);
void MRISdtr(MRIS *mris);
    //
    // These are only way to create a surface and destroy a surface.
    // The destruction frees any child structures.
    //
    // There are functions below for editing the surface by adding vertex positions, edges, face information, etc.
    // There is even one for creating one similar to a subset of another's vertices and faces - MRIScreateWithSimilarTopologyAsSubset

MRI_SURFACE* MRISoverAlloc              (                   int max_vertices, int max_faces, int nvertices, int nfaces) ;
MRI_SURFACE* MRISalloc                  (                                                    int nvertices, int nfaces) ;
    //
    // Allocates an MRIS then calls MRISctr

void MRISfree(MRIS **pmris) ;
    //
    // The only way to delete a surface.  All the substructures are also freed, and the *pmris set to nullptr
    
void MRISreallocVerticesAndFaces(MRIS *mris, int nvertices, int nfaces) ;
    //
    // Used by code that is deforming the surface

void MRISacquireNverticesFrozen(MRIS *mris);
void MRISreleaseNverticesFrozen(MRIS *mris);
    //
    // Used by some code that depends on this number not changing
    // but lots of such dependent code have not been changed to use this
    
// Make and free temp properties for all the vertices
// whilest making sure they are the right size by stopping the nvertices changing when any such exist
//
float* MRISmakeFloatPerVertex(MRIS *mris);
void   MRISfreeFloatPerVertex(MRIS *mris, float** pp);

// The following create a copy of a surface, whilest deleteing some vertices and some faces
// The faces that are kept must not reference any vertices which are not kept.
// There is a function to help generate such a mapping...
//
MRIS* MRIScreateWithSimilarTopologyAsSubset(
    MRIS_const * src,
    size_t       nvertices,     // the mapToDstVno entries must each be less than this
    int const*   mapToDstVno,   // src->nvertices entries, with the entries being -1 (vertex should be ignored) or the vno within the dst surface the face maps to
    size_t       nfaces,        // the mapToDstFno entries must each be less than this
    int const*   mapToDstFno);  // src->nfaces entries, with the entries being -1 (face should be ignored) or the fno within the dst surface the face maps to
    //
    // Used to extract some of the vertices, some of the faces, and to renumber them.
    // mrisurf_deform uses this for several purposes.

MRIS* MRIScreateWithSimilarXYZAsSubset(
    MRIS_const * src,
    size_t       nvertices,     // the mapToDstVno entries must each be less than this
    int const*   mapToDstVno,   // src->nvertices entries, with the entries being -1 (vertex should be ignored) or the vno within the dst surface the face maps to
    size_t       nfaces,        // the mapToDstFno entries must each be less than this
    int const*   mapToDstFno);  // src->nfaces entries, with the entries being -1 (face should be ignored) or the fno within the dst surface the face maps to

MRIS* MRIScreateWithSimilarPropertiesAsSubset(
    MRIS_const * src,
    size_t       nvertices,     // the mapToDstVno entries must each be less than this
    int const*   mapToDstVno,   // src->nvertices entries, with the entries being -1 (vertex should be ignored) or the vno within the dst surface the face maps to
    size_t       nfaces,        // the mapToDstFno entries must each be less than this
    int const*   mapToDstFno);  // src->nfaces entries, with the entries being -1 (face should be ignored) or the fno within the dst surface the face maps to

void MRIScreateSimilarTopologyMapsForNonripped(
    MRIS_const * src,
    size_t     * pnvertices,     // the mapToDstVno entries must each be less than this
    int const* * pmapToDstVno,   // src->nvertices entries, with the entries being -1 (vertex should be ignored) or the vno within the dst surface the face maps to
    size_t     * pnfaces,        // the mapToDstFno entries must each be less than this
    int const* * pmapToDstFno);  // src->nfaces entries, with the entries being -1 (face should be ignored) or the fno within the dst surface the face maps to
    // The caller should free the * pmapToDstVno and * pmapToDstFno

// There are various fields in the VERTEX and FACE and others that are used for many purposes at different
// times, and clashing uses could cause a big problem.  Start working towards a reservation system.
//
typedef enum MRIS_TempAssigned {
    MRIS_TempAssigned_Vertex_marked,
    MRIS_TempAssigned_Vertex_marked2,
    MRIS_TempAssigned__end
} MRIS_TempAssigned;

int  MRIS_acquireTemp      (MRIS* mris, MRIS_TempAssigned temp);                               // save the result to use later to ...
void MRIS_checkAcquiredTemp(MRIS* mris, MRIS_TempAssigned temp, int MRIS_acquireTemp_result);  // ... check that you own it
void MRIS_releaseTemp      (MRIS* mris, MRIS_TempAssigned temp, int MRIS_acquireTemp_result);  // ... be allowed to release it


FaceNormCacheEntry const * getFaceNorm(MRIS const * const mris, int fno);
void setFaceNorm(MRIS const * const mris, int fno, float nx, float ny, float nz);


// Support for writing traces that can be compared across test runs to help find where differences got introduced  
//
typedef struct {
    unsigned long hash;
} MRIS_HASH;

void mrisVertexHash(MRIS_HASH* hash, MRIS const * mris, int vno);

void mris_hash_init (MRIS_HASH* hash, MRIS const * mris);
void mris_hash_add  (MRIS_HASH* hash, MRIS const * mris);
void mris_hash_print(MRIS_HASH const* hash, FILE* file);
void mris_print_hash(FILE* file, MRIS const * mris, const char* prefix, const char* suffix);
void mris_print_diff(FILE* file, MRIS const * lhs, MRIS const * rhs);

// This structs are used with the TESS functions
typedef struct tface_type_
{
  int imnr,i,j,f;
  int num;
  int v[4];
}
tface_type;
typedef struct tvertex_type_
{
  int imnr,i,j;
  int num;
  int f[9];
}
tvertex_type;

typedef struct {
  // Structure for computing info at an equal number of hops from
  // a given vertex (not exactly the same as equal mm distance)
  int cvtx; // center vertex
  int nhops;
  char *hit; // vector to indicate whether vertex has been hit
  int *nperhop; // number of vertices per hop
  int **vtxlist; // list of vertices for each hop
  int *nperhop_alloced; // number of vertices alloced per hop
} SURFHOPLIST;


#define IPFLAG_HVARIABLE                0x0001 /* for parms->flags */
#define IPFLAG_NO_SELF_INT_TEST         0x0002
#define IPFLAG_QUICK                    0x0004 /* sacrifice qualty for speed */
#define IPFLAG_ADD_VERTICES             0x0008
#define IP_USE_CURVATURE                0x0010         /* was 0x0001 !!! */
#define IP_NO_RIGID_ALIGN               0x0020         /* was 0x0002 !!! */
#define IP_RETRY_INTEGRATION            0x0040         /* was 0x0004 !!! */
/* VECTORIAL_REGISTRATION*/
#define IP_USE_MULTIFRAMES              0x0080
#define IP_NO_SULC                      0x0100
#define IP_USE_INFLATED                 0x0200
#define IP_NO_FOLD_REMOVAL              0x0400
#define IP_DONT_ALLOW_FOLDS             0x0800
/* MRIScorrectTopology : topology preserving patch deformation */
#define IPFLAG_PRESERVE_TOPOLOGY_CONVEXHULL 0x1000 /* apply topology
preserving gradient */
#define IPFLAG_PRESERVE_SPHERICAL_POSITIVE_AREA 0x2000 /* apply gradients
that preserve
positive areas */
#define IPFLAG_MAXIMIZE_SPHERICAL_POSITIVE_AREA 0x4000  /* apply  gradients
that will
maximize the
positive areas */
#define IPFLAG_NOSCALE_TOL            0x8000   // don't scale tol with navgs
#define IPFLAG_FORCE_GRADIENT_OUT    0x10000
#define IPFLAG_FORCE_GRADIENT_IN     0x20000
#define IPFLAG_FIND_FIRST_WM_PEAK    0x40000  // for Matt Glasser/David Van Essen

#define INTEGRATE_LINE_MINIMIZE    0  /* use quadratic fit */
#define INTEGRATE_MOMENTUM         1
#define INTEGRATE_ADAPTIVE         2
#define INTEGRATE_LM_SEARCH        3  /* binary search for minimum */

#define MRIS_SURFACE               0
#define MRIS_PATCH                 1
#define MRIS_CUT                   MRIS_PATCH
#define MRIS_PLANE                 2
#define MRIS_ELLIPSOID             3
#define MRIS_SPHERE                4
#define MRIS_PARAMETERIZED_SPHERE  5
#define MRIS_RIGID_BODY            6
#define MRIS_SPHERICAL_PATCH       7
#define MRIS_UNORIENTED_SPHERE     8
#define MRIS_PIAL_SURFACE          9

// different Hausdorff distance modes
#define HDIST_MODE_SYMMETRIC_MEAN 0
double MRIScomputeHausdorffDistance(MRI_SURFACE *mris, int mode) ;


/*
  the following structure is built after the surface has been morphed
  into a parameterizable surface (such as a sphere or an ellipsoid). It
  contains relevant data (such as curvature or sulc) which represented
  in the plane of the parameterization. This then allows operations such
  as convolution to be carried out directly in the parameterization space
*/

/*
   define it for a sphere, the first axis (x or u) goes from 0 to pi,
   the second axis (y or v) goes from 0 to 2pi.
*/
#define PHI_MAX                   M_PI
#define U_MAX                     PHI_MAX
#define X_DIM(mrisp)              (mrisp->Ip->cols)
#define U_DIM(mrisp)              X_DIM(mrisp)
#define PHI_DIM(mrisp)            X_DIM(mrisp)
#define PHI_MAX_INDEX(mrisp)      (X_DIM(mrisp)-1)
#define U_MAX_INDEX(mrisp)        (PHI_MAX_INDEX(mrisp))

#define THETA_MAX                 (2.0*M_PI)
#define Y_DIM(mrisp)              (mrisp->Ip->rows)
#define V_DIM(mrisp)              (Y_DIM(mrisp))
#define THETA_DIM(mrisp)          (Y_DIM(mrisp))
#define THETA_MAX_INDEX(mrisp)    (Y_DIM(mrisp)-1)
#define V_MAX_INDEX(mrisp)        (THETA_MAX_INDEX(mrisp))

#include "image.h" // IMAGE
typedef struct
{
  MRI_SURFACE  *mris ;        /* surface it came from (if any) */
  IMAGE        *Ip ;
#if 0
  /* 2-d array of curvature, or sulc in parms */
  float        data[X_DIM][Y_DIM] ;
  float        distances[X_DIM][Y_DIM] ;
  int          vertices[X_DIM][Y_DIM] ;   /* vertex numbers */
#endif
  float        sigma ;                    /* blurring scale */
  float        radius ;
  float        scale ;
}
MRI_SURFACE_PARAMETERIZATION, MRI_SP ;

#define L_ANGLE              0.25f /*was 0.01*/ /* coefficient of angle term */
#define L_AREA               1.0f    /* coefficient of angle term */
#define N_AVERAGES           4096
#define WRITE_ITERATIONS     10
#define NITERATIONS          1
#define NO_PROJECTION        0
#define PROJECT_ELLIPSOID    1
#define ELLIPSOID_PROJECTION PROJECT_ELLIPSOID
#define PROJECT_PLANE        2
#define PLANAR_PROJECTION    PROJECT_PLANE
#define PROJECT_SPHERE       3

#define MOMENTUM             0.8
#define EPSILON              0.25

#define TOL                  1e-6  /* minimum error tolerance for unfolding */
#define DELTA_T              0.1

/* VECTORIAL_REGISTRATION */
#include "field_code.h"

typedef struct INTEGRATION_PARMS
{
  double  tol ;               /* tolerance for terminating a step */
  float   l_angle ;           /* coefficient of angle term */
  float   l_pangle ;          /* coefficient of "positive" angle term - only penalize narrower angles */
  float   l_area ;            /* coefficient of (negative) area term */
  float   l_parea ;           /* coefficient of (all) area term */
  float   l_nlarea ;          /* coefficient of nonlinear area term */
  float   l_nldist ;          /* coefficient of nonlinear distance term */
  float   l_thick_normal ;    // coefficient to keep vector field close to surface normal
  float   l_thick_spring ;    // coefficient to keep vector field spaced on pial surface
  float   l_ashburner_triangle ;   // coefficient for (Ashburner, 1999) invertibility term
  float   l_ashburner_lambda ;  // amount of regularization
  float   l_corr ;            /* coefficient of correlation term */
  float   l_ocorr ;           // overlay correlation weight
  float   l_pcorr ;           /* polar correlation for rigid body */
  float   l_curv ;            /* coefficient of curvature term */
  float   l_norm ;            /* coefficient of surface self-intersection */
  float   l_scurv ;           /* coefficient of curvature term */
  float   l_lap ;             // coefficient of laplacian term
  float   l_link ;            /* coefficient of link term to keep
                                 white and pial vertices approximately
                                 along normal */
  float   l_spring ;          /* coefficient of spring term */
  float   l_nlspring ;        /* nonlinear spring term */
  float   l_max_spring ;      /* draws vertices towards the most distant neighbor */
  float   l_spring_norm ;     /* coefficient of normalize spring term */
  float   l_tspring ;         /* coefficient of tangential spring term */
  float   l_nltspring ;       /* coefficient of nonlinear tangential spring term */
  float   l_nspring ;         /* coefficient of normal spring term */
  float   l_repulse ;         /* repulsize force on tessellation */
  float   l_repulse_ratio ;   /* repulsize force on tessellation */
  float   l_boundary ;        /* coefficient of boundary term */
  float   l_dist ;            /* coefficient of distance term */
  float   l_location ;        // target location term
  float   l_neg ;
  float   l_intensity ;       /* for settling surface at a specified val */
  float   l_sphere ;          /* for expanding the surface to a sphere */
  float   l_expand ;          /* for uniformly expanding the surface */
  float   l_grad ;            /* gradient term */
  float   l_convex ;          /* convexity term */
  float   l_tsmooth ;         /* thickness smoothness term */
  float   l_surf_repulse ;    /* repulsive orig surface (for white->pial) */
  float   l_osurf_repulse ;   /* repulsive outer surface (for layer IV) */
  float   l_external ;        /* external (user-defined) coefficient */
  float   l_thick_parallel ;  // term that encourages thickness vectors to be parallel
  float   l_thick_min ;       // term that encourages thickness vectors to be minimal
  float   l_shrinkwrap ;      /* move in if MRI=0 and out otherwise */
  float   l_expandwrap ;      /* move out */
  float   l_unfold ;          /* move inwards along normal */
  float   l_dura ;            // move away from dura
  float   l_histo ;           // increase the likelihood of the entire volume given the surfaces
  double  l_map ;             // for MAP deformation
  double  l_map2d ;             // for 2D MAP deformation (intensity x distance)
  double  dura_thresh ;
  MRI     *mri_dura ;         /* ratio of early to late echo -
                                         dura shows up bright */
  int     n_averages ;        /* # of averages */
  int     min_averages ;
  int     first_pass_averages;// # of averages to use in the first pass
  int     nbhd_size ;
  int     max_nbrs ;
  int     write_iterations ;  /* # of iterations between saving movies */
  char    base_name[STRLEN] ; /* base name of movie files */
  int     projection ;        /* what kind of projection to do */
  int     niterations ;       /* max # of time steps */
  float   a ;
  float   b ;
  float   c ;                 /* ellipsoid parameters */
  int     start_t ;           /* starting time step */
  int     t ;                 /* current time */
  
  FILE    * const fp ;        /* for logging results, write by calling INTEGRATION_PARMS_<various> functions */
  
  float   Hdesired ;          /* desired (mean) curvature */
  int     integration_type ;  /* line minimation or momentum */
  double  momentum ;
  double  dt ;                /* current time step (for momentum only) */
  double  base_dt ;           /* base time step (for momentum only) */
  int     flags ;
  double  dt_increase ;       /* rate at which time step increases */
  double  dt_decrease ;       /* rate at which time step decreases */
  double  error_ratio ;       /* ratio at which to undo step */
  double  epsilon ;           /* constant in Sethian inflation */
  double  desired_rms_height; /* desired height above tangent plane */
  double  starting_sse ;
  double  ending_sse ;
  double  scale ;             /* scale current distances to mimic spring */
  MRI_SP  *mrisp ;            /* parameterization  of this surface */
  int     frame_no ;          /* current frame in template parameterization */
  MRI_SP  *mrisp_template ;   /* parameterization of canonical surface */
  MRI_SP  *mrisp_blurred_template ; /* parameterization of canonical
                                       surface convolve with Gaussian */
  double  area_coef_scale ;
  float   sigma ;             /* blurring scale */

  /* VECTORIAL_REGISTRATION
     The average template mrisp is assumed to be composed of several
     different fields (look in 'field_code.h').
     MRISvectorRegistration will use 'nfields' fields,
     with their corresponding
     location at fields[n].frames, 0 <= n< nfields in the mrisp structure.
     The field code is in fields[n].field (look in 'field_code.h').
     Corresponding correlation terms are in
     fields[n].l_corrs and fields[n].l_pcorrs.
     MRISvectorRegistration will use the structure VALS_VP in v->vp
  */
#define MAX_NUMBER_OF_FIELDS_IN_VECTORIAL_REGISTRATION  50
#define MNOFIV MAX_NUMBER_OF_FIELDS_IN_VECTORIAL_REGISTRATION
  int nfields;                  /* the number of fields in mrisp */
  FIELD_LABEL fields[MNOFIV];   /* information for each field */
  /*    int     ncorrs;            /\* the number of fields in mrisp *\/ */
  /*    int     corrfields[MNOFIV];/\* field code (see below) *\/ */
  /*    int     frames[MNOFIV];    /\* corresponding frame in mrisp *\/  */
  /*    int     types[MNOFIV];     /\* the field type
    (default,distance field...) *\/ */
  /*    float   l_corrs[MNOFIV];   /\* correlation coefficient *\/ */
  /*    float   l_pcorrs[MNOFIV];  /\* polar correlation coefficient *\/ */
  /*    float   sses[MNOFIV];      /\* corresponding sse *\/ */

  MRI     *mri_brain ;        /* for settling surfaace to e.g. g/w border */
  MRI     *mri_smooth ;       /* smoothed version of mri_brain */
  void    *user_parms ;       /* arbitrary spot for user to put stuff */
  MRI     *mri_dist ;         /* distance map for repulsion term */
  float   target_radius ;
  int     ignore_energy ;     // when no valid energy func availb...integrate
  int     check_tol ;         // to avoid changing mris_make_surfaces
  char    *overlay_dir;       // subject/overlay_dir/parms->fields[n].fname
  int     nsurfaces ;         // if 0 use default
  MRI     *mri_ll ;           // log-likelihood image
  double  rmin ;
  double  rmax ;              // for nonlinear spring term
  int     var_smoothness ;    // for space/time varying weights on
                              // metric distortion/likelihood
  float   *vsmoothness ;      // variable smoothness coefficients
                              // (one per vertex)
  float   *dist_error ;       // the values for each vertex of the
                              // various error type
  float   *area_error ;       // needed for the variable smoothness
                              // gradient calculation
  float   *geometry_error ;
  int     which_norm ;        // mean or median normalization
  int     abs_norm ;
  int     grad_dir ;          // use this instead of gradient direction
  int     fill_interior ;     // use filled interior to constrain gradient to not leave surface
  double  rms ;
  int     complete_dist_mat ; //whether to sample or use complete dist mat
  int     nsubjects ;
  int     nlabels ;
  void       **mht_array;
  MRI_SURFACE **mris_array ;
  MRI_SURFACE *mris_ico ;     // for sampling from a spherical template
  void         *mht ;        // hash table for surface vertex/face lookup
  int          smooth_averages ;
  int          ico_order ;  // which icosahedron to use
  int          remove_neg ;
  MRI          *mri_hires ;
  MRI          *mri_hires_smooth ;
  MRI          *mri_vno ;
  MRI          *mri_template ;
  int          which_surface ;
  float        trinarize_thresh;   // for trinarizing curvature in registration
  int          smooth_intersections ;  // run soap bubble smoothing during surface positioning
  int          uncompress ;            // run code to remove compressions in tessellation
  double       min_dist ;
  HISTOGRAM    *h_wm ;
  HISTOGRAM    *h_gm ;
  HISTOGRAM    *h_nonbrain ;
  MRI          *mri_labels ;   // hires labeling of interior of WM, GM and nonbrain
  MRI          *mri_white ;
  MRI          *mri_aseg ;
  HISTOGRAM    *hwm ;
  HISTOGRAM    *hgm ;
  HISTOGRAM    *hout ;

  HISTOGRAM2D  *h2d_wm ;
  HISTOGRAM2D  *h2d_gm ;
  HISTOGRAM2D  *h2d_out ;
  HISTOGRAM2D  *h2d ;
  MRI          *mri_volume_fractions ;  // the partial volume fractions associated with the boundaries in this mris
  MRI          *mri_dtrans ;   // distance to surface
  float        resolution ;  // at which to compute distance transforms and such
  double       target_intensity ;
  double       stressthresh ;
  int          explode_flag ;
  
#ifdef __cplusplus
  INTEGRATION_PARMS() : fp(NULL) {}
  INTEGRATION_PARMS(FILE* file) : fp(file) {}
#endif 
  
}
INTEGRATION_PARMS ;

void INTEGRATION_PARMS_copy   (INTEGRATION_PARMS* dst, INTEGRATION_PARMS const * src);

void INTEGRATION_PARMS_setFp  (INTEGRATION_PARMS* parms, FILE* file);
void INTEGRATION_PARMS_openFp (INTEGRATION_PARMS* parms, const char* name, const char* mode);
void INTEGRATION_PARMS_closeFp(INTEGRATION_PARMS* parms);
void INTEGRATION_PARMS_copyFp (INTEGRATION_PARMS* dst, INTEGRATION_PARMS const * src);


extern double (*gMRISexternalGradient)                   (MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
extern double (*gMRISexternalSSE)                        (MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
extern double (*gMRISexternalRMS)                        (MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
extern int    (*gMRISexternalTimestep)                   (MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
extern int    (*gMRISexternalRipVertices)                (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
extern int    (*gMRISexternalClearSSEStatus)             (MRI_SURFACE *mris) ;
extern int    (*gMRISexternalReduceSSEIncreasedGradients)(MRI_SURFACE *mris, double pct) ;

// These are for backwards compatibility for when we don't want to fix
// the vertex area. Default is to fix, but this can be changed by setting to 0.
//
int MRISsetFixVertexAreaValue(int value);
int MRISgetFixVertexAreaValue(void);

/* The following structure is used in MRISvectorRegistration
   The original values are loaded in orig_vals
   The used values are in vals (for instance the  blurred orig_vals)
*/
typedef struct
{
  float *vals;
  float *orig_vals;
  int nvals;
}
VALS_VP;

const char *MRISurfSrcVersion(void);

/* new functions */
int MRISvectorRegister(MRI_SURFACE *mris,
                       MRI_SP *mrisp_template,
                       INTEGRATION_PARMS *parms, int max_passes,
                       float min_degrees, float max_degrees,
                       int nangles);
int MRISrigidBodyAlignVectorLocal(MRI_SURFACE *mris,
                                  INTEGRATION_PARMS *old_parms);
int MRISrigidBodyAlignVectorGlobal(MRI_SURFACE *mris,
                                   INTEGRATION_PARMS *parms,
                                   float min_degrees,
                                   float max_degrees, int nangles);
void MRISnormalizeField(MRIS *mris , int distance_field, int which_norm);
int  MRISsmoothCurvatures(MRI_SURFACE *mris, int niterations) ;
void MRISsetCurvaturesToValues(MRIS *mris,int fno);
void MRISsetCurvaturesToOrigValues(MRIS *mris,int fno);
void MRISsetOrigValuesToCurvatures(MRIS *mris,int fno);
void MRISsetOrigValuesToValues(MRIS *mris,int fno);

MRI *MRISbinarizeVolume(MRI_SURFACE *mris,
                        MRI_REGION *region,
                        float resolution,
                        float distance_from_surface);

int MRISfindClosestCanonicalVertex(MRI_SURFACE *mris, float x, float y,
                                   float z) ;
int MRISfindClosestOriginalVertex(MRI_SURFACE *mris, float x, float y,
                                  float z) ;
double MRIScomputeWhiteVolume(MRI_SURFACE *mris, 
                              MRI *mri_aseg, 
                              double resolution) ;
MRI *MRISfillWhiteMatterInterior(MRI_SURFACE *mris,
                                 MRI *mri_aseg,
                                 MRI *mri_filled,
                                 double resolution,
                                 int wm_val, int gm_val, int csf_val);
int MRISfindClosestWhiteVertex(MRI_SURFACE *mris, float x, float y,
                               float z) ;
int MRISfindClosestVertex(MRI_SURFACE *mris,
                          float x, float y, float z,
                          float *dmin, int which_vertices);
double MRIScomputeSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
double MRIScomputeSSEExternal(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                              double *ext_sse) ;
double       MRIScomputeCorrelationError(MRI_SURFACE *mris,
    MRI_SP *mrisp_template, int fno) ;
int          MRISallocExtraGradients(MRI_SURFACE *mris) ;
MRI_SURFACE  *MRISread(const char *fname) ;
MRI_SURFACE  *MRISreadOverAlloc(const char *fname, double nVFMultiplier) ;
MRI_SURFACE  *MRISfastRead(const char *fname) ;
int          MRISreadOriginalProperties(MRI_SURFACE *mris,const  char *sname) ;
int          MRISreadCanonicalCoordinates(MRI_SURFACE *mris,
                                          const  char *sname) ;
int          MRISreadInflatedCoordinates(MRI_SURFACE *mris,
                                         const  char *sname) ;
int          MRISreadFlattenedCoordinates(MRI_SURFACE *mris,
                                          const  char *sname) ;
int          MRISreadPialCoordinates(MRI_SURFACE *mris, const char *sname);
int          MRISreadWhiteCoordinates(MRI_SURFACE *mris,const  char *sname) ;
int          MRIScomputeCanonicalCoordinates(MRI_SURFACE *mris) ;
int          MRIScanonicalToWorld(MRI_SURFACE *mris, double phi, double theta,
                                  double *pxw, double *pyw, double *pzw) ;
int          MRISreadPatch(MRI_SURFACE *mris,const  char *pname) ;
int          MRISreadPatchNoRemove(MRI_SURFACE *mris,const  char *pname) ;
int          MRISreadTriangleProperties(MRI_SURFACE *mris,
                                        const  char *mris_fname) ;
int          MRISreadBinaryCurvature(MRI_SURFACE *mris,
                                     const  char *mris_fname) ;
int          MRISreadCurvatureFile(MRI_SURFACE *mris,const char *fname) ;
float        *MRISreadNewCurvatureVector(MRI_SURFACE *mris,
                                         const  char *sname) ;
int          MRISreadNewCurvatureIntoArray(const char *fname,
                                           int in_array_size,
                                           float** out_array) ;
float        *MRISreadCurvatureVector(MRI_SURFACE *mris,const  char *sname) ;
int          MRISreadCurvatureIntoArray(const char *fname,
                                        int in_array_size,
                                        float** out_array) ;
int          MRISreadFloatFile(MRI_SURFACE *mris,const char *fname) ;
#define MRISreadCurvature MRISreadCurvatureFile

MRI *MRISloadSurfVals(const char *srcvalfile,
                      const char *typestring,
                      MRI_SURFACE *Surf,
                      const char *subject,
                      const char *hemi,
                      const char *subjectsdir);
int          MRISreadValues(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadValuesIntoArray(const char *fname,
                                     int in_array_size,
                                     float** out_array) ;
int          MRISreadAnnotation(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteVertexLocations(MRI_SURFACE *mris, char *fname, int which_vertices) ;
int          MRISimportVertexCoords(MRI_SURFACE *mris, float *locations[3], int which_vertices);
int          MRISwriteAnnotation(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadAnnotationIntoArray(const char *fname,
                                         int in_array_size,
                                         int** out_array);
int          MRISreadCTABFromAnnotationIfPresent(const char *fname,
                                                 COLOR_TABLE** out_table);
int          MRISisCTABPresentInAnnotation(const char *fname, int* present);
int          MRISreadValuesBak(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadImagValues(MRI_SURFACE *mris,const  char *fname) ;
int          MRIScopyImagValuesToValues(MRI_SURFACE *mris) ;
int          MRIScopyMarksToAnnotation(MRI_SURFACE *mris) ;
int          MRIScopyValsToAnnotations(MRI_SURFACE *mris) ;
int          MRIScopyValuesToImagValues(MRI_SURFACE *mris) ;
int          MRIScopyStatsToValues(MRI_SURFACE *mris) ;
int          MRIScopyStatsFromValues(MRI_SURFACE *mris) ;


int MRISsetCroppedToZero(MRI_SURFACE *mris) ;
int MRIScopyFromCropped(MRI_SURFACE *mris, int which) ;
int MRIScopyToCropped(MRI_SURFACE *mris, int which) ;
int MRISwriteCropped(MRI_SURFACE *mris, char *fname) ;

int          MRIScopyValToVal2(MRI_SURFACE *mris) ;
int          MRIScopyValToVal2Bak(MRI_SURFACE *mris) ;
int          MRIScopyValToValBak(MRI_SURFACE *mris) ;
int          MRISsqrtVal(MRI_SURFACE *mris) ;
int          MRISmulVal(MRI_SURFACE *mris, float mul) ;

int          MRISwrite(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteAscii(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteWhiteNormals(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteNormalsAscii(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadNormals(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteNormals(MRI_SURFACE *mris,const  char *fname) ;
int mrisNormalFace(MRIS *mris, int fac, int n, float norm[]);
int          MRISwritePrincipalDirection(MRI_SURFACE *mris, int dir_index, const  char *fname) ;
int          MRISwriteVTK(MRI_SURFACE *mris,const  char *fname);
int          MRISwriteCurvVTK(MRI_SURFACE *mris, const char *fname);
int          MRISwriteGeo(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteICO(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteSTL(MRI_SURFACE *mris, const char *fname) ;
int          MRISwritePatchAscii(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteDists(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteCurvature(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadNewCurvatureFile(MRI_SURFACE *mris,const  char *fname) ;
int          MRISrectifyCurvature(MRI_SURFACE *mris) ;
#define NORM_NONE  -1
#define NORM_MEAN   0
#define NORM_MEDIAN 1
#define NORM_MAX    2

int          MRISnormalizeCurvature(MRI_SURFACE *mris, int norm_type) ;
int          MRISnormalizeCurvatureVariance(MRI_SURFACE *mris) ;
int          MRISzeroMeanCurvature(MRI_SURFACE *mris) ;
int          MRISnonmaxSuppress(MRI_SURFACE *mris) ;
int          MRISscaleCurvatures(MRI_SURFACE *mris,
                                 float min_curv, float max_curv) ;
int          MRISwriteAreaError(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteAngleError(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwritePatch(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteValues(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteD(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteCurvatureToWFile(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteTriangleProperties(MRI_SURFACE *mris,
                                         const char *mris_fname);
int          MRISaverageCurvatures(MRI_SURFACE *mris, int navgs) ;
int          MRISminFilterCurvatures(MRI_SURFACE *mris, int niter) ;
int          MRISmaxFilterCurvatures(MRI_SURFACE *mris, int niter) ;
int          MRISaverageMarkedCurvatures(MRI_SURFACE *mris, int navgs) ;
double       MRIScomputeAverageCurvature(MRI_SURFACE *mris, double *psigma) ;
int          MRISaverageVertexPositions(MRI_SURFACE *mris, int navgs) ;
int          MRIScomputeNormal(MRIS *mris, int which, int vno,
                               double *pnx, double *pny, double *pnz) ;

int   MRISintegrate(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int n_avgs);
int   mrisLogIntegrationParms(FILE *fp, MRI_SURFACE *mris,
			      INTEGRATION_PARMS *parms) ;
MRI_SURFACE  *MRISprojectOntoSphere(MRI_SURFACE *mris_src,
                                    MRI_SURFACE *mris_dst, double r) ;
MRI_SURFACE  *MRISprojectOntoEllipsoid(MRI_SURFACE *mris_src,
                                       MRI_SURFACE *mris_dst,
                                       float a, float b, float c) ;
int          MRISsampleDistances(MRI_SURFACE *mris, int *nbr_count,int n_nbrs);
int          MRISsampleAtEachDistance(MRI_SURFACE *mris, int nbhd_size,
                                      int nbrs_per_distance) ;
int          MRISscaleDistances(MRI_SURFACE *mris, float scale) ;
MRI_SURFACE  *MRISradialProjectOntoEllipsoid(MRI_SURFACE *mris_src,
    MRI_SURFACE *mris_dst,
    float a, float b, float c);
    
int MRISorigVertexToVoxel(MRI_SURFACE *,
                          VERTEX *v,
                          MRI *mri,
                          double *pxv, double *pyv, double *pzv) ;
int MRISwhiteVertexToVoxel(MRI_SURFACE *,
                           VERTEX *v,
                           MRI *mri,
                           double *pxv, double *pyv, double *pzv) ;
int          MRISvertexToVoxel(MRI_SURFACE *, VERTEX *v, MRI *mri,
                               double *pxv, double *pyv,
                               double *pzv) ;
int          MRISvertexCoordToVoxel(MRI_SURFACE *, VERTEX *v, MRI *mri,
                                    int coord,
                                    double *pxv, double *pyv,
                                    double *pzv) ;
#if 0
int          MRISworldToTalairachVoxel(MRI_SURFACE *mris, MRI *mri,
                                       double xw, double yw, double zw,
                                       double *pxv, double *pyv, double *pzv) ;
#endif

int          MRISsurfaceRASToVoxel(MRI_SURFACE *mris, MRI *mri, double r, 
                                   double a, double s, 
                                   double *px, double *py, double *pz) ;
                                   
// THE FOLLOWING IS NOT THREAD SAFE!
int          MRISsurfaceRASToVoxelCached(MRI_SURFACE *mris,
                                         MRI *mri,
                                         double r, double a, double s, 
                                         double *px, double *py, double *pz) ;


typedef struct MRIS_SurfRAS2VoxelMap {
    MRI*    mri;                            // was held in mris->mri_sras2vox
    MATRIX* sras2vox;                       // was held in mris->m_sras2vox
    VECTOR * volatile v1[_MAX_FS_THREADS];  // used to avoid repeated allocations
    VECTOR * volatile v2[_MAX_FS_THREADS];  // used to avoid repeated allocations
} MRIS_SurfRAS2VoxelMap;

void MRIS_useRAS2VoxelMap(MRIS_SurfRAS2VoxelMap * cache_nonconst,   // accesses cache thread safely
        MRI const * const mri,
        double r, double a, double s, double *px, double *py, double *pz);
    
void MRIS_loadRAS2VoxelMap(MRIS_SurfRAS2VoxelMap* cache,            // not thread safe
        MRI const * const mri, MRI_SURFACE const * const mris);
void MRIS_unloadRAS2VoxelMap(MRIS_SurfRAS2VoxelMap* cache);         // not thread safe

MRIS_SurfRAS2VoxelMap* MRIS_makeRAS2VoxelMap(                       // not thread safe
        MRI const * const mri, MRI_SURFACE const * const mris);
void MRIS_freeRAS2VoxelMap(MRIS_SurfRAS2VoxelMap** const cachePtr); // not thread safe

// these are the inverse of the previous two
int          MRISsurfaceRASFromVoxel(MRI_SURFACE *mris, MRI *mri, 
                                     double x, double y, double z, 
                                     double *pr, double *pa, double *ps) ;
int          MRISsurfaceRASFromVoxelCached(MRI_SURFACE *mris, MRI *mri, 
                                           double x, double y, double z, 
                                           double *pr, double *pa, double *ps);

int          MRISsurfaceRASToTalairachVoxel(MRI_SURFACE *mris, MRI *mri,
    double xw, double yw, double zw,
    double *pxv, double *pyv, double *pzv) ;

int          MRIStalairachToVertex(MRI_SURFACE *mris,
                                   double xt, double yt, double zt) ;
int           MRIScanonicalToVertex(MRI_SURFACE *mris,
                                    double phi,
                                    double theta) ;
MRI_SURFACE  *MRIStalairachTransform(MRI_SURFACE *mris_src,
                                     MRI_SURFACE *mris_dst);
MRI_SURFACE  *MRISunfold(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                         int max_passes) ;
MRI_SURFACE  *MRISquickSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                              int max_passes) ;
int          MRISregister(MRI_SURFACE *mris,
                          MRI_SP *mrisp_template,
                          INTEGRATION_PARMS *parms,
                          int max_passes,
                          float min_degrees,
                          float max_degrees,
                          int nangles) ;
int          MRISrigidBodyAlignLocal(MRI_SURFACE *mris,
                                     INTEGRATION_PARMS *parms) ;
int          MRISrigidBodyAlignGlobal(MRI_SURFACE *mris,
                                      INTEGRATION_PARMS *parms,
                                      float min_degrees,
                                      float max_degrees,
                                      int nangles) ;
MRI_SURFACE  *MRISunfoldOnSphere(MRI_SURFACE *mris,
                                 INTEGRATION_PARMS *parms,
                                 int max_passes);
MRI_SURFACE  *MRISflatten(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
MRI_SURFACE  *MRISremoveNegativeVertices(MRI_SURFACE *mris,
    INTEGRATION_PARMS *parms,
    int min_neg, float min_neg_pct) ;
int          MRIScomputeFaceAreas(MRI_SURFACE *mris) ;
int          MRISupdateEllipsoidSurface(MRI_SURFACE *mris) ;
MRI_SURFACE  *MRISrotate(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst,
                         float alpha, float beta, float gamma) ;

MRI          *MRISwriteIntoVolume(MRI_SURFACE *mris, MRI *mri, int which) ;
MRI_SURFACE  *MRISreadFromVolume(MRI *mri, MRI_SURFACE *mris, int which) ;



int          MRIScomputeTriangleProperties(MRI_SURFACE *mris) ;
int          MRISpaintVolume(MRI_SURFACE *mris, LTA *lta, MRI *mri) ;
int          MRISsampleStatVolume(MRI_SURFACE *mris, void *sv,int time,
                                  int use_talairach_xform);

int          MRISflattenPatch(MRI_SURFACE *mris) ;
int          MRISflattenPatchRandomly(MRI_SURFACE *mris) ;
int          MRIScomputeMeanCurvature(MRI_SURFACE *mris) ;

int          MRIScomputeEulerNumber(MRI_SURFACE *mris, int *pnvertices,
                                    int *pnfaces, int *pnedges) ;
int          MRIStopologicalDefectIndex(MRI_SURFACE *mris) ;
int          MRISremoveTopologicalDefects(MRI_SURFACE *mris,
                                          float curv_thresh);

int          MRIScomputeSecondFundamentalFormAtVertex(MRI_SURFACE *mris,
                                                      int vno,
                                                      int *vertices,
                                                      int vnum) ;
int          MRIScomputeSecondFundamentalFormThresholded(MRI_SURFACE *mris,
                                                         double thresh) ;
int          MRIScomputeSecondFundamentalForm(MRI_SURFACE *mris) ;
int          MRISuseCurvatureDifference(MRI_SURFACE *mris) ;
int          MRISuseCurvatureStretch(MRI_SURFACE *mris) ;
int          MRISuseCurvatureMax(MRI_SURFACE *mris) ;
int          MRISuseCurvatureMin(MRI_SURFACE *mris) ;
int          MRISuseNegCurvature(MRI_SURFACE *mris) ;
int          MRISuseAreaErrors(MRI_SURFACE *mris) ;
int          MRISuseGaussianCurvature(MRI_SURFACE *mris) ;
int          MRISclearCurvature(MRI_SURFACE *mris) ;
void         MRISclearD(MRIS *mris) ;
int          MRISusePrincipalCurvature(MRI_SURFACE *mris) ;
int          MRISuseMeanCurvature(MRI_SURFACE *mris) ;
int          MRISuseK1Curvature(MRI_SURFACE *mris);
int          MRISuseK2Curvature(MRI_SURFACE *mris);
int          MRISusePrincipalCurvatureFunction(MRI_SURFACE*		pmris, 
                                               float 			(*f)(float k1, float k2));


int          MRIScomputeCurvatureIndices(MRI_SURFACE *mris,
                                         double *pici, double *pfi) ;
int          MRISuseCurvatureRatio(MRI_SURFACE *mris) ;
int          MRISuseCurvatureContrast(MRI_SURFACE *mris) ;

double MRIScomputeFolding(MRI_SURFACE *mris) ;
double MRISavgVetexRadius(MRIS *Surf, double *StdDev);

int          MRISprojectOntoCylinder(MRI_SURFACE *mris, float radius) ;
double       MRISaverageRadius(MRI_SURFACE *mris) ;
double       MRISaverageCanonicalRadius(MRI_SURFACE *mris) ;
double       MRISmaxRadius(MRI_SURFACE *mris) ;
int          MRISinflateToSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
int          MRISinflateBrain(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
double       MRISrmsTPHeight(MRI_SURFACE *mris) ;
double       MRIStotalVariation(MRI_SURFACE *mris) ;
double       MRIStotalVariationDifference(MRI_SURFACE *mris) ;
double       MRIScurvatureError(MRI_SURFACE *mris, double Kd) ;
MRI_SURFACE  *MRISscaleBrain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst,
                             float scale) ;
int          MRISstoreMetricProperties(MRI_SURFACE *mris) ;
int          MRISrestoreMetricProperties(MRI_SURFACE *mris) ;
double       MRIScomputeAnalyticDistanceError(MRI_SURFACE *mris, int which,
    FILE *fp);
int          MRISzeroNegativeAreas(MRI_SURFACE *mris) ;
int          MRIScountNegativeTriangles(MRI_SURFACE *mris) ;
int          MRIScountMarked(MRI_SURFACE *mris, int mark_threshold) ;
int MRIScountRipped(MRIS *mris);
int MRIScountAllMarked(MRIS *mris);
int          MRIScountTotalNeighbors(MRI_SURFACE *mris, int nsize) ;
int          MRISstoreMeanCurvature(MRI_SURFACE *mris) ;
int          MRISreadTetherFile(MRI_SURFACE *mris,
                                const char *fname,
                                float radius) ;
int          MRISreadVertexPositions(MRI_SURFACE *mris,const  char *fname) ;
int          MRISspringTermWithGaussianCurvature(MRI_SURFACE *mris,
                                                 double gaussian_norm,
                                                 double l_spring) ;
int          MRISmarkedSpringTerm(MRI_SURFACE *mris, double l_spring) ;
double       MRISmomentumTimeStep(MRI_SURFACE *mris,
                                  float momentum,
                                  float dt,
                                  float tol,
                                  float n_averages) ;
int          MRISapplyGradient(MRI_SURFACE *mris, double dt) ;
int          MRIScomputeNormals(MRI_SURFACE *mris) ;
int          MRIScomputeSurfaceNormals(MRI_SURFACE *mris,
                                       int which,
                                       int navgs) ;
int          MRIScomputeMetricProperties(MRI_SURFACE *mris) ;
double       MRISrescaleMetricProperties(MRIS *surf);
int          MRISrestoreOldPositions(MRI_SURFACE *mris) ;
int          MRISstoreCurrentPositions(MRI_SURFACE *mris) ;
int          MRISupdateSurface(MRI_SURFACE *mris) ;
double       MRISpercentDistanceError(MRI_SURFACE *mris) ;
int          MRISscaleBrainArea(MRI_SURFACE *mris) ;


int          MRISPsetFrameVal(MRI_SP *mrisp, int frame, float val) ;
MRI_SP       *MRISPcombine(MRI_SP *mrisp, MRI_SP *mrisp_template, int fno);
MRI_SP       *MRISPaccumulate(MRI_SP *mrisp, MRI_SP *mrisp_template, int fno);
int          MRISPcoordinate(MRI_SP *mrisp, float x, float y, float z,
                             int *pu, int *pv) ;

typedef struct MRISPfunctionValResultForAlpha {
    double curr;
    double next;
} MRISPfunctionValResultForAlpha;

void MRISPfunctionVal_radiusR(                                                      // returns the value that would be stored in resultsForEachFno[0] for fnoLo
                              MRI_SURFACE_PARAMETERIZATION *mrisp,  
                              MRISPfunctionValResultForAlpha* resultsForEachAlpha,  // must be numAlphas elements
                              MRI_SURFACE *mris,
                              float r, float x, float y, float z, 
                              int fno, bool getNextAlso,                            // always fills in resultsForEachAlpha.curr for fno, optionally fills in .next for fno+1
                              const float* alphas, float numAlphas,                 // rotate x,y,z by these alphas (radians) and get the values
                              bool trace) ;                                         // note: this rotation is around the z axis, hence z does not change
                             
double       MRISPfunctionValTraceable(MRI_SURFACE_PARAMETERIZATION *mrisp,
                              MRI_SURFACE *mris,
                              float x, float y, float z, int fno, bool trace) ;
double       MRISPfunctionVal(MRI_SURFACE_PARAMETERIZATION *mrisp,
                              MRI_SURFACE *mris,
                              float x, float y, float z, int fno) ;
                              
MRI_SP       *MRIStoParameterizationBarycentric(MRI_SURFACE *mris, MRI_SP *mrisp,
						float scale, int fno) ;
MRI_SURFACE  *MRISfromParameterizationBarycentric(MRI_SP *mrisp, MRI_SURFACE *mris,
						  int fno) ;

MRI_SP       *MRIStoParameterization(MRI_SURFACE *mris, MRI_SP *mrisp,
                                     float scale, int fno) ;
MRI_SURFACE  *MRISfromParameterization(MRI_SP *mrisp, MRI_SURFACE *mris,
                                       int fno) ;
MRI_SURFACE  *MRISnormalizeFromParameterization(MRI_SP *mrisp,
                                                MRI_SURFACE *mris, int fno) ;
MRI_SP       *MRISgradientToParameterization(MRI_SURFACE *mris, MRI_SP *mrisp,
                                             float scale) ;
MRI_SURFACE  *MRISgradientFromParameterization(MRI_SP*mrisp,MRI_SURFACE *mris);

MRI_SP       *MRIScoordsToParameterization(MRI_SURFACE *mris, MRI_SP *mrisp,
                                           float scale, int which_vertices) ;
MRI_SURFACE  *MRIScoordsFromParameterization(MRI_SP *mrisp, MRI_SURFACE *mris, int which_vertices);


float         MRISPsample(MRI_SP *mrisp, float x, float y, float z, int fno) ;
MRI_SP       *MRISPblur(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma,
                        int fno) ;
MRI_SP       *MRISPconvolveGaussian(MRI_SP *mrisp_src, MRI_SP *mrisp_dst,
                                    float sigma, float radius, int fno) ;
MRI_SP       *MRISPalign(MRI_SP *mrisp_orig, MRI_SP *mrisp_src,
                         MRI_SP *mrisp_tmp, MRI_SP *mrisp_dst) ;
MRI_SP       *MRISPtranslate(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, int du,
                             int dv) ;
MRI_SP       *MRISPclone(MRI_SP *mrisp_src) ;
MRI_SP       *MRISPalloc(float scale, int nfuncs) ;
int          MRISPfree(MRI_SP **pmrisp) ;
MRI_SP       *MRISPread(char *fname) ;
int          MRISPwrite(MRI_SP *mrisp, char *fname) ;

int          MRISwriteArea(MRI_SURFACE *mris,const  char *sname) ;
int          MRISwriteMarked(MRI_SURFACE *mris,const  char *sname) ;
int          MRISmarkedToCurv(MRI_SURFACE *mris) ;
int          MRISdToCurv(MRI_SURFACE *mris) ;

double       MRISParea(MRI_SP *mrisp) ;


#define ORIGINAL_VERTICES   0
#define ORIG_VERTICES       ORIGINAL_VERTICES
#define GOOD_VERTICES       1
#define TMP_VERTICES        2
#define CANONICAL_VERTICES  3
#define CURRENT_VERTICES    4
#define INFLATED_VERTICES   5
#define FLATTENED_VERTICES  6
#define PIAL_VERTICES       7
#define TMP2_VERTICES       8
#define WHITE_VERTICES      9
#define TARGET_VERTICES     10
#define LAYERIV_VERTICES    11
#define LAYER4_VERTICES     LAYERIV_VERTICES
#define VERTEX_NORMALS      12
#define PIAL_NORMALS        13
#define WHITE_NORMALS       14

// Two ways of saving and restoring the VERTEX:xyz values
//
// The push/pop is preferable to using VERTEX members because it uses all the entries in a cache line
//
typedef struct MRISsavedXYZ MRISsavedXYZ;
MRISsavedXYZ* MRISsaveXYZ(MRIS *mris);
void          MRISloadXYZ(MRIS *mris, MRISsavedXYZ const *  pMRISsavedXYZ);  // does set
void          MRISfreeXYZ(MRIS *mris, MRISsavedXYZ      ** ppMRISsavedXYZ);  // does not set
void          MRISpopXYZ (MRIS *mris, MRISsavedXYZ      ** ppMRISsavedXYZ);  // set then free

int MRISsaveVertexPositions   (MRIS *mris, int which) ;
int MRISrestoreVertexPositions(MRIS *mris, int which) ;

int MRISrestoreNormals(MRIS *mris, int which) ;
int MRISsaveNormals   (MRIS *mris, int which) ;

/* constants for vertex->tethered */
#define TETHERED_NONE           0
#define TETHERED_FRONTAL_POLE   1
#define TETHERED_OCCIPITAL_POLE 2
#define TETHERED_TEMPORAL_POLE  3

#define MAX_TALAIRACH_Y         100.0f

/* constants for mris->hemisphere */
#define LEFT_HEMISPHERE         0
#define RIGHT_HEMISPHERE        1
#define NO_HEMISPHERE           2
#define BOTH_HEMISPHERES        3

#if 0
#define DEFAULT_A  44.0f
#define DEFAULT_B  122.0f
#define DEFAULT_C  70.0f
#else
#define DEFAULT_A  122.0f
#define DEFAULT_B  122.0f
#define DEFAULT_C  122.0f
#endif
#define DEFAULT_RADIUS  100.0f

#define MAX_DIM    DEFAULT_B

#if 1
#define DT_INCREASE  1.1 /* 1.03*/
#define DT_DECREASE  0.5
#else
#define DT_INCREASE  1.0 /* 1.03*/
#define DT_DECREASE  1.0
#endif
#define DT_MIN       0.01
#ifdef ERROR_RATIO
#undef ERROR_RATIO
#endif
#define ERROR_RATIO  1.03  /* 1.01 then 1.03 */

#define NO_LABEL     -1


#define WHITE_SURFACE 0
#define GRAY_SURFACE  1
#define GRAY_MID      2

#define REVERSE_X     0
#define REVERSE_Y     1
#define REVERSE_Z     2

int   MRISpositionSurface(MRI_SURFACE *mris, MRI *mri_brain,
                          MRI *mri_smooth, INTEGRATION_PARMS *parms);
int   MRISpositionSurfaces(MRI_SURFACE *mris, MRI **mri_flash,
                           int nvolumes, INTEGRATION_PARMS *parms);

int MRISpositionSurface_mef(MRI_SURFACE *mris,
                            MRI *mri_30,
                            MRI *mri_5,
                            INTEGRATION_PARMS *parms,
                            float weight30,
                            float weight5);

int   MRISmoveSurface(MRI_SURFACE *mris, MRI *mri_brain,
                      MRI *mri_smooth, INTEGRATION_PARMS *parms);
int   MRISscaleVals(MRI_SURFACE *mris, float scale) ;
int   MRISsetVals(MRI_SURFACE *mris, float val) ;
int   MRISsetValBaks(MRI_SURFACE *mris, float val) ;
int   MRISaverageD(MRI_SURFACE *mris, int navgs) ;
int   MRISgaussianFilterD(MRI_SURFACE *mris, double wt) ;
int   MRISaverageVals(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageVal2baks(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageVal2s(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageMarkedVals(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageMarkedStats(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageEveryOtherVertexPositions(MRI_SURFACE *mris,
                                           int navgs,
                                           int which) ;
int   MRISsoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs) ;
int   MRISsoapBubbleOrigVertexPositions(MRI_SURFACE *mris, int navgs) ;
int   MRISsoapBubbleTargetVertexPositions(MRI_SURFACE *mris, int navgs) ;
MRI   *MRISwriteSurfaceIntoVolume(MRI_SURFACE *mris, MRI *mri_template,
                                  MRI *mri) ;
#if 0
int   MRISmeasureCorticalThickness(MRI_SURFACE *mris, MRI *mri_brain,
                                   MRI *mri_wm, float nsigma) ;
#else
int   MRISmeasureCorticalThickness(MRI_SURFACE *mris, int nbhd_size,
                                   float max_thickness) ;
#endif

#include "mrishash.h"
int  MRISmeasureThicknessFromCorrespondence(MRI_SURFACE *mris, MHT *mht, float max_thick) ;
int MRISfindClosestOrigVertices(MRI_SURFACE *mris, int nbhd_size) ;
int MRISfindClosestPialVerticesCanonicalCoords(MRI_SURFACE *mris, int nbhd_size) ;

int   MRISmarkRandomVertices(MRI_SURFACE *mris, float prob_marked) ;
int   MRISmarkNegativeVertices(MRI_SURFACE *mris, int mark) ;
int   MRISripNegativeVertices(MRI_SURFACE *mris) ;
int   MRISclearGradient(MRI_SURFACE *mris) ;
int   MRISclearMarks(MRI_SURFACE *mris) ;
int   MRISclearMark2s(MRI_SURFACE *mris) ;
int   MRISclearFaceMarks(MRI_SURFACE *mris) ;
int   MRISclearFixedValFlags(MRI_SURFACE *mris) ;
int   MRIScopyFixedValFlagsToMarks(MRI_SURFACE *mris) ;
int   MRISclearAnnotations(MRI_SURFACE *mris) ;
int   MRISsetMarks(MRI_SURFACE *mris, int mark) ;
int   MRISsequentialAverageVertexPositions(MRI_SURFACE *mris, int navgs) ;
int   MRISreverseCoords(MRI_SURFACE *mris, int which_direction, int reverse_face_order, int which_coords) ;
int   MRISreverse(MRI_SURFACE *mris, int which, int reverse_face_order) ;
int   MRISdisturbOriginalDistances(MRI_SURFACE *mris, double max_pct) ;
double MRISstoreAnalyticDistances(MRI_SURFACE *mris, int which) ;
int   MRISnegateValues(MRI_SURFACE *mris) ;
int   MRIScopyMeansToValues(MRI_SURFACE *mris) ;
int   MRIScopyCurvatureToValues(MRI_SURFACE *mris) ;
int   MRIScopyCurvatureToImagValues(MRI_SURFACE *mris) ;
int   MRIScopyCurvatureFromValues(MRI_SURFACE *mris) ;
int   MRIScopyVal2ToVal(MRI_SURFACE *mris) ;
int   MRIScopyVal2BakToVal(MRI_SURFACE *mris) ;
int   MRIScopyCurvatureFromImagValues(MRI_SURFACE *mris) ;
int   MRIScopyImaginaryMeansToValues(MRI_SURFACE *mris) ;
int   MRIScopyStandardErrorsToValues(MRI_SURFACE *mris) ;
int   MRIScopyImagValuesToCurvature(MRI_SURFACE *mris) ;
int   MRIScopyValuesToCurvature(MRI_SURFACE *mris) ;
int   MRISaccumulateMeansInVolume(MRI_SURFACE *mris, MRI *mri, int mris_dof,
                                  int mri_dof, int coordinate_system,int sno);
int   MRISaccumulateStandardErrorsInVolume(MRI_SURFACE *mris, MRI *mri,
    int mris_dof, int mri_dof,
    int coordinate_system, int sno) ;
int   MRISaccumulateMeansOnSurface(MRI_SURFACE *mris, int total_dof,
                                   int new_dof);
int   MRISaccumulateImaginaryMeansOnSurface(MRI_SURFACE *mris, int total_dof,
    int new_dof);
int   MRISaccumulateStandardErrorsOnSurface(MRI_SURFACE *mris,
    int total_dof,int new_dof);
#define GRAY_WHITE     1
#define GRAY_CSF       2
int
MRIScomputeMaxGradBorderValuesPial(MRI_SURFACE *mris,MRI *mri_brain,
                                   MRI *mri_smooth, double sigma,
                                   float max_thickness,
                                   float dir, FILE *log_fp,
                                   int callno, MRI *mri_mask) ;
int
MRIScomputeMaxGradBorderValues(MRI_SURFACE *mris,MRI *mri_brain,
                               MRI *mri_smooth, double sigma,
                               float max_thickness, float dir, FILE *log_fp,
                               MRI *mri_wm, int callno) ;
int MRIScomputeInvertedGrayWhiteBorderValues(MRI_SURFACE *mris,
                                             MRI *mri_brain,
                                             MRI *mri_smooth,
                                             double inside_hi,
                                             double border_hi,
                                             double border_low,
                                             double outside_low,
                                             double outside_hi,
                                             double sigma,
                                             float max_thickness,
                                             FILE *log_fp);
int MRIScomputeInvertedPialBorderValues(MRI_SURFACE *mris,
                                        MRI *mri_brain,
                                        MRI *mri_smooth,
                                        double inside_hi,
                                        double border_hi,
                                        double border_low,
                                        double outside_low,
                                        double outside_hi,
                                        double sigma,
                                        float max_thickness,
                                        FILE *log_fp);
int   MRIScomputeBorderValues(MRI_SURFACE *mris,
                              MRI *mri_brain,
                              MRI *mri_smooth,
                              double inside_hi,
                              double border_hi,
                              double border_low,
                              double outside_low,
                              double outside_hi,
                              double sigma,
                              float max_dist,
                              FILE *log_fp,
                              int white,
                              MRI *mri_mask, double thresh, int flags, MRI *mri_aseg);
int  MRIScomputeWhiteSurfaceValues(MRI_SURFACE *mris, MRI *mri_brain,
                                   MRI *mri_smooth);
int  MRIScomputeGraySurfaceValues(MRI_SURFACE *mris, MRI *mri_brain,
                                  MRI *mri_smooth, float gray_surface);
int  MRIScomputeDistanceErrors(MRI_SURFACE *mris, int nbhd_size,int max_nbrs);
int  MRISsetAllMarks(MRI_SURFACE *mris, int mark) ;
int  MRISscaleCurvature(MRI_SURFACE *mris, float scale) ;
int  MRISwriteTriangularSurface(MRI_SURFACE *mris,const char *fname) ;

int  MRISbuildFileName(MRI_SURFACE *mris, const char *sname, char *fname) ;
int  MRISsmoothSurfaceNormals(MRI_SURFACE *mris, int niter) ;
int  MRISsoapBubbleD(MRI_SURFACE *mris, int niter) ;
int  MRISsoapBubbleVals(MRI_SURFACE *mris, int niter) ;
int  MRISmodeFilterVals(MRI_SURFACE *mris, int niter) ;
int  MRISmodeFilterAnnotations(MRI_SURFACE *mris, int niter) ;
int  MRISmodeFilterZeroVals(MRI_SURFACE *mris) ;
int  MRISreadBinaryAreas(MRI_SURFACE *mris,const  char *mris_fname) ;
int  MRISwriteAreaErrorToValFile(MRI_SURFACE *mris,const  char *name) ;
int  MRIStransform(MRI_SURFACE *mris, MRI *mri,
                   TRANSFORM *transform, MRI *mri_dst) ;
int  MRISmatrixMultiply(MRIS *mris, MATRIX *M);
int MRISltaMultiply(MRIS *surf, const LTA *lta);

int  MRISanisotropicScale(MRI_SURFACE *mris, float sx, float sy, float sz) ;
double MRIScomputeVertexSpacingStats(MRI_SURFACE *mris, double *psigma,
                                     double *pmin, double *pmax, int *pvno,
                                     int *pvno2, int which_vertices);
double MRIScomputeTotalVertexSpacingStats(MRI_SURFACE *mris, double *psigma,
                                          double *pmin, double *pmax,
                                          int *pvno,
                                          int *pvno2);
double MRIScomputeFaceAreaStats(MRI_SURFACE *mris, double *psigma,
                                double *pmin, double *pmax);
int MRISprintTessellationStats(MRI_SURFACE *mris, FILE *fp) ;
int MRISprintVertexStats(MRI_SURFACE *mris, int vno, FILE *fp, int which_vertices) ;
int MRISmergeIcosahedrons(MRI_SURFACE *mri_src, MRI_SURFACE *mri_dst) ;
int MRISinverseSphericalMap(MRI_SURFACE *mris, MRI_SURFACE *mris_ico) ;

////////////////////////////////////////////////////
/* for mris_topo_fixer */
typedef struct
{
  MRIS *mris;
  int n_vertices,n_faces; //number of vertices & faces before top. correction
  int *vtrans_to,*ftrans_to,*vtrans_from,*ftrans_from;
  MRIS *mris_source;
}
MRI_PATCH, MRIP;

/* Different verbose mode for mris_fix_topology and mris_topo_fix */
#define VERBOSE_MODE_DEFAULT 1
#define VERBOSE_MODE_LOW 2
#define VERBOSE_MODE_MEDIUM 3
#define VERBOSE_MODE_HIGH 4

#include "histo.h" // HISTOGRAM

typedef struct
{
  int verbose; // verbose mode
  int smooth;  // smoothing defect
  int match;   // using local intensity to match surface
  int defect_number; // the defect_number
  int mode; // which mode to use (not used so far)
  int no_self_intersections; //to prevent self-intersection
  int contrast; //direction of the contrast
  int write; //writing out temporary surfaces using fname
  char fname[STRLEN] ;
  double nattempts_percent;
  int nattempts;
  int minimal_mode;
  int nminattempts;
  double minimal_loop_percent;
  int nminimal_attempts;
  int max_face;
  void* patchdisk; // the patching surfaces
  // loglikelihood
  double  l_mri ;
  double  l_curv ;
  double  l_qcurv;
  double  l_unmri ;
  int volume_resolution; /* used if l_unmri is on */
  MRI *mri; //brain volume
  MRI *mri_wm; //wm volume
  HISTOGRAM *h_k1, *h_k2,*h_gray,*h_white,*h_dot,*h_border, *h_grad;
  MRI *mri_gray_white, *mri_k1_k2;
  MATRIX *transformation_matrix;
  //defect info
  MRIP *mrip;
  MRIS *mris_defect;
  void *defect_list;
  void   *dp;
  //statistics
  float fitness,initial_fitness;
  int ngeneratedpatches,nselfintersectingpatches;
}
TOPOFIX_PARMS;

void MRIScomputeInitialFitness(MRIS *mris, TOPOFIX_PARMS *parms);
void MRIScomputeDistanceVolume(TOPOFIX_PARMS *parms,
                               float distance_to_surface);
void MRISsaveLocal(MRIS *mris, TOPOFIX_PARMS *parms, char *name); //floflo
int MRISisSurfaceValid(MRIS *mris, int patch,int verbose);
void MRISinitTopoFixParameters(MRIS *mris, TOPOFIX_PARMS *parms);
void MRISinitDefectParameters(MRIS *mris, TOPOFIX_PARMS *parms);
void TOPOFIXfreeDP(TOPOFIX_PARMS *parms);
void MRISinitDefectPatch(MRIS *mris, TOPOFIX_PARMS *parms);
void MRISdefectMatch(MRIS *mris, TOPOFIX_PARMS *parms);
void MRISprintInfo(TOPOFIX_PARMS *parms);
double MRIScomputeFitness(MRIS* mris,TOPOFIX_PARMS *parms,int verbose);
int IsMRISselfIntersecting(MRI_SURFACE *mris);
void MRISmapOntoSphere(MRIS *mris);
void MRISidentifyDefects(MRIS *mris, TOPOFIX_PARMS *parms);
void MRISmarkPatchVertices(MRIS *mris,
                           TOPOFIX_PARMS *parms,
                           int ninitvertices);
void MRISmarkBorderVertices(MRIS *mris,
                            TOPOFIX_PARMS *parms,
                            int mark);

#define GREEDY_SEARCH 0
#define GENETIC_SEARCH 1
#define RANDOM_SEARCH 2

typedef struct
{
  int     max_patches ;
  int     max_unchanged ;
  int     niters;  /* stop the genetic algorithm after n iterations */
  int     genetic; /* to use the genetic algorithm */
  int     search_mode; /* which topology correction (greedy,genetic,random)*/
  int     retessellation_mode; /* which retessellation (ordering dependent)*/
  int     vertex_eliminate; /* kill less used vertices */
  int     initial_selection; /* find set of initial chromosomes */
  double  l_mri ;
  double  l_curv ;
  double  l_qcurv;
  double  l_unmri ;
  int volume_resolution; /* used if l_unmri is on */
  int keep; /* keep every vertex in the defect */
  int edge_table; /*using edge table or not */
  char *save_fname; /* save results into folder name save_fname */
  int defect_number; /* save only for one specific defect */
  int verbose; /* outputs information */
  int smooth; /* smoothing the patch */
  int match; /* matching the patch onto the surface using local parameters */
  int movie; /* save interesting files for movie purpose */
  int correct_defect; /* correct only one single defect */
  int check_surface_intersection; /* check if self-intersection happens */
  int optimal_mapping; /* find optmal mapping by genrating sevral mappings */
}
TOPOLOGY_PARMS ;

int MRISmarkOrientationChanges(MRI_SURFACE *mris);
MRIS* MRISextractMainComponent(MRI_SURFACE *mris,
                               int do_not_extract,
                               int verbose,
                               int *ncpts);
MRIS* MRISextractMarkedVertices(MRIS *mris);
MRIS* MRISremoveRippedSurfaceElements(MRIS *mris);

MRI_SURFACE *MRIScorrectTopology(MRI_SURFACE *mris, MRI *mri, 
   MRI *mri_wm, int nsmooth, TOPOLOGY_PARMS *parms, char *defectbasename);

int MRISsmoothOnSphere(MRIS* mris, int niters);
int mrisCountIntersectingFaces(MRIS *mris, int*flist , int nfaces);
int MRIScountNegativeFaces(MRI_SURFACE *mris) ;
int MRISevertSurface(MRI_SURFACE *mris) ;
int MRISripDefectiveFaces(MRI_SURFACE *mris) ;
int MRISunrip(MRI_SURFACE *mris) ;
int MRISdivideLongEdges(MRI_SURFACE *mris, double thresh) ;
int MRISdivideEdges(MRI_SURFACE *mris, int npoints) ;
int MRISremoveTriangleLinks(MRI_SURFACE *mris) ;
int MRISsetOriginalFileName(char *orig_name) ;
int MRISsetSulcFileName(const char *sulc_name) ;
int MRISsetCurvatureName(int nth, const char *name);
int MRISprintCurvatureNames(FILE *fp);
int MRISsetInflatedFileName(char *inflated_name) ;
int MRISsetRegistrationSigmas(float *sigmas, int nsigmas) ;

int MRISextractVertexCoords(MRI_SURFACE *mris, float *locations[3], int which_vertices) ;
int MRISimporttVertexCoords(MRI_SURFACE *mris, float *locations[3], int which_vertices) ;
int MRISextractCurvatureVector(MRI_SURFACE *mris, float *curvs) ;
int MRISextractCurvatureDoubleVector(MRI_SURFACE *mris, double *curvs) ;
int MRISextractCurvatureVectorDouble(MRI_SURFACE *mris, double *curvs, int offset) ;
#define MRISexportCurvatureVector  MRISextractCurvatureVector

int MRISimportCurvatureVector(MRI_SURFACE *mris, float *curvs) ;
int MRISimportValVector(MRI_SURFACE *mris, float *vals) ;
int MRISexportValVector(MRI_SURFACE *mris, float *vals) ;
int MRISexportValVectorDouble(MRI_SURFACE *mris, double *vals, int offset) ;
int MRISimportValFromMatrixColumn(MRI_SURFACE *mris, MATRIX *m, int col) ;
int MRISimportValFromMRI(MRI_SURFACE *mris, MRI *mri, int frame) ;
MRI *MRISreadParameterizationToSurface(MRI_SURFACE *mris, char *fname) ;


/* multi-timepoint (or stc) files */
int  MRISwriteStc(char *fname, MATRIX *m_data, float epoch_begin_lat,
                  float sample_period, int *vertices) ;

#define MAX_LINKS 10
typedef struct
{
  double     tx, ty, tz ;   // target coords
  int        l_in ;         // inside label
  int        l_out ;        // outside label
  HISTOGRAM *h_in ;         // inside label histogram
  HISTOGRAM *h_out ;        // inside label histogram
  double     mag ;          // directional derivative in normal dir initially
  double     p ;            // p-value for user to fill in
}
VERTEX_INFO ;

typedef struct
{
  MRI_SURFACE *mris[MAX_SURFACES] ;
  int         nsurfaces ;
  int         nvertices ;   // sum of individual surface nvertices
  int         nfaces ;      // sum of individual surface faces
  int         labels[MAX_SURFACES] ;  // interior label of the structure
  MRI_SURFACE *mris_total ; // single surface with all other surfs in it
  VERTEX_INFO *vi ;                   // one/vertex in mris_total
  int         vstart[MAX_SURFACES] ;  /* starting index of this
                                         surface in mris_total vertices */
}
MRI_SURFACE_ARRAY, MSA ;
#if 1
float  MRISdistanceToSurface(MRI_SURFACE *mris, MHT *mht,
                             float x0, float y0, float z0,
                             float nx, float ny, float nz) ;
int    MRISexpandSurface(MRI_SURFACE *mris,
                         float distance,
                         INTEGRATION_PARMS *parms, int use_thickness, int nsurfaces) ;
int MRISripZeroThicknessRegions(MRI_SURFACE *mris) ;

#endif

/* cortical ribbon */
MRI   *MRISribbon(MRI_SURFACE *inner_mris,
                  MRI_SURFACE *outer_mris,
                  MRI *mri_src,
                  MRI *mri_dst);
MRI   *MRISaccentuate(MRI *mri_src,
                      MRI *mri_dst,
                      int lo_thresh,
                      int hi_thresh);
MRI *MRISfillInterior(MRI_SURFACE *mris,
                      double resolution,
                      MRI *mri_interior) ;
int MRISfillInteriorRibbonTest(char *subject, int UseNew, FILE *fp);
MRI   *MRISshell(MRI *mri_src,
                 MRI_SURFACE *mris,
                 MRI *mri_dst,
                 int clearflag);
MRI   *MRISfloodoutside(MRI *mri_src,MRI *mri_dst);
// src never used. dst is modified

/* high resolution cortical ribbon */
MRI   *MRISpartialribbon(MRI_SURFACE *inner_mris_lh,
                         MRI_SURFACE *outer_mris_lh,
                         MRI_SURFACE *inner_mris_rh,
                         MRI_SURFACE *outer_mris_rh,
                         MRI *mri_src,
                         MRI *mri_dst,
                         MRI *mri_mask);
MRI   *MRISpartialaccentuate(MRI *mri_src,
                             MRI *mri_dst,
                             int lo_thresh,
                             int hi_thresh);
MRI   *MRISpartialshell(MRI *mri_src,
                        MRI_SURFACE *mris,
                        MRI *mri_dst,
                        int clearflag);
MRI   *MRISpartialfloodoutside(MRI *mri_src,
                               MRI *mri_dst);
// src is used. dst must be initialized

#define MRIS_BINARY_QUADRANGLE_FILE    0    /* homegrown */
#define MRIS_ASCII_TRIANGLE_FILE       1    /* homegrown */
#define MRIS_ASCII_FILE MRIS_ASCII_TRIANGLE_FILE
#define MRIS_GEO_TRIANGLE_FILE         2    /* movie.byu format */
#define MRIS_ICO_SURFACE               3
#define MRIS_TRIANGULAR_SURFACE        MRIS_ICO_SURFACE
#define MRIS_ICO_FILE                  4
#define MRIS_VTK_FILE                  5
#define MRIS_STL_FILE                  6
#define MRIS_GIFTI_FILE                7
#define MRIS_ANNOT_FILE                8
#define MRIS_VOLUME_FILE               9


#define IS_QUADRANGULAR(mris)                   \
  (mris->type == MRIS_BINARY_QUADRANGLE_FILE)


/* actual constants are in mri.h */
#define RH_LABEL           MRI_RIGHT_HEMISPHERE
#define RH_LABEL2          MRI_RIGHT_HEMISPHERE2
#define LH_LABEL           MRI_LEFT_HEMISPHERE

#define MRISPvox(m,u,v)   (*IMAGEFpix(m->Ip,u,v))


/* VECTORIAL_REGISTRATION*/
int MRISPfunctionVectorVals(MRI_SURFACE_PARAMETERIZATION *mrisp,
                            MRI_SURFACE *mris,
                            float x, float y, float z,
                            int *frames, int nframes, double *vals);
MRI_SURFACE * MRISfromParameterizations(MRI_SP *mrisp,
                                        MRI_SURFACE *mris, int *frames,
                                        int *indices,int nframes);
MRI_SP * MRIStoParameterizations(MRI_SURFACE *mris,
                                 MRI_SP *mrisp, float scale,int *frames,
                                 int *indices,int nframes);
MRI_SP * MRISPblurFrames(MRI_SP *mrisp_src,
                         MRI_SP *mrisp_dst, float sigma, int *frames,
                         int nframes);

typedef struct
{
  float   x ;
  float   y ;
  float   z ;
}
SMALL_VERTEX ;

typedef struct
{
  int   nvertices ;
  SMALL_VERTEX *vertices ;
}
SMALL_SURFACE ;


int   MRISpositionSurfaceArray(MRI_SURFACE_ARRAY *msa, MRI *mri_brain,
                               MRI *mri_smooth, INTEGRATION_PARMS *parms);

#define MRISSread  MRISreadVerticesOnly
SMALL_SURFACE *MRISreadVerticesOnly(char *fname) ;
int           MRISSfree(SMALL_SURFACE **pmriss) ;


int  MRISsubsampleDist(MRI_SURFACE *mris, float spacing) ;
int  MRISwriteDecimation(MRI_SURFACE *mris, char *fname) ;
int  MRISreadDecimation(MRI_SURFACE *mris, char *fname) ;


#define VERTEX_COORDS      0
#define VERTEX_VALS        1
#define VERTEX_VAL         VERTEX_VALS
#define VERTEX_AREA        2
#define VERTEX_CURV        3
#define VERTEX_CURVATURE   VERTEX_CURV
#define VERTEX_LABEL       4
#define VERTEX_ANNOTATION  5
#define VERTEX_DX          6
#define VERTEX_DY          7
#define VERTEX_DZ          8
#define VERTEX_STATS       9
#define VERTEX_LOGODDS     10
#define VERTEX_MARKS       11


int MRISclearOrigArea(MRI_SURFACE *mris) ;
int MRISclearOrigDistances(MRI_SURFACE *mris) ;
int MRIScombine(MRI_SURFACE *mris_src, MRI_SURFACE *mris_total,
                MRIS_HASH_TABLE *mht, int which) ;
int MRISsphericalCopy(MRI_SURFACE *mris_src, MRI_SURFACE *mris_total,
                      MRIS_HASH_TABLE *mht, int which) ;
int   MRISorigAreaToCurv(MRI_SURFACE *mris) ;
int   MRISareaToCurv(MRI_SURFACE *mris) ;
int   MRISclear(MRI_SURFACE *mris, int which) ;
int   MRISnormalize(MRI_SURFACE *mris, int dof, int which) ;

int  MRIScopyMRI(MRIS *Surf, MRI *Src, int Frame, char *Field);
MRI *MRIcopyMRIS(MRI *mri, MRIS *surf, int Frame, char *Field);

MRI *MRISsmoothMRI(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *IncMask, MRI *Targ);
MRI *MRISsmoothMRIFast(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *IncMask,  MRI *Targ);
MRI *MRISsmoothMRIFastD(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *IncMask,  MRI *Targ);
int MRISsmoothMRIFastCheck(int nSmoothSteps);
int MRISsmoothMRIFastFrame(MRIS *Surf, MRI *Src, int frame, int nSmoothSteps, MRI *IncMask);


int  MRISclearFlags(MRI_SURFACE *mris, int flags) ;
int  MRISsetCurvature(MRI_SURFACE *mris, float val) ;
int  MRISsetFlags(MRI_SURFACE *mris, int flags) ;

int MRISgaussianFilterD(MRI_SURFACE *mris, double wt) ;
int MRISmedianFilterD(MRI_SURFACE *mris, int nmedians, int vtotal) ;
int MRISmedianFilterVals(MRI_SURFACE *mris, int nmedians) ;
int MRISmedianFilterVal2s(MRI_SURFACE *mris, int nmedians) ;
int MRISmedianFilterVal2baks(MRI_SURFACE *mris, int nmedians) ;
int MRISmedianFilterCurvature(MRI_SURFACE *mris, int nmedians);
int MRISfileNameType(const char *fname) ;
unsigned long MRISeraseOutsideOfSurface(float h,
                                        MRI* mri_dst,
                                        MRIS *mris,
                                        unsigned char val) ;

/* Some utility functions to handle reading and writing annotation
   values. MRISRGBToAnnot stuffs an r,g,b tuple into an annotation
   value and MRISAnnotToRGB separates an annotation value into an
   r,g,b tuple. */
#define MRISAnnotToRGB(annot,r,g,b)             \
  r = annot & 0xff ;                            \
  g = (annot >> 8) & 0xff ;                     \
  b = (annot >> 16) & 0xff ;
#define MRISRGBToAnnot(r,g,b,annot)                                     \
  annot = ((r) & 0xff) | (((g) & 0xff) << 8) | (((b) & 0xff) << 16);

int MRISextendedNeighbors(MRIS *SphSurf,int TargVtxNo, int CurVtxNo,
                          double DotProdThresh, int *XNbrVtxNo,
                          double *XNbrDotProd, int *nXNbrs,
                          int nXNbrsMax, int DistType);
MRI *MRISgaussianSmooth(MRIS *Surf, MRI *Src, double GStd, MRI *Targ,
                        double TruncFactor);
MRI *MRISdistSphere(MRIS *surf, double dmax);
int MRISgaussianWeights(MRIS *surf, MRI *dist, double GStd);
MRI *MRISspatialFilter(MRI *Src, MRI *wdist, MRI *Targ);

MATRIX *surfaceRASToSurfaceRAS_(MRI *src, MRI *dst, LTA *lta);

int MRISsurf2surf(MRIS *mris, MRI *dst, LTA *lta);
// convert all vertex positions
int MRISsurf2surfAll(MRIS *mris, MRI *dst, LTA *lta);
void MRISsetReadFrame(int frame);
int MRISaddCommandLine(MRI_SURFACE *mris, char *cmdline) ;
int MRISgetReadFrame(void);
int MRISabsCurvature(MRI_SURFACE *mris) ;
int MRISabsVals(MRI_SURFACE *mris) ;
int MRISsmoothFrames(MRI_SURFACE *mris, MRI *mri, int navgs) ;
int MRISwriteFrameToValues(MRI_SURFACE *mris, MRI *mri, int frame) ;
int MRISreadFrameFromValues(MRI_SURFACE *mris, MRI *mri, int frame) ;
MRI *MRISar1(MRIS *surf, MRI *src, MRI *mask, MRI *ar1);
int **MRIScrsLUT(MRIS *surf, MRI *src);
int MRIScrsLUTFree(int **crslut);
int MRISremoveOverlapWithSmoothing(MRI_SURFACE *mris,
                                   INTEGRATION_PARMS *parms) ;
int MRISupsampleIco(MRI_SURFACE *mris, MRI_SURFACE *mris_new) ;
int MRIScopyVolGeomFromMRI(MRI_SURFACE *mris, MRI *mri) ;

// the next function assumes the surface vertex positions are 
//    set in the voxel coordinates of srcMri
// and it constructs a valid RAS for the surface
void MRISsetVolumeForSurface(MRI_SURFACE *mris, MRI *srcMri);

MRI *MRISremoveRippedFromMask(MRIS *surf, MRI *mask, MRI *outmask);
int MRISremoveIntersections(MRI_SURFACE *mris) ;
int MRIScopyMarkedToMarked2(MRI_SURFACE *mris) ;
int MRIScopyMarked2ToMarked(MRI_SURFACE *mris) ;
int MRIScopyMarkedToMarked3(MRI_SURFACE *mris) ;
int MRIScopyMarked3ToMarked(MRI_SURFACE *mris) ;
int MRISexpandMarked(MRI_SURFACE *mris) ;
int MRISnotMarked(MRI_SURFACE *mri) ;
double MRISsmoothingArea(MRIS *mris, int vtxno, int niters);

int MRISopenMarked(MRI_SURFACE *mris, int order) ;
int MRIScloseMarked(MRI_SURFACE *mris, int order) ;
int MRISerodeMarked(MRI_SURFACE *mris, int ndil) ;
int MRISdilateMarked(MRI_SURFACE *mris, int ndil) ;
int MRISerodeRipped(MRI_SURFACE *mris, int ndil) ;
int MRISdilateRipped(MRI_SURFACE *mris, int ndil) ;
MRI *MRISdilateConfined(MRIS *surf,
                        MRI *mask,
                        int annotidmask,
                        int niters,
                        int newid);
MRI *MRISfbirnMask_SFG_Cing(MRIS *surf);
MRI *MRISfbirnMask_MOF_RACing(MRIS *surf);

int   MRISvalidVertices(MRI_SURFACE *mris) ;
int MRISmarkedVertices(MRI_SURFACE *mris) ;
int MRISmarkVerticesWithValOverThresh(MRI_SURFACE *mris, float thresh) ;

int MRIScomputeClassStatistics(MRI_SURFACE *mris,
                               MRI *mri,
                               float *pwhite_mean,
                               float *pwhite_std,
                               float *pgray_mean,
                               float *pgray_std) ;
int MRIScomputeClassModes(MRI_SURFACE *mris,
                          MRI *mri,
                          float *pwhite_mode,
                          float *pgray_mode,
                          float *pcsf_mode,
			  float *pwhite_std,
			  float *pgray_std,
			  float *pcsf_std);
int MRISrasToVoxel(MRI_SURFACE *mris,
                   MRI *mri,
                   double xs, double ys, double zs,
                   double *pxv, double *pyv, double *pzv) ;
int MRISrestoreRipFlags(MRI_SURFACE *mris) ;
int MRISstoreRipFlags(MRI_SURFACE *mris) ;
int MRISripMedialWall(MRI_SURFACE *mris) ;
int MRISripMarked(MRI_SURFACE *mris) ;
int MRISripUnmarked(MRI_SURFACE *mris) ;
int MRISzeroMedialWallCurvature(MRI_SURFACE *mris) ;
int MRISvertexNormalInVoxelCoords(MRI_SURFACE *mris,
                                  MRI *mri,
                                  int vno,
                                  double *pnx, double *pny, double *pnz) ;

#define MRISgetCoords(v,c,vx,vy,vz) \
 switch(c) { \
   case ORIGINAL_VERTICES:  (*vx) = (v)->origx;  (*vy) = (v)->origy;  (*vz) = (v)->origz; break; \
   case TMP_VERTICES:       (*vx) = (v)->tx;     (*vy) = (v)->ty;     (*vz) = (v)->tz; break; \
   case CANONICAL_VERTICES: (*vx) = (v)->cx;     (*vy) = (v)->cy;     (*vz) = (v)->cz; break; \
   case CURRENT_VERTICES:   (*vx) = (v)->x;      (*vy) = (v)->y;      (*vz) = (v)->z; break; \
   case TARGET_VERTICES:   (*vx) = (v)->targx;   (*vy) = (v)->targy;  (*vz) = (v)->targz; break; \
   case INFLATED_VERTICES:  (*vx) = (v)->infx;   (*vy) = (v)->infy;   (*vz) = (v)->infz; break; \
   case FLATTENED_VERTICES: (*vx) = (v)->fx;     (*vy) = (v)->fy;     (*vz) = (v)->fz; break; \
   case PIAL_VERTICES:      (*vx) = (v)->pialx;  (*vy) = (v)->pialy;  (*vz) = (v)->pialz; break; \
   case TMP2_VERTICES:      (*vx) = (v)->tx2;    (*vy) = (v)->ty2;    (*vz) = (v)->tz2; break; \
   case WHITE_VERTICES:     (*vx) = (v)->whitex; (*vy) = (v)->whitey; (*vz) = (v)->whitez; break; \
   default: break; \
 }

#include "label.h" // LABEL
int MRISlogOdds(MRI_SURFACE *mris, LABEL *area, double slope)  ;
MRI_SP  *MRISPorLabel(MRI_SP *mrisp, MRI_SURFACE *mris, LABEL *area) ;
MRI_SP  *MRISPandLabel(MRI_SP *mrisp, MRI_SURFACE *mris, LABEL *area) ;
MRI *MRISlabel2Mask(MRIS *surf, LABEL *lb, MRI *mask);
int   MRIScomputeAverageCircularPhaseGradient(MRI_SURFACE *mris,
    LABEL *area,
    float *pdx,
    float *pdy,
    float *pdz);

int MRISmaskLabel(MRI_SURFACE *mris, LABEL *area) ;
int MRISmaskNotLabel(MRI_SURFACE *mris, LABEL *area) ;
int MRISripLabel(MRI_SURFACE *mris, LABEL *area) ;
int MRISripNotLabel(MRI_SURFACE *mris, LABEL *area) ;
int MRISsegmentMarked(MRI_SURFACE *mris,
                      LABEL ***plabel_array,
                      int *pnlabels,
                      float min_label_area) ;
int MRISsegmentAnnotated(MRI_SURFACE *mris,
                         LABEL ***plabel_array,
                         int *pnlabels,
                         float min_label_area) ;

int MRISaverageGradients(MRI_SURFACE *mris, int num_avgs) ;
int MRISaverageGradientsFast(MRI_SURFACE *mris, int num_avgs);
int MRISaverageGradientsFastCheck(int num_avgs);

int MRISnormalTermWithGaussianCurvature(MRI_SURFACE *mris,double l_lambda) ;
int MRISnormalSpringTermWithGaussianCurvature(MRI_SURFACE *mris,
                                              double gaussian_norm,
                                              double l_spring) ;
int MRISmakeDensityMap(MRI_SURFACE *mris,
                       double resolution,
                       double radius,
                       int diag_no,
                       MRI **pmri);
double MRIScomputeWhiteVolume(MRI_SURFACE *mris,
                              MRI *mri_aseg,
                              double resolution);
int MRIShistoThresholdCurvature(MRI_SURFACE *mris, float thresh_pct);
int MRISthresholdCurvature(MRI_SURFACE *mris, float thresh, int use_abs);
int MRISbinarizeCurvature(MRI_SURFACE *mris,
                          float thresh,
                          float low,
                          float high,
                          int use_abs);
int MRISsetVal2(MRI_SURFACE *mris, float val);
MRI *MRIScomputeDistanceToSurface(MRI_SURFACE *mris,
                                  MRI *mri,
                                  float resolution) ;
int MRISdistanceTransform(MRI_SURFACE *mris,LABEL *area, int mode) ;
int MRISinvertMarks(MRI_SURFACE *mris) ;


/* value field labels for swap_vertex_fields(), also used in tcl script */
#define FIELD_CURV      0
#define FIELD_CURVBAK   1
#define FIELD_VAL       2
#define FIELD_IMAG_VAL  3
#define FIELD_VAL2      4
#define FIELD_VALBAK    5
#define FIELD_VAL2BAK   6
#define FIELD_STAT      7

#define SCLV_VAL          0
#define SCLV_VAL2         1
#define SCLV_VALBAK       2
#define SCLV_VAL2BAK      3
#define SCLV_VALSTAT      4
#define SCLV_IMAG_VAL     5
#define SCLV_MEAN         6
#define SCLV_MEAN_IMAG    7
#define SCLV_STD_ERROR    8
#define NUM_SCALAR_VALUES 9
HISTOGRAM *MRISgetHistogram(MRI_SURFACE *mris, int nbins, int field);

// Discrete Principal Curvature and Related  vvvvvvvvvvvvvvvvvv
// Column / width for formatted output
#define		G_LC				50
#define		G_RC				30

int		slprints(
    char*		apch_txt
);

void	cprintf(
	char*		apch_left,
	float		af_right
);

void	cprints(
	char*		apch_left,
	char*		apch_right
);

void	cprintd(
	char*		apch_left,
	int		a_right
);

short	FACES_aroundVertex_reorder(
	MRIS*			apmris,
    	int			avertex,
    	VECTOR*			pv_geometricOrder
);

short	FACE_vertexIndex_find(
    	FACE*			pFace,
    	int 			avertex 
);

short	VECTOR_elementIndex_findNotEqual(
	VECTOR*			apV,
	float			af_searchTerm
);

short	VECTOR_elementIndex_find(
	VECTOR*			apV,
	float			af_searchTerm
);

short	MRIS_vertexProgress_print(
    	MRIS*			apmris,
    	int			avertex,
    	char*			apch_message
);

int	FACE_vertexIndexAtMask_find(
	FACE*			apFACE_I,
	VECTOR*			apv_verticesCommon
);

short	VERTICES_commonInFaces_find(
	FACE*			apFACE_I,
	FACE*			apFACE_J,
	VECTOR*			apv_verticesCommon
);

short	FACES_Hcurvature_determineSign(
    	MRIS*			apmris,
    	int			apFACE_O_fno,
    	int			apFACE_I_fno
);

int	VERTEX_faceAngles_determine(
    	MRIS*			apmris,
    	int			avertex,
    	VECTOR*			apv_angle
);

int	VERTEX_faceMinMaxAngles_determine(
    	MRIS*			apmris,
    	int			avertex,
    	int*			ap_minIndex,
    	float*			apf_minAngle,
    	int*			ap_maxIndex,
    	float*			apf_maxAngle
);

int	MRIS_facesAtVertices_reorder(
    	MRIS*			apmris
);

int	MRIScomputeGeometricProperties(
    	MRIS*			apmris
);

short	FACES_aroundVertex_reorder(
    	MRIS*			apmris,
    	int			avertex,
    	VECTOR*			pv_geometricOrder
);

float	FACES_angleNormal_find(
    	MRIS*			apmris,
    	int			apFACE_I_fno,
    	int			apFACE_J_fno
);

float	FACES_commonEdgeLength_find(
    	MRIS*			apmris,
    	FACE*			apFACE_I,
    	FACE*			apFACE_J
);

short	MRIS_discreteKH_compute(
	MRIS*			apmris
);

short	MRIS_discretek1k2_compute(
	MRIS*			apmris,
	short			ab_signedPrincipals
);

// The discrete curvature calculations are based on the Gauss-Bonnet Scheme.
// 
// For 'K' and 'H', see
// 
// @article{1280456,
//  author = {Evgeni Magid and Octavian Soldea and Ehud Rivlin},
//  title = {A comparison of Gaussian and mean curvature estimation 
// 	  methods on triangular meshes of range image data},
//  journal = {Comput. Vis. Image Underst.},
//  volume = {107},
//  number = {3},
//  year = {2007},
//  issn = {1077-3142},
//  pages = {139--159},
//  doi = {http://dx.doi.org/10.1016/j.cviu.2006.09.007},
//  publisher = {Elsevier Science Inc.},
//  address = {New York, NY, USA},
//  }
// 
// and for k1 and k2 see:
// 
// @misc{ meyer02discrete,
//   author = "M. MEYER and M. DESBRUN and P. SCHR and A. BARR",
//   title = "Discrete DifferentialGeometry Operators for Triangulated 2-Manifolds",
//   text = "MEYER, M., DESBRUN, M., SCHR ODER, P., AND BARR, 
// 	  A. H. Discrete DifferentialGeometry
//     	  Operators for Triangulated 2-Manifolds, 2002. VisMath.",
//   year = "2002",
//   url = "citeseer.ist.psu.edu/meyer02discrete.html" }
// 

short	MRIScomputeSecondFundamentalFormDiscrete(
	MRIS*			apmris,
	short			ab_signedPrincipals
);

int  	MRISminMaxCurvatureIndicesLookup(
  	MRI_SURFACE*  		apmris,
  	int*   			ap_vertexMin,
  	int*   			ap_vertexMax
);

int  	MRISvertexCurvature_set(
  	MRI_SURFACE*  		apmris,
  	int   			aindex,
  	float   		af_val
);

int	MRISzeroCurvature(
  	MRI_SURFACE*  		apmris
);

int	MRISuseK1Curvature(
  	MRI_SURFACE*  		mris
);

int	MRISuseK2Curvature(
  	MRI_SURFACE*  		mris
);

// end Discrete Principal Curvature and Related ^^^^^^^^^^^^^^^^^^

void UpdateMRIS(MRI_SURFACE *mris,const char *fname);
int  MRISreadTransform(MRIS *mris,const char *fname);
int MRISaddToValues(MRI_SURFACE *mris, float val) ;
int MRISsetValues(MRI_SURFACE *mris, float val) ;
LABEL *MRISannotation_to_label(MRI_SURFACE *mris, int annot_index) ;


// thickness stuff
int MRISminimizeThicknessFunctional(MRI_SURFACE *mris,
                                    INTEGRATION_PARMS *parms,
                                    float max_thick) ;
int MRISmeasureDistanceBetweenSurfaces(MRI_SURFACE *mris,
                                       MRI_SURFACE *mris2,
                                       int signed_dist) ;
int MRISwriteCoordsToIco(MRI_SURFACE *mris,
                         MRI_SURFACE *mris_ico,
                         int which_vertices);
int MRISvertexCoord2XYZ_float (VERTEX * v,
                               int which,
                               float  *x, float  *y, float  *z);
int MRISvertexCoord2XYZ_double (VERTEX * v,
                               int which,
                               double  *x, double  *y, double  *z);
int MRISsampleFaceNormal(MRI_SURFACE *mris, int fno, double x, double y, double z, 
                         float *px, float *py, float *pz) ;
int
MRISsampleFaceCoordsCanonical(MHT *mht, MRI_SURFACE *mris, float x, float y, float z, int which, 
                              float *px, float *py, float *pz) ;
MRI *MRISmapToSurface(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, MRI *mri_src_features, MRI *mri_dst_features) ;
int MRISsampleFaceCoords(MRI_SURFACE *mris,
                         int fno,
                         double x, double y, double z,
                         int which_coords,
                         int which_barcyentric,
                         float *px, float *py, float *pz);
MRI *MRISlaplacian(MRI_SURFACE *mris,
                   MRI *mri_cmatrix,
                   double inner_width,
                   double outer_width);
double MRISsampleValue(MRI_SURFACE *mris,
                       FACE *f,
                       double xp, double yp, double zp, 
                       int which, MRI *mri_vals) ;
int MRIScopyAnnotationsToMarkedIndex(MRI_SURFACE *mris) ;
int MRISmaxMarked(MRI_SURFACE *mris) ;
int MRISscaleVertexCoordinates(MRI_SURFACE *mris, double scale) ;
int MRIScurvToMarked(MRI_SURFACE *mris) ;
int MRIScurvToD(MRI_SURFACE *mris) ;
int MRISreadMarked(MRI_SURFACE *mris, const char *sname) ;
int MRISstoreTangentPlanes(MRI_SURFACE *mris, int which_vertices) ;
double MRISsampleFace(MRI_SURFACE *mris, int fno, int which, double x, double y, double z, double val0, double val1, double val2);
int MRISrepositionSurface(MRI_SURFACE *mris, MRI *mri, int *target_vnos, float *target_vals, 
                          int nv, int nsize, double sigma, int flags)  ;
int MRISrepositionSurfaceToCoordinate(MRI_SURFACE *mris, MRI *mri, int target_vno, 
                                      float tx, 
                                      float ty, 
                                      float tz, 
                                      int nsize, double sigma, int flags)  ;
int face_barycentric_coords(MRI_SURFACE const *mris, int fno, int which_vertices,
                            double cx, double cy, double cz, double *pl1, double *pl2, double *pl3) ;

MRI *MRIScomputeFlattenedVolume(MRI_SURFACE *mris,
                                MRI *mri,
                                double res,
                                int nsamples,
                                int normalize,
                                MRI **pmri_vertices,
                                int smooth_iters,
                                double wm_dist,
                                double outside_dist);
int MRIStrinarizeCurvature(MRI_SURFACE *mris, float binarize_thresh) ;
int MRISthresholdValIntoMarked(MRI_SURFACE *mris, float thresh) ;
int MRISremoveCompressedRegions(MRI_SURFACE *mris, double min_dist) ;
int  MRISweightedSoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs) ;
int MRIStaubinSmooth(MRI_SURFACE *mris, int niters, double lambda, double mu, int which) ;
MRI_SURFACE *MRISconcat(MRI_SURFACE *mris1, MRI_SURFACE *mris2, MRI_SURFACE *mris) ;

#define TAUBIN_UNIFORM_WEIGHTS   0
#define TAUBIN_INVERSE_WEIGHTS   1
#define TAUBIN_EDGE_WEIGHTS      2

int MRISvertexNormalToVoxelScaled(MRI_SURFACE *mris,
				  VERTEX *v,
				  MRI *mri,
				  double *pnx, double *pny, double *pnz) ;
int
MRISvertexNormalToVoxel(MRI_SURFACE *mris,
			VERTEX *v,
			MRI *mri,
			double *pnx, double *pny, double *pnz) ;
MRI *MRIcomputeLaminarVolumeFractions(MRI_SURFACE *mris, double res, MRI *mri_src, MRI *mri_vfracs) ;
MRIS *MRIStessellate(MRI *mri,  int value, int all_flag);
void TESSaddFace(MRI *mri, int imnr, int i, int j, int f, int prev_flag, int *pface_index, 
		 tface_type *face, int *face_index_table0, int *face_index_table1);
void TESScheckFace(MRI *mri, int im0, int i0, int j0, int im1, int i1,int j1,
		int f, int n, int v_ind, int prev_flag, int all_flag, int value,
		tface_type *face, int *pface_index, int *face_index_table0, 
		   int *face_index_table1,	tvertex_type *vertex);
int TESSaddVertex(MRI *mri, int imnr, int i, int j, int *pvertex_index,  int *vertex_index_table, tvertex_type *vertex);
int TESSfacep(MRI *mri, int im0, int i0, int j0, int im1, int i1, int j1, int value, int all_flag);

#define SURFACE_SMOOTH_STEPS_TO_SIGMA(iter)   (sqrt((double)iter) * M_PI / 2.0)
#define SIGMA_TO_SURFACE_SMOOTH_STEPS(sigma)  SQR(2.0*sigma/M_PI)

SURFHOPLIST *SetSurfHopListAlloc(MRI_SURFACE *Surf, int nHops);
SURFHOPLIST *SetSurfHopList(int CenterVtx, MRI_SURFACE *Surf, int nHops);
int SurfHopListFree(SURFHOPLIST **shl0);
MRI *MRISarN(MRIS *surf, MRI *src, MRI *mask, MRI *arN, int N);
MRI *MRISsmoothKernel(MRIS *surf, MRI *src, MRI *mask, MRI *mrikern, MATRIX *globkern, SURFHOPLIST ***pshl, MRI *out);
int MRISmeasureLaplaceStreamlines(MRI_SURFACE *mris, MRI *mri_laplace, MRI *mri_intensity, MRI *mri_profiles) ;
MRI *MRISsolveLaplaceEquation(MRI_SURFACE *mris, MRI *mri, double res) ;
int MRIScountEdges(MRIS *surf);
int MRISedges(MRIS *surf);
int MRISfixAverageSurf7(MRIS *surf7);
double mrisRmsValError(MRI_SURFACE *mris, MRI *mri);

// for sorting vertices using qsort
typedef struct{
  int vtxno;
  float x,y,z;
} VERTEX_SORT; 
int CompareVertexCoords(const void *v1, const void *v2);
// for sorting faces using qsort
typedef struct{
  int faceno;
  int v0, v1, v2;
} FACE_SORT; 
int CompareFaceVertices(const void *vf1, const void *vf2);
// for making the surface deterministic after decimation
MRIS *MRISsortVertices(MRIS *mris0);


// Create a surface
//
MRIS *MRISclone(MRIS const * mris_src) ;
MRIS* MRISunion(MRIS const * mris, MRIS const * mris2);


// mrisurf_topology needed by more
//
//  Vertices, like Faces, come into existence when the surface is created with a vertex and face count.
//  Edges are implicit (MRI_EDGE is more than just an edge), and are created by telling each of the end vertices that they are neighbors.
//  Faces get associated with three edges associated with three vertices (VERTICES_PER_FACE is 3)
//
#define mrisCheckVertexVertexTopology(_MRIS) true // mrisCheckVertexVertexTopologyWkr(__FILE__,__LINE__,_MRIS,false)
#define mrisCheckVertexFaceTopology(_MRIS)   true // mrisCheckVertexFaceTopologyWkr  (__FILE__,__LINE__,_MRIS,false)
bool mrisCheckVertexVertexTopologyWkr(const char* file, int line, MRIS const * mris, bool always);
bool mrisCheckVertexFaceTopologyWkr  (const char* file, int line, MRIS const * mris, bool always);
                                            // includes a mrisCheckVertexVertexTopology check

//  Vertices
//
int mrisVertexVSize                 (MRIS const * mris, int vno);
static int  mrisVertexNeighborIndex (MRIS const * mris, int vno1, int vno2);
static bool mrisVerticesAreNeighbors(MRIS const * mris, int vno1, int vno2);

void mrisAddEdge   (MRIS* mris, int vno1, int vno2);
void mrisRemoveEdge(MRIS *mris, int vno1, int vno2);

// Neighbourhoods
//
#define MAX_NEIGHBORS (400)
void mrisVertexReplacingNeighbors(MRIS * mris, int vno, int vnum);
void mrisForgetNeighborhoods     (MRIS * mris);

int  MRISresetNeighborhoodSize      (MRIS *mris, int nsize) ;
void MRISsetNeighborhoodSize        (MRIS *mris, int nsize) ;   // Doesn't compute distances if not already present
void MRISsetNeighborhoodSizeAndDist (MRIS *mris, int nsize) ;   // Always computes distances

#define MRIS_MAX_NEIGHBORHOOD_LINKS 50  // bound on nlinks
int  MRISfindNeighborsAtVertex      (MRIS *mris, int vno, int nlinks, size_t listCapacity, int* vlist, int* hops);
    // sets vlist[*] to the neighboring vno
    // sets hops [*] to -1 for [vno] and the number of hops for all entries returned in the vlist
    // returns the number of neighbors

// dist and dist_orig
//      can be freed at any time
// dist is created by calls to MRISsetNeighborhoodSizeAndDist
// dist_orig must be explicitly created
//
void MRISfreeDistsButNotOrig(MRIS *mris);

void MRISmakeDistOrig (MRIS *mris, int vno);                        // makes it the same size as the current VERTEX.dist
void MRISgrowDistOrig (MRIS *mris, int vno, int minimumCapacity);   // same size as current or bigger
void MRISfreeDistOrigs(MRIS *mris);

//  Faces
//
void mrisSetVertexFaceIndex(MRIS *mris, int vno, int fno);
    // is being used outside mrissurf_topology but shouldn't be
    
void setFaceAttachmentDeferred(MRIS* mris, bool to);                                // for performance reasons, defer adding them; or do all the deferred ones
void mrisAttachFaceToEdges   (MRIS* mris, int fno, int vno1, int vno2, int vno3);   // edges must already exist
void mrisAttachFaceToVertices(MRIS* mris, int fno, int vno1, int vno2, int vno3);   // adds any needed edges

void  MRISflipFaceAroundV1(MRIS *mris, int fno);
void  MRISreverseFaceOrder(MRIS *mris);


// Adding missing parts of the topology representation
//
void mrisCompleteTopology(MRIS *mris);
    //
    // The best way to build a tesselation is to simply create all the vertices then ...
    //      setFaceAttachmentDeferred(true)
    //      lots of calls to mrisAttachFaceToVertices()  
    //      setFaceAttachmentDeferred(false)
    //
    // However lots of existing code just adds a misc set of edges and faces by mechanisms
    // then they call this function to add any missing edges and to calculate the mris->avg_nbrs


// Vertices and Faces can have their ripflag set
//
void MRISsetRipInFacesWithRippedVertices(MRIS* mris);


// Removing 
//
void MRISrenumberRemovingRippedFacesAndVertices(MRIS *mris);

void MRISremoveRipped(MRIS *mris);
    // cleans up the neighbours etc. but does not renumber


// Marked
//
bool mrisAnyVertexOfFaceMarked(MRIS *mris, int fno);


// Static function implementations
//
static int mrisVertexNeighborIndex(MRIS const *mris, int vno1, int vno2) {
  cheapAssert(0 <= vno1 && vno1 < mris->nvertices);
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno1];
  if (!vt->v) return -1;    // happens in invalid surfaces 
  int n;
  for (n = 0; n < vt->vnum; n++) {
    if (vt->v[n] == vno2) return n;
  }
  return -1;
}


static bool mrisVerticesAreNeighbors(MRIS const * const mris, int const vno1, int const vno2)
{
  return 0 <= mrisVertexNeighborIndex(mris, vno1, vno2);
}


// Inputs to the metric properties
//
void MRISsetXYZwkr(MRIS *mris, int vno, float x, float y, float z, const char * file, int line, bool* laterTime);
#define MRISsetXYZ(_MRIS,_VNO, _X,_Y,_Z) { \
    static bool _laterTime; \
    MRISsetXYZwkr((_MRIS),(_VNO),(_X),(_Y),(_Z), __FILE__, __LINE__, &_laterTime); \
  }
    //
    // VERTEX:xyz can be set directly, one at a time, or via one of the following operations
    // that iterate across many vertices.
    //
    // This results in all the derived properties (distances, face areas, normals, angles, ...) being invalid
    // until recomputed.  However the use of invalid properties is not yet detected.

void mrisFindMiddleOfGray(MRIS *mris);

int  MRIStranslate (MRIS *mris, float dx, float dy, float dz);
void MRISmoveOrigin(MRIS *mris, float x0, float y0, float z0);
int  MRISscale     (MRIS *mris, double scale);

void MRISblendXYZandTXYZ(MRIS* mris, float xyzScale, float txyzScale);  // x = x*xyzScale + tx*txyzScale  etc.

void mrisDisturbVertices(MRIS *mris, double amount);

MRIS* MRIScenter(MRIS *mris_src, MRIS *mris_dst) ;
void  MRIScenterSphere(MRIS *mris);
void  MRISrecenter(MRIS *mris, int which_move, int which_target) ;

MRIS* MRISprojectOntoTranslatedSphere(MRIS* mris_src, MRIS* mris_dst, 
    double r,
    double x0, double y0, double z0);


// Vals
//
void MRISsetOriginalXYZ(MRIS *mris, int vno, float x, float y, float z);
void MRISsetOriginalXYZfromXYZ(MRIS *mris);

void MRIScheckForNans(MRIS *mris);
