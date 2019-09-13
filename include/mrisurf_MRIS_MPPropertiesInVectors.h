
#pragma once
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
#include "mrisurf_aaa.h"
struct MRIS_MP {
                            MRIS* underlyingMRIS          ;  //  for properties that are read from the underlying MRIS
    MRIS_MP*                      in_src                  ;  //  since the in are not written, they can be shared by copies
    int                           in_ref_count            ;  //  check the src doesn't go away
    VERTEX_TOPOLOGY const *       vertices_topology       ;  //  pointer copied from MRIS
    FACE_TOPOLOGY   const *       faces_topology          ;  //  pointer copied from MRIS
    int*                          v_VSize                 ;
    float**                       v_dist_buffer           ;
    const float*                  f_norm_orig_area        ;
    char*                         f_normSet               ;
    vertices_per_face_t*          f_v                     ;
    float*                        f_area                  ;
    angles_per_triangle_t*        f_angle                 ;
    angles_per_triangle_t*        f_orig_angle            ;
    char*                         f_ripflag               ;
    FloatXYZ*                     f_norm                  ;
    pSeveralFloat*                v_dist                  ;  // size() is vtotal.    distance to neighboring vertices based on  xyz   
    pSeveralFloat*                v_dist_orig             ;  // size() is vtotal.    distance to neighboring vertices based on origxyz
    int*                          v_dist_capacity         ;  //  -- should contain at least vtx_vtotal elements   
    float*                        v_x                     ;  //  current coordinates	
    float*                        v_y                     ;  //  use MRISsetXYZ() to set
    float*                        v_z                     ;
    float*                        v_origx                 ;  //  original coordinates, see also MRIS::origxyz_status
    float*                        v_origy                 ;  //  use MRISsetOriginalXYZ(, 
    float*                        v_origz                 ;  //  or MRISsetOriginalXYZfromXYZ to set
    float*                        v_nx                    ;
    float*                        v_ny                    ;
    float*                        v_nz                    ;  //  curr normal
    float*                        v_wnx                   ;
    float*                        v_wny                   ;
    float*                        v_wnz                   ;  //  white normal
    float*                        v_dx                    ;
    float*                        v_dy                    ;
    float*                        v_dz                    ;  //  current change in position
    float*                        v_curv                  ;  //  curr curvature
    float*                        v_cx                    ;
    float*                        v_cy                    ;
    float*                        v_cz                    ;  //  coordinates in canonical coordinate system
    float*                        v_pialx                 ;
    float*                        v_pialy                 ;
    float*                        v_pialz                 ;  //  pial surface coordinates
    float*                        v_whitex                ;
    float*                        v_whitey                ;
    float*                        v_whitez                ;  //  white surface coordinates
    float*                        v_area                  ;
    float*                        v_origarea              ;
    int*                          v_fno                   ;  //  face that this vertex is in 
    char*                         v_neg                   ;  //  1 if the normal vector is inverted 
    char*                         v_border                ;  //  flag 
    char*                         v_ripflag               ;  //  vertex no longer exists - placed last to load the next vertex into cache
    int                           nvertices               ;  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
    int                           nfaces                  ;  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
    pSeveralVERTEX                vertices                ;
    pSeveralFACE                  faces                   ;
    pSeveralFaceNormCacheEntry    faceNormCacheEntries    ;
    pSeveralFaceNormDeferredEntry faceNormDeferredEntries ;
    float                         xctr                    ;
    float                         yctr                    ;
    float                         zctr                    ;
    float                         xlo                     ;
    float                         ylo                     ;
    float                         zlo                     ;
    float                         xhi                     ;
    float                         yhi                     ;
    float                         zhi                     ;
    float                         total_area              ;
    double                        avg_vertex_area         ;
    double                        avg_vertex_dist         ;  //  set by MRIScomputeAvgInterVertexDist
    double                        std_vertex_dist         ;
    float                         orig_area               ;
    float                         neg_area                ;
    float                         neg_orig_area           ;  //  amount of original surface in folds
    double                        radius                  ;  //  radius (if status==MRIS_SPHERE)
    MRIS_Status                   status                  ;  //  type of surface (e.g. sphere, plane)
    MRIS_Status                   origxyz_status          ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
    int                           patch                   ;  //  if a patch of the surface
    char                          nsize                   ;  //  size of neighborhoods or -1
    char                          dist_nsize              ;  //  max mrisComputeVertexDistances has computed distances out to
    int                           noscale                 ;  //  don't scale by surface area if true
};		// MRIS_MP

