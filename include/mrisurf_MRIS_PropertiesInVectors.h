
#pragma once
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
#include "mrisurf_aaa.h"
struct MRISPV {
             vertices_per_face_t* f_v                      ;
    float*                        f_area                   ;
    angles_per_triangle_t*        f_angle                  ;
    angles_per_triangle_t*        f_orig_angle             ;
    char*                         f_ripflag                ;
    char*                         f_oripflag               ;
    int*                          f_marked                 ;
    PDMATRIX*                     f_norm                   ;
    A3PDMATRIX*                   f_gradNorm               ;
    //  put the pointers before the ints, before the shorts, before uchars, to reduce size
    //  the whole fits in much less than one cache line, so further ordering is no use
    pSeveralInt*                  v_f                      ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
    pSeveralUchar*                v_n                      ;  // size() is num.    array[v->num] the face.v[*] index for this vertex        
    pSeveralInt*                  v_e                      ;  //  edge state for neighboring vertices                      
    pSeveralInt*                  v_v                      ;  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
    short*                        v_vnum                   ;  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
    short*                        v_v2num                  ;  //  number of 1, or 2-hop neighbors                          
    short*                        v_v3num                  ;  //  number of 1,2,or 3-hop neighbors                         
    short*                        v_vtotal                 ;  //  total # of neighbors. copy of vnum.nsizeCur              
    short*                        v_nsizeMaxClock          ;  //  copy of mris->nsizeMaxClock when v#num                   
    uchar*                        v_nsizeMax               ;  //  the max nsize that was used to fill in vnum etc          
    uchar*                        v_nsizeCur               ;  //  index of the current v#num in vtotal                     
    uchar*                        v_num                    ;  //  number of neighboring faces                              
    //  managed by MRISfreeDists[_orig] and MRISmakeDists[_orig]
    pSeveralFloat*                v_dist                   ;  // size() is vtotal.    distance to neighboring vertices based on  xyz   
    pSeveralFloat*                v_dist_orig              ;  // size() is vtotal.    distance to neighboring vertices based on origxyz
    int*                          v_dist_capacity          ;  //  -- should contain at least vtx_vtotal elements   
    int*                          v_dist_orig_capacity     ;  //  -- should contain at least vtx_vtotal elements   
    float*                        v_x                      ;  //  current coordinates	
    float*                        v_y                      ;  //  use MRISsetXYZ() to set
    float*                        v_z                      ;
    float*                        v_origx                  ;  //  original coordinates, see also MRIS::origxyz_status
    float*                        v_origy                  ;  //  use MRISsetOriginalXYZ(, 
    float*                        v_origz                  ;  //  or MRISsetOriginalXYZfromXYZ to set
    float*                        v_nx                     ;
    float*                        v_ny                     ;
    float*                        v_nz                     ;  //  curr normal
    float*                        v_pnx                    ;
    float*                        v_pny                    ;
    float*                        v_pnz                    ;  //  pial normal
    float*                        v_wnx                    ;
    float*                        v_wny                    ;
    float*                        v_wnz                    ;  //  white normal
    float*                        v_onx                    ;
    float*                        v_ony                    ;
    float*                        v_onz                    ;  //  original normal
    float*                        v_dx                     ;
    float*                        v_dy                     ;
    float*                        v_dz                     ;  //  current change in position
    float*                        v_odx                    ;
    float*                        v_ody                    ;
    float*                        v_odz                    ;  //  last change of position (for momentum, 
    float*                        v_tdx                    ;
    float*                        v_tdy                    ;
    float*                        v_tdz                    ;  //  temporary storage for averaging gradient
    float*                        v_curv                   ;  //  curr curvature
    float*                        v_curvbak                ;
    float*                        v_val                    ;  //  scalar data value (file: rh.val, sig2-rh.w)
    float*                        v_imag_val               ;  //  imaginary part of complex data value
    float*                        v_cx                     ;
    float*                        v_cy                     ;
    float*                        v_cz                     ;  //  coordinates in canonical coordinate system
    float*                        v_tx                     ;
    float*                        v_ty                     ;
    float*                        v_tz                     ;  //  tmp coordinate storage
    float*                        v_t2x                    ;
    float*                        v_t2y                    ;
    float*                        v_t2z                    ;  //  another tmp coordinate storage
    float*                        v_targx                  ;
    float*                        v_targy                  ;
    float*                        v_targz                  ;  //  target coordinates
    float*                        v_pialx                  ;
    float*                        v_pialy                  ;
    float*                        v_pialz                  ;  //  pial surface coordinates
    float*                        v_whitex                 ;
    float*                        v_whitey                 ;
    float*                        v_whitez                 ;  //  white surface coordinates
    float*                        v_l4x                    ;
    float*                        v_l4y                    ;
    float*                        v_l4z                    ;  //  layerIV surface coordinates
    float*                        v_infx                   ;
    float*                        v_infy                   ;
    float*                        v_infz                   ;  //  inflated coordinates
    float*                        v_fx                     ;
    float*                        v_fy                     ;
    float*                        v_fz                     ;  //  flattened coordinates
    int*                          v_px                     ;
    int*                          v_qx                     ;
    int*                          v_py                     ;
    int*                          v_qy                     ;
    int*                          v_pz                     ;
    int*                          v_qz                     ;  //  rational coordinates for exact calculations
    float*                        v_e1x                    ;
    float*                        v_e1y                    ;
    float*                        v_e1z                    ;  //  1st basis vector for the local tangent plane
    float*                        v_e2x                    ;
    float*                        v_e2y                    ;
    float*                        v_e2z                    ;  //  2nd basis vector for the local tangent plane
    float*                        v_pe1x                   ;
    float*                        v_pe1y                   ;
    float*                        v_pe1z                   ;  //  1st basis vector for the local tangent plane
    float*                        v_pe2x                   ;
    float*                        v_pe2y                   ;
    float*                        v_pe2z                   ;  //  2nd basis vector for the local tangent plane
    float*                        v_nc                     ;  //  curr length normal comp 
    float*                        v_val2                   ;  //  complex comp data value (file: sig3-rh.w) 
    float*                        v_valbak                 ;  //  scalar data stack 
    float*                        v_val2bak                ;  //  complex comp data stack 
    float*                        v_stat                   ;  //  statistic 
    int*                          v_undefval               ;  //  [previously dist=0] 
    int*                          v_old_undefval           ;  //  for smooth_val_sparse 
    int*                          v_fixedval               ;  //  [previously val=0] 
    float*                        v_fieldsign              ;  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
    float*                        v_fsmask                 ;  //  significance mask (file: rh.fm) 
    float*                        v_d                      ;  //  for distance calculations 
    int*                          v_annotation             ;  //  area label (defunct--now from label file name!) 
    char*                         v_oripflag               ;
    char*                         v_origripflag            ;  //  cuts flags 
    p_void*                       v_vp                     ;  //  to store user's information 
    float*                        v_theta                  ;
    float*                        v_phi                    ;  //  parameterization 
    float*                        v_area                   ;
    float*                        v_origarea               ;
    float*                        v_group_avg_area         ;
    float*                        v_K                      ;  //  Gaussian curvature 
    float*                        v_H                      ;  //  mean curvature 
    float*                        v_k1                     ;
    float*                        v_k2                     ;  //  the principal curvatures 
    float*                        v_mean                   ;
    float*                        v_mean_imag              ;  //  imaginary part of complex statistic 
    float*                        v_std_error              ;
    uint*                         v_flags                  ;
    int*                          v_fno                    ;  //  face that this vertex is in 
    int*                          v_cropped                ;
    short*                        v_marked                 ;  //  for a variety of uses 
    short*                        v_marked2                ;
    short*                        v_marked3                ;
    char*                         v_neg                    ;  //  1 if the normal vector is inverted 
    char*                         v_border                 ;  //  flag 
    char*                         v_ripflag                ;  //  vertex no longer exists - placed last to load the next vertex into cache
    //  Fields being maintained by specialist functions
    int                           nverticesFrozen          ;  //  # of vertices on surface is frozen
    int                           nvertices                ;  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
    int                           nfaces                   ;  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
    bool                          faceAttachmentDeferred   ;  //  defer connecting faces to vertices for performance reasons
    int                           nedges                   ;  //  # of edges on surface
    int                           ncorners                 ;  //  # of triangle corners
    int                           nstrips                  ;
    pSeveralVERTEX_TOPOLOGY       vertices_topology        ;
    pSeveralVERTEX                vertices                 ;
    p_p_void                      dist_storage             ;  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
    p_p_void                      dist_orig_storage        ;  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
    int                           tempsAssigned            ;  //  State of various temp fields that can be borrowed if not already in use
    pSeveralFACE                  faces                    ;
    pSeveralMRI_EDGE              edges                    ;
    pSeveralMRI_CORNER            corners                  ;
    pSeveralFaceNormCacheEntry    faceNormCacheEntries     ;
    pSeveralFaceNormDeferredEntry faceNormDeferredEntries  ;
    pSeveralSTRIP                 strips                   ;
    float                         xctr                     ;
    float                         yctr                     ;
    float                         zctr                     ;
    float                         xlo                      ;
    float                         ylo                      ;
    float                         zlo                      ;
    float                         xhi                      ;
    float                         yhi                      ;
    float                         zhi                      ;
    float                         x0                       ;  //  center of spherical expansion
    float                         y0                       ;
    float                         z0                       ;
    //  v_temporal_pole, v_frontal_pole, and v_occipital_pole don't appear to be used, and are unusual being pointers to vertices
    float                         max_curv                 ;
    float                         min_curv                 ;
    float                         total_area               ;
    double                        avg_vertex_area          ;
    double                        avg_vertex_dist          ;  //  set by MRIScomputeAvgInterVertexDist
    double                        std_vertex_dist          ;
    float                         orig_area                ;
    float                         neg_area                 ;
    float                         neg_orig_area            ;  //  amount of original surface in folds
    int                           zeros                    ;
    int                           hemisphere               ;  //  which hemisphere
    int                           initialized              ;
    PLTA                          lta                      ;
    PMATRIX                       SRASToTalSRAS_           ;
    PMATRIX                       TalSRASToSRAS_           ;
    int                           free_transform           ;
    double                        radius                   ;  //  radius (if status==MRIS_SPHERE)
    float                         a                        ;
    float                         b                        ;
    float                         c                        ;  //  ellipsoid parameters
    MRIS_fname_t                  fname                    ;  //  file it was originally loaded from
    float                         Hmin                     ;  //  min mean curvature
    float                         Hmax                     ;  //  max mean curvature
    float                         Kmin                     ;  //  min Gaussian curvature
    float                         Kmax                     ;  //  max Gaussian curvature
    double                        Ktotal                   ;  //  total Gaussian curvature
    MRIS_Status                   status                   ;  //  type of surface (e.g. sphere, plane)
    MRIS_Status                   origxyz_status           ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
    int                           patch                    ;  //  if a patch of the surface
    int                           nlabels                  ;
    PMRIS_AREA_LABEL              labels                   ;  //  nlabels of these (may be null)
    char                          nsize                    ;  //  size of neighborhoods or -1
    uchar                         vtotalsMightBeTooBig     ;  //  MRISsampleDistances sets this
    short                         nsizeMaxClock            ;  //  changed whenever an edge is added or removed, which invalidates the vertex v#num values
    char                          max_nsize                ;  //  max the neighborhood size has been set to (typically 3)
    char                          dist_nsize               ;  //  max mrisComputeVertexDistances has computed distances out to
    char                          dist_orig_nsize          ;  //  max mrisComputeOriginalVertexDistances has computed distances out to
    char                          dist_alloced_flags       ;  //  two flags, set when any dist(1) or dist_orig(2) allocated
    float                         avg_nbrs                 ;  //  mean # of vertex neighbors
    p_void                        vp                       ;  //  for misc. use
    float                         alpha                    ;  //  rotation around z-axis
    float                         beta                     ;  //  rotation around y-axis
    float                         gamma                    ;  //  rotation around x-axis
    float                         da                       ;
    float                         db                       ;
    float                         dg                       ;  //  old deltas
    int                           type                     ;  //  what type of surface was this initially
    int                           max_vertices             ;  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
    int                           max_faces                ;  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
    MRIS_subject_name_t           subject_name             ;  //  name of the subject
    float                         canon_area               ;
    int                           noscale                  ;  //  don't scale by surface area if true
    pSeveralFloat                 dx2                      ;  //  an extra set of gradient (not always alloced)
    pSeveralFloat                 dy2                      ;
    pSeveralFloat                 dz2                      ;
    PCOLOR_TABLE                  ct                       ;
    int                           orig_xyzspace            ;  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
    int                           useRealRAS               ;  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
    VOL_GEOM                      vg                       ;  //  volume info from which this surface is created. valid iff vg.valid = 1
    MRIS_cmdlines_t               cmdlines                 ;
    int                           ncmds                    ;
    float                         group_avg_surface_area   ;  //  average of total surface area for group
    int                           group_avg_vtxarea_loaded ;  //  average vertex area for group at each vertex
    int                           triangle_links_removed   ;  //  for quad surfaces
    p_void                        user_parms               ;  //  for whatever the user wants to hang here
    PMATRIX                       m_sras2vox               ;  //  for converting surface ras to voxel
    PMRI                          mri_sras2vox             ;  //  volume that the above matrix is for
    p_void                        mht                      ;
    p_void                        temps                    ;
};		// MRISPV

