
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
struct MRISPV {
      vertices_per_face_t * f_v       ;
    float                 * f_area    ;
    angles_per_triangle_t * f_angle   ;
    char                  * f_ripflag ;
    PDMATRIX              * f_norm    ;
    //  put the pointers before the ints, before the shorts, before uchars, to reduce size
    //  the whole fits in much less than one cache line, so further ordering is no use
    pSeveralInt   * v_f             ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
    uchar         * v_num           ;  //  number of neighboring faces                              
    //  managed by MRISfreeDists[_orig] and MRISmakeDists[_orig]
    pSeveralFloat * v_dist          ;  // size() is vtotal.    distance to neighboring vertices based on  xyz   
    int           * v_dist_capacity ;  //  -- should contain at least vtx_vtotal elements   
    float         * v_x             ;  //  current coordinates	
    float         * v_y             ;  //  use MRISsetXYZ() to set
    float         * v_z             ;
    float         * v_nx            ;
    float         * v_ny            ;
    float         * v_nz            ;  //  curr normal
    float         * v_area          ;
    float         * v_origarea      ;
    char          * v_neg           ;  //  1 if the normal vector is inverted 
    char          * v_border        ;  //  flag 
    char          * v_ripflag       ;  //  vertex no longer exists - placed last to load the next vertex into cache
    //  Fields being maintained by specialist functions
    int            nvertices       ;  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
    int            nfaces          ;  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
    float          xctr            ;
    float          yctr            ;
    float          zctr            ;
    float          xlo             ;
    float          ylo             ;
    float          zlo             ;
    float          xhi             ;
    float          yhi             ;
    float          zhi             ;
    float          total_area      ;
    double         avg_vertex_area ;
    double         avg_vertex_dist ;  //  set by MRIScomputeAvgInterVertexDist
    double         std_vertex_dist ;
    float          neg_area        ;
    float          neg_orig_area   ;  //  amount of original surface in folds
    double         radius          ;  //  radius (if status==MRIS_SPHERE)
    MRIS_Status    status          ;  //  type of surface (e.g. sphere, plane)
    MRIS_Status    origxyz_status  ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
    char           nsize           ;  //  size of neighborhoods or -1
    char           dist_nsize      ;  //  max mrisComputeVertexDistances has computed distances out to
};
