    namespace AllM {
    struct Face : public Repr_Elt {
        typedef AllM::Surface Surface;
        typedef AllM::Vertex  Vertex;
        inline Face                        (                                            );
        inline Face (                        Face const & src                           );
        inline Face (                        Representation* representation, size_t idx );
        int fno     () const { return idx; }

        inline Vertex                v              ( size_t i                           ) const ;
        inline float                 area           (                                    ) const ;
        inline angles_per_triangle_t angle          (                                    ) const ;
        inline angles_per_triangle_t orig_angle     (                                    ) const ;
        inline char                  ripflag        (                                    ) const ;
        inline char                  oripflag       (                                    ) const ;
        inline int                   marked         (                                    ) const ;
        inline PDMATRIX              norm           (                                    ) const ;
        inline A3PDMATRIX            gradNorm       (                                    ) const ;
        
        inline void                  set_v          ( size_t i,                Vertex to       ) ;
        inline void                  set_area       (                           float to       ) ;
        inline void                  set_angle      (           angles_per_triangle_t to       ) ;
        inline void                  set_orig_angle (           angles_per_triangle_t to       ) ;
        inline void                  set_ripflag    (                            char to       ) ;
        inline void                  set_oripflag   (                            char to       ) ;
        inline void                  set_marked     (                             int to       ) ;
        inline void                  set_norm       (                        PDMATRIX to       ) ;
        inline void                  set_gradNorm   (                      A3PDMATRIX to       ) ;
    }; // Face

    struct Vertex : public Repr_Elt {
        typedef AllM::Surface Surface;
        typedef AllM::Face    Face;
        inline Vertex (                                                                   );
        inline Vertex (                        Vertex const & src                         );
        inline Vertex (                        Representation* representation, size_t idx );
        int vno       () const { return idx; }

        // put the pointers before the ints, before the shorts, before uchars, to reduce size
        // the whole fits in much less than one cache line, so further ordering is no use
        inline Face   f                  ( size_t i            ) const ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
        inline size_t n                  ( size_t i            ) const ;  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        inline int    e                  ( size_t i            ) const ;  //  edge state for neighboring vertices                      
        inline Vertex v                  ( size_t i            ) const ;  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        inline short  vnum               (                     ) const ;  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        inline short  v2num              (                     ) const ;  //  number of 1, or 2-hop neighbors                          
        inline short  v3num              (                     ) const ;  //  number of 1,2,or 3-hop neighbors                         
        inline short  vtotal             (                     ) const ;  //  total # of neighbors. copy of vnum.nsizeCur              
        inline short  nsizeMaxClock      (                     ) const ;  //  copy of mris->nsizeMaxClock when v#num                   
        inline uchar  nsizeMax           (                     ) const ;  //  the max nsize that was used to fill in vnum etc          
        inline uchar  nsizeCur           (                     ) const ;  //  index of the current v#num in vtotal                     
        inline uchar  num                (                     ) const ;  //  number of neighboring faces                              
        // managed by MRISfreeDists[_orig] and MRISmakeDists[_orig]
        inline float  dist               ( size_t i            ) const ;  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        inline float  dist_orig          ( size_t i            ) const ;  // size() is vtotal.    distance to neighboring vertices based on origxyz
        inline int    dist_capacity      (                     ) const ;  //  -- should contain at least vtx_vtotal elements   
        inline int    dist_orig_capacity (                     ) const ;  //  -- should contain at least vtx_vtotal elements   
        // 
        inline float  x                  (                     ) const ;  //  current coordinates	
        inline float  y                  (                     ) const ;  //  use MRISsetXYZ() to set
        inline float  z                  (                     ) const ;
        // 
        inline float  origx              (                     ) const ;  //  original coordinates, see also MRIS::origxyz_status
        inline float  origy              (                     ) const ;  //  use MRISsetOriginalXYZ(, 
        inline float  origz              (                     ) const ;  //  or MRISsetOriginalXYZfromXYZ to set
        // 
        inline float  nx                 (                     ) const ;
        inline float  ny                 (                     ) const ;
        inline float  nz                 (                     ) const ;  //  curr normal
        inline float  pnx                (                     ) const ;
        inline float  pny                (                     ) const ;
        inline float  pnz                (                     ) const ;  //  pial normal
        // 
        inline float  wnx                (                     ) const ;
        inline float  wny                (                     ) const ;
        inline float  wnz                (                     ) const ;  //  white normal
        inline float  onx                (                     ) const ;
        inline float  ony                (                     ) const ;
        inline float  onz                (                     ) const ;  //  original normal
        inline float  dx                 (                     ) const ;
        inline float  dy                 (                     ) const ;
        inline float  dz                 (                     ) const ;  //  current change in position
        inline float  odx                (                     ) const ;
        inline float  ody                (                     ) const ;
        inline float  odz                (                     ) const ;  //  last change of position (for momentum, 
        inline float  tdx                (                     ) const ;
        inline float  tdy                (                     ) const ;
        inline float  tdz                (                     ) const ;  //  temporary storage for averaging gradient
        inline float  curv               (                     ) const ;  //  curr curvature
        inline float  curvbak            (                     ) const ;
        inline float  val                (                     ) const ;  //  scalar data value (file: rh.val, sig2-rh.w)
        inline float  imag_val           (                     ) const ;  //  imaginary part of complex data value
        inline float  cx                 (                     ) const ;
        inline float  cy                 (                     ) const ;
        inline float  cz                 (                     ) const ;  //  coordinates in canonical coordinate system
        inline float  tx                 (                     ) const ;
        inline float  ty                 (                     ) const ;
        inline float  tz                 (                     ) const ;  //  tmp coordinate storage
        inline float  t2x                (                     ) const ;
        inline float  t2y                (                     ) const ;
        inline float  t2z                (                     ) const ;  //  another tmp coordinate storage
        inline float  targx              (                     ) const ;
        inline float  targy              (                     ) const ;
        inline float  targz              (                     ) const ;  //  target coordinates
        inline float  pialx              (                     ) const ;
        inline float  pialy              (                     ) const ;
        inline float  pialz              (                     ) const ;  //  pial surface coordinates
        inline float  whitex             (                     ) const ;
        inline float  whitey             (                     ) const ;
        inline float  whitez             (                     ) const ;  //  white surface coordinates
        inline float  l4x                (                     ) const ;
        inline float  l4y                (                     ) const ;
        inline float  l4z                (                     ) const ;  //  layerIV surface coordinates
        inline float  infx               (                     ) const ;
        inline float  infy               (                     ) const ;
        inline float  infz               (                     ) const ;  //  inflated coordinates
        inline float  fx                 (                     ) const ;
        inline float  fy                 (                     ) const ;
        inline float  fz                 (                     ) const ;  //  flattened coordinates
        inline int    px                 (                     ) const ;
        inline int    qx                 (                     ) const ;
        inline int    py                 (                     ) const ;
        inline int    qy                 (                     ) const ;
        inline int    pz                 (                     ) const ;
        inline int    qz                 (                     ) const ;  //  rational coordinates for exact calculations
        inline float  e1x                (                     ) const ;
        inline float  e1y                (                     ) const ;
        inline float  e1z                (                     ) const ;  //  1st basis vector for the local tangent plane
        inline float  e2x                (                     ) const ;
        inline float  e2y                (                     ) const ;
        inline float  e2z                (                     ) const ;  //  2nd basis vector for the local tangent plane
        inline float  pe1x               (                     ) const ;
        inline float  pe1y               (                     ) const ;
        inline float  pe1z               (                     ) const ;  //  1st basis vector for the local tangent plane
        inline float  pe2x               (                     ) const ;
        inline float  pe2y               (                     ) const ;
        inline float  pe2z               (                     ) const ;  //  2nd basis vector for the local tangent plane
        inline float  nc                 (                     ) const ;  //  curr length normal comp 
        inline float  val2               (                     ) const ;  //  complex comp data value (file: sig3-rh.w) 
        inline float  valbak             (                     ) const ;  //  scalar data stack 
        inline float  val2bak            (                     ) const ;  //  complex comp data stack 
        inline float  stat               (                     ) const ;  //  statistic 
        // 
        inline int    undefval           (                     ) const ;  //  [previously dist=0] 
        inline int    old_undefval       (                     ) const ;  //  for smooth_val_sparse 
        inline int    fixedval           (                     ) const ;  //  [previously val=0] 
        // 
        inline float  fieldsign          (                     ) const ;  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        inline float  fsmask             (                     ) const ;  //  significance mask (file: rh.fm) 
        inline float  d                  (                     ) const ;  //  for distance calculations 
        inline int    annotation         (                     ) const ;  //  area label (defunct--now from label file name!) 
        inline char   oripflag           (                     ) const ;
        inline char   origripflag        (                     ) const ;  //  cuts flags 
        inline p_void vp                 (                     ) const ;  //  to store user's information 
        inline float  theta              (                     ) const ;
        inline float  phi                (                     ) const ;  //  parameterization 
        inline float  area               (                     ) const ;
        inline float  origarea           (                     ) const ;
        inline float  group_avg_area     (                     ) const ;
        inline float  K                  (                     ) const ;  //  Gaussian curvature 
        inline float  H                  (                     ) const ;  //  mean curvature 
        inline float  k1                 (                     ) const ;
        inline float  k2                 (                     ) const ;  //  the principal curvatures 
        inline float  mean               (                     ) const ;
        inline float  mean_imag          (                     ) const ;  //  imaginary part of complex statistic 
        inline float  std_error          (                     ) const ;
        inline uint   flags              (                     ) const ;
        inline int    fno                (                     ) const ;  //  face that this vertex is in 
        inline int    cropped            (                     ) const ;
        inline short  marked             (                     ) const ;  //  for a variety of uses 
        inline short  marked2            (                     ) const ;
        inline short  marked3            (                     ) const ;
        inline char   neg                (                     ) const ;  //  1 if the normal vector is inverted 
        inline char   border             (                     ) const ;  //  flag 
        inline char   ripflag            (                     ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache
        inline void   which_coords       (int which, float *x, float *y, float *z) const ;
        // put the pointers before the ints, before the shorts, before uchars, to reduce size
        // the whole fits in much less than one cache line, so further ordering is no use
        
        inline void   set_f              ( size_t i,   Face to       ) ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
        inline void   set_n              ( size_t i, size_t to       ) ;  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        inline void   set_e              ( size_t i,    int to       ) ;  //  edge state for neighboring vertices                      
        inline void   set_v              ( size_t i, Vertex to       ) ;  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        inline void   set_vnum           (            short to       ) ;  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        inline void   set_v2num          (            short to       ) ;  //  number of 1, or 2-hop neighbors                          
        inline void   set_v3num          (            short to       ) ;  //  number of 1,2,or 3-hop neighbors                         
        inline void   set_vtotal         (            short to       ) ;  //  total # of neighbors. copy of vnum.nsizeCur              
        inline void   set_nsizeMaxClock  (            short to       ) ;  //  copy of mris->nsizeMaxClock when v#num                   
        inline void   set_nsizeMax       (            uchar to       ) ;  //  the max nsize that was used to fill in vnum etc          
        inline void   set_nsizeCur       (            uchar to       ) ;  //  index of the current v#num in vtotal                     
        inline void   set_num            (            uchar to       ) ;  //  number of neighboring faces                              
        // 
        inline void   set_x              (            float to       ) ;  //  current coordinates	
        inline void   set_y              (            float to       ) ;  //  use MRISsetXYZ() to set
        inline void   set_z              (            float to       ) ;
        // 
        inline void   set_nx             (            float to       ) ;
        inline void   set_ny             (            float to       ) ;
        inline void   set_nz             (            float to       ) ;  //  curr normal
        inline void   set_pnx            (            float to       ) ;
        inline void   set_pny            (            float to       ) ;
        inline void   set_pnz            (            float to       ) ;  //  pial normal
        // 
        inline void   set_wnx            (            float to       ) ;
        inline void   set_wny            (            float to       ) ;
        inline void   set_wnz            (            float to       ) ;  //  white normal
        inline void   set_onx            (            float to       ) ;
        inline void   set_ony            (            float to       ) ;
        inline void   set_onz            (            float to       ) ;  //  original normal
        inline void   set_dx             (            float to       ) ;
        inline void   set_dy             (            float to       ) ;
        inline void   set_dz             (            float to       ) ;  //  current change in position
        inline void   set_odx            (            float to       ) ;
        inline void   set_ody            (            float to       ) ;
        inline void   set_odz            (            float to       ) ;  //  last change of position (for momentum, 
        inline void   set_tdx            (            float to       ) ;
        inline void   set_tdy            (            float to       ) ;
        inline void   set_tdz            (            float to       ) ;  //  temporary storage for averaging gradient
        inline void   set_curv           (            float to       ) ;  //  curr curvature
        inline void   set_curvbak        (            float to       ) ;
        inline void   set_val            (            float to       ) ;  //  scalar data value (file: rh.val, sig2-rh.w)
        inline void   set_imag_val       (            float to       ) ;  //  imaginary part of complex data value
        inline void   set_cx             (            float to       ) ;
        inline void   set_cy             (            float to       ) ;
        inline void   set_cz             (            float to       ) ;  //  coordinates in canonical coordinate system
        inline void   set_tx             (            float to       ) ;
        inline void   set_ty             (            float to       ) ;
        inline void   set_tz             (            float to       ) ;  //  tmp coordinate storage
        inline void   set_t2x            (            float to       ) ;
        inline void   set_t2y            (            float to       ) ;
        inline void   set_t2z            (            float to       ) ;  //  another tmp coordinate storage
        inline void   set_targx          (            float to       ) ;
        inline void   set_targy          (            float to       ) ;
        inline void   set_targz          (            float to       ) ;  //  target coordinates
        inline void   set_pialx          (            float to       ) ;
        inline void   set_pialy          (            float to       ) ;
        inline void   set_pialz          (            float to       ) ;  //  pial surface coordinates
        inline void   set_whitex         (            float to       ) ;
        inline void   set_whitey         (            float to       ) ;
        inline void   set_whitez         (            float to       ) ;  //  white surface coordinates
        inline void   set_l4x            (            float to       ) ;
        inline void   set_l4y            (            float to       ) ;
        inline void   set_l4z            (            float to       ) ;  //  layerIV surface coordinates
        inline void   set_infx           (            float to       ) ;
        inline void   set_infy           (            float to       ) ;
        inline void   set_infz           (            float to       ) ;  //  inflated coordinates
        inline void   set_fx             (            float to       ) ;
        inline void   set_fy             (            float to       ) ;
        inline void   set_fz             (            float to       ) ;  //  flattened coordinates
        inline void   set_px             (              int to       ) ;
        inline void   set_qx             (              int to       ) ;
        inline void   set_py             (              int to       ) ;
        inline void   set_qy             (              int to       ) ;
        inline void   set_pz             (              int to       ) ;
        inline void   set_qz             (              int to       ) ;  //  rational coordinates for exact calculations
        inline void   set_e1x            (            float to       ) ;
        inline void   set_e1y            (            float to       ) ;
        inline void   set_e1z            (            float to       ) ;  //  1st basis vector for the local tangent plane
        inline void   set_e2x            (            float to       ) ;
        inline void   set_e2y            (            float to       ) ;
        inline void   set_e2z            (            float to       ) ;  //  2nd basis vector for the local tangent plane
        inline void   set_pe1x           (            float to       ) ;
        inline void   set_pe1y           (            float to       ) ;
        inline void   set_pe1z           (            float to       ) ;  //  1st basis vector for the local tangent plane
        inline void   set_pe2x           (            float to       ) ;
        inline void   set_pe2y           (            float to       ) ;
        inline void   set_pe2z           (            float to       ) ;  //  2nd basis vector for the local tangent plane
        inline void   set_nc             (            float to       ) ;  //  curr length normal comp 
        inline void   set_val2           (            float to       ) ;  //  complex comp data value (file: sig3-rh.w) 
        inline void   set_valbak         (            float to       ) ;  //  scalar data stack 
        inline void   set_val2bak        (            float to       ) ;  //  complex comp data stack 
        inline void   set_stat           (            float to       ) ;  //  statistic 
        // 
        inline void   set_undefval       (              int to       ) ;  //  [previously dist=0] 
        inline void   set_old_undefval   (              int to       ) ;  //  for smooth_val_sparse 
        inline void   set_fixedval       (              int to       ) ;  //  [previously val=0] 
        // 
        inline void   set_fieldsign      (            float to       ) ;  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        inline void   set_fsmask         (            float to       ) ;  //  significance mask (file: rh.fm) 
        inline void   set_d              (            float to       ) ;  //  for distance calculations 
        inline void   set_annotation     (              int to       ) ;  //  area label (defunct--now from label file name!) 
        inline void   set_oripflag       (             char to       ) ;
        inline void   set_origripflag    (             char to       ) ;  //  cuts flags 
        inline void   set_vp             (           p_void to       ) ;  //  to store user's information 
        inline void   set_theta          (            float to       ) ;
        inline void   set_phi            (            float to       ) ;  //  parameterization 
        inline void   set_area           (            float to       ) ;
        inline void   set_origarea       (            float to       ) ;
        inline void   set_group_avg_area (            float to       ) ;
        inline void   set_K              (            float to       ) ;  //  Gaussian curvature 
        inline void   set_H              (            float to       ) ;  //  mean curvature 
        inline void   set_k1             (            float to       ) ;
        inline void   set_k2             (            float to       ) ;  //  the principal curvatures 
        inline void   set_mean           (            float to       ) ;
        inline void   set_mean_imag      (            float to       ) ;  //  imaginary part of complex statistic 
        inline void   set_std_error      (            float to       ) ;
        inline void   set_flags          (             uint to       ) ;
        inline void   set_fno            (              int to       ) ;  //  face that this vertex is in 
        inline void   set_cropped        (              int to       ) ;
        inline void   set_marked         (            short to       ) ;  //  for a variety of uses 
        inline void   set_marked2        (            short to       ) ;
        inline void   set_marked3        (            short to       ) ;
        inline void   set_neg            (             char to       ) ;  //  1 if the normal vector is inverted 
        inline void   set_border         (             char to       ) ;  //  flag 
        inline void   set_ripflag        (             char to       ) ;  //  vertex no longer exists - placed last to load the next vertex into cache
    }; // Vertex

    struct Surface : public Repr_Elt {
        typedef AllM::Face    Face;
        typedef AllM::Vertex  Vertex;
        inline Surface (                                );
        inline Surface ( Surface const & src            );
        inline Surface ( Representation* representation );

        // Fields being maintained by specialist functions
        inline int                   nverticesFrozen          (                               ) const ;  //  # of vertices on surface is frozen
        inline int                   nvertices                (                               ) const ;  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        inline int                   nfaces                   (                               ) const ;  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        inline bool                  faceAttachmentDeferred   (                               ) const ;  //  defer connecting faces to vertices for performance reasons
        inline int                   nedges                   (                               ) const ;  //  # of edges on surface
        inline int                   ncorners                 (                               ) const ;  //  # of triangle corners
        inline int                   nstrips                  (                               ) const ;
        inline Vertex                vertices                 ( size_t i                      ) const ;
        inline p_p_void              dist_storage             (                               ) const ;  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        inline p_p_void              dist_orig_storage        (                               ) const ;  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        inline int                   tempsAssigned            (                               ) const ;  //  State of various temp fields that can be borrowed if not already in use
        inline Face                  faces                    ( size_t i                      ) const ;
        inline MRI_EDGE              edges                    ( size_t i                      ) const ;
        inline MRI_CORNER            corners                  ( size_t i                      ) const ;
        inline FaceNormCacheEntry    faceNormCacheEntries     ( size_t i                      ) const ;
        inline FaceNormDeferredEntry faceNormDeferredEntries  ( size_t i                      ) const ;
        inline STRIP                 strips                   ( size_t i                      ) const ;
        inline float                 xctr                     (                               ) const ;
        inline float                 yctr                     (                               ) const ;
        inline float                 zctr                     (                               ) const ;
        inline float                 xlo                      (                               ) const ;
        inline float                 ylo                      (                               ) const ;
        inline float                 zlo                      (                               ) const ;
        inline float                 xhi                      (                               ) const ;
        inline float                 yhi                      (                               ) const ;
        inline float                 zhi                      (                               ) const ;
        inline float                 x0                       (                               ) const ;  //  center of spherical expansion
        inline float                 y0                       (                               ) const ;
        inline float                 z0                       (                               ) const ;
        // v_temporal_pole, v_frontal_pole, and v_occipital_pole don't appear to be used, and are unusual being pointers to vertices
        // 
        inline float                 max_curv                 (                               ) const ;
        inline float                 min_curv                 (                               ) const ;
        inline float                 total_area               (                               ) const ;
        inline double                avg_vertex_area          (                               ) const ;
        inline double                avg_vertex_dist          (                               ) const ;  //  set by MRIScomputeAvgInterVertexDist
        inline double                std_vertex_dist          (                               ) const ;
        inline float                 orig_area                (                               ) const ;
        inline float                 neg_area                 (                               ) const ;
        inline float                 neg_orig_area            (                               ) const ;  //  amount of original surface in folds
        inline int                   zeros                    (                               ) const ;
        inline int                   hemisphere               (                               ) const ;  //  which hemisphere
        inline int                   initialized              (                               ) const ;
        inline PLTA                  lta                      (                               ) const ;
        inline PMATRIX               SRASToTalSRAS_           (                               ) const ;
        inline PMATRIX               TalSRASToSRAS_           (                               ) const ;
        inline int                   free_transform           (                               ) const ;
        inline double                radius                   (                               ) const ;  //  radius (if status==MRIS_SPHERE)
        inline float                 a                        (                               ) const ;
        inline float                 b                        (                               ) const ;
        inline float                 c                        (                               ) const ;  //  ellipsoid parameters
        inline MRIS_fname_t          fname                    (                               ) const ;  //  file it was originally loaded from
        inline float                 Hmin                     (                               ) const ;  //  min mean curvature
        inline float                 Hmax                     (                               ) const ;  //  max mean curvature
        inline float                 Kmin                     (                               ) const ;  //  min Gaussian curvature
        inline float                 Kmax                     (                               ) const ;  //  max Gaussian curvature
        inline double                Ktotal                   (                               ) const ;  //  total Gaussian curvature
        inline MRIS_Status           status                   (                               ) const ;  //  type of surface (e.g. sphere, plane)
        inline MRIS_Status           origxyz_status           (                               ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline int                   patch                    (                               ) const ;  //  if a patch of the surface
        inline int                   nlabels                  (                               ) const ;
        inline PMRIS_AREA_LABEL      labels                   (                               ) const ;  //  nlabels of these (may be null)
        inline p_void                vp                       (                               ) const ;  //  for misc. use
        inline float                 alpha                    (                               ) const ;  //  rotation around z-axis
        inline float                 beta                     (                               ) const ;  //  rotation around y-axis
        inline float                 gamma                    (                               ) const ;  //  rotation around x-axis
        inline float                 da                       (                               ) const ;
        inline float                 db                       (                               ) const ;
        inline float                 dg                       (                               ) const ;  //  old deltas
        inline int                   type                     (                               ) const ;  //  what type of surface was this initially
        inline int                   max_vertices             (                               ) const ;  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        inline int                   max_faces                (                               ) const ;  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        inline MRIS_subject_name_t   subject_name             (                               ) const ;  //  name of the subject
        inline float                 canon_area               (                               ) const ;
        inline int                   noscale                  (                               ) const ;  //  don't scale by surface area if true
        inline float                 dx2                      ( size_t i                      ) const ;  //  an extra set of gradient (not always alloced)
        inline float                 dy2                      ( size_t i                      ) const ;
        inline float                 dz2                      ( size_t i                      ) const ;
        inline PCOLOR_TABLE          ct                       (                               ) const ;
        inline int                   orig_xyzspace            (                               ) const ;  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        inline int                   useRealRAS               (                               ) const ;  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        inline VOL_GEOM              vg                       (                               ) const ;  //  volume info from which this surface is created. valid iff vg.valid = 1
        inline MRIS_cmdlines_t       cmdlines                 (                               ) const ;
        inline int                   ncmds                    (                               ) const ;
        inline float                 group_avg_surface_area   (                               ) const ;  //  average of total surface area for group
        inline int                   group_avg_vtxarea_loaded (                               ) const ;  //  average vertex area for group at each vertex
        inline int                   triangle_links_removed   (                               ) const ;  //  for quad surfaces
        inline p_void                user_parms               (                               ) const ;  //  for whatever the user wants to hang here
        inline PMATRIX               m_sras2vox               (                               ) const ;  //  for converting surface ras to voxel
        inline PMRI                  mri_sras2vox             (                               ) const ;  //  volume that the above matrix is for
        inline p_void                mht                      (                               ) const ;
        inline p_void                temps                    (                               ) const ;
        
        inline void                  set_strips               ( size_t i,            STRIP to       ) ;
        inline void                  set_xctr                 (                      float to       ) ;
        inline void                  set_yctr                 (                      float to       ) ;
        inline void                  set_zctr                 (                      float to       ) ;
        inline void                  set_xlo                  (                      float to       ) ;
        inline void                  set_ylo                  (                      float to       ) ;
        inline void                  set_zlo                  (                      float to       ) ;
        inline void                  set_xhi                  (                      float to       ) ;
        inline void                  set_yhi                  (                      float to       ) ;
        inline void                  set_zhi                  (                      float to       ) ;
        inline void                  set_x0                   (                      float to       ) ;  //  center of spherical expansion
        inline void                  set_y0                   (                      float to       ) ;
        inline void                  set_z0                   (                      float to       ) ;
        // v_temporal_pole, v_frontal_pole, and v_occipital_pole don't appear to be used, and are unusual being pointers to vertices
        // 
        inline void                  set_max_curv             (                      float to       ) ;
        inline void                  set_min_curv             (                      float to       ) ;
        inline void                  set_total_area           (                      float to       ) ;
        inline void                  set_avg_vertex_area      (                     double to       ) ;
        inline void                  set_avg_vertex_dist      (                     double to       ) ;  //  set by MRIScomputeAvgInterVertexDist
        inline void                  set_std_vertex_dist      (                     double to       ) ;
        inline void                  set_orig_area            (                      float to       ) ;
        inline void                  set_neg_area             (                      float to       ) ;
        inline void                  set_neg_orig_area        (                      float to       ) ;  //  amount of original surface in folds
        inline void                  set_zeros                (                        int to       ) ;
        inline void                  set_hemisphere           (                        int to       ) ;  //  which hemisphere
        inline void                  set_fname                (               MRIS_fname_t to       ) ;  //  file it was originally loaded from
        inline void                  set_Hmin                 (                      float to       ) ;  //  min mean curvature
        inline void                  set_Hmax                 (                      float to       ) ;  //  max mean curvature
        inline void                  set_Kmin                 (                      float to       ) ;  //  min Gaussian curvature
        inline void                  set_Kmax                 (                      float to       ) ;  //  max Gaussian curvature
        inline void                  set_Ktotal               (                     double to       ) ;  //  total Gaussian curvature
        inline void                  set_status               (                MRIS_Status to       ) ;  //  type of surface (e.g. sphere, plane)
        inline void                  set_origxyz_status       (                MRIS_Status to       ) ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline void                  set_patch                (                        int to       ) ;  //  if a patch of the surface
        inline void                  set_nlabels              (                        int to       ) ;
        inline void                  set_labels               (           PMRIS_AREA_LABEL to       ) ;  //  nlabels of these (may be null)
        inline void                  set_vp                   (                     p_void to       ) ;  //  for misc. use
        inline void                  set_alpha                (                      float to       ) ;  //  rotation around z-axis
        inline void                  set_beta                 (                      float to       ) ;  //  rotation around y-axis
        inline void                  set_gamma                (                      float to       ) ;  //  rotation around x-axis
        inline void                  set_da                   (                      float to       ) ;
        inline void                  set_db                   (                      float to       ) ;
        inline void                  set_dg                   (                      float to       ) ;  //  old deltas
        inline void                  set_type                 (                        int to       ) ;  //  what type of surface was this initially
    }; // Surface

    } // namespace AllM
