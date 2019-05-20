    namespace XYZPositionConsequences {
    struct Face : public MRIS_Elt {
        inline Face (                             );
        inline Face ( Face const & src            );
        inline Face ( MRIS* mris, size_t idx      );
        inline Face ( DistortM::Face const & src  );
        inline Face ( Distort::Face const & src   );
        inline Face ( AnalysisM::Face const & src );
        inline Face ( Analysis::Face const & src  );
        inline Face ( AllM::Face const & src      );

        inline vertices_per_face_t   v          (   ) const ;
        inline float                 area       (   ) const ;
        inline angles_per_triangle_t angle      (   ) const ;
        inline angles_per_triangle_t orig_angle (   ) const ;
        inline char                  ripflag    (   ) const ;
        inline char                  oripflag   (   ) const ;
        inline int                   marked     (   ) const ;
        inline PDMATRIX              norm       (   ) const ;
        inline A3PDMATRIX            gradNorm   (   ) const ;
                   
        inline void set_orig_angle (  angles_per_triangle_t to ) ;
        inline void set_ripflag    (                   char to ) ;
        inline void set_oripflag   (                   char to ) ;
        inline void set_marked     (                    int to ) ;
    };

    struct Vertex : public MRIS_Elt {
        inline Vertex (                               );
        inline Vertex ( Vertex const & src            );
        inline Vertex ( MRIS* mris, size_t idx        );
        inline Vertex ( DistortM::Vertex const & src  );
        inline Vertex ( Distort::Vertex const & src   );
        inline Vertex ( AnalysisM::Vertex const & src );
        inline Vertex ( Analysis::Vertex const & src  );
        inline Vertex ( AllM::Vertex const & src      );

        // put the pointers before the ints, before the shorts, before uchars, to reduce size
        // the whole fits in much less than one cache line, so further ordering is no use
        inline Face   f             ( size_t i  ) const ;  // size() is num.    array[v->num] the fno's of the neighboring faces            
        inline size_t n             ( size_t i  ) const ;  // size() is num.    array[v->num] the face.v[*] index for this vertex           
        inline int    e             ( size_t i  ) const ;  //  edge state for neighboring vertices                                          
        inline Vertex v             ( size_t i  ) const ;  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        inline short  vnum          (           ) const ;  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i,                     
        inline short  v2num         (           ) const ;  //  number of 1, or 2-hop neighbors                                              
        inline short  v3num         (           ) const ;  //  number of 1,2,or 3-hop neighbors                                             
        inline short  vtotal        (           ) const ;  //  total # of neighbors. copy of vnum.nsizeCur                                  
        inline short  nsizeMaxClock (           ) const ;  //  copy of mris->nsizeMaxClock when v#num                                       
        inline uchar  nsizeMax      (           ) const ;  //  the max nsize that was used to fill in vnum etc                              
        inline uchar  nsizeCur      (           ) const ;  //  index of the current v#num in vtotal                                         
        inline uchar  num           (           ) const ;  //  number of neighboring faces                                                  
        //           
        inline float  x             (           ) const ;  //  current coordinates	                                                         
        inline float  y             (           ) const ;  //  use MRISsetXYZ() to set                                                      
        inline float  z             (           ) const ;                                                                                   
        //           
        inline float  nx            (           ) const ;                                                                                   
        inline float  ny            (           ) const ;                                                                                   
        inline float  nz            (           ) const ;  //  curr normal                                                                  
        inline float  cx            (           ) const ;                                                                                   
        inline float  cy            (           ) const ;                                                                                   
        inline float  cz            (           ) const ;  //  coordinates in canonical coordinate system                                   
        inline float  area          (           ) const ;                                                                                   
        inline char   ripflag       (           ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache     
        // put the pointers before the ints, before the shorts, before uchars, to reduce size
        // the whole fits in much less than one cache line, so further ordering is no use
                   
        inline void set_f       ( size_t i,   Face to ) ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
        inline void set_n       ( size_t i, size_t to ) ;  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        inline void set_e       ( size_t i,    int to )                  ;  //  edge state for neighboring vertices                      
        //         
        inline void set_nx      (            float to )                                                                                 ;
        inline void set_ny      (            float to )                                                                                 ;
        inline void set_nz      (            float to )                                                                ;  //  curr normal
        inline void set_cx      (            float to )                                                                                 ;
        inline void set_cy      (            float to )                                                                                 ;
        inline void set_cz      (            float to )                                 ;  //  coordinates in canonical coordinate system
        inline void set_ripflag (             char to )   ;  //  vertex no longer exists - placed last to load the next vertex into cache
    };

    struct Surface : public MRIS_Elt {
        inline Surface (                                );
        inline Surface ( Surface const & src            );
        inline Surface ( MRIS* mris, size_t idx         );
        inline Surface ( DistortM::Surface const & src  );
        inline Surface ( Distort::Surface const & src   );
        inline Surface ( AnalysisM::Surface const & src );
        inline Surface ( Analysis::Surface const & src  );
        inline Surface ( AllM::Surface const & src      );

        inline STRIP            strips           ( size_t i  ) const ;                                                                               
        inline float            xctr             (           ) const ;                                                                               
        inline float            yctr             (           ) const ;                                                                               
        inline float            zctr             (           ) const ;                                                                               
        inline float            xlo              (           ) const ;                                                                               
        inline float            ylo              (           ) const ;                                                                               
        inline float            zlo              (           ) const ;                                                                               
        inline float            xhi              (           ) const ;                                                                               
        inline float            yhi              (           ) const ;                                                                               
        inline float            zhi              (           ) const ;                                                                               
        inline float            x0               (           ) const ;  //  center of spherical expansion                                            
        inline float            y0               (           ) const ;                                                                               
        inline float            z0               (           ) const ;                                                                               
        inline PVERTEX          v_temporal_pole  (           ) const ;                                                                               
        inline PVERTEX          v_frontal_pole   (           ) const ;                                                                               
        inline PVERTEX          v_occipital_pole (           ) const ;                                                                               
        inline float            max_curv         (           ) const ;                                                                               
        inline float            min_curv         (           ) const ;                                                                               
        inline float            total_area       (           ) const ;                                                                               
        inline double           avg_vertex_area  (           ) const ;                                                                               
        inline double           avg_vertex_dist  (           ) const ;  //  set by MRIScomputeAvgInterVertexDist                                     
        inline double           std_vertex_dist  (           ) const ;                                                                               
        inline float            orig_area        (           ) const ;                                                                               
        inline float            neg_area         (           ) const ;                                                                               
        inline float            neg_orig_area    (           ) const ;  //  amount of original surface in folds                                      
        inline int              zeros            (           ) const ;                                                                               
        inline int              hemisphere       (           ) const ;  //  which hemisphere                                                         
        inline MRIS_fname_t     fname            (           ) const ;  //  file it was originally loaded from                                       
        inline float            Hmin             (           ) const ;  //  min mean curvature                                                       
        inline float            Hmax             (           ) const ;  //  max mean curvature                                                       
        inline float            Kmin             (           ) const ;  //  min Gaussian curvature                                                   
        inline float            Kmax             (           ) const ;  //  max Gaussian curvature                                                   
        inline double           Ktotal           (           ) const ;  //  total Gaussian curvature                                                 
        inline MRIS_Status      status           (           ) const ;  //  type of surface (e.g. sphere, plane)                                     
        inline MRIS_Status      origxyz_status   (           ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline int              patch            (           ) const ;  //  if a patch of the surface                                                
        inline int              nlabels          (           ) const ;                                                                               
        inline PMRIS_AREA_LABEL labels           (           ) const ;  //  nlabels of these (may be null)                                           
        inline p_void           vp               (           ) const ;  //  for misc. use                                                            
        inline float            alpha            (           ) const ;  //  rotation around z-axis                                                   
        inline float            beta             (           ) const ;  //  rotation around y-axis                                                   
        inline float            gamma            (           ) const ;  //  rotation around x-axis                                                   
        inline float            da               (           ) const ;                                                                               
        inline float            db               (           ) const ;                                                                               
        inline float            dg               (           ) const ;  //  old deltas                                                               
        inline int              type             (           ) const ;  //  what type of surface was this initially                                  
                   
        inline void set_strips           ( size_t i,            STRIP to )                                              ;
        inline void set_xctr             (                      float to )                                              ;
        inline void set_yctr             (                      float to )                                              ;
        inline void set_zctr             (                      float to )                                              ;
        inline void set_xlo              (                      float to )                                              ;
        inline void set_ylo              (                      float to )                                              ;
        inline void set_zlo              (                      float to )                                              ;
        inline void set_xhi              (                      float to )                                              ;
        inline void set_yhi              (                      float to )                                              ;
        inline void set_zhi              (                      float to )                                              ;
        inline void set_x0               (                      float to )           ;  //  center of spherical expansion
        inline void set_y0               (                      float to )                                              ;
        inline void set_z0               (                      float to )                                              ;
        inline void set_v_temporal_pole  (                    PVERTEX to )                                              ;
        inline void set_v_frontal_pole   (                    PVERTEX to )                                              ;
        inline void set_v_occipital_pole (                    PVERTEX to )                                              ;
        inline void set_max_curv         (                      float to )                                              ;
        inline void set_min_curv         (                      float to )                                              ;
        inline void set_total_area       (                      float to )                                              ;
        inline void set_avg_vertex_area  (                     double to )                                              ;
        inline void set_zeros            (                        int to )                                              ;
        inline void set_hemisphere       (                        int to )                        ;  //  which hemisphere
        inline void set_Hmin             (                      float to )                      ;  //  min mean curvature
        inline void set_Hmax             (                      float to )                      ;  //  max mean curvature
        inline void set_Kmin             (                      float to )                  ;  //  min Gaussian curvature
        inline void set_Kmax             (                      float to )                  ;  //  max Gaussian curvature
        inline void set_Ktotal           (                     double to )                ;  //  total Gaussian curvature
        inline void set_nlabels          (                        int to )                                              ;
        inline void set_labels           (           PMRIS_AREA_LABEL to )          ;  //  nlabels of these (may be null)
        inline void set_vp               (                     p_void to )                           ;  //  for misc. use
        inline void set_alpha            (                      float to )                  ;  //  rotation around z-axis
        inline void set_beta             (                      float to )                  ;  //  rotation around y-axis
        inline void set_gamma            (                      float to )                  ;  //  rotation around x-axis
        inline void set_da               (                      float to )                                              ;
        inline void set_db               (                      float to )                                              ;
        inline void set_dg               (                      float to )                              ;  //  old deltas
        inline void set_type             (                        int to ) ;  //  what type of surface was this initially
    };

    } // namespace XYZPositionConsequences
