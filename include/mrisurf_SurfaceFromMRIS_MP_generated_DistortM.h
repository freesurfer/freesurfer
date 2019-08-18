    namespace DistortM {
    struct Face : public Repr_Elt {
        typedef DistortM::Surface Surface;
        typedef DistortM::Vertex  Vertex;
        inline Face                        (                                            );
        inline Face (                        Face const & src                           );
        inline Face (                        Representation* representation, size_t idx );
        inline Face (                        AllM::Face const & src                     );
        int fno     () const { return idx; }

        inline Vertex                v              ( size_t i                          ) const ;
        inline float                 area           (                                   ) const ;
        inline angles_per_triangle_t angle          (                                   ) const ;
        inline angles_per_triangle_t orig_angle     (                                   ) const ;
        inline char                  ripflag        (                                   ) const ;
        inline FloatXYZ              norm           (                                   ) const ;
        
        inline void                  set_orig_angle (          angles_per_triangle_t to       ) ;
        inline void                  set_ripflag    (                           char to       ) ;
    }; // Face

    struct Vertex : public Repr_Elt {
        typedef DistortM::Surface Surface;
        typedef DistortM::Face    Face;
        inline Vertex (                                                                   );
        inline Vertex (                        Vertex const & src                         );
        inline Vertex (                        Representation* representation, size_t idx );
        inline Vertex (                        AllM::Vertex const & src                   );
        int vno       () const { return idx; }

        inline float dist          ( size_t i          ) const ;  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        inline float dist_orig     ( size_t i          ) const ;  // size() is vtotal.    distance to neighboring vertices based on origxyz
        inline int   dist_capacity (                   ) const ;  //  -- should contain at least vtx_vtotal elements   
        inline float x             (                   ) const ;  //  current coordinates	
        inline float y             (                   ) const ;  //  use MRISsetXYZ() to set
        inline float z             (                   ) const ;
        inline float origx         (                   ) const ;  //  original coordinates, see also MRIS::origxyz_status
        inline float origy         (                   ) const ;  //  use MRISsetOriginalXYZ(, 
        inline float origz         (                   ) const ;  //  or MRISsetOriginalXYZfromXYZ to set
        inline float nx            (                   ) const ;
        inline float ny            (                   ) const ;
        inline float nz            (                   ) const ;  //  curr normal
        inline float wnx           (                   ) const ;
        inline float wny           (                   ) const ;
        inline float wnz           (                   ) const ;  //  white normal
        inline float dx            (                   ) const ;
        inline float dy            (                   ) const ;
        inline float dz            (                   ) const ;  //  current change in position
        inline float curv          (                   ) const ;  //  curr curvature
        inline float cx            (                   ) const ;
        inline float cy            (                   ) const ;
        inline float cz            (                   ) const ;  //  coordinates in canonical coordinate system
        inline float pialx         (                   ) const ;
        inline float pialy         (                   ) const ;
        inline float pialz         (                   ) const ;  //  pial surface coordinates
        inline float whitex        (                   ) const ;
        inline float whitey        (                   ) const ;
        inline float whitez        (                   ) const ;  //  white surface coordinates
        inline float area          (                   ) const ;
        inline float origarea      (                   ) const ;
        inline int   fno           (                   ) const ;  //  face that this vertex is in 
        inline char  neg           (                   ) const ;  //  1 if the normal vector is inverted 
        inline char  border        (                   ) const ;  //  flag 
        inline char  ripflag       (                   ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache
        inline void  which_coords  (int which, float *x, float *y, float *z) const ;
        
        inline void  set_nx        (          float to       ) ;
        inline void  set_ny        (          float to       ) ;
        inline void  set_nz        (          float to       ) ;  //  curr normal
        inline void  set_wnx       (          float to       ) ;
        inline void  set_wny       (          float to       ) ;
        inline void  set_wnz       (          float to       ) ;  //  white normal
        inline void  set_dx        (          float to       ) ;
        inline void  set_dy        (          float to       ) ;
        inline void  set_dz        (          float to       ) ;  //  current change in position
        inline void  set_curv      (          float to       ) ;  //  curr curvature
        inline void  set_cx        (          float to       ) ;
        inline void  set_cy        (          float to       ) ;
        inline void  set_cz        (          float to       ) ;  //  coordinates in canonical coordinate system
        inline void  set_pialx     (          float to       ) ;
        inline void  set_pialy     (          float to       ) ;
        inline void  set_pialz     (          float to       ) ;  //  pial surface coordinates
        inline void  set_whitex    (          float to       ) ;
        inline void  set_whitey    (          float to       ) ;
        inline void  set_whitez    (          float to       ) ;  //  white surface coordinates
        inline void  set_origarea  (          float to       ) ;
        inline void  set_fno       (            int to       ) ;  //  face that this vertex is in 
        inline void  set_neg       (           char to       ) ;  //  1 if the normal vector is inverted 
        inline void  set_border    (           char to       ) ;  //  flag 
        inline void  set_ripflag   (           char to       ) ;  //  vertex no longer exists - placed last to load the next vertex into cache
    }; // Vertex

    struct MRIS_MP : public Repr_Elt {
        typedef DistortM::Surface Surface;
        typedef DistortM::Face    Face;
        typedef DistortM::Vertex  Vertex;
        inline MRIS_MP (                                            );
        inline MRIS_MP ( MRIS_MP const & src                        );
        inline MRIS_MP ( Representation* representation, size_t idx );
        inline MRIS_MP ( AllM::MRIS_MP const & src                  );

    }; // MRIS_MP

    struct Surface : public Repr_Elt {
        typedef DistortM::Face    Face;
        typedef DistortM::Vertex  Vertex;
        inline Surface                                                (                                );
        inline Surface                                                ( Surface const & src            );
        inline Surface                                                ( Representation* representation );
        inline Surface                                                ( AllM::Surface const & src      );
        void freeDistsButNotOrig() { MRISfreeDistsButNotOrig(repr); }

        inline int                   nvertices               (                    ) const ;  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        inline int                   nfaces                  (                    ) const ;  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        inline Vertex                vertices                ( size_t i           ) const ;
        inline Face                  faces                   ( size_t i           ) const ;
        inline FaceNormCacheEntry    faceNormCacheEntries    ( size_t i           ) const ;
        inline FaceNormDeferredEntry faceNormDeferredEntries ( size_t i           ) const ;
        inline float                 xctr                    (                    ) const ;
        inline float                 yctr                    (                    ) const ;
        inline float                 zctr                    (                    ) const ;
        inline float                 xlo                     (                    ) const ;
        inline float                 ylo                     (                    ) const ;
        inline float                 zlo                     (                    ) const ;
        inline float                 xhi                     (                    ) const ;
        inline float                 yhi                     (                    ) const ;
        inline float                 zhi                     (                    ) const ;
        inline float                 total_area              (                    ) const ;
        inline double                avg_vertex_area         (                    ) const ;
        inline double                avg_vertex_dist         (                    ) const ;  //  set by MRIScomputeAvgInterVertexDist
        inline double                std_vertex_dist         (                    ) const ;
        inline float                 orig_area               (                    ) const ;
        inline float                 neg_area                (                    ) const ;
        inline float                 neg_orig_area           (                    ) const ;  //  amount of original surface in folds
        inline double                radius                  (                    ) const ;  //  radius (if status==MRIS_SPHERE)
        inline MRIS_Status           status                  (                    ) const ;  //  type of surface (e.g. sphere, plane)
        inline MRIS_Status           origxyz_status          (                    ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline int                   patch                   (                    ) const ;  //  if a patch of the surface
        inline int                   noscale                 (                    ) const ;  //  don't scale by surface area if true
        
        inline void                  set_xctr                (           float to       ) ;
        inline void                  set_yctr                (           float to       ) ;
        inline void                  set_zctr                (           float to       ) ;
        inline void                  set_xlo                 (           float to       ) ;
        inline void                  set_ylo                 (           float to       ) ;
        inline void                  set_zlo                 (           float to       ) ;
        inline void                  set_xhi                 (           float to       ) ;
        inline void                  set_yhi                 (           float to       ) ;
        inline void                  set_zhi                 (           float to       ) ;
        inline void                  set_total_area          (           float to       ) ;
        inline void                  set_avg_vertex_area     (          double to       ) ;
    }; // Surface

    } // namespace DistortM
