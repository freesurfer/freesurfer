    namespace TopologyM {
    struct Face : public Repr_Elt {
        typedef TopologyM::Surface Surface;
        typedef TopologyM::Vertex  Vertex;
        inline Face                        (                                            );
        inline Face (                        Face const & src                           );
        inline Face (                        Representation* representation, size_t idx );
        inline Face (                        AllM::Face const & src                     );
        int fno     () const { return idx; }

        inline Vertex v           ( size_t i            ) const ;
        inline char   ripflag     (                     ) const ;
        
        inline void   set_v       ( size_t i, Vertex to       ) ;
        inline void   set_ripflag (             char to       ) ;
    }; // Face

    struct Vertex : public Repr_Elt {
        typedef TopologyM::Surface Surface;
        typedef TopologyM::Face    Face;
        inline Vertex (                                                                   );
        inline Vertex (                        Vertex const & src                         );
        inline Vertex (                        Representation* representation, size_t idx );
        inline Vertex (                        AllM::Vertex const & src                   );
        int vno       () const { return idx; }

        inline char ripflag      (           ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache
        inline void which_coords (int which, float *x, float *y, float *z) const ;
        
        inline void set_ripflag  (   char to       ) ;  //  vertex no longer exists - placed last to load the next vertex into cache
    }; // Vertex

    struct MRIS_MP : public Repr_Elt {
        typedef TopologyM::Surface Surface;
        typedef TopologyM::Face    Face;
        typedef TopologyM::Vertex  Vertex;
        inline MRIS_MP (                                            );
        inline MRIS_MP ( MRIS_MP const & src                        );
        inline MRIS_MP ( Representation* representation, size_t idx );
        inline MRIS_MP ( AllM::MRIS_MP const & src                  );

    }; // MRIS_MP

    struct Surface : public Repr_Elt {
        typedef TopologyM::Face    Face;
        typedef TopologyM::Vertex  Vertex;
        inline Surface (                                );
        inline Surface ( Surface const & src            );
        inline Surface ( Representation* representation );
        inline Surface ( AllM::Surface const & src      );

        inline int                   nvertices               (           ) const ;  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        inline int                   nfaces                  (           ) const ;  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        inline Vertex                vertices                ( size_t i  ) const ;
        inline Face                  faces                   ( size_t i  ) const ;
        inline FaceNormCacheEntry    faceNormCacheEntries    ( size_t i  ) const ;
        inline FaceNormDeferredEntry faceNormDeferredEntries ( size_t i  ) const ;
        inline double                radius                  (           ) const ;  //  radius (if status==MRIS_SPHERE)
        inline MRIS_Status           status                  (           ) const ;  //  type of surface (e.g. sphere, plane)
        inline MRIS_Status           origxyz_status          (           ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline int                   patch                   (           ) const ;  //  if a patch of the surface
        inline int                   noscale                 (           ) const ;  //  don't scale by surface area if true
    }; // Surface

    } // namespace TopologyM
