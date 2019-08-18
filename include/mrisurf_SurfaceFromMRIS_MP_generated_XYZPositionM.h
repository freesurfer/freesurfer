    namespace XYZPositionM {
    struct Face : public Repr_Elt {
        typedef XYZPositionM::Surface Surface;
        typedef XYZPositionM::Vertex  Vertex;
        inline Face                        (                                            );
        inline Face (                        Face const & src                           );
        inline Face (                        Representation* representation, size_t idx );
        inline Face (                        AllM::Face const & src                     );
        int fno     () const { return idx; }

        inline Vertex v           ( size_t i         ) const ;
        inline char   ripflag     (                  ) const ;
        
        inline void   set_ripflag (          char to       ) ;
    }; // Face

    struct Vertex : public Repr_Elt {
        typedef XYZPositionM::Surface Surface;
        typedef XYZPositionM::Face    Face;
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
        inline float cx            (                   ) const ;
        inline float cy            (                   ) const ;
        inline float cz            (                   ) const ;  //  coordinates in canonical coordinate system
        inline char  ripflag       (                   ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache
        inline void  which_coords  (int which, float *x, float *y, float *z) const ;
        
        inline void  set_x         (          float to       ) ;  //  current coordinates	
        inline void  set_y         (          float to       ) ;  //  use MRISsetXYZ() to set
        inline void  set_z         (          float to       ) ;
        inline void  set_cx        (          float to       ) ;
        inline void  set_cy        (          float to       ) ;
        inline void  set_cz        (          float to       ) ;  //  coordinates in canonical coordinate system
        inline void  set_ripflag   (           char to       ) ;  //  vertex no longer exists - placed last to load the next vertex into cache
    }; // Vertex

    struct MRIS_MP : public Repr_Elt {
        typedef XYZPositionM::Surface Surface;
        typedef XYZPositionM::Face    Face;
        typedef XYZPositionM::Vertex  Vertex;
        inline MRIS_MP (                                            );
        inline MRIS_MP ( MRIS_MP const & src                        );
        inline MRIS_MP ( Representation* representation, size_t idx );
        inline MRIS_MP ( AllM::MRIS_MP const & src                  );

    }; // MRIS_MP

    struct Surface : public Repr_Elt {
        typedef XYZPositionM::Face    Face;
        typedef XYZPositionM::Vertex  Vertex;
        inline Surface                                                (                                );
        inline Surface                                                ( Surface const & src            );
        inline Surface                                                ( Representation* representation );
        inline Surface                                                ( AllM::Surface const & src      );
        void freeDistsButNotOrig() { MRISfreeDistsButNotOrig(repr); }

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

    } // namespace XYZPositionM
