    namespace ExistenceM {
    struct Face : public Repr_Elt {
        typedef ExistenceM::Surface Surface;
        typedef ExistenceM::Vertex  Vertex;
        inline Face                        (                                            );
        inline Face (                        Face const & src                           );
        inline Face (                        Representation* representation, size_t idx );
        inline Face (                        AllM::Face const & src                     );
        int fno     () const { return idx; }

        inline char ripflag     (          ) const ;
        
        inline void set_ripflag (  char to       ) ;
    }; // Face

    struct Vertex : public Repr_Elt {
        typedef ExistenceM::Surface Surface;
        typedef ExistenceM::Face    Face;
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
        typedef ExistenceM::Surface Surface;
        typedef ExistenceM::Face    Face;
        typedef ExistenceM::Vertex  Vertex;
        inline MRIS_MP (                                            );
        inline MRIS_MP ( MRIS_MP const & src                        );
        inline MRIS_MP ( Representation* representation, size_t idx );
        inline MRIS_MP ( AllM::MRIS_MP const & src                  );

    }; // MRIS_MP

    struct Surface : public Repr_Elt {
        typedef ExistenceM::Face    Face;
        typedef ExistenceM::Vertex  Vertex;
        inline Surface (                                );
        inline Surface ( Surface const & src            );
        inline Surface ( Representation* representation );
        inline Surface ( AllM::Surface const & src      );

        inline double      radius             (                 ) const ;  //  radius (if status==MRIS_SPHERE)
        inline MRIS_Status status             (                 ) const ;  //  type of surface (e.g. sphere, plane)
        inline MRIS_Status origxyz_status     (                 ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline int         patch              (                 ) const ;  //  if a patch of the surface
        inline int         noscale            (                 ) const ;  //  don't scale by surface area if true
        
        inline void        set_status         (  MRIS_Status to       ) ;  //  type of surface (e.g. sphere, plane)
        inline void        set_origxyz_status (  MRIS_Status to       ) ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline void        set_patch          (          int to       ) ;  //  if a patch of the surface
    }; // Surface

    } // namespace ExistenceM
