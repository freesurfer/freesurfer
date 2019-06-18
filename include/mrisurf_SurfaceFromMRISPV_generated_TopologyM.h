    namespace TopologyM {
    struct Face : public Repr_Elt {
        inline Face                        (                                            );
        inline Face (                        Face const & src                           );
        inline Face (                        Representation* representation, size_t idx );
        inline Face (                        AllM::Face const & src                     );
        int fno     () const { return idx; }

        inline Vertex v       ( size_t i  ) const ;
        inline char   ripflag (           ) const ;
        
        inline void set_v       ( size_t i, Vertex to ) ;
        inline void set_ripflag (             char to ) ;
    };

    struct Vertex : public Repr_Elt {
        inline Vertex (                                                                   );
        inline Vertex (                        Vertex const & src                         );
        inline Vertex (                        Representation* representation, size_t idx );
        inline Vertex (                        AllM::Vertex const & src                   );
        int vno       () const { return idx; }

        // put the pointers before the ints, before the shorts, before uchars, to reduce size
        // the whole fits in much less than one cache line, so further ordering is no use
        inline Face  f            ( size_t i  ) const ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
        inline uchar num          (           ) const ;  //  number of neighboring faces                              
        inline char  ripflag      (           ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache
        inline void  which_coords (int which, float *x, float *y, float *z) const ;
        // put the pointers before the ints, before the shorts, before uchars, to reduce size
        // the whole fits in much less than one cache line, so further ordering is no use
        
        inline void set_f       ( size_t i,  Face to ) ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
        inline void set_num     (           uchar to ) ;  //  number of neighboring faces                              
        inline void set_ripflag (            char to ) ;  //  vertex no longer exists - placed last to load the next vertex into cache
    };

    struct Surface : public Repr_Elt {
        inline Surface (                                );
        inline Surface ( Surface const & src            );
        inline Surface ( Representation* representation );
        inline Surface ( AllM::Surface const & src      );

        // Fields being maintained by specialist functions
        inline int         nvertices      (           ) const ;  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        inline int         nfaces         (           ) const ;  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        inline Vertex      vertices       ( size_t i  ) const ;
        inline Face        faces          ( size_t i  ) const ;
        inline double      radius         (           ) const ;  //  radius (if status==MRIS_SPHERE)
        inline MRIS_Status status         (           ) const ;  //  type of surface (e.g. sphere, plane)
        inline MRIS_Status origxyz_status (           ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
    };

    } // namespace TopologyM
