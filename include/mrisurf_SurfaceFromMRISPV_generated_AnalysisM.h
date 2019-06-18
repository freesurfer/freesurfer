    namespace AnalysisM {
    struct Face : public Repr_Elt {
        inline Face                        (                                            );
        inline Face (                        Face const & src                           );
        inline Face (                        Representation* representation, size_t idx );
        inline Face (                        AllM::Face const & src                     );
        int fno     () const { return idx; }

        inline Vertex                v       ( size_t i  ) const ;
        inline float                 area    (           ) const ;
        inline angles_per_triangle_t angle   (           ) const ;
        inline char                  ripflag (           ) const ;
        inline PDMATRIX              norm    (           ) const ;
        
        inline void set_ripflag (  char to ) ;
    };

    struct Vertex : public Repr_Elt {
        inline Vertex (                                                                   );
        inline Vertex (                        Vertex const & src                         );
        inline Vertex (                        Representation* representation, size_t idx );
        inline Vertex (                        AllM::Vertex const & src                   );
        int vno       () const { return idx; }

        // put the pointers before the ints, before the shorts, before uchars, to reduce size
        // the whole fits in much less than one cache line, so further ordering is no use
        inline Face  f             ( size_t i  ) const ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
        inline uchar num           (           ) const ;  //  number of neighboring faces                              
        // managed by MRISfreeDists[_orig] and MRISmakeDists[_orig]
        inline float dist          ( size_t i  ) const ;  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        inline int   dist_capacity (           ) const ;  //  -- should contain at least vtx_vtotal elements   
        // 
        inline float x             (           ) const ;  //  current coordinates	
        inline float y             (           ) const ;  //  use MRISsetXYZ() to set
        inline float z             (           ) const ;
        // 
        // 
        inline float nx            (           ) const ;
        inline float ny            (           ) const ;
        inline float nz            (           ) const ;  //  curr normal
        // 
        // 
        // 
        inline float area          (           ) const ;
        inline float origarea      (           ) const ;
        inline char  neg           (           ) const ;  //  1 if the normal vector is inverted 
        inline char  border        (           ) const ;  //  flag 
        inline char  ripflag       (           ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache
        inline void  which_coords  (int which, float *x, float *y, float *z) const ;
        // put the pointers before the ints, before the shorts, before uchars, to reduce size
        // the whole fits in much less than one cache line, so further ordering is no use
        
        inline void set_f        ( size_t i,  Face to ) ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
        // 
        // 
        inline void set_nx       (           float to ) ;
        inline void set_ny       (           float to ) ;
        inline void set_nz       (           float to ) ;  //  curr normal
        inline void set_origarea (           float to ) ;
        inline void set_neg      (            char to ) ;  //  1 if the normal vector is inverted 
        inline void set_border   (            char to ) ;  //  flag 
        inline void set_ripflag  (            char to ) ;  //  vertex no longer exists - placed last to load the next vertex into cache
    };

    struct Surface : public Repr_Elt {
        inline Surface (                                );
        inline Surface ( Surface const & src            );
        inline Surface ( Representation* representation );
        inline Surface ( AllM::Surface const & src      );

        // Fields being maintained by specialist functions
        inline int         nvertices       (           ) const ;  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        inline int         nfaces          (           ) const ;  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        inline Vertex      vertices        ( size_t i  ) const ;
        inline Face        faces           ( size_t i  ) const ;
        inline float       xctr            (           ) const ;
        inline float       yctr            (           ) const ;
        inline float       zctr            (           ) const ;
        inline float       xlo             (           ) const ;
        inline float       ylo             (           ) const ;
        inline float       zlo             (           ) const ;
        inline float       xhi             (           ) const ;
        inline float       yhi             (           ) const ;
        inline float       zhi             (           ) const ;
        inline float       total_area      (           ) const ;
        inline double      avg_vertex_area (           ) const ;
        inline double      avg_vertex_dist (           ) const ;  //  set by MRIScomputeAvgInterVertexDist
        inline double      std_vertex_dist (           ) const ;
        inline float       neg_area        (           ) const ;
        inline float       neg_orig_area   (           ) const ;  //  amount of original surface in folds
        inline double      radius          (           ) const ;  //  radius (if status==MRIS_SPHERE)
        inline MRIS_Status status          (           ) const ;  //  type of surface (e.g. sphere, plane)
        inline MRIS_Status origxyz_status  (           ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        
        inline void set_xctr            (   float to ) ;
        inline void set_yctr            (   float to ) ;
        inline void set_zctr            (   float to ) ;
        inline void set_xlo             (   float to ) ;
        inline void set_ylo             (   float to ) ;
        inline void set_zlo             (   float to ) ;
        inline void set_xhi             (   float to ) ;
        inline void set_yhi             (   float to ) ;
        inline void set_zhi             (   float to ) ;
        inline void set_total_area      (   float to ) ;
        inline void set_avg_vertex_area (  double to ) ;
    };

    } // namespace AnalysisM
