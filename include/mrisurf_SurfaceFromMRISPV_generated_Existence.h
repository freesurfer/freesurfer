    namespace Existence {
    struct Face : public Repr_Elt {
        inline Face                        (                                            );
        inline Face (                        Face const & src                           );
        inline Face (                        Representation* representation, size_t idx );
        inline Face (                        TopologyM::Face const & src                );
        inline Face (                        Topology::Face const & src                 );
        inline Face (                        XYZPositionM::Face const & src             );
        inline Face (                        XYZPosition::Face const & src              );
        inline Face (                        XYZPositionConsequencesM::Face const & src );
        inline Face (                        XYZPositionConsequences::Face const & src  );
        inline Face (                        DistortM::Face const & src                 );
        inline Face (                        Distort::Face const & src                  );
        inline Face (                        AnalysisM::Face const & src                );
        inline Face (                        Analysis::Face const & src                 );
        inline Face (                        AllM::Face const & src                     );
        int fno     () const { return idx; }

        inline char ripflag (   ) const ;
        
        inline void set_ripflag (  char to ) ;
    };

    struct Vertex : public Repr_Elt {
        inline Vertex (                                                                     );
        inline Vertex (                        Vertex const & src                           );
        inline Vertex (                        Representation* representation, size_t idx   );
        inline Vertex (                        TopologyM::Vertex const & src                );
        inline Vertex (                        Topology::Vertex const & src                 );
        inline Vertex (                        XYZPositionM::Vertex const & src             );
        inline Vertex (                        XYZPosition::Vertex const & src              );
        inline Vertex (                        XYZPositionConsequencesM::Vertex const & src );
        inline Vertex (                        XYZPositionConsequences::Vertex const & src  );
        inline Vertex (                        DistortM::Vertex const & src                 );
        inline Vertex (                        Distort::Vertex const & src                  );
        inline Vertex (                        AnalysisM::Vertex const & src                );
        inline Vertex (                        Analysis::Vertex const & src                 );
        inline Vertex (                        AllM::Vertex const & src                     );
        int vno       () const { return idx; }

        inline char ripflag      (    ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache
        inline void which_coords (int which, float *x, float *y, float *z) const ;
        
        inline void set_ripflag (  char to ) ;  //  vertex no longer exists - placed last to load the next vertex into cache
    };

    struct Surface : public Repr_Elt {
        inline Surface (                                               );
        inline Surface ( Surface const & src                           );
        inline Surface ( Representation* representation                );
        inline Surface ( TopologyM::Surface const & src                );
        inline Surface ( Topology::Surface const & src                 );
        inline Surface ( XYZPositionM::Surface const & src             );
        inline Surface ( XYZPosition::Surface const & src              );
        inline Surface ( XYZPositionConsequencesM::Surface const & src );
        inline Surface ( XYZPositionConsequences::Surface const & src  );
        inline Surface ( DistortM::Surface const & src                 );
        inline Surface ( Distort::Surface const & src                  );
        inline Surface ( AnalysisM::Surface const & src                );
        inline Surface ( Analysis::Surface const & src                 );
        inline Surface ( AllM::Surface const & src                     );

        inline double      radius         (   ) const ;  //  radius (if status==MRIS_SPHERE)
        inline MRIS_Status status         (   ) const ;  //  type of surface (e.g. sphere, plane)
        inline MRIS_Status origxyz_status (   ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
    };

    } // namespace Existence
