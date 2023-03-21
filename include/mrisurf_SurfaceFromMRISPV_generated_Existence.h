    namespace Existence {
    struct Face : public Repr_Elt {
        typedef Existence::Surface Surface;
        typedef Existence::Vertex  Vertex;
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

        inline char ripflag      (          ) const ;
        inline char oripflag     (          ) const ;
        inline int  marked       (          ) const ;
        
        inline void set_ripflag  (  char to       ) ;
        inline void set_oripflag (  char to       ) ;
        inline void set_marked   (   int to       ) ;
    }; // Face

    struct Vertex : public Repr_Elt {
        typedef Existence::Surface Surface;
        typedef Existence::Face    Face;
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

        inline char ripflag      (           ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache
        inline void which_coords (int which, float *x, float *y, float *z) const ;
        
        inline void set_ripflag  (   char to       ) ;  //  vertex no longer exists - placed last to load the next vertex into cache
    }; // Vertex

    struct Surface : public Repr_Elt {
        typedef Existence::Face    Face;
        typedef Existence::Vertex  Vertex;
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

        inline int                 initialized              (           ) const ;
        inline PLTA                lta                      (           ) const ;
        inline PMATRIX             SRASToTalSRAS_           (           ) const ;
        inline PMATRIX             TalSRASToSRAS_           (           ) const ;
        inline int                 free_transform           (           ) const ;
        inline double              radius                   (           ) const ;  //  radius (if status==MRIS_SPHERE)
        inline float               a                        (           ) const ;
        inline float               b                        (           ) const ;
        inline float               c                        (           ) const ;  //  ellipsoid parameters
        inline MRIS_fname_t        fname                    (           ) const ;  //  file it was originally loaded from
        inline MRIS_Status         status                   (           ) const ;  //  type of surface (e.g. sphere, plane)
        inline MRIS_Status         origxyz_status           (           ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline int                 patch                    (           ) const ;  //  if a patch of the surface
        inline int                 max_vertices             (           ) const ;  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        inline int                 max_faces                (           ) const ;  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        inline MRIS_subject_name_t subject_name             (           ) const ;  //  name of the subject
        inline float               canon_area               (           ) const ;
        inline int                 noscale                  (           ) const ;  //  don't scale by surface area if true
        inline float               dx2                      ( size_t i  ) const ;  //  an extra set of gradient (not always alloced)
        inline float               dy2                      ( size_t i  ) const ;
        inline float               dz2                      ( size_t i  ) const ;
        inline PCOLOR_TABLE        ct                       (           ) const ;
        inline int                 orig_xyzspace            (           ) const ;  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        inline int                 useRealRAS               (           ) const ;  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        inline VOL_GEOM            vg                       (           ) const ;  //  volume info from which this surface is created. valid iff vg.valid = 1
        inline MRIS_cmdlines_t     cmdlines                 (           ) const ;
        inline int                 ncmds                    (           ) const ;
        inline float               group_avg_surface_area   (           ) const ;  //  average of total surface area for group
        inline int                 group_avg_vtxarea_loaded (           ) const ;  //  average vertex area for group at each vertex
        inline int                 triangle_links_removed   (           ) const ;  //  for quad surfaces
        inline p_void              user_parms               (           ) const ;  //  for whatever the user wants to hang here
        inline PMATRIX             m_sras2vox               (           ) const ;  //  for converting surface ras to voxel
        inline PMRI                mri_sras2vox             (           ) const ;  //  volume that the above matrix is for
        inline p_void              mht                      (           ) const ;
        inline p_void              temps                    (           ) const ;
    }; // Surface

    } // namespace Existence
