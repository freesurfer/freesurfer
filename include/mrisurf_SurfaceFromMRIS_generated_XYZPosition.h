    namespace XYZPosition {
    struct Face : public MRIS_Elt {
        inline Face (                                            );
        inline Face ( Face const & src                           );
        inline Face ( MRIS* mris, size_t idx                     );
        inline Face ( XYZPositionConsequencesM::Face const & src );
        inline Face ( XYZPositionConsequences::Face const & src  );
        inline Face ( DistortM::Face const & src                 );
        inline Face ( Distort::Face const & src                  );
        inline Face ( AnalysisM::Face const & src                );
        inline Face ( Analysis::Face const & src                 );
        inline Face ( AllM::Face const & src                     );

        inline vertices_per_face_t v        (   ) const ;
        inline char                ripflag  (   ) const ;
        inline char                oripflag (   ) const ;
        inline int                 marked   (   ) const ;
                   
        inline void set_ripflag  (  char to ) ;
        inline void set_oripflag (  char to ) ;
        inline void set_marked   (   int to ) ;
    };

    struct Vertex : public MRIS_Elt {
        inline Vertex (                                              );
        inline Vertex ( Vertex const & src                           );
        inline Vertex ( MRIS* mris, size_t idx                       );
        inline Vertex ( XYZPositionConsequencesM::Vertex const & src );
        inline Vertex ( XYZPositionConsequences::Vertex const & src  );
        inline Vertex ( DistortM::Vertex const & src                 );
        inline Vertex ( Distort::Vertex const & src                  );
        inline Vertex ( AnalysisM::Vertex const & src                );
        inline Vertex ( Analysis::Vertex const & src                 );
        inline Vertex ( AllM::Vertex const & src                     );

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
        inline float  cx            (           ) const ;                                                                                   
        inline float  cy            (           ) const ;                                                                                   
        inline float  cz            (           ) const ;  //  coordinates in canonical coordinate system                                   
        inline char   ripflag       (           ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache     
        // put the pointers before the ints, before the shorts, before uchars, to reduce size
        // the whole fits in much less than one cache line, so further ordering is no use
                   
        inline void set_f       ( size_t i,   Face to ) ;  // size() is num.    array[v->num] the fno's of the neighboring faces         
        inline void set_n       ( size_t i, size_t to ) ;  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        inline void set_e       ( size_t i,    int to )                  ;  //  edge state for neighboring vertices                      
        inline void set_cx      (            float to )                                                                                 ;
        inline void set_cy      (            float to )                                                                                 ;
        inline void set_cz      (            float to )                                 ;  //  coordinates in canonical coordinate system
        inline void set_ripflag (             char to )   ;  //  vertex no longer exists - placed last to load the next vertex into cache
    };

    struct Surface : public MRIS_Elt {
        inline Surface (                                               );
        inline Surface ( Surface const & src                           );
        inline Surface ( MRIS* mris, size_t idx                        );
        inline Surface ( XYZPositionConsequencesM::Surface const & src );
        inline Surface ( XYZPositionConsequences::Surface const & src  );
        inline Surface ( DistortM::Surface const & src                 );
        inline Surface ( Distort::Surface const & src                  );
        inline Surface ( AnalysisM::Surface const & src                );
        inline Surface ( Analysis::Surface const & src                 );
        inline Surface ( AllM::Surface const & src                     );

        inline MRIS_fname_t fname          (   ) const ;  //  file it was originally loaded from                                       
        inline MRIS_Status  status         (   ) const ;  //  type of surface (e.g. sphere, plane)                                     
        inline MRIS_Status  origxyz_status (   ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline int          patch          (   ) const ;  //  if a patch of the surface                                                
        inline p_void       vp             (   ) const ;  //  for misc. use                                                            
        inline float        alpha          (   ) const ;  //  rotation around z-axis                                                   
        inline float        beta           (   ) const ;  //  rotation around y-axis                                                   
        inline float        gamma          (   ) const ;  //  rotation around x-axis                                                   
        inline float        da             (   ) const ;                                                                               
        inline float        db             (   ) const ;                                                                               
        inline float        dg             (   ) const ;  //  old deltas                                                               
        inline int          type           (   ) const ;  //  what type of surface was this initially                                  
                   
        inline void set_vp    (  p_void to )                           ;  //  for misc. use
        inline void set_alpha (   float to )                  ;  //  rotation around z-axis
        inline void set_beta  (   float to )                  ;  //  rotation around y-axis
        inline void set_gamma (   float to )                  ;  //  rotation around x-axis
        inline void set_da    (   float to )                                              ;
        inline void set_db    (   float to )                                              ;
        inline void set_dg    (   float to )                              ;  //  old deltas
        inline void set_type  (     int to ) ;  //  what type of surface was this initially
    };

    } // namespace XYZPosition
