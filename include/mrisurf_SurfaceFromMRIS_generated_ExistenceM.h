    namespace ExistenceM {
    struct Face : public MRIS_Elt {
        inline Face (                        );
        inline Face ( Face const & src       );
        inline Face ( MRIS* mris, size_t idx );
        inline Face ( AllM::Face const & src );

        inline char ripflag  (   ) const ;
        inline char oripflag (   ) const ;
        inline int  marked   (   ) const ;
                   
        inline void set_ripflag  (  char to ) ;
        inline void set_oripflag (  char to ) ;
        inline void set_marked   (   int to ) ;
    };

    struct Vertex : public MRIS_Elt {
        inline Vertex (                          );
        inline Vertex ( Vertex const & src       );
        inline Vertex ( MRIS* mris, size_t idx   );
        inline Vertex ( AllM::Vertex const & src );

        inline char ripflag (   ) const ;  //  vertex no longer exists - placed last to load the next vertex into cache
                   
        inline void set_ripflag (  char to ) ;  //  vertex no longer exists - placed last to load the next vertex into cache
    };

    struct Surface : public MRIS_Elt {
        inline Surface (                           );
        inline Surface ( Surface const & src       );
        inline Surface ( MRIS* mris, size_t idx    );
        inline Surface ( AllM::Surface const & src );

        inline MRIS_fname_t fname          (   ) const ;  //  file it was originally loaded from                                       
        inline MRIS_Status  status         (   ) const ;  //  type of surface (e.g. sphere, plane)                                     
        inline MRIS_Status  origxyz_status (   ) const ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline int          patch          (   ) const ;  //  if a patch of the surface                                                
                   
        inline void set_fname          (  MRIS_fname_t to )                                        ;  //  file it was originally loaded from
        inline void set_status         (   MRIS_Status to )                                      ;  //  type of surface (e.g. sphere, plane)
        inline void set_origxyz_status (   MRIS_Status to ) ;  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        inline void set_patch          (           int to )                                                 ;  //  if a patch of the surface
    };

    } // namespace ExistenceM
