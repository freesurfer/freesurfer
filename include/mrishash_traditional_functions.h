#ifndef MHT_ONLY_STATIC
    
// The specifications of the virtual functions.
// In the traditional c-style code, they are called passing the MHT as the first argument,
// and they often passed the MRIS as the second parameter, even though that should have been captured in the MHT
//
// They should be called c++-style,   p->foo(...) rather than MHTfoo(p,...)

// Add/remove the faces of which vertex vno is a part
//
MHT_VIRTUAL int  MHT_FUNCTION(addAllFaces)                  (MHT_THIS_PARAMETER MHT_MRIS_PARAMETER  int vno) MHT_ABSTRACT;
MHT_VIRTUAL int  MHT_FUNCTION(removeAllFaces)               (MHT_THIS_PARAMETER MHT_MRIS_PARAMETER  int vno) MHT_ABSTRACT;

// Surface self-intersection (Uses MHT initialized with FACES)
//
MHT_VIRTUAL int MHT_FUNCTION(doesFaceIntersect)             (MHT_THIS_PARAMETER MHT_MRIS_PARAMETER  int fno                                 )                MHT_ABSTRACT;
MHT_VIRTUAL int MHT_FUNCTION(isVectorFilled)                (MHT_CONST_THIS_PARAMETER               int vno, float dx, float dy, float dz   ) MHT_CONST_THIS MHT_ABSTRACT;

// Find nearest vertex/vertices (Uses MHT initialized with VERTICES)
//
MHT_VIRTUAL int MHT_FUNCTION(findClosestVertexNoXYZ)        (MHT_THIS_PARAMETER MHT_MRIS_PARAMETER float x, float y, float z, float *min_dist       ) MHT_ABSTRACT;
MHT_VIRTUAL int MHT_FUNCTION(findClosestSetVertexNo)        (MHT_THIS_PARAMETER MHT_MRIS_PARAMETER float x, float y, float z                        ) MHT_ABSTRACT;
MHT_VIRTUAL int MHT_FUNCTION(findVnoOfClosestVertexInTable) (MHT_THIS_PARAMETER MHT_MRIS_PARAMETER float x, float y, float z, int do_global_search  ) MHT_ABSTRACT;
MHT_VIRTUAL int MHT_FUNCTION(findClosestVertexGeneric)      (MHT_THIS_PARAMETER                    double probex, double probey, double probez, 
                                                                                                   double in_max_distance_mm, int in_max_halfmhts, 
                                                                                                   int *vtxnum,  double *vtx_distance               ) MHT_ABSTRACT;

// Find closest face
//
MHT_VIRTUAL void MHT_FUNCTION(findClosestFaceNoGeneric)(MHT_THIS_PARAMETER MHT_MRIS_PARAMETER 
                              //---------- inputs --------------
                              double probex, double probey, double probez,
                              // How far to search: set one or both
                              double in_max_distance_mm, /* Use large number 
                                                            to ignore */
                              int    in_max_mhts,  /* Use -1 to ignore */
                              // only faces that projection is interior to (Use -1 to ignore )
                              int    project_into_face, 
                              //---------- outputs -------------
                              int *pfno, 
                              double *pface_distance) MHT_ABSTRACT;

#endif

#ifndef MHT_ONLY_VIRTUAL

// The static functions that need to be given an MRIS or other Surface as a parameter
//
MHT_STATIC_MEMBER int MHT_FUNCTION(BruteForceClosestFace)(MHT_MRIS_PARAMETER  
                             float x, float y, float z, 
                             int which,                  // which surface within mris to search
                             float *dmin);

MHT_STATIC_MEMBER int MHT_FUNCTION(testIsMRISselfIntersecting)(MHT_MRIS_PARAMETER  float res);


#endif  // non virtual






#undef MHT_MRIS_PARAMETER
#undef MHT_THIS_PARAMETER
#undef MHT_MRIS_PARAMETER_NOCOMMA
#undef MHT_THIS_PARAMETER_NOCOMMA
#undef MHT_CONST_THIS_PARAMETER
#undef MHT_CONST_THIS
#undef MHT_FUNCTION
#undef MHT_STATIC_MEMBER
#undef MHT_VIRTUAL
#undef MHT_ABSTRACT
