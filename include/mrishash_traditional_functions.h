#ifndef MHT_ONLY_STATIC


// Add/remove the faces of which vertex vno is a part
//
MHT_VIRTUAL int  MHT_FUNCTION(addAllFaces)   (MHT_THIS_PARAMETER MHT_MRIS_PARAMETER int vno) MHT_ABSTRACT
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithFaces();
  return mht->addAllFaces(vno); }
#endif

MHT_VIRTUAL int  MHT_FUNCTION(removeAllFaces)(MHT_THIS_PARAMETER MHT_MRIS_PARAMETER int vno) MHT_ABSTRACT
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithFaces();
  return mht->removeAllFaces(vno); }
#endif



// Surface self-intersection (Uses MHT_FUNCTION() initialized with FACES)
//
MHT_VIRTUAL int MHT_FUNCTION(doesFaceIntersect)(MHT_THIS_PARAMETER MHT_MRIS_PARAMETER int fno) MHT_ABSTRACT
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithFaces();
  return mht->doesFaceIntersect(fno); }
#endif



MHT_VIRTUAL int MHT_FUNCTION(isVectorFilled)(MHT_CONST_THIS_PARAMETER  int vno, 
                                                float dx, float dy, float dz) MHT_CONST_THIS MHT_ABSTRACT
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ return mht->isVectorFilled(vno,dx,dy,dz); }
#endif


// Find nearest vertex/vertices (Uses MHT_FUNCTION() initialized with VERTICES)
//
MHT_VIRTUAL int MHT_FUNCTION(findClosestVertexGeneric)(MHT_THIS_PARAMETER
                                double probex, double probey, double probez,
                                double in_max_distance_mm, 
                                int in_max_halfmhts,
                                int *vtxnum, 
                                double *vtx_distance) MHT_ABSTRACT
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ return mht->findClosestVertexGeneric(
                                probex, probey, probez,
                                in_max_distance_mm, 
                                in_max_halfmhts,
                                vtxnum, 
                                vtx_distance); }
#endif


MHT_VIRTUAL int     MHT_FUNCTION(findClosestVertexNoXYZ)(MHT_THIS_PARAMETER
                               MHT_MRIS_PARAMETER  
                               float x, float y, float z, 
                               float *min_dist) MHT_ABSTRACT
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithVertices();
  return mht->findClosestVertexNoXYZ(x,y,z,min_dist); }
#endif


                             
MHT_VIRTUAL int     MHT_FUNCTION(findClosestSetVertexNo)(MHT_THIS_PARAMETER
                                MHT_MRIS_PARAMETER  
                                float x, float y, float z) MHT_ABSTRACT
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithVertices();
  return mht->findClosestSetVertexNo(x,y,z); }
#endif


                                            
MHT_VIRTUAL int MHT_FUNCTION(findVnoOfClosestVertexInTable)(MHT_THIS_PARAMETER
                                MHT_MRIS_PARAMETER
                                float x, float y, float z, int do_global_search) MHT_ABSTRACT
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithVertices();
  return mht->findVnoOfClosestVertexInTable(x,y,z,do_global_search); }
#endif



// utilities for finding closest face
//
MHT_VIRTUAL void MHT_FUNCTION(findClosestFaceNoGeneric)(MHT_THIS_PARAMETER
                              MHT_MRIS_PARAMETER 
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
                              double *pface_distance) MHT_ABSTRACT
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithFaces();    // seen to fail when Vertices
  mht->findClosestFaceNoGeneric(probex, probey, probez,
                              in_max_distance_mm,
                              in_max_mhts,
                              project_into_face,
                              pfno, 
                              pface_distance); }
#endif


#endif  // not MHT_STATIC_ONLY
                     

// The static functions that need to be given an MRIS or other Surface as a parameter
//
#ifndef MHT_ONLY_VIRTUAL

        
MHT_STATIC_MEMBER int MHT_FUNCTION(BruteForceClosestFace)(MHT_MRIS_PARAMETER  
                             float x, float y, float z, 
                             int which,                  // which surface within mris to search
                             float *dmin)
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ return MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::BruteForceClosestFace(mris,x,y,z,which,dmin); }
#endif


MHT_STATIC_MEMBER int MHT_FUNCTION(testIsMRISselfIntersecting)(MHT_MRIS_PARAMETER  float res)
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ return MRIS_HASH_TABLE::testIsMRISselfIntersecting(mris, res); }
#endif

// Construction and destruction
// Given that the mris is stored in the MHT_FUNCTION(), it is very unclear why the following functions require it to be passed in again!
//
MHT_STATIC_MEMBER MRIS_HASH_TABLE* MHT_FUNCTION(createFaceTable)(
    MHT_MRIS_PARAMETER_NOCOMMA)
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ return MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::createFaceTable(mris); }
#endif


MHT_STATIC_MEMBER MRIS_HASH_TABLE* MHT_FUNCTION(createFaceTable_Resolution)(
    MHT_MRIS_PARAMETER  
    int   which, 
    float res)
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ return MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::createFaceTable_Resolution(mris,which,res); }
#endif


MHT_STATIC_MEMBER MRIS_HASH_TABLE* MHT_FUNCTION(createVertexTable)(
    MHT_MRIS_PARAMETER  
    int which)
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ return MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::createVertexTable(mris,which); }
#endif

                                    
MHT_STATIC_MEMBER MRIS_HASH_TABLE* MHT_FUNCTION(createVertexTable_Resolution)(
    MHT_MRIS_PARAMETER 
    int which,
    float res)
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ return MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::createVertexTable_Resolution(mris,which,res); }
#endif

MHT_STATIC_MEMBER void MHT_FUNCTION(free)(MRIS_HASH_TABLE**mht)
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::free(mht); }
#endif

int MHT_FUNCTION(which)(MHT_THIS_PARAMETER_NOCOMMA)     // Whether uses the ORIGINAL, CANONICAL, CURRENT, ... ###xyz values
#ifndef MHT_TRADITIONAL_IMPL
    ;
#else
{ return mht->which(); }
#endif


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
