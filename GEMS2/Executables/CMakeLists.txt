OPTION( MAKE_SPARSE_INITIAL_MESHES "Make sparse initial meshes" OFF )

#
IF( MAKE_SPARSE_INITIAL_MESHES )
  #
  ADD_DEFINITIONS( -DUSE_TETGEN )
  
  # 
  MESSAGE( WARNING 
           "Since you're opting to make sparse (i.e., irregular) initial meshes, the TetGen library will be included, which has a restrictive license." )
           
  # Find TetGen 
  FIND_PATH( TETGEN_INCLUDE_DIR tetgen.h ../../tetgen1.4.2
            DOC "Include directory of TetGen (i.e. containing tetgen.h)" )
  FIND_LIBRARY( TETGEN_LIB tet ../../tetgen1.4.2 )

  SET( TETGEN_FOUND "NO" )
  IF ( TETGEN_INCLUDE_DIR AND TETGEN_LIB )
    SET( TETGEN_FOUND "YES" )
  ENDIF ( TETGEN_INCLUDE_DIR AND TETGEN_LIB )

  IF ( TETGEN_FOUND )
    INCLUDE_DIRECTORIES( ${TETGEN_INCLUDE_DIR} )
  ELSE ( TETGEN_FOUND )
    MESSAGE( FATAL_ERROR "Cannot build without TETGEN." )
  ENDIF ( TETGEN_FOUND )           
           
ENDIF( MAKE_SPARSE_INITIAL_MESHES )


# # Add an executable
# ADD_EXECUTABLE( kvlAverageMeshes
#                 kvlAverageMeshes.cxx )
# TARGET_LINK_LIBRARIES( kvlAverageMeshes
#                        kvlGEMSCommon
#                       )

# Add an executable
ADD_EXECUTABLE( kvlBuildAtlasMesh
                kvlBuildAtlasMesh.cxx
                kvlAtlasMeshBuilder.cxx
                kvlAtlasMeshCollectionFastReferencePositionCost.cxx
                kvlAtlasMeshCollectionModelLikelihoodCalculator.cxx
                kvlMultiResolutionAtlasMesher.cxx
                )
TARGET_LINK_LIBRARIES( kvlBuildAtlasMesh
                       kvlGEMSCommon
                       ${TETGEN_LIB} 
                       )

