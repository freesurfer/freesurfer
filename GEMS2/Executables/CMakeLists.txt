# Find TetGen (needed for building atlases)
FIND_PATH( TETGEN_INCLUDE_DIR tetgen.h ../tetgen1.4.2
           DOC "Include directory of TetGen (i.e. containing tetgen.h)" )
FIND_LIBRARY( TETGEN_LIB tet ../tetgen1.4.2 )

SET( TETGEN_FOUND "NO" )
IF ( TETGEN_INCLUDE_DIR AND TETGEN_LIB )
  SET( TETGEN_FOUND "YES" )
ENDIF ( TETGEN_INCLUDE_DIR AND TETGEN_LIB )

IF ( TETGEN_FOUND )
  INCLUDE_DIRECTORIES( ${TETGEN_INCLUDE_DIR} )
ELSE ( TETGEN_FOUND )
  MESSAGE( FATAL_ERROR "Cannot build without TETGEN." )
ENDIF ( TETGEN_FOUND )


# Add an executable
ADD_EXECUTABLE( kvlAverageMeshes
                kvlAverageMeshes.cxx )
TARGET_LINK_LIBRARIES( kvlAverageMeshes
                       kvlGEMSCommon
                      )

# Add an executable
ADD_EXECUTABLE( kvlCompressImage
                kvlCompressImage.cxx )
TARGET_LINK_LIBRARIES( kvlCompressImage
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlBuildAtlasMesh
                kvlBuildAtlasMesh.cxx
                kvlAtlasMeshBuilder.cxx
                kvlAtlasMeshCollectionFastReferencePositionCost.cxx
                kvlAtlasMeshCollectionModelLikelihoodCalculator.cxx
                kvlMultiResolutionAtlasMesher.cxx
                kvlAtlasMeshVertexProcessor2.cxx
                kvlAtlasMeshCollectionPositionCostCalculator2.cxx
                )
TARGET_LINK_LIBRARIES( kvlBuildAtlasMesh
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlConvert
                kvlConvert.cxx )
TARGET_LINK_LIBRARIES( kvlConvert
                       kvlGEMSCommon
                       )

