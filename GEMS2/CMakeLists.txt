CMAKE_MINIMUM_REQUIRED ( VERSION 3.5 )

PROJECT(GEMS)

# Output directories.
SET (LIBRARY_OUTPUT_PATH ${GEMS_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${GEMS_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
SET(GEMS_LIBRARY_PATH "${LIBRARY_OUTPUT_PATH}")
SET(GEMS_EXECUTABLE_PATH "${EXECUTABLE_OUTPUT_PATH}")

# 
INCLUDE_DIRECTORIES(${GEMS_SOURCE_DIR})


# Find ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR
           "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# IF ( NOT ITK_SOURCE_DIR )
#   FIND_PATH( ITK_SOURCE_DIR InsightLogo.gif
#              DOC "Source directory of ITK (needed for including some vxl header files" )
# ENDIF( NOT ITK_SOURCE_DIR )
# 
# IF ( ITK_SOURCE_DIR )
#   INCLUDE_DIRECTORIES( ${ITK_SOURCE_DIR}/Utilities/vxl/v3p/netlib )
# ELSE ( ITK_SOURCE_DIR )
#   MESSAGE( FATAL_ERROR "Cannot build without ITK_SOURCE_DIR." )
# ENDIF ( ITK_SOURCE_DIR )


OPTION( BUILD_SHARED_LIBS "Build GEMS with shared libraries." ON )

OPTION( BUILD_EXECUTABLES "Build command line executables" OFF )

OPTION( BUILD_TESTING "Build tests" OFF )

OPTION( BUILD_GUI "Build GUI components (requires FLTK and VTK)" OFF )

OPTION( BUILD_MATLAB "Build Matlab wrappers" OFF )



# Find TetGen
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


# Find Gmm++
FIND_PATH( GMM_INCLUDE_DIR config.h ../gmm-5.0
           DOC "Include directory of Gmm++ (i.e. containing config.h)" )

SET( GMM_FOUND "NO" )
IF ( GMM_INCLUDE_DIR )
  SET( GMM_FOUND "YES" )
ENDIF ( GMM_INCLUDE_DIR )

IF ( GMM_FOUND )
  INCLUDE_DIRECTORIES( ${GMM_INCLUDE_DIR}
                       ${GMM_INCLUDE_DIR}/include )
ELSE ( GMM_FOUND )
  MESSAGE( FATAL_ERROR "Cannot build without Gmm++." )
ENDIF ( GMM_FOUND )


# Add libary
ADD_LIBRARY( kvlGEMSCommon
             kvlAtlasMeshCollection.cxx
             gzstream.C
             kvlAtlasMeshCollectionValidator.cxx
             kvlAtlasParameterEstimator.cxx
             kvlAtlasMeshSmoother.cxx
             kvlCompressionLookupTable.cxx
             kvlAtlasMeshBuilder.cxx
             kvlAtlasMeshCollectionReferencePositionCost.cxx
             kvlAtlasMeshCollectionFastReferencePositionCost.cxx
             kvlAtlasMeshCollectionModelLikelihoodCalculator.cxx
             kvlAtlasMeshSegmenter.cxx
             kvlAtlasMeshSegmentationDriver.cxx
             kvlEMSegmenter.cxx
             kvlRegisterer.cxx
             kvlAtlasMeshToIntensityImagePartialVolumeGradientCalculator.cxx
             kvlMultiResolutionAtlasMesher.cxx
             kvlCroppedImageReader.cxx
             kvlAtlasMeshVertexProcessor2.cxx
             kvlAtlasMeshCollectionPositionCostCalculator2.cxx
             kvlImageSmoother.cxx
             kvlAtlasMeshDeformationLevenbergMarquardt.cxx
             kvlAtlasMeshToIntensityImageGradientCalculator.cxx
             #kvlAtlasMeshHamiltonianPositionSampler.cxx
             kvlAtlasMeshDeformationOptimizer.cxx
             kvlAtlasMeshDeformationLevenbergMarquardtOptimizer.cxx
             kvlAtlasMeshDeformationConjugateGradientOptimizer.cxx
             kvlAtlasMeshDeformationConjugateGradientOptimizerCPU.cxx
             kvlAtlasMeshAveragingConjugateGradientOptimizer.cxx
             kvlAtlasMeshRasterizorCPU.cxx
             kvlAtlasMeshAlphaDrawerCPU.cxx
             kvlAtlasMeshVisitCounterCPU.cxx
             kvlAtlasMeshToIntensityImageGradientCalculatorCPU.cxx
             itkMGHImageIO.cxx
             itkMGHImageIOFactory.cxx )
TARGET_LINK_LIBRARIES( kvlGEMSCommon
                       ${ITK_LIBRARIES}
                       ${TETGEN_LIB} 
                       ${ZLIB_LIBRARIES} 
                     )


# Build command line executables                       
IF(BUILD_EXECUTABLES)

  SUBDIRS( Executables )
  
ENDIF(BUILD_EXECUTABLES)

# Build tests                      
IF(BUILD_TESTING)

  SUBDIRS( Testing )
  
  SET(Boost_USE_STATIC_LIBS ON)
  FIND_PACKAGE( Boost 1.59.0 )
  IF(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  ENDIF(Boost_FOUND)

ENDIF(BUILD_TESTING)


# Build GUI
IF(BUILD_GUI)

  # Find FLTK
  FIND_PACKAGE(FLTK)
  IF(FLTK_FOUND)
    INCLUDE_DIRECTORIES(${FLTK_INCLUDE_DIR})
  ELSE(FLTK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build GUI components without FLTK.  Please set FLTK_DIR.")
  ENDIF(FLTK_FOUND)


  # Find VTK
  FIND_PACKAGE(VTK)
  IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
  ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR
           "Cannot build GUI components without VTK.  Please set VTK_DIR.")
  ENDIF(VTK_FOUND)

  SUBDIRS( GUI )

ENDIF(BUILD_GUI)


# Build Matlab wrappers
IF(BUILD_MATLAB)
  FIND_PACKAGE( Matlab REQUIRED )

  SUBDIRS( Matlab )

ENDIF(BUILD_MATLAB)

