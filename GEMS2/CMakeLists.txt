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
INCLUDE_DIRECTORIES(interfaces)


# Find ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR
           "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

OPTION( BUILD_SHARED_LIBS "Build GEMS with shared libraries." ON )

OPTION( BUILD_EXECUTABLES "Build command line executables" OFF )

OPTION( BUILD_TESTING "Build tests" OFF )

OPTION( BUILD_GUI "Build GUI components (requires FLTK and VTK)" OFF )

OPTION( BUILD_MATLAB "Build Matlab wrappers" OFF )

#
IF( BUILD_EXECUTABLES OR BUILD_GUI )
  #
  ADD_DEFINITIONS( -DUSE_DYNAMIC_MESH )
  
  # 
  MESSAGE( WARNING 
           "Since you're building executables to compute meshes, ITK dynamic meshes will be used (internally using std::map instead of std::vector to store lists of points, cells, etc). Everything will still work, but anything that only uses existing meshes will be slower than it could be (up to 100% slower)" )
           
ENDIF( BUILD_EXECUTABLES OR BUILD_GUI )




# Find CUDA
FIND_PACKAGE( CUDA )
IF( CUDA_FOUND )
  # Enable C++11 syntax
  SET(CMAKE_CXX_FLAGS "-g -std=c++11 ${CMAKE_CXX_FLAGS}")

  INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
  SET(CMAKE_CXX_FLAGS "-DCUDA_FOUND ${CMAKE_CXX_FLAGS}")

  INCLUDE_DIRECTORIES(cuda)
  SUBDIRS(cuda)
ENDIF( CUDA_FOUND )


# Add libary
ADD_LIBRARY( kvlGEMSCommon
             kvlAtlasMeshCollection.cxx
             gzstream.C
             kvlAtlasMeshCollectionValidator.cxx
             kvlAtlasParameterEstimator.cxx
             kvlAtlasMeshSmoother.cxx
             kvlCompressionLookupTable.cxx
             kvlCroppedImageReader.cxx
             #kvlAtlasMeshAveragingConjugateGradientOptimizer.cxx
             kvlAtlasMeshRasterizor.cxx
             kvlAtlasMeshAlphaDrawer.cxx
             kvlAtlasMeshMultiAlphaDrawer.cxx
             kvlAtlasMeshVisitCounter.cxx
             kvlAtlasMeshStatisticsCollector.cxx
             kvlAtlasMeshLabelImageStatisticsCollector.cxx
             kvlAtlasMeshProbabilityImageStatisticsCollector.cxx
             kvlAtlasMeshDeformationOptimizer.cxx
             kvlAtlasMeshDeformationFixedStepGradientDescentOptimizer.cxx
             kvlAtlasMeshDeformationGradientDescentOptimizer.cxx
             kvlAtlasMeshDeformationConjugateGradientOptimizer.cxx
             kvlAtlasMeshDeformationLBFGSOptimizer.cxx
             kvlAtlasMeshPositionCostAndGradientCalculator.cxx
             kvlAtlasMeshToLabelImageCostAndGradientCalculator.cxx
             kvlAtlasMeshToIntensityImageCostAndGradientCalculator.cxx
             kvlConditionalGaussianEntropyCostAndGradientCalculator.cxx
             kvlMutualInformationCostAndGradientCalculator.cxx
             kvlAtlasMeshToPointSetCostAndGradientCalculator.cxx
             kvlHistogrammer.cxx
             kvlAtlasMeshSummaryDrawer.cxx
             itkMGHImageIO.cxx
             itkMGHImageIOFactory.cxx )
TARGET_LINK_LIBRARIES( kvlGEMSCommon
                       ${ITK_LIBRARIES}
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

