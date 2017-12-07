CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

PROJECT(GEMS)

INCLUDE_DIRECTORIES(${GEMS_SOURCE_DIR})

# Find ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR
           "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Enable C++11 syntax
set(CMAKE_CXX_STANDARD 11)


add_subdirectory(pybind11)
pybind11_add_module(GEMS2Python
        GEMS2Python.cxx
        pyKvlCalculator.cxx
        pyKvlImage.cxx
        pyKvlMesh.cxx
        pyKvlTransform.cxx
        pyKvlOptimizer.cxx
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
        kvlAverageAtlasMeshPositionCostAndGradientCalculator.cxx
        kvlHistogrammer.cxx
        kvlAtlasMeshSummaryDrawer.cxx
        itkMGHImageIO.cxx
        itkMGHImageIOFactory.cxx
        )

TARGET_LINK_LIBRARIES(GEMS2Python
        PRIVATE ${ITK_LIBRARIES}
        PRIVATE ${ZLIB_LIBRARIES}
        )

