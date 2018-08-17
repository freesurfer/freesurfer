find_package(Matlab REQUIRED)

# build the mex file containing all the real functionality
matlab_add_mex(
  NAME kvlGEMSMatlab
  SRC kvlGEMSMatlab.cxx kvlMatlabRunnerArray.cxx kvlMatlabObjectArray.cxx
  LINK_TO ${ITK_LIBRARIES} kvlGEMSCommon
)

# copy some scripts
configure_file(kvlReadCompressionLookupTable.m
  ${GEMS_LIBRARY_PATH}/kvlReadCompressionLookupTable.m COPYONLY)
configure_file(kvlMergeAlphas.m
  ${GEMS_LIBRARY_PATH}/kvlMergeAlphas.m COPYONLY)
configure_file(showImage.m
  ${GEMS_LIBRARY_PATH}/showImage.m COPYONLY)
configure_file(kvlColorCodeProbabilityImages.m
  ${GEMS_LIBRARY_PATH}/kvlColorCodeProbabilityImages.m COPYONLY)
configure_file(mosaicImages.m
  ${GEMS_LIBRARY_PATH}/mosaicImages.m COPYONLY)
configure_file(kvlEvaluateMeshPositionInVectorFormat.m
  ${GEMS_LIBRARY_PATH}/kvlEvaluateMeshPositionInVectorFormat.m COPYONLY)
configure_file(backprojectKroneckerProductBasisFunctions.m
  ${GEMS_LIBRARY_PATH}/backprojectKroneckerProductBasisFunctions.m COPYONLY)
configure_file(computePrecisionOfKroneckerProductBasisFunctions.m
  ${GEMS_LIBRARY_PATH}/computePrecisionOfKroneckerProductBasisFunctions.m COPYONLY)
configure_file(projectKroneckerProductBasisFunctions.m
  ${GEMS_LIBRARY_PATH}/projectKroneckerProductBasisFunctions.m COPYONLY)
configure_file(kvlReadSharedGMMParameters.m
  ${GEMS_LIBRARY_PATH}/kvlReadSharedGMMParameters.m COPYONLY)
configure_file(kvlReadCompressionLookupTable.m
  ${GEMS_LIBRARY_PATH}/kvlReadCompressionLookupTable.m COPYONLY)
configure_file(kvlWriteCompressionLookupTable.m
  ${GEMS_LIBRARY_PATH}/kvlWriteCompressionLookupTable.m COPYONLY)
configure_file(kvlWriteSharedGMMParameters.m
  ${GEMS_LIBRARY_PATH}/kvlWriteSharedGMMParameters.m COPYONLY)
configure_file(kvlWarpMesh.m
  ${GEMS_LIBRARY_PATH}/kvlWarpMesh.m COPYONLY)

# for each function, also create a .m file calling the mex file with the correct arguments
set(runnerName "GetImageBuffer")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "SetImageBuffer")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "ReadCroppedImage")
set(runnerNumberOfOutputs "4")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "ReadMeshCollection")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetMesh")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "RasterizeAtlasMesh")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "SmoothImageBuffer")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetAlphasInMeshNodes")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "SetAlphasInMeshNodes")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "SmoothMeshCollection")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "ScaleMeshCollection")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "BiasFieldCorrectImage")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "WriteImage")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "CreateImage")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "DeformMesh")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetLevenbergMarquardtOptimizer")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "Clear")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "SetMaximumNumberOfThreads")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "ScaleMesh")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "SmoothMesh")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "ReadImage")
set(runnerNumberOfOutputs "2")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetTransformMatrix")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetMeshNodePositions")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "SetMeshCollectionPositions")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "WriteMeshCollection")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "TransformMeshCollection")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "SetKOfMeshCollection")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "CreateTransform")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "EvaluateMeshPosition")
set(runnerNumberOfOutputs "2")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "EvaluateMeshPositionWithEntropy")
set(runnerNumberOfOutputs "2")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "SetMeshNodePositions")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "ChangeK")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m 
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetImageFromOptimizer")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m 
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "CreateRGBImage")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m 
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "WriteRGBImage")
set(runnerNumberOfOutputs "0")
configure_file(kvlRunnerWrapperTemplate.m 
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetCroppedRegion")
set(runnerNumberOfOutputs "4")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetCostAndGradientCalculator")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetAverageAtlasMeshPositionCostAndGradientCalculator")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "GetOptimizer")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "StepOptimizer")
set(runnerNumberOfOutputs "2")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)

set(runnerName "CreateMeshCollection")
set(runnerNumberOfOutputs "1")
configure_file(kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m)
