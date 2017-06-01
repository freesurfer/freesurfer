# Build the mex file containing all the real functionality
MATLAB_ADD_MEX( NAME kvlGEMSMatlab
                SRC kvlGEMSMatlab.cxx kvlMatlabRunnerArray.cxx kvlMatlabObjectArray.cxx
                LINK_TO ${ITK_LIBRARIES} kvlGEMSCommon )


# Copy some scripts
CONFIGURE_FILE( kvlReadCompressionLookupTable.m
  ${GEMS_LIBRARY_PATH}/kvlReadCompressionLookupTable.m COPYONLY )
CONFIGURE_FILE( kvlReduceAlphas.m
  ${GEMS_LIBRARY_PATH}/kvlReduceAlphas.m COPYONLY )
CONFIGURE_FILE( showImage.m
  ${GEMS_LIBRARY_PATH}/showImage.m COPYONLY )
CONFIGURE_FILE( kvlColorCodeProbabilityImages.m
  ${GEMS_LIBRARY_PATH}/kvlColorCodeProbabilityImages.m COPYONLY )
CONFIGURE_FILE( smoothWithSeparableBasisFunctions.m
  ${GEMS_LIBRARY_PATH}/smoothWithSeparableBasisFunctions.m COPYONLY )
CONFIGURE_FILE( smoothWithSeparableBasisFunctionsWithWeights.m
  ${GEMS_LIBRARY_PATH}/smoothWithSeparableBasisFunctionsWithWeights.m COPYONLY )
CONFIGURE_FILE( expandBasisFunctions.m
  ${GEMS_LIBRARY_PATH}/expandBasisFunctions.m COPYONLY )
CONFIGURE_FILE( mosaicImages.m
  ${GEMS_LIBRARY_PATH}/mosaicImages.m COPYONLY )
CONFIGURE_FILE( kvlEvaluateMeshPositionInVectorFormat.m
  ${GEMS_LIBRARY_PATH}/kvlEvaluateMeshPositionInVectorFormat.m COPYONLY )
CONFIGURE_FILE( kvlMinFunc.m
  ${GEMS_LIBRARY_PATH}/kvlMinFunc.m COPYONLY )
CONFIGURE_FILE( kvlBackTrack.m
  ${GEMS_LIBRARY_PATH}/kvlBackTrack.m COPYONLY )
CONFIGURE_FILE( computeNumberOfSeparableBasisFunctions.m
  ${GEMS_LIBRARY_PATH}/computeNumberOfSeparableBasisFunctions.m COPYONLY )
CONFIGURE_FILE( computeRegularizers.m
  ${GEMS_LIBRARY_PATH}/computeRegularizers.m COPYONLY )
CONFIGURE_FILE( computeSeparableBasisFunctions.m
  ${GEMS_LIBRARY_PATH}/computeSeparableBasisFunctions.m COPYONLY )
CONFIGURE_FILE( constructAtWA.m
  ${GEMS_LIBRARY_PATH}/constructAtWA.m COPYONLY )
CONFIGURE_FILE( constructAtWr.m
  ${GEMS_LIBRARY_PATH}/constructAtWr.m COPYONLY )
CONFIGURE_FILE( expandSeparableBasisFunctions.m
  ${GEMS_LIBRARY_PATH}/expandSeparableBasisFunctions.m COPYONLY )
CONFIGURE_FILE( smoothWithSeparableBasisFunctionsWithWeightsMultiSpec.m
  ${GEMS_LIBRARY_PATH}/smoothWithSeparableBasisFunctionsWithWeightsMultiSpec.m COPYONLY )


# For each function, also create a .m file calling the mex file
# with the correct arguments. 
SET( runnerName "GetImageBuffer" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "SetImageBuffer" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "ReadCroppedImage" )
SET( runnerNumberOfOutputs "2" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "ReadMeshCollection" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "GetMesh" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "RasterizeAtlasMesh" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "SmoothImageBuffer" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "GetAlphasInMeshNodes" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "SetAlphasInMeshNodes" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "SmoothMeshCollection" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "ScaleMeshCollection" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "BiasFieldCorrectImage" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "WriteImage" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "CreateImage" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "DeformMesh" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "GetLevenbergMarquardtOptimizer" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "Clear" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "SetMaximumNumberOfThreads" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "ScaleMesh" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "SmoothMesh" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "ReadImage" )
SET( runnerNumberOfOutputs "2" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "GetTransformMatrix" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "GetMeshNodePositions" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "SetMeshCollectionPositions" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "WriteMeshCollection" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "TransformMeshCollection" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "SetKOfMeshCollection" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "CreateTransform" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "EvaluateMeshPosition" )
SET( runnerNumberOfOutputs "2" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "EvaluateMeshPositionWithEntropy" )
SET( runnerNumberOfOutputs "2" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "SetMeshNodePositions" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "ChangeK" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m 
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "GetImageFromOptimizer" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m 
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "CreateRGBImage" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m 
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "WriteRGBImage" )
SET( runnerNumberOfOutputs "0" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m 
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "GetCroppedRegion" )
SET( runnerNumberOfOutputs "4" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "GetCostAndGradientCalculator" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "GetOptimizer" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "StepOptimizer" )
SET( runnerNumberOfOutputs "2" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

SET( runnerName "CreateMeshCollection" )
SET( runnerNumberOfOutputs "1" )
CONFIGURE_FILE( kvlRunnerWrapperTemplate.m
  ${GEMS_LIBRARY_PATH}/kvl${runnerName}.m )

  