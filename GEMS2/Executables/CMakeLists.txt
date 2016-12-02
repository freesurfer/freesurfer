# Add an executable
ADD_EXECUTABLE( kvlAverageMeshes
                kvlAverageMeshes.cxx )
TARGET_LINK_LIBRARIES( kvlAverageMeshes
                       kvlGEMSCommon
                      )

# Add an executable
ADD_EXECUTABLE( kvlSmoothMesh
                kvlSmoothMesh.cxx )
TARGET_LINK_LIBRARIES( kvlSmoothMesh
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlCrop
                kvlCrop.cxx )
TARGET_LINK_LIBRARIES( kvlCrop
                       ${ITK_LIBRARIES} )

# Add an executable
ADD_EXECUTABLE( kvlChangeKOfMeshCollection
                kvlChangeKOfMeshCollection.cxx )
TARGET_LINK_LIBRARIES( kvlChangeKOfMeshCollection
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlChangePointParametersOfMeshCollection
                kvlChangePointParametersOfMeshCollection.cxx )
TARGET_LINK_LIBRARIES( kvlChangePointParametersOfMeshCollection
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlMapDeformationField
                kvlMapDeformationField.cxx )
TARGET_LINK_LIBRARIES( kvlMapDeformationField
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlUpsampleMesh
                kvlUpsampleMesh.cxx )
TARGET_LINK_LIBRARIES( kvlUpsampleMesh
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlSetReferencePosition
                kvlSetReferencePosition.cxx )
TARGET_LINK_LIBRARIES( kvlSetReferencePosition
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlCompressImage
                kvlCompressImage.cxx )
TARGET_LINK_LIBRARIES( kvlCompressImage
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlBuildAtlasMeshWithoutGUI
                kvlBuildAtlasMeshWithoutGUI.cxx )
TARGET_LINK_LIBRARIES( kvlBuildAtlasMeshWithoutGUI
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlCropImage
                kvlCropImage.cxx )
TARGET_LINK_LIBRARIES( kvlCropImage
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlCombineManualSegmentationProtocols
                kvlCombineManualSegmentationProtocols.cxx )
TARGET_LINK_LIBRARIES( kvlCombineManualSegmentationProtocols
                       ${ITK_LIBRARIES} )

# Add an executable
ADD_EXECUTABLE( kvlReduceNumberOfIntensityLevels
                kvlReduceNumberOfIntensityLevels.cxx )
TARGET_LINK_LIBRARIES( kvlReduceNumberOfIntensityLevels
                       ${ITK_LIBRARIES} )

# Add an executable
ADD_EXECUTABLE( kvlAugmentManualSegmentationWithAutomatedSegmentation
                kvlAugmentManualSegmentationWithAutomatedSegmentation.cxx )
TARGET_LINK_LIBRARIES( kvlAugmentManualSegmentationWithAutomatedSegmentation
                       ${ITK_LIBRARIES} )

# Add an executable
ADD_EXECUTABLE( kvlScaleMeshCollection
                kvlScaleMeshCollection.cxx )
TARGET_LINK_LIBRARIES( kvlScaleMeshCollection
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlMapLabelsOfAtlasMesh
                kvlMapLabelsOfAtlasMesh.cxx )
TARGET_LINK_LIBRARIES( kvlMapLabelsOfAtlasMesh
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlSegmentWithoutGUI
                kvlSegmentWithoutGUI.cxx )
TARGET_LINK_LIBRARIES( kvlSegmentWithoutGUI
                       kvlGEMSCommon
                       )


# Add an executable
ADD_EXECUTABLE( kvlEvaluateQuantitatively
                kvlEvaluateQuantitatively.cxx )
TARGET_LINK_LIBRARIES( kvlEvaluateQuantitatively
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlEvaluateQuantitativelyWithHausdorffDistance
                kvlEvaluateQuantitativelyWithHausdorffDistance.cxx )
TARGET_LINK_LIBRARIES( kvlEvaluateQuantitativelyWithHausdorffDistance
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlCorrectLeftSideLabelsToRightSideLabels
                kvlCorrectLeftSideLabelsToRightSideLabels.cxx )
TARGET_LINK_LIBRARIES( kvlCorrectLeftSideLabelsToRightSideLabels
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlOverwriteEmptySlices
                kvlOverwriteEmptySlices.cxx )
TARGET_LINK_LIBRARIES( kvlOverwriteEmptySlices
                       kvlGEMSCommon
                       )


# Add an executable
ADD_EXECUTABLE( kvlOptimizePositionsWithICM
                kvlOptimizePositionsWithICM.cxx )
TARGET_LINK_LIBRARIES( kvlOptimizePositionsWithICM
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlCalculatePositionCost
                kvlCalculatePositionCost.cxx )
TARGET_LINK_LIBRARIES( kvlCalculatePositionCost
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlCalculateAtlasMeshCollectionTotalCost
                kvlCalculateAtlasMeshCollectionTotalCost.cxx )
TARGET_LINK_LIBRARIES( kvlCalculateAtlasMeshCollectionTotalCost
                       kvlGEMSCommon
                       )


# Add an executable
ADD_EXECUTABLE( kvlGenerateWeightedColorChannels
                kvlGenerateWeightedColorChannels.cxx )
TARGET_LINK_LIBRARIES( kvlGenerateWeightedColorChannels
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlColorCodeProbabilityImages
                kvlColorCodeProbabilityImages.cxx )
TARGET_LINK_LIBRARIES( kvlColorCodeProbabilityImages
                       kvlGEMSCommon
                       )

# Add an executable
# ADD_EXECUTABLE( kvlSamplePositionsFromMeshCollection
#                 kvlSamplePositionsFromMeshCollection.cxx )
# TARGET_LINK_LIBRARIES( kvlSamplePositionsFromMeshCollection
#                        kvlGEMSCommon
#                        ${ITK_LIBRARIES} )

# Add an executable
ADD_EXECUTABLE( kvlAddMeshToMeshCollection
                kvlAddMeshToMeshCollection.cxx )
TARGET_LINK_LIBRARIES( kvlAddMeshToMeshCollection
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlMergeMeshCollections
                kvlMergeMeshCollections.cxx )
TARGET_LINK_LIBRARIES( kvlMergeMeshCollections
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlAddCSFToBucknerData
                kvlAddCSFToBucknerData.cxx )
TARGET_LINK_LIBRARIES( kvlAddCSFToBucknerData
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlAddBorder
                kvlAddBorder.cxx )
TARGET_LINK_LIBRARIES( kvlAddBorder
                       ${ITK_LIBRARIES} )

# Add an executable
ADD_EXECUTABLE( kvlQuantifyPosteriorProbabilityImages
                kvlQuantifyPosteriorProbabilityImages.cxx )
TARGET_LINK_LIBRARIES( kvlQuantifyPosteriorProbabilityImages
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlEMSegment
                kvlEMSegment.cxx )
TARGET_LINK_LIBRARIES( kvlEMSegment
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlAddProbabilityImages
                kvlAddProbabilityImages.cxx )
TARGET_LINK_LIBRARIES( kvlAddProbabilityImages
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlMathematicalMorphology
                kvlMathematicalMorphology.cxx )
TARGET_LINK_LIBRARIES( kvlMathematicalMorphology
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlConvert
                kvlConvert.cxx )
TARGET_LINK_LIBRARIES( kvlConvert
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlRegister
                kvlRegister.cxx )
TARGET_LINK_LIBRARIES( kvlRegister
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlResample
                kvlResample.cxx )
TARGET_LINK_LIBRARIES( kvlResample
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlBinaryThresholdImage
                kvlBinaryThresholdImage.cxx )
TARGET_LINK_LIBRARIES( kvlBinaryThresholdImage
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlUpsample
                kvlUpsample.cxx )
TARGET_LINK_LIBRARIES( kvlUpsample
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlRasterizeAtlasMesh
                kvlRasterizeAtlasMesh.cxx )
TARGET_LINK_LIBRARIES( kvlRasterizeAtlasMesh
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlAutoCrop
                kvlAutoCrop.cxx )
TARGET_LINK_LIBRARIES( kvlAutoCrop
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlMaskImage
                kvlMaskImage.cxx )
TARGET_LINK_LIBRARIES( kvlMaskImage
                       kvlGEMSCommon
                       )

# Add an executable
ADD_EXECUTABLE( kvlApplyTransform
                kvlApplyTransform.cxx )
TARGET_LINK_LIBRARIES( kvlApplyTransform
                       kvlGEMSCommon
                       )
                       
                  
