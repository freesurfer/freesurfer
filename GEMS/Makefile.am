if HAVE_ITK_LIBS

SUBDIRS=scripts data

ITK_THINGS= $(ITK_CFLAGS) -DHAVE_ITK_LIBS

AM_CXXFLAGS=-I$(top_srcdir)/include	\
	-fno-strict-aliasing \
	-Wno-uninitialized \
	-Wno-unused \
	-Wno-sign-compare \
	-Wno-return-type \
	-Wno-deprecated \
	-DTETLIBRARY \
	-DSELF_CHECK \
	$(GL_CFLAGS) \
	$(VXL_CFLAGS) \
	-I$(top_srcdir)/GEMS \
	-I$(top_srcdir)/tetgen \
	-I$(ITK_DIR)/Utilities/zlib \
	-I$(ITK_DIR)/Utilities/vxl/v3p/netlib \
	-I$(top_srcdir)/itkio \
	-I$(top_srcdir)/itkutils \
	$(ITK_THINGS) \
	$(FLTK_CXXFLAGS)\
	$(VTK_CFLAGS)

AM_LDFLAGS=

noinst_LIBRARIES=libkvlAtlas3DCommon.a

libkvlAtlas3DCommon_a_SOURCES= \
	kvlAtlasMeshCollection.cxx \
	gzstream.cxx \
	kvlAtlasMeshCollectionValidator.cxx \
	kvlAtlasParameterEstimator.cxx \
	kvlAtlasMeshSmoother.cxx \
	kvlCompressionLookupTable.cxx \
	kvlAtlasMeshBuilder.cxx \
	kvlAtlasMeshCollectionReferencePositionCost.cxx \
	kvlAtlasMeshCollectionFastReferencePositionCost.cxx \
	kvlAtlasMeshCollectionModelLikelihoodCalculator.cxx \
	kvlAtlasMeshSegmenter.cxx \
	kvlAtlasMeshSegmentationDriver.cxx \
	kvlEMSegmenter.cxx \
	kvlRegisterer.cxx \
	kvlAtlasMeshToIntensityImagePartialVolumeGradientCalculator.cxx \
	kvlMultiResolutionAtlasMesher.cxx \
	kvlCroppedImageReader.cxx \
	kvlAtlasMeshVertexProcessor2.cxx \
	kvlAtlasMeshCollectionPositionCostCalculator2.cxx \
	kvlImageSmoother.cxx \
	kvlAtlasMeshDeformationLevenbergMarquardt.cxx \
	kvlAtlasMeshHamiltonianPositionSampler.cxx

am_libkvlAtlas3DCommon_a_OBJECTS=$(libkvlAtlas3DCommon_a_SOURCES:%.cxx=%.o)
libkvlAtlas3DCommon_a_LIBADD= \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a	 \
	$(LIBS_ITK) \
	$(top_builddir)/tetgen/libtetgen.a \
	$(addprefix $(top_builddir)/, $(LIBS_MGH))

bin_PROGRAMS=\
	kvlSmoothMesh \
	kvlCrop \
	kvlChangeKOfMeshCollection \
	kvlChangePointParametersOfMeshCollection \
	kvlMapDeformationField \
	kvlUpsampleMesh \
	kvlSetReferencePosition \
	kvlCompressImage \
	kvlCropImage \
	kvlCombineManualSegmentationProtocols \
	kvlReduceNumberOfIntensityLevels \
	kvlAugmentManualSegmentationWithAutomatedSegmentation \
	kvlScaleMeshCollection \
	kvlMapLabelsOfAtlasMesh \
	kvlSegmentWithoutGUI \
	kvlEvaluateQuantitatively \
	kvlEvaluateQuantitativelyWithHausdorffDistance \
	kvlCorrectLeftSideLabelsToRightSideLabels \
	kvlOverwriteEmptySlices \
	kvlOptimizePositionsWithICM \
	kvlCalculatePositionCost \
	kvlCalculateAtlasMeshCollectionTotalCost \
	kvlColorCodeProbabilityImages \
	kvlSamplePositionsFromMeshCollection \
	kvlAddMeshToMeshCollection \
	kvlMergeMeshCollections \
	kvlAddBorder \
	kvlQuantifyPosteriorProbabilityImages \
	kvlEMSegment \
	kvlAddProbabilityImages \
	kvlMathematicalMorphology \
	kvlConvert \
	kvlRegister \
	kvlResample \
	kvlBinaryThresholdImage \
	kvlUpsample \
	kvlRasterizeAtlasMesh \
	kvlAutoCrop \
	kvlMaskImage \
	kvlApplyTransform

#kvlAddCSFToBucknerData
#kvlBuildAtlasMeshWithoutGUI
#kvlAtlasMeshRasterizorTest
#kvlAtlasMeshCollectionTest
#kvlLabelImageToCompactMeshRepresentation
#kvlValidateMeshCollection
#kvlGenerateToyExampleData
#kvlConstructSparseMeshCollection
#kvlConstructMultiResolutionMeshCollection
#kvlTestImageSmoother
#kvlGenerateWeightedColorChannels


# kvlAtlasMeshRasterizorTest_SOURCES=kvlAtlasMeshRasterizorTest.cxx
# kvlAtlasMeshRasterizorTest_LDADD=$(LIBS_ITK) \
#					 libkvlAtlas3DCommon.a \
#					 $(top_builddir)/itkio/libitkio.a \
#					 $(top_builddir)/itkutils/libitkutils.a \
#					 -lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
#					 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlAtlasMeshRasterizorTest_LDFLAGS= $(OS_LDFLAGS)

# kvlAtlasMeshCollectionTest_SOURCES=kvlAtlasMeshCollectionTest.cxx
# kvlAtlasMeshCollectionTest_LDADD=$(LIBS_ITK) \
#					 libkvlAtlas3DCommon.a \
#					 $(top_builddir)/itkio/libitkio.a \
#					 $(top_builddir)/itkutils/libitkutils.a \
#					 -lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
#					 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlAtlasMeshCollectionTest_LDFLAGS= $(OS_LDFLAGS)

# kvlLabelImageToCompactMeshRepresentation_SOURCES=kvlLabelImageToCompactMeshRepresentation.cxx
# kvlLabelImageToCompactMeshRepresentation_LDADD=$(LIBS_ITK) \
#					 libkvlAtlas3DCommon.a \
#					 $(top_builddir)/itkio/libitkio.a \
#					 $(top_builddir)/itkutils/libitkutils.a \
#					 -lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
#					 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlLabelImageToCompactMeshRepresentation_LDFLAGS= $(OS_LDFLAGS)

# kvlValidateMeshCollection_SOURCES=kvlValidateMeshCollection.cxx
# kvlValidateMeshCollection_LDADD=$(LIBS_ITK) \
#					 libkvlAtlas3DCommon.a \
#					 $(top_builddir)/itkio/libitkio.a \
#					 $(top_builddir)/itkutils/libitkutils.a \
#					 -lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
#					 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlValidateMeshCollection_LDFLAGS= $(OS_LDFLAGS)

# kvlGenerateToyExampleData_SOURCES=kvlGenerateToyExampleData.cxx
# kvlGenerateToyExampleData_LDADD=$(LIBS_ITK) \
#					 libkvlAtlas3DCommon.a \
#					 $(top_builddir)/itkio/libitkio.a \
#					 $(top_builddir)/itkutils/libitkutils.a \
#					 -lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
#					 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlGenerateToyExampleData_LDFLAGS= $(OS_LDFLAGS)

kvlSmoothMesh_SOURCES=kvlSmoothMesh.cxx
kvlSmoothMesh_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlSmoothMesh_LDFLAGS= $(OS_LDFLAGS)

kvlCrop_SOURCES=kvlCrop.cxx
kvlCrop_LDADD=$(LIBS_ITK) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlCrop_LDFLAGS= $(OS_LDFLAGS)

kvlChangeKOfMeshCollection_SOURCES=kvlChangeKOfMeshCollection.cxx
kvlChangeKOfMeshCollection_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlChangeKOfMeshCollection_LDFLAGS= $(OS_LDFLAGS)

kvlChangePointParametersOfMeshCollection_SOURCES=kvlChangePointParametersOfMeshCollection.cxx
kvlChangePointParametersOfMeshCollection_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlChangePointParametersOfMeshCollection_LDFLAGS= $(OS_LDFLAGS)

kvlMapDeformationField_SOURCES=kvlMapDeformationField.cxx
kvlMapDeformationField_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlMapDeformationField_LDFLAGS= $(OS_LDFLAGS)

kvlUpsampleMesh_SOURCES=kvlUpsampleMesh.cxx
kvlUpsampleMesh_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlUpsampleMesh_LDFLAGS= $(OS_LDFLAGS)

kvlSetReferencePosition_SOURCES=kvlSetReferencePosition.cxx
kvlSetReferencePosition_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlSetReferencePosition_LDFLAGS= $(OS_LDFLAGS)

kvlCompressImage_SOURCES=kvlCompressImage.cxx
kvlCompressImage_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlCompressImage_LDFLAGS= $(OS_LDFLAGS)

# kvlBuildAtlasMeshWithoutGUI_SOURCES=kvlBuildAtlasMeshWithoutGUI.cxx
# kvlBuildAtlasMeshWithoutGUI_LDADD=$(LIBS_ITK) \
# 				 libkvlAtlas3DCommon.a \
# 				 $(top_builddir)/itkio/libitkio.a \
# 				 $(top_builddir)/itkutils/libitkutils.a \
# 				 -lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
# 				 $(top_builddir)/tetgen/libtetgen.a \
# 				 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlBuildAtlasMeshWithoutGUI_LDFLAGS= $(OS_LDFLAGS)

# kvlConstructSparseMeshCollection_SOURCES=kvlConstructSparseMeshCollection.cxx
# kvlConstructSparseMeshCollection_LDADD=$(LIBS_ITK) \
#					 libkvlAtlas3DCommon.a \
#					 $(top_builddir)/itkio/libitkio.a \
#					 $(top_builddir)/itkutils/libitkutils.a \
#					 -lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
#					 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlConstructSparseMeshCollection_LDFLAGS= $(OS_LDFLAGS)

# kvlConstructMultiResolutionMeshCollection_SOURCES=kvlConstructMultiResolutionMeshCollection.cxx
# kvlConstructMultiResolutionMeshCollection_LDADD=$(LIBS_ITK) \
#					 libkvlAtlas3DCommon.a \
#					 $(top_builddir)/itkio/libitkio.a \
#					 $(top_builddir)/itkutils/libitkutils.a \
#					 -lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
#					 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlConstructMultiResolutionMeshCollection_LDFLAGS= $(OS_LDFLAGS)

kvlCropImage_SOURCES=kvlCropImage.cxx
kvlCropImage_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlCropImage_LDFLAGS= $(OS_LDFLAGS)

kvlCombineManualSegmentationProtocols_SOURCES=kvlCombineManualSegmentationProtocols.cxx
kvlCombineManualSegmentationProtocols_LDADD=$(LIBS_ITK) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlCombineManualSegmentationProtocols_LDFLAGS= $(OS_LDFLAGS)

kvlReduceNumberOfIntensityLevels_SOURCES=kvlReduceNumberOfIntensityLevels.cxx
kvlReduceNumberOfIntensityLevels_LDADD=$(LIBS_ITK) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKStatistics -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlReduceNumberOfIntensityLevels_LDFLAGS= $(OS_LDFLAGS)

kvlAugmentManualSegmentationWithAutomatedSegmentation_SOURCES=kvlAugmentManualSegmentationWithAutomatedSegmentation.cxx
kvlAugmentManualSegmentationWithAutomatedSegmentation_LDADD=$(LIBS_ITK) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlAugmentManualSegmentationWithAutomatedSegmentation_LDFLAGS= $(OS_LDFLAGS)

kvlScaleMeshCollection_SOURCES=kvlScaleMeshCollection.cxx
kvlScaleMeshCollection_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlScaleMeshCollection_LDFLAGS= $(OS_LDFLAGS)

kvlMapLabelsOfAtlasMesh_SOURCES=kvlMapLabelsOfAtlasMesh.cxx
kvlMapLabelsOfAtlasMesh_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlMapLabelsOfAtlasMesh_LDFLAGS= $(OS_LDFLAGS)

kvlSegmentWithoutGUI_SOURCES=kvlSegmentWithoutGUI.cxx
kvlSegmentWithoutGUI_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlSegmentWithoutGUI_LDFLAGS= $(OS_LDFLAGS)

kvlEvaluateQuantitatively_SOURCES=kvlEvaluateQuantitatively.cxx
kvlEvaluateQuantitatively_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlEvaluateQuantitatively_LDFLAGS= $(OS_LDFLAGS)

kvlEvaluateQuantitativelyWithHausdorffDistance_SOURCES=kvlEvaluateQuantitativelyWithHausdorffDistance.cxx
kvlEvaluateQuantitativelyWithHausdorffDistance_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlEvaluateQuantitativelyWithHausdorffDistance_LDFLAGS= $(OS_LDFLAGS)

kvlCorrectLeftSideLabelsToRightSideLabels_SOURCES=kvlCorrectLeftSideLabelsToRightSideLabels.cxx
kvlCorrectLeftSideLabelsToRightSideLabels_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlCorrectLeftSideLabelsToRightSideLabels_LDFLAGS= $(OS_LDFLAGS)

kvlOverwriteEmptySlices_SOURCES=kvlOverwriteEmptySlices.cxx
kvlOverwriteEmptySlices_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlOverwriteEmptySlices_LDFLAGS= $(OS_LDFLAGS)

kvlOptimizePositionsWithICM_SOURCES=kvlOptimizePositionsWithICM.cxx
kvlOptimizePositionsWithICM_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlOptimizePositionsWithICM_LDFLAGS= $(OS_LDFLAGS)

kvlCalculatePositionCost_SOURCES=kvlCalculatePositionCost.cxx
kvlCalculatePositionCost_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlCalculatePositionCost_LDFLAGS= $(OS_LDFLAGS)

kvlCalculateAtlasMeshCollectionTotalCost_SOURCES=kvlCalculateAtlasMeshCollectionTotalCost.cxx
kvlCalculateAtlasMeshCollectionTotalCost_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKAlgorithms -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlCalculateAtlasMeshCollectionTotalCost_LDFLAGS= $(OS_LDFLAGS)

# kvlTestImageSmoother_SOURCES=kvlTestImageSmoother.cxx
# kvlTestImageSmoother_LDADD=$(LIBS_ITK) \
#					 libkvlAtlas3DCommon.a \
#					 $(top_builddir)/itkio/libitkio.a \
#					 $(top_builddir)/itkutils/libitkutils.a \
#					 -lITKCommon -lITKIO -lITKznz \
#					 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlTestImageSmoother_LDFLAGS= $(OS_LDFLAGS)

# kvlGenerateWeightedColorChannels_SOURCES=kvlGenerateWeightedColorChannels.cxx
# kvlGenerateWeightedColorChannels_LDADD=$(LIBS_ITK) \
#					 libkvlAtlas3DCommon.a \
#					 $(top_builddir)/itkio/libitkio.a \
#					 $(top_builddir)/itkutils/libitkutils.a \
#					 -lITKCommon -lITKIO -lITKznz \
#					 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlGenerateWeightedColorChannels_LDFLAGS= $(OS_LDFLAGS)

kvlColorCodeProbabilityImages_SOURCES=kvlColorCodeProbabilityImages.cxx
kvlColorCodeProbabilityImages_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlColorCodeProbabilityImages_LDFLAGS= $(OS_LDFLAGS)

kvlSamplePositionsFromMeshCollection_SOURCES=kvlSamplePositionsFromMeshCollection.cxx
kvlSamplePositionsFromMeshCollection_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlSamplePositionsFromMeshCollection_LDFLAGS= $(OS_LDFLAGS)

kvlAddMeshToMeshCollection_SOURCES=kvlAddMeshToMeshCollection.cxx
kvlAddMeshToMeshCollection_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlAddMeshToMeshCollection_LDFLAGS= $(OS_LDFLAGS)

kvlMergeMeshCollections_SOURCES=kvlMergeMeshCollections.cxx
kvlMergeMeshCollections_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlMergeMeshCollections_LDFLAGS= $(OS_LDFLAGS)

# kvlAddCSFToBucknerData_SOURCES=kvlAddCSFToBucknerData.cxx
# kvlAddCSFToBucknerData_LDADD=$(LIBS_ITK) \
# 				 libkvlAtlas3DCommon.a \
# 				 $(top_builddir)/itkio/libitkio.a \
# 				 $(top_builddir)/itkutils/libitkutils.a \
# 				 -lITKCommon -lITKIO -lITKznz \
# 				 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
# kvlAddCSFToBucknerData_LDFLAGS= $(OS_LDFLAGS)

kvlAddBorder_SOURCES=kvlAddBorder.cxx
kvlAddBorder_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlAddBorder_LDFLAGS= $(OS_LDFLAGS)

kvlQuantifyPosteriorProbabilityImages_SOURCES=kvlQuantifyPosteriorProbabilityImages.cxx
kvlQuantifyPosteriorProbabilityImages_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlQuantifyPosteriorProbabilityImages_LDFLAGS= $(OS_LDFLAGS)

kvlEMSegment_SOURCES=kvlEMSegment.cxx
kvlEMSegment_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlEMSegment_LDFLAGS= $(OS_LDFLAGS)

kvlAddProbabilityImages_SOURCES=kvlAddProbabilityImages.cxx
kvlAddProbabilityImages_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlAddProbabilityImages_LDFLAGS= $(OS_LDFLAGS)

kvlMathematicalMorphology_SOURCES=kvlMathematicalMorphology.cxx
kvlMathematicalMorphology_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlMathematicalMorphology_LDFLAGS= $(OS_LDFLAGS)

kvlConvert_SOURCES=kvlConvert.cxx
kvlConvert_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlConvert_LDFLAGS= $(OS_LDFLAGS)

kvlRegister_SOURCES=kvlRegister.cxx
kvlRegister_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlRegister_LDFLAGS= $(OS_LDFLAGS)

kvlResample_SOURCES=kvlResample.cxx
kvlResample_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlResample_LDFLAGS= $(OS_LDFLAGS)

kvlBinaryThresholdImage_SOURCES=kvlBinaryThresholdImage.cxx
kvlBinaryThresholdImage_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlBinaryThresholdImage_LDFLAGS= $(OS_LDFLAGS)

kvlUpsample_SOURCES=kvlUpsample.cxx
kvlUpsample_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlUpsample_LDFLAGS= $(OS_LDFLAGS)

kvlRasterizeAtlasMesh_SOURCES=kvlRasterizeAtlasMesh.cxx
kvlRasterizeAtlasMesh_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlRasterizeAtlasMesh_LDFLAGS= $(OS_LDFLAGS)

kvlAutoCrop_SOURCES=kvlAutoCrop.cxx
kvlAutoCrop_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlAutoCrop_LDFLAGS= $(OS_LDFLAGS)

kvlMaskImage_SOURCES=kvlMaskImage.cxx
kvlMaskImage_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlMaskImage_LDFLAGS= $(OS_LDFLAGS)

kvlApplyTransform_SOURCES=kvlApplyTransform.cxx
kvlApplyTransform_LDADD=$(LIBS_ITK) \
	libkvlAtlas3DCommon.a \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	-lITKCommon -lITKIO -lITKznz \
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlApplyTransform_LDFLAGS= $(OS_LDFLAGS)

if HAVE_FLTK_LIBS
noinst_LIBRARIES+=libkvlAtlas3DGUI.a
libkvlAtlas3DGUI_a_SOURCES = vtkFlRenderWindowInteractor.cxx \
	kvlImageViewer.cxx
libkvlAtlas3DGUI_a_LIBADD=$(LIBS_ITK) $(LIBS_VTK)
SUFFIXES: .fl .cxx .h
.fl.h .fl.cxx:
	$(FLTK_DIR)/bin/fluid -c $<

bin_PROGRAMS+=kvlEstimateAtlasParameters

# Add GUI executable
kvlEstimateAtlasParameters_SOURCES = \
	kvlAtlasParameterEstimationConsoleGUI.fl \
	kvlEstimateAtlasParameters.cxx \
	kvlAtlasParameterEstimationConsole.cxx
kvlEstimateAtlasParameters_LDADD=$(LIBS_ITK) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	libkvlAtlas3DCommon.a \
	libkvlAtlas3DGUI.a \
	-lITKCommon -lITKBasicFilters -lITKAlgorithms -lITKIO -lITKznz \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
	$(VTK_SQLITE_LIB) \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB) \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lvtkWidgets -lvtkHybrid $(VTK_HDF5_LIB) \
	-lvtkIO	$(VTK_LSDYNA_LIB) -lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	$(FLTK_LIBS) $(FLTK_LDFLAGS)\
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlEstimateAtlasParameters_LDFLAGS = $(OS_LDFLAGS)


# Add GUI executable 
bin_PROGRAMS+=kvlThresholdImage

kvlThresholdImage_SOURCES = \
	kvlImageThresholdingConsoleGUI.fl \
	kvlThresholdImage.cxx \
	kvlImageThresholdingConsole.cxx \
	kvlImageThresholdingConsole.h

kvlThresholdImage_LDADD=$(LIBS_ITK) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	libkvlAtlas3DCommon.a \
	libkvlAtlas3DGUI.a \
	-lITKCommon -lITKBasicFilters -lITKAlgorithms -lITKIO -lITKznz \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
	$(VTK_SQLITE_LIB) \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB) \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lvtkWidgets -lvtkHybrid $(VTK_HDF5_LIB) \
	-lvtkIO	$(VTK_LSDYNA_LIB) -lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	$(FLTK_LIBS) $(FLTK_LDFLAGS)\
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlThresholdImage_LDFLAGS = $(OS_LDFLAGS)


# Add GUI executable 
bin_PROGRAMS+=kvlViewImage

kvlViewImage_SOURCES = \
	kvlImageViewingConsoleGUI.fl \
	kvlViewImage.cxx \
	kvlImageViewingConsole.cxx \
	kvlImageViewingConsole.h

kvlViewImage_LDADD=$(LIBS_ITK) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	libkvlAtlas3DCommon.a \
	libkvlAtlas3DGUI.a \
	-lITKCommon -lITKBasicFilters -lITKAlgorithms -lITKIO -lITKznz \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
	$(VTK_SQLITE_LIB) \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB) \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lvtkWidgets -lvtkHybrid $(VTK_HDF5_LIB) \
	-lvtkIO	$(VTK_LSDYNA_LIB) -lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	$(FLTK_LIBS) $(FLTK_LDFLAGS)\
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlViewImage_LDFLAGS = $(OS_LDFLAGS)


# Add GUI executable 
bin_PROGRAMS+=kvlViewMeshCollectionWithGUI

kvlViewMeshCollectionWithGUI_SOURCES = \
	kvlAtlasMeshViewingConsoleGUI.fl \
	kvlViewMeshCollectionWithGUI.cxx \
	kvlAtlasMeshViewingConsole.cxx
kvlViewMeshCollectionWithGUI_LDADD=$(LIBS_ITK) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	libkvlAtlas3DCommon.a \
	libkvlAtlas3DGUI.a \
	-lITKCommon -lITKBasicFilters -lITKAlgorithms -lITKIO -lITKznz \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
	$(VTK_SQLITE_LIB) \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB) \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lvtkWidgets -lvtkHybrid $(VTK_HDF5_LIB) \
	-lvtkIO	$(VTK_LSDYNA_LIB) -lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	$(FLTK_LIBS) $(FLTK_LDFLAGS)\
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlViewMeshCollectionWithGUI_LDFLAGS = $(OS_LDFLAGS)


## Add GUI executable
#bin_PROGRAMS+=kvlSegment
#
#kvlSegment_SOURCES = kvlAtlasMeshSegmenterConsoleGUI.fl \
#                     kvlSegment.cxx \
#                     kvlAtlasMeshSegmenterConsole.cxx
#kvlSegment_LDADD=$(LIBS_ITK) \
#                                 $(top_builddir)/itkio/libitkio.a \
#                                 $(top_builddir)/itkutils/libitkutils.a \
#                                 libkvlAtlas3DCommon.a \
#                                 libkvlAtlas3DGUI.a \
#                                 -lITKCommon -lITKBasicFilters -lITKAlgorithms -lITKIO -lITKznz \
#                                $(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
#                                $(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
#                                $(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
#                                $(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
#                                $(VTK_SQLITE_LIB) \
#                                -lvtkImaging -lvtkFiltering \
#                                -lvtkCommon -lvtksys \
#                                -lvtkexoIIc -lvtkNetCDF \
#                                -lvtkVolumeRendering -lvtkRendering -lvtkftgl \
#                                -lvtkWidgets -lvtkHybrid \
#                                -lvtkIO -lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
#                                 $(FLTK_LIBS) $(FLTK_LDFLAGS)\
#                                 $(addprefix $(top_builddir)/, $(LIBS_MGH) )
#kvlSegment_LDFLAGS = $(OS_LDFLAGS)


# Add GUI executable
bin_PROGRAMS+=kvlEvaluateSegmentation

kvlEvaluateSegmentation_SOURCES = \
	kvlSegmentationEvaluationConsoleGUI.fl \
	kvlEvaluateSegmentation.cxx \
	kvlSegmentationEvaluationConsole.cxx
kvlEvaluateSegmentation_LDADD=$(LIBS_ITK) \
	$(top_builddir)/itkio/libitkio.a \
	$(top_builddir)/itkutils/libitkutils.a \
	libkvlAtlas3DCommon.a \
	libkvlAtlas3DGUI.a \
	-lITKCommon -lITKBasicFilters -lITKAlgorithms -lITKIO -lITKznz \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
	$(VTK_SQLITE_LIB) \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB) \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lvtkWidgets -lvtkHybrid $(VTK_HDF5_LIB) \
	-lvtkIO	$(VTK_LSDYNA_LIB) -lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	$(FLTK_LIBS) $(FLTK_LDFLAGS)\
	$(addprefix $(top_builddir)/, $(LIBS_MGH) )
kvlEvaluateSegmentation_LDFLAGS = $(OS_LDFLAGS)


# put a wrapper around the gui bins, used to setup vtk enviro vars
install-exec-hook:kvlThresholdImage kvlEstimateAtlasParameters
	cp $(top_builddir)/GEMS/kvlThresholdImage \
		$(DESTDIR)$(bindir)/kvlThresholdImage.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/kvlThresholdImage
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/kvlThresholdImage
	echo "if (-e /bin/pwd) then" >> $(DESTDIR)$(bindir)/kvlThresholdImage
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlThresholdImage
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlThresholdImage
	echo "endif" >> $(DESTDIR)$(bindir)/kvlThresholdImage
	echo "kvlThresholdImage.bin \$$argv" >> $(DESTDIR)$(bindir)/kvlThresholdImage
	chmod a+x $(DESTDIR)$(bindir)/kvlThresholdImage
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/kvlThresholdImage.bin; fi
	cp $(top_builddir)/GEMS/kvlEstimateAtlasParameters \
		$(DESTDIR)$(bindir)/kvlEstimateAtlasParameters.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/kvlEstimateAtlasParameters
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/kvlEstimateAtlasParameters
	echo "if (-e /bin/pwd) then" >> \
		$(DESTDIR)$(bindir)/kvlEstimateAtlasParameters
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlEstimateAtlasParameters
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlEstimateAtlasParameters
	echo "endif" >> $(DESTDIR)$(bindir)/kvlEstimateAtlasParameters
	echo "kvlEstimateAtlasParameters.bin \$$argv" >> \
		$(DESTDIR)$(bindir)/kvlEstimateAtlasParameters
	chmod a+x $(DESTDIR)$(bindir)/kvlEstimateAtlasParameters
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d \
		$(DESTDIR)$(bindir)/kvlEstimateAtlasParameters.bin; fi
	cp $(top_builddir)/GEMS/kvlViewImage \
		$(DESTDIR)$(bindir)/kvlViewImage.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/kvlViewImage
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/kvlViewImage
	echo "if (-e /bin/pwd) then" >> $(DESTDIR)$(bindir)/kvlViewImage
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlViewImage
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlViewImage
	echo "endif" >> $(DESTDIR)$(bindir)/kvlViewImage
	echo "kvlViewImage.bin \$$argv" >> $(DESTDIR)$(bindir)/kvlViewImage
	chmod a+x $(DESTDIR)$(bindir)/kvlViewImage
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/kvlViewImage.bin; fi
	cp $(top_builddir)/GEMS/kvlViewMeshCollectionWithGUI \
		$(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI
	echo "if (-e /bin/pwd) then" >> $(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI
	echo "endif" >> $(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI
	echo "kvlViewMeshCollectionWithGUI.bin \$$argv" >> $(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI
	chmod a+x $(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI.bin; fi
	cp $(top_builddir)/GEMS/kvlEvaluateSegmentation \
		$(DESTDIR)$(bindir)/kvlEvaluateSegmentation.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/kvlEvaluateSegmentation
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/kvlEvaluateSegmentation
	echo "if (-e /bin/pwd) then" >> $(DESTDIR)$(bindir)/kvlEvaluateSegmentation
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlEvaluateSegmentation
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/kvlEvaluateSegmentation
	echo "endif" >> $(DESTDIR)$(bindir)/kvlEvaluateSegmentation
	echo "kvlEvaluateSegmentation.bin \$$argv" >> $(DESTDIR)$(bindir)/kvlEvaluateSegmentation
	chmod a+x $(DESTDIR)$(bindir)/kvlEvaluateSegmentation
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/kvlEvaluateSegmentation.bin; fi 
	cp $(top_builddir)/GEMS/kvlApplyTransform $(DESTDIR)$(bindir)/kvlApplyTransform.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/kvlApplyTransform
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/kvlApplyTransform
	echo "kvlApplyTransform.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/kvlApplyTransform
	cp $(top_builddir)/GEMS/kvlAutoCrop $(DESTDIR)$(bindir)/kvlAutoCrop.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/kvlAutoCrop
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/kvlAutoCrop
	echo "kvlAutoCrop.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/kvlAutoCrop
	cp $(top_builddir)/GEMS/kvlRegister $(DESTDIR)$(bindir)/kvlRegister.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/kvlRegister
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/kvlRegister
	echo "kvlRegister.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/kvlRegister

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/kvlThresholdImage.bin \
				$(DESTDIR)$(bindir)/kvlEstimateAtlasParameters.bin \
				$(DESTDIR)$(bindir)/kvlViewImage.bin \
				$(DESTDIR)$(bindir)/kvlViewMeshCollectionWithGUI.bin \
				$(DESTDIR)$(bindir)/kvlEvaluateSegmentation.bin

EXTRA_DIST=\
	gzstream.h \
	itkAffine3DTransform.h \
	kvlParameterOrderPowellOptimizer.h \
	kvlAtlasMeshAlphaDrawer.h \
	kvlAtlasMeshBuilder.h \
	kvlAtlasMeshCollectionFastReferencePositionCost.h \
	kvlAtlasMeshCollection.h \
	kvlAtlasMeshCollectionModelLikelihoodCalculator.h \
	kvlAtlasMeshCollectionPositionCostCalculator2.h \
	kvlAtlasMeshCollectionPositionCostCalculator.h \
	kvlAtlasMeshCollectionReferencePositionCost.h \
	kvlAtlasMeshCollectionValidator.h \
	kvlAtlasMeshDeformationLevenbergMarquardt.h \
	kvlAtlasMesh.h \
	kvlAtlasMeshHamiltonianPositionSampler.h \
	kvlAtlasMeshLabelStatisticsCollector.h \
	kvlAtlasMeshMinLogLikelihoodCalculator.h \
	kvlAtlasMeshMultiAlphaDrawer.h \
	kvlAtlasMeshPositionGradientCalculator.h \
	kvlAtlasMeshRasterizor.h \
	kvlAtlasMeshSegmentationDriver.h \
	kvlAtlasMeshSegmenterConsole.h \
	kvlAtlasMeshSegmenter.h \
	kvlAtlasMeshSmoother.h \
	kvlAtlasMeshSummaryDrawer.h \
	kvlAtlasMeshToIntensityImageGradientCalculator.h \
	kvlAtlasMeshToIntensityImagePartialVolumeGradientCalculator.h \
	kvlAtlasMeshToProbabilityImageGradientCalculator.h \
	kvlAtlasMeshVertexProcessor2.h \
	kvlAtlasMeshVertexProcessor.h \
	kvlAtlasMeshViewingConsole.h \
	kvlAtlasMeshVisitCounter.h \
	kvlAtlasParameterEstimationConsoleGUI.h \
	kvlAtlasParameterEstimationConsole.h \
	kvlAtlasParameterEstimator.h \
	kvlCompressionLookupTable.h \
	kvlCroppedImageReader.h \
	kvlEMSegmenter.h \
	kvlImageSmoother.h \
	kvlImageThresholdingConsole.h \
	kvlImageViewer.h \
	kvlImageViewingConsole.h \
	kvlMultiResolutionAtlasMesher.h \
	kvlProgressReporter.h \
	kvlRegisterer.h \
	kvlResampler.h \
	kvlSegmentationEvaluationConsole.h \
	kvlTetrahedronAspectRatio.h \
	vnl_kvl_qr.h \
	vtkFlRenderWindowInteractor.h \
	itkAffine3DTransform.txx \
	kvlAtlasMeshRasterizor.txx \
	vnl_kvl_qr.txx

endif
endif

clean-local:
	rm -f kvlAtlasMeshViewingConsoleGUI.cxx kvlImageThresholdingConsoleGUI.cxx \
	kvlImageViewingConsoleGUI.cxx kvlSegmentationEvaluationConsoleGUI.cxx \
	rm kvlImageThresholdingConsoleGUI.h kvlImageViewingConsoleGUI.h

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

