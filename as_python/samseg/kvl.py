class KVL:
    def __init__(self, number_of_threads):
        self.number_of_threads = number_of_threads
        # % Clean up
        # kvlClear; % Clear all the wrapped C++ stuff
        # close all;
        #
        # % Specify the maximum number of threads the C++ stuff will use. The more threads you can use
        # % the faster, provided you have a matching amount of cores on your system - up to
        # % the point where memory becomes the bottle neck.
        # % If the following command is not provided, the number of cores on your system will be used
        # kvlSetMaximumNumberOfThreads( numberOfThreads );
        pass

    ## Example uses:

    # As used in run_samsegment:
    #   [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
    #   meshCollection = kvlReadMeshCollection( meshCollectionFileName );
    #   mesh = kvlGetMesh( meshCollection, -1 );
    #   [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
    #   alphas = kvlGetAlphasInMeshNodes( mesh );
    #   [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );
    #   [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, sharedGMMParameters, FreeSurferLabels, colors );
    #   [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );
    #   % Create tailored atlas
    #   kvlSetAlphasInMeshNodes( mesh, alphas );
    #   [ template, transform ] = kvlReadImage( templateFileName );
    #   templateImageBuffer = kvlGetImageBuffer( template );
    #   priors = kvlRasterizeAtlasMesh( mesh, size( templateImageBuffer ) );
    #   transformMatrix = kvlGetTransformMatrix( transform );
    # kvlWriteImage( origTemplate, transformedTemplateFileName, ...
    #                kvlCreateTransform( double( worldToWorldTransformMatrix * kvlGetTransformMatrix( origTemplateTransform ) ) ) );

    # As used in samsseg_registerAtlas:
    # % Read in image, and figure out where the mesh is located
    # [ image, imageToWorldTransform ] = kvlReadImage( imageFileName );
    # imageToWorldTransformMatrix = double( kvlGetTransformMatrix( imageToWorldTransform ) );
    #
    # [ template, templateImageToWorldTransform ] = kvlReadImage( templateFileName );
    # templateImageToWorldTransformMatrix = double( kvlGetTransformMatrix( templateImageToWorldTransform ) );
    #   meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
    #                                           kvlCreateTransform( initialImageToImageTransformMatrix ), ...
    #                                           K * prod( downSamplingFactors ) );
    #   mesh = kvlGetMesh( meshCollection, -1 );
    #   meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
    #                                           kvlCreateTransform( imageToWorldTransformMatrix \ ...
    #                                                               templateImageToWorldTransformMatrix ), ...
    #                                           K * prod( downSamplingFactors ) );
    #   mesh = kvlGetMesh( meshCollection, -1 );
    #   nodePositions = kvlGetMeshNodePositions( mesh );
    #   kvlSetMeshNodePositions( mesh, nodePositions );
    # imageBuffer = kvlGetImageBuffer( image );
    # image = kvlCreateImage( imageBuffer );
    # kvlScaleMesh( mesh, 1 ./ downSamplingFactors );
    # alphas = kvlGetAlphasInMeshNodes( mesh );
    #   kvlSetAlphasInMeshNodes( mesh, alphas )
    #   priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    # calculator = kvlGetCostAndGradientCalculator( 'MutualInformation', ...
    #                                                image, 'Affine' );
    # [ cost gradient ] = kvlEvaluateMeshPosition( calculator, mesh );
    #   priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    #   nodePositions = kvlGetMeshNodePositions( mesh );
    #   kvlSetMeshNodePositions( mesh, trialNodePositions );
    #   [ trialCost trialGradient ] = kvlEvaluateMeshPosition( calculator, mesh );
    #     kvlSetMeshNodePositions( mesh, nodePositions );
    # originalNodePositions = kvlGetMeshNodePositions( mesh );
    #   priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    #   colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
    # optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
    #                                 'Verbose', 1, ...
    #                                 'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ...
    #                                 'LineSearchMaximalDeformationIntervalStopCriterion', ...
    #                                 lineSearchMaximalDeformationIntervalStopCriterion, ...
    #                                 'BFGS-MaximumMemoryLength', 12 ); % Affine registration only has 12 DOF
    #   [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer )
    #     priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    #     colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
    #     nodePositions = kvlGetMeshNodePositions( mesh );
    # nodePositions = kvlGetMeshNodePositions( mesh );
    # kvlWriteImage( template, transformedTemplateFileName, ...
    #                kvlCreateTransform( desiredTemplateImageToWorldTransformMatrix ) );
    # priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    # colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );

    # As used in createAtlasMeshForAffineRegistration
    #   priors( :, :, :, 1 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/background.nii' ) );
    #   priors( :, :, :, 2 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/wm.nii' ) );
    #   priors( :, :, :, 3 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/gm.nii' ) );
    #   priors( :, :, :, 4 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/csf.nii' ) );
    #   if ( numberOfClasses > 4 )
    #     priors( :, :, :, 5 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/bone.nii' ) );
    #     priors( :, :, :, 6 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/soft.nii' ) );
    #   end
    #
    #   % Also get the image-to-world transform of this template
    #   [ ~, transform ] = kvlReadImage( '/data/testing/atlas/wm.nii' );
    #   transformMatrix = kvlGetTransformMatrix( transform );
    # meshCollection = kvlCreateMeshCollection( priors, meshSize );
    # kvlWriteMeshCollection( meshCollection, affineRegistrationMeshCollectionFileName );
    # template = kvlCreateImage( single( templateImageBuffer ) );
    # affineRegistrationTemplateFileName = [ outputFileNameBase '_template.nii' ];
    # transform = kvlCreateTransform( double( transformMatrix ) );
    # kvlWriteImage( template, affineRegistrationTemplateFileName, transform );

    # As used in samsegment:
    #   [ images( contrastNumber ), transform, nonCroppedImageSize, croppingOffset ] = ...
    #                                     kvlReadCroppedImage( imageFileNames{ contrastNumber }, transformedTemplateFileName );
    #   imageBuffers( :, :, :, contrastNumber ) = kvlGetImageBuffer( images( contrastNumber ) ); % Get the actual imageBuffer
    # [ ~, imageToWorldTransform ] = kvlReadImage( imageFileNames{1} );
    # imageToWorldTransformMatrix = kvlGetTransformMatrix( imageToWorldTransform );
    # meshCollection = kvlReadMeshCollection( meshCollectionFileName, transform, modelSpecifications.K );
    # mesh = kvlGetMesh( meshCollection, -1 );
    # backgroundPrior = kvlRasterizeAtlasMesh( mesh, imageSize, labelNumber );
    # smoothedBackgroundPrior = kvlSmoothImageBuffer( backgroundPrior, modelSpecifications.brainMaskingSmoothingSigma ./ voxelSpacing );
    #   % kvlSetImageBuffer( images( contrastNumber ), imageBuffers( :, :, :, contrastNumber ) );
    # [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
    # alphas = kvlGetAlphasInMeshNodes( mesh );
    # [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );
    # kvlSetAlphasInMeshNodes( mesh, alphas );
    # [ reducedAlphas, reducedNames, reducedFreeSurferLabels, reducedColors, reducingLookupTable ] = ...
    #                             kvlMergeAlphas( alphas, names, modelSpecifications.sharedGMMParameters, FreeSurferLabels, colors );
    #   priors = kvlRasterizeAtlasMesh( mesh, imageSize ); % Without specifying a specific label, will rasterize all simultaneously, rasterize everything
    #   rgbBuffer = kvlColorCodeProbabilityImages( priors, colors );
    #   % Setting priors in mesh to the inital values that also happen to
    #   % be reduced
    #   kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
    #       smoothedBuffer = kvlSmoothImageBuffer( single( buffer ), smoothingSigmas );
    #       smoothedMask = kvlSmoothImageBuffer( single( downSampledMask ), smoothingSigmas );
    #   kvlScaleMesh( mesh, 1./downSamplingFactors );
    #   kvlSmoothMesh( mesh, smoothingSigmas );
    #   smoothedReducedAlphas = kvlGetAlphasInMeshNodes( mesh );
    #     oldColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, downSampledImageSize ), reducedColors );
    #     tmp = reshape( kvlRasterizeAtlasMesh( mesh, downSampledImageSize ), [ prod( downSampledImageSize ) numberOfClasses ] );
    #         kvlClear( downSampledBiasCorrectedImages( contrastNumber ) );
    #              kvlCreateImage( single( downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) ) );
    #     calculator = kvlGetCostAndGradientCalculator( 'AtlasMeshToIntensityImage', ...
    #                                                    downSampledBiasCorrectedImages, ...
    #                                                    'Sliding', ...
    #                                                    transform, ...
    #                                                    means, variances, mixtureWeights, numberOfGaussiansPerClass );
    #     optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
    #                                     'Verbose', optimizationOptions.verbose, ...
    #                                     'MaximalDeformationStopCriterion', optimizationOptions.maximalDeformationStopCriterion, ...
    #                                     'LineSearchMaximalDeformationIntervalStopCriterion', ...
    #                                       optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion, ...
    #                                     'MaximumNumberOfIterations', optimizationOptions.maximumNumberOfDeformationIterations, ...
    #                                     'BFGS-MaximumMemoryLength', optimizationOptions.BFGSMaximumMemoryLength );
    #     nodePositionsBeforeDeformation = kvlGetMeshNodePositions( mesh );
    #       [ minLogLikelihoodTimesDeformationPrior, maximalDeformation ] = kvlStepOptimizer( optimizer );
    #     kvlClear( calculator );
    #     kvlClear( optimizer );
    #     nodePositionsAfterDeformation = kvlGetMeshNodePositions( mesh );
    #       newColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, downSampledImageSize ), reducedColors );
    #   kvlScaleMesh( mesh, downSamplingFactors );
    # kvlSetAlphasInMeshNodes( mesh, alphas )
    # priors = kvlRasterizeAtlasMesh( mesh, imageSize );
    # kvlWriteImage( kvlCreateImage( uncroppedFreeSurferSegmentation ), ...
    #                fullfile( savePath, 'crispSegmentation.nii' ), ...
    #                imageToWorldTransform );
    #   kvlWriteImage( kvlCreateImage( biasField ), outputFileName, imageToWorldTransform );
    #   kvlWriteImage( kvlCreateImage( biasCorrected ), outputFileName, imageToWorldTransform );
