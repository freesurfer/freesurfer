# ﻿
# % We do the optimization in a multi-resolution type of scheme, where large
# % deformations are quickly found using smoothed versions of the atlas mesh, and the fine
# % details are then found on gradually less smoothed versions until the original atlas mesh is used for the optimization.
# % Multi-resolution is a standard procedure in image registration: the initial
# % blurring avoids getting stuck in the first local optimum of the search space, and get the rough major
# % deformations right instead.
# numberOfMultiResolutionLevels = length( optimizationOptions.multiResolutionSpecification );
# historyWithinEachMultiResolutionLevel = struct( [] );
# for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
#
#   %
#   maximumNumberOfIterations = optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).maximumNumberOfIterations;
#   estimateBiasField = optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).estimateBiasField;
#   historyOfCost = [ 1/eps ];
#   historyOfMaximalDeformationApplied = [];
#   historyOfTimeTakenIntensityParameterUpdating = [];
#   historyOfTimeTakenDeformationUpdating = [];
#   fprintf('maximumNumberOfIterations %d\n',maximumNumberOfIterations);
#
#
#   % Downsample the images, the mask, the mesh, and the bias field basis functions
#   % Must be integer
#   downSamplingFactors = max( round( optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).targetDownsampledVoxelSpacing ...
#                                     ./ voxelSpacing ), [ 1 1 1 ] )
#   downSampledMask = mask(  1 : downSamplingFactors( 1 ) : end, ...
#                            1 : downSamplingFactors( 2 ) : end, ...
#                            1 : downSamplingFactors( 3 ) : end );
#   downSampledMaskIndices = find( downSampledMask );
#   downSampledImageBuffers = [];
#   for contrastNumber = 1 : numberOfContrasts
#     % if ( multiResolutionLevel == numberOfMultiResolutionLevels )
#     if true
#       % No image smoothing
#       downSampledImageBuffers( :, :, :, contrastNumber ) = imageBuffers( 1 : downSamplingFactors( 1 ) : end, ...
#                                                                          1 : downSamplingFactors( 2 ) : end, ...
#                                                                          1 : downSamplingFactors( 3 ) : end, ...
#                                                                          contrastNumber );
#     else
#       % Try image smoothing
#       buffer = imageBuffers( 1 : downSamplingFactors( 1 ) : end, ...
#                              1 : downSamplingFactors( 2 ) : end, ...
#                              1 : downSamplingFactors( 3 ) : end, ...
#                              contrastNumber );
#       smoothingSigmas = downSamplingFactors / 2 / sqrt( 2 * log( 2 ) ); % Variance chosen to approximately
#                                                                         % match normalized binomial filter
#                                                                         % (1/4, 1/2, 1/4) for downsampling
#                                                                         % factor of 2
#       smoothingSigmas( find( downSamplingFactors == 1 ) ) = 0.0;
#       smoothedBuffer = kvlSmoothImageBuffer( single( buffer ), smoothingSigmas );
#       smoothedMask = kvlSmoothImageBuffer( single( downSampledMask ), smoothingSigmas );
#       downSampledImageBuffers( :, :, :, contrastNumber ) = downSampledMask .* ( smoothedBuffer ./ ( smoothedMask + eps ) );
#     end
#
#   end
#
#   %
#   downSampledKroneckerProductBasisFunctions = cell( 0, 0 );
#   for dimensionNumber = 1 : 3
#     A = kroneckerProductBasisFunctions{ dimensionNumber };
#     downSampledKroneckerProductBasisFunctions{ dimensionNumber } = A( 1 : downSamplingFactors( dimensionNumber ) : end, : );
#   end
#   downSampledImageSize = size( downSampledImageBuffers( :, :, :, 1 ) );
#
#
#   % Read the atlas mesh to be used for this multi-resolution level, taking into account the downsampling to position it
#   % correctly
#   downSamplingTransformMatrix = diag( [ 1./downSamplingFactors 1 ] );
#   totalTransformationMatrix = downSamplingTransformMatrix * double( kvlGetTransformMatrix( transform ) );
#   meshCollection = ...
#         kvlReadMeshCollection( optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).atlasFileName, ...
#                                 kvlCreateTransform( totalTransformationMatrix ), modelSpecifications.K );
#   mesh = kvlGetMesh( meshCollection, -1 );
#
#   % Get the initial mesh node positions, also transforming them back into template space
#   % (i.e., undoing the affine registration that we applied) for later usage
#   initialNodePositions = kvlGetMeshNodePositions( mesh );
#   numberOfNodes = size( initialNodePositions, 1 );
#   tmp = ( totalTransformationMatrix \ [ initialNodePositions ones( numberOfNodes, 1 ) ]' )';
#   initialNodePositionsInTemplateSpace = tmp( :, 1 : 3 );
#
#
#   % If this is not the first multi-resolution level, apply the warp computed during the previous level
#   if ( multiResolutionLevel > 1 )
#     % Get the warp in template space
#     nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel = ...
#             historyWithinEachMultiResolutionLevel( multiResolutionLevel-1 ).finalNodePositionsInTemplateSpace - ...
#             historyWithinEachMultiResolutionLevel( multiResolutionLevel-1 ).initialNodePositionsInTemplateSpace;
#     initialNodeDeformationInTemplateSpace = kvlWarpMesh( ...
#                   optimizationOptions.multiResolutionSpecification( multiResolutionLevel-1 ).atlasFileName, ...
#                   nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel, ...
#                   optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).atlasFileName );
#
#     % Apply this warp on the mesh node positions in template space, and transform into current space
#     desiredNodePositionsInTemplateSpace = initialNodePositionsInTemplateSpace + initialNodeDeformationInTemplateSpace;
#     tmp = ( totalTransformationMatrix * ...
#             [ desiredNodePositionsInTemplateSpace ones( numberOfNodes, 1 ) ]' )';
#     desiredNodePositions = tmp( :, 1 : 3 );
#
#     %
#     kvlSetMeshNodePositions( mesh, desiredNodePositions );
#
#   end
#
#
#
#   % Set priors in mesh to the reduced (super-structure) ones
#   alphas = kvlGetAlphasInMeshNodes( mesh );
#   reducedAlphas = kvlMergeAlphas( alphas, names, modelSpecifications.sharedGMMParameters, FreeSurferLabels, colors );
#   kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
#
#
#
#
#   % Algorithm-wise, we're just estimating sets of parameters for one given data (MR scan) that is
#   % known and fixed throughout. However, in terms of bias field correction it will be computationally
#   % more efficient to pre-compute the bias field corrected version of the scan ("corrected" with
#   % the current estimate of the bias field) once and pass that on to different routines instead of the
#   % original data.
#   % For convenience (although potentially a recipe for future bug introduction), I'm also keeping a
#   % vectorized form of that around -- this will be useful in various places in the EM-parts. So
#   % effectively I have two redundant variables "downSampledBiasCorrectedImageBuffers" and "biasCorrectedData"
#   % that really just encode the variable "biasFieldCoefficients" and so need to be meticiously updated each time
#   % "biasFieldCoefficients" is updated (!)
#   downSampledBiasCorrectedImageBuffers = zeros( [ downSampledImageSize numberOfContrasts ] );
#   biasCorrectedData = zeros( [ length( downSampledMaskIndices ) numberOfContrasts ] );
#   for contrastNumber = 1 : numberOfContrasts
#     downSampledBiasField = backprojectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
#     tmp = downSampledImageBuffers( :, :, :, contrastNumber ) - downSampledBiasField .* downSampledMask;
#     downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) = tmp;
#     biasCorrectedData( :, contrastNumber ) = tmp( downSampledMaskIndices );
#   end
#
#
#   % Compute a color coded version of the atlas prior in the atlas's current pose, i.e., *before*
#   % we start deforming. We'll use this just for visualization purposes
#   if ( showFigures )
#     oldColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, downSampledImageSize ), reducedColors );
#   end
#
#
#
#   historyWithinEachIteration = struct( [] );
#   priors = zeros( length( downSampledMaskIndices ), numberOfClasses );
#   posteriors = zeros( length( downSampledMaskIndices ), numberOfGaussians ); % Gaussian mixture models burst out into
#                                                                              % individual Gaussian components
#
#   % Easier to work with vector notation in the EM computations
#   % reshape into a matrix
#   data = zeros( [ length( downSampledMaskIndices ) numberOfContrasts ] );
#   for contrastNumber = 1:numberOfContrasts
#     tmp = reshape( downSampledImageBuffers( :, :, :, contrastNumber ), [ prod(downSampledImageSize) 1 ] );
#     data( :, contrastNumber ) = tmp( downSampledMaskIndices );
#   end
#
#   % Main iteration loop over both EM and deformation
#   for iterationNumber = 1 : maximumNumberOfIterations
#
#     %
#     startTimeIntensityParameterUpdating = tic;
#
#     %
#     % Part I: estimate Gaussian mixture model parameters, as well as bias field parameters using EM.
#     %
#
#     % Get the priors at the current mesh position
#     tmp = reshape( kvlRasterizeAtlasMesh( mesh, downSampledImageSize ), [ prod( downSampledImageSize ) numberOfClasses ] );
#     priors( : ) = double( tmp( downSampledMaskIndices, : ) ) / 65535;
#
#     if ( iterationNumber == 1 )
#       historyWithinEachMultiResolutionLevel( multiResolutionLevel ).priorsAtStart = priors;
#     end
#
#
#     % Start EM iterations.
#     if ( ( multiResolutionLevel == 1 ) && ( iterationNumber == 1 ) )
#
#       % Initialize the mixture parameters if this is the first time ever you run this
#       means = zeros( numberOfGaussians, numberOfContrasts );
#       variances = zeros( numberOfGaussians, numberOfContrasts, numberOfContrasts );
#       mixtureWeights = zeros( numberOfGaussians, 1 );
#       for classNumber = 1 : numberOfClasses
#
#         % Calculate the global weighted mean and variance of this class, where the weights are given by the prior
#         prior = priors( :, classNumber );
#         mean = data' * prior / sum( prior );
#         tmp = data - repmat( mean', [ size( data, 1 ) 1 ] );
#         variance = tmp' * ( tmp .* repmat( prior, [ 1 numberOfContrasts ] ) ) / sum( prior );
#         if modelSpecifications.useDiagonalCovarianceMatrices
#           % Force diagonal covariance matrices
#           variance = diag( diag( variance ) );
#         end
#
#
#         % Based on this, initialize the mean and variance of the individual Gaussian components in this class'
#         % mixture model: variances are simply copied from the global class variance, whereas the means are
#         % determined by splitting the [ mean-sqrt( variance ) mean+sqrt( variance ) ] domain into equal intervals,
#         % the middle of which are taken to be the means of the Gaussians. Mixture weights are initialized to be
#         % all equal.
#         %
#         % This actually creates a mixture model that mimics the single Gaussian quite OK-ish: to visualize this
#         % do e.g.,
#         %
#         %  for numberOfComponents = 1 : 7
#         %    intervalSize = 2 / numberOfComponents;
#         %    means = -1 + intervalSize/2 + [ 0 : numberOfComponents-1 ] * intervalSize;
#         %    figure
#         %    x = [ -6 : .1 : 6 ];
#         %    gauss = exp( -x.^2/2 );
#         %    plot( gauss )
#         %    hold on
#         %    mixture = zeros( size( x ) );
#         %    for i = 1 : numberOfComponents
#         %      gauss = exp( -( x - means( i ) ).^2/2 );
#         %      plot( gauss / numberOfComponents, 'g' )
#         %      mixture = mixture + gauss / numberOfComponents;
#         %    end
#         %    plot( mixture, 'r' )
#         %    grid
#         %    title( [ num2str( numberOfComponents ) ' components' ] )
#         %  end
#         %
#         numberOfComponents = numberOfGaussiansPerClass( classNumber );
#         for componentNumber = 1 : numberOfComponents
#           gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
#
#           variances( gaussianNumber, :, : ) = variance;
#           intervalSize = 2 * sqrt( diag( variance ) ) / numberOfComponents;
#           means( gaussianNumber, : ) = ( mean - sqrt( diag( variance ) ) + intervalSize/2 + ( componentNumber - 1 ) * intervalSize )';
#           mixtureWeights( gaussianNumber ) = 1 / numberOfComponents;
#         end
#
#       end % End loop over classes
#
#
#       % Also remember the overall data variance for later usage in a conjugate prior on the variances
#       dataMean = sum( data )' / size( data, 1 );
#       tmp = data - repmat( dataMean', [ size( data, 1 ) 1 ] );
#       dataVariance = diag( diag( tmp' * tmp ) ) / size( data, 1 );
#       numberOfPseudoMeasurementsOfWishartPrior = 1; % In Oula's code this was effectively 2 * ( numberOfContrasts + 2 )
#                                                     % although I have no clue why
#       pseudoVarianceOfWishartPrior = dataVariance / numberOfPseudoMeasurementsOfWishartPrior;
#
#     end % End test need for initialization
#
#     % stopCriterionEM = 1e-5;
#     historyOfEMCost = [ 1/eps ];
#     for EMIterationNumber = 1 : 100
#       %
#       % E-step: compute the posteriors based on the current parameters.
#       %
#       for classNumber = 1 : numberOfClasses
#         prior = priors( :, classNumber );
#
#         numberOfComponents = numberOfGaussiansPerClass( classNumber );
#         for componentNumber = 1 : numberOfComponents
#           gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
#
#           mean = means( gaussianNumber, : )';
#           variance = squeeze( variances( gaussianNumber, :, : ) );
#           L = chol( variance, 'lower' );  % variance = L * L'
#           tmp = L \ ( biasCorrectedData' - repmat( mean, [ 1 size( biasCorrectedData, 1 ) ] ) );
#           squaredMahalanobisDistances = ( sum( tmp.^2, 1 ) )';
#           sqrtDeterminantOfVariance = prod( diag( L ) ); % Same as sqrt( det( variance ) )
#           gaussianLikelihoods = exp( -squaredMahalanobisDistances / 2 ) / ( 2 * pi )^( numberOfContrasts / 2 ) / sqrtDeterminantOfVariance;
#
#           posteriors( :, gaussianNumber ) = gaussianLikelihoods .* ( mixtureWeights( gaussianNumber ) * prior );
#
#         end % End loop over mixture components
#
#       end % End loop over classes
#       normalizer = sum( posteriors, 2 ) + eps;
#       if 0
#         x = zeros( downSampledImageSize );
#         x( downSampledMaskIndices ) = -log( normalizer );
#         figure
#         showImage( x )
#       end
#       posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfGaussians ] );
#       minLogLikelihood = -sum( log( normalizer ) );
#       intensityModelParameterCost = 0;
#       for gaussianNumber = 1 : numberOfGaussians
#         variance = squeeze( variances( gaussianNumber, :, : ) );
#
#         % Evaluate unnormalized Wishart distribution (conjugate prior on precisions) with parameters
#         %
#         %   scale matrix V = inv( pseudoVarianceOfWishartPrior * numberOfPseudoMeasurementsOfWishartPrior )
#         %
#         % and
#         %
#         %   degrees of freedom n = numberOfPseudoMeasurementsOfWishartPrior + numberOfContrasts + 1
#         %
#         % which has pseudoVarianceOfWishartPrior as the MAP solution in the absence of any data
#         %
#         minLogUnnormalizedWishart = ...
#             trace( variance \ pseudoVarianceOfWishartPrior ) * numberOfPseudoMeasurementsOfWishartPrior / 2 + ...
#             numberOfPseudoMeasurementsOfWishartPrior / 2 * log( det( variance ) );
#         intensityModelParameterCost = intensityModelParameterCost + minLogUnnormalizedWishart;
#       end
#       historyOfEMCost = [ historyOfEMCost; minLogLikelihood + intensityModelParameterCost ];
#
#       % Show some figures
#       if ( showFigures )
#         for classNumber = 1 : numberOfClasses
#           posterior = zeros( downSampledImageSize );
#           numberOfComponents = numberOfGaussiansPerClass( classNumber );
#           for componentNumber = 1 : numberOfComponents
#             gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
#             posterior( downSampledMaskIndices ) = posterior( downSampledMaskIndices ) + ...
#                                                   posteriors( :, gaussianNumber );
#           end
#           figure( posteriorFigure )
#           subplot( floor( sqrt( numberOfClasses ) ), ...
#                    ceil( numberOfClasses / floor( sqrt( numberOfClasses ) ) ), ...
#                    classNumber )
#           showImage( posterior )
#         end
#         clear posterior
#
#         figure( costFigure )
#         subplot( 2, 1, 1 )
#         plot( historyOfEMCost( 2 : end ) )
#         title( 'EM cost' )
#         subplot(2, 1, 2 )
#         plot( historyOfCost( 2 : end ) )
#         title( 'Cost' )
#
#         figure( biasFieldFigure )
#         for contrastNumber = 1 : numberOfContrasts
#           subplot( numberOfContrasts, 2, ( contrastNumber - 1 ) * numberOfContrasts + 1 )
#           showImage( exp( downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) ) );
#           subplot( numberOfContrasts, 2, ( contrastNumber - 1 ) * numberOfContrasts + 2 )
#           downSampledBiasField = backprojectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, ...
#                                                                             biasFieldCoefficients( :, contrastNumber ) );
#           showImage( exp( downSampledBiasField ) .* downSampledMask )
#         end
#         drawnow
#
#       end % End test if we need to show some figures
#
#
#       % Check for convergence
#       % relativeChangeCost = ( historyOfEMCost(end-1) - historyOfEMCost(end) ) /  historyOfEMCost(end)
#       % if ( relativeChangeCost < stopCriterionEM )
#       changeCostPerVoxel = ( historyOfEMCost(end-1) - historyOfEMCost(end) ) / length( downSampledMaskIndices );
#       if ( changeCostPerVoxel < optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion )
#         % Converged
#         disp( 'EM converged!' )
#         break;
#       end
#
#
#
#
#       %
#       % M-step: update the model parameters based on the current posterior
#       %
#       % First the mixture model parameters
#       for gaussianNumber = 1 : numberOfGaussians
#         posterior = posteriors( :, gaussianNumber );
#
#         mean = biasCorrectedData' * posterior ./ sum( posterior );
#         tmp = biasCorrectedData - repmat( mean', [ size( biasCorrectedData, 1 ) 1 ] );
#         %variance = ( tmp' * ( tmp .* repmat( posterior, [ 1 numberOfContrasts ] ) ) + dataVariance ) ...
#         %            / ( 2 * ( numberOfContrasts + 2 ) + sum( posterior ) );
#         variance = ( tmp' * ( tmp .* repmat( posterior, [ 1 numberOfContrasts ] ) ) + ...
#                                 pseudoVarianceOfWishartPrior * numberOfPseudoMeasurementsOfWishartPrior ) ...
#                     / ( sum( posterior ) + numberOfPseudoMeasurementsOfWishartPrior );
#         if modelSpecifications.useDiagonalCovarianceMatrices
#           % Force diagonal covariance matrices
#           variance = diag( diag( variance ) );
#         end
#
#         variances( gaussianNumber, :, : ) = variance;
#         means( gaussianNumber, : ) = mean';
#
#       end
#       mixtureWeights = sum( posteriors + eps )';
#       for classNumber = 1 : numberOfClasses
#         % mixture weights are normalized (those belonging to one mixture sum to one)
#         numberOfComponents = numberOfGaussiansPerClass( classNumber );
#         gaussianNumbers = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + [ 1 : numberOfComponents ];
#
#         mixtureWeights( gaussianNumbers ) = mixtureWeights( gaussianNumbers ) / sum( mixtureWeights( gaussianNumbers ) );
#       end
#
#
#       % Now update the parameters of the bias field model.
#       %  if ( ( multiResolutionLevel == 1 ) && ( iterationNumber ~= 1 ) ) % Don't attempt bias field correction until
#       %                                                                    % decent mixture model parameters are available
#       if ( estimateBiasField && ( iterationNumber > 1 ) ) % Don't attempt bias field correction until
#                                                           % decent mixture model parameters are available
#
#
#         %
#         % Bias field correction: implements Eq. 8 in the paper
#         %
#         %    Van Leemput, "Automated Model-based Bias Field Correction of MR Images of the Brain", IEEE TMI 1999
#         %
#         precisions = zeros( size( variances ) );
#         for classNumber = 1 : numberOfGaussians
#           precisions( classNumber, :, : ) = reshape( inv( squeeze( variances( classNumber, :, : ) ) ), ...
#                                                      [ 1 numberOfContrasts numberOfContrasts ] );
#         end
#
#         lhs = zeros( prod( numberOfBasisFunctions ) * numberOfContrasts ); % left-hand side of linear system
#         rhs = zeros( prod( numberOfBasisFunctions ) * numberOfContrasts, 1 ); % right-hand side of linear system
#         weightsImageBuffer = zeros( downSampledImageSize );
#         tmpImageBuffer = zeros( downSampledImageSize );
#         for contrastNumber1 = 1 : numberOfContrasts
#           tmp = zeros( size( data, 1 ), 1 );
#           for contrastNumber2 = 1 : numberOfContrasts
#             classSpecificWeights = posteriors .* repmat( squeeze( precisions( :, contrastNumber1, contrastNumber2 ) )', ...
#                                                          [ size( posteriors, 1 ) 1 ] );
#             weights = sum( classSpecificWeights, 2 );
#
#             % Build up stuff needed for rhs
#             predicted = sum( classSpecificWeights .* repmat( means( :, contrastNumber2 )', [ size( posteriors, 1 ) 1 ] ), 2 ) ...
#                         ./ ( weights + eps );
#             residue = data( :, contrastNumber2 ) - predicted;
#             tmp = tmp + weights .* residue;
#
#             % Fill in submatrix of lhs
#             weightsImageBuffer( downSampledMaskIndices ) = weights;
#             lhs( ( contrastNumber1 - 1 ) * prod( numberOfBasisFunctions ) + [ 1 : prod( numberOfBasisFunctions ) ], ...
#                  ( contrastNumber2 - 1 ) * prod( numberOfBasisFunctions ) + [ 1 : prod( numberOfBasisFunctions ) ] ) = ...
#                   computePrecisionOfKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, weightsImageBuffer );
#
#           end % End loop over contrastNumber2
#
#           tmpImageBuffer( downSampledMaskIndices ) = tmp;
#           rhs( ( contrastNumber1 - 1 ) * prod( numberOfBasisFunctions ) + [ 1 : prod( numberOfBasisFunctions ) ] ) = ...
#                           projectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, tmpImageBuffer );
#
#         end % End loop over contrastNumber1
#
#         % lhs = lhs + diag( 0.001 * diag( lhs ) );
#
#         biasFieldCoefficients = reshape( lhs \ rhs, [ prod( numberOfBasisFunctions ) numberOfContrasts ] );
#         for contrastNumber = 1 : numberOfContrasts
#           downSampledBiasField = backprojectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
#           tmp = downSampledImageBuffers( :, :, :, contrastNumber ) - downSampledBiasField .* downSampledMask;
#           downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) = tmp;
#           biasCorrectedData( :, contrastNumber ) = tmp( downSampledMaskIndices );
#         end
#
#       end % End test if multiResolutionLevel == 1
#
#
#     end % End EM iterations
#     historyOfEMCost = historyOfEMCost( 2 : end );
#     timeTakenIntensityParameterUpdating = toc( startTimeIntensityParameterUpdating );
#     historyOfTimeTakenIntensityParameterUpdating = [ historyOfTimeTakenIntensityParameterUpdating; ...
#                                                      timeTakenIntensityParameterUpdating ];
#
#
#     %
#     % Part II: update the position of the mesh nodes for the current mixture model and bias field parameter estimates
#     %
#
#     %
#     startTimeDeformationUpdating = tic;
#
#     % Create ITK images to pass on to the mesh node position cost calculator
#     if ( exist( 'downSampledBiasCorrectedImages' ) == 1 )
#       % Clean up mess from any previous iteration
#       for contrastNumber = 1 : numberOfContrasts
#         kvlClear( downSampledBiasCorrectedImages( contrastNumber ) );
#       end
#     end
#     for contrastNumber = 1 : numberOfContrasts
#       downSampledBiasCorrectedImages( contrastNumber ) = ...
#              kvlCreateImage( single( downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) ) );
#     end
#
#     % Set up cost calculator
#     calculator = kvlGetCostAndGradientCalculator( 'AtlasMeshToIntensityImage', ...
#                                                    downSampledBiasCorrectedImages, ...
#                                                    'Sliding', ...
#                                                    transform, ...
#                                                    means, variances, mixtureWeights, numberOfGaussiansPerClass );
#
#     %optimizerType = 'ConjugateGradient';
#     optimizerType = 'L-BFGS';
#     optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
#                                     'Verbose', optimizationOptions.verbose, ...
#                                     'MaximalDeformationStopCriterion', optimizationOptions.maximalDeformationStopCriterion, ...
#                                     'LineSearchMaximalDeformationIntervalStopCriterion', ...
#                                       optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion, ...
#                                     'MaximumNumberOfIterations', optimizationOptions.maximumNumberOfDeformationIterations, ...
#                                     'BFGS-MaximumMemoryLength', optimizationOptions.BFGSMaximumMemoryLength );
#
#     historyOfDeformationCost = [];
#     historyOfMaximalDeformation = [];
#     nodePositionsBeforeDeformation = kvlGetMeshNodePositions( mesh );
#     deformationStartTime = tic;
#     while true
#       %
#       stepStartTime = tic;
#       [ minLogLikelihoodTimesDeformationPrior, maximalDeformation ] = kvlStepOptimizer( optimizer );
#       disp( [ 'maximalDeformation ' num2str( maximalDeformation ) ' took ' num2str( toc( stepStartTime ) ) ' sec' ] )
#
#       if ( maximalDeformation == 0 )
#         break;
#       end
#
#       %
#       historyOfDeformationCost = [ historyOfDeformationCost; minLogLikelihoodTimesDeformationPrior ];
#       historyOfMaximalDeformation = [ historyOfMaximalDeformation; maximalDeformation ];
#
#     end % End loop over iterations
#     kvlClear( calculator );
#     kvlClear( optimizer );
#     % haveMoved = ( length( historyOfDeformationCost ) > 0 );
#     nodePositionsAfterDeformation = kvlGetMeshNodePositions( mesh );
#     maximalDeformationApplied = sqrt( max( sum( ...
#                 ( nodePositionsAfterDeformation - nodePositionsBeforeDeformation ).^2, 2 ) ) );
#     disp( '==============================' )
#     disp( [ 'iterationNumber: ' num2str( iterationNumber ) ] )
#     disp( [ '    maximalDeformationApplied: ' num2str( maximalDeformationApplied ) ] )
#     disp( [ '  ' num2str( toc( deformationStartTime ) ) ' sec' ] )
#     disp( '==============================' )
#
#
#     % Show a little movie comparing before and after deformation so far...
#     if ( showFigures )
#       figure( deformationMovieFigure )
#       newColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, downSampledImageSize ), reducedColors );
#
#       set( deformationMovieFigure, 'position', get( 0, 'ScreenSize' ) );
#       for i = 1 : 10
#         priorVisualizationAlpha = 0.4;
#         backgroundImage = exp( downSampledImageBuffers( :, :, :, 1 ) );
#         backgroundImage = backgroundImage - min( backgroundImage(:) );
#         backgroundImage = backgroundImage / max( backgroundImage(:) );
#
#         % showImage( oldColorCodedPriors )
#         imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( backgroundImage, [ 1 1 1 3 ] ) + ...
#                       priorVisualizationAlpha * oldColorCodedPriors;
#         showImage( imageToShow )
#         drawnow
#         pause( 0.1 )
#         % showImage( newColorCodedPriors )
#         imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( backgroundImage, [ 1 1 1 3 ] ) + ...
#                       priorVisualizationAlpha * newColorCodedPriors;
#         showImage( imageToShow )
#         drawnow
#         pause( 0.1 )
#       end
#     end
#
#     % Keep track of the cost function we're optimizing
#     historyOfCost = [ historyOfCost; minLogLikelihoodTimesDeformationPrior + intensityModelParameterCost ];
#     historyOfMaximalDeformationApplied = [ historyOfMaximalDeformationApplied; maximalDeformationApplied ];
#     timeTakenDeformationUpdating = toc( startTimeDeformationUpdating );
#     historyOfTimeTakenDeformationUpdating = [ historyOfTimeTakenDeformationUpdating; ...
#                                               timeTakenDeformationUpdating ];
#
#
#     % Save something about how the estimation proceeded
#     %historyWithinEachIteration( iterationNumber ).priors = priors;
#     %historyWithinEachIteration( iterationNumber ).posteriors = posteriors;
#     historyWithinEachIteration( iterationNumber ).historyOfEMCost = historyOfEMCost;
#     historyWithinEachIteration( iterationNumber ).mixtureWeights = mixtureWeights;
#     historyWithinEachIteration( iterationNumber ).means = means;
#     historyWithinEachIteration( iterationNumber ).variances = variances;
#     historyWithinEachIteration( iterationNumber ).biasFieldCoefficients = biasFieldCoefficients;
#     historyWithinEachIteration( iterationNumber ).historyOfDeformationCost = historyOfDeformationCost;
#     historyWithinEachIteration( iterationNumber ).historyOfMaximalDeformation = historyOfMaximalDeformation;
#     historyWithinEachIteration( iterationNumber ).maximalDeformationApplied = maximalDeformationApplied;
#
#     % Determine if we should stop the overall iterations over the two set of parameters
#     %  if ( ( ~haveMoved ) || ...
#     %        ( ( ( historyOfCost( end-1 ) - historyOfCost( end ) ) / historyOfCost( end ) ) ...
#     %          < relativeCostDecreaseStopCriterion ) || ...
#     %        ( maximalDeformationApplied < maximalDeformationAppliedStopCriterion ) )
#     if ( ( ( ( historyOfCost( end-1 ) - historyOfCost( end ) ) / length( downSampledMaskIndices ) ) ...
#            < optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion ) ) % If EM converges in one iteration and mesh node optimization doesn't do anything
#
#       % Converged
#       break;
#     end
#
#
#   end % End looping over global iterations for this multiresolution level
#   historyOfCost = historyOfCost( 2 : end );
#
#   % Get the final node positions
#   finalNodePositions = kvlGetMeshNodePositions( mesh );
#
#   % Transform back in template space (i.e., undoing the affine registration
#   % that we applied), and save for later usage
#   tmp = ( totalTransformationMatrix \ [ finalNodePositions ones( numberOfNodes, 1 ) ]' )';
#   finalNodePositionsInTemplateSpace = tmp( :, 1 : 3 );
#
#
#   % Save something about how the estimation proceeded
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSamplingFactors = downSamplingFactors;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSampledImageBuffers = downSampledImageBuffers;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSampledMask = downSampledMask;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).initialNodePositions = initialNodePositions;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).finalNodePositions = finalNodePositions;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).initialNodePositionsInTemplateSpace = ...
#                                                                              initialNodePositionsInTemplateSpace;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).finalNodePositionsInTemplateSpace = ...
#                                                                              finalNodePositionsInTemplateSpace;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration = ...
#                                                                       historyWithinEachIteration;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfCost = historyOfCost;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfMaximalDeformationApplied = ...
#                                                                       historyOfMaximalDeformationApplied;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenIntensityParameterUpdating = ...
#                                                                       historyOfTimeTakenIntensityParameterUpdating;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenDeformationUpdating = ...
#                                                                       historyOfTimeTakenDeformationUpdating;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).priorsAtEnd = priors;
#   historyWithinEachMultiResolutionLevel( multiResolutionLevel ).posteriorsAtEnd = posteriors;
#
# end % End loop over multiresolution levels
#
