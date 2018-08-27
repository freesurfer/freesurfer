function [ FreeSurferLabels, names, volumesInCubicMm ] = samsegment( imageFileNames, transformedTemplateFileName, ...
                                                                     modelSpecifications, optimizationOptions, ...
                                                                     savePath, showFigures )
%
%



% Print input options
disp( '==========================' );
disp( 'imageFileNames: ' )
for i = 1 : length( imageFileNames )
  disp( imageFileNames{ i } )
end
fprintf( '-----\n' )

disp( 'transformedTemplateFileName: ' )
disp( transformedTemplateFileName )
fprintf( '-----\n' )

disp( 'modelSpecifications: ' )
fieldNames = fieldnames( modelSpecifications );
for i = 1 : length( fieldNames )
  fieldName = fieldNames{ i };
  disp( fieldName )
  fields = getfield( modelSpecifications, fieldName );
  for j = 1 : length( fields )
    disp( fields(j) )
  end
  disp( ' ' )
end
fprintf( '-----\n' )  

disp( 'optimizationOptions:' )
fieldNames = fieldnames( optimizationOptions );
for i = 1 : length( fieldNames )
  fieldName = fieldNames{ i };
  disp( fieldName )
  fields = getfield( optimizationOptions, fieldName );
  for j = 1 : length( fields )
    disp( fields(j) )
  end
  disp( ' ' )
end
fprintf( '-----\n' )  

disp( 'savePath: ' )
disp( savePath )
fprintf( '-----\n' )

disp( 'showFigures: ' )
disp( showFigures )
fprintf( '-----\n' )
  

% Save input variables in a "history" structure
history = struct;
history.input = struct;
history.input.imageFileNames = imageFileNames;
history.input.transformedTemplateFileName = transformedTemplateFileName;
history.input.modelSpecifications = modelSpecifications;
history.input.optimizationOptions = optimizationOptions;
history.input.savePath = savePath;
history.input.showFigures = showFigures;



%
numberOfContrasts = length( imageFileNames );

% Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
% translation, rotation, scaling, and skewing) as well - this transformation will later be used
% to initially transform the location of the atlas mesh's nodes into the coordinate system of
% the image.
%
imageBuffers = [];
for contrastNumber = 1 : numberOfContrasts
  % Get the pointers to image and the corresponding transform
  [ images( contrastNumber ), transform, nonCroppedImageSize, croppingOffset ] = ...
                                    kvlReadCroppedImage( imageFileNames{ contrastNumber }, transformedTemplateFileName ); 
  imageBuffers( :, :, :, contrastNumber ) = kvlGetImageBuffer( images( contrastNumber ) ); % Get the actual imageBuffer
end
imageSize = [ size( imageBuffers, 1 ) size( imageBuffers, 2 ) size( imageBuffers, 3 ) ];

if ( showFigures )
  for contrastNumber = 1 : numberOfContrasts
    figure
    showImage( imageBuffers( :, :, :, contrastNumber ) ); % Automatically displays middle slices in each direction
  end
end

% Also read in the voxel spacing -- this is needed since we'll be specifying bias field smoothing kernels, downsampling 
% steps etc in mm.
[ ~, imageToWorldTransform ] = kvlReadImage( imageFileNames{1} );
imageToWorldTransformMatrix = kvlGetTransformMatrix( imageToWorldTransform );
voxelSpacing = sum( imageToWorldTransformMatrix( 1 : 3, 1 : 3 ).^2 ).^( 1/2 );


% Read the atlas mesh from file, immediately applying the previously determined transform to the location
% of its nodes. Rather than one mesh, the atlas consists of a so-called "collection" of meshes: they're
% all the same except for the location of their mesh nodes. The first mesh, so-called "reference mesh"
% has index -1: it represents the average shape (location of the mesh nodes) of all the other meshes. If
% the atlas was built from N training subjects, there will be N other meshes as well, each one warping
% the atlas optimally to each of the training subjects
%
% You also need to provide a value for K, which determines the flexibility of the atlas mesh, i.e., how
% much it will typically deform. Higher values correspond to stiffer meshes.
%
meshCollection = kvlReadMeshCollection( modelSpecifications.atlasFileName, transform, modelSpecifications.K );

% Retrieve the reference mesh, i.e., the mesh representing the average shape.
mesh = kvlGetMesh( meshCollection, -1 ); 

% Get a Matlab matrix containing a copy of the probability vectors in each mesh node (size numberOfNodes x
% numberOfLabels ). 
alphas = kvlGetAlphasInMeshNodes( mesh );

% Mask away uninteresting voxels. This is done by a poor man's implementation of a dilation operation on
% a non-background class mask; followed by a cropping to the area covered by the mesh (needed because 
% otherwise there will be voxels in the data with prior probability zero of belonging to any class)
labelNumber = 0; % background label
backgroundPrior = kvlRasterizeAtlasMesh( mesh, imageSize, labelNumber ); % Volume with background probs
if 1
  % Threshold background prior at 0.5 - this helps for atlases built from imperfect (i.e., automatic) 
  % segmentations, whereas background areas don't have zero probability for non-background structures
  backgroundPrior( backgroundPrior > 2^8 ) = 2^16-1;
end
if( showFigures )
  figure
  subplot( 2, 2, 1 )
  showImage( backgroundPrior )
  subplot( 2, 2, 2 )
  showImage( mosaicImages( 2^16 - 1 - double( backgroundPrior ), double( imageBuffers(:,:,:,1) ), 10 ) )
end
smoothedBackgroundPrior = kvlSmoothImageBuffer( backgroundPrior, modelSpecifications.brainMaskingSmoothingSigma ./ voxelSpacing );
if( showFigures )
  subplot( 2, 2, 3 )
  showImage( smoothedBackgroundPrior )
end
% 65535 = 2^16 - 1. priors are stored as 16bit ints
% To put the threshold in perspective: for Gaussian smoothing with a 3D isotropic kernel with variance 
% diag( sigma^2, sigma^2, sigma^2 ) a single binary "on" voxel at distance sigma results in a value of
% 1/( sqrt(2*pi)*sigma )^3 * exp( -1/2 ).
% More generally, a single binary "on" voxel at some Eucledian distance d results in a value of
% 1/( sqrt(2*pi)*sigma )^3 * exp( -1/2*d^2/sigma^2 ). Turning this around, if we threshold this at some
% value "t", a single binary "on" voxel will cause every voxel within Eucledian distance 
%
%   d = sqrt( -2*log( t * ( sqrt(2*pi)*sigma )^3 ) * sigma^2 )
%
% of it to be included in the mask.
%
% As an example, for 1mm isotropic data, the choice of sigma=3 and t=0.01 yields ... complex value -> 
% actually a single "on" voxel will then not make any voxel survive, as the normalizing constant (achieved
% at Mahalanobis distance zero) is already < 0.01 
brainMask = ( 1 - single( smoothedBackgroundPrior ) / 65535 ) > modelSpecifications.brainMaskingThreshold;

% Crop to area covered by the mesh
areaCoveredAlphas = [ zeros( size( alphas, 1 ), 1, 'single' ) ones( size( alphas, 1 ), 1, 'single' ) ];
kvlSetAlphasInMeshNodes( mesh, areaCoveredAlphas );
areaCoveredByMesh = kvlRasterizeAtlasMesh( mesh, imageSize, 1 );
kvlSetAlphasInMeshNodes( mesh, alphas );
brainMask = brainMask & ( areaCoveredByMesh > 0 );


% Mask each of the inputs
for contrastNumber = 1 : numberOfContrasts
  imageBuffer = imageBuffers( :, :, :, contrastNumber );
  imageBuffer( find( ~brainMask ) ) = 0;
  imageBuffers( :, :, :, contrastNumber ) = imageBuffer;
  % kvlSetImageBuffer( images( contrastNumber ), imageBuffers( :, :, :, contrastNumber ) );
end

if( showFigures )
  subplot( 2, 2, 4 )
  showImage( imageBuffers( :, :, :, 1 ) )
end



% Let's prepare for the bias field correction that is part of the imaging model. It assumes
% an additive effect, whereas the MR physics indicate it's a multiplicative one - so we log
% transform the data first. In order to do so, mask out zeros from
% the images.
% This removes any voxel where any contrast has a zero value
% (messes up log)
mask = true( imageSize ); % volume of ones within the mask
for contrastNumber = 1 : numberOfContrasts
  mask = mask .* ( imageBuffers( :, :, :, contrastNumber ) > 0 );
end
maskIndices = find( mask );
for contrastNumber = 1 : numberOfContrasts
  buffer = imageBuffers( :, :, :, contrastNumber );
  buffer( maskIndices ) = log( buffer( maskIndices ) );
  buffer = buffer .* mask;
  imageBuffers( :, :, :, contrastNumber ) = buffer;
end



% Merge classes into "super-structures" that define a single Gaussian mixture model shared between the classes belonging
% to the same super-structure
FreeSurferLabels = modelSpecifications.FreeSurferLabels;
names = modelSpecifications.names;
colors = modelSpecifications.colors;
[ reducedAlphas, reducedNames, reducedFreeSurferLabels, reducedColors, translationTable ] = ...
                        kvlMergeAlphas( alphas, names, modelSpecifications.sharedGMMParameters, FreeSurferLabels, colors );


if ( showFigures )
  % Rasterizing a color-code prior can potentially be done a lot more efficient using the kvl::AtlasMeshSummaryDrawer class,
  % which first determines the color in each mesh vertex and then interpolates that (only four colors to interpolate, so
  % much more efficient). However, this requires yet another Matlab wrapper with a kvl::CompressionLookupTable to be
  % created from the colors (or read from file), which is a bit tedious so I'm leaving it at that for now...
  fprintf( 'Visualizing the atlas mesh; this takes quite some time and is only here for tutorial purposes...' )
  priors = kvlRasterizeAtlasMesh( mesh, imageSize ); % Without specifying a specific label, will rasterize all simultaneously, rasterize everything
  rgbBuffer = kvlColorCodeProbabilityImages( priors, colors );
  figure
  showImage( rgbBuffer )
  clear priors rgbBuffer
  fprintf( 'done\n' )
  drawnow;
end




% The fact that we merge several neuroanatomical structures into "super"-structures for the purpose of model
% parameter estimaton, but at the same time represent each of these super-structures with a mixture of Gaussians,
% creates something of a messy situation when implementing this stuff. To avoid confusion, let's define a few
% conventions that we'll closely follow in the code as follows:
%
%   - classNumber = 1 ... numberOfClasses  -> indexes a specific super-structure (there are numberOfClasses superstructures)
%   - numberOfGaussiansPerClass            -> a numberOfClasses-dimensional vector that indicates the number of components
%                                             in the Gaussian mixture model associated with each class 
%   - gaussianNumber = 1 .... numberOfGaussians  -> indexes a specific Gaussian distribution; there are 
%                                                   numberOfGaussians = sum( numberOfGaussiansPerClass ) of those in total
%
% In certain situations it will be easier to index things using "classNumber" (which is the only relevant thing when 
% estimating mesh deformations) than using "gaussianNumber" (which is the only relevant thing when estimating mixture 
% model parameters) and vice versa, so we need to have a standardized way of going back and forth between those. For a 
% given classNumber, there are numberOfComponents = numberOfGaussiansPerClass( classNumber ) components in its mixture 
% model. By convention we convert the pair ( classNumber, componentNumber ) into gaussianNumber as follows:
%
%    gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
%
numberOfGaussiansPerClass = [ modelSpecifications.sharedGMMParameters.numberOfComponents ];
numberOfClasses = length( numberOfGaussiansPerClass );
numberOfGaussians = sum( numberOfGaussiansPerClass );



% Our bias model is a linear combination of a set of basis functions. We are using so-called
% "DCT-II" basis functions, i.e., the lowest few frequency components of the Discrete Cosine
% Transform.
%
kroneckerProductBasisFunctions = cell( 0, 0 );
numberOfBasisFunctions = zeros( 1, 3 );
for dimensionNumber = 1 : 3
  N = imageSize( dimensionNumber ); % Number of data points
  delta =  modelSpecifications.biasFieldSmoothingKernelSize / voxelSpacing( dimensionNumber ); % Measured in voxels
  M = ceil( N / delta ) + 1; % Number of frequencies to use
  Nvirtual = ( M - 1 ) * delta; % Virtual (possibly non-integer) number of data points to make sure 
                                % we reach the targeted smoothing kernel size
  js = ( [ 0 : N-1 ]' + 1/2 ) * pi / Nvirtual;
  A = cos( js * [ 0 : M-1 ] ) * sqrt( 2 / Nvirtual );
  A( :, 1 ) = A( :, 1 ) / sqrt( 2 );
  
  if showFigures
    % Show smoothing kernel
    figure
    smootherMatrix = A * ( ( A' * A ) \ A' );
    subplot( 2, 2, 1 )
    imshow( smootherMatrix, [] )
    subplotCounter = 2;
    for rowNumber = round( [ N/4 N/2 3*N/4 ] )
      subplot( 2, 2, subplotCounter )
      plot( smootherMatrix( rowNumber, : ) )
      set( gca, 'xlim', [ 1 N ] )
      grid
      hold on
      ylim = get( gca, 'ylim' );
      line( ( rowNumber ) * [ 1 1 ], ylim, 'color', 'k' )
      line( ( rowNumber + delta ) * [ 1 1 ], ylim, 'color', 'r', 'linestyle', '--' )
      line( ( rowNumber - delta ) * [ 1 1 ], ylim, 'color', 'r', 'linestyle', '--' )
      subplotCounter = subplotCounter + 1;
    end
  end
    
  kroneckerProductBasisFunctions{ dimensionNumber } = A;
  numberOfBasisFunctions( dimensionNumber ) = M;    
end
biasFieldCoefficients = zeros( prod( numberOfBasisFunctions ), numberOfContrasts ); % No bias field to start with


if ( showFigures )
  posteriorFigure = figure;
  costFigure = figure;
  deformationMovieFigure = figure;
  biasFieldFigure = figure;
end


% We do the optimization in a multi-resolution type of scheme, where large
% deformations are quickly found using smoothed versions of the atlas mesh, and the fine
% details are then found on gradually less smoothed versions until the original atlas mesh is used for the optimization.
% Multi-resolution is a standard procedure in image registration: the initial
% blurring avoids getting stuck in the first local optimum of the search space, and get the rough major
% deformations right instead.
numberOfMultiResolutionLevels = length( optimizationOptions.multiResolutionSpecification );
historyWithinEachMultiResolutionLevel = struct( [] );
for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
    
  %  
  maximumNumberOfIterations = optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).maximumNumberOfIterations;
  estimateBiasField = optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).estimateBiasField; 
  historyOfCost = [ 1/eps ];
  historyOfMaximalDeformationApplied = [];
  historyOfTimeTakenIntensityParameterUpdating = [];
  historyOfTimeTakenDeformationUpdating = [];
  fprintf('maximumNumberOfIterations %d\n',maximumNumberOfIterations);

  
  % Downsample the images, the mask, the mesh, and the bias field basis functions
  % Must be integer
  downSamplingFactors = max( round( optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).targetDownsampledVoxelSpacing ...
                                    ./ voxelSpacing ), [ 1 1 1 ] )
  downSampledMask = mask(  1 : downSamplingFactors( 1 ) : end, ...
                           1 : downSamplingFactors( 2 ) : end, ...
                           1 : downSamplingFactors( 3 ) : end );
  downSampledMaskIndices = find( downSampledMask );
  downSampledImageBuffers = [];
  for contrastNumber = 1 : numberOfContrasts
    % if ( multiResolutionLevel == numberOfMultiResolutionLevels )
    if true
      % No image smoothing
      downSampledImageBuffers( :, :, :, contrastNumber ) = imageBuffers( 1 : downSamplingFactors( 1 ) : end, ...
                                                                         1 : downSamplingFactors( 2 ) : end, ...
                                                                         1 : downSamplingFactors( 3 ) : end, ...
                                                                         contrastNumber );
    else
      % Try image smoothing
      buffer = imageBuffers( 1 : downSamplingFactors( 1 ) : end, ...
                             1 : downSamplingFactors( 2 ) : end, ...
                             1 : downSamplingFactors( 3 ) : end, ...
                             contrastNumber );
      smoothingSigmas = downSamplingFactors / 2 / sqrt( 2 * log( 2 ) ); % Variance chosen to approximately 
                                                                        % match normalized binomial filter 
                                                                        % (1/4, 1/2, 1/4) for downsampling 
                                                                        % factor of 2
      smoothingSigmas( find( downSamplingFactors == 1 ) ) = 0.0;
      smoothedBuffer = kvlSmoothImageBuffer( single( buffer ), smoothingSigmas );
      smoothedMask = kvlSmoothImageBuffer( single( downSampledMask ), smoothingSigmas );     
      downSampledImageBuffers( :, :, :, contrastNumber ) = downSampledMask .* ( smoothedBuffer ./ ( smoothedMask + eps ) );                                                                     
    end
    
  end

  % 
  downSampledKroneckerProductBasisFunctions = cell( 0, 0 );
  for dimensionNumber = 1 : 3
    A = kroneckerProductBasisFunctions{ dimensionNumber };
    downSampledKroneckerProductBasisFunctions{ dimensionNumber } = A( 1 : downSamplingFactors( dimensionNumber ) : end, : );
  end
  downSampledImageSize = size( downSampledImageBuffers( :, :, :, 1 ) );
  
  
  % Read the atlas mesh to be used for this multi-resolution level, taking into account the downsampling to position it 
  % correctly
  downSamplingTransformMatrix = diag( [ 1./downSamplingFactors 1 ] );
  totalTransformationMatrix = downSamplingTransformMatrix * double( kvlGetTransformMatrix( transform ) );
  meshCollection = ...
        kvlReadMeshCollection( optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).atlasFileName, ...
                                kvlCreateTransform( totalTransformationMatrix ), modelSpecifications.K );
  mesh = kvlGetMesh( meshCollection, -1 );
  
  % Get the initial mesh node positions, also transforming them back into template space 
  % (i.e., undoing the affine registration that we applied) for later usage 
  initialNodePositions = kvlGetMeshNodePositions( mesh );  
  numberOfNodes = size( initialNodePositions, 1 );
  tmp = ( totalTransformationMatrix \ [ initialNodePositions ones( numberOfNodes, 1 ) ]' )';
  initialNodePositionsInTemplateSpace = tmp( :, 1 : 3 );


  % If this is not the first multi-resolution level, apply the warp computed during the previous level
  if ( multiResolutionLevel > 1 )
    % Get the warp in template space
    nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel = ...
            historyWithinEachMultiResolutionLevel( multiResolutionLevel-1 ).finalNodePositionsInTemplateSpace - ...
            historyWithinEachMultiResolutionLevel( multiResolutionLevel-1 ).initialNodePositionsInTemplateSpace;
    initialNodeDeformationInTemplateSpace = kvlWarpMesh( ...
                  optimizationOptions.multiResolutionSpecification( multiResolutionLevel-1 ).atlasFileName, ...
                  nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel, ...
                  optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).atlasFileName );

    % Apply this warp on the mesh node positions in template space, and transform into current space  
    desiredNodePositionsInTemplateSpace = initialNodePositionsInTemplateSpace + initialNodeDeformationInTemplateSpace;
    tmp = ( totalTransformationMatrix * ...
            [ desiredNodePositionsInTemplateSpace ones( numberOfNodes, 1 ) ]' )';
    desiredNodePositions = tmp( :, 1 : 3 );

    %
    kvlSetMeshNodePositions( mesh, desiredNodePositions );

  end    

    
  
  % Set priors in mesh to the reduced (super-structure) ones
  alphas = kvlGetAlphasInMeshNodes( mesh );
  reducedAlphas = kvlMergeAlphas( alphas, names, modelSpecifications.sharedGMMParameters, FreeSurferLabels, colors );
  kvlSetAlphasInMeshNodes( mesh, reducedAlphas )

  
    
  
  % Algorithm-wise, we're just estimating sets of parameters for one given data (MR scan) that is
  % known and fixed throughout. However, in terms of bias field correction it will be computationally
  % more efficient to pre-compute the bias field corrected version of the scan ("corrected" with
  % the current estimate of the bias field) once and pass that on to different routines instead of the
  % original data. 
  % For convenience (although potentially a recipe for future bug introduction), I'm also keeping a
  % vectorized form of that around -- this will be useful in various places in the EM-parts. So 
  % effectively I have two redundant variables "downSampledBiasCorrectedImageBuffers" and "biasCorrectedData"
  % that really just encode the variable "biasFieldCoefficients" and so need to be meticiously updated each time 
  % "biasFieldCoefficients" is updated (!)
  downSampledBiasCorrectedImageBuffers = zeros( [ downSampledImageSize numberOfContrasts ] );
  biasCorrectedData = zeros( [ length( downSampledMaskIndices ) numberOfContrasts ] );
  for contrastNumber = 1 : numberOfContrasts
    downSampledBiasField = backprojectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
    tmp = downSampledImageBuffers( :, :, :, contrastNumber ) - downSampledBiasField .* downSampledMask;
    downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) = tmp;
    biasCorrectedData( :, contrastNumber ) = tmp( downSampledMaskIndices );
  end
  
  
  % Compute a color coded version of the atlas prior in the atlas's current pose, i.e., *before*
  % we start deforming. We'll use this just for visualization purposes
  if ( showFigures )
    oldColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, downSampledImageSize ), reducedColors );
  end

  
  
  historyWithinEachIteration = struct( [] );
  priors = zeros( length( downSampledMaskIndices ), numberOfClasses );
  posteriors = zeros( length( downSampledMaskIndices ), numberOfGaussians ); % Gaussian mixture models burst out into
                                                                             % individual Gaussian components

  % Easier to work with vector notation in the EM computations
  % reshape into a matrix
  data = zeros( [ length( downSampledMaskIndices ) numberOfContrasts ] );  
  for contrastNumber = 1:numberOfContrasts
    tmp = reshape( downSampledImageBuffers( :, :, :, contrastNumber ), [ prod(downSampledImageSize) 1 ] ); 
    data( :, contrastNumber ) = tmp( downSampledMaskIndices );
  end

  % Main iteration loop over both EM and deformation
  for iterationNumber = 1 : maximumNumberOfIterations
    
    %
    startTimeIntensityParameterUpdating = tic;
    
    %
    % Part I: estimate Gaussian mixture model parameters, as well as bias field parameters using EM.
    %
    
    % Get the priors at the current mesh position
    tmp = reshape( kvlRasterizeAtlasMesh( mesh, downSampledImageSize ), [ prod( downSampledImageSize ) numberOfClasses ] );
    priors( : ) = double( tmp( downSampledMaskIndices, : ) ) / 65535;
    
    if ( iterationNumber == 1 )
      historyWithinEachMultiResolutionLevel( multiResolutionLevel ).priorsAtStart = priors;                                                                   
    end
    

    % Start EM iterations.
    if ( ( multiResolutionLevel == 1 ) && ( iterationNumber == 1 ) )

      % Initialize the mixture parameters if this is the first time ever you run this
      means = zeros( numberOfGaussians, numberOfContrasts );
      variances = zeros( numberOfGaussians, numberOfContrasts, numberOfContrasts );
      mixtureWeights = zeros( numberOfGaussians, 1 );
      for classNumber = 1 : numberOfClasses
      
        % Calculate the global weighted mean and variance of this class, where the weights are given by the prior
        prior = priors( :, classNumber );
        mean = data' * prior / sum( prior );
        tmp = data - repmat( mean', [ size( data, 1 ) 1 ] );
        variance = tmp' * ( tmp .* repmat( prior, [ 1 numberOfContrasts ] ) ) / sum( prior );
        if modelSpecifications.useDiagonalCovarianceMatrices
          % Force diagonal covariance matrices
          variance = diag( diag( variance ) );
        end
       
       
        % Based on this, initialize the mean and variance of the individual Gaussian components in this class' 
        % mixture model: variances are simply copied from the global class variance, whereas the means are
        % determined by splitting the [ mean-sqrt( variance ) mean+sqrt( variance ) ] domain into equal intervals, 
        % the middle of which are taken to be the means of the Gaussians. Mixture weights are initialized to be 
        % all equal.
        %
        % This actually creates a mixture model that mimics the single Gaussian quite OK-ish: to visualize this
        % do e.g.,
        %
        %  for numberOfComponents = 1 : 7
        %    intervalSize = 2 / numberOfComponents;
        %    means = -1 + intervalSize/2 + [ 0 : numberOfComponents-1 ] * intervalSize;
        %    figure
        %    x = [ -6 : .1 : 6 ];
        %    gauss = exp( -x.^2/2 );
        %    plot( gauss )
        %    hold on 
        %    mixture = zeros( size( x ) );
        %    for i = 1 : numberOfComponents
        %      gauss = exp( -( x - means( i ) ).^2/2 );
        %      plot( gauss / numberOfComponents, 'g' )
        %      mixture = mixture + gauss / numberOfComponents;
        %    end
        %    plot( mixture, 'r' )
        %    grid
        %    title( [ num2str( numberOfComponents ) ' components' ] )
        %  end          
        %
        numberOfComponents = numberOfGaussiansPerClass( classNumber );
        for componentNumber = 1 : numberOfComponents
          gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
        
          variances( gaussianNumber, :, : ) = variance;
          intervalSize = 2 * sqrt( diag( variance ) ) / numberOfComponents;
          means( gaussianNumber, : ) = ( mean - sqrt( diag( variance ) ) + intervalSize/2 + ( componentNumber - 1 ) * intervalSize )';
          mixtureWeights( gaussianNumber ) = 1 / numberOfComponents;
        end        
        
      end % End loop over classes
      
      
      % Also remember the overall data variance for later usage in a conjugate prior on the variances      
      dataMean = sum( data )' / size( data, 1 );
      tmp = data - repmat( dataMean', [ size( data, 1 ) 1 ] );
      dataVariance = diag( diag( tmp' * tmp ) ) / size( data, 1 );
      numberOfPseudoMeasurementsOfWishartPrior = 1; % In Oula's code this was effectively 2 * ( numberOfContrasts + 2 )
                                                    % although I have no clue why
      pseudoVarianceOfWishartPrior = dataVariance / numberOfPseudoMeasurementsOfWishartPrior;
      
    end % End test need for initialization

    % stopCriterionEM = 1e-5;
    historyOfEMCost = [ 1/eps ];
    for EMIterationNumber = 1 : 100
      %
      % E-step: compute the posteriors based on the current parameters.
      %
      for classNumber = 1 : numberOfClasses
        prior = priors( :, classNumber );
        
        numberOfComponents = numberOfGaussiansPerClass( classNumber );
        for componentNumber = 1 : numberOfComponents
          gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
          
          mean = means( gaussianNumber, : )';
          variance = squeeze( variances( gaussianNumber, :, : ) );
          L = chol( variance, 'lower' );  % variance = L * L'
          tmp = L \ ( biasCorrectedData' - repmat( mean, [ 1 size( biasCorrectedData, 1 ) ] ) );
          squaredMahalanobisDistances = ( sum( tmp.^2, 1 ) )';
          sqrtDeterminantOfVariance = prod( diag( L ) ); % Same as sqrt( det( variance ) )
          gaussianLikelihoods = exp( -squaredMahalanobisDistances / 2 ) / ( 2 * pi )^( numberOfContrasts / 2 ) / sqrtDeterminantOfVariance;
          
          posteriors( :, gaussianNumber ) = gaussianLikelihoods .* ( mixtureWeights( gaussianNumber ) * prior );
          
        end % End loop over mixture components

      end % End loop over classes
      normalizer = sum( posteriors, 2 ) + eps;
      if 0
        x = zeros( downSampledImageSize );
        x( downSampledMaskIndices ) = -log( normalizer );
        figure
        showImage( x )
      end
      posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfGaussians ] );
      minLogLikelihood = -sum( log( normalizer ) );
      intensityModelParameterCost = 0;
      useRestrictedGMMs = false;
      if ~useRestrictedGMMs
        for gaussianNumber = 1 : numberOfGaussians
          variance = squeeze( variances( gaussianNumber, :, : ) );
          
          % Evaluate unnormalized Wishart distribution (conjugate prior on precisions) with parameters 
          %
          %   scale matrix V = inv( pseudoVarianceOfWishartPrior * numberOfPseudoMeasurementsOfWishartPrior )
          %
          % and 
          %
          %   degrees of freedom n = numberOfPseudoMeasurementsOfWishartPrior + numberOfContrasts + 1
          %
          % which has pseudoVarianceOfWishartPrior as the MAP solution in the absence of any data
          %
          minLogUnnormalizedWishart = ...
              trace( variance \ pseudoVarianceOfWishartPrior ) * numberOfPseudoMeasurementsOfWishartPrior / 2 + ...
              numberOfPseudoMeasurementsOfWishartPrior / 2 * log( det( variance ) );
          intensityModelParameterCost = intensityModelParameterCost + minLogUnnormalizedWishart;
        end
      end  
      historyOfEMCost = [ historyOfEMCost; minLogLikelihood + intensityModelParameterCost ];

      % Show some figures
      if ( showFigures )
        for classNumber = 1 : numberOfClasses
          posterior = zeros( downSampledImageSize );
          numberOfComponents = numberOfGaussiansPerClass( classNumber );
          for componentNumber = 1 : numberOfComponents
            gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
            posterior( downSampledMaskIndices ) = posterior( downSampledMaskIndices ) + ...
                                                  posteriors( :, gaussianNumber );
          end                                        
          figure( posteriorFigure )
          subplot( floor( sqrt( numberOfClasses ) ), ...
                   ceil( numberOfClasses / floor( sqrt( numberOfClasses ) ) ), ...
                   classNumber )
          showImage( posterior )
        end
        clear posterior
        
        figure( costFigure )
        subplot( 2, 1, 1 )
        plot( historyOfEMCost( 2 : end ) )
        title( 'EM cost' )
        subplot(2, 1, 2 )
        plot( historyOfCost( 2 : end ) )
        title( 'Cost' )
        
        figure( biasFieldFigure )
        for contrastNumber = 1 : numberOfContrasts
          subplot( numberOfContrasts, 2, ( contrastNumber - 1 ) * numberOfContrasts + 1 )
          showImage( exp( downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) ) );
          subplot( numberOfContrasts, 2, ( contrastNumber - 1 ) * numberOfContrasts + 2 )
          downSampledBiasField = backprojectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, ...
                                                                            biasFieldCoefficients( :, contrastNumber ) );
          showImage( exp( downSampledBiasField ) .* downSampledMask )
        end
        drawnow

      end % End test if we need to show some figures

      
      % Check for convergence
      % relativeChangeCost = ( historyOfEMCost(end-1) - historyOfEMCost(end) ) /  historyOfEMCost(end)
      % if ( relativeChangeCost < stopCriterionEM )
      changeCostPerVoxel = ( historyOfEMCost(end-1) - historyOfEMCost(end) ) / length( downSampledMaskIndices );
      if ( changeCostPerVoxel < optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion )
        % Converged
        disp( 'EM converged!' )
        break;
      end 
      

           
           
      %
      % M-step: update the model parameters based on the current posterior 
      %
      % First the mixture model parameters
      if ~useRestrictedGMMs
        for gaussianNumber = 1 : numberOfGaussians
          posterior = posteriors( :, gaussianNumber );

          mean = biasCorrectedData' * posterior ./ sum( posterior );
          tmp = biasCorrectedData - repmat( mean', [ size( biasCorrectedData, 1 ) 1 ] );
          %variance = ( tmp' * ( tmp .* repmat( posterior, [ 1 numberOfContrasts ] ) ) + dataVariance ) ...
          %            / ( 2 * ( numberOfContrasts + 2 ) + sum( posterior ) );
          variance = ( tmp' * ( tmp .* repmat( posterior, [ 1 numberOfContrasts ] ) ) + ...
                                  pseudoVarianceOfWishartPrior * numberOfPseudoMeasurementsOfWishartPrior ) ...
                      / ( sum( posterior ) + numberOfPseudoMeasurementsOfWishartPrior );
          if modelSpecifications.useDiagonalCovarianceMatrices
            % Force diagonal covariance matrices
            variance = diag( diag( variance ) );
          end

          variances( gaussianNumber, :, : ) = variance;
          means( gaussianNumber, : ) = mean';
        end
      else
        %
        for classNumber = 1 : numberOfClasses
          numberOfComponents = numberOfGaussiansPerClass( classNumber );
          gaussianNumbers = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + [ 1 : numberOfComponents ];

          if 0
            save myxxx.mat biasCorrectedData posteriors gaussianNumbers      
          end
          [ mixtureMeans, mixtureSigma, ~ ] = kvlFitRestrictedGMM( biasCorrectedData, posteriors( :, gaussianNumbers ), 2 );
          
          variances( gaussianNumbers, :, : ) = diag( mixtureSigma.^2 );
          if ( ~modelSpecifications.useDiagonalCovarianceMatrices & ( length( mixtureSigma ) > 1 ) )
            % Issue an error if non-diagonal covariances are asked for 
            error( 'Non-diagonal covariance matrices for multicontrast models not implemented for restricted GMMs' )            
          end
          means( gaussianNumbers, : ) = mixtureMeans;
        end % End loop over classNumber
        
      end

      mixtureWeights = sum( posteriors + eps )';
      for classNumber = 1 : numberOfClasses
        % mixture weights are normalized (those belonging to one mixture sum to one)
        numberOfComponents = numberOfGaussiansPerClass( classNumber );
        gaussianNumbers = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + [ 1 : numberOfComponents ];

        mixtureWeights( gaussianNumbers ) = mixtureWeights( gaussianNumbers ) / sum( mixtureWeights( gaussianNumbers ) );
      end      
      

      % Now update the parameters of the bias field model. 
      %  if ( ( multiResolutionLevel == 1 ) && ( iterationNumber ~= 1 ) ) % Don't attempt bias field correction until 
      %                                                                    % decent mixture model parameters are available  
      if ( estimateBiasField && ( iterationNumber > 1 ) ) % Don't attempt bias field correction until 
                                                          % decent mixture model parameters are available  

        
        %
        % Bias field correction: implements Eq. 8 in the paper
        %
        %    Van Leemput, "Automated Model-based Bias Field Correction of MR Images of the Brain", IEEE TMI 1999
        %                       
        precisions = zeros( size( variances ) );
        for classNumber = 1 : numberOfGaussians
          precisions( classNumber, :, : ) = reshape( inv( squeeze( variances( classNumber, :, : ) ) ), ...
                                                     [ 1 numberOfContrasts numberOfContrasts ] );
        end  
          
        lhs = zeros( prod( numberOfBasisFunctions ) * numberOfContrasts ); % left-hand side of linear system
        rhs = zeros( prod( numberOfBasisFunctions ) * numberOfContrasts, 1 ); % right-hand side of linear system
        weightsImageBuffer = zeros( downSampledImageSize );
        tmpImageBuffer = zeros( downSampledImageSize );
        for contrastNumber1 = 1 : numberOfContrasts
          tmp = zeros( size( data, 1 ), 1 );
          for contrastNumber2 = 1 : numberOfContrasts
            classSpecificWeights = posteriors .* repmat( squeeze( precisions( :, contrastNumber1, contrastNumber2 ) )', ...
                                                         [ size( posteriors, 1 ) 1 ] );
            weights = sum( classSpecificWeights, 2 );
            
            % Build up stuff needed for rhs
            predicted = sum( classSpecificWeights .* repmat( means( :, contrastNumber2 )', [ size( posteriors, 1 ) 1 ] ), 2 ) ...
                        ./ ( weights + eps );
            residue = data( :, contrastNumber2 ) - predicted;
            tmp = tmp + weights .* residue;
            
            % Fill in submatrix of lhs
            weightsImageBuffer( downSampledMaskIndices ) = weights;
            lhs( ( contrastNumber1 - 1 ) * prod( numberOfBasisFunctions ) + [ 1 : prod( numberOfBasisFunctions ) ], ...
                 ( contrastNumber2 - 1 ) * prod( numberOfBasisFunctions ) + [ 1 : prod( numberOfBasisFunctions ) ] ) = ... 
                  computePrecisionOfKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, weightsImageBuffer );

          end % End loop over contrastNumber2
          
          tmpImageBuffer( downSampledMaskIndices ) = tmp;
          rhs( ( contrastNumber1 - 1 ) * prod( numberOfBasisFunctions ) + [ 1 : prod( numberOfBasisFunctions ) ] ) = ...
                          projectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, tmpImageBuffer );
                          
        end % End loop over contrastNumber1
        
        % lhs = lhs + diag( 0.001 * diag( lhs ) );
        
        biasFieldCoefficients = reshape( lhs \ rhs, [ prod( numberOfBasisFunctions ) numberOfContrasts ] );
        for contrastNumber = 1 : numberOfContrasts
          downSampledBiasField = backprojectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
          tmp = downSampledImageBuffers( :, :, :, contrastNumber ) - downSampledBiasField .* downSampledMask;
          downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) = tmp;
          biasCorrectedData( :, contrastNumber ) = tmp( downSampledMaskIndices );
        end

      end % End test if multiResolutionLevel == 1
    
      
    end % End EM iterations
    historyOfEMCost = historyOfEMCost( 2 : end );
    timeTakenIntensityParameterUpdating = toc( startTimeIntensityParameterUpdating );  
    historyOfTimeTakenIntensityParameterUpdating = [ historyOfTimeTakenIntensityParameterUpdating; ...
                                                     timeTakenIntensityParameterUpdating ];

                                                     
    %
    % Part II: update the position of the mesh nodes for the current mixture model and bias field parameter estimates
    %
    
    %
    startTimeDeformationUpdating = tic;

    % Create ITK images to pass on to the mesh node position cost calculator
    if ( exist( 'downSampledBiasCorrectedImages' ) == 1 )
      % Clean up mess from any previous iteration
      for contrastNumber = 1 : numberOfContrasts
        kvlClear( downSampledBiasCorrectedImages( contrastNumber ) );
      end
    end
    for contrastNumber = 1 : numberOfContrasts
      downSampledBiasCorrectedImages( contrastNumber ) = ...
             kvlCreateImage( single( downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) ) );      
    end

    % Set up cost calculator
    calculator = kvlGetCostAndGradientCalculator( 'AtlasMeshToIntensityImage', ...
                                                   downSampledBiasCorrectedImages, ...
                                                   'Sliding', ...
                                                   transform, ...
                                                   means, variances, mixtureWeights, numberOfGaussiansPerClass );

    %optimizerType = 'ConjugateGradient';
    optimizerType = 'L-BFGS';
    optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
                                    'Verbose', optimizationOptions.verbose, ...
                                    'MaximalDeformationStopCriterion', optimizationOptions.maximalDeformationStopCriterion, ... 
                                    'LineSearchMaximalDeformationIntervalStopCriterion', ...
                                      optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion, ...
                                    'MaximumNumberOfIterations', optimizationOptions.maximumNumberOfDeformationIterations, ... 
                                    'BFGS-MaximumMemoryLength', optimizationOptions.BFGSMaximumMemoryLength );

    historyOfDeformationCost = [];
    historyOfMaximalDeformation = [];
    nodePositionsBeforeDeformation = kvlGetMeshNodePositions( mesh );
    deformationStartTime = tic;
    while true
      %
      stepStartTime = tic;
      [ minLogLikelihoodTimesDeformationPrior, maximalDeformation ] = kvlStepOptimizer( optimizer );
      disp( [ 'maximalDeformation ' num2str( maximalDeformation ) ' took ' num2str( toc( stepStartTime ) ) ' sec' ] )
      
      if ( maximalDeformation == 0 )
        break;
      end

      %
      historyOfDeformationCost = [ historyOfDeformationCost; minLogLikelihoodTimesDeformationPrior ];
      historyOfMaximalDeformation = [ historyOfMaximalDeformation; maximalDeformation ];
    
    end % End loop over iterations
    kvlClear( calculator );
    kvlClear( optimizer );
    % haveMoved = ( length( historyOfDeformationCost ) > 0 );
    nodePositionsAfterDeformation = kvlGetMeshNodePositions( mesh );
    maximalDeformationApplied = sqrt( max( sum( ...
                ( nodePositionsAfterDeformation - nodePositionsBeforeDeformation ).^2, 2 ) ) );
    disp( '==============================' )
    disp( [ 'iterationNumber: ' num2str( iterationNumber ) ] )
    disp( [ '    maximalDeformationApplied: ' num2str( maximalDeformationApplied ) ] )
    disp( [ '  ' num2str( toc( deformationStartTime ) ) ' sec' ] )
    disp( '==============================' )

    
    % Show a little movie comparing before and after deformation so far...
    if ( showFigures )
      figure( deformationMovieFigure )
      newColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, downSampledImageSize ), reducedColors );
      
      set( deformationMovieFigure, 'position', get( 0, 'ScreenSize' ) );
      for i = 1 : 10
        priorVisualizationAlpha = 0.4;
        backgroundImage = exp( downSampledImageBuffers( :, :, :, 1 ) );
        backgroundImage = backgroundImage - min( backgroundImage(:) );
        backgroundImage = backgroundImage / max( backgroundImage(:) );

        % showImage( oldColorCodedPriors )
        imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( backgroundImage, [ 1 1 1 3 ] ) + ...
                      priorVisualizationAlpha * oldColorCodedPriors;
        showImage( imageToShow )
        drawnow
        pause( 0.1 )
        % showImage( newColorCodedPriors )
        imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( backgroundImage, [ 1 1 1 3 ] ) + ...
                      priorVisualizationAlpha * newColorCodedPriors;
        showImage( imageToShow )
        drawnow
        pause( 0.1 )
      end
    end
    
    % Keep track of the cost function we're optimizing
    historyOfCost = [ historyOfCost; minLogLikelihoodTimesDeformationPrior + intensityModelParameterCost ];
    historyOfMaximalDeformationApplied = [ historyOfMaximalDeformationApplied; maximalDeformationApplied ];
    timeTakenDeformationUpdating = toc( startTimeDeformationUpdating );  
    historyOfTimeTakenDeformationUpdating = [ historyOfTimeTakenDeformationUpdating; ...
                                              timeTakenDeformationUpdating ];
    
    
    % Save something about how the estimation proceeded
    %historyWithinEachIteration( iterationNumber ).priors = priors;
    %historyWithinEachIteration( iterationNumber ).posteriors = posteriors;
    historyWithinEachIteration( iterationNumber ).historyOfEMCost = historyOfEMCost;
    historyWithinEachIteration( iterationNumber ).mixtureWeights = mixtureWeights;
    historyWithinEachIteration( iterationNumber ).means = means;
    historyWithinEachIteration( iterationNumber ).variances = variances;
    historyWithinEachIteration( iterationNumber ).biasFieldCoefficients = biasFieldCoefficients;
    historyWithinEachIteration( iterationNumber ).historyOfDeformationCost = historyOfDeformationCost;
    historyWithinEachIteration( iterationNumber ).historyOfMaximalDeformation = historyOfMaximalDeformation;
    historyWithinEachIteration( iterationNumber ).maximalDeformationApplied = maximalDeformationApplied;

    % Determine if we should stop the overall iterations over the two set of parameters
    %  if ( ( ~haveMoved ) || ...
    %        ( ( ( historyOfCost( end-1 ) - historyOfCost( end ) ) / historyOfCost( end ) ) ...
    %          < relativeCostDecreaseStopCriterion ) || ...
    %        ( maximalDeformationApplied < maximalDeformationAppliedStopCriterion ) )
    if ( ( ( ( historyOfCost( end-1 ) - historyOfCost( end ) ) / length( downSampledMaskIndices ) ) ...
           < optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion ) ) % If EM converges in one iteration and mesh node optimization doesn't do anything
         
      % Converged
      break;
    end
    
    
  end % End looping over global iterations for this multiresolution level
  historyOfCost = historyOfCost( 2 : end );
    
  % Get the final node positions 
  finalNodePositions = kvlGetMeshNodePositions( mesh );
  
  % Transform back in template space (i.e., undoing the affine registration
  % that we applied), and save for later usage 
  tmp = ( totalTransformationMatrix \ [ finalNodePositions ones( numberOfNodes, 1 ) ]' )';
  finalNodePositionsInTemplateSpace = tmp( :, 1 : 3 );
  
  
  % Save something about how the estimation proceeded
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSamplingFactors = downSamplingFactors;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSampledImageBuffers = downSampledImageBuffers;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSampledMask = downSampledMask;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).initialNodePositions = initialNodePositions;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).finalNodePositions = finalNodePositions;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).initialNodePositionsInTemplateSpace = ...
                                                                             initialNodePositionsInTemplateSpace;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).finalNodePositionsInTemplateSpace = ...
                                                                             finalNodePositionsInTemplateSpace;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration = ...
                                                                      historyWithinEachIteration;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfCost = historyOfCost;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfMaximalDeformationApplied = ...
                                                                      historyOfMaximalDeformationApplied;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenIntensityParameterUpdating = ...
                                                                      historyOfTimeTakenIntensityParameterUpdating;  
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenDeformationUpdating = ...
                                                                      historyOfTimeTakenDeformationUpdating;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).priorsAtEnd = priors;                                                                   
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).posteriorsAtEnd = posteriors;
    
end % End loop over multiresolution levels


% Save something about how the estimation proceeded
history.imageBuffers = imageBuffers;
history.mask = mask;
history.historyWithinEachMultiResolutionLevel = historyWithinEachMultiResolutionLevel;
eval( [ 'save ' savePath '/history.mat history -v7.3' ] );


% OK, now that all the parameters have been estimated, try to segment the original, full resolution image
% with all the original labels (instead of the reduced "super"-structure labels we created).

% Get bias field corrected images
biasCorrectedImageBuffers = zeros( [ imageSize numberOfContrasts ] );
biasFields = zeros( [ imageSize numberOfContrasts ] );
for contrastNumber = 1 : numberOfContrasts
  biasField = backprojectKroneckerProductBasisFunctions( kroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
  biasCorrectedImageBuffers( :, :, :, contrastNumber ) = ...
              imageBuffers( :, :, :, contrastNumber ) - biasField .* mask;
  biasFields( :, :, :, contrastNumber ) = biasField;
end



% Read the atlas, applying the affine registration transform
meshCollection = kvlReadMeshCollection( modelSpecifications.atlasFileName, transform, modelSpecifications.K );
mesh = kvlGetMesh( meshCollection, -1 );

% Get the mesh node positions transformed back into template space (i.e., undoing the affine registration that we applied)
nodePositions = kvlGetMeshNodePositions( mesh );  
numberOfNodes = size( nodePositions, 1 );
transformMatrix = double( kvlGetTransformMatrix( transform ) );
tmp = ( transformMatrix \ [ nodePositions ones( numberOfNodes, 1 ) ]' )';
nodePositionsInTemplateSpace = tmp( :, 1 : 3 );

% Get the estimated warp in template space
estimatedNodeDeformationInTemplateSpace = ...
          kvlWarpMesh( optimizationOptions.multiResolutionSpecification( end ).atlasFileName, ...
                        historyWithinEachMultiResolutionLevel( end ).finalNodePositionsInTemplateSpace ...
                        - historyWithinEachMultiResolutionLevel( end ).initialNodePositionsInTemplateSpace, ...
                        modelSpecifications.atlasFileName );

% Apply this warp on the mesh node positions in template space, and transform into current space  
desiredNodePositionsInTemplateSpace = nodePositionsInTemplateSpace + estimatedNodeDeformationInTemplateSpace;
tmp = ( transformMatrix * [ desiredNodePositionsInTemplateSpace ones( numberOfNodes, 1 ) ]' )';
desiredNodePositions = tmp( :, 1 : 3 );

%
kvlSetMeshNodePositions( mesh, desiredNodePositions );

%
alphas = kvlGetAlphasInMeshNodes( mesh );
numberOfStructures = size( alphas, 2 );



% Get the priors as dictated by the current mesh position
data = reshape( biasCorrectedImageBuffers, [ prod( imageSize ) numberOfContrasts ] ); 
priors = kvlRasterizeAtlasMesh( mesh, imageSize );
priors = reshape( priors, [ prod( imageSize ) numberOfStructures ] );  


% Ignore everything that's has zero intensity
priors = priors( maskIndices, : );
data = data( maskIndices, : );


% Calculate the posteriors
posteriors = zeros( size( priors ), 'double' );
for structureNumber = 1 : numberOfStructures

  prior = single( priors( :, structureNumber ) ) / 65535;
  
  mixedLikelihoods = zeros( length( maskIndices ), 1 );
  for classNumber = 1 : numberOfClasses
    %
    fraction = translationTable( classNumber, structureNumber );
    if ( fraction < 1e-10 )
      continue;
    end 
  
    % Compute likelihood of this class (aka mixture model)
    likelihoods = zeros( length( maskIndices ), 1 );
    numberOfComponents = numberOfGaussiansPerClass( classNumber );
    for componentNumber = 1 : numberOfComponents
      gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;

      mean = means( gaussianNumber, : )';
      variance = squeeze( variances( gaussianNumber, :, : ) );
      mixtureWeight = mixtureWeights( gaussianNumber );
      
      L = chol( variance, 'lower' );  % variance = L * L'
      tmp = L \ ( data' - repmat( mean, [ 1 size( data, 1 ) ] ) );
      squaredMahalanobisDistances = ( sum( tmp.^2, 1 ) )';
      sqrtDeterminantOfVariance = prod( diag( L ) ); % Same as sqrt( det( variance ) )
      gaussianLikelihoods = exp( -squaredMahalanobisDistances / 2 ) / ( 2 * pi )^( numberOfContrasts / 2 ) / sqrtDeterminantOfVariance;

      likelihoods = likelihoods + gaussianLikelihoods * mixtureWeight;
    end

    % 
    mixedLikelihoods = mixedLikelihoods + likelihoods * fraction;
  end 
  
  posteriors( :, structureNumber ) = mixedLikelihoods .* prior;

end % End loop over structures

normalizer = sum( posteriors, 2 ) + eps;
posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfStructures ] );


% Compute volumes in mm^3
volumeOfOneVoxel = abs( det( imageToWorldTransformMatrix( 1:3, 1:3 ) ) );
volumesInCubicMm = ( sum( posteriors ) )' * volumeOfOneVoxel;


% Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention 
[ ~, structureNumbers ] = max( posteriors, [], 2 );
freeSurferSegmentation = zeros( imageSize, 'uint16' );
for structureNumber = 1 : numberOfStructures
  freeSurferSegmentation( maskIndices( find( structureNumbers == structureNumber ) ) ) = FreeSurferLabels( structureNumber );
end

% Write to file, remembering to un-crop the segmentation to the original image size
uncroppedFreeSurferSegmentation = zeros( nonCroppedImageSize, 'single' );
uncroppedFreeSurferSegmentation( croppingOffset( 1 ) + [ 1 : imageSize( 1 ) ], ...
                                 croppingOffset( 2 ) + [ 1 : imageSize( 2 ) ], ...
                                 croppingOffset( 3 ) + [ 1 : imageSize( 3 ) ] ) = freeSurferSegmentation;
fprintf( 'Writing out freesurfer segmentation\n' );
kvlWriteImage( kvlCreateImage( uncroppedFreeSurferSegmentation ), ...
               fullfile( savePath, 'crispSegmentation.nii' ), ...
               imageToWorldTransform );


% Also write out the bias field and the bias corrected image, each time remembering to un-crop the images 
for contrastNumber = 1 : numberOfContrasts
  [ dataPath, scanName, ext ] = fileparts( imageFileNames{ contrastNumber } );
  
  % First bias field
  biasField = zeros( nonCroppedImageSize, 'single' );
  biasField( croppingOffset( 1 ) + [ 1 : imageSize( 1 ) ], ...
             croppingOffset( 2 ) + [ 1 : imageSize( 2 ) ], ...
             croppingOffset( 3 ) + [ 1 : imageSize( 3 ) ] ) = exp( biasFields( :, :, :, contrastNumber ) ) .* mask;
  outputFileName = fullfile( savePath, [ scanName '_biasField.nii' ] );
  kvlWriteImage( kvlCreateImage( biasField ), outputFileName, imageToWorldTransform );


  % Then bias field corrected data
  biasCorrected = zeros( nonCroppedImageSize, 'single' );
  biasCorrected( croppingOffset( 1 ) + [ 1 : imageSize( 1 ) ], ...
                 croppingOffset( 2 ) + [ 1 : imageSize( 2 ) ], ...
                 croppingOffset( 3 ) + [ 1 : imageSize( 3 ) ] ) = exp( biasCorrectedImageBuffers( :, :, :, contrastNumber ) );
  outputFileName = fullfile( savePath, [ scanName '_biasCorrected.nii' ] );
  kvlWriteImage( kvlCreateImage( biasCorrected ), outputFileName, imageToWorldTransform );

end


