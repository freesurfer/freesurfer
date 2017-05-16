%
%  imageFileName
%  templateFileName = 'mni305_masked_autoCropped.mgz';
%  meshCollectionFileName = 'CurrentMeshCollection30New.txt.gz';
%  compressionLookupTableFileName = 'namedCompressionLookupTable.txt';
%


%
downSamplingFactor = 3; % Use "1" for no downsampling
%  K = 1e-7; % Mesh stiffness -- compared to normal models, the entropy cost function is normalized 
%            % (i.e., measures an average *per voxel*), so that this needs to be scaled down by the
%            % number of voxels that are covered
tissueTypesToUse = [ 1 2 3 ]; % 1=WM, 2=GM, 3=CSF
maximalDeformationStopCriterion = 0.05; % Measured in voxels
lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion; % Doesn't seem to matter very much
%  initializeUsingCenterOfGravityAlignment = false;
%  showFigures = true;



kvlClear; % Clear all the wrapped C++ stuff
close all;

% Read in image, and figure out where the mesh is located
[ image, imageToWorldTransform ] = kvlReadImage( imageFileName );
[ template, templateImageToWorldTransform ] = kvlReadImage( templateFileName );

imageToWorldTransformMatrix = double( kvlGetTransformMatrix( imageToWorldTransform ) );
templateImageToWorldTransformMatrix = double( kvlGetTransformMatrix( templateImageToWorldTransform ) );
initialWorldToWorldTransformMatrix = eye( 4 );
if true
  % Provide an initial (non-identity) affine transform guestimate
  
  % Rotation around X-axis (direction from left to right ear)
  theta = pi/180*10.0;
  rotationMatrix = eye( 4 );
  rotationMatrix( 2 : 3, 2 : 3 ) = [ cos( theta ) -sin(theta); sin(theta) cos( theta ) ]; 
  initialWorldToWorldTransformMatrix = rotationMatrix * initialWorldToWorldTransformMatrix;
  
  % Isotropic scaling
  scaling = 0.9;
  scalingMatrix = diag( [ scaling scaling scaling 1 ] );
  initialWorldToWorldTransformMatrix = scalingMatrix * initialWorldToWorldTransformMatrix;
  
  K = K / scaling^3;
end
initialImageToImageTransformMatrix = imageToWorldTransformMatrix \ ...
                  ( initialWorldToWorldTransformMatrix * templateImageToWorldTransformMatrix );

if 1
  % Use initial transform to define the reference (rest) position of the mesh (i.e., the one
  % where the log-prior term is zero)
  meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
                                          kvlCreateTransform( initialImageToImageTransformMatrix ), ...
                                          K * downSamplingFactor^3 );
  mesh = kvlGetMesh( meshCollection, -1 );
else
  % "Proper" initialization: apply the initial transform but don't let it affect the deformation
  % prior
  meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
                                          kvlCreateTransform( imageToWorldTransformMatrix \ ...
                                                              templateImageToWorldTransformMatrix ), ...
                                          K * downSamplingFactor^3 );
  mesh = kvlGetMesh( meshCollection, -1 );
  nodePositions = kvlGetMeshNodePositions( mesh );
  tmp = [ nodePositions ones( size( nodePositions, 1 ), 1 ) ];
  tmp = ( imageToWorldTransformMatrix \ ...
          ( initialWorldToWorldTransformMatrix * imageToWorldTransformMatrix ) ) * tmp';
  nodePositions = tmp( 1:3, : )';
  kvlSetMeshNodePositions( mesh, nodePositions );

end

  
% Get image data
imageBuffer = kvlGetImageBuffer( image );
if showFigures
  figure
  showImage( imageBuffer );
end

% Downsample
if ( downSamplingFactor ~= 1 )
  imageBuffer = imageBuffer( 1 : downSamplingFactor : end, ...
                             1 : downSamplingFactor : end, ...
                             1 : downSamplingFactor : end );
  image = kvlCreateImage( imageBuffer );
  kvlScaleMesh( mesh, 1/downSamplingFactor );
end




% Compute the "reduced" alphas - those referring to the "super"-structures 
originalAlphas = kvlGetAlphasInMeshNodes( mesh );
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
sameGaussianParameters = cell(0,0);
sameGaussianParameters{1} = [ 0 ];  % Background is a separate class
sameGaussianParameters{2} = [ 2 7 16 41 46 28 60 ]; % Force all white matter structures to have the same intensity model
sameGaussianParameters{3} = [ 3 8 42 47 11 50 17 53 18 54 26 58 77 80 10 49 12 51 13 52 ]; % Same for several gray matter structures
sameGaussianParameters{4} = [ 4 5 14 15 24 43 44 72 30 62 31 63 ]; % Same for CSF
[ reducedAlphas, reducingLookupTable ] = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );
size( reducedAlphas )
max( abs( sum( reducedAlphas, 2 ) - 1 ) )  % Make sure these vectors really sum to 1

reducedColors = [ 0 255 0 255; ...
                  0 0 255 255; ...
                  255 0 0 255 ];



% Set the reduced alphas to be the alphas of the mesh
kvlSetAlphasInMeshNodes( mesh, reducedAlphas )


% Get rid of the background class
reducedAlphas = kvlGetAlphasInMeshNodes( mesh );
reducedAlphas = reducedAlphas( :, 2 : end );
kvlSetAlphasInMeshNodes( mesh, reducedAlphas )

% Get rid of some other classes
reducedAlphas = kvlGetAlphasInMeshNodes( mesh );
reducedAlphas = reducedAlphas( :, tissueTypesToUse );
kvlSetAlphasInMeshNodes( mesh, reducedAlphas )

  
if showFigures
  figure
  priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
  for reducedLabel = 1 : size( reducedAlphas, 2 )
    subplot( 2, 2, reducedLabel )
    showImage( priors( :, :, :, reducedLabel ) )
  end
end


% 
if initializeUsingCenterOfGravityAlignment
  %
  [ xtmp, ytmp, ztmp ] = ndgrid( 1 : size( imageBuffer, 1 ), ...
                                 1 : size( imageBuffer, 2 ), ...
                                 1 : size( imageBuffer, 3 ) );
  centerOfGravityImage = [ xtmp(:) ytmp(:) ztmp(:) ]' * imageBuffer(:) / sum( imageBuffer(:) );
  
  priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
  tmp = sum( priors, 4 );
  centerOfGravityAtlas = [ xtmp(:) ytmp(:) ztmp(:) ]' * tmp(:) / sum( tmp(:) ); 

  %
  initialTranslation = double( centerOfGravityImage - centerOfGravityAtlas );
  nodePositions = kvlGetMeshNodePositions( mesh );
  nodePositions = nodePositions + repmat( initialTranslation', [ size( nodePositions, 1 ) 1 ] );
  kvlSetMeshNodePositions( mesh, nodePositions );
  
  if showFigures
    figure
    priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    for reducedLabel = 1 : size( reducedAlphas, 2 )
      subplot( 2, 2, reducedLabel )
      showImage( priors( :, :, :, reducedLabel ) )
    end
  end

end




%
numberOfClasses = size( reducedAlphas, 2 );
numberOfNodes = size( reducedAlphas, 1 );
originalNodePositions = kvlGetMeshNodePositions( mesh );


% Visualize starting situation
priorVisualizationAlpha = 0.4;
if showFigures
  figure
  priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
  colorCodedPriors = kvlColorCodeProbabilityImages( priors, reducedColors );
  imageToShow = ( imageBuffer - min( imageBuffer(:) ) ) / ( max( imageBuffer(:) ) - min( imageBuffer(:) ) );
  imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( imageToShow, [ 1 1 1 3 ] ) + ...
                priorVisualizationAlpha * colorCodedPriors;
  showImage( imageToShow )
  drawnow
end


% Get a registration cost and stick into an optimizer
calculator = kvlGetCostAndGradientCalculator( 'ConditionalGaussianEntropy', ...
                                                  image, 'Affine' );
optimizerType = 'L-BFGS';
optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
                                'Verbose', 1, ...
                                'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ... 
                                'LineSearchMaximalDeformationIntervalStopCriterion', ...
                                lineSearchMaximalDeformationIntervalStopCriterion, ...
                                'BFGS-MaximumMemoryLength', 12 ); % Affine registration only has 12 DOF
                                
numberOfIterations = 0;
tic
while true
  %
  [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer )
  if ( maximalDeformation == 0 )
    break;
  end
  numberOfIterations = numberOfIterations + 1;

  %
  % Visualize progress
  %
  if showFigures
    % Show figure
    priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    colorCodedPriors = kvlColorCodeProbabilityImages( priors, reducedColors );
    imageToShow = ( imageBuffer - min( imageBuffer(:) ) ) / ( max( imageBuffer(:) ) - min( imageBuffer(:) ) );
    alpha = .5;
    imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( imageToShow, [ 1 1 1 3 ] ) + ...
                  priorVisualizationAlpha * colorCodedPriors;
    showImage( imageToShow )
    drawnow

    % Show affine matrix, retrieved from any four non-colinear points before and after registration
    nodePositions = kvlGetMeshNodePositions( mesh );
    pointNumbers = [ 1 101 201 301 ];
    originalY = [ originalNodePositions( pointNumbers, : )'; 1 1 1 1 ];
    Y = [ nodePositions( pointNumbers, : )'; 1 1 1 1 ];
    extraImageToImageTransformMatrix = Y * inv( originalY );
    scaling = svd( extraImageToImageTransformMatrix( 1 : 3, 1 : 3 ) ); 
    disp( [ 'scaling: ' num2str( scaling' ) ] )
  end  

end % End loop over iterations
numberOfIterations
toc

  
% Retrieve the implicitly applied affine matrix from any four non-colinear points before and after registration,
% taking into account the downsampling that we applied
nodePositions = kvlGetMeshNodePositions( mesh );
pointNumbers = [ 1 101 201 301 ];
originalY = [ downSamplingFactor * originalNodePositions( pointNumbers, : )'; 1 1 1 1 ];
Y = [ downSamplingFactor * nodePositions( pointNumbers, : )'; 1 1 1 1 ];
extraImageToImageTransformMatrix = Y * inv( originalY );

% Final result: the image-to-image (from template to image) as well as the world-to-world transform that
% we computed (the latter would be the identity matrix if we didn't move the image at all)
imageToImageTransformMatrix = extraImageToImageTransformMatrix * initialImageToImageTransformMatrix;
worldToWorldTransformMatrix = imageToWorldTransformMatrix * imageToImageTransformMatrix * ...
                              inv( templateImageToWorldTransformMatrix );
[ dummy, templateFileNameBase, templateFileNameExtension ] = fileparts( templateFileName );
transformationMatricesFileName = fullfile( savePath, ...
                                           [ templateFileNameBase '_coregistrationMatrices.mat' ] );
eval( [ 'save ' transformationMatricesFileName ' imageToImageTransformMatrix worldToWorldTransformMatrix;' ] )


% For historical reasons, we applied the estimated transformation to the template; let's do that now
desiredTemplateImageToWorldTransformMatrix = imageToWorldTransformMatrix * imageToImageTransformMatrix                   
transformedTemplateFileName = fullfile( savePath, ...
                                        [ templateFileNameBase '_coregistered' templateFileNameExtension ] );
kvlWriteImage( template, transformedTemplateFileName, ...
               kvlCreateTransform( desiredTemplateImageToWorldTransformMatrix ) );



