function [ worldToWorldTransformMatrix, transformedTemplateFileName ] = samseg_registerAtlas( imageFileName, meshCollectionFileName, templateFileName, savePath, showFigures, worldToWorldTransformMatrix, InitLTAFile )
%

% For converting from RAS to LPS. Itk/GEMS/SAMSEG uses LPS internally
RAS2LPS = diag([-1 -1 1 1]);

if ( nargin < 6 )
  worldToWorldTransformMatrix = [];
end

% Print out the input
fprintf('entering registerAtlas\n');
imageFileName
meshCollectionFileName
templateFileName
savePath
showFigures
worldToWorldTransformMatrix
InitLTAFile



% Read in image and template, as well as their coordinates in world (mm) space
[ image, imageToWorldTransform ] = kvlReadImage( imageFileName );
imageToWorldTransformMatrix = double( kvlGetTransformMatrix( imageToWorldTransform ) );

[ template, templateImageToWorldTransform ] = kvlReadImage( templateFileName );
templateImageToWorldTransformMatrix = double( kvlGetTransformMatrix( templateImageToWorldTransform ) );
[ ~, templateFileNameBase, templateFileNameExtension ] = fileparts( templateFileName );


%
if ( isempty( worldToWorldTransformMatrix ) )
  %
  % The solution is not externally (secretly) given, so we need to compute it
  %

  % Some hard-coded parameter settings first
  targetDownsampledVoxelSpacing = 3.0; % In mm
  K = 1e-7; % Mesh stiffness -- compared to normal models, the entropy cost function is normalized
            % (i.e., measures an average *per voxel*), so that this needs to be scaled down by the
            % number of voxels that are covered
  maximalDeformationStopCriterion = 0.005; % Measured in voxels
  lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion; % Doesn't seem to matter very much


  % Initialization
  if(isempty(InitLTAFile))
    initialWorldToWorldTransformMatrix = eye( 4 );
  else
    [~,InitLTA] = lta_read(InitLTAFile);
    if(InitLTA.type == 0)
      % The LTA is vox2vox
      initialWorldToWorldTransformMatrix = ...
	  RAS2LPS * InitLTA.dstmri.vox2ras0 * InitLTA.xform * inv(InitLTA.srcmri.vox2ras0) * inv(RAS2LPS);
    else
      % The LTA is ras2ras
      initialWorldToWorldTransformMatrix = RAS2LPS * InitLTA.xform * inv(RAS2LPS);      
    end
  end

  if true
    % Provide an initial (non-identity) affine transform guestimate

    % Rotation around X-axis (direction from left to right ear)
    theta = pi/180 * -10.0;
    rotationMatrix = eye( 4 );
    rotationMatrix( 2 : 3, 2 : 3 ) = [ cos( theta ) -sin(theta); sin(theta) cos( theta ) ];
    initialWorldToWorldTransformMatrix = rotationMatrix * initialWorldToWorldTransformMatrix;

    % Isotropic scaling
    scaling = 0.9 * ones( 1, 3 );
    scalingMatrix = diag( [ scaling 1 ] );

    initialWorldToWorldTransformMatrix = scalingMatrix * initialWorldToWorldTransformMatrix;

    K = K / prod( scaling );
  end
  initialImageToImageTransformMatrix = imageToWorldTransformMatrix \ ...
                    ( initialWorldToWorldTransformMatrix * templateImageToWorldTransformMatrix );

  % Figure out how much to downsample (depends on voxel size)
  voxelSpacing = sum( imageToWorldTransformMatrix( 1 : 3, 1 : 3 ).^2 ).^( 1/2 );
  downSamplingFactors = max( round( targetDownsampledVoxelSpacing ./ voxelSpacing ), [ 1 1 1 ] )
  if 1
    % Use initial transform to define the reference (rest) position of the mesh (i.e., the one
    % where the log-prior term is zero)
    meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
                                            kvlCreateTransform( initialImageToImageTransformMatrix ), ...
                                            K * prod( downSamplingFactors ) );
    mesh = kvlGetMesh( meshCollection, -1 );
  else
    % "Proper" initialization: apply the initial transform but don't let it affect the deformation
    % prior
    meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
                                            kvlCreateTransform( imageToWorldTransformMatrix \ ...
                                                                templateImageToWorldTransformMatrix ), ...
                                            K * prod( downSamplingFactors ) );
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
  imageBuffer = imageBuffer( 1 : downSamplingFactors( 1 ) : end, ...
                             1 : downSamplingFactors( 2 ) : end, ...
                             1 : downSamplingFactors( 3 ) : end );
  image = kvlCreateImage( imageBuffer );
  kvlScaleMesh( mesh, 1 ./ downSamplingFactors );
  alphas = kvlGetAlphasInMeshNodes( mesh );
  gmClassNumber = 3;  % Needed for displaying purposes
  if 0
    % Get rid of the background class
    alphas = alphas( :, 2 : end );
    kvlSetAlphasInMeshNodes( mesh, alphas )
    gmClassNumber = gmClassNumber-1;
  end
  numberOfClasses = size( alphas, 2 );
  colors = 255 * [ hsv( numberOfClasses ) ones( numberOfClasses, 1 ) ];

  if showFigures
    figure
    priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    for classNumber = 1 : numberOfClasses
      subplot( 2, 3, classNumber )
      showImage( priors( :, :, :, classNumber ) )
    end
  end

  %
  % Get a registration cost and use it to evaluate some promising starting point proposals
  calculator = kvlGetCostAndGradientCalculator( 'MutualInformation', ...
                                                image, 'Affine' );

  [ cost gradient ] = kvlEvaluateMeshPosition( calculator, mesh );
  if true
    %
    [ xtmp, ytmp, ztmp ] = ndgrid( 1 : size( imageBuffer, 1 ), ...
                                   1 : size( imageBuffer, 2 ), ...
                                   1 : size( imageBuffer, 3 ) );
    centerOfGravityImage = [ xtmp(:) ytmp(:) ztmp(:) ]' * imageBuffer(:) / sum( imageBuffer(:) );
    priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    %tmp = sum( priors, 4 );
    tmp = sum( priors( :, :, :, 2 : end ), 4 );
    centerOfGravityAtlas = [ xtmp(:) ytmp(:) ztmp(:) ]' * tmp(:) / sum( tmp(:) );
    %
    initialTranslation = double( centerOfGravityImage - centerOfGravityAtlas );
    nodePositions = kvlGetMeshNodePositions( mesh );
    trialNodePositions = nodePositions + repmat( initialTranslation', [ size( nodePositions, 1 ) 1 ] );
    kvlSetMeshNodePositions( mesh, trialNodePositions );
    [ trialCost trialGradient ] = kvlEvaluateMeshPosition( calculator, mesh );
    if ( trialCost >= cost )
      % Center of gravity was not a success; revert to what we had before
      kvlSetMeshNodePositions( mesh, nodePositions );
    else
      % This is better starting position; remember that we applied it
      initialImageToImageTransformMatrix( 1 : 3, 4 ) = ...
              initialImageToImageTransformMatrix( 1 : 3, 4 ) + diag( downSamplingFactors ) * initialTranslation;
    end

  end
  %
  originalNodePositions = kvlGetMeshNodePositions( mesh );
  % Visualize starting situation
  priorVisualizationAlpha = 0.4;
  if showFigures
    figure

    priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
    colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
    mask = ( sum( double( priors ), 4 ) / (2^16-1) ) > .5;
    subplot( 2, 2, 1 )
    showImage( imageBuffer );
    subplot( 2, 2, 2 )
    showImage( imageBuffer .* mask );
    subplot( 2, 2, 3 )
    imageToShow = ( imageBuffer - min( imageBuffer(:) ) ) / ( max( imageBuffer(:) ) - min( imageBuffer(:) ) );
    imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( imageToShow, [ 1 1 1 3 ] ) + ...
                  priorVisualizationAlpha * colorCodedPriors;
    showImage( imageToShow )
    subplot( 2, 2, 4 )
    tmp = double( priors( :, :, :, gmClassNumber ) ) / ( 2^16-1 );
    showImage( mosaicImages( tmp, imageBuffer .* mask, 2 ) );

    drawnow
  end

  % Get an optimizer, and stick the cost function into it
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
    %return
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
      colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
      mask = ( sum( double( priors ), 4 ) / (2^16-1) ) > .5;
      subplot( 2, 2, 1 )
      showImage( imageBuffer );
      subplot( 2, 2, 2 )
      showImage( imageBuffer .* mask );
      subplot( 2, 2, 3 )
      imageToShow = ( imageBuffer - min( imageBuffer(:) ) ) / ( max( imageBuffer(:) ) - min( imageBuffer(:) ) );
      imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( imageToShow, [ 1 1 1 3 ] ) + ...
                    priorVisualizationAlpha * colorCodedPriors;
      showImage( imageToShow )
      subplot( 2, 2, 4 )
      tmp = double( priors( :, :, :, gmClassNumber ) ) / ( 2^16-1 );
      showImage( mosaicImages( tmp, imageBuffer .* mask, 2 ) );
      drawnow

      % Show affine matrix, retrieved from any four non-colinear points before and after registration
      nodePositions = kvlGetMeshNodePositions( mesh );
      pointNumbers = [ 1 111 202 303 ];
      originalY = [ originalNodePositions( pointNumbers, : )'; 1 1 1 1 ];
      Y = [ nodePositions( pointNumbers, : )'; 1 1 1 1 ];
      extraImageToImageTransformMatrix = Y * inv( originalY );
      scaling = svd( extraImageToImageTransformMatrix( 1 : 3, 1 : 3 ) );
      disp( [ 'scaling: ' num2str( scaling' ) ] )
    end

  end % End loop over iterations
  numberOfIterations
  toc


  % For debugging and/or quality control purposes, save a picture of the registration result to file
  priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
  colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
  mask = ( sum( double( priors ), 4 ) / (2^16-1) ) > .5;
  overlayQcImage = ( imageBuffer - min( imageBuffer(:) ) ) / ( max( imageBuffer(:) ) - min( imageBuffer(:) ) );
  overlayQcImage = ( 1 - priorVisualizationAlpha ) * repmat( overlayQcImage, [ 1 1 1 3 ] ) + ...
                priorVisualizationAlpha * colorCodedPriors;

  tmp = double( priors( :, :, :, gmClassNumber ) ) / ( 2^16-1 );
  mosaicQcImage = mosaicImages( tmp, imageBuffer .* mask, 2 );

  overlayCollage = getCollage( overlayQcImage, 10 );
  mosaicCollage = getCollage( mosaicQcImage, 10 );

  borderSize = 20;
  DIM = [ size( overlayCollage, 1 ) size( overlayCollage, 2 ) ];
  qcFigure = zeros( [ DIM( 1 )  2 * DIM( 2 ) 3 ]  + [ 2*borderSize 3*borderSize 0 ] ) + .5;
  qcFigure( borderSize + [ 1 : DIM( 1 ) ], borderSize + [ 1 : DIM( 2 ) ], : ) = overlayCollage;
  qcFigure( borderSize + [ 1 : DIM( 1 ) ], ...
            2 * borderSize + DIM( 2 ) + [ 1 : DIM( 2 ) ], : ) = mosaicCollage;
  qcFigureFileName = fullfile( savePath, ...
                              [ templateFileNameBase '_coregistrationCqFigure.png' ] );
  imwrite( qcFigure, qcFigureFileName )


  % Retrieve the implicitly applied affine matrix from any four non-colinear points before and after registration,
  % taking into account the downsampling that we applied
  nodePositions = kvlGetMeshNodePositions( mesh );
  pointNumbers = [ 1 111 202 303 ];
  originalY = [ diag( downSamplingFactors ) * originalNodePositions( pointNumbers, : )'; 1 1 1 1 ];
  Y = [ diag( downSamplingFactors ) * nodePositions( pointNumbers, : )'; 1 1 1 1 ];
  extraImageToImageTransformMatrix = Y * inv( originalY );

  % Final result: the image-to-image (from template to image) as well as the world-to-world transform that
  % we computed (the latter would be the identity matrix if we didn't move the image at all)
  imageToImageTransformMatrix = extraImageToImageTransformMatrix * initialImageToImageTransformMatrix;
  worldToWorldTransformMatrix = imageToWorldTransformMatrix * imageToImageTransformMatrix * ...
                                inv( templateImageToWorldTransformMatrix );

else
  % The world-to-world transfrom is externally given, so let's just compute the corresponding image-to-image
  % transform (needed for subsequent computations) and be done
  imageToImageTransformMatrix = inv( imageToWorldTransformMatrix ) * worldToWorldTransformMatrix * ...
                                templateImageToWorldTransformMatrix
end % End test if the solution is externally given


% Save the image-to-image and the world-to-world affine registration matrices
transformationMatricesFileName = fullfile( savePath, ...
                                           [ templateFileNameBase '_coregistrationMatrices.mat' ] );
eval( [ 'save ' transformationMatricesFileName ' imageToImageTransformMatrix worldToWorldTransformMatrix;' ] )


% Compute the talairach.xfm
% Load fsaverage orig.mgz -- this is the ultimate target/destination
fshome = getenv('FREESURFER_HOME');
fnamedst = sprintf('%s/subjects/fsaverage/mri/orig.mgz',fshome);
fsaorig = MRIread(fnamedst,1);
% Compute the vox2vox from the template to fsaverage assuming they
%   share world RAS space
RAS2LPS = diag([-1 -1 1 1]);
M = inv(RAS2LPS*fsaorig.vox2ras)*(templateImageToWorldTransformMatrix);
% Compute the input to fsaverage vox2vox by combining the
% input-template vox2vox and the template-fsaverage vox2vox
X = M*inv(imageToImageTransformMatrix);
% Now write out the LTA. This can be used as the talairach.lta in recon-all
invol = MRIread(imageFileName,1); % have to reread to get header info
lta.type = 0;
lta.xform = X;
lta.srcfile = imageFileName;
lta.srcmri = invol;
lta.srcmri.vol = [];
lta.dstfile = fnamedst;
lta.dstmri = fsaorig;
lta.dstmri.vol = [];
lta.subject = 'fsaverage';
ltaFileName = sprintf('%s/samseg.talairach.lta',savePath);
lta_write(ltaFileName,lta);
fprintf('Done computng and writing out LTA %s\n',ltaFileName);


% For historical reasons, we applied the estimated transformation to the template; let's do that now
desiredTemplateImageToWorldTransformMatrix = imageToWorldTransformMatrix * imageToImageTransformMatrix
transformedTemplateFileName = fullfile( savePath, ...
                                        [ templateFileNameBase '_coregistered' templateFileNameExtension ] );
kvlWriteImage( template, transformedTemplateFileName, ...
               kvlCreateTransform( desiredTemplateImageToWorldTransformMatrix ) );

