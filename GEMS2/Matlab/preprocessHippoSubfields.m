%
% This script will initialize the position of our atlas mesh nodes to something that more or less matches the outer boundaries
% of the whole hippocampus as determined by FreeSurfer. It's really only meant to provide a parameter intialization for the
% optimization of a parametric model, and you could certainly experiment with leaving it out altogether...
%
% You'll notice that there are not much comments here - that's because the real code is in the file "processHippoSubfields.m"
% where I've also tried to document every little thing I do. So once you understand how things work there, please come back
% here to see what's going (if you decide this type of hackish preprocessing is useful).
%
% There are some section with more comments - these are parts that have no equivalent in "processHippoSubfields.m" but for
% which I still wanted to provide you with documentation.
%


% Clean up the Matlab work space
kvlClear % Clear all the wrapped C++ stuff
close all
clear all


% 
asegFileName = 'aseg.mgz';  % FreeSurfer's volumetric segmentation results. This are non-probabilistic, "crisp" definitions
boundingFileName = 'imageDump_coregistered.mgz'; % Bounding box
meshCollectionFileName = 'CurrentMeshCollection30.gz'; % The tetrahedral atlas mesh
compressionLookupTableFileName = 'compressionLookupTable.txt'; % Look-up table belonging to the atlas


%
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );


%
[ aseg, transform ] = kvlReadCroppedImage( asegFileName, boundingFileName );
asegBuffer = kvlGetImageBuffer( aseg );
asegSize = size( asegBuffer );
figure
showImage( aseg )


%
meshCollection = kvlReadMeshCollection( meshCollectionFileName );
kvlTransformMeshCollection( meshCollection, transform );
K = 0.01;
kvlSetKOfMeshCollection( meshCollection, K );


% Retrieve the reference mesh, i.e., the mesh representing the average shape.
mesh = kvlGetMesh( meshCollection, -1 );
originalNodePositions = kvlGetMeshNodePositions( mesh );


% Just for illustrative purposes, let's also display this mesh warped onto each
% of the training subjects, as computed during the group-wise registration during
% the atlas building
for meshNumber = 0 : 9  % C-style indexing
  pause( .1 )
  showImage( kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( kvlGetMesh( meshCollection, meshNumber ), asegSize ), colors ) )
end



% In this first step, we're cheating in that we're going to segment a "fake" intensity image obtained by
% giving FreeSurfer's ASEG hippocampus segmentation an artificial, high intensity, and all the rest an 
% artificial, low (but not zero, for reasons that will soon become clear) intensity. Let's generate this
% artificial image here.
% 
% First, we look up what intensity value corresponds to hippocampus in FreeSurfer's ASEG result, and then
% we use that information
hippoLabel = FreeSurferLabels( find( strcmp( 'Right-Hippocampus', cellstr( names ) ) ) );
hippoIndices = find( asegBuffer == hippoLabel );
cheatingImageBuffer = ones( asegSize, 'single' );
cheatingImageBuffer( hippoIndices ) = 255;
cheatingImage = kvlCreateImage( cheatingImageBuffer );
figure
showImage( cheatingImage, [], [ 0 5 ] )  % Using dynamic range for display that shows what's going on



%
sameGaussianParameters = cell(0,0);
sameGaussianParameters{1} = [];
for FreeSurferLabel = { 'CSF', ...
                        'Right-Cerebral-White-Matter', ...
                        'Right-Cerebral-Cortex' }
  sameGaussianParameters{1} = [ sameGaussianParameters{1} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
sameGaussianParameters{2} = [];
for FreeSurferLabel = { 'Right-Hippocampus', ...
                        'right_CA2-3', ...
                        'right_CA1', ...
                        'right_fimbria', ...
                        'right_presubiculum', ...
                        'right_hippocampal_fissure', ...
                        'right_CA4-DG', ...
                        'right_subiculum' }
  sameGaussianParameters{2} = [ sameGaussianParameters{2} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end



%
originalAlphas = kvlGetAlphasInMeshNodes( mesh );
cheatingAlphas = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );
kvlSetAlphasInMeshNodes( mesh, cheatingAlphas )
figure
priors = kvlRasterizeAtlasMesh( mesh, asegSize );
for cheatingLabel = 1 : size( cheatingAlphas, 2 )
  subplot( 1, 2, cheatingLabel )
  showImage( priors( :, :, :, cheatingLabel ) )
end


%
meshSmoothingSigma = 3.0; 
fprintf( 'Smoothing mesh with kernel size %f ...', meshSmoothingSigma )
kvlSmoothMesh( mesh, meshSmoothingSigma )
fprintf( 'done\n' )
figure
priors = kvlRasterizeAtlasMesh( mesh, asegSize );
for cheatingLabel = 1 : size( cheatingAlphas, 2 )
  subplot( 1, 2, cheatingLabel )
  showImage( priors( :, :, :, cheatingLabel ) )
end


%
cheatingOptimizer = kvlGetLevenbergMarquardtOptimizer( mesh, cheatingImage, transform );
cheatingMeans = [ 1 255 ]';  % We know the mean intensity for each of the two classes in the image since we created it
                             % accordingly - no need to estimate these.
cheatingVariances = ( [ 1 1 ]' ).^2; % I guess much smaller values could be used, but anyway the variance is already
                                     % tiny compared to the distance between the means, so that we're really optimizing
                                     % the warp that best explains class labels rather than intensities

figure
subplot( 2, 2, 1 )
oldPrior = kvlRasterizeAtlasMesh( mesh, asegSize, 1 );
showImage( oldPrior );
subplot( 2, 2, 2 )
showImage( cheatingImage )
  
kvlSetOptimizerProperties( cheatingOptimizer, cheatingMeans, cheatingVariances );
historyOfMinLogLikelihoodTimesPrior = [ 1/eps ];
for positionUpdatingIterationNumber = 1 : 200
  % Calculate a good step. The first one is very slow because of various set-up issues
  tic
  [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlDeformOneStep( cheatingOptimizer );
  elapsedTime = toc;
  disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
  minLogLikelihoodTimesPrior
  historyOfMinLogLikelihoodTimesPrior = [ historyOfMinLogLikelihoodTimesPrior; minLogLikelihoodTimesPrior ];
    
  % Test if we need to stop
  if ( ( maximalDeformation <= 0.05 ) | ...
        ( ( historyOfMinLogLikelihoodTimesPrior( end-1 ) - historyOfMinLogLikelihoodTimesPrior( end ) ) ...
          / historyOfMinLogLikelihoodTimesPrior( end ) < 1e-5 ) )
    break;
  end

  % Show what we have
  subplot( 2, 2, 3 )
  showImage( kvlRasterizeAtlasMesh( mesh, asegSize, 1 ) );
  subplot( 2, 2, 4 )
  plot( historyOfMinLogLikelihoodTimesPrior( 2 : end ) )
  drawnow
end
kvlClear( cheatingOptimizer )


% OK, we're done. let's modify the mesh atlas in such a way that our computed mesh node positions are
% assigned to what was originally the mesh warp corresponding to the first training subject.
kvlSetAlphasInMeshNodes( mesh, originalAlphas )
updatedNodePositions = kvlGetMeshNodePositions( mesh );
kvlSetMeshCollectionPositions( meshCollection, ... %
                               originalNodePositions, ... % reference position (average "shape")
                               updatedNodePositions );


% Compare the average shape we started with, with the shape we have computed now in a little movie
figure
originalPositionColorCodedPriors = ...
      kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( kvlGetMesh( meshCollection, -1 ), asegSize ), colors );
updatedPositionColorCodedPriors = ...
      kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( kvlGetMesh( meshCollection, 0 ), asegSize ), colors );
for i=1:10
  showImage( originalPositionColorCodedPriors )
  drawnow
  pause( 0.1 )
  showImage( updatedPositionColorCodedPriors )
  drawnow
  pause( 0.1 )
end


% Write the resulting atlas mesh to file
transformMatrix = kvlGetTransformMatrix( transform );
inverseTransform = kvlCreateTransform( inv( transformMatrix ) );
kvlTransformMeshCollection( meshCollection, inverseTransform );
kvlWriteMeshCollection( meshCollection, 'warpedOriginalMesh.txt' );



