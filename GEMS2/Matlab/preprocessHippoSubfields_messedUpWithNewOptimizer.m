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
figure
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
if 0
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
end



%
cheatingMeans = [ 1 255 ]';  % We know the mean intensity for each of the two classes in the image since we created it
                             % accordingly - no need to estimate these.
cheatingVariances = ( [ 1 1 ]' ).^2; % I guess much smaller values could be used, but anyway the variance is already
                                     % tiny compared to the distance between the means, so that we're really optimizing
                                     % the warp that best explains class labels rather than intensities

if 0
  % Check repeatability of multi-threaded gradient computation. The contributions of different threads are
  % added, and depending on the order of summation and the exact tetrahedra assigned to the different threads
  % the results will be slightly different
  kvlSetMaximumNumberOfThreads( 16 )

  x = originalNodePositions(:);
  costs = zeros( 10, 1 );
  gradients = zeros( length(x), 10 );
  for i=1:10
    [ cost gradient ] = kvlEvaluateMeshPositionInVectorFormat( x, mesh, cheatingImage, transform, cheatingMeans, cheatingVariances );
    costs( i ) = cost;
    gradients( :, i ) = gradient;
  end

  meanGradient = sum( gradients, 2 ) / size( gradients, 2 );
  [ maxDiff maxInd ] = max( abs( repmat( meanGradient, [ 1 10 ] ) - gradients ) ./ ...
                            ( repmat( abs( meanGradient ), [ 1 10 ] ) + eps ) * 100 );
  gradients( maxInd( 1 ), : )
  gradients( maxInd( 1 ), 2 : end ) - gradients( maxInd( 1 ), 1 )



    kvlSetMaximumNumberOfThreads( 1 )
  [ cost gradient ] = kvlEvaluateMeshPositionInVectorFormat( x, mesh, cheatingImage, transform, cheatingMeans, cheatingVariances );
  [ maxDiff maxInd ] = max( abs( repmat( gradient, [ 1 10 ] ) - gradients ) ./ ...
                            ( repmat( abs( gradient ), [ 1 10 ] ) + eps ) * 100 );
  gradients( maxInd( 1 ), : )
  gradients( maxInd( 1 ), : ) - gradient( maxInd( 1 ) )
  return
end


if 0

% Let's optimize the mesh node locations using Mark Schmidt's minFunc
costFunctionHandle = @(x)kvlEvaluateMeshPositionInVectorFormat( x, mesh, cheatingImage, transform, cheatingMeans, cheatingVariances );

if 0
  addpath /home/koen/software/minFunc

  options = [];
  options.display = 'full'; % [ off | final | (iter) | full | excessive ]
  %  options.maxFunEvals = 2000; % Maximum number of function evaluations allowed (1000). It's about 0.1 sec for one evaluation
  %  options.MaxIter = 2000; % Maximum number of iterations allowed (500)
  options.Method = 'cg';
  options.LS = 0; % 0...2 recommended for objective functions that return NaN (Armijo Bactracking instead of Wolfe stuff)
                  %   - 0: Backtrack w/ Step Size Halving
                  %   - 1: Backtrack w/ Quadratic/Cubic Interpolation from new function values
                  %   - 2: Backtrack w/ Cubic Interpolation from new function + gradient
                  %   values (default for 'bb' and 'sd')
                  %
                  % 1 doens't work; simple 0 works best
  %  options.Method = 'scg';
  %  options.Method = 'newton0';
  %  options.Method = 'lbfgs';
  %  options.Method = 'pnewton0';
  %  reshapedPosition = minFunc( costFunctionHandle, reshapedInitialPosition, options );
  tic
  [ reshapedPosition cost exitflag output ] = minFunc( costFunctionHandle, originalNodePositions(:), options );
  elapsedTime = toc;

  output

  figure
  subplot( 2, 1, 1 )
  costMinFunc = output.trace.fval;
  timeMinFunc = [ 0 : length( costMinFunc )-1 ] /  ( length( costMinFunc ) - 1 ) * elapsedTime;
  plot( timeMinFunc, costMinFunc )
  title( 'cost' )
  subplot( 2, 1, 2 )
  plot( timeMinFunc, output.trace.funcCount )
  title( 'number of function evaluations' )
  disp( [ 'elapsedTime: ' num2str( elapsedTime ) ] );
  updatedNodePositions = reshape( reshapedPosition, [ length( reshapedPosition ) / 3 3 ] );

else

  tic
  [ reshapedPosition cost exitflag output ] = kvlMinFunc( costFunctionHandle, originalNodePositions(:) );
  elapsedTime = toc;

  output

  figure
  subplot( 2, 1, 1 )
  costMinFunc = output.trace.fval;
  timeMinFunc = [ 0 : length( costMinFunc )-1 ] /  ( length( costMinFunc ) - 1 ) * elapsedTime;
  plot( timeMinFunc, costMinFunc )
  title( 'cost' )
  subplot( 2, 1, 2 )
  plot( timeMinFunc, output.trace.funcCount )
  title( 'number of function evaluations' )
  disp( [ 'elapsedTime: ' num2str( elapsedTime ) ] );
  updatedNodePositions = reshape( reshapedPosition, [ length( reshapedPosition ) / 3 3 ] );

end


end


%kvlSetMaximumNumberOfThreads( 1 );

% For comparision, let's also use our own LM-based optimizer and see how it fares
kvlSetMeshNodePositions( mesh, originalNodePositions );
if 0
  cheatingOptimizer = kvlGetLevenbergMarquardtOptimizer( mesh, cheatingImage, transform );
  %maximalDeformationStopCriterion = 0.05;
  maximalDeformationStopCriterion = 0;
  relativeChangeInCostStopCriterion = 1e-5;
  maximumNumberOfPositionUpdatingIterations = 200;
else
  cheatingOptimizer = kvlGetConjugateGradientOptimizer( mesh, cheatingImage, transform );
  %maximalDeformationStopCriterion = 0;
  %relativeChangeInCostStopCriterion = 1e-9;
  maximalDeformationStopCriterion = 0.00001;
  relativeChangeInCostStopCriterion = 0;
  maximumNumberOfPositionUpdatingIterations = 500;
  %cheatingVariances = cheatingVariances / ( 2 * pi )
end


figure
subplot( 2, 2, 1 )
oldPrior = kvlRasterizeAtlasMesh( mesh, asegSize, 1 );
showImage( oldPrior );
subplot( 2, 2, 2 )
showImage( cheatingImage )
  
kvlSetOptimizerProperties( cheatingOptimizer, cheatingMeans, cheatingVariances );
historyOfMinLogLikelihoodTimesPrior = [ 1/eps ];
tStart = tic;
for positionUpdatingIterationNumber = 1 : maximumNumberOfPositionUpdatingIterations
  % Calculate a good step. The first one is very slow because of various set-up issues
  tic
  [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlDeformOneStep( cheatingOptimizer );
  elapsedTime = toc;
  disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
  minLogLikelihoodTimesPrior
  historyOfMinLogLikelihoodTimesPrior = [ historyOfMinLogLikelihoodTimesPrior; minLogLikelihoodTimesPrior ];
  %pause

  % Test if we need to stop
  if ( ( maximalDeformation <= maximalDeformationStopCriterion ) | ...
        ( ( historyOfMinLogLikelihoodTimesPrior( end-1 ) - historyOfMinLogLikelihoodTimesPrior( end ) ) ...
          / historyOfMinLogLikelihoodTimesPrior( end ) < relativeChangeInCostStopCriterion ) )
    break;
  end

  % Show what we have
  if 0
    subplot( 2, 2, 3 )
    showImage( kvlRasterizeAtlasMesh( mesh, asegSize, 1 ) );
    subplot( 2, 2, 4 )
    plot( historyOfMinLogLikelihoodTimesPrior( 2 : end ) )
    grid
    drawnow
  end
end
subplot( 2, 2, 3 )
showImage( kvlRasterizeAtlasMesh( mesh, asegSize, 1 ) );
subplot( 2, 2, 4 )
plot( historyOfMinLogLikelihoodTimesPrior( 2 : end ) )
grid
drawnow
kvlClear( cheatingOptimizer )
updatedNodePositions = kvlGetMeshNodePositions( mesh );
elapsedTime = toc( tStart );
disp( [ 'elapsedTime: ' num2str( elapsedTime ) ] );

% Take a snapshot
snapshotFileName = [ 'snaphot' num2str( tic ) '.png' ];
tmp = getframe( gcf );
imwrite( tmp.cdata, snapshotFileName );

% Save output
saveNumber = 1;
saveSessionFileName = [ 'saveSession_' sprintf( '%02d', saveNumber ) '.mat' ];
while exist( saveSessionFileName )
  saveNumber = saveNumber + 1;
  saveSessionFileName = [ 'saveSession_' sprintf( '%02d', saveNumber ) '.mat' ];
end
historyOfMinLogLikelihoodTimesPriorName = [ 'historyOfMinLogLikelihoodTimesPrior_' sprintf( '%02d', saveNumber ) ];
updatedNodePositionsName = [ 'updatedNodePositions_' sprintf( '%02d', saveNumber ) ];
elapsedTimeName = [ 'elapsedTime_' sprintf( '%02d', saveNumber ) ];
eval( [ historyOfMinLogLikelihoodTimesPriorName ' = historyOfMinLogLikelihoodTimesPrior;' ] );
eval( [ updatedNodePositionsName ' = updatedNodePositions;' ] );
eval( [ elapsedTimeName ' = elapsedTime;' ] );
eval( [ 'save ' saveSessionFileName ' ' historyOfMinLogLikelihoodTimesPriorName ' ' updatedNodePositionsName ' ' elapsedTimeName ] )



if 0
  load saveSession_01.mat            
  load saveSession_02.mat

  %  for i=2:5
  %    i
  %    historyOfMinLogLikelihoodTimesPrior_01( i )
  %    historyOfMinLogLikelihoodTimesPrior_02( i )
  %  
  %    historyOfMinLogLikelihoodTimesPrior_01( i ) - historyOfMinLogLikelihoodTimesPrior_02( i )
  %    pause
  %  end

  cost_01 = historyOfMinLogLikelihoodTimesPrior_01( 2 : end );
  time_01 = [ 0 : length( cost_01 )-1 ] /  ( length( cost_01 ) - 1 ) * elapsedTime_01;
  cost_02 = historyOfMinLogLikelihoodTimesPrior_02( 2 : end );
  time_02 = [ 0 : length( cost_02 )-1 ] /  ( length( cost_02 ) - 1 ) * elapsedTime_02;

  figure
  plot( time_01, cost_01 )
  hold on
  plot( time_02, cost_02, 'r' )
  legend( '01', '02' )


  

end


if 0
  figure
  plot( timeMinFunc, costMinFunc )
  costOwn = historyOfMinLogLikelihoodTimesPrior( 2 : end );
  timeOwn = [ 0 : length( costOwn )-1 ] /  ( length( costOwn ) - 1 ) * elapsedTime;
  hold on
  plot( timeOwn, costOwn, 'r' )
  legend( 'minFunc', 'own' )
  grid
end


return


% OK, we're done. let's modify the mesh atlas in such a way that our computed mesh node positions are
% assigned to what was originally the mesh warp corresponding to the first training subject.
kvlSetAlphasInMeshNodes( mesh, originalAlphas )
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



