function [ targetDeformation, averageDistance, maximumDistance ] = kvlWarpMesh( sourceMeshCollectionFileName, sourceDeformation, targetMeshCollectionFileName, showFigures )
%
% Applies sourceDeformation to the reference node positions of the sourceMeshCollection given by 
% sourceMeshCollectionFileName, and computes a targetDeformation on the reference node positions of the 
% targetMeshCollection given by targetMeshCollectionFileName that mimicks the same overall dense 3D deformation.
% 


if ( nargin == 0 )
  % 
  % Test run this function
  %

  %
  if 1
    sourceMeshCollectionFileName = ...
        '/data/tmp/tmpBuckner/tmp3/scratch/koenLogDir10SubjectAtlas3/CurrentMeshCollection30_multires6.gz';
  else
    % Test special case of matching mesh connectivity
    sourceMeshCollectionFileName = ...
        '/data/tmp/tmpBuckner/tmp3/scratch/koenLogDir10SubjectAtlas3/CurrentMeshCollection30.gz';
  end
  targetMeshCollectionFileName = ...
      '/data/tmp/tmpBuckner/tmp3/scratch/koenLogDir10SubjectAtlas3/CurrentMeshCollection30.gz';

      
  %
  sourceMeshCollection = kvlReadMeshCollection( sourceMeshCollectionFileName );
  sourceReferenceMesh = kvlGetMesh( sourceMeshCollection, -1 ); 
  sourceReferencePosition = kvlGetMeshNodePositions( sourceReferenceMesh );
  sourceDeformedMesh = kvlGetMesh( sourceMeshCollection, 8 ); % Purposefully deformed
  sourceDeformedPosition = kvlGetMeshNodePositions( sourceDeformedMesh );
  sourceDeformation = sourceDeformedPosition - sourceReferencePosition;
   
  %
  showFigures = true;
  targetDeformation = kvlWarpMesh( sourceMeshCollectionFileName, sourceDeformation, ...
                                   targetMeshCollectionFileName, showFigures );
                                
  % Write deformed target mesh collection to file
  targetMeshCollection = kvlReadMeshCollection( targetMeshCollectionFileName );
  targetReferenceMesh = kvlGetMesh( targetMeshCollection, -1 ); 
  targetReferencePosition = kvlGetMeshNodePositions( targetReferenceMesh );
  targetPosition = targetReferencePosition + targetDeformation;
  kvlSetMeshCollectionPositions( targetMeshCollection, targetReferencePosition, targetPosition );
  [a b c ] = fileparts( targetMeshCollectionFileName );
  resultFileName = fullfile( a, [ b '_debugDeformed' c ] );
  kvlWriteMeshCollection( targetMeshCollection, resultFileName );
  disp( [ 'Wrote result in file: ' resultFileName ] )
                                
  return;  
end

if ( nargin < 4 )
  showFigures = false;
end
      

% Set flexibility of mesh deformation      
K = .001;
      
      
%
sourceMeshCollection = kvlReadMeshCollection( sourceMeshCollectionFileName );
sourceReferenceMesh = kvlGetMesh( sourceMeshCollection, -1 ); 
sourceReferencePosition = kvlGetMeshNodePositions( sourceReferenceMesh );
sourceNumberOfNodes = size( sourceReferencePosition, 1 );

%
targetMeshCollection = kvlReadMeshCollection( targetMeshCollectionFileName, kvlCreateTransform( eye( 4 ) ), K );
targetReferenceMesh = kvlGetMesh( targetMeshCollection, -1 ); 
targetReferencePosition = kvlGetMeshNodePositions( targetReferenceMesh );
targetNumberOfNodes = size( targetReferencePosition, 1 );

% In the special case of identical mesh connectivity, no need to do any optimization
if ( targetNumberOfNodes == sourceNumberOfNodes )
  if ( max( abs( sourceReferencePosition(:) - targetReferencePosition(:) ) ) < 1e-2 )
    % The reference meshes seem to be the same - therefore simply copy the deformation
    targetDeformation = sourceDeformation;
    averageDistance = 0.0;
    maximumDistance = 0.0;
    
    return
  end
end


%
imageSize = max( sourceReferencePosition ) + 1;


% Rasterize the deformation by abusing alpha drawer
deformation = sourceDeformation;
denseDeformation = zeros( [ imageSize 3 ] );
if ( max( abs( deformation(:) ) ) > 0 )
  %
  maxDeformation = max( deformation(:) );
  minDeformation = min( deformation(:) );

  kvlSetAlphasInMeshNodes( sourceReferenceMesh, single( deformation - minDeformation ) ./ ...
                                                repmat( maxDeformation - minDeformation, [ sourceNumberOfNodes 3 ] ) );
  % denseDeformation = double( kvlRasterizeAtlasMesh( sourceReferenceMesh, imageSize ) ) / ( 2^16 -1 ) ...
  %                               * ( maxDeformation - minDeformation ) + minDeformation;
  tmp = kvlRasterizeAtlasMesh( sourceReferenceMesh, imageSize );
  
  % Unvisited voxels are marked by zeroes in all three coordinates - except possibly for the origin
  % which has all three coordinates zero as its natural state
  validMask = ( sum( tmp ~= 0, 4 ) ~= 0 );
  validMask( 1, 1, 1 ) = 1;
  if showFigures
    figure
    showImage( validMask );
  end  

  %
  denseDeformation = double( tmp ) / ( 2^16 -1 ) * ( maxDeformation - minDeformation ) + minDeformation;
                                 
  % Due to tetrahedral inside/outside checking performed in the rasterizer, some voxels on the boundary of the 
  % image grid are never visited. There also seems to be non-visited voxels inside the image grid (bug in rasterizer).
  % We can resolve both issues by giving to each non-visited voxel the deformation field value of the closest voxel 
  % that has been visited, augmented by knowledge that certain deformation components of voxels on the image grid
  % boundaries are zero by theory (sliding boundary conditions)
  
  [ distances, closestIndices ] = bwdist( validMask );
  for directionNumber = 1 : 3
    corrected = denseDeformation( :, :, :, directionNumber );
    corrected( find( ~validMask ) ) =  corrected( closestIndices( ~validMask ) );
    
    % Boundaries are known a priori
    if ( directionNumber == 1 )
      corrected( 1, :, : ) = 0;
      corrected( end, :, : ) = 0;
    elseif ( directionNumber == 2 )
      corrected( :, 1, : ) = 0;
      corrected( :, end, : ) = 0;
    else  
      corrected( :, :, 1 ) = 0;
      corrected( :, :, end ) = 0;
    end
    
    denseDeformation( :, :, :, directionNumber ) = corrected;
    
  end % End loop over directions
  
end % End test if there is deformation at all

% 
if showFigures
  for dimensionNumber = 1 : 3
    figure
    showImage( denseDeformation( :, :, :, dimensionNumber ) )
    drawnow
  end
end


%
if showFigures
  [ x y z ] = ndgrid( [ 0 : imageSize(1)-1 ], [ 0 : imageSize(2)-1 ], [ 0 : imageSize(3)-1 ] );
  densePositions = zeros( [ imageSize 3 ] );
  densePositions( :, :, :, 1 ) = x;
  densePositions( :, :, :, 2 ) = y;
  densePositions( :, :, :, 3 ) = z;

  figure
  q = quiver3( densePositions( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 1 ), ...
               densePositions( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 2 ), ...
               densePositions( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 3 ), ...
               denseDeformation( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 1 ), ...
               denseDeformation( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 2 ), ...
               denseDeformation( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 3 ) );
   title( 'Dense deformation' )  
   
end          
         
% OK this seems to work. Now get (interpolate) the deformation at the target mesh nodes (in reference position)
desiredTargetNodePosition = zeros( targetNumberOfNodes, 3 );
for directionNumber = 1 : 3
  desiredTargetNodePosition( :, directionNumber ) = ...
            targetReferencePosition( :, directionNumber ) + ...
            interpn( denseDeformation( :, :, :, directionNumber ), ...
                     targetReferencePosition( :, 1 )+1, ...
                     targetReferencePosition( :, 2 )+1, ...
                     targetReferencePosition( :, 3 )+1 );
end

%
if showFigures
  figure
  desiredDeformation = desiredTargetNodePosition - targetReferencePosition;
  quiver3( targetReferencePosition( 1 : 10 : end, 1 ), ...
          targetReferencePosition( 1 : 10 : end, 2 ), ...
          targetReferencePosition( 1 : 10 : end, 3 ), ...
          desiredDeformation( 1 : 10 : end, 1 ), ...
          desiredDeformation( 1 : 10 : end, 2 ), ...
          desiredDeformation( 1 : 10 : end, 3 ) );
  title( 'desired deformation' )        
end

          
% Now deform the target mesh to try and bring the position of its mesh nodes close(r) to their target positions         
image = kvlCreateImage( zeros( 10, 10, 10, 'single' ) ); % Dummy but we need it with the current interface
transform = kvlCreateTransform( eye( 4 ) );
calculator = kvlGetCostAndGradientCalculator( 'PointSet', image, 'Sliding', transform, [], [], [], [], ...
                                              desiredTargetNodePosition );

% Get an optimizer, and stick the cost function into it
optimizerType = 'L-BFGS'; % 'FixedStepGradientDescent','GradientDescent','ConjugateGradient', or 'L-BFGS'
maximalDeformationStopCriterion = 0.005; % Measured in voxels
lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion; % Doesn't seem to matter very much
optimizer = kvlGetOptimizer( optimizerType, targetReferenceMesh, calculator, ...
                                'Verbose', 1, ...
                                'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ... 
                                'LineSearchMaximalDeformationIntervalStopCriterion', ...
                                lineSearchMaximalDeformationIntervalStopCriterion, ...
                                'BFGS-MaximumMemoryLength', 12 ); % Affine registration only has 12 DOF
                                
numberOfIterations = 0;
startTime = tic;
if showFigures
  historyOfAverageDistance = [];
  historyOfMaximumDistance = [];
  iterationFigure = figure;
end
while true

  %
  if showFigures
    targetNodePositions = kvlGetMeshNodePositions( targetReferenceMesh );
    distances = sqrt( sum( ( targetNodePositions - desiredTargetNodePosition ).^2, 2 ) );
    averageDistance = sum( distances ) / length( distances )
    maximumDistance = max( distances );
    historyOfAverageDistance = [ historyOfAverageDistance; averageDistance ];
    historyOfMaximumDistance = [ historyOfMaximumDistance; maximumDistance ];
    
    figure( iterationFigure )
    subplot( 2, 1, 1 )
    plot( historyOfAverageDistance )
    grid
    title( [ 'average distance: ' num2str( historyOfAverageDistance( end ) ) ] )
    subplot( 2, 1, 2 )
    plot( historyOfMaximumDistance )
    grid
    title( [ 'max distance: ' num2str( historyOfMaximumDistance( end ) ) ] )
    drawnow
  end  

  %
  [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer )
  %return
  if ( maximalDeformation == 0 )
    break;
  end
  numberOfIterations = numberOfIterations + 1;

end
toc( startTime )

%
targetNodePositions = kvlGetMeshNodePositions( targetReferenceMesh );
targetDeformation = targetNodePositions - targetReferencePosition;
distances = sqrt( sum( ( targetNodePositions - desiredTargetNodePosition ).^2, 2 ) );
averageDistance = sum( distances ) / length( distances );
maximumDistance = max( distances );

if showFigures
  figure
  quiver3( targetReferencePosition( 1 : 10 : end, 1 ), ...
           targetReferencePosition( 1 : 10 : end, 2 ), ...
           targetReferencePosition( 1 : 10 : end, 3 ), ...
           targetDeformation( 1 : 10 : end, 1 ), ...
           targetDeformation( 1 : 10 : end, 2 ), ...
           targetDeformation( 1 : 10 : end, 3 ) );
  title( 'Achieved deformation' )
  
  figure
  nonAchievedDeformationComponent = targetDeformation - desiredDeformation;
  quiver3( targetReferencePosition( 1 : 10 : end, 1 ), ...
           targetReferencePosition( 1 : 10 : end, 2 ), ...
           targetReferencePosition( 1 : 10 : end, 3 ), ...
           nonAchievedDeformationComponent( 1 : 10 : end, 1 ), ...
           nonAchievedDeformationComponent( 1 : 10 : end, 2 ), ...
           nonAchievedDeformationComponent( 1 : 10 : end, 3 ), ...
           0 ); % Don't scale arrows
  title( 'Difference between achieved and desired deformation (true scale)' )
  
  
end



