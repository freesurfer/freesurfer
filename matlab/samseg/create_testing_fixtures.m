﻿
% Read in image, and figure out where the mesh is located
% test cases: read in the same image and assert transform and image buffer
% is the same
[ image, imageToWorldTransform ] = kvlReadImage( '/home/ys/freesurfer/GEMS2/Testing/test.nii' );
imageToWorldTransformMatrix = double( kvlGetTransformMatrix( imageToWorldTransform ) );
createdTransform = kvlCreateTransform( imageToWorldTransformMatrix );
imageBuffer = kvlGetImageBuffer( image );
image = kvlCreateImage( imageBuffer );
save('test_kvlImage.mat', 'imageToWorldTransformMatrix', 'imageBuffer');

% test case: write the same image out to disk and asser the file is created
% and when read in again we get the same thing.
kvlWriteImage( image, '/home/ys/freesurfer/GEMS2/Testing/test_written.nii', imageToWorldTransform );

% given a mesh collection file assert we get the same node positions and
% alphas and rasterization result.
meshCollection = kvlReadMeshCollection( '/usr/local/freesurfer/data/GEMS/CurrentMeshCollection30.gz' );
mesh = kvlGetMesh( meshCollection, -1 );
nodePositions = kvlGetMeshNodePositions( mesh );
alphas = kvlGetAlphasInMeshNodes( mesh );
priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
save('test_kvlMesh.mat', 'nodePositions', 'alphas', 'priors');

% Test mesh mutator kvlScaleMesh
meshCollection = kvlReadMeshCollection( '/usr/local/freesurfer/data/GEMS/CurrentMeshCollection30.gz' );
mesh = kvlGetMesh( meshCollection, -1 );
kvlScaleMesh(mesh, 1/10);
nodePositions = kvlGetMeshNodePositions( mesh );
save('test_kvlScaleMesh.mat', 'nodePositions');

% Test mesh mutator kvlSetAlphasInMeshNodes
meshCollection = kvlReadMeshCollection( '/usr/local/freesurfer/data/GEMS/CurrentMeshCollection30.gz' );
mesh = kvlGetMesh( meshCollection, -1 );
alphas = kvlGetAlphasInMeshNodes( mesh );
kvlSetAlphasInMeshNodes(mesh, alphas+1);
alphas = kvlGetAlphasInMeshNodes( mesh );
save('test_kvlSetAlphasInMeshNodes.mat', 'alphas');

% Test mesh mutator kvlSetMeshNodePositions
meshCollection = kvlReadMeshCollection( '/usr/local/freesurfer/data/GEMS/CurrentMeshCollection30.gz' );
mesh = kvlGetMesh( meshCollection, -1 );
nodePositions = kvlGetMeshNodePositions( mesh );
kvlSetMeshNodePositions(mesh, nodePositions+1);
nodePositions = kvlGetMeshNodePositions( mesh );
save('test_kvlSetMeshNodePositions.mat', 'nodePositions');


% cost and gradient calculator tests. Assert we get the same cost and
% gradient values.
calculator = kvlGetCostAndGradientCalculator( 'MutualInformation', image, 'Affine' );
[ cost gradient ] = kvlEvaluateMeshPosition( calculator, mesh );
save('test_kvlCostAndGradientMutualInformation.mat', 'cost', 'gradient')

% optimizer test. Run a few steps and write intermediate results to disk
history = []
optimizerType = 'L-BFGS';
maximalDeformationStopCriterion = 0.005; % Measured in voxels
lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion; % Doesn't seem to matter very much
optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
                                'Verbose', 1, ...
                                'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ...
                                'LineSearchMaximalDeformationIntervalStopCriterion', ...
                                lineSearchMaximalDeformationIntervalStopCriterion, ...
                                'BFGS-MaximumMemoryLength', 12 ); % Affine registration only has 12 DOF
[ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer );
save('test_kvlGetOptimizer.mat', 'minLogLikelihoodTimesPrior', 'maximalDeformation')
