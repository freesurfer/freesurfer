﻿
% Read in image, and figure out where the mesh is located
% test cases: read in the same image and assert transform and image buffer
% is the same
[ image, imageToWorldTransform ] = kvlReadImage( '/home/ys/freesurfer/GEMS2/Testing/test.nii' );
imageToWorldTransformMatrix = double( kvlGetTransformMatrix( imageToWorldTransform ) );
createdTransform = kvlCreateTransform( imageToWorldTransformMatrix );
imageBuffer = kvlGetImageBuffer( image );
image = kvlCreateImage( imageBuffer );

% given a mesh collection file assert we get the same node positions and
% alphas and rasterization result.
meshCollection = kvlReadMeshCollection( '/home/ys/freesurfer/GEMS2/Testing/test.txt.gz' );
mesh = kvlGetMesh( meshCollection, -1 );
nodePositions = kvlGetMeshNodePositions( mesh );
alphas = kvlGetAlphasInMeshNodes( mesh );
priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );


kvlScaleMesh(mesh, 1/10);
nodePositions_after_scaling = kvlGetMeshNodePositions( mesh );

kvlSetAlphasInMeshNodes(mesh, alphas+1);
alphas_after_modification = kvlGetAlphasInMeshNodes( mesh );


% Test mesh mutator kvlSetMeshNodePositions
kvlSetMeshNodePositions(mesh, nodePositions+1);
nodePositions_after_modification = kvlGetMeshNodePositions( mesh );

% cost and gradient calculator tests. Assert we get the same cost and
% gradient values.
calculator = kvlGetCostAndGradientCalculator( 'MutualInformation', image, 'Affine' );
[ mutualInformation_cost mutualInformation_gradient ] = kvlEvaluateMeshPosition( calculator, mesh );

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

save('/media/sf_matlab_data/integration.mat', ...
    'imageToWorldTransformMatrix', ...
    'imageBuffer', ...
    'nodePositions', ...
    'alphas', ...
    'priors', ...
    'nodePositions_after_scaling', ...
    'nodePositions_after_modification', ...
    'mutualInformation_cost', ...
    'mutualInformation_gradient', ...
    'minLogLikelihoodTimesPrior', ...
    'maximalDeformation')


