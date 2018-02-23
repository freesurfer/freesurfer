meshCollectionFileName = '/home/willy/work/cm_m/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2/atlasForAffineRegistration.txt.gz';
imageFileName = '/home/willy/work/cm_m/innolitics_testing/buckner40/004/orig.mgz';
scaling = 0.9 * ones( 1, 3 );
scalingMatrix = diag( [ scaling 1 ] );
initialWorldToWorldTransformMatrix = scalingMatrix;
K = 1e-7;
K = K / (0.9 * 0.9 * 0.9);
meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
                                        kvlCreateTransform( initialWorldToWorldTransformMatrix), ...
                                        K  );
mesh = kvlGetMesh( meshCollection, -1 );
repeatCount = 3;
for thread_count = 1: 6
    kvlSetMaximumNumberOfThreads( thread_count );

    [ image, imageToWorldTransform ] = kvlReadImage( imageFileName );
  calculator = kvlGetCostAndGradientCalculator( 'MutualInformation', ...
                                                image, 'Affine' );

  maximalDeformationStopCriterion = 0.005; % Measured in voxels
  lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion; % Doesn't seem to matter very much
  optimizerType = 'L-BFGS';
  optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
                                  'Verbose', 1, ...
                                  'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ...
                                  'LineSearchMaximalDeformationIntervalStopCriterion', ...
                                  lineSearchMaximalDeformationIntervalStopCriterion, ...
                                  'BFGS-MaximumMemoryLength', 12 ); % Affine registration only has 12 DOF

    startingTime = tic;
    for i = 1: repeatCount
        [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer )
    end
    elapsedTime = toc(startingTime) / repeatCount;
    productivity = 1.0 / elapsedTime;
    thread_productivity = productivity / thread_count;
    fprintf('threads=%d time = %f thread_throughput=%f\n', thread_count, elapsedTime,thread_productivity);
end