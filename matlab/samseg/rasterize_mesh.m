meshCollectionFileName = '/home/willy/work/cm_p/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2/atlasForAffineRegistration.txt.gz';
scaling = 0.9 * ones( 1, 3 );
scalingMatrix = diag( [ scaling 1 ] );
initialWorldToWorldTransformMatrix = scalingMatrix;
K = 1e-7;
K = K / (0.9 * 0.9 * 0.9);
meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
                                        kvlCreateTransform( initialWorldToWorldTransformMatrix), ...
                                        K  );
mesh = kvlGetMesh( meshCollection, -1 );
imageSize = [200, 200, 200];

repeatCount = 10;
for thread_count = 1: 6
    kvlSetMaximumNumberOfThreads( thread_count );
    startingTime = tic;
    for i = 1: repeatCount
        priors = kvlRasterizeAtlasMesh( mesh, imageSize );
    end
    elapsedTime = toc(startingTime) / repeatCount;
    productivity = 1.0 / elapsedTime;
    thread_productivity = productivity / thread_count;
    fprintf('threads=%d time = %f thread_throughput=%f\n', thread_count, elapsedTime,thread_productivity);
end