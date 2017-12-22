import numpy as np
import scipy.ndimage
import scipy.io
import GEMS2Python
import os
def require_np_array(np_array):
    return np.require(np_array, requirements=['F_CONTIGUOUS', 'ALIGNED'])

ASSERTIOINS_ON = False
def assert_close(golden, trial, **kwargs):
    if ASSERTIOINS_ON:
        np.testing.assert_allclose(golden, trial, **kwargs)

def samseg_registerAtlas(imageFileName,
                         meshCollectionFileName,
                         templateFileName,
                         savePath,
                         showFigures=False,
                         worldToWorldTransformMatrix=None):
    # function [ worldToWorldTransformMatrix, transformedTemplateFileName ] = samseg_registerAtlas( imageFileName, meshCollectionFileName, templateFileName, savePath, showFigures, worldToWorldTransformMatrix )
    # %
    # fixture = struct;
    # if ( nargin < 6 )
    #   worldToWorldTransformMatrix = [];
    # end
    image_dir = os.path.dirname(imageFileName)
    fixture = scipy.io.loadmat(os.path.join(image_dir, 'register_atlas_fixture.mat'),
                               struct_as_record=False, squeeze_me=True)['fixture']
    # % Print out the input
    # fprintf('entering registerAtlas\n');
    print('entering registerAtlas')
    # imageFileName
    print(imageFileName)
    # meshCollectionFileName
    print(meshCollectionFileName)
    # templateFileName
    print(templateFileName)
    # savePath
    print(savePath)
    # showFigures
    print(showFigures)
    # worldToWorldTransformMatrix
    print(worldToWorldTransformMatrix)
    #
    #
    # % Read in image and template, as well as their coordinates in world (mm) space
    # [ image, imageToWorldTransform ] = kvlReadImage( imageFileName );
    image = GEMS2Python.KvlImage(imageFileName)
    # imageToWorldTransformMatrix = double( kvlGetTransformMatrix( imageToWorldTransform ) );
    imageToWorldTransformMatrix = image.transform_matrix.as_numpy_array
    assert_close(fixture.sourceImageToWorldTransformMatrix, imageToWorldTransformMatrix)
    # [ template, templateImageToWorldTransform ] = kvlReadImage( templateFileName );
    template = GEMS2Python.KvlImage(templateFileName)

    # templateImageToWorldTransformMatrix = double( kvlGetTransformMatrix( templateImageToWorldTransform ) );
    templateImageToWorldTransformMatrix = template.transform_matrix.as_numpy_array
    assert_close(fixture.templateImageToWorldTransformMatrix, templateImageToWorldTransformMatrix)
    basepath, templateFileNameExtension = os.path.splitext(templateFileName)
    _, templateFileNameBase = os.path.split(basepath)
    #
    # [ ~, templateFileNameBase, templateFileNameExtension ] = fileparts( templateFileName );
    #
    #
    # %
    if worldToWorldTransformMatrix is None:
        # if ( isempty( worldToWorldTransformMatrix ) )
        #   %
        #   % The solution is not externally (secretly) given, so we need to compute it
        #   %
        #
        #   % Some hard-coded parameter settings first
        #   targetDownsampledVoxelSpacing = 3.0; % In mm
        targetDownsampledVoxelSpacing = 3.0
        #   K = 1e-7; % Mesh stiffness -- compared to normal models, the entropy cost function is normalized
        K = 1e-7
        #             % (i.e., measures an average *per voxel*), so that this needs to be scaled down by the
        #             % number of voxels that are covered
        #   maximalDeformationStopCriterion = 0.005; % Measured in voxels
        maximalDeformationStopCriterion = 0.005
        #   lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion; % Doesn't seem to matter very much
        lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion
        #
        #
        #   % Initialization
        #   initialWorldToWorldTransformMatrix = eye( 4 );
        #   if true
        #     % Provide an initial (non-identity) affine transform guestimate
        #
        #     % Rotation around X-axis (direction from left to right ear)
        #     theta = pi/180 * -10.0;
        theta = np.pi / 180 * -10.0
        #     rotationMatrix = eye( 4 );
        rotationMatrix = np.identity(4, dtype=np.double)
        #     rotationMatrix( 2 : 3, 2 : 3 ) = [ cos( theta ) -sin(theta); sin(theta) cos( theta ) ];
        #     initialWorldToWorldTransformMatrix = rotationMatrix * initialWorldToWorldTransformMatrix;
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        rotationMatrix[1, 1] = cos_theta
        rotationMatrix[1, 2] = -sin_theta
        rotationMatrix[2, 1] = sin_theta
        rotationMatrix[2, 2] = cos_theta
        #
        #     % Isotropic scaling
        #     scaling = 0.9 * ones( 1, 3 );
        scaling = 0.9
        #     scalingMatrix = diag( [ scaling 1 ] );
        scalingMatrix = np.diag([scaling, scaling, scaling, 1.0])
        #
        #     initialWorldToWorldTransformMatrix = scalingMatrix * initialWorldToWorldTransformMatrix;
        initialWorldToWorldTransformMatrix = scalingMatrix @ rotationMatrix

        #     K = K / prod( scaling );
        K = K / np.prod(scaling)
        #   end
        assert_close(fixture.initialWorldToWorldTransformMatrix, initialWorldToWorldTransformMatrix)


        #   initialImageToImageTransformMatrix = imageToWorldTransformMatrix \ ...
        #                     ( initialWorldToWorldTransformMatrix * templateImageToWorldTransformMatrix );
        multiplied = (initialWorldToWorldTransformMatrix @ templateImageToWorldTransformMatrix)
        initialImageToImageTransformMatrix = np.linalg.solve(imageToWorldTransformMatrix, multiplied)
        assert_close(fixture.imageToWorldTransformMatrix, imageToWorldTransformMatrix)
        assert_close(fixture.multiplied, initialWorldToWorldTransformMatrix @ templateImageToWorldTransformMatrix)
        assert_close(fixture.templateImageToWorldTransformMatrix, templateImageToWorldTransformMatrix)
        assert_close(fixture.initialWorldToWorldTransformMatrix, initialWorldToWorldTransformMatrix)
        assert_close(fixture.initialImageToImageTransformMatrix, initialImageToImageTransformMatrix)

        #
        #   % Figure out how much to downsample (depends on voxel size)
        #   voxelSpacing = sum( imageToWorldTransformMatrix( 1 : 3, 1 : 3 ).^2 ).^( 1/2 );
        voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3,0:3]**2, axis=0)**(1/2)
        #   downSamplingFactors = max( round( targetDownsampledVoxelSpacing ./ voxelSpacing ), [ 1 1 1 ] )
        downSamplingFactors = np.round(targetDownsampledVoxelSpacing / voxelSpacing)
        downSamplingFactors[downSamplingFactors < 1] = 1
        assert_close(fixture.downSamplingFactors, downSamplingFactors)
        #
        #   if 1
        #     % Use initial transform to define the reference (rest) position of the mesh (i.e., the one
        #     % where the log-prior term is zero)
        #     meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
        #                                             kvlCreateTransform( initialImageToImageTransformMatrix ), ...
        #                                             K * prod( downSamplingFactors ) );
        #     mesh = kvlGetMesh( meshCollection, -1 );
        mesh_collection = GEMS2Python.KvlMeshCollection()
        mesh_collection.read(meshCollectionFileName)
        mesh_collection.k = K * np.prod(downSamplingFactors)
        mesh_collection.transform(GEMS2Python.KvlTransform(require_np_array(initialImageToImageTransformMatrix)))
        mesh = mesh_collection.reference_mesh
        #   else
        #     % "Proper" initialization: apply the initial transform but don't let it affect the deformation
        #     % prior
        #     meshCollection = kvlReadMeshCollection( meshCollectionFileName, ...
        #                                             kvlCreateTransform( imageToWorldTransformMatrix \ ...
        #                                                                 templateImageToWorldTransformMatrix ), ...
        #                                             K * prod( downSamplingFactors ) );
        #     mesh = kvlGetMesh( meshCollection, -1 );
        #     nodePositions = kvlGetMeshNodePositions( mesh );
        #     tmp = [ nodePositions ones( size( nodePositions, 1 ), 1 ) ];
        #     tmp = ( imageToWorldTransformMatrix \ ...
        #             ( initialWorldToWorldTransformMatrix * imageToWorldTransformMatrix ) ) * tmp';
        #     nodePositions = tmp( 1:3, : )';
        #     kvlSetMeshNodePositions( mesh, nodePositions );
        #
        #   end
        #
        #   % Get image data
        #   imageBuffer = kvlGetImageBuffer( image );
        imageBuffer = image.getImageBuffer()
        #   assert_close(fixture.initialImageBuffer, imageBuffer;)
        assert_close(fixture.initialImageBuffer, imageBuffer)
        #   if showFigures
        #     figure
        #     showImage( imageBuffer );
        #   end
        #
        #   % Downsample
        #   imageBuffer = imageBuffer( 1 : downSamplingFactors( 1 ) : end, ...
        #                              1 : downSamplingFactors( 2 ) : end, ...
        #                              1 : downSamplingFactors( 3 ) : end );
        imageBuffer = imageBuffer[
                      ::int(downSamplingFactors[0]),
                      ::int(downSamplingFactors[1]),
                      ::int(downSamplingFactors[2])]
        #   assert_close(fixture.downsampledImageBuffer, imageBuffer;)
        assert_close(fixture.downsampledImageBuffer, imageBuffer)
        #   image = kvlCreateImage( imageBuffer );
        image = GEMS2Python.KvlImage(require_np_array(imageBuffer))
        #   kvlScaleMesh( mesh, 1 ./ downSamplingFactors );
        assert_close(fixture.preScaleMesh, mesh.points)
        mesh.scale(1 / downSamplingFactors)
        assert_close(fixture.postScaleMesh, mesh.points)
        #
        #   alphas = kvlGetAlphasInMeshNodes( mesh );
        alphas = mesh.alphas
        assert_close(fixture.downsampledAlphas, alphas)
        #   gmClassNumber = 3;  % Needed for displaying purposes
        gmClassNumber = 3  # Needed for displaying purposes
        #   if 0
        #     % Get rid of the background class
        #     alphas = alphas( :, 2 : end );
        #     kvlSetAlphasInMeshNodes( mesh, alphas )
        #     gmClassNumber = gmClassNumber-1;
        #   end
        #   numberOfClasses = size( alphas, 2 );
        numberOfClasses = alphas.shape[1]
        #   colors = 255 * [ hsv( numberOfClasses ) ones( numberOfClasses, 1 ) ];

        #
        #   if showFigures
        #     figure
        #     priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
        #     for classNumber = 1 : numberOfClasses
        #       subplot( 2, 3, classNumber )
        #       showImage( priors( :, :, :, classNumber ) )
        #     end
        #   end
        #
        #   %
        #   % Get a registration cost and use it to evaluate some promising starting point proposals
        #   calculator = kvlGetCostAndGradientCalculator( 'MutualInformation', ...
        #                                                 image, 'Affine' );

        calculator = GEMS2Python.KvlCostAndGradientCalculator('MutualInformation', [image], 'Affine')
        #   [ cost gradient ] = kvlEvaluateMeshPosition( calculator, mesh );
        cost, gradient = calculator.evaluate_mesh_position(mesh)
        #   assert_close(fixture.initialCost, cost;)
        #   assert_close(fixture.initialGradient, gradient;)
        assert_close(fixture.initialCost, cost)
        assert_close(fixture.initialGradient, gradient)
        #   if true
        #     %
        #     [ xtmp, ytmp, ztmp ] = ndgrid( 1 : size( imageBuffer, 1 ), ...
        #                                    1 : size( imageBuffer, 2 ), ...
        #                                    1 : size( imageBuffer, 3 ) );
        #     centerOfGravityImage = [ xtmp(:) ytmp(:) ztmp(:) ]' * imageBuffer(:) / sum( imageBuffer(:) );

        # the + 1 is to account for MATLAB being 1 indexed
        centerOfGravityImage = np.array(scipy.ndimage.measurements.center_of_mass(imageBuffer))

        assert_close(fixture.centerOfGravityImage, np.array(centerOfGravityImage) + 1, rtol=1e-5)

        #     priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
        priors = mesh.rasterize(imageBuffer.shape)
        #     %tmp = sum( priors, 4 );
        #     tmp = sum( priors( :, :, :, 2 : end ), 4 );
        tmp = np.sum(priors[:, :, :, 1:], axis=3)
        #     centerOfGravityAtlas = [ xtmp(:) ytmp(:) ztmp(:) ]' * tmp(:) / sum( tmp(:) );
        centerOfGravityAtlas = np.array(scipy.ndimage.measurements.center_of_mass(tmp))
        #     %
        #     initialTranslation = double( centerOfGravityImage - centerOfGravityAtlas );
        initialTranslation = centerOfGravityImage - centerOfGravityAtlas
        #     nodePositions = kvlGetMeshNodePositions( mesh );
        nodePositions = mesh.points
        #     trialNodePositions = nodePositions + repmat( initialTranslation', [ size( nodePositions, 1 ) 1 ] );
        trialNodePositions = nodePositions + initialTranslation
        assert_close(fixture.trialNodePositions, trialNodePositions, rtol=1e-5)
        #     kvlSetMeshNodePositions( mesh, trialNodePositions );
        mesh.points = trialNodePositions
        #     [ trialCost trialGradient ] = kvlEvaluateMeshPosition( calculator, mesh );
        trialCost, trialGradient = calculator.evaluate_mesh_position(mesh)
        #     if ( trialCost >= cost )
        if trialCost >= cost:
            #       % Center of gravity was not a success; revert to what we had before
            #       kvlSetMeshNodePositions( mesh, nodePositions );
            mesh.points = nodePositions
        #     else
        else:
            #       % This is better starting position; remember that we applied it
            #       initialImageToImageTransformMatrix( 1 : 3, 4 ) = ...
            #               initialImageToImageTransformMatrix( 1 : 3, 4 ) + diag( downSamplingFactors ) * initialTranslation;
            initialImageToImageTransformMatrix[0:3, 3] = initialImageToImageTransformMatrix[0:3, 3] + (np.diag( downSamplingFactors ) @ initialTranslation)
        #   end
        #
        #   end
        #
        #   %
        #   originalNodePositions = kvlGetMeshNodePositions( mesh );
        # assert_close(fixture.initialImageToImageTransformMatrix, initialImageToImageTransformMatrix)
        originalNodePositions = mesh.points
        #   assert_close(fixture.originalNodePositions, originalNodePositions;)
        assert_close(fixture.originalNodePositions, originalNodePositions, rtol=1e-5)

        #   % Visualize starting situation
        #   priorVisualizationAlpha = 0.4;
        #   if showFigures
        #     figure
        #
        #     priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
        #     colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
        #     mask = ( sum( double( priors ), 4 ) / (2^16-1) ) > .5;
        #     subplot( 2, 2, 1 )
        #     showImage( imageBuffer );
        #     subplot( 2, 2, 2 )
        #     showImage( imageBuffer .* mask );
        #     subplot( 2, 2, 3 )
        #     imageToShow = ( imageBuffer - min( imageBuffer(:) ) ) / ( max( imageBuffer(:) ) - min( imageBuffer(:) ) );
        #     imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( imageToShow, [ 1 1 1 3 ] ) + ...
        #                   priorVisualizationAlpha * colorCodedPriors;
        #     showImage( imageToShow )
        #     subplot( 2, 2, 4 )
        #     tmp = double( priors( :, :, :, gmClassNumber ) ) / ( 2^16-1 );
        #     showImage( mosaicImages( tmp, imageBuffer .* mask, 2 ) );
        #
        #     drawnow
        #   end
        #
        #   % Get an optimizer, and stick the cost function into it
        #   optimizerType = 'L-BFGS';
        #   optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
        #                                   'Verbose', 1, ...
        #                                   'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ...
        #                                   'LineSearchMaximalDeformationIntervalStopCriterion', ...
        #                                   lineSearchMaximalDeformationIntervalStopCriterion, ...
        #                                   'BFGS-MaximumMemoryLength', 12 ); % Affine registration only has 12 DOF
        optimization_parameters = {
            'Verbose': 1.0,
            'MaximalDeformationStopCriterion': maximalDeformationStopCriterion,
            'LineSearchMaximalDeformationIntervalStopCriterion': lineSearchMaximalDeformationIntervalStopCriterion,
            'BFGS-MaximumMemoryLength': 12.0  # Affine registration only has 12 DOF
        }
        optimizer = GEMS2Python.KvlOptimizer(
            'L-BFGS',
            mesh,
            calculator,
            optimization_parameters
        )

        #   numberOfIterations = 0;
        numberOfIterations = 0
        #   tic
        minLogLikelihoodTimesPriors = []
        maximalDeformations = []
        costs = []
        gradients = []
        #   while truek
        while True:
            cost, gradient = calculator.evaluate_mesh_position(mesh)
            print("cost = {}".format(cost))
            costs.append(cost)
            gradients.append(gradient)
            #     %
            #     [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer )
            minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer()
            #     assert_close(fixture.minLogLikelihoodTimesPriors, [fixture.minLogLikelihoodTimesPriors, minLogLikelihoodTimesPrior];)
            minLogLikelihoodTimesPriors.append(minLogLikelihoodTimesPrior)
            #     assert_close(fixture.maximalDeformations, [fixture.maximalDeformations, maximalDeformation];)
            maximalDeformations.append(maximalDeformation)
            #     %return
            #     if ( maximalDeformation == 0 )
            #       break;
            #     end
            if maximalDeformation == 0:
                break
            numberOfIterations += 1
            # scipy.io.savemat(
            #     '/Users/ys/work/freesurfer/GEMS2/Testing/matlab_data/integration_cost.mat',
            #         {
            #             'points': mesh.points,
            #             'pythonCost': cost,
            #         }
            # )
            #
            #
            #     %
            #     % Visualize progress
            #     %
            #     if showFigures
            #       % Show figure
            #       priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
            #       colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
            #       mask = ( sum( double( priors ), 4 ) / (2^16-1) ) > .5;
            #       subplot( 2, 2, 1 )
            #       showImage( imageBuffer );
            #       subplot( 2, 2, 2 )
            #       showImage( imageBuffer .* mask );
            #       subplot( 2, 2, 3 )
            #       imageToShow = ( imageBuffer - min( imageBuffer(:) ) ) / ( max( imageBuffer(:) ) - min( imageBuffer(:) ) );
            #       imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( imageToShow, [ 1 1 1 3 ] ) + ...
            #                     priorVisualizationAlpha * colorCodedPriors;
            #       showImage( imageToShow )
            #       subplot( 2, 2, 4 )
            #       tmp = double( priors( :, :, :, gmClassNumber ) ) / ( 2^16-1 );
            #       showImage( mosaicImages( tmp, imageBuffer .* mask, 2 ) );
            #       drawnow
            #
            #       % Show affine matrix, retrieved from any four non-colinear points before and after registration
            #       nodePositions = kvlGetMeshNodePositions( mesh );
            #       pointNumbers = [ 1 111 202 303 ];
            #       originalY = [ originalNodePositions( pointNumbers, : )'; 1 1 1 1 ];
            #       Y = [ nodePositions( pointNumbers, : )'; 1 1 1 1 ];
            #       extraImageToImageTransformMatrix = Y * inv( originalY );
            #       scaling = svd( extraImageToImageTransformMatrix( 1 : 3, 1 : 3 ) );
            #       disp( [ 'scaling: ' num2str( scaling' ) ] )
            #     end
            #
            #   end % End loop over iterations
            #   numberOfIterations
            #   toc

        # assert_close(fixture.costs, costs)
        # assert_close(fixture.gradients, gradients)
        # assert_close(fixture.minLogLikelihoodTimesPriors, minLogLikelihoodTimesPriors)
        # assert_close(fixture.maximalDeformations, maximalDeformations)

        #
        #   % For debugging and/or quality control purposes, save a picture of the registration result to file
        #   priors = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ) );
        #   colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
        #   mask = ( sum( double( priors ), 4 ) / (2^16-1) ) > .5;
        #   overlayQcImage = ( imageBuffer - min( imageBuffer(:) ) ) / ( max( imageBuffer(:) ) - min( imageBuffer(:) ) );
        #   overlayQcImage = ( 1 - priorVisualizationAlpha ) * repmat( overlayQcImage, [ 1 1 1 3 ] ) + ...
        #                 priorVisualizationAlpha * colorCodedPriors;
        #
        #   tmp = double( priors( :, :, :, gmClassNumber ) ) / ( 2^16-1 );
        #   mosaicQcImage = mosaicImages( tmp, imageBuffer .* mask, 2 );
        #
        #   overlayCollage = getCollage( overlayQcImage, 10 );
        #   mosaicCollage = getCollage( mosaicQcImage, 10 );
        #
        #   borderSize = 20;
        #   DIM = [ size( overlayCollage, 1 ) size( overlayCollage, 2 ) ];
        #   qcFigure = zeros( [ DIM( 1 )  2 * DIM( 2 ) 3 ]  + [ 2*borderSize 3*borderSize 0 ] ) + .5;
        #   qcFigure( borderSize + [ 1 : DIM( 1 ) ], borderSize + [ 1 : DIM( 2 ) ], : ) = overlayCollage;
        #   qcFigure( borderSize + [ 1 : DIM( 1 ) ], ...
        #             2 * borderSize + DIM( 2 ) + [ 1 : DIM( 2 ) ], : ) = mosaicCollage;
        #   qcFigureFileName = fullfile( savePath, ...
        #                               [ templateFileNameBase '_coregistrationCqFigure.png' ] );
        #   imwrite( qcFigure, qcFigureFileName )
        #
        #
        #   % Retrieve the implicitly applied affine matrix from any four non-colinear points before and after registration,
        #   % taking into account the downsampling that we applied
        #   nodePositions = kvlGetMeshNodePositions( mesh );
        nodePositions = mesh.points
        # assert_close(fixture.finalNodePositions, nodePositions;)
        # assert_close(fixture.finalNodePositions, nodePositions, atol=2)
        #   pointNumbers = [ 1 111 202 303 ];
        pointNumbers = [0, 110, 201, 302 ]
        #   originalY = [ diag( downSamplingFactors ) * originalNodePositions( pointNumbers, : )'; 1 1 1 1 ];
        originalY = np.vstack((np.diag( downSamplingFactors ) @ originalNodePositions[pointNumbers].T, [1, 1, 1, 1]))

        #   Y = [ diag( downSamplingFactors ) * nodePositions( pointNumbers, : )'; 1 1 1 1 ];
        Y = np.vstack((np.diag( downSamplingFactors ) @ nodePositions[pointNumbers].T, [1, 1, 1, 1]))
        #   extraImageToImageTransformMatrix = Y * inv( originalY );
        extraImageToImageTransformMatrix = Y @ np.linalg.inv(originalY)
        #
        #   % Final result: the image-to-image (from template to image) as well as the world-to-world transform that
        #   % we computed (the latter would be the identity matrix if we didn't move the image at all)
        #   imageToImageTransformMatrix = extraImageToImageTransformMatrix * initialImageToImageTransformMatrix;
        imageToImageTransformMatrix = extraImageToImageTransformMatrix @ initialImageToImageTransformMatrix
        #   worldToWorldTransformMatrix = imageToWorldTransformMatrix * imageToImageTransformMatrix * ...
        #                                 inv( templateImageToWorldTransformMatrix );
        worldToWorldTransformMatrix = imageToWorldTransformMatrix @ imageToImageTransformMatrix @ np.linalg.inv( templateImageToWorldTransformMatrix )


        # assert_close(fixture.imageToImageTransformMatrix[:, :3], imageToImageTransformMatrix[:, :3], atol=2e-2)
        # assert_close(fixture.imageToImageTransformMatrix[:, :4], imageToImageTransformMatrix[:, :4], atol=1)
        # assert_close(fixture.worldToWorldTransformMatrix[:, :3], worldToWorldTransformMatrix[:, :3], atol=2e-2)
        # assert_close(fixture.worldToWorldTransformMatrix[:, :4], worldToWorldTransformMatrix[:, :4], atol=1)
        # #
    # else
    else:
        #   % The world-to-world transfrom is externally given, so let's just compute the corresponding image-to-image
        #   % transform (needed for subsequent computations) and be done
        #   imageToImageTransformMatrix = inv( imageToWorldTransformMatrix ) * worldToWorldTransformMatrix * ...
        #                                 templateImageToWorldTransformMatrix
        imageToImageTransformMatrix = np.linalg.inv(imageToWorldTransformMatrix) * worldToWorldTransformMatrix @ templateImageToWorldTransformMatrix

        # end % End test if the solution is externally given
    #
    #
    # % Save the image-to-image and the world-to-world affine registration matrices

    scipy.io.savemat(os.path.join(image_dir, templateFileNameBase+'_coregistrationMatrices.mat'),
                     {'imageToImageTransformMatrix': imageToImageTransformMatrix,
                      'worldToWorldTransformMatrix': worldToWorldTransformMatrix,
                      'costs': costs,
                      }
                     )
    print(fixture.imageToImageTransformMatrix - imageToImageTransformMatrix)

    # transformationMatricesFileName = fullfile( savePath, ...
    #                                            [ templateFileNameBase '_coregistrationMatrices.mat' ] );
    # eval( [ 'save ' transformationMatricesFileName ' imageToImageTransformMatrix worldToWorldTransformMatrix;' ] )
    #
    #
    # % Compute the talairach.xfm
    # % Load fsaverage orig.mgz -- this is the ultimate target/destination
    # fshome = getenv('FREESURFER_HOME');
    # fnamedst = sprintf('%s/subjects/fsaverage/mri/orig.mgz',fshome);
    # fsaorig = MRIread(fnamedst,1);
    # % Compute the vox2vox from the template to fsaverage assuming they
    # %   share world RAS space
    # RAS2LPS = diag([-1 -1 1 1]);
    # M = inv(RAS2LPS*fsaorig.vox2ras)*(templateImageToWorldTransformMatrix);
    # % Compute the input to fsaverage vox2vox by combining the
    # % input-template vox2vox and the template-fsaverage vox2vox
    # X = M*inv(imageToImageTransformMatrix);
    # % Now write out the LTA. This can be used as the talairach.lta in recon-all
    # invol = MRIread(imageFileName,1); % have to reread to get header info
    # lta.type = 0;
    # lta.xform = X;
    # lta.srcfile = imageFileName;
    # lta.srcmri = invol;
    # lta.srcmri.vol = [];
    # lta.dstfile = fnamedst;
    # lta.dstmri = fsaorig;
    # lta.dstmri.vol = [];
    # lta.subject = 'fsaverage';
    # ltaFileName = sprintf('%s/samseg.talairach.lta',savePath);
    # lta_write(ltaFileName,lta);
    # fprintf('Done computng and writing out LTA %s\n',ltaFileName);
    #
    #
    # % For historical reasons, we applied the estimated transformation to the template; let's do that now
    # desiredTemplateImageToWorldTransformMatrix = imageToWorldTransformMatrix * imageToImageTransformMatrix
    desiredTemplateImageToWorldTransformMatrix = imageToWorldTransformMatrix @ imageToImageTransformMatrix

    # transformedTemplateFileName = fullfile( savePath, ...
    #                                         [ templateFileNameBase '_coregistered' templateFileNameExtension ] );
    transformedTemplateFileName = os.path.join(savePath, templateFileNameBase, '_coregistered' + templateFileNameExtension)
    # kvlWriteImage( template, transformedTemplateFileName, ...
    #                kvlCreateTransform( desiredTemplateImageToWorldTransformMatrix ) );
    template.write(transformedTemplateFileName, GEMS2Python.KvlTransform(desiredTemplateImageToWorldTransformMatrix))
    # assert_close(fixture.desiredTemplateImageToWorldTransformMatrix, desiredTemplateImageToWorldTransformMatrix;)
    #
    # save('/media/sf_matlab_data/register_atlas_fixture.mat', 'fixture')
    #
    return worldToWorldTransformMatrix, transformedTemplateFileName

if __name__ == '__main__':
    import os
    affineRegistrationMeshCollectionFileName = '/Users/ys/Downloads/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2/atlasForAffineRegistration.txt.gz'
    templateFileName = '/Users/ys/Downloads/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2/template.nii'
    rootPath = '/Users/ys/Downloads/innolitics_testing/buckner40/'
    for f in os.listdir(rootPath):
        fullpath = os.path.join(rootPath, f)
        if os.path.isdir(fullpath):
            imageFileName = os.path.join(fullpath, 'orig.mgz')
            savePath = fullpath
            samseg_registerAtlas(imageFileName=imageFileName,
                                 meshCollectionFileName=affineRegistrationMeshCollectionFileName,
                                 templateFileName=templateFileName,
                                 savePath=savePath,
                                 showFigures=False,
                                 worldToWorldTransformMatrix=None)


    # fullpath = os.path.join(rootPath, '008')
    # if os.path.isdir(fullpath):
    #     imageFileName = os.path.join(fullpath, 'orig.mgz')
    #     savePath = fullpath
    #     samseg_registerAtlas(imageFileName=imageFileName,
    #                          meshCollectionFileName=affineRegistrationMeshCollectionFileName,
    #                          templateFileName=templateFileName,
    #                          savePath=savePath,
    #                          showFigures=False,
    #                          worldToWorldTransformMatrix=None)
