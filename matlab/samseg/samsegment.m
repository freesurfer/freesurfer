function samsegment( imageFileNames, transformedTemplateFileName, meshCollectionFileName, compressionLookupTableFileName, ...
                     modelSpecifications, optimizationOptions, savePath, showFigures )
%
%

% TODO: 
%   - the separate cases for uni- vs. multi-contrast Gaussian model parameters is mind-boggling
%     and utterly unnecessary: not only do the dimensions swap, also their names are inexplicably 
%     different ("variances" vs. "sigmas"). For the bias field coefficients EM update and mesh node
%     position optimizer I already introduced "todoMeans" and "todoVariances" that handle both cases 
%     transparently in the same way -- this should be done in all the rest of the code also (removing
%     all the tests if numberOfContrasts > 1 from the code). Also: should we really keep variances 
%     around, are precisions instead (I think the latter)?
%   - The separate function for writing out results as non-cropped images at the very end should 
%     be removed and its functionality (returning information of how large the original volume is
%     and where the cropped region we're segmenting is located within that volume) should be folded 
%     into the cropped reader. Right now it's prone to bugs because a manually-set extra margin 
%     needs to be consistent between two separate Mex functions (!)
%   - make a distinction between variables "numberOfGaussians" (e.g., 17 individual Gaussians) 
%     vs. "numberOfClasses" (e.g., 7 classes, each represented by a Gaussian mixture model) 
%     throughout the code 
%   - allow for priors on intensity model (e.g., ridge regression on diagonal of lhs bias field correction
%     matrix). Right now there is a tiny conjugate prior on the Gaussian variances without taking it
%     into account in cost function evaluation.
%      



% Print input options
disp( '==========================' );
disp( 'imageFileNames: ' )
for i = 1 : length( imageFileNames )
  disp( imageFileNames{ i } )
end
fprintf( '-----\n' )

disp( 'transformedTemplateFileName: ' )
disp( transformedTemplateFileName )
fprintf( '-----\n' )

disp( 'meshCollectionFileName: ' )
disp( meshCollectionFileName )
fprintf( '-----\n' )

disp( 'compressionLookupTableFileName: ' )
disp( compressionLookupTableFileName )
fprintf( '-----\n' )

disp( 'modelSpecifications: ' )
fieldNames = fieldnames( modelSpecifications );
for i = 1 : length( fieldNames )
  fieldName = fieldNames{ i };
  disp( fieldName )
  fields = getfield( modelSpecifications, fieldName );
  for j = 1 : length( fields )
    disp( fields(j) )
  end
  disp( ' ' )
end
fprintf( '-----\n' )  

disp( 'optimizationOptions:' )
fieldNames = fieldnames( optimizationOptions );
for i = 1 : length( fieldNames )
  fieldName = fieldNames{ i };
  disp( fieldName )
  fields = getfield( optimizationOptions, fieldName );
  for j = 1 : length( fields )
    disp( fields(j) )
  end
  disp( ' ' )
end
fprintf( '-----\n' )  

disp( 'savePath: ' )
disp( savePath )
fprintf( '-----\n' )

disp( 'showFigures: ' )
disp( showFigures )
fprintf( '-----\n' )
  



%
numberOfContrasts = length( imageFileNames );

% Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
% translation, rotation, scaling, and skewing) as well - this transformation will later be used
% to initially transform the location of the atlas mesh's nodes into the coordinate system of
% the image.
%
imageBuffers = [];
for contrastNumber = 1 : numberOfContrasts
  % Get the pointers to image and the corresponding transform
  [ images( contrastNumber ), transform ] = kvlReadCroppedImage( imageFileNames{ contrastNumber }, transformedTemplateFileName ); 
  imageBuffers( :, :, :, contrastNumber ) = kvlGetImageBuffer( images( contrastNumber ) ); % Get the actual imageBuffer
end

if ( showFigures )
  for contrastNumber = 1 : numberOfContrasts
    figure
    showImage( imageBuffers( :, :, :, contrastNumber ) ); % Automatically displays middle slices in each direction
  end
end

% Also read in the voxel spacing -- this is needed since we'll be specifying bias field smoothing kernels, downsampling 
% steps etc in mm.
[ ~, at ] = kvlReadImage( imageFileNames{1} );
atm = double( kvlGetTransformMatrix( at ) );
voxelSpacing = sum( atm( 1 : 3, 1 : 3 ).^2 ).^( 1/2 );



% Read the atlas mesh from file, immediately applying the previously determined transform to the location
% of its nodes. Rather than one mesh, the atlas consists of a so-called "collection" of meshes: they're
% all the same except for the location of their mesh nodes. The first mesh, so-called "reference mesh"
% has index -1: it represents the average shape (location of the mesh nodes) of all the other meshes. If
% the atlas was built from N training subjects, there will be N other meshes as well, each one warping
% the atlas optimally to each of the training subjects
%
% You also need to provide a value for K, which determines the flexibility of the atlas mesh, i.e., how
% much it will typically deform. Higher values correspond to stiffer meshes.
%
meshCollection = kvlReadMeshCollection( meshCollectionFileName, transform, modelSpecifications.K );

% Retrieve the reference mesh, i.e., the mesh representing the average shape.
mesh = kvlGetMesh( meshCollection, -1 );



% Skull strip the images
labelNumber = 0;
sz = [size(imageBuffers,1) size(imageBuffers,2) size(imageBuffers,3)];
backgroundPrior = kvlRasterizeAtlasMesh( mesh, sz, labelNumber );

if( showFigures )
  figure
  subplot( 2, 2, 1 )
  showImage( backgroundPrior )
  subplot( 2, 2, 2 )
  showImage( mosaicImages( 2^16 - 1 - double( backgroundPrior ), double( imageBuffers(:,:,:,1) ), 10 ) )
end
smoothedBackgroundPrior = kvlSmoothImageBuffer( backgroundPrior, modelSpecifications.brainMaskingSmoothingSigma ./ voxelSpacing );
if( showFigures )
  subplot( 2, 2, 3 )
  showImage( smoothedBackgroundPrior )
end
brainMask = ( 1 - single( smoothedBackgroundPrior ) / 65535 ) > modelSpecifications.brainMaskingThreshold;

for contrastNumber = 1 : numberOfContrasts
  imageBuffer = imageBuffers( :, :, :, contrastNumber );
  imageBuffer( find( ~brainMask ) ) = 0;
  imageBuffers( :, :, :, contrastNumber ) = imageBuffer;
  % kvlSetImageBuffer( images( contrastNumber ), imageBuffers( :, :, :, contrastNumber ) );
end

if( showFigures )
  subplot( 2, 2, 4 )
  showImage( imageBuffers( :, :, :, 1 ) )
end



% Let's prepare for the bias field correction that is part of the imaging model. It assumes
% an additive effect, whereas the MR physics indicate it's a multiplicative one - so we log
% transform the data first. In order to do so, mask out zeros from the images.
mask = ( imageBuffers( :, :, :, 1 ) > 0 );
for contrastNumber = 2 : numberOfContrasts
  mask = mask .* ( imageBuffers( :, :, :, contrastNumber ) > 0 );
end
maskIndices = find( mask );
for contrastNumber = 1 : numberOfContrasts
  buffer = imageBuffers( :, :, :, contrastNumber );
  buffer( maskIndices ) = 1000 * log( buffer( maskIndices ) ); % The 1000 factor isn't really necessary; it's just
                                                              % easier to have Gaussian means of 100 350 than
                                                              % 0.1 and 0.35 for my human brain
  buffer = buffer .* mask;
  imageBuffers( :, :, :, contrastNumber ) = buffer;
end



% FreeSurfer (http://surfer.nmr.mgh.harvard.edu) has a standardized way of representation segmentations,
% both manual and automated, as images in which certain intensity levels correspond to well-defined
% anatomical structures - for instance an intensity value 17 always corresponds to the left hippocampus.
% The text file FreeSurferColorLUT.txt distributed with FreeSurfer contains all these definitions, as well 
% as a color for each structure in
% RGBA (Red-Green-Blue and Alpha (opaqueness)) format with which FreeSurfer will always represent segmented
% structures in its visualization tools
%
% Let's read the contents of this "compressionLookupTable.txt" file, and show the names of the structures
% being considered. The results are automatically sorted according to their "compressed label", i.e., the
% first result corresponds to the first entry in the vector of probabilities associated with each node in
% our atlas mesh.
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );

% Get a Matlab matrix containing a copy of the probability vectors in each mesh node (size numberOfNodes x
% numberOfLabels ). 
alphas = kvlGetAlphasInMeshNodes( mesh );

% Remove any structures that don't exist in the data (e.g., half a brain missing in ex vivo images)
mergeOptions = struct;
mergeOptions( 1 ).mergedName = 'Unknown';
mergeOptions( 1 ).searchStrings = modelSpecifications.missingStructureSearchStrings;
[ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );
kvlSetAlphasInMeshNodes( mesh, alphas );



% Although we have many labels to segment, and each of these labels has its own Gaussian mixture model
% whose parameters (mean, variance, mixture weight) we have to estimate from the data, it may makes sense to restrict
% the degrees of freedom in the model somewhat by specifying that some of these labels have the same parameters
% governing their Gaussian mixture model. For example, we'd expect no intensity differences between the left and right 
% part of each structure.
%
% Compute the "reduced" alphas - those referring to the "super"-structures. At the same
% time also build an inverse lookup table (mapping from original class number onto a
% reduced class number) that we will need to compute the final segmentation.
[ reducedAlphas, reducedNames, reducedFreeSurferLabels, reducedColors, reducingLookupTable ] = ...
                            kvlMergeAlphas( alphas, names, modelSpecifications.sharedGMMParameters, FreeSurferLabels, colors );


if ( showFigures )
  fprintf( 'Visualizing the atlas mesh; this takes quite some time and is only here for tutorial purposes...' )
  priors = kvlRasterizeAtlasMesh( mesh, [size(imageBuffers,1) size(imageBuffers,2) size(imageBuffers,3)]); % Without specifying a specific label, will rasterize all simultaneously, rasterize everything
  rgbBuffer = kvlColorCodeProbabilityImages( priors, colors );
  figure
  showImage( rgbBuffer )
  clear priors rgbBuffer
  fprintf( 'done\n' )
  drawnow;
end




% Initialize classNumber weights, which are at first split evenly
numberOfGaussiansPerClass = [ modelSpecifications.sharedGMMParameters.numberOfComponents ];
numberOfClasses = length( numberOfGaussiansPerClass );
EMweights=zeros( 1, sum( numberOfGaussiansPerClass ) );
shift = 0;
for classNumber = 1 : numberOfClasses
  EMweights( classNumber + shift : classNumber + shift + numberOfGaussiansPerClass( classNumber ) - 1 ) = ...
                                                                            1 / numberOfGaussiansPerClass( classNumber );
  shift = shift + numberOfGaussiansPerClass( classNumber ) - 1;
end



% Our bias model is a linear combination of a set of basis functions. We are using so-called
% "DCT-II" basis functions, i.e., the lowest few frequency components of the Discrete Cosine
% Transform.
%
biasFieldSmoothingKernelSize = 50.0;  % Distance in mm of sinc function center to first zero crossing
kroneckerProductBasisFunctions = cell( 0, 0 );
numberOfBasisFunctions = zeros( 1, 3 );
for dimensionNumber = 1 : 3
  N = size( imageBuffers, dimensionNumber ); % Number of data points
  delta =  biasFieldSmoothingKernelSize / voxelSpacing( dimensionNumber ); % Measured in voxels
  M = ceil( N / delta ) + 1; % Number of frequencies to use
  Nvirtual = ( M - 1 ) * delta; % Virtual (possibly non-integer) number of data points to make sure 
                                % we reach the targeted smoothing kernel size
  js = ( [ 0 : N-1 ]' + 1/2 ) * pi / Nvirtual;
  A = cos( js * [ 0 : M-1 ] ) * sqrt( 2 / Nvirtual );
  A( :, 1 ) = A( :, 1 ) / sqrt( 2 );
  
  if showFigures
    % Show smoothing kernel
    figure
    smootherMatrix = A * ( ( A' * A ) \ A' );
    subplot( 2, 2, 1 )
    imshow( smootherMatrix, [] )
    subplotCounter = 2;
    for rowNumber = round( [ N/4 N/2 3*N/4 ] )
      subplot( 2, 2, subplotCounter )
      plot( smootherMatrix( rowNumber, : ) )
      set( gca, 'xlim', [ 1 N ] )
      grid
      hold on
      ylim = get( gca, 'ylim' );
      line( ( rowNumber ) * [ 1 1 ], ylim, 'color', 'k' )
      line( ( rowNumber + delta ) * [ 1 1 ], ylim, 'color', 'r', 'linestyle', '--' )
      line( ( rowNumber - delta ) * [ 1 1 ], ylim, 'color', 'r', 'linestyle', '--' )
      subplotCounter = subplotCounter + 1;
    end
  end
    
  kroneckerProductBasisFunctions{ dimensionNumber } = A;
  numberOfBasisFunctions( dimensionNumber ) = M;    
end
biasFieldCoefficients = zeros( prod( numberOfBasisFunctions ), numberOfContrasts ); % No bias field to start with


if ( showFigures )
  posteriorFigure = figure;
  costFigure = figure;
  deformationMovieFigure = figure;
  biasFieldFigure = figure;
end


% We do the optimization in a multi-resolution type of scheme, where large
% deformations are quickly found using smoothed versions of the atlas mesh, and the fine
% details are then found on gradually less smoothed versions until the original atlas mesh is used for the optimization.
% Multi-resolution is a standard procedure in image registration: the initial
% blurring avoids getting stuck in the first local optimum of the search space, and get the rough major
% deformations right instead.
numberOfMultiResolutionLevels = length( optimizationOptions.multiResolutionSpecification );
historyWithinEachMultiResolutionLevel = struct( [] );
for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
    
  %  
  maximumNumberOfIterations = optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).maximumNumberOfIterations;
  historyOfCost = [ 1/eps ];
  historyOfMaximalDeformationApplied = [];
  historyOfTimeTakenIntensityParameterUpdating = [];
  historyOfTimeTakenDeformationUpdating = [];
  fprintf('maximumNumberOfIterations %d\n',maximumNumberOfIterations);

  % 
  kvlSetAlphasInMeshNodes( mesh, reducedAlphas )

  % Downsample the images, the mask, the mesh, and the bias field basis functions
  downSamplingFactors = max( round( optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).targetDownsampledVoxelSpacing ...
                                    ./ voxelSpacing ), [ 1 1 1 ] )
  downSampledMask = mask(  1 : downSamplingFactors( 1 ) : end, ...
                           1 : downSamplingFactors( 2 ) : end, ...
                           1 : downSamplingFactors( 3 ) : end );
  downSampledMaskIndices = find( downSampledMask );                         
  downSampledImageBuffers = [];
  for contrastNumber = 1 : numberOfContrasts
    % if ( multiResolutionLevel == numberOfMultiResolutionLevels )
    if true
      % No image smoothing
      downSampledImageBuffers( :, :, :, contrastNumber ) = imageBuffers( 1 : downSamplingFactors( 1 ) : end, ...
                                                                         1 : downSamplingFactors( 2 ) : end, ...
                                                                         1 : downSamplingFactors( 3 ) : end, ...
                                                                         contrastNumber );
    else
      % Try image smoothing
      buffer = imageBuffers( 1 : downSamplingFactors( 1 ) : end, ...
                             1 : downSamplingFactors( 2 ) : end, ...
                             1 : downSamplingFactors( 3 ) : end, ...
                             contrastNumber );
      smoothingSigmas = downSamplingFactors / 2 / sqrt( 2 * log( 2 ) ); % Variance chosen to approximately 
                                                                        % match normalized binomial filter 
                                                                        % (1/4, 1/2, 1/4) for downsampling 
                                                                        % factor of 2
      smoothingSigmas( find( downSamplingFactors == 1 ) ) = 0.0;
      smoothedBuffer = kvlSmoothImageBuffer( single( buffer ), smoothingSigmas );
      smoothedMask = kvlSmoothImageBuffer( single( downSampledMask ), smoothingSigmas );     
      downSampledImageBuffers( :, :, :, contrastNumber ) = downSampledMask .* ( smoothedBuffer ./ ( smoothedMask + eps ) );                                                                     
    end
    
  end
  kvlScaleMesh( mesh, 1./downSamplingFactors );
  downSampledKroneckerProductBasisFunctions = cell( 0, 0 );
  for dimensionNumber = 1 : 3
    A = kroneckerProductBasisFunctions{ dimensionNumber };
    downSampledKroneckerProductBasisFunctions{ dimensionNumber } = A( 1 : downSamplingFactors( dimensionNumber ) : end, : );
  end
  DIM = size( downSampledImageBuffers( :, :, :, 1 ) );
  
  
  % Smooth the mesh using a Gaussian kernel.
  smoothingSigmas = optimizationOptions.multiResolutionSpecification( multiResolutionLevel ).meshSmoothingSigma ./ ...
                    voxelSpacing ./ downSamplingFactors; % In voxels
  fprintf( 'Smoothing mesh with kernel size %f %f %f...', smoothingSigmas( 1 ), smoothingSigmas( 2 ), smoothingSigmas( 3 ) )
  kvlSmoothMesh( mesh, smoothingSigmas );
  smoothedReducedAlphas = kvlGetAlphasInMeshNodes( mesh );
  fprintf( 'done\n' )
  
  
  % Algorithm-wise, we're just estimating sets of parameters for one given data (MR scan) that is
  % known and fixed throughout. However, in terms of bias field correction it will be computationally
  % more efficient to pre-compute the bias field corrected version of the scan ("corrected" with
  % the current estimate of the bias field) once and pass that on to different routines instead of the
  % original data. 
  % For convenience (although potentially a recipe for future bug introduction), I'm also keeping a
  % vectorized form of that around -- this will be useful in various places in the EM-parts. So 
  % effectively I have two redundant variables "downSampledBiasCorrectedImageBuffers" and "biasCorrectedData"
  % that really just encode the variable "biasFieldCoefficients" and so need to be meticiously updated each time 
  % "biasFieldCoefficients" is updated (!)
  downSampledBiasCorrectedImageBuffers = zeros( [ DIM numberOfContrasts ] );
  biasCorrectedData = zeros( [ length( downSampledMaskIndices ) numberOfContrasts ] );
  for contrastNumber = 1 : numberOfContrasts
    downSampledBiasField = backprojectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
    tmp = downSampledImageBuffers( :, :, :, contrastNumber ) - downSampledBiasField .* downSampledMask;
    downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) = tmp;
    biasCorrectedData( :, contrastNumber ) = tmp( downSampledMaskIndices );
  end
  
  
  % Compute a color coded version of the atlas prior in the atlas's current pose, i.e., *before*
  % we start deforming. We'll use this just for visualization purposes
  if ( showFigures )
    oldColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, DIM ), reducedColors );
  end

  
  
  historyWithinEachIteration = struct( [] );
  priors = zeros( length( downSampledMaskIndices ), length( numberOfGaussiansPerClass ) );
  posteriors = zeros( length( downSampledMaskIndices ), sum( numberOfGaussiansPerClass ) ); % Gaussian mixture models burst out into
                                                                                            % individual Gaussian components
  data = zeros( [ length( downSampledMaskIndices ) numberOfContrasts ] );  % Easier to work with vector notation in the EM computations
  for contrastNumber = 1:numberOfContrasts
    tmp = reshape( downSampledImageBuffers( :, :, :, contrastNumber ), [ prod(DIM) 1 ] ); 
    data( :, contrastNumber ) = tmp( downSampledMaskIndices );
  end
  for iterationNumber = 1 : maximumNumberOfIterations
    
    %
    startTimeIntensityParameterUpdating = tic;
    
    % Part I: estimate Gaussian mean and variances as well as bias
    % field parameters using EM.

    % Get the priors at the current mesh position
    tmp = reshape( kvlRasterizeAtlasMesh( mesh, DIM ), [ prod( DIM ) length( numberOfGaussiansPerClass ) ] );
    priors( : ) = double( tmp( downSampledMaskIndices, : ) ) / 65535;
    
    if ( iterationNumber == 1 )
      historyWithinEachMultiResolutionLevel( multiResolutionLevel ).priorsAtStart = priors;                                                                   
    end
    

    % Start EM iterations. Initialize the parameters if this is the
    % first time ever you run this
    if ( ( multiResolutionLevel == 1 ) && ( iterationNumber == 1 ) )
      
      if ( numberOfContrasts > 1 ) 
        % Initialize multi-contrast
        meansPerClass = zeros( numberOfContrasts, numberOfClasses );
        sigmasPerClass = zeros( numberOfContrasts, numberOfContrasts, numberOfClasses );
        
        means = zeros( numberOfContrasts, sum( numberOfGaussiansPerClass ) );
        sigmas = zeros( numberOfContrasts, numberOfContrasts, sum( numberOfGaussiansPerClass ) );
        
        
        shift = 0;
        for classNumber = 1 : numberOfClasses  
          prior = priors(:,classNumber);
          mean = data' * prior / ( sum( prior ));
          
          X = data - repmat(mean',size(data,1),1);
          X = X.*sqrt(repmat(prior,[1 numberOfContrasts]));
                    
          sigma = X'*X/( sum( prior ));
          if modelSpecifications.useDiagonalCovarianceMatrices
            % Force diagonal covariance matrices
            sigma = diag( diag( sigma ) );
          end
          
          meansPerClass( :, classNumber ) = mean;
          sigmasPerClass(:,:,classNumber) = sigma;
        end
        
        % Set the means to different values within a classNumber
        shift = 0;
        for classNumber = 1:numberOfClasses
          for i = 1:floor(numberOfGaussiansPerClass(classNumber)/2)
            meansSpread(i) = (floor(numberOfGaussiansPerClass(classNumber)/2))+(i-1);
          end
          
          if(mod(numberOfGaussiansPerClass(classNumber),2)==0)
            meansSpread = [meansSpread -meansSpread];
          elseif(mod(numberOfGaussiansPerClass(classNumber),2)~=0 && numberOfGaussiansPerClass(classNumber)~=1)
            meansSpread = [meansSpread 0 -meansSpread];
          elseif(numberOfGaussiansPerClass(classNumber)==1)
            meansSpread=0;
          end
          
          means(:,shift+classNumber : shift+classNumber+numberOfGaussiansPerClass(classNumber)-1) = (repmat(meansPerClass(:,classNumber),[1 numberOfGaussiansPerClass(classNumber)])) + ...
              (repmat(meansSpread,[numberOfContrasts 1]).*repmat(sqrt(diag(sigmasPerClass(:,:,classNumber))),[1 numberOfGaussiansPerClass(classNumber)]));
          
          
          shift=shift+numberOfGaussiansPerClass(classNumber)-1;
          clear meansSpread;
        end
        
        shift = 0;
        for classNumber = 1:numberOfClasses
          sigmas(:,:,classNumber+shift : classNumber+shift+numberOfGaussiansPerClass(classNumber) - 1) = repmat(sigmasPerClass(:,:,classNumber),[1 1 numberOfGaussiansPerClass(classNumber)]);
          shift = shift + numberOfGaussiansPerClass(classNumber) -1;
        end
        clear meansPerClass variancesPerClass meansPerImage variancesPerImage;
        
        %Finally compute a prior sigma which is used to get rid of
        %numerical singularities
        priorSigma = eye(numberOfContrasts).*cov(data);
        
      else 
        % Initialize uni-contrast
        means = zeros(sum(numberOfGaussiansPerClass), 1);
        variances = zeros(sum(numberOfGaussiansPerClass), 1);
        
        meansPerClass = zeros( numberOfClasses, 1 );
        variancesPerClass = zeros( numberOfClasses, 1);
        
        
        shift = 0;
        for classNumber = 1 : numberOfClasses  %loop over the classNumberes we want to extract
          
          prior = priors(:,classNumber);
          meansPerClass(classNumber) = data' * prior / ( sum( prior ));
          variancesPerClass(classNumber) = ((data - meansPerClass(classNumber)).^2 )' * prior / ( sum( prior ));
          shift=shift+numberOfGaussiansPerClass(classNumber) -1;
          
        end
                
        % Set the means to different values within a classNumber
        shift = 0;
        for classNumber = 1:numberOfClasses
          
          for i = 1:floor(numberOfGaussiansPerClass(classNumber)/2)
            meansSpread(i) = -(floor(numberOfGaussiansPerClass(classNumber)/2))+(i-1);
          end
          
          if(mod(numberOfGaussiansPerClass(classNumber),2)==0)
            meansSpread = [meansSpread -meansSpread];
          elseif(mod(numberOfGaussiansPerClass(classNumber),2)~=0 && numberOfGaussiansPerClass(classNumber)~=1)
            meansSpread = [meansSpread 0 -meansSpread];
          elseif(numberOfGaussiansPerClass(classNumber)==1)
            meansSpread=0;
          end
                    
                    
          means(shift+classNumber : shift+classNumber+numberOfGaussiansPerClass(classNumber)-1) = (meansPerClass(classNumber)) + (meansSpread.*sqrt(variancesPerClass(classNumber)));
                    
                    
          shift=shift+numberOfGaussiansPerClass(classNumber)-1;
          clear meansSpread;
        end
        
        shift = 0;
        for classNumber = 1:numberOfClasses
          
          variances(classNumber+shift : classNumber+shift+numberOfGaussiansPerClass(classNumber) - 1) = variancesPerClass(classNumber);
          shift = shift + numberOfGaussiansPerClass(classNumber) -1;
        end
        
        clear meansPerClass;
        clear variancesPerClass;
        
        % Finally compute a prior variance which is used to get rid of
        % numerical singularities
        priorVar = var( data );

      end % End test multi vs. unicontrast
      
    end % End test need for initialization

    % stopCriterionEM = 1e-5;
    historyOfEMCost = [ 1/eps ];
    for EMIterationNumber = 1 : 100
      %
      % E-step: compute the posteriors based on the current parameters.
      %
      shift = 0;
      for classNumber = 1 : numberOfClasses
        prior = priors( :, classNumber );
        for gaussian = 1 : numberOfGaussiansPerClass( classNumber )
          if ( numberOfContrasts > 1 )
            mean = means( :, classNumber + shift + gaussian - 1 );
            sigma = sigmas( :, :, classNumber + shift + gaussian - 1 );
            posteriors( :, classNumber + shift + gaussian - 1 ) = ...
                EMweights( classNumber + shift + gaussian - 1 ) .* mvnpdf( biasCorrectedData, mean', sigma ) .* prior;
          else
            mean = means( classNumber + shift + gaussian - 1 );
            variance = variances( classNumber + shift + gaussian - 1 );
            posteriors( :, classNumber + shift + gaussian - 1 ) = ...
                EMweights( classNumber + shift + gaussian - 1 ) .* normpdf( biasCorrectedData, mean, sqrt(variance) ) .* prior;
          end
        end
        
        shift = shift + numberOfGaussiansPerClass( classNumber ) - 1;
      end
      normalizer = sum( posteriors, 2 ) + eps;
      posteriors = posteriors ./ repmat( normalizer, [ 1 sum( numberOfGaussiansPerClass ) ] );
      minLogLikelihood = -sum( log( normalizer ) ); 
      historyOfEMCost = [ historyOfEMCost; minLogLikelihood ];
      
      % Show some figures
      if ( showFigures )
        shift = 0;
        for classNumber = 1 : numberOfClasses
          % if ( numberOfGaussiansPerClass( classNumber ) == 1 )
          %  continue;
          % end
          posterior = zeros( DIM );
          
          posterior( downSampledMaskIndices ) = ...
                   sum( posteriors( :, classNumber + shift : classNumber + shift + numberOfGaussiansPerClass ( classNumber ) - 1 ), ...
                        2 );
          shift = shift + numberOfGaussiansPerClass( classNumber ) - 1;
          figure( posteriorFigure )
          subplot( 2, 4, classNumber )
          showImage( posterior )
        end
        clear posterior
        
        figure( costFigure )
        subplot( 2, 1, 1 )
        plot( historyOfEMCost( 2 : end ) )
        title( 'EM cost' )
        subplot(2, 1, 2 )
        plot( historyOfCost( 2 : end ) )
        title( 'Cost' )
        
        figure( biasFieldFigure )
        for contrastNumber = 1 : numberOfContrasts
          subplot( numberOfContrasts, 2, ( contrastNumber - 1 ) * numberOfContrasts + 1 )
          showImage( exp( downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) / 1000 ) );
          subplot( numberOfContrasts, 2, ( contrastNumber - 1 ) * numberOfContrasts + 2 )
          downSampledBiasField = backprojectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, ...
                                                                            biasFieldCoefficients( :, contrastNumber ) );
          showImage( exp( downSampledBiasField / 1000 ) .* downSampledMask )
        end
        drawnow

      end % End test if we need to show some figures

      
      % Check for convergence
      % relativeChangeCost = ( historyOfEMCost(end-1) - historyOfEMCost(end) ) /  historyOfEMCost(end)
      % if ( relativeChangeCost < stopCriterionEM )
      changeCostPerVoxel = ( historyOfEMCost(end-1) - historyOfEMCost(end) ) / length( downSampledMaskIndices );
      if ( changeCostPerVoxel < optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion )
        % Converged
        disp( 'EM converged!' )
        break;
      end 
      

           
           
      %
      % M-step: update the model parameters based on the current
      % posterior estimate.
      %
      
      % First the weights
      shift = 0;
      for classNumber = 1 : numberOfClasses
        if ( numberOfGaussiansPerClass( classNumber ) == 1 )
          continue;
        end
        
        weightNorm = sum( sum( posteriors( :, ...
                                           classNumber + shift : classNumber + shift + numberOfGaussiansPerClass( classNumber ) - 1 ), 2 ) );
        for gaussian = 1 : numberOfGaussiansPerClass( classNumber )
          EMweights( classNumber + shift + gaussian - 1 ) = sum( posteriors( :, classNumber + shift + gaussian - 1 ) ) / weightNorm;
        end
        shift = shift + numberOfGaussiansPerClass( classNumber ) - 1;
      end
      if ( optimizationOptions.verbose )
        EMweights
      end
      
      %Then the means and (co)variances
      shift = 0;
      for classNumber = 1 : numberOfClasses
        for gaussian = 1 : numberOfGaussiansPerClass( classNumber )
          posterior = posteriors( :, classNumber + shift + gaussian - 1 );
          mean = biasCorrectedData' * posterior ./ sum( posterior );
          
          if ( numberOfContrasts > 1 )
            X = biasCorrectedData - repmat( mean', size( biasCorrectedData, 1 ), 1 );
            X = X .* sqrt( repmat( posterior, [ 1 numberOfContrasts ] ) );
            % Update using the prior sigma
            sigma = ( X'*X + priorSigma ) / ( 2 * ( numberOfContrasts + 2 ) + sum( posterior ) );
            if modelSpecifications.useDiagonalCovarianceMatrices
              % Force diagonal covariance matrices
              sigma = diag( diag( sigma ) );
            end
            
            means( :, classNumber + shift + gaussian - 1 ) = mean;
            sigmas( :, :, classNumber + shift + gaussian -1 ) = sigma;
          else
            variance = ( ( ( ( biasCorrectedData - mean ).^2 )' * posterior ) + priorVar ) / ...
                       ( 2 * ( numberOfContrasts + 2 ) + sum( posterior ) );
            means( classNumber + shift + gaussian - 1 ) = mean;
            variances( classNumber + shift + gaussian - 1 ) = variance;
          end
        end
        shift = shift + numberOfGaussiansPerClass( classNumber ) - 1;

      end % End loop over classNumberes
      
            
      % Update parameters of the bias field model. I'm only doing this for the very first
      % multi-resolution level, where the atlas is really smooth; for subsequent multi-resolution
      % levels the bias field parameters are kept fixed to their values estimated at the very first
      % level. 
      %  if ( ( multiResolutionLevel == 1 ) && ( iterationNumber ~= 1 ) ) % Don't attempt bias field correction until 
      %                                                                    % decent mixture model parameters are available  
      if ( iterationNumber ~= 1 ) % Don't attempt bias field correction until 
                                  % decent mixture model parameters are available  

        % TODO: the fact that means and variances are not consistent across cases numberOfContrasts ~= 1
        % is driving me nuts: not only "variances" vs. "sigmas", but also dimensions! For now I'll work
        % around it - the "todoMeans" and "todoVariances" should in the future replace the means and 
        % variances/sigmas currently in use throughout the entire code
        todoMeans = zeros( sum( numberOfGaussiansPerClass ), numberOfContrasts );
        todoVariances = zeros( sum( numberOfGaussiansPerClass ), numberOfContrasts, numberOfContrasts );
        for classNumber = 1 : sum( numberOfGaussiansPerClass )
          if ( numberOfContrasts > 1 )
            todoMeans( classNumber, : ) = means( :, classNumber )';
            todoVariances( classNumber, :, : ) = sigmas( :, :, classNumber );
          else
            todoMeans( classNumber ) = means( classNumber );
            todoVariances( classNumber ) = variances( classNumber );
          end
        end
        
        % Real stuff starts here
        
        %
        % Bias field correction: implements Eq. 8 in the paper
        %
        %    Van Leemput, "Automated Model-based Bias Field Correction of MR Images of the Brain", IEEE TMI 1999
        %                       
        precisions = zeros( size( todoVariances ) );
        for classNumber = 1 : sum( numberOfGaussiansPerClass )
          precisions( classNumber, :, : ) = reshape( inv( squeeze( todoVariances( classNumber, :, : ) ) ), ...
                                                     [ 1 numberOfContrasts numberOfContrasts ] );
        end  
          
        lhs = zeros( prod( numberOfBasisFunctions ) * numberOfContrasts ); % left-hand side of linear system
        rhs = zeros( prod( numberOfBasisFunctions ) * numberOfContrasts, 1 ); % right-hand side of linear system
        weightsImageBuffer = zeros( DIM );
        tmpImageBuffer = zeros( DIM );
        for contrastNumber1 = 1 : numberOfContrasts
          tmp = zeros( size( data, 1 ), 1 );
          for contrastNumber2 = 1 : numberOfContrasts
            classSpecificWeights = posteriors .* repmat( squeeze( precisions( :, contrastNumber1, contrastNumber2 ) )', ...
                                                         [ size( posteriors, 1 ) 1 ] );
            weights = sum( classSpecificWeights, 2 );
            
            % Build up stuff needed for rhs
            predicted = sum( classSpecificWeights .* repmat( todoMeans( :, contrastNumber2 )', [ size( posteriors, 1 ) 1 ] ), 2 ) ...
                        ./ ( weights + eps );
            residue = data( :, contrastNumber2 ) - predicted;
            tmp = tmp + weights .* residue;
            
            % Fill in submatrix of lhs
            weightsImageBuffer( downSampledMaskIndices ) = weights;
            lhs( ( contrastNumber1 - 1 ) * prod( numberOfBasisFunctions ) + [ 1 : prod( numberOfBasisFunctions ) ], ...
                 ( contrastNumber2 - 1 ) * prod( numberOfBasisFunctions ) + [ 1 : prod( numberOfBasisFunctions ) ] ) = ... 
                  computePrecisionOfKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, weightsImageBuffer );
            
          end % End loop over contrastNumber2
          
          tmpImageBuffer( downSampledMaskIndices ) = tmp;
          rhs( ( contrastNumber1 - 1 ) * prod( numberOfBasisFunctions ) + [ 1 : prod( numberOfBasisFunctions ) ] ) = ...
                          projectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, tmpImageBuffer );
                          
        end % End loop over contrastNumber1
        
        % lhs = lhs + diag( 0.001 * diag( lhs ) );
        
        biasFieldCoefficients = reshape( lhs \ rhs, [ prod( numberOfBasisFunctions ) numberOfContrasts ] );
        for contrastNumber = 1 : numberOfContrasts
          downSampledBiasField = backprojectKroneckerProductBasisFunctions( downSampledKroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
          tmp = downSampledImageBuffers( :, :, :, contrastNumber ) - downSampledBiasField .* downSampledMask;
          downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) = tmp;
          biasCorrectedData( :, contrastNumber ) = tmp( downSampledMaskIndices );
        end

      end % End test if multiResolutionLevel == 1
    
      
    end % End EM iterations
    historyOfEMCost = historyOfEMCost( 2 : end );
    timeTakenIntensityParameterUpdating = toc( startTimeIntensityParameterUpdating );  
    historyOfTimeTakenIntensityParameterUpdating = [ historyOfTimeTakenIntensityParameterUpdating; ...
                                                     timeTakenIntensityParameterUpdating ];


    %
    % Part II: update the position of the mesh nodes for the current set of Gaussian parameters
    %
    
    %
    startTimeDeformationUpdating = tic;

    % Before calling the mesh node position optimizer, which at this stage can only handle Gaussian intensity distributions
    % rather than the Gaussian mixture models we have, let's burst out our small number of classes into a larger number
    % of individual Gaussians
    optimizerAlphas = zeros( size( reducedAlphas, 1 ), sum( numberOfGaussiansPerClass ), 'single' );
    shift = 0;
    for classNumber = 1 : numberOfClasses
      alpha = smoothedReducedAlphas( :, classNumber );
      for gaussian = 1 : numberOfGaussiansPerClass( classNumber )
        optimizerAlphas( :, shift + classNumber + gaussian - 1 ) = EMweights( shift + classNumber + gaussian - 1 ) .* alpha;
      end
      shift = shift + numberOfGaussiansPerClass( classNumber ) - 1;
    end
    kvlSetAlphasInMeshNodes( mesh, optimizerAlphas );
    

    % Create ITK images to pass on to the mesh node position optimizer
    if ( exist( 'downSampledBiasCorrectedImages' ) == 1 )
      % Clean up mess from any previous iteration
      for contrastNumber = 1 : numberOfContrasts
        kvlClear( downSampledBiasCorrectedImages( contrastNumber ) );
      end
    end
    for contrastNumber = 1 : numberOfContrasts
      downSampledBiasCorrectedImages( contrastNumber ) = ...
             kvlCreateImage( single( downSampledBiasCorrectedImageBuffers( :, :, :, contrastNumber ) ) );      
    end

    
    if ( numberOfContrasts > 1 )
      precisions = zeros(size(sigmas));
      for classNumber = 1 : sum(numberOfGaussiansPerClass)
          precisions(:,:,classNumber) = inv(sigmas(:,:,classNumber));
      end
        
    else
      precisions = zeros(size(variances));
      for classNumber = 1 : sum(numberOfGaussiansPerClass)
        precisions(classNumber) = 1./variances(classNumber);            
      end
    end
    
    % TODO: the fact that means and variances are not consistent across cases numberOfContrasts ~= 1
    % is driving me nuts: not only "variances" vs. "sigmas", but also dimensions! For now I'll work
    % around it - the "todoMeans" and "todoVariances" should in the future replace the means and 
    % variances/sigmas currently in use throughout the entire code
    todoMeans = zeros( sum( numberOfGaussiansPerClass ), numberOfContrasts );
    todoVariances = zeros( sum( numberOfGaussiansPerClass ), numberOfContrasts, numberOfContrasts );
    for classNumber = 1 : sum( numberOfGaussiansPerClass )
      if ( numberOfContrasts > 1 )
        todoMeans( classNumber, : ) = means( :, classNumber )';
        todoVariances( classNumber, :, : ) = sigmas( :, :, classNumber );
      else
        todoMeans( classNumber ) = means( classNumber );
        todoVariances( classNumber ) = variances( classNumber );
      end
    end
    
    % Real stuff starts here
    precisions = zeros( size( todoVariances ) );
    for classNumber = 1 : sum( numberOfGaussiansPerClass )
      precisions( classNumber, :, : ) = reshape( inv( squeeze( todoVariances( classNumber, :, : ) ) ), ...
                                                  [ 1 numberOfContrasts numberOfContrasts ] );
    end  
    calculator = kvlGetCostAndGradientCalculator( 'AtlasMeshToIntensityImage', ...
                                                   downSampledBiasCorrectedImages, ...
                                                   'Sliding', ...
                                                   transform, ...
                                                   todoMeans, ...
                                                   precisions );

    %optimizerType = 'ConjugateGradient';
    optimizerType = 'L-BFGS';
    optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
                                    'Verbose', optimizationOptions.verbose, ...
                                    'MaximalDeformationStopCriterion', optimizationOptions.maximalDeformationStopCriterion, ... 
                                    'LineSearchMaximalDeformationIntervalStopCriterion', ...
                                      optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion, ...
                                    'MaximumNumberOfIterations', optimizationOptions.maximumNumberOfDeformationIterations, ... 
                                    'BFGS-MaximumMemoryLength', optimizationOptions.BFGSMaximumMemoryLength );

    historyOfDeformationCost = [];
    historyOfMaximalDeformation = [];
    nodePositionsBeforeDeformation = kvlGetMeshNodePositions( mesh );
    deformationStartTime = tic;
    while true
      %
      stepStartTime = tic;
      [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer );
      disp( [ 'maximalDeformation ' num2str( maximalDeformation ) ' took ' num2str( toc( stepStartTime ) ) ' sec' ] )
      
      if ( maximalDeformation == 0 )
        break;
      end

      %
      historyOfDeformationCost = [ historyOfDeformationCost; minLogLikelihoodTimesPrior ];
      historyOfMaximalDeformation = [ historyOfMaximalDeformation; maximalDeformation ];
    
    end % End loop over iterations
    kvlClear( calculator );
    kvlClear( optimizer );
    % haveMoved = ( length( historyOfDeformationCost ) > 0 );
    nodePositionsAfterDeformation = kvlGetMeshNodePositions( mesh );
    maximalDeformationApplied = sqrt( max( sum( ...
                ( nodePositionsAfterDeformation - nodePositionsBeforeDeformation ).^2, 2 ) ) );
    disp( '==============================' )
    disp( [ 'iterationNumber: ' num2str( iterationNumber ) ] )
    disp( [ '    maximalDeformationApplied: ' num2str( maximalDeformationApplied ) ] )
    disp( [ '  ' num2str( toc( deformationStartTime ) ) ' sec' ] )
    disp( '==============================' )
    
    % Undo the artificial bursting out of mixture models into individual Gaussian components
    kvlSetAlphasInMeshNodes( mesh, smoothedReducedAlphas );

    
    % Show a little movie comparing before and after deformation so far...
    if ( showFigures )
      figure( deformationMovieFigure )
      kvlSetAlphasInMeshNodes( mesh, smoothedReducedAlphas );
      newColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, DIM ), reducedColors );
      
      set( deformationMovieFigure, 'position', get( 0, 'ScreenSize' ) );
      for i = 1 : 10
        priorVisualizationAlpha = 0.4;
        backgroundImage = exp( downSampledImageBuffers( :, :, :, 1 ) / 1000 );
        backgroundImage = backgroundImage - min( backgroundImage(:) );
        backgroundImage = backgroundImage / max( backgroundImage(:) );

        % showImage( oldColorCodedPriors )
        imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( backgroundImage, [ 1 1 1 3 ] ) + ...
                      priorVisualizationAlpha * oldColorCodedPriors;
        showImage( imageToShow )
        drawnow
        pause( 0.1 )
        % showImage( newColorCodedPriors )
        imageToShow = ( 1 - priorVisualizationAlpha ) * repmat( backgroundImage, [ 1 1 1 3 ] ) + ...
                      priorVisualizationAlpha * newColorCodedPriors;
        showImage( imageToShow )
        drawnow
        pause( 0.1 )
      end
    end
    
    % Keep track of the cost function we're optimizing
    historyOfCost = [ historyOfCost; minLogLikelihoodTimesPrior ];
    historyOfMaximalDeformationApplied = [ historyOfMaximalDeformationApplied; maximalDeformationApplied ];
    timeTakenDeformationUpdating = toc( startTimeDeformationUpdating );  
    historyOfTimeTakenDeformationUpdating = [ historyOfTimeTakenDeformationUpdating; ...
                                              timeTakenDeformationUpdating ];
    
    
    % Save something about how the estimation proceeded
    %historyWithinEachIteration( iterationNumber ).priors = priors;
    %historyWithinEachIteration( iterationNumber ).posteriors = posteriors;
    historyWithinEachIteration( iterationNumber ).historyOfEMCost = historyOfEMCost;
    historyWithinEachIteration( iterationNumber ).EMweights = EMweights;
    historyWithinEachIteration( iterationNumber ).means = means;
    if ( numberOfContrasts > 1 )
      historyWithinEachIteration( iterationNumber ).sigmas = sigmas;
    else
      historyWithinEachIteration( iterationNumber ).variances = variances;
    end
    historyWithinEachIteration( iterationNumber ).biasFieldCoefficients = biasFieldCoefficients;
    historyWithinEachIteration( iterationNumber ).historyOfDeformationCost = historyOfDeformationCost;
    historyWithinEachIteration( iterationNumber ).historyOfMaximalDeformation = historyOfMaximalDeformation;
    historyWithinEachIteration( iterationNumber ).maximalDeformationApplied = maximalDeformationApplied;

    % Determine if we should stop the overall iterations over the two set of parameters
    %  if ( ( ~haveMoved ) || ...
    %        ( ( ( historyOfCost( end-1 ) - historyOfCost( end ) ) / historyOfCost( end ) ) ...
    %          < relativeCostDecreaseStopCriterion ) || ...
    %        ( maximalDeformationApplied < maximalDeformationAppliedStopCriterion ) )
    if ( ( ( ( historyOfCost( end-1 ) - historyOfCost( end ) ) / length( downSampledMaskIndices ) ) ...
           < optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion ) ) % If EM converges in one iteration and mesh node optimization doesn't do anything
         
      % Converged
      break;
    end
    
  end % End looping over global iterations for this multiresolution level
  historyOfCost = historyOfCost( 2 : end );
    
  % Undo the down sampling on the mesh
  kvlScaleMesh( mesh, downSamplingFactors );    
    
  % Save something about how the estimation proceeded
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSamplingFactors = downSamplingFactors;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSampledImageBuffers = downSampledImageBuffers;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSampledMask = downSampledMask;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration = ...
                                                                      historyWithinEachIteration;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfCost = historyOfCost;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfMaximalDeformationApplied = ...
                                                                      historyOfMaximalDeformationApplied;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenIntensityParameterUpdating = ...
                                                                      historyOfTimeTakenIntensityParameterUpdating;  
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenDeformationUpdating = ...
                                                                      historyOfTimeTakenDeformationUpdating;
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).priorsAtEnd = priors;                                                                   
  historyWithinEachMultiResolutionLevel( multiResolutionLevel ).posteriorsAtEnd = posteriors;
    
end % End loop over multiresolution levels



% Save something about how the estimation proceeded
history.imageBuffers = imageBuffers;
history.mask = mask;
history.historyWithinEachMultiResolutionLevel = historyWithinEachMultiResolutionLevel;
eval( [ 'save ' savePath '/history.mat history -v7.3' ] );


% OK, now that all the parameters have been estimated, try to segment the original, full resolution image
% with all the original labels (instead of the reduced "super"-structure labels we created).

% Get bias field corrected images
DIM = size( imageBuffers( :, :, :, 1 ) );
biasCorrectedImageBuffers = zeros( [ DIM numberOfContrasts ] );
biasFields = zeros( [ DIM numberOfContrasts ] );
for contrastNumber = 1 : numberOfContrasts
  biasField = backprojectKroneckerProductBasisFunctions( kroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
  biasCorrectedImageBuffers( :, :, :, contrastNumber ) = ...
              imageBuffers( :, :, :, contrastNumber ) - biasField .* mask;
  biasFields( :, :, :, contrastNumber ) = biasField;
end


% Undo the collapsing of several structures into "super"-structures
kvlSetAlphasInMeshNodes( mesh, alphas )
numberOfClasses = size( alphas, 2 );


% Get the priors as dictated by the current mesh position
% should kvlRasterizeAtlasMesh be GPUed here?
data = reshape( biasCorrectedImageBuffers, [ prod( DIM ) numberOfContrasts ] ); 
priors = kvlRasterizeAtlasMesh( mesh, DIM );
priors = reshape( priors, [ prod( DIM ) numberOfClasses ] );  


% Ignore everything that's has zero intensity
priors = priors( maskIndices, : );
data = data( maskIndices, : );


% Calculate the posteriors
posteriors = zeros( size( priors ), 'double' );
shift = 0;
for classNumber = 1 : numberOfClasses
    
  reducedClass = reducingLookupTable( classNumber );
  if ( reducedClass == 1 )
    shift = 1;
  else
    shift = sum( numberOfGaussiansPerClass( 1 : reducedClass-1 ) ) + 1;
  end
  
  prior = single( priors( :, classNumber ) ) / 65535;
  for gaussian = 1 : numberOfGaussiansPerClass( reducedClass )
    if ( numberOfContrasts > 1 )
      mean = means( :, shift + gaussian - 1 );
      sigma = sigmas( :, :, shift + gaussian - 1 );
      EMweight = EMweights( shift + gaussian - 1 );
      posteriors( :, classNumber ) = posteriors( :, classNumber ) + mvnpdf( data, mean', sigma ) .* EMweight;
    else
      mean = means( shift + gaussian - 1 );
      variance = variances( shift + gaussian - 1 );
      EMweight = EMweights( shift + gaussian - 1 );
      posteriors( :, classNumber ) = posteriors( :, classNumber ) + normpdf( data, mean, sqrt(variance) ) .* EMweight;
    end
  end
  posteriors( :, classNumber ) = posteriors( :, classNumber ) .* prior;

end % End loop over classNumberes

normalizer = sum( posteriors, 2 ) + eps;
posteriors = posteriors ./ repmat(normalizer, [ 1 numberOfClasses ] );
posteriors = uint16( round( posteriors * 65535 ) );

% Write out FreeSurfer segmentation
seg = zeros( prod( DIM ), numberOfClasses, 'uint16' );
seg( maskIndices, : ) = posteriors;
seg = reshape( seg ,[ DIM numberOfClasses ] );
[ maxProbability maxIndices ] = max( seg, [], 4 );
% Put into FreeSurfer format
labeling = zeros( DIM, 'single' );
for label = 1 : length( FreeSurferLabels )
  ind = find( maxIndices == label );
  labeling( ind ) = FreeSurferLabels( label );
end
%  [pathstr,name,ext] = fileparts(imageFileNames{1});
fprintf( 'Writing out freesurfer seg\n' );
samseg_writeOutFreeSurferSeg( imageFileNames{1}, transformedTemplateFileName, labeling, savePath, 'segSubSpace.mgz' );


% write out the bias field and the bias corrected image:
for contrastNumber = 1 : numberOfContrasts
  [ dataPath, scanName, ext ] = fileparts( imageFileNames{ contrastNumber } );
  outputfile = [ scanName '_biasField.mgz' ];
  samseg_writeOutFreeSurferSeg( imageFileNames{ contrastNumber }, transformedTemplateFileName, ...
                                exp( biasFields( :, :, :, contrastNumber ) / 1000 ) .* mask, savePath, outputfile );
  
  outputfile = [ scanName '_biasCorrected.mgz' ];
  samseg_writeOutFreeSurferSeg( imageFileNames{ contrastNumber }, transformedTemplateFileName, ...
                                exp( biasCorrectedImageBuffers( :, :, :, contrastNumber ) / 1000 ), savePath, outputfile );
end
