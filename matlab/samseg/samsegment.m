% This is a different version of Oula's segment.m script. It has been
% modified so that inputs and some parameters must be set prior to the
% script being run (as in run_samseg.m). It also allows for "GPU"
% versions of some functions to be run. This does not actually make
% calls to a GPU, they are just functions that are GPU-able.
% 
% $Id: samsegment.m,v 1.1 2017/01/26 00:21:25 greve Exp $
%
% This script provides a basic idea of how to do whole-brain parcellation
% using tetrahedral mesh based atlases without having to deal with sophisticated
% C++ code.
%
% The overall algorithm implemented here is described very briefly in section VI.B of the paper
%
%      Encoding Probabilistic Brain Atlases Using Bayesian Inference
%      K. Van Leemput
%      IEEE Transactions on Medical Imaging, vol. 28, no. 6, pp. 822-837, June 2009
%
% except that there a gradient-descent optimization of the atlas deformation was used, and
% a much more powerful method is used here.
%
% The code is quite memory-hungry; in order to avoid hanging your system when
% Matlab tries to allocate more memory, limit the amount it is allowed to use. If it's
% not sufficient, Matlab will give an error but your system will never become non-responsive.
%
% On my system Linux system with 8G RAM, I use
%  ulimit -m 6500000
%  ulimit -v 8000000
%
% Before starting Matlab, you should also make sure the correct C++ library on your
% system is loaded, rather than the one that is shipped with Matlab:
%  export LD_PRELOAD=/usr/lib/libstdc++.so.6

% These are inputs that must be set prior to running
fprintf('input file %s\n',imageFileName);
fprintf('output path %s\n',savePath);
fprintf('nThreads = %d\n',nThreads);
fprintf('UseGPU = %d\n',UseGPU);
fprintf('downSamplingFactor = %d\n',downSamplingFactor);
fprintf('maxNuberOfIterationPerMultiResolutionLevel(1) = %d\n',maxNuberOfIterationPerMultiResolutionLevel(1));
fprintf('maxNuberOfIterationPerMultiResolutionLevel(2) = %d\n',maxNuberOfIterationPerMultiResolutionLevel(2));

% These are specified in run_samseg.m
%meshCollectionFileName = '/autofs/cluster/koen/koen/GEMSapplications/wholeBrain/CurrentMeshCollection30New.txt.gz';
%compressionLookupTableFileName = '/autofs/cluster/koen/koen/GEMSapplications/wholeBrain/namedCompressionLookupTable.txt'; % Look-up table belonging to the atlas (see below)

% Specify the maximum number of threads the C++ stuff will use. The more threads you can use
% the faster, provided you have a matching amount of cores on your system - up to
% the point where memory becomes the bottle neck.
% If the following command is not provided, the number of cores on your system will be used
kvlSetMaximumNumberOfThreads(nThreads);

% Clean up the Matlab work space
kvlClear % Clear all the wrapped C++ stuff
close all

%Set this to true if you want to see some figures during the run.
showImages=false;

% Provide the location of the image to be segmented, as well as the atlas that has been
% pre-registered affinely (i.e., 12 degrees of freedom) to the image.
transformedTemplateFileName = sprintf('%s/mni305_masked_autoCropped_coregistered.mgz',savePath);

numberOfImages = length(imageFileNames);

% Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
% translation, rotation, scaling, and skewing) as well - this transformation will later be used
% to initially transform the location of the atlas mesh's nodes into the coordinate system of
% the image.
%

for nima = 1:numberOfImages
  % Get the pointers to image and the corresponding transform
  [ images(nima), transform ] = kvlReadCroppedImage( imageFileNames{nima}, transformedTemplateFileName ); 
  imageBuffers(:,:,:,nima) = kvlGetImageBuffer( images(nima) ); %Get the actual imageBuffer
end

if(showImages)
  for nima = 1 : numberOfImages
    figure
    showImage( imageBuffers(:,:,:,nima) ); % Automatically displays middle slices in each direction
  end
end

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

K = 0.1;
meshCollection = kvlReadMeshCollection( meshCollectionFileName, transform, K );

% Retrieve the reference mesh, i.e., the mesh representing the average shape.
mesh = kvlGetMesh( meshCollection, -1 );

%Skull strip the images

labelNumber = 0;
sz = [size(imageBuffers,1) size(imageBuffers,2) size(imageBuffers,3)];
if(UseGPU ==0)
  backgroundPrior = kvlRasterizeAtlasMesh( mesh, sz, labelNumber );
else
  fprintf('Running kvlRasterizeAtlasMeshGPU()\n');
  backgroundPrior = kvlRasterizeAtlasMeshGPU( mesh, sz, labelNumber );
end

if(showImages)
  figure
  subplot( 2, 2, 1 )
  showImage( backgroundPrior )
  subplot( 2, 2, 2 )
  showImage( mosaicImages( 2^16 - 1 - double( backgroundPrior ), double( imageBuffers(:,:,:,1) ), 10 ) )
end

smoothingSigma = 2; % sqrt of the variance of a Gaussian blurring kernel
smoothedBackgroundPrior = kvlSmoothImageBuffer( backgroundPrior, smoothingSigma );

if(showImages)
    subplot( 2, 2, 3 )
    showImage( smoothedBackgroundPrior )
end

brainMask = ( 1 - single( smoothedBackgroundPrior ) / 65535 ) > .01;

for nima = 1 : numberOfImages
  imageBuffer = kvlGetImageBuffer( images(nima) );
  imageBuffer( find( ~brainMask ) ) = 0;
  kvlSetImageBuffer( images(nima), imageBuffer );
  imageBuffers(:,:,:,nima) = imageBuffer;
  clear imageBuffer;
end

if(showImages)
  subplot(2,2,4)
  showImage(imageBuffers(:,:,:,1))
end


% FreeSurfer (http://surfer.nmr.mgh.harvard.edu) has a standardized way of representation segmentations,
% both manual and automated, as images in which certain intensity levels correspond to well-defined
% anatomical structures - for instance an intensity value 17 always corresponds to the left hippocampus.
% The text file FreeSurferColorLUT.txt distributed with FreeSurfer (a copy of which is found in the
% Atlas3D source directory) contains all these definitions, as well as a color for each structure in
% RGBA (Red-Green-Blue and Alpha (opaqueness)) format with which FreeSurfer will always represent segmented
% structures in its visualization tools
%

%
% Let's read the contents of this "compressionLookupTable.txt" file, and show the names of the structures
% being considered. The results are automatically sorted according to their "compressed label", i.e., the
% first result corresponds to the first entry in the vector of probabilities associated with each node in
% our atlas mesh.
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
names


if(showImages)
  fprintf( 'Visualizing the atlas mesh; this takes quite some time and is only here for tutorial purposes...' )
  priors = kvlRasterizeAtlasMesh( mesh, [size(imageBuffers,1) size(imageBuffers,2) size(imageBuffers,3)]); % Without specifying a specific label, will rasterize all simultaneously, rasterize everything
  rgbBuffer = kvlColorCodeProbabilityImages( priors, colors );
  figure
  showImage( rgbBuffer )
  clear priors rgbBuffer
  fprintf( 'done\n' )
  drawnow;
end


% Although we now have about 41 labels to segment or so, and each of these labels has its own Gaussian distribution
% whose parameters (mean + variance) we have to estimate from the data, it may makes sense to restrict
% the degrees of freedom in the model somewhat by specifying that some of these labels have the same parameters
% governing their Gaussian. For example, we'd expect no intensity differences between the left and right part of
% each structure.
%


% Note: in order to know for sure what we're doing, the structures to merge into "super"-structures are specified
% by their FreeSurfer intensity values - you'd effectively have to look up in "compressionLookupTable.txt" what they mean.
sameGaussianParameters = cell(0,0);
sameGaussianParameters{1} = [ 0 ];  % Background is a separate class
sameGaussianParameters{2} = [ 2 7 46 16 41 28 60 85 ]; % Force all white matter structures to have the same intensity model
sameGaussianParameters{3} = [ 3 8 42 47 11 50 17 53 18 54 26 58 77 80 ]; % Same for several gray matter structures
sameGaussianParameters{4} = [ 4 5 14 24 15 43 44 72 30 62 31 63 ]; % Same for CSF
sameGaussianParameters{5} = [ 10 49 ]; % Force left and right thalamus to  have the same intensity model
sameGaussianParameters{6} = [ 12 51 ]; % Same for putamen
sameGaussianParameters{7} = [ 13 52 ]; % Same for pallidum

%Decides how many Gaussians you want to use to model each class, this set-up has typically worked for me.
numberOfGaussiansPerClass=[3 2 3 3 2 2 2]; 
numberOfClasses = length(numberOfGaussiansPerClass);

% Get a Matlab matrix containing a copy of the probability vectors in each mesh node (size numberOfNodes x
% numberOfLabels ).
originalAlphas = kvlGetAlphasInMeshNodes( mesh );

%Check that these sum to one.
max( abs( sum( originalAlphas, 2 ) - 1 ) )


% Compute the "reduced" alphas - those referring to the "super"-structures
[ reducedAlphas, reducingLookupTable ] = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );

% Set the reduced alphas to be the alphas of the mesh
kvlSetAlphasInMeshNodes( mesh, reducedAlphas );


%Initialize class weights, which are at first split evenly
EMweights=zeros(1,sum(numberOfGaussiansPerClass));

shift=0;
for class = 1:numberOfClasses
  EMweights(class + shift : class + shift + numberOfGaussiansPerClass(class) - 1) = 1/numberOfGaussiansPerClass(class);
  shift = shift + numberOfGaussiansPerClass(class)-1;
end

%You can choose to downsample the image to make the code even faster, however this makes the segmentations less accurate.
%This is generally not needed, because the algorithm runs quite fast
%anyway.
if(downSamplingFactor ~= 1)
  for nima = 1 : numberOfImages
    
    buffer = kvlGetImageBuffer( images(:,nima) );
    buffer = buffer( 1 : downSamplingFactor : end, ...
		     1 : downSamplingFactor : end, ...
		     1 : downSamplingFactor : end );
    fullResolutionImages(:,nima) = images(:,nima); % Save for later use
    images(:,nima) = kvlCreateImage( buffer );
    downSampledImageBuffers(:,:,:,nima) = buffer;
    clear buffer;
  end
  imageBuffers = downSampledImageBuffers;
  clear downSampledImageBuffers;
  kvlScaleMesh( mesh, 1/downSamplingFactor );
end


%Mask zeros out from the images.
buffer1 = kvlGetImageBuffer(images(1));
mask = (buffer1 > 0);
if(numberOfImages > 1)
    for nima = 2:numberOfImages
        bufferTemp = kvlGetImageBuffer(images(nima));
        maskIndicesTemp = (bufferTemp>0);
        mask = mask.*maskIndicesTemp;
    end
end
maskIndices = find(mask);
clear buffer1 bufferTemp maskIndicesTemp;

% Let's prepare for the bias field correction that is part of the imaging model. It assumes
% an additive effect, whereas the MR physics indicate it's a multiplicative one - so we log
% transform the data first
for nima = 1:numberOfImages
    buffer = kvlGetImageBuffer(images(nima));
    DIM = size(buffer);
    mask = zeros( DIM );
    mask( maskIndices ) = 1;
    buffer( maskIndices ) = 1000* log( buffer( maskIndices ) ); % The 1000 factor isn't really necessary; it's just
                                                                % easier to have Gaussian means of 100 350 than
                                                                % 0.1 and 0.35 for my human brain
    
    buffer = buffer .* mask;
    kvlSetImageBuffer( images(nima), buffer );
    imageBuffers(:,:,:,nima) = buffer;
    clear buffer
end

% Our bias model is a linear combination of a set of basis functions. We are using so-called
% "DCT-II" basis functions, i.e., the lowest few frequency components of the Discrete Cosine
% Transform.

distanceToFirstZeroCrossing = 50; % Larger kernels yield smoother bias field estimates
DIM = size(imageBuffers(:,:,:,1));
maximumNumberOfBasisFunctions = round( ( 2 * DIM(1) / distanceToFirstZeroCrossing + 1 ) / 2 ); 



% Algorithm-wise, we're just estimating sets of parameters for one given data (MR scan) that is
% known and fixed throughout. However, in terms of bias field correction it will be computationally
% much more efficient to pre-compute the bias field corrected version of the scan ("corrected" with
% the current estimate of the bias field) once and pass that on to different routines instead of the
% original data. Initially the bias field estimate is flat everywhere, so it's the same as the original
% image


for nima = 1:numberOfImages
    buffer = kvlGetImageBuffer(images(nima));
    biasCorrectedImages(nima) = kvlCreateImage( buffer );
end

% In our implementation, we gradually increase the number of basis functions used to represent the
% bias field, starting from very low-frequency components only (i.e., very wide smoothing kernels)
% and then moving slowly into more and more higher-frequency components (i.e., tighter kernels). This
% way massive global bias fields with one image area clearly much darker/brigher than the rest are
% removed quickly, before the fine details are filled in, helping avoid getting stuck in local optima

numberOfBasisFunctions = 1;

% We do the optimization in a multi-resolution type of scheme, where large
% deformations are quickly found using smoothed versions of the atlas mesh, and the fine
% details are then found on gradually less smoothed versions until the original atlas mesh is used for the optimization.
% Multi-resolution is a standard procedure in image registration: the initial
% blurring avoids getting stuck in the first local optimum of the search space, and get the rough major
% deformations right instead.
% Specify the size (in terms of voxels in the original image *before* downsampling) of the standard deviation
% of the Gaussian kernel used to smooth the priors/mesh. Use
%
%  meshSmoothingSigmas = [ 0 ]';
%
% if you don't want to use multi-resolution

meshSmoothingSigmas = [2.0 0]';

% Precalculate the spatially smoothed vector of prior probabilities in the current mesh.
numberOfMultiResolutionLevels = length( meshSmoothingSigmas );
smoothedReducedAlphas = zeros( size( reducedAlphas, 1 ), ...
    size( reducedAlphas, 2 ), ...
    numberOfMultiResolutionLevels, ...
    'single' );
for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
    
  % Smooth the mesh using a Gaussian kernel.
  meshSmoothingSigma = meshSmoothingSigmas( multiResolutionLevel ); 
  fprintf( 'Smoothing mesh with kernel size %f ...', meshSmoothingSigma)
  kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
  kvlSmoothMesh( mesh, meshSmoothingSigma);
  %kvlSmoothMesh( mesh, meshSmoothingSigma)
  fprintf( 'done\n' )
  
  % Store the result
  smoothedReducedAlphas( :, :, multiResolutionLevel ) = kvlGetAlphasInMeshNodes( mesh );
  
end


if(showImages)
  posteriorFigure = figure;
  costFigure = figure;
  deformationMovieFigure = figure;
  biasFieldFigure = figure;
  colorsTMP = zeros(length(numberOfGaussiansPerClass),4);
  
  colorsTMP(1,:) = [0 0 0 0];
  colorsTMP(2,:) = [0 255 0 255];
  colorsTMP(3,:) = [0 0 255 255];
  colorsTMP(4,:) = [255 0 0 255];
  colorsTMP(5,:) = [0 118 14 255];
  colorsTMP(6,:) = [236 13 176 255];
  colorsTMP(7,:) = [12 48 255 255];
  
end


for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
    
  % Set up the black box optimizer for the mesh nodes
  if ( exist( 'optimizer' ) == 1 )
    % The optimizer is very memory hungry when run in multithreaded mode.
    % Let's clear any old ones we may have lying around
    kvlClear( optimizer );
  end
  disp('Getting optimizer...');
  if(UseGPU==0)
    optimizer = kvlGetConjugateGradientOptimizer( mesh, biasCorrectedImages, transform );
  else
    fprintf('Running kvlGetConjugateGradientOptimizerGPU()\n');
    optimizer = kvlGetConjugateGradientOptimizerGPU( mesh, biasCorrectedImages, transform );
  end
  disp('Optimizer got.');
    
  maximumNumberOfIterations = maxNuberOfIterationPerMultiResolutionLevel(multiResolutionLevel);
  historyOfCost = [ 1/eps ];
  fprintf('maximumNumberOfIterations %d\n',maximumNumberOfIterations);
  
  % Compute a color coded version of the atlas prior in the atlas's current pose, i.e., *before*
  % we start deforming. We'll use this just for visualization purposes
  if(showImages)
    kvlSetAlphasInMeshNodes( mesh, smoothedReducedAlphas( :, :, multiResolutionLevel ) );
    oldColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, DIM ),colorsTMP );
  end
  
  
  for iterationNumber = 1 : maximumNumberOfIterations
    
    % Part I: estimate Gaussian mean and variances as well as bias
    % field parameters using EM.
    
    % Get the priors as dictated by the current mesh position, as well as the original and
    % the bias-field corrected intensities
    clear data;
    clear biasCorrectedData;
    data = zeros([prod(DIM) numberOfImages]);
    biasCorrectedData = zeros([prod(DIM) numberOfImages]);
    
    for nima = 1:numberOfImages
      data(:,nima) = reshape(double( kvlGetImageBuffer( images(nima) ) ), [prod(DIM) 1]); % Easier to work with vector notation in the computations
      biasCorrectedData(:,nima) = reshape(double(kvlGetImageBuffer( biasCorrectedImages(nima) ) ),[prod(DIM) 1]);
    end
    
    %Set the smoothed probability vectors into the mesh nodes
    kvlSetAlphasInMeshNodes( mesh, smoothedReducedAlphas( :, :, multiResolutionLevel ) );
    priors = kvlRasterizeAtlasMesh( mesh, DIM );
    priors = reshape( priors, [ prod( DIM ) length(numberOfGaussiansPerClass) ] );
    
    % Ignore everything that has zero intensity
    priors = priors( maskIndices, : );
    data = data( maskIndices, : );
    biasCorrectedData = biasCorrectedData( maskIndices, : );
    
    priors = double(priors) / 65535; %back to probabilities between 0 and 1
    
    
    
    % Start EM iterations. Initialize the parameters if this is the
    % first time ever you run this
    if ( ( multiResolutionLevel == 1 ) && ( iterationNumber == 1 ) )
      
      if(numberOfImages>1) %Initialize multi-contrast
	meansPerClass = zeros( numberOfImages, numberOfClasses );
	sigmasPerClass = zeros( numberOfImages, numberOfImages, numberOfClasses );
	
	means = zeros(numberOfImages, sum(numberOfGaussiansPerClass));
	sigmas = zeros(numberOfImages, numberOfImages, sum(numberOfGaussiansPerClass));
	
	
	shift = 0;
	for class = 1 : numberOfClasses  
	  prior = priors(:,class);
	  mean = biasCorrectedData' * prior / ( sum( prior ));
	  
	  X = biasCorrectedData - repmat(mean',size(biasCorrectedData,1),1);
	  X = X.*sqrt(repmat(prior,[1 numberOfImages]));
                    
	  sigma = X'*X/( sum( prior ));
	  meansPerClass( :, class ) = mean;
	  sigmasPerClass(:,:,class) = sigma;
	end
	
	%Set the means to different values within a class
	shift = 0;
	for class = 1:numberOfClasses
	  for i = 1:floor(numberOfGaussiansPerClass(class)/2)
	    meansSpread(i) = (floor(numberOfGaussiansPerClass(class)/2))+(i-1);
	  end
	  
	  if(mod(numberOfGaussiansPerClass(class),2)==0)
	    meansSpread = [meansSpread -meansSpread];
	  elseif(mod(numberOfGaussiansPerClass(class),2)~=0 && numberOfGaussiansPerClass(class)~=1)
	    meansSpread = [meansSpread 0 -meansSpread];
	  elseif(numberOfGaussiansPerClass(class)==1)
	    meansSpread=0;
	  end
	  
	  means(:,shift+class : shift+class+numberOfGaussiansPerClass(class)-1) = (repmat(meansPerClass(:,class),[1 numberOfGaussiansPerClass(class)])) + ...
	      (repmat(meansSpread,[numberOfImages 1]).*repmat(sqrt(diag(sigmasPerClass(:,:,class))),[1 numberOfGaussiansPerClass(class)]));
	  
	  
	  shift=shift+numberOfGaussiansPerClass(class)-1;
	  clear meansSpread;
	end
	
	shift = 0;
	for class = 1:numberOfClasses
	  sigmas(:,:,class+shift : class+shift+numberOfGaussiansPerClass(class) - 1) = repmat(sigmasPerClass(:,:,class),[1 1 numberOfGaussiansPerClass(class)]);
	  shift = shift + numberOfGaussiansPerClass(class) -1;
	end
	clear meansPerClass variancesPerClass meansPerImage variancesPerImage;
	
	%Finally compute a prior sigma which is used to get rid of
	%numerical singularities
	priorSigma = eye(numberOfImages).*cov(biasCorrectedData);
      else %initialize uni-contrast
	means = zeros(sum(numberOfGaussiansPerClass), 1);
	variances = zeros(sum(numberOfGaussiansPerClass), 1);
	
	meansPerClass = zeros( numberOfClasses, 1 );
	variancesPerClass = zeros( numberOfClasses, 1);
	
	
	shift = 0;
	for class = 1 : numberOfClasses  %loop over the classes we want to extract
	  
	  prior = priors(:,class);
	  meansPerClass(class) = biasCorrectedData' * prior / ( sum( prior ));
	  variancesPerClass(class) = ((biasCorrectedData - meansPerClass(class)).^2 )' * prior / ( sum( prior ));
	  shift=shift+numberOfGaussiansPerClass(class) -1;
	  
	end
                
	%Set the means to different values within a class
	
	shift = 0;
	
	for class = 1:numberOfClasses
	  
	  for i = 1:floor(numberOfGaussiansPerClass(class)/2)
	    meansSpread(i) = -(floor(numberOfGaussiansPerClass(class)/2))+(i-1);
	  end
	  
	  if(mod(numberOfGaussiansPerClass(class),2)==0)
	    meansSpread = [meansSpread -meansSpread];
	  elseif(mod(numberOfGaussiansPerClass(class),2)~=0 && numberOfGaussiansPerClass(class)~=1)
	    meansSpread = [meansSpread 0 -meansSpread];
	  elseif(numberOfGaussiansPerClass(class)==1)
	    meansSpread=0;
	  end
                    
                    
	  means(shift+class : shift+class+numberOfGaussiansPerClass(class)-1) = (meansPerClass(class)) + (meansSpread.*sqrt(variancesPerClass(class)));
                    
                    
	  shift=shift+numberOfGaussiansPerClass(class)-1;
	  clear meansSpread;
	end
	
	shift = 0;
	for class = 1:numberOfClasses
	  
	  variances(class+shift : class+shift+numberOfGaussiansPerClass(class) - 1) = variancesPerClass(class);
	  shift = shift + numberOfGaussiansPerClass(class) -1;
	end
	
	
	
	
	clear meansPerClass;
	clear variancesPerClass;
	
	%Finally compute a prior variance which is used to get rid of
	%numerical singularities
	priorVar = var(biasCorrectedData);
      end
      
      
    end % End test need for initialization
    
        
    stopCriterionEM = 1e-5;
    historyOfEMCost = [ 1/eps ];
    for EMIterationNumber = 1 : 100
      % E-step: compute the posteriors based on the current
      % parameters.
      
      shift = 0;
      for class = 1:numberOfClasses
	prior = priors(:,class);
	for gaussian = 1:numberOfGaussiansPerClass(class)
	  if(numberOfImages>1)
	    mean = means( :, class + shift + gaussian -1 );
	    sigma = sigmas( :, :, class + shift + gaussian -1 );
	    posteriors(:,class + shift + gaussian -1) = EMweights(class + shift + gaussian -1).*mvnpdf(biasCorrectedData, mean', sigma).*prior;
	  else
	    mean = means(class + shift + gaussian -1 );
	    variance = variances(class + shift + gaussian -1 );
	    posteriors(:,class + shift + gaussian -1) = EMweights(class + shift + gaussian -1).*normpdf(biasCorrectedData,mean,sqrt(variance)).*prior;
	  end
	end
	
	shift=shift+numberOfGaussiansPerClass(class)-1;
      end
      
      normalizer = sum( posteriors, 2 ) + eps;
      posteriors = posteriors ./ repmat( normalizer, [ 1 sum(numberOfGaussiansPerClass) ] );
            
            
      %
      % M-step: update the model parameters based on the current
      % posterior estimate.
      %
      
      %First the weights
      shift = 0;
      for class = 1:numberOfClasses
	if(numberOfGaussiansPerClass(class)==1)
	  continue;
	end
	
	weightNorm = sum(sum(posteriors(:,class+shift : class+shift+numberOfGaussiansPerClass(class) - 1),2));
	for gaussian = 1:numberOfGaussiansPerClass(class)
	  EMweights(class + shift + gaussian - 1) = (sum(posteriors(:,class+shift+gaussian-1)))/(weightNorm);
	end
	shift = shift + numberOfGaussiansPerClass(class) -1;
      end
      % Probably want to print this out only for verbose
      EMweights
      
      %Then the means and (co)variances
      shift = 0;
      for class = 1 : numberOfClasses
	for gaussian = 1 : numberOfGaussiansPerClass(class)
	  posterior = posteriors( :, class + shift + gaussian - 1);
	  mean = biasCorrectedData' * posterior ./ (sum( posterior ));
	  
	  if(numberOfImages>1)
	    X = biasCorrectedData - repmat(mean',size(biasCorrectedData,1),1);
	    X = X.*sqrt(repmat(posterior,[1 numberOfImages]));
	    %Update using the prior sigma
	    sigma = (X'*X + priorSigma)/(2*(numberOfImages + 2) + sum( posterior ));
	    
	    means(:,class + shift + gaussian - 1) = mean;
	    sigmas(:,:,class + shift + gaussian -1) = sigma;
	  else
	    variance = ((((biasCorrectedData - mean).^2)' * posterior) + priorVar) / (2*(numberOfImages + 2) + sum( posterior));
	    means(class + shift + gaussian -1) = mean;
	    variances(class + shift + gaussian -1) = variance;
	  end
	end
	shift = shift + numberOfGaussiansPerClass(class) -1;
      end
      
            
            % Update parameters of the bias field model. I'm only doing this for the very first
            % multi-resolution level, where the atlas is really smooth; for subsequent multi-resolution
            % levels the bias field parameters are kept fixed to their values estimated at the very first
            % level. 
            if ( multiResolutionLevel == 1)
                
                if(numberOfImages>1) 
                    
                    precisions = zeros(size(sigmas));
                    
                    for class = 1 : sum(numberOfGaussiansPerClass)
                        precisions(:,:,class) = inv(sigmas(:,:,class));
                    end
                    predictedImage = zeros( DIM(1), DIM(2), DIM(3), numberOfImages*numberOfImages );
                    weightImage = zeros( DIM(1), DIM(2), DIM(3), numberOfImages*numberOfImages  );
                    dataImage = zeros( DIM(1), DIM(2), DIM(3), numberOfImages*numberOfImages );
                    counter = 1;
                    for channel1 = 1:numberOfImages
                        for channel2 = 1:numberOfImages
                            tmp = posteriors .* repmat( squeeze(precisions(channel1,channel2,:))', [ size( posteriors, 1 ) 1 ] );
                            predicted = sum( tmp .* repmat( means(channel2,:), [ size( posteriors, 1 ) 1 ] ), 2 ) ./ ( sum(tmp,2) + eps );
                            predictedTMP = zeros(DIM);
                            predictedTMP(maskIndices) = predicted;
                            weightTMP = zeros(DIM);
                            weightTMP(maskIndices) = sum(tmp,2);
                            dataTMP = zeros(DIM);
                            dataTMP(maskIndices) = data(:,channel2);
                            
                            predictedImage(:,:,:,counter) =  predictedTMP;
                            weightImage(:,:,:,counter) = weightTMP;
                            dataImage(:,:,:, counter) = dataTMP;
                            
                            counter = counter + 1;
                        end
                    end
                    
                    
                    residueImage = dataImage - predictedImage;
                    
                    clear tmp predicted predictedTMP weightTMP dataTMP dataImage predictedImage;
                    
                    
                    if ( EMIterationNumber == 1 )
                        [ estimatedBiasField, coefficients ] = smoothWithSeparableBasisFunctionsWithWeightsMultiSpec( residueImage, weightImage, numberOfBasisFunctions, maximumNumberOfBasisFunctions);
                    else
                        [ estimatedBiasField, coefficients ] = smoothWithSeparableBasisFunctionsWithWeightsMultiSpec( residueImage, weightImage, numberOfBasisFunctions); % This is just a faster version of the above, but needs the above to have been
                        % called at least once...
                    end
                    
                    for nima = 1:numberOfImages
                        biasField = estimatedBiasField(:,:,:,nima);
                        biasCorrectedData(:,nima) = data(:,nima) - biasField(maskIndices); % Correct the data for the estimated bias field
                    end
                    
                else 
                    tmp = posteriors ./ repmat( variances', [ size( posteriors, 1 ) 1 ] );
                    weights = sum( tmp, 2 );
                    predicted = sum( tmp .* repmat( means', [ size( posteriors, 1 ) 1 ] ), 2 ) ./ ( weights + eps);
                    predictedImage = zeros( DIM );
                    predictedImage( maskIndices ) = predicted;
                    weightImage = zeros( DIM  );
                    weightImage( maskIndices ) = weights;
                    dataImage = zeros( DIM );
                    dataImage( maskIndices ) = data;
                    residueImage = dataImage - predictedImage; % The bias field is simply a (weighted) smoothed version of this
                    if ( EMIterationNumber == 1 )
                        [ estimatedBiasField, coefficients ] = smoothWithSeparableBasisFunctionsWithWeights( residueImage, weightImage, numberOfBasisFunctions, maximumNumberOfBasisFunctions );
                    else
                        [ estimatedBiasField, coefficients ] = smoothWithSeparableBasisFunctionsWithWeights( residueImage, weightImage, numberOfBasisFunctions ); % This is just a faster version of the above, but needs the above to have been
                        % called at least once...
                    end
                    biasCorrectedData = data - estimatedBiasField( maskIndices ); % Correct the data for the estimated bias field
                    
                end
            end
         
            minLogLikelihood = -sum( log( normalizer ) ); 
            historyOfEMCost = [ historyOfEMCost; minLogLikelihood ];
            
            
            % Check for convergence
            relativeChangeCost = ( historyOfEMCost(end-1) - historyOfEMCost(end) ) /  historyOfEMCost(end)
            
            
            if ( relativeChangeCost < 1e-4 )
                if ( numberOfBasisFunctions < maximumNumberOfBasisFunctions )
                    %We've converged, but we're not yet using the full bias field
                    %model. Increase flexibility of bias field
                    numberOfBasisFunctions = numberOfBasisFunctions + 1;
                    disp( [ 'Increasing numberOfBasisFunctions to ' num2str( numberOfBasisFunctions ) ] )
                    
                else
                    if ( relativeChangeCost < stopCriterionEM )
                        %Converged
                        disp( 'EM converged!' )
                        break;
                    end
                end
            end
            
   
            % Show some figures
            if(showImages)
                shift = 0;
                for class = 1:numberOfClasses
                    if(numberOfGaussiansPerClass(class)==1)
                        continue;
                    end
                    posterior = zeros( DIM );
                    
                    posterior(maskIndices) = sum(posteriors(:,class+shift : class+shift+numberOfGaussiansPerClass(class) - 1),2);
                    shift = shift + numberOfGaussiansPerClass(class) -1;
                    figure(posteriorFigure)
                    subplot(2,4,class)
                    showImage(posterior)
                end
                clear posterior
                
                figure( costFigure )
                subplot(2,1,1)
                plot( historyOfEMCost( 2 : end ) )
                title( 'EM cost' )
                subplot(2,1,2)
                plot( historyOfCost( 2 : end ) )
                title( 'Cost' )
                
                figure(biasFieldFigure)
                for n = 1:numberOfImages
                    subplot(numberOfImages,2,(n-1)*numberOfImages+1)
                    tmp = zeros(DIM);
                    tmp(maskIndices) = biasCorrectedData(:,n);
                    showImage(exp(tmp/1000));
                    subplot(numberOfImages,2,(n-1)*numberOfImages+2)
                    showImage(exp( estimatedBiasField(:,:,:,n) / 1000 ) .* mask)
                end
                drawnow
            end
        end
        
        
        
        %Before calling the optimizer, let's split the alphas in the mesh
        %nodes according to the EM weights.
        optimizerAlphas = zeros(size(reducedAlphas,1),sum(numberOfGaussiansPerClass),'single');
        shift = 0;
        for class = 1:numberOfClasses
            alpha = smoothedReducedAlphas(:,class,multiResolutionLevel);
            for gaussian = 1:numberOfGaussiansPerClass(class)
                optimizerAlphas(:,shift + class + gaussian -1) = EMweights(shift + class + gaussian -1).*alpha;
            end
            shift = shift + numberOfGaussiansPerClass(class) -1;
        end
        kvlSetAlphasInMeshNodes(mesh, optimizerAlphas);
        
        
        % Make sure the mesh node optimizer is made aware of the current bias field parameter estimate
        for nima = 1 : numberOfImages
            buffer = zeros(DIM,'single');
            buffer(maskIndices) = biasCorrectedData(:,nima);
            kvlSetImageBuffer( biasCorrectedImages(nima), buffer );
            clear buffer;
        end
        
        
        if(numberOfImages > 1)
            precisions = zeros(size(sigmas));
            for class = 1 : sum(numberOfGaussiansPerClass)
                precisions(:,:,class) = inv(sigmas(:,:,class));
            end
            
        else
            precisions = zeros(size(variances));
            for class = 1 : sum(numberOfGaussiansPerClass)
                precisions(class) = 1./variances(class);            
            end
        end
        
        %
        % Part II: update the position of the mesh nodes for the current set of Gaussian parameters
        %
        
        haveMoved = false; % Keep track if we've ever moved or not
        disp(['Setting up means and precisions']);
	if(UseGPU ==0) 
	  if(numberOfImages>1)
            kvlSetOptimizerProperties( optimizer, means', precisions );
	  else
            kvlSetOptimizerProperties( optimizer, means, precisions );
	  end
	else
	  if(numberOfImages>1)
	    fprintf('Running kvlSetOptimizerPropertiesGPU()\n');	    
            kvlSetOptimizerPropertiesGPU( optimizer, means', precisions );
	  else
	    fprintf('Running kvlSetOptimizerPropertiesGPU()\n');	    
            kvlSetOptimizerPropertiesGPU( optimizer, means, precisions );
	  end
	end
	
        disp(['Means and precisions set']);
        %disp(['means set']);
        for positionUpdatingIterationNumber = 1 : positionUpdatingMaximumNumberOfIterations
            % Calculate a good step. The first one is very slow because of various set-up issues
            tic
	    if(UseGPU ==0)
	      [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlDeformOneStep( optimizer );
	    else
	      fprintf('Running kvlDeformOneStepGPU()\n');	    
	      [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlDeformOneStepGPU( optimizer );
	    end
	    elapsedTime = toc;
            %disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
            %minLogLikelihoodTimesPrior
            if ( maximalDeformation > 0 )
                haveMoved = true;
            end
            
            % Test if we need to stop
            if ( maximalDeformation <= maximalDeformationStopCriterion )
                disp( 'maximalDeformation is too small; stopping' );
                break;
            end
        end
        
        
        % Show a little movie comparing before and after deformation so far...
        if(showImages)
            figure( deformationMovieFigure )
            kvlSetAlphasInMeshNodes( mesh, smoothedReducedAlphas( :, :, multiResolutionLevel ) );
            newColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, DIM ),colorsTMP );
            
            set( deformationMovieFigure, 'position', get(0,'ScreenSize') );
            for i=1:10
                showImage( oldColorCodedPriors )
                drawnow
                pause( 0.1 )
                showImage( newColorCodedPriors )
                drawnow
                pause( 0.1 )
            end
        end
        
        % Keep track of the cost function we're optimizing
        historyOfCost = [ historyOfCost; minLogLikelihoodTimesPrior ];
        
        
        % Determine if we should stop the overall iterations over the two set of parameters
        if ( ( ~haveMoved ) || ( ( ( historyOfCost( end-1 ) - historyOfCost( end ) ) / historyOfCost( end ) ) <  1e-6 ) )
            
            % Converged
            break;
        end
        
        
    end % End looping over global iterations
    
    
    
end % End loop over multiresolution levels




% OK, now that all the parameters have been estimated, try to segment the original, full resolution image
% with all the original labels (instead of the reduced "super"-structure labels we created).

% Clear some memory
kvlClear( optimizer )


for nima = 1 : numberOfImages
    biasCorrectedImageBuffers(:,:,:,nima) = kvlGetImageBuffer(biasCorrectedImages(nima));
end

% Undo the downsampling, if it was used.
if(downSamplingFactor ~= 1)
    
    for nima = 1 : numberOfImages
        buffer = kvlGetImageBuffer( fullResolutionImages(:,nima) );
        DIM = size( buffer );
        maskIndices = find( buffer > 0 );
        mask = zeros( DIM );
        mask( maskIndices ) = 1;
        buffer( maskIndices ) = 1000 * log( buffer( maskIndices ) );
        buffer = buffer .* mask;
        originalImageBuffers(:,:,:,nima) = buffer;
        clear buffer;
    end
    kvlScaleMesh( mesh, downSamplingFactor );
    
    biasCorrectedImageBuffers = zeros([DIM numberOfImages]);
    % Also reconstruct the estimated bias field at the original resolution
    disp( 'Reconstructing the bias field at the original resolution' )
    for nima = 1 : numberOfImages
        biasFields(:,:,:,nima) = expandBasisFunctions( coefficients(:,:,:,nima), ...
            DIM( 1 ), ...
            DIM( 2 ), ...
            DIM( 3 ), ...
            downSamplingFactor ) .* mask;
        
        %Correct the data for the estimated bias field
        biasCorrectedImageBuffers(:,:,:,nima) = originalImageBuffers(:,:,:,nima) - biasFields(:,:,:,nima);
    end
    
end


% Undo the collapsing of several structures into "super"-structures
kvlSetAlphasInMeshNodes( mesh, originalAlphas )
numberOfClasses = size( originalAlphas, 2 );


% Get the priors as dictated by the current mesh position
% should kvlRasterizeAtlasMesh be GPUed here?
data = double( reshape( biasCorrectedImageBuffers, [ prod( DIM ) numberOfImages ] ) ); 
if 1
    priors = kvlRasterizeAtlasMesh( mesh, DIM );
    priors = reshape( priors, [ prod( DIM ) numberOfClasses ] );  
else
    % This is the same as the above; it's slower but requires less memory
    priors = zeros( [ prod( DIM ) numberOfClasses ], 'uint16' ); 
    for classNumber = 1 : numberOfClasses
        prior = kvlRasterizeAtlasMesh( mesh, DIM, classNumber-1 ); % Remember! Non-Matlab programming start with index 0 rather than 1!
        priors( :, classNumber ) = prior(:);
    end
end


% Ignore everything that's has zero intensity
priors = priors( maskIndices, : );
data = data( maskIndices, : );





% Calculate the posteriors
posteriors = zeros( size(priors), 'double' );
shift = 0;
for class = 1:numberOfClasses
    
    reducedClass = reducingLookupTable(class);
    if(reducedClass==1)
        shift=1;
    else
        shift=sum(numberOfGaussiansPerClass(1:reducedClass-1))+1;
    end
    
    
    prior = single(priors(:,class))/65535;
    for gaussian = 1:numberOfGaussiansPerClass(reducedClass)
        if(numberOfImages>1)
            mean = means(:,shift + gaussian -1);
            sigma = sigmas(:,:,shift + gaussian -1);
            EMweight = EMweights(shift + gaussian -1);
            posteriors(:,class) = posteriors(:,class) + mvnpdf(data, mean', sigma).*EMweight;
        else
            mean = means(shift + gaussian -1);
            variance = variances(shift + gaussian -1);
            EMweight = EMweights(shift + gaussian -1);
            posteriors(:,class) = posteriors(:,class) + normpdf(data, mean, sqrt(variance)).*EMweight;
        end
    end
    posteriors(:,class) = posteriors(:,class).*prior;
end

normalizer = sum( posteriors, 2 ) + eps;
posteriors = posteriors ./ repmat(normalizer, [1 numberOfClasses]);
posteriors = uint16( round( posteriors * 65535 ) );

%Write out FreeSurfer segmentation
seg = zeros(prod(DIM),numberOfClasses,'uint16');
seg(maskIndices,:)=posteriors;
seg=reshape(seg ,[DIM numberOfClasses]);
[maxProbability maxIndices]=max(seg, [], 4);
%Put into FreeSurfer format
labeling = zeros(DIM,'single');

for label = 1:length(FreeSurferLabels)
    ind = find(maxIndices == label);
    labeling(ind) = FreeSurferLabels(label);
end
[pathstr,name,ext] = fileparts(imageFileNames{1});
fprintf('Writing out freesurfer seg\n');
samseg_writeOutFreeSurferSeg(imageFileNames{1},transformedTemplateFileName,labeling,savePath,'segSubSpace.mgz');
%writeOutFreeSurferSeg(imageFileNames{1},transformedTemplateFileName,labeling,pathstr);
%im = kvlReadCroppedImage(imageFileNames{1},transformedTemplateFileName);
%kvlSetImageBuffer(im,labeling);
%kvlWriteImage(im,[pathstr '/' 'seg.mgz']);

% write out the bias field and the bias corrected image:
for n = 1:numberOfImages
  [dataPath, scanName, ext] = fileparts(imageFileNames{n});
  outputfile = [scanName '_biasField.mgz'];
  samseg_writeOutFreeSurferSeg(imageFileNames{n},transformedTemplateFileName,exp(estimatedBiasField(:,:,:,n)/1000),savePath, outputfile);
  
  outputfile = [scanName '_biasCorrected.mgz'];
  samseg_writeOutFreeSurferSeg(imageFileNames{n},transformedTemplateFileName,exp(biasCorrectedImageBuffers(:,:,:,n)/1000),savePath, outputfile);
end