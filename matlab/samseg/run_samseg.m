function run_samseg( imageFileName1, savePath, nThreadsStr, UseGPUStr, imageFileName2, imageFileName3, imageFileName4, imageFileName5, imageFileName6, exvivoStr, missingStructureSearchString1, missingStructureSearchString2, missingStructureSearchString3, missingStructureSearchString4, missingStructureSearchString5, missingStructureSearchString6, missingStructureSearchString7, missingStructureSearchString8 )
% This function is a wrapper for running the samseg matlab
% scripts. This wrapper can be compiled (meaning that all the
% inputs are strings).
% 
% $Id: run_samseg.m,v 1.1 2017/01/26 00:22:47 greve Exp $
%

fprintf('Matlab version %s\n',version);


% Parse input -- filling in defaults if the entire long list of input strings is not provided
missingStructureSearchStrings = {};
for missingStructureNumber = 1 : 8
  if ( nargin >= ( 10 + missingStructureNumber ) )
    eval( [ 'missingStructureSearchString = missingStructureSearchString' num2str( missingStructureNumber ) ';' ] )
    missingStructureSearchStrings = [ missingStructureSearchStrings missingStructureSearchString ];
  end
end

if ( nargin < 10 )
  exvivo = 0;
else
  exvivo = str2num( exvivoStr );
end

imageFileNames = { imageFileName1 };
for imageFileNameNumber = 2 : 6
  if ( nargin >= ( 3 + imageFileNameNumber ) )
    eval( [ 'imageFileName = imageFileName' num2str( imageFileNameNumber ) ';' ] )
    imageFileNames = [ imageFileNames imageFileName ];
  end
end

if ( nargin < 3 )
  error( 'At least three input arguments expected' )
else
  numberOfThreads = str2num( nThreadsStr );
end



% Display input
imageFileNames
savePath
numberOfThreads
exvivo
missingStructureSearchStrings



% set SAMSEG_DATA_DIR as an environment variable, eg,
%  setenv SAMSEG_DATA_DIR /autofs/cluster/koen/koen/GEMSapplications/wholeBrain
AvgDataDir = getenv( 'SAMSEG_DATA_DIR' );
templateFileName = sprintf('%s/mni305_masked_autoCropped.mgz',AvgDataDir);
meshCollectionFileName = sprintf('%s/CurrentMeshCollection30New.txt.gz',AvgDataDir);
compressionLookupTableFileName = sprintf('%s/namedCompressionLookupTable.txt',AvgDataDir);
 
 
samsegStartTime = tic;

% Clean up
kvlClear; % Clear all the wrapped C++ stuff
close all;

% Specify the maximum number of threads the C++ stuff will use. The more threads you can use
% the faster, provided you have a matching amount of cores on your system - up to
% the point where memory becomes the bottle neck.
% If the following command is not provided, the number of cores on your system will be used
kvlSetMaximumNumberOfThreads( numberOfThreads );



% 
if exvivo

  % Specify which classes share the same intensity Gaussian mixture model (and the number of components within each model)
  sharedGMMParameters = struct;
  sharedGMMParameters( 1 ).mergedName = 'Unknown'; % Background and what is normally CSF
  sharedGMMParameters( 1 ).searchStrings = { 'Unknown', 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus' };
  sharedGMMParameters( 1 ).numberOfComponents = 1;
  sharedGMMParameters( 2 ).mergedName = 'Global WM'; % WM
  sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
  sharedGMMParameters( 2 ).numberOfComponents = 1;
  sharedGMMParameters( 3 ).mergedName = 'Global GM'; % GM
  sharedGMMParameters( 3 ).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities', 'Putamen' };
  sharedGMMParameters( 3 ).numberOfComponents = 1;
  sharedGMMParameters( 4 ).mergedName = 'Thalamus'; % Thalamus
  sharedGMMParameters( 4 ).searchStrings = { 'Thalamus' };
  sharedGMMParameters( 4 ).numberOfComponents = 1;
  sharedGMMParameters( 5 ).mergedName = 'Pallidum'; % Pallidum
  sharedGMMParameters( 5 ).searchStrings = { 'Pallidum' };
  sharedGMMParameters( 5 ).numberOfComponents = 1;

  buildTailoredAffineRegistrationAtlas = true;
  
else
  
  % Specify which classes share the same intensity Gaussian mixture model (and the number of components within each model)
  sharedGMMParameters = struct;
  sharedGMMParameters( 1 ).mergedName = 'Unknown'; % Background
  sharedGMMParameters( 1 ).searchStrings = { 'Unknown'};
  sharedGMMParameters( 1 ).numberOfComponents = 3;
  sharedGMMParameters( 2 ).mergedName = 'Global WM'; % WM
  sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
  sharedGMMParameters( 2 ).numberOfComponents = 2;
  sharedGMMParameters( 3 ).mergedName = 'Global GM'; % GM
  sharedGMMParameters( 3 ).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities' };
  sharedGMMParameters( 3 ).numberOfComponents = 3;
  sharedGMMParameters( 4 ).mergedName = 'Global CSF'; % CSF
  sharedGMMParameters( 4 ).searchStrings = { 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus' };
  sharedGMMParameters( 4 ).numberOfComponents = 3;
  sharedGMMParameters( 5 ).mergedName = 'Thalamus'; % Thalamus
  sharedGMMParameters( 5 ).searchStrings = { 'Thalamus' };
  sharedGMMParameters( 5 ).numberOfComponents = 2;
  sharedGMMParameters( 6 ).mergedName = 'Pallidum'; % Pallidum
  sharedGMMParameters( 6 ).searchStrings = { 'Pallidum' };
  sharedGMMParameters( 6 ).numberOfComponents = 2;
  sharedGMMParameters( 7 ).mergedName = 'Putamen'; % Putamen
  sharedGMMParameters( 7 ).searchStrings = { 'Putamen' };
  sharedGMMParameters( 7 ).numberOfComponents = 2;
  
  buildTailoredAffineRegistrationAtlas = false;

end


%
if buildTailoredAffineRegistrationAtlas
  % Create a tailor-made atlas for affine registration purposes

  % Read mesh 
  meshCollection = kvlReadMeshCollection( meshCollectionFileName );
  mesh = kvlGetMesh( meshCollection, -1 );
  [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );

  % Get a Matlab matrix containing a copy of the probability vectors in each mesh node (size numberOfNodes x
  % numberOfLabels ). 
  alphas = kvlGetAlphasInMeshNodes( mesh );

  % Remove non-existing structures
  mergeOptions = struct;
  mergeOptions( 1 ).mergedName = 'Unknown';
  mergeOptions( 1 ).searchStrings = missingStructureSearchStrings;
  [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );

  % Get global tissue types
  [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, sharedGMMParameters, FreeSurferLabels, colors );
  
  % Additionally move some deep gray matter structures into global GM
  mergeOptions = struct;
  mergeOptions( 1 ).mergedName = 'Global GM';
  mergeOptions( 1 ).searchStrings = { 'Thalamus', 'Pallidum', 'Putamen' };
  [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );
  
  % Create tailored atlas
  kvlSetAlphasInMeshNodes( mesh, alphas ); 
  [ template, transform ] = kvlReadImage( templateFileName ); 
  templateImageBuffer = kvlGetImageBuffer( template );
  priors = kvlRasterizeAtlasMesh( mesh, size( templateImageBuffer ) );
  transformMatrix = kvlGetTransformMatrix( transform );
  [ affineRegistrationMeshCollectionFileName, affineRegistrationTemplateFileName ] = ...
                                                createAtlasMeshForAffineRegistration( priors, transformMatrix, savePath );

else

  %
  affineRegistrationMeshCollectionFileName = sprintf( '%s/SPM12_6classes_30x30x30_meshCollection.txt.gz', AvgDataDir );
  affineRegistrationTemplateFileName = sprintf( '%s/SPM12_6classes_30x30x30_template.nii', AvgDataDir );


end % End test build tailored atlas for affine registration



%  
fprintf('entering registerAtlas\n');
if true
  showFigures = true;
  worldToWorldTransformMatrix = samseg_registerAtlas( imageFileNames{ 1 }, ...
                                                      affineRegistrationMeshCollectionFileName, ...
                                                      affineRegistrationTemplateFileName, ...
                                                      savePath, ...
                                                      showFigures );
else
  % No registration
  worldToWorldTransformMatrix = eye( 4 );
end



% For historical reasons the samsegment script figures out the affine transformation from
% a transformed MNI template (where before transformation this template defines the segmentation
% mesh atlas domain). This is a bit silly really, but for now let's just play along and make
% sure we generate it
[ origTemplate, origTemplateTransform ] = kvlReadImage( templateFileName );  
transformedTemplateFileName = sprintf( '%s/mni305_masked_autoCropped_coregistered.mgz', savePath );
kvlWriteImage( origTemplate, transformedTemplateFileName, ...
               kvlCreateTransform( double( worldToWorldTransformMatrix * kvlGetTransformMatrix( origTemplateTransform ) ) ) );


% 
fprintf('entering samsegment \n');

% Set various model specifications
modelSpecifications = struct;
modelSpecifications.missingStructureSearchStrings = missingStructureSearchStrings;
modelSpecifications.sharedGMMParameters = sharedGMMParameters;
modelSpecifications.useDiagonalCovarianceMatrices = false;
modelSpecifications.brainMaskingSmoothingSigma = 3; % sqrt of the variance of a Gaussian blurring kernel 
modelSpecifications.brainMaskingThreshold = 0.01;
modelSpecifications.K = 0.1; % Stiffness of the mesh
if exvivo
  modelSpecifications.brainMaskingThreshold = -Inf; % Disable brain masking
  modelSpecifications.useDiagonalCovarianceMatrices = true;
end


% Set various optimization options
optimizationOptions = struct;
optimizationOptions.multiResolutionSpecification = struct;
optimizationOptions.multiResolutionSpecification( 1 ).meshSmoothingSigma = 2.0; % In mm
optimizationOptions.multiResolutionSpecification( 1 ).targetDownsampledVoxelSpacing = 2.0; % In mm
optimizationOptions.multiResolutionSpecification( 1 ).maximumNumberOfIterations = 100;
optimizationOptions.multiResolutionSpecification( 2 ).meshSmoothingSigma = 0.0; % In mm
optimizationOptions.multiResolutionSpecification( 2 ).targetDownsampledVoxelSpacing = 1.0; % In mm
optimizationOptions.multiResolutionSpecification( 2 ).maximumNumberOfIterations = 100;
optimizationOptions.maximumNumberOfDeformationIterations = 20;
optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion = 1e-4;
optimizationOptions.verbose = 0;
optimizationOptions.maximalDeformationStopCriterion = 0.001; % Measured in pixels
optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion = optimizationOptions.maximalDeformationStopCriterion; % Idem
% optimizationOptions.relativeCostDecreaseStopCriterion = 1e-6;
optimizationOptions.maximalDeformationAppliedStopCriterion = 0.0;
optimizationOptions.BFGSMaximumMemoryLength = 12;

showFigures = false; % Set this to true if you want to see some figures during the run.
if exvivo
  showFigures = true;
end

samsegment( imageFileNames, transformedTemplateFileName, meshCollectionFileName, compressionLookupTableFileName, ...
            modelSpecifications, optimizationOptions, savePath, showFigures );

fprintf('#@# samseg done %6.4f min\n',toc( samsegStartTime )/60);

