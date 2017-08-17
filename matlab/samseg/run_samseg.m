function retval = run_samseg(varargin)
% Run with no arguments to get help
retval = 1;
mfileversion = '$Id$';
 
% This is a structure to handle reading in of command-line args.  When
% adding a new arg, create a new field with the default value, then
% add a case in parse_args(). If the arg needs to be checked, add a
% case in check_params. Maybe add something to dump_args to print out
% the arg value. For clarity, do not reference this structure
% outside of this mfile.
cmdargs.involfiles = '';
cmdargs.missingStructures = '';
cmdargs.outdir = '';
cmdargs.regmatfile = '';
cmdargs.nthreads = 1;
cmdargs.debug = 0;
cmdargs.exvivo = 0;
%% Print useage if there are no arguments %%
if(nargin == 0)
  print_usage(cmdargs)
  return;
end
%% Parse the command-line arguments %%
cmdargs = parse_args(cmdargs,varargin);
if(isempty(cmdargs)) return; end
cmdargs = check_params(cmdargs);
if(isempty(cmdargs)) return; end
dump_args(cmdargs);

fprintf('%s\n',mfileversion);
fprintf('Matlab version %s\n',version);

savePath = cmdargs.outdir;
numberOfThreads = cmdargs.nthreads;
RegMatFile = cmdargs.regmatfile;
exvivo = cmdargs.exvivo;

ninvolfiles = size(cmdargs.involfiles,1);
imageFileNames = cell(0,0);
for n = 1:ninvolfiles
  imageFileNames{n} = deblank(cmdargs.involfiles(n,:));
end

numMissing = size(cmdargs.missingStructures,1);
missingStructureSearchStrings = cell(0,0);
for n = 1:numMissing
  missingStructureSearchStrings{n} = deblank(cmdargs.missingStructures(n,:));
end

% Display input
imageFileNames
savePath
numberOfThreads
exvivo
missingStructureSearchStrings

% Create the output folder
if ~exist(savePath, 'dir')
  mkdir(savePath);
end

% set SAMSEG_DATA_DIR as an environment variable, eg,
% setenv SAMSEG_DATA_DIR /autofs/cluster/koen/koen/GEMSapplications/wholeBrain
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
if(isempty(RegMatFile))
  fprintf('entering registerAtlas\n');
  showFigures = false;
  worldToWorldTransformMatrix = samseg_registerAtlas( imageFileNames{ 1 }, ...
                                                      affineRegistrationMeshCollectionFileName, ...
                                                      affineRegistrationTemplateFileName, ...
                                                      savePath, ...
                                                      showFigures );
else
  fprintf('Not performing registration:\n');
  fprintf('  Loading reg file %s\n',RegMatFile);
  load(RegMatFile);
  fname = sprintf('%s/SPM12_6classes_30x30x30_template_coregistrationMatrices.mat',savePath);
  save(fname,'worldToWorldTransformMatrix','imageToImageTransformMatrix');
end

% For historical reasons the samsegment script figures out the affine transformation from
% a transformed MNI template (where before transformation this template defines the segmentation
% mesh atlas domain). This is a bit silly really, but for now let's just play along and make
% sure we generate it
[ origTemplate, origTemplateTransform ] = kvlReadImage( templateFileName );  
transformedTemplateFileName = sprintf( '%s/mni305_masked_autoCropped_coregistered.mgz', savePath );
kvlWriteImage( origTemplate, transformedTemplateFileName, ...
               kvlCreateTransform( double( worldToWorldTransformMatrix * kvlGetTransformMatrix( origTemplateTransform ) ) ) );

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
retval = 0;
return


%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function cmdargs = parse_args(cmdargs,varargin)

inputargs = varargin{1};
ninputargs = length(inputargs);

narg = 1;
while(narg <= ninputargs)
  
  flag = deblank(inputargs{narg});
  narg = narg + 1;
  if(cmdargs.debug) fprintf(1,'Argument: %s\n',flag); end
  if(~isstr(flag))
    flag
    fprintf(1,'ERROR: All Arguments must be a string\n');
    error;
  end
  
  switch(flag)
    
   case '--i',
    arg1check(flag,narg,ninputargs);
    cmdargs.involfiles = strvcat(cmdargs.involfiles,inputargs{narg});
    narg = narg + 1;
    
   case '--threads',
    arg1check(flag,narg,ninputargs);
    cmdargs.nthreads = sscanf(inputargs{narg},'%d');
    narg = narg + 1;
    
   case '--o',
    arg1check(flag,narg,ninputargs);
    cmdargs.outdir = inputargs{narg};
    narg = narg + 1;
    
   case '--regmat',
    arg1check(flag,narg,ninputargs);
    cmdargs.regmatfile = inputargs{narg};
    narg = narg + 1;
    
   case '--missing',
    arg1check(flag,narg,ninputargs);
    cmdargs.missingStructures = strvcat(cmdargs.missingStructures,inputargs{narg});
    narg = narg + 1;

   case '--debug',
    cmdargs.debug = 1;
   
   case '--exvivo',
    cmdargs.exvivo = 1;
    
   otherwise
    fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
    cmdargs = [];
    return;
    
  end % --- switch(flag) ----- %

end % while(narg <= ninputargs)

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function cmdargs = check_params(cmdargs)
  if(isempty(cmdargs.involfiles))
    fprintf('ERROR: must specify at least one input with --i\n');
    error;
  end
  if(isempty(cmdargs.outdir))
    fprintf('ERROR: must specify an output dir with --o\n');
    error;
  end
return;

%--------------------------------------------------%
%% Check that there is at least one more argument %%
function arg1check(flag,nflag,nmax)
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    error;
  end
return;
%--------------------------------------------------%
%% Check that there are at least two more arguments %%
function arg2check(flag,nflag,nmax)
  if(nflag > nmax-1 ) 
    fprintf(1,'ERROR: Flag %s needs two arguments',flag);
    error;
  end
return;

%------------- Print Usage ---------------------%
function print_usage(cmdargs)
  fprintf('USAGE:\n');
  fprintf('run_samseg\n');
  fprintf(' --o <outdir>          : output folder\n');
  fprintf(' --i <input>           : input volume (add --i input for more)\n');
  fprintf(' --threads <nthreads>  : number of threads\n');
  fprintf(' --regmat <regmat>     : regmat from previous run\n');
  fprintf(' --missing <struct>    : specify a missing structure\n');
  fprintf(' --exvivo              : run samseg exvivo');
return

%------------- Dump Args ---------------------%
function dump_args(cmdargs)
  fprintf(' outdir %s\n',cmdargs.outdir);
  for n = 1:size(cmdargs.involfiles)
    fprintf(' %d %s\n',n,cmdargs.involfiles(n,:));
  end
  fprintf('nthreads %d\n',cmdargs.nthreads);
  if(~isempty(cmdargs.regmatfile))
    fprintf('regmatfile %s\n',cmdargs.regmatfile)
  end
    
return
