function retval = run_samseg(varargin)
% Run with no arguments to get help
retval = 1;
 
% This is a structure to handle reading in of command-line args.  When
% adding a new arg, create a new field with the default value, then
% add a case in parse_args(). If the arg needs to be checked, add a
% case in check_params. Maybe add something to dump_args to print out
% the arg value. For clarity, do not reference this structure
% outside of this mfile.
cmdargs.involfiles = '';
cmdargs.outdir = '';
cmdargs.regmatfile = '';
cmdargs.initltafile = '';
cmdargs.regonly = 0;
cmdargs.nthreads = 1;
cmdargs.stiffness = 0.1;
cmdargs.showfigs = 0;
cmdargs.nobrainmask = 0;
cmdargs.diagcovs = 0;
cmdargs.debug = 0;


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

fprintf('Matlab version %s\n',version);
fprintf('2 cmd Stiffness %g\n',cmdargs.stiffness);
    
savePath = cmdargs.outdir;
numberOfThreads = cmdargs.nthreads;
showFigures = cmdargs.showfigs;
noBrainMasking = cmdargs.nobrainmask;
useDiagonalCovarianceMatrices = cmdargs.diagcovs;
RegMatFile = cmdargs.regmatfile;
InitLTAFile = cmdargs.initltafile;

ninvolfiles = size(cmdargs.involfiles,1);
imageFileNames = cell(0,0);
for n = 1:ninvolfiles
  imageFileNames{n} = deblank(cmdargs.involfiles(n,:));
end

% Display input
imageFileNames
savePath
numberOfThreads

% Create the output folder
if ~exist(savePath, 'dir')
  mkdir(savePath);
end

% set SAMSEG_DATA_DIR as an environment variable, eg,
% setenv SAMSEG_DATA_DIR /autofs/cluster/koen/koen/GEMSapplications/wholeBrain
samsegDataDir = getenv( 'SAMSEG_DATA_DIR' );


samsegStartTime = tic;

% Clean up
kvlClear; % Clear all the wrapped C++ stuff
close all;

% Specify the maximum number of threads the C++ stuff will use. The more threads you can use
% the faster, provided you have a matching amount of cores on your system - up to
% the point where memory becomes the bottle neck.
% If the following command is not provided, the number of cores on your system will be used
kvlSetMaximumNumberOfThreads( numberOfThreads );



% Affine registration
templateFileName = fullfile( samsegDataDir, 'template.nii' );
affineRegistrationMeshCollectionFileName = fullfile( samsegDataDir, 'atlasForAffineRegistration.txt.gz' );
worldToWorldTransformMatrix = [];
if ( ~isempty( RegMatFile ) )
  load( RegMatFile, 'worldToWorldTransformMatrix' );
end

[ worldToWorldTransformMatrix, transformedTemplateFileName ] = ...
          samseg_registerAtlas( imageFileNames{ 1 }, ...
                                affineRegistrationMeshCollectionFileName, ...
                                templateFileName, ...
                                savePath, ...
                                showFigures, ...
                                worldToWorldTransformMatrix,...
				InitLTAFile);
                                                    
if(cmdargs.regonly) 
  fprintf('Registration-only requested, so quiting now)');
  retval = 0;
  return;
end
               
% FreeSurfer (http://surfer.nmr.mgh.harvard.edu) has a standardized way of representation segmentations,
% both manual and automated, as images in which certain intensity levels correspond to well-defined
% anatomical structures - for instance an intensity value 17 always corresponds to the left hippocampus.
% The text file FreeSurferColorLUT.txt distributed with FreeSurfer contains all these definitions, as well 
% as a color for each structure in
% RGBA (Red-Green-Blue and Alpha (opaqueness)) format with which FreeSurfer will always represent segmented
% structures in its visualization tools
%
% Let's read the contents of the "compressionLookupTable.txt" file, and show the names of the structures
% being considered. The results are automatically sorted according to their "compressed label", i.e., the
% first result corresponds to the first entry in the vector of probabilities associated with each node in
% our atlas mesh.
compressionLookupTableFileName = fullfile( samsegDataDir, 'compressionLookupTable.txt' );
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
               

% Because we have many labels to segment, and each of these labels has its own Gaussian mixture model
% whose parameters (mean, variance, mixture weight) we have to estimate from the data, it may makes sense to restrict
% the degrees of freedom in the model somewhat by specifying that some of these labels have the same parameters
% governing their Gaussian mixture model. For example, we'd expect no intensity differences between the left and right 
% part of each structure.
% The way we implement this is by defining "super-structures" (i.e., a global white matter tissue class), and therefore
% work with a simplied ("reduced") model during the entire parameter estimation phase. At the same time we also build 
% an inverse lookup table (mapping from original class number onto a reduced class number (super-structure)) that we 
% will need to compute the final segmentation.
%  
sharedGMMParametersFileName = fullfile( samsegDataDir, 'sharedGMMParameters.txt' );
sharedGMMParameters = kvlReadSharedGMMParameters( sharedGMMParametersFileName );
               
               

% Set various model specifications
modelSpecifications = struct;
modelSpecifications.atlasFileName = fullfile( samsegDataDir, 'atlas_level2.txt.gz' );
modelSpecifications.FreeSurferLabels = FreeSurferLabels;
modelSpecifications.names = names;
modelSpecifications.colors = colors;
modelSpecifications.sharedGMMParameters = sharedGMMParameters;
modelSpecifications.useDiagonalCovarianceMatrices = useDiagonalCovarianceMatrices;
modelSpecifications.brainMaskingSmoothingSigma = 3; % sqrt of the variance of a Gaussian blurring kernel
if noBrainMasking
  modelSpecifications.brainMaskingThreshold = -Inf;
else
  modelSpecifications.brainMaskingThreshold = 0.01;
end
modelSpecifications.K = cmdargs.stiffness; % Stiffness of the mesh
fprintf('Stiffness %g\n',modelSpecifications.K);
modelSpecifications.biasFieldSmoothingKernelSize = 50.0;  % Distance in mm of sinc function center to first zero crossing

% Set various optimization options
optimizationOptions = struct;
optimizationOptions.multiResolutionSpecification = struct;
optimizationOptions.multiResolutionSpecification( 1 ).atlasFileName = fullfile( samsegDataDir, 'atlas_level1.txt.gz' );
optimizationOptions.multiResolutionSpecification( 1 ).targetDownsampledVoxelSpacing = 2.0; % In mm
optimizationOptions.multiResolutionSpecification( 1 ).maximumNumberOfIterations = 100;
optimizationOptions.multiResolutionSpecification( 1 ).estimateBiasField = true;
optimizationOptions.multiResolutionSpecification( 2 ).atlasFileName = fullfile( samsegDataDir, 'atlas_level2.txt.gz' );
optimizationOptions.multiResolutionSpecification( 2 ).targetDownsampledVoxelSpacing = 1.0; % In mm
optimizationOptions.multiResolutionSpecification( 2 ).maximumNumberOfIterations = 100;
optimizationOptions.multiResolutionSpecification( 2 ).estimateBiasField = true; % Switching this off will use the bias field estimated 
                                                                                % at lower resolution(s)
optimizationOptions.maximumNumberOfDeformationIterations = 20;
optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion = 1e-4;
optimizationOptions.verbose = 0;
optimizationOptions.maximalDeformationStopCriterion = 0.001; % Measured in pixels
optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion = optimizationOptions.maximalDeformationStopCriterion; % Idem
% optimizationOptions.relativeCostDecreaseStopCriterion = 1e-6;
optimizationOptions.maximalDeformationAppliedStopCriterion = 0.0;
optimizationOptions.BFGSMaximumMemoryLength = 12;


[ FreeSurferLabels, names, volumesInCubicMm ] = samsegment( imageFileNames, transformedTemplateFileName, ...
                                                            modelSpecifications, optimizationOptions, ...
                                                            savePath, showFigures );
names
volumesInCubicMm

fprintf('#@# run_samseg done %6.4f min\n',toc( samsegStartTime )/60);
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
    
   case '--stiffness',
    arg1check(flag,narg,ninputargs);
    cmdargs.stiffness = sscanf(inputargs{narg},'%f');
    fprintf('cmd Stiffness %g\n',cmdargs.stiffness);
    narg = narg + 1;
    
   case '--o',
    arg1check(flag,narg,ninputargs);
    cmdargs.outdir = inputargs{narg};
    narg = narg + 1;
    
   case '--regmat',
    arg1check(flag,narg,ninputargs);
    cmdargs.regmatfile = inputargs{narg};
    narg = narg + 1;
    
   case '--regonly',
    cmdargs.regonly = 1;
   case '--reg-only',
    cmdargs.regonly = 1;
    
   case '--initlta',
    arg1check(flag,narg,ninputargs);
    cmdargs.initltafile = inputargs{narg};
    narg = narg + 1;
    
   case '--showfigs',
    cmdargs.showfigs = 1;

   case '--nobrainmask'
    cmdargs.nobrainmask = 1;
  
   case 'diagcovs'
    cmdarg.diagcovs = 1;
    
   case '--debug',
    cmdargs.debug = 1;
   
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
  fprintf(' --showfigs            : show figures during run\n');
  fprintf(' --nobrainmask         : no initial brain masking based on affine atlas registration\n' );
  fprintf(' --diagcovs            : use diagonal covariance matrices (only affect multi-contrast case)\n' );
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
  fprintf('showfigs %d\n',cmdargs.showfigs );
  fprintf('nobrainmask %d\n',cmdargs.nobrainmask );
  fprintf('diagcovs %d\n',cmdargs.diagcovs );
    
return
