function retval = run_samseg(varargin)
% Run with no arguments to get help
retval = 1;
mfileversion = '$Id$';
AvgDataDir = getenv('SAMSEG_DATA_DIR');
 
% This is a structure to handle reading in of command-line args.  When
% adding a new arg, create a new field with the default value, then
% add a case in parse_args(). If the arg needs to be checked, add a
% case in check_params. Maybe add something to dump_args to print out
% the arg value. For clarity, do not reference this structure
% outside of this mfile.
cmdargs.involfiles = '';
cmdargs.outdir = '';
cmdargs.regmatfile = '';
cmdargs.nthreads = 1;
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

fprintf('%s\n',mfileversion);
fprintf('Matlab version %s\n',version);

savePath = cmdargs.outdir;
nThreads = cmdargs.nthreads;
RegMatFile = cmdargs.regmatfile;

ninvolfiles = size(cmdargs.involfiles,1);
imageFileNames = cell(0,0);
for n = 1:ninvolfiles
  imageFileNames{n} = deblank(cmdargs.involfiles(n,:));
end

% Create the output folder
mkdirp(savePath);
imageFileName = deblank(cmdargs.involfiles(1,:));

samsegStartTime = tic;

fprintf('entering kvlClear\n');
kvlClear; % Clear all the wrapped C++ stuff
close all;

if(isempty(RegMatFile))
  fprintf('entering registerAtlas\n');
  % Switch on if you want to initialize the registration by matching
  % (translation) the centers  of gravity 
  initializeUsingCenterOfGravityAlignment = false;
  showFigures = false;
  samseg_registerAtlas
else
  fprintf('Not performing registration\n');
  fprintf('  Loading reg %s\n',RegMatFile);
  load(RegMatFile);
  fname = sprintf('%s/SPM12_6classes_30x30x30_template_coregistrationMatrices.mat',savePath);
  save(fname,'worldToWorldTransformMatrix','imageToImageTransformMatrix');
end

% set SAMSEG_DATA_DIR as an environment variable, eg,
%  setenv SAMSEG_DATA_DIR /autofs/cluster/koen/koen/GEMSapplications/wholeBrain
templateFileName = sprintf('%s/mni305_masked_autoCropped.mgz',AvgDataDir);
meshCollectionFileName = sprintf('%s/CurrentMeshCollection30New.txt.gz',AvgDataDir);
% This is bascially an LUT
compressionLookupTableFileName = sprintf('%s/namedCompressionLookupTable.txt',AvgDataDir);

% For historical reasons the samsegment script figures out the affine transformation from
% a transformed MNI template (where before transformation this template defines the segmentation
% mesh atlas domain). This is a bit silly really, but for now let's just play along and make
% sure we generate it
[ origTemplate, origTemplateTransform ] = kvlReadImage( templateFileName );  
transformedTemplateFileName = sprintf('%s/mni305_masked_autoCropped_coregistered.mgz',savePath);
kvlWriteImage( origTemplate, transformedTemplateFileName, ...
               kvlCreateTransform( double( worldToWorldTransformMatrix * kvlGetTransformMatrix( origTemplateTransform ) ) ) );





showFigures = false; % Set this to true if you want to see some figures during the run.

multiResolutionSpecification = struct( [] );
multiResolutionSpecification{ 1 }.meshSmoothingSigma = 2.0; % In mm
multiResolutionSpecification{ 1 }.targetDownsampledVoxelSpacing = 2.0; % In mm
multiResolutionSpecification{ 1 }.maximumNumberOfIterations = 100;
multiResolutionSpecification{ 2 }.meshSmoothingSigma = 0.0; % In mm
multiResolutionSpecification{ 2 }.targetDownsampledVoxelSpacing = 1.0; % In mm
multiResolutionSpecification{ 2 }.maximumNumberOfIterations = 100;

maximumNumberOfDeformationIterations = 20;
absoluteCostPerVoxelDecreaseStopCriterion = 1e-4;
verbose = 0;

maximalDeformationStopCriterion = 0.001; % Measured in pixels
lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion; % Idem
% relativeCostDecreaseStopCriterion = 1e-6;
maximalDeformationAppliedStopCriterion = 0.0;
BFGSMaximumMemoryLength = 12;
K = 0.1; % Stiffness of the mesh
brainMaskingSmoothingSigma = 3; % sqrt of the variance of a Gaussian blurring kernel 
brainMaskingThreshold = 0.01;

samsegment;

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
  fprintf(' --o outdir : output folder\n');
  fprintf(' --i input : input volume (add --i input for more)\n');
  fprintf(' --threads nthreads\n');
  fprintf(' --regmat regmat : regmat from previous run\n');
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
