function run_samseg(imageFileName1,savePath,nThreadsStr,UseGPUStr,imageFileName2,imageFileName3,imageFileName4,imageFileName5,imageFileName6)
% This function is a wrapper for running the samseg matlab
% scripts. This wrapper can be compiled (meaning that all the
% inputs are strings).
% 
% $Id: run_samseg.m,v 1.1 2017/01/26 00:22:47 greve Exp $
%

fprintf('Matlab version %s\n',version);

nThreads = sscanf(nThreadsStr,'%d');
fprintf('input file1 %s\n',imageFileName1);
fprintf('input file2 %s\n',imageFileName2);
fprintf('input file3 %s\n',imageFileName3);
fprintf('input file4 %s\n',imageFileName4);
fprintf('input file5 %s\n',imageFileName5);
fprintf('input file6 %s\n',imageFileName6);
fprintf('output path %s\n',savePath);
fprintf('nThreads = %d\n',nThreads);

if(strcmp(imageFileName2,'none')) imageFileName2 = ''; end
if(strcmp(imageFileName3,'none')) imageFileName3 = ''; end
if(strcmp(imageFileName4,'none')) imageFileName4 = ''; end
if(strcmp(imageFileName5,'none')) imageFileName5 = ''; end
if(strcmp(imageFileName6,'none')) imageFileName6 = ''; end

imageFileName = imageFileName1;

AvgDataDir = getenv('SAMSEG_DATA_DIR');
 
samsegStartTime = tic;

fprintf('entering kvlClear\n');
kvlClear; % Clear all the wrapped C++ stuff
close all;

fprintf('entering registerAtlas\n');
% Switch on if you want to initialize the registration by matching
% (translation) the centers  of gravity 
initializeUsingCenterOfGravityAlignment = false;
showFigures = true;
samseg_registerAtlas


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



fprintf('entering samsegment \n');
%subject with different contrasts you can feed those in as well. Make sure
%though that the scans are coregistered
%imageFileNames{2} = '...'; %If you have more scans of the same
imageFileNames = cell(0,0);
imageFileNames{1} =  imageFileName1;
if(strlen(imageFileName2)>0) 
  imageFileNames{2} =  imageFileName2;
end
if(strlen(imageFileName3)>0) 
  imageFileNames{3} =  imageFileName3;
end
if(strlen(imageFileName4)>0) 
  imageFileNames{4} =  imageFileName4;
end
if(strlen(imageFileName5)>0) 
  imageFileNames{5} =  imageFileName5;
end
if(strlen(imageFileName6)>0) 
  imageFileNames{6} =  imageFileName6;
end


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

