function run_samseg(imageFileName1,savePath,nThreadsStr,imageFileName2,imageFileName3)
% run_samseg(imageFileName1,savePath,nThreadsStr,imageFileName2,imageFileName3)
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
fprintf('output path %s\n',savePath);
fprintf('nThreads = %d\n',nThreads);

if(strcmp(imageFileName2,'none')) imageFileName2 = ''; end
if(strcmp(imageFileName3,'none')) imageFileName3 = ''; end

imageFileName = imageFileName1;

% set SAMSEG_DATA_DIR as an environment variable, eg,
%  setenv SAMSEG_DATA_DIR /autofs/cluster/koen/koen/GEMSapplications/wholeBrain
AvgDataDir = getenv('SAMSEG_DATA_DIR');
templateFileName = sprintf('%s/mni305_masked_autoCropped.mgz',AvgDataDir);
meshCollectionFileName = sprintf('%s/CurrentMeshCollection30New.txt.gz',AvgDataDir);
% This is bascially an LUT
compressionLookupTableFileName = sprintf('%s/namedCompressionLookupTable.txt',AvgDataDir);
 
samsegStartTime = tic;

fprintf('entering kvlClear\n');
kvlClear; % Clear all the wrapped C++ stuff
close all;

fprintf('entering registerAtlas\n');
% Switch on if you want to initialize the registration by matching
% (translation) the centers  of gravity 
initializeUsingCenterOfGravityAlignment = false;
useSPMForAffineRegistration = true;
if useSPMForAffineRegistration
  samseg_registerToAtlas
else
  K = 1e-7; % Mesh stiffness -- compared to normal models, the entropy cost function is normalized 
            % (i.e., measures an average *per voxel*), so that this needs to be scaled down by the
            % number of voxels that are covered
  showFigures = true;
  samseg_registerAtlas
end
  

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
downSamplingFactor = 1;  % Use 1 for no downsampling
maxNuberOfIterationPerMultiResolutionLevel(1) = 5; % default 5
maxNuberOfIterationPerMultiResolutionLevel(2) = 20; % default 20
maximumNumberOfDeformationIterations = 500;
maximalDeformationStopCriterion = 0.001; % Measured in pixels
lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion; % Idem
meshSmoothingSigmas = [ 2.0 0 ]'; % UsemeshSmoothingSigmas = [ 0 ]' if you don't want to use multi-resolution
relativeCostDecreaseStopCriterion = 1e-6;
maximalDeformationAppliedStopCriterion = 0.0;
BFGSMaximumMemoryLength = 12;
K = 0.1; % Stiffness of the mesh
if useSPMForAffineRegistration
  brainMaskingSmoothingSigma = 2; % sqrt of the variance of a Gaussian blurring kernel 
  brainMaskingThreshold = 0.01;
else
  brainMaskingSmoothingSigma = 5; % 2; % sqrt of the variance of a Gaussian blurring kernel 
  brainMaskingThreshold = 0.01;
end


samsegment;

fprintf('#@# samseg done %6.4f min\n',toc( samsegStartTime )/60);

