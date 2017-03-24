function run_samseg(imageFileName1,savePath,nThreadsStr,UseGPUStr,imageFileName2,imageFileName3,imageFileName4,imageFileName5,imageFileName6)
% run_samseg(imageFileName1,savePath,nThreadsStr,UseGPUStr,imageFileName2,imageFileName3)
% This function is a wrapper for running the samseg matlab
% scripts. This wrapper can be compiled (meaning that all the
% inputs are strings).
% 
% $Id: run_samseg.m,v 1.1 2017/01/26 00:22:47 greve Exp $
%

fprintf('Matlab version %s\n',version);

nThreads = sscanf(nThreadsStr,'%d');
UseGPU = sscanf(UseGPUStr,'%d');
fprintf('input file1 %s\n',imageFileName1);
fprintf('input file2 %s\n',imageFileName2);
fprintf('input file3 %s\n',imageFileName3);
fprintf('input file4 %s\n',imageFileName4);
fprintf('input file5 %s\n',imageFileName5);
fprintf('input file6 %s\n',imageFileName6);
fprintf('output path %s\n',savePath);
fprintf('nThreads = %d\n',nThreads);
fprintf('UseGPU = %d\n',UseGPU);

if(strcmp(imageFileName2,'none')) imageFileName2 = ''; end
if(strcmp(imageFileName3,'none')) imageFileName3 = ''; end
if(strcmp(imageFileName2,'none')) imageFileName4 = ''; end
if(strcmp(imageFileName3,'none')) imageFileName5 = ''; end
if(strcmp(imageFileName3,'none')) imageFileName6 = ''; end

imageFileName = imageFileName1;

% set SAMSEG_DATA_DIR as an environment variable, eg,
%  setenv SAMSEG_DATA_DIR /autofs/cluster/koen/koen/GEMSapplications/wholeBrain
AvgDataDir = getenv('SAMSEG_DATA_DIR');
templateFileName = sprintf('%s/mni305_masked_autoCropped.mgz',AvgDataDir);
meshCollectionFileName = sprintf('%s/CurrentMeshCollection30New.txt.gz',AvgDataDir);
% This is bascially an LUT
compressionLookupTableFileName = sprintf('%s/namedCompressionLookupTable.txt',AvgDataDir);
 
tic;

fprintf('entering kvlClear\n');
kvlClear; % Clear all the wrapped C++ stuff
close all;

fprintf('entering registerToAtlas\n');
% Switch on if you want to initialize the registration by matching
% (translation) the centers  of gravity 
initializeUsingCenterOfGravityAlignment = 0; 
samseg_registerToAtlas

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
downSamplingFactor = 1;  % Use 1 for no downsampling
maxNuberOfIterationPerMultiResolutionLevel(1) = 5; % default 5
maxNuberOfIterationPerMultiResolutionLevel(2) = 20; % default 20
positionUpdatingMaximumNumberOfIterations = 500;
maximalDeformationStopCriterion = 1e-3;
samsegment;

fprintf('#@# samseg done %6.4f min\n',toc/60);