% SAXSlice.m
% Dec 7, 1998
% 
% Selective Averaging assuming uncorrelated noise.  This uses a
% matrix formuluation to solve the problem (as apposed to SelAvg.m
% which is simply selective averaging with correlations in the
% paradigm assumed to be cancelled by proper counter balancing).
%
% This is for multiple runs in a single session.
%
% This differs from SAX in that the inputs and par files
% are specified explicitly (not using SubjDir).  Also, the
% naming conventions are different in that this script requires
% a prefix for all names which is prepended to the mnSAX name.
% 
% fmri toolbox functions called:
% 1. LdBFile
% 2. SvBFile
% 3. HanKernel
% 4. fMRIPreProc
% 5. SelAvgX.m
%
% See: runSAXSlice
%
% global TimeWindow, TimeOffset, TR;
% global InputFiles, ParFiles;
% global hAvgFile, datFile, dofFile;
% global HanRadius;
% global tTPExclude
% global QuitOnError;
%
%


%
% fmri_saxslice.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

fprintf(1,'\n');
fprintf(1,'  --- SAXSlice: Starting ------\n');
fprintf(1,'fmri_saxslice.m @FS_VERSION@\n');

if( ~exist('QuitOnError') ) QuitOnError = 1; end

%%% ---- Check that all the variables are defined --- %%%
VarNameList = [];
VarNameList = strvcat(VarNameList,'TimeWindow');
VarNameList = strvcat(VarNameList,'TimeOffset');
VarNameList = strvcat(VarNameList,'RemoveBaseline');
VarNameList = strvcat(VarNameList,'RemoveTrend');
VarNameList = strvcat(VarNameList,'UsePercent');
VarNameList = strvcat(VarNameList,'TR');
VarNameList = strvcat(VarNameList,'nSkip');
VarNameList = strvcat(VarNameList,'RescaleTarget');
VarNameList = strvcat(VarNameList,'tPreStim');
VarNameList = strvcat(VarNameList,'InputFiles');
VarNameList = strvcat(VarNameList,'ParFiles');
VarNameList = strvcat(VarNameList,'TPExclFiles');
VarNameList = strvcat(VarNameList,'hAvgFile');
VarNameList = strvcat(VarNameList,'dofFile');
VarNameList = strvcat(VarNameList,'datFile');
VarNameList = strvcat(VarNameList,'VxStatFile');

nVar = size(VarNameList,1);
for n = 1:nVar,
  if( exist(deblank(VarNameList(n,:))) ~= 1)
    fprintf(2,'Error: Variable %s does not exist\n',VarNameList(n,:));
    if(QuitOnError) quit; 
    else return
    end
  end
end

% nRee = 20;
nRuns = size(InputFiles,1);

if( size(ParFiles,1) ~= nRuns )
  fprintf(2,'Incorrect number of ParFiles\n');
  if(QuitOnError) quit; 
  else            return;
  end
end

if( size(hAvgFile,1) ~= 1 )
  fprintf(2,'Incorrect number of hAvgFile\n');
  if(QuitOnError) quit; 
  else            return;
  end
end

nHEst = floor(TimeWindow/TR);
if(nHEst < 1)
  fprintf('TimeWindow too small, but be > TR\n');
  if(QuitOnError) quit; 
  else            return;
  end
end
nCond = [];

for n = 1:nRuns,
  fprintf('- Loading ParFile %d/%d %s\n',n,nRuns,ParFiles(n,:));
  Par(:,:,n) = LdPar(ParFiles(n,:));
  Par(:,1,n) = Par(:,1,n) + TimeOffset;
  nCond(n) = max(Par(:,2,n)) - min(Par(:,2,n))+1;
end

% Check that all runs have the same number of conditions:
if(length( find(diff(nCond)~=0)) ~= 0 )
  fprintf(2,'Error:All runs do not have the same number of conditions\n');
  fprintf(2,'nCond = %d\n',nCond);
  if(QuitOnError) quit; 
  else            return;
  end
end

fprintf(1,'  Found %d Conditions\n',nCond(1));

%%% ---- Load the fMRI Runs ------- %%%
for n = 1:nRuns,
    fprintf(1,'----- Loading %s ----\n',InputFiles(n,:));
    y = LdBFile(InputFiles(n,:));
    nRows = size(y,1);
    nCols = size(y,2);
    nTP   = size(y,3);
    nV    = nRows*nCols;
    Slice(:,:,:,n) = y;
end

%%% -- Spatial Filtering Setup-- %%
if( ~exist('HanRadius') ) HanFilter = [];
else
  if(HanRadius < 1)
    fprintf(1,'Error: HanRadius = %g, must be >= 1\n',HanRadius);
    if(QuitOnError) quit; 
    else            return;
    end
  end
  fprintf(1,'Using Spatial Filter, HanRad = %g\n',HanRadius);
  HanFilter = HanKernel(HanRadius);
end

%%%% --- Check for exclusions (including skips) ---%%%%
TPExclFiles
TPExclude = zeros(nTP,nRuns);
for n = 1:nRuns,
  tpexcl = ldtpexcl(TPExclFiles(n,:));
  if(nSkip > 0)
    tSkip = TR*[0:(nSkip-1)]';
    tpexcl = [tpexcl; tSkip;];
    tpexcl = unique(tpexcl);
  end
  fprintf('Run %2d Exclusions: ',n);
  fprintf('%2d ',tpexcl');
  fprintf('\n',tpexcl');
  
  l = length(tpexcl);
  if(l ~= 0)
    TPExclude(floor(tpexcl/TR)+1,n) = ones(l,1);
  end
end

%fprintf('------- Exclusion Vectors ------');
%TPExclude
%fprintf('------- oooooooooooooooo ------');

%%%% ----- Implement Skipping ------ %%%%%
%if(nSkip > 0)
%   fprintf(1,'Skipping %d observations \n',nSkip);
%   tSkip = [0:(nSkip-1)*TR];
%   tTPExclude = [tTPExclude tSkip];
%   %[Slice Par] = clipskip(Slice,Par,nSkip,TR);
%   %nTP = nTP - nSkip;
%end

fprintf(1,'--- Preprocessing Runs ---- \n');
tic;
[fSlice fMean fSlope DOFAdj] = ...
   fMRIPreProc(Slice,HanFilter,0,0);
%   fMRIPreProc(Slice,HanFilter,RemoveBaseline,RemoveTrend);
fprintf(1,'Preprocessing Time: %g\n',toc);

%%%--- PreStim --- %%%%
nPreStim = floor(tPreStim/TR);
fprintf('nPreStim = %d\n',nPreStim);

%%%%%-------------------------------%%%%%%%%%%
fprintf(1,' --- Selectively Avg Runs ------\n');
tic;
[sHAvg sEVar sumXtX DOF Base Trend ] = ...
   SelAvgX(fSlice,Par,TR,nHEst,nPreStim,...
           RemoveBaseline,RemoveTrend,TPExclude);
fprintf(1,'Selective Averaging Time: %g\n',toc);
%%%%%-------------------------------%%%%%%%%%%

if(UsePercent)
  fprintf(1,'Computing Percent Signal Change\n');
  if(nRuns == 1) b = Base;
  else           b = mean(Base,3);
  end
  ind0 = find(b==0);
  l0 = length(ind0);
  if(l0 ~= 0)
    indnz = find(b ~= 0);
    fprintf(1,'Found %d voxels with mean zero\n',length(ind0));
    bmin = min(abs(b(indnz)));
    b(ind0) = bmin;
    fprintf(1,'Resetting zero mean voxels with %g\n',bmin);
  end
  nsH = size(sHAvg,3);
  sHAvg = sHAvg ./ repmat(b, [1 1 nsH]);
  clear b;
end

fprintf(1,'  --- Saving Slice SelXAvg ------\n');
fprintf(1,'     %s\n',hAvgFile);
[ySA dofSA] = sxa2sa(sEVar,sumXtX,sHAvg,nHEst+nPreStim,nTP);
SvBFile(ySA, hAvgFile); % save in selavg format %

%SvBFile(avgEEt,'/homes/nmrnew/home/greve/greve4/rot/eet.bfloat');

% This is now saved in hAvgFile %
%fprintf(1,'  --- Saving  Error Variance ------\n');
%SvBFile(sEVar, eVarFile);

% This is saved in the datFile %
%fprintf(1,'  --- Saving sumXtX File ------\n');
%SvBFile(sumXtX,sumXtXFile);

fprintf(1,'  Saving dat File to \n');
fprintf(1,'     %s\n',datFile);
nbins = nCond(1);
svsxadat(datFile,TR,TimeWindow,tPreStim,nbins,DOF,nRuns,nTP,...
         nRuns,nCols,nSkip,RescaleTarget,sumXtX);

fprintf(1,'  --- Saving dof File ------\n');
fprintf(1,'     %s\n',dofFile);
fid=fopen(deblank(dofFile),'w');
if( fid == -1 )
  msg = sprintf('Could not open dof file %s\n',dofFile);
  qoe(msg);
  error(msg);
end
fprintf(fid,'%d %d %d\n',[ [0:nbins-1]; dofSA; dofSA-1]);
fclose(fid);

fprintf(1,'  --- Saving Voxel Stats ------\n');
fprintf(1,'     %s\n',VxStatFile);
vxstats = [reshape(Base, [nRows*nCols nRuns])'; ...
           reshape(Trend,[nRows*nCols nRuns])'; ...
           reshape(std(fSlice,0,3),[nRows*nCols nRuns])';];

vxstats = reshape(vxstats', [nRows nCols (3*nRuns)]);
SvBFile(vxstats,VxStatFile);

fprintf(1,'-------- SAXSlice: Done ----------\n');
