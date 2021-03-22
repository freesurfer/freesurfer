% isxgslice
% Intersubject Statistical Grinder
%
% global hAvgFile;
% global StatFile;
% global ActiveCond, ControlCond;
% global TestType, HDelMin, HDelMax;
% global OutputFormat 
%  0 = log(p) (natural log)
%  1 = log10(p)
%  2 = p
%  3 = test value
%
% global CmpIdeal
% global QuitOnError;
%
%
%


%
% fmri_isxgslice.m
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

fprintf(1,' $Id: fmri_isxgslice.m,v 1.3 2011/03/02 00:04:06 nicks Exp $ \n');

if( ~exist('QuitOnError') ) QuitOnError = 1; end

%%% ---- Check that all the variables are defined --- %%%
VarNameList = [];
VarNameList = strvcat(VarNameList,'hAvgFile');
VarNameList = strvcat(VarNameList,'datFile');
VarNameList = strvcat(VarNameList,'StatFile');
VarNameList = strvcat(VarNameList,'ActiveCond');
VarNameList = strvcat(VarNameList,'ControlCond');
VarNameList = strvcat(VarNameList,'TestType');
VarNameList = strvcat(VarNameList,'UsePMin');
VarNameList = strvcat(VarNameList,'OutputFormat');
VarNameList = strvcat(VarNameList,'CorrIdeal');
VarNameList = strvcat(VarNameList,'CmpIdeal');

nVar = size(VarNameList,1);
for n = 1:nVar,
  if( exist(deblank(VarNameList(n,:))) ~= 1)
    msg = sprintf('Variable %s does not exist\n',VarNameList(n,:));
    qoe(msg);error(msg);
  end
end
fprintf('UsePMin = %d\n',UsePMin );

if(CmpIdeal)
  VarNameList = [];
  VarNameList = strvcat(VarNameList,'ihMag');
  VarNameList = strvcat(VarNameList,'ihDelta');
  VarNameList = strvcat(VarNameList,'ihTau');
  nVar = size(VarNameList,1);
  for n = 1:nVar,
    if( exist(deblank(VarNameList(n,:))) ~= 1)
      msg = sprintf('Variable %s does not exist\n',VarNameList(n,:));
      qoe(msg);error(msg);
    end
  end
end

if(CorrIdeal)
  VarNameList = [];
  VarNameList = strvcat(VarNameList,'ihDelta');
  VarNameList = strvcat(VarNameList,'ihTau');
  nVar = size(VarNameList,1);
  for n = 1:nVar,
    if( exist(deblank(VarNameList(n,:))) ~= 1)
      msg = sprintf('Variable %s does not exist\n',VarNameList(n,:));
      qoe(msg);error(msg);
    end
  end
end

%%%% ----- Check that variables are the proper values ---- %%%%%
% Check the Test Type %
if( isempty( strmatch(upper(TestType),{'T ','TM','FM','F0','FD','FC','FCD','FDC'},'exact')))
  fprintf(2,'Error: Unkown TestType %s',TestType);
  if(QuitOnError) quit;
  else            return;
  end
end

if( OutputFormat < 0 | OutputFormat > 3)
  fprintf(2,'Error: OutputFormat=%d, must be within 0 to 3\n',OutputFormat);
  if(QuitOnError) quit;
  else            return;
  end
end

%%% -------- read in the dat file ---------- %%
[nNNCond nHEst DOF  TR nRuns nTP nRows nCols ...
   nSkip DTOrder Rescale TW tPreStim HanRad BASeg ...
   GammaFit gfDelta gfTau NullCondId SumXtX] =   fmri_lddat(datFile);
nVoxels = nRows*nCols;

if( exist('HDelMin') & GammaFit > 0)
  msg = 'Cannot specify a delay range with gamma-fit average files';
  qoe(msg); error(msg);
end

if(CorrIdeal & GammaFit > 0)
  msg = 'Cannot correlate ideal HDR with gamma-fit average files';
  qoe(msg); error(msg);
end

fprintf(1,'nNNCond = %d\n',nNNCond);
fprintf(1,'nHEst   = %d\n',nHEst);
fprintf(1,'DOF = %d\n',DOF);
fprintf(1,'nRuns = %d\n',nRuns);
fprintf(1,'nTP = %d\n',nTP);
fprintf(1,'nSkip = %d\n',nSkip);
fprintf(1,'DTOrder = %d\n',DTOrder);
fprintf(1,'TR = %g\n',TR);
fprintf(1,'TW = %g\n',TW);
fprintf(1,'tPreStim= %g\n',tPreStim);
fprintf(1,'NullCondId = %d\n',NullCondId);

fprintf(1,'TestType = %s\n',TestType);

nCond = nNNCond + 1;
Nch = nHEst*nNNCond;
nPreStim = floor(tPreStim/TR);

%% --- Read the hAvg File ---- %%%
fprintf(1,'Reading hAvg File \n');
ysa = fmri_ldbfile(hAvgFile);
ysa = permute(ysa, [3 1 2]);
ysa = reshape(ysa, [nHEst 2 nCond nRows nCols]);
hAvg = squeeze(ysa(:,1,[2:nCond],:,:));
hAvg = reshape(hAvg, [Nch nRows nCols]);
hAvg = permute(hAvg, [2 3 1]);
hStd = squeeze(ysa(:,2,[2:nCond],:,:));
hStd = reshape(hStd, [Nch nRows nCols]);
hStd = permute(hStd, [2 3 1]);
clear ySA;


if(GammaFit > 0) HDelTest = 1;
else
%%%% ---- Set Defaults ---- %%%%
if( ~exist('HDelMin') )
  HDelMin = 0;
  fprintf(2,'Info: Setting HDelMin to 0\n');
end
if( ~exist('HDelMax') )
  HDelMax = (nHEst - 1) * TR - tPreStim;
  fprintf(2,'Info: Setting HDelMax to %g\n',HDelMax);
end
if(HDelMin > HDelMax)
  fprintf(2,'Error: HDelMin=%g > HDelMax=%g\n',HDelMin,HDelMax);
   if(QuitOnError) quit;
   else            return;
   end
end
if(HDelMin < 0)
  fprintf(2,'Error: HDelMin=%g must be > 0\n',HDelMin);
   if(QuitOnError) quit;
   else            return;
   end
end
if(HDelMax > TR*(nHEst-1) )
  fprintf(2,'Error: HDelMax=%g must be < TimeWindow %g\n',...
          HDelMax,TR*(nHEst-1));
   if(QuitOnError) quit;
   else            return;
   end
 end
 newHDelMin = TR*round(HDelMin/TR)+tPreStim;
 newHDelMax = TR*round(HDelMax/TR)+tPreStim;
 fprintf(1,'Info: HDelMin = %g, HDelMax = %g\n',newHDelMin,newHDelMax);
 HDelTest = ([newHDelMin/TR:newHDelMax/TR]') + 1; %'
end %% if(GammaFit)else %%

fprintf(1,'-------------------------\n');
fprintf(1,'Original Conditions Ids \n');
fprintf(1,'Active Condition Ids: ');
fprintf(1,'%d ',ActiveCond);
fprintf(1,'\n');
fprintf(1,'Control Condition Ids: ');
fprintf(1,'%d ',ControlCond);
fprintf(1,'\n');

% Adjust Condition Ids to Account for Null Condition %
ind = find(ActiveCond <= NullCondId);
ActiveCond(ind) = ActiveCond(ind) - NullCondId;
ind = find(ControlCond <= NullCondId);
ControlCond(ind) = ControlCond(ind) - NullCondId;
fprintf(1,'-------------------------\n');
fprintf(1,'Conditions Ids Ajusted for NullCondId = %d\n',NullCondId);
fprintf(1,'Active Condition Ids: ');
fprintf(1,'%d ',ActiveCond);
fprintf(1,'\n');
fprintf(1,'Control Condition Ids: ');
fprintf(1,'%d ',ControlCond);
fprintf(1,'\n');
fprintf(1,'-------------------------\n');

fprintf(1,'Creating Restriction Matrix ------------\n');
RM = fmri_mrestriction(TestType, nHEst, nCond, ...
                       ActiveCond, ControlCond, HDelTest);
if(CorrIdeal)
  fprintf(1,'   Adding Correlation Constraint\n');
  fprintf(1,'   Delta = %g, Tau = %g\n',ihDelta,ihTau);
  RM = fmri_mcorrrestriction(RM,TR,nHEst,ihDelta,ihTau);
end

if(CmpIdeal)
  fprintf(1,'Comparing to Ideal HDIR');
  t = TR*[0:nHEst-1]';%'
  hHDIR = ihMag*fmri_hemodyn(t,ihDelta,ihTau);
  q = repmat(hHDIR,nCond-1,1);
  %%% q = [q; 0]; %%% Mean no longer there
else
  q = zeros(size(hAvg,3),1);
end

fprintf(1,'--- Intersubject Grinding --- \n');
[valStat sigStat polStat] = fmri_isxgrinder(hAvg,hStd,DOF,RM);

if(OutputFormat == 3)
  fprintf(1,' Saving Stat Map to \n');
  fprintf(1,'    %s\n',StatFile);   
  fmri_svbfile(valSig, StatFile);
  fprintf(1,' ---------- fmri_isxgslice.m : Done ------------\n');
  return;
end

if(CmpIdeal)
  sigStat = 1 - sigStat ;
end

if(UsePMin)
  fprintf(1,'Searching for min p-value at each voxel\n');
  [sigStat polStat] = fmri_pmin(sigStat,polStat);
end

if(OutputFormat == 0)     sigStat = -log(sigStat);
elseif(OutputFormat == 1) sigStat = -log10(sigStat);
end

%% Assign Polarity to value %%
sigStat = sigStat .* polStat;

fprintf(1,'Saving Stat Map to \n');
fprintf(1,'    %s\n',StatFile);   
fmri_svbfile(sigStat, StatFile);
fprintf(1,' ---------- fmri_isxgslice.m : Done ------------\n');
