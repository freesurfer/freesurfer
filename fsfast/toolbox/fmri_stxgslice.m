% StXGSlice.m
% StX Slice Grinder
%
% global hAvgFile, eVarFile, sumXtXFile, dofFile;
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
% fmri_stxgslice.m
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

fprintf(1,' ---------- StXGSlice.m : Starting ------------\n');

if( ~exist('QuitOnError') ) QuitOnError = 1; end

%%% ---- Check that all the variables are defined --- %%%
VarNameList = [];
VarNameList = strvcat(VarNameList,'hAvgFile');
%VarNameList = strvcat(VarNameList,'hcovFile');
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
    fprintf(2,'Error: Variable %s does not exist\n',VarNameList(n,:));
    if(QuitOnError) quit; 
    else return
    end
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
      fprintf(2,'Error: Variable %s does not exist\n',VarNameList(n,:));
      if(QuitOnError) quit; 
      else return
      end
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
      fprintf(2,'Error: Variable %s does not exist\n',VarNameList(n,:));
      if(QuitOnError) quit; 
      else return
      end
    end
  end
end

if(exist('cesstem') ~= 1) cesstem = ''; end


%%%% ----- Check that all input files exist ------ %%%%%%%%
InFileName = [];
InFileName = strvcat(InFileName,hAvgFile);
%InFileName = strvcat(InFileName,eVarFile);
%InFileName = strvcat(InFileName,sumXtXFile);
%InFileName = strvcat(InFileName,dofFile);
InFileName = strvcat(InFileName,datFile);
nFile = size(InFileName,1);
for n = 1:nFile,
  if( isempty( dir( deblank(InFileName(n,:)) ) ) )
    fprintf(2,'Error: Cannot find file %s\n',InFileName(n,:) );
    if(QuitOnError) quit; 
    else return
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
%[nNNCond nHEst DOF  TR nRuns nTP nRows nCols ...
%   nSkip DTOrder Rescale TW tPreStim HanRad BASeg ...
%   GammaFit gfDelta gfTau NullCondId SumXtX] =   fmri_lddat(datFile);
%hd = fmri_lddat2(datFile);
%ySA = fmri_ldbfile(hAvgFile);
%[hAvg eVar] = fmri_sa2sxa(ySA,hd.Nh);
%clear ySA;

%% --- Read the hAvg File ---- %%%
fprintf(1,'Reading hAvg File \n');
[hAvg eVar hd] = fast_ldsxabfile(hAvgFile);

if( exist('HDelMin') & hd.GammaFit > 0)
  msg = 'Cannot specify a delay range with gamma-fit average files';
  qoe(msg); error(msg);
end

if(CorrIdeal & hd.GammaFit > 0)
  msg = 'Cannot correlate ideal HDR with gamma-fit average files';
  qoe(msg); error(msg);
end

nVoxels = hd.Nrows*hd.Ncols;

fprintf(1,'nNNCond = %d\n',hd.Nnnc);
fprintf(1,'nHEst   = %d\n',hd.Nh);
fprintf(1,'DOF = %d\n',hd.DOF);
fprintf(1,'nRuns = %d\n',hd.Nruns);
fprintf(1,'nTP = %d\n',hd.Ntp);
fprintf(1,'nSkip = %d\n',hd.Nskip);
fprintf(1,'DTOrder = %d\n',hd.DTOrder);
fprintf(1,'TR = %g\n',hd.TR);
fprintf(1,'TER = %g\n',hd.TER);
fprintf(1,'TW = %g\n',hd.TimeWindow);
fprintf(1,'tPreStim= %g\n',hd.TPreStim);
fprintf(1,'NullCondId = %d\n',hd.NullCondId);

fprintf(1,'TestType = %s\n',TestType);

nPreStim = floor(hd.TPreStim/hd.TER);

%Ch = fmri_ldbfile(hcovFile);
Ch = hd.hCovMtx;

if(~isempty(CMtxFile))
  fprintf('Loading %s\n',CMtxFile);
  tmp = load(CMtxFile);
  RM = tmp.ContrastMtx_0;
  nRM = size(RM,2);
  nh  = size(hAvg,3);
  if(nRM ~= nh)
    msg = sprintf('hAvg size (%d) is inconsistent with CMtx (%d)',nh,nRM);
    qoe(msg);error(msg);
  end
else
  if(hd.GammaFit > 0) HDelTest = 1;
  else
  %%%% ---- Set Defaults ---- %%%%
  if( ~exist('HDelMin') )
    HDelMin = -hd.TPreStim;
    fprintf(2,'Info: Setting HDelMin to %g\n',HDelMin);
  end
  if( ~exist('HDelMax') )
    HDelMax = (hd.Nh - 1) * hd.TER - hd.TPreStim;
    fprintf(2,'Info: Setting HDelMax to %g\n',HDelMax);
  end
  if(HDelMin > HDelMax)
    fprintf(2,'Error: HDelMin=%g > HDelMax=%g\n',HDelMin,HDelMax);
     if(QuitOnError) quit;
     else            return;
     end
  end
  if(HDelMin < -hd.TPreStim)
    fprintf(2,'Error: HDelMin=%g must be > %g\n',HDelMin,-hd.TPreStim);
     if(QuitOnError) quit;
     else            return;
     end
  end
  if(HDelMax > hd.TER*(hd.Nh-1) )
    fprintf(2,'Error: HDelMax=%g must be < TimeWindow %g\n',...
            HDelMax,hd.TER*(hd.Nh-1));
     if(QuitOnError) quit;
     else            return;
     end
   end
   newHDelMin = hd.TER*round(HDelMin/hd.TER)+hd.TER*round(hd.TPreStim/hd.TER);
   newHDelMax = hd.TER*round(HDelMax/hd.TER)+hd.TER*round(hd.TPreStim/hd.TER);
   fprintf(1,'Info: newHDelMin = %g, newHDelMax = %g\n',newHDelMin,newHDelMax);
   HDelTest = round([newHDelMin/hd.TER:newHDelMax/hd.TER]') + 1; %'
   fprintf(1,'Info: HDelTest '); 
   fprintf(1,' %d',HDelTest);
   fprintf(1,'\n');
  end %% if(hd.GammaFit)else %%

  fprintf(1,'-------------------------\n');
  fprintf(1,'Original Conditions Ids \n');
  fprintf(1,'Active Condition Ids: ');
  fprintf(1,'%d ',ActiveCond);
  fprintf(1,'\n');
  fprintf(1,'Control Condition Ids: ');
  fprintf(1,'%d ',ControlCond);
  fprintf(1,'\n');

  % Adjust Condition Ids to Account for Null Condition %
  if(hd.NullCondId ~= 0)
    ind0    = find(ActiveCond == 0);
    indNull = find(ActiveCond == hd.NullCondId);
    ActiveCond(ind0)    = hd.NullCondId;
    ActiveCond(indNull) = 0;
  
    ind0    = find(ControlCond == 0);
    indNull = find(ControlCond == hd.NullCondId);
    ControlCond(ind0)    = hd.NullCondId;
    ControlCond(indNull) = 0;
  end

  fprintf(1,'-------------------------\n');
  fprintf(1,'Conditions Ids Ajusted for NullCondId = %d\n',hd.NullCondId);
  fprintf(1,'Active Condition Ids: ');
  fprintf(1,'%d ',ActiveCond);
  fprintf(1,'\n');
  fprintf(1,'Control Condition Ids: ');
  fprintf(1,'%d ',ControlCond);
  fprintf(1,'\n');
  fprintf(1,'-------------------------\n');

  fprintf(1,'Creating Restriction Matrix ------------\n');
  RM = fmri_mrestriction(TestType, hd.Nh, hd.Nc, ...
                       ActiveCond, ControlCond, HDelTest);

  if(CorrIdeal)
    fprintf(1,'   Adding Correlation Constraint\n');
    fprintf(1,'   Delta = %g, Tau = %g\n',ihDelta,ihTau);
    RM = fmri_mcorrrestriction(RM,hd.TER,hd.Nh,ihDelta,ihTau);
  end
end

q = zeros(size(hAvg,3),1);

fprintf(1,'Grinding  ------------\n');
[vSig pSig ces] = fmri_stxgrinder(TestType,hAvg,eVar,Ch,hd.DOF,RM,q);

if(~isempty(statstem))
  fname = sprintf('%s_%03d.bfloat',statstem,sliceno);
  fmri_svbfile(vSig, fname); 
end

if(~isempty(cesstem))
  % Compute ces as percent of baseline 
  instem = fmri_getstem(hAvgFile);
  hoffsetname = sprintf('%s-offset_%03d.bfloat',instem,sliceno);
  hoffset = fmri_ldbfile(hoffsetname);
  if(isempty(hoffset))
    fprintf('ERROR: could not load %s\n',hoffsetname);
    return;
  end
  indz = find(hoffset==0);
  hoffset(indz) = 10^10;
  nces = size(ces,3);
  cespct = 100*ces./repmat(hoffset,[1 1 nces]);

  fname = sprintf('%s_%03d.bfloat',cesstem,sliceno);
  fmri_svbfile(ces, fname); 
  fname = sprintf('%spct_%03d.bfloat',cesstem,sliceno);
  fmri_svbfile(cespct, fname); 
end

if(CmpIdeal)
  pSig = 1 - pSig ;
end

if(strncmp('T',upper(TestType),1))
  SignMask = ((vSig>=0) - (vSig<0));
  pSig = pSig .* SignMask;
else
  SignMask = [];
end

if( ~isempty(pminstem) | ~isempty(ipminstem) )
  [pSigMin ipSigMin] = min(abs(pSig),[],3);

  %% Go through some gymnastics to get the signed pmin %%
  [nrows ncols nplanes] = size(pSig);
  xx1 = [1:nrows]' * ones(1,ncols); %'
  I1 = reshape1d(xx1);
  clear xx1;
  xx2 = ones(nrows,1) * [1:ncols]; 
  I2 = reshape1d(xx2);
  clear xx2;
  itmp = sub2ind(size(pSig),I1,I2,reshape1d(ipSigMin));
  pSigMin = reshape(pSig(itmp), [nrows ncols]);      
  clear I1 I2 itmp;

  % Bonferoni Correction
  pSigMin = pSigMin * size(pSig,3);  

  if(~isempty(pminstem))
    if(OutputFormat == 0)     pSigMin = -log(abs(pSigMin))  .*sign(pSigMin);
    elseif(OutputFormat == 1) pSigMin = -log10(abs(pSigMin)).*sign(pSigMin);
    end
    fname = sprintf('%s_%03d.bfloat',pminstem,sliceno);
    fmri_svbfile(pSigMin, fname); 
  end

  if(~isempty(ipminstem))
    ind = find(abs(pSigMin) < 4);
    ipSigMin(ind) = 0;
    fname = sprintf('%s_%03d.bfloat',ipminstem,sliceno);
    fmri_svbfile(ipSigMin, fname);  
  end

  SignMask = ((pSigMin>=0) - (pSigMin<0));

end

% Dont let pSig = zero (messes up log10)
iz = find(abs(pSig) < 10^-300);
pSig(iz) = sign(pSig(iz)) * 10^-300;
iz = find(pSig == 0);
pSig(iz) = 10^-300; % Have to do this because sign(0) = 0

if(OutputFormat == 0)     pSig = -log(abs(pSig))  .*sign(pSig);
elseif(OutputFormat == 1) pSig = -log10(abs(pSig)).*sign(pSig);
end
fmri_svbfile(pSig, StatFile);

if(~isempty(fstatstem) | ~isempty(fsigstem) )
  fprintf(1,'Performing F-Test \n');
  [vSig pSig] = fmri_stxgrinder('F',hAvg,eVar,Ch,hd.DOF,RM,q);

  if(size(RM,1) > 1)
    tmp = RM * reshape(hAvg,[size(hAvg,1)*size(hAvg,2) size(hAvg,3)])'; %'
    tmp = mean(tmp);
    SignMask = sign(tmp);
    SignMask = reshape(SignMask',[size(hAvg,1) size(hAvg,2)]); %'
    clear hAvg tmp;
    pSig = pSig.*SignMask;
  end

  if(OutputFormat == 0)     pSig = -log(abs(pSig))  .*sign(pSig);
  elseif(OutputFormat == 1) pSig = fast_log10p(pSig);
  end

  if(~isempty(fstatstem))
    fname = sprintf('%s_%03d.bfloat',fstatstem,sliceno);
    fmri_svbfile(vSig, fname);
  end

  if(~isempty(fsigstem))
    fname = sprintf('%s_%03d.bfloat',fsigstem,sliceno);
    fmri_svbfile(pSig, fname);
  end

end

