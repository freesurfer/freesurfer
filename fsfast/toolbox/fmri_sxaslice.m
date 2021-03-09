% fmri_sxaslice2.m
% March 31, 1999
% Modification of fmri_sxaslice to include fitting for the amplitude
% of a specified hemodyn response.
% Douglas N. Greve, MGH-NMR 
%
% Performs selective averaging using a matrix formulation
% (this is the same as deconvolution) assuming white noise.
%
% This is for multiple runs in a single session.
%
%
%
%


%
% fmri_sxaslice.m
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
fprintf(1,'  --- fmri_sxaslice: Starting ------\n');
fprintf(1,'fmri_sxaslice.m @FS_VERSION@\n');

if( ~exist('QuitOnError') ) QuitOnError = 1; end

%%% ---- Check that all the variables are defined --- %%%
VarNameList = [];
VarNameList = strvcat(VarNameList,'TimeWindow');
VarNameList = strvcat(VarNameList,'TimeOffset');
VarNameList = strvcat(VarNameList,'RemoveBaseline');
VarNameList = strvcat(VarNameList,'RemoveTrend');
VarNameList = strvcat(VarNameList,'TR');
VarNameList = strvcat(VarNameList,'TER');
VarNameList = strvcat(VarNameList,'nSkip');
VarNameList = strvcat(VarNameList,'NullCondId');
VarNameList = strvcat(VarNameList,'Nconds');
VarNameList = strvcat(VarNameList,'GammaFit');
VarNameList = strvcat(VarNameList,'gfDelta');
VarNameList = strvcat(VarNameList,'gfTau');
VarNameList = strvcat(VarNameList,'RescaleTarget');
VarNameList = strvcat(VarNameList,'tPreStim');
VarNameList = strvcat(VarNameList,'InputFiles');
VarNameList = strvcat(VarNameList,'ParFiles');
VarNameList = strvcat(VarNameList,'hAvgFile');
VarNameList = strvcat(VarNameList,'dofFile');
%VarNameList = strvcat(VarNameList,'hcovFile');
VarNameList = strvcat(VarNameList,'segFile');
VarNameList = strvcat(VarNameList,'ReeFile');
VarNameList = strvcat(VarNameList,'datFile');
VarNameList = strvcat(VarNameList,'VxStatFile');
nVar = size(VarNameList,1);
for n = 1:nVar,
  if( exist(deblank(VarNameList(n,:))) ~= 1)
    msg = sprintf('Error: Variable %s does not exist',VarNameList(n,:));
    qoe(msg);error(msg);
  end
end

if(exist('NoiseType') ~= 1) NoiseType = 0; end
if(NoiseType ~= 0) 
  if(exist('BASegment') ~= 1) 
    msg = 'Variable BASegment does not exist';
    qoe(msg);error(msg);
  end
  if(exist('nNACorr') ~= 1) 
    msg = 'Variable nNACorr does not exist';
    qoe(msg);error(msg);
  end
  if(nNACorr < 2) NoiseType = 0; end
end

if(NoiseType == 0) fprintf(1,'Assuming white noise\n'); end

if(exist('SegStem') == 0) SegStem = ''; end

%%% Get the number of runs from the number of input files %%%
nRuns = size(InputFiles,1);
if( size(ParFiles,1) ~= nRuns )
  msg = 'Incorrect number of ParFiles';
  qoe(msg);error(msg);
end

if( size(hAvgFile,1) ~= 1 )
  msg = 'Incorrect number of hAvgFile';
  qoe(msg);error(msg);
end

%% Compute the number of delay-points in the HDR to estimate
nHEst = floor(TimeWindow/TER);
if(nHEst < 1)
  msg = 'TimeWindow too small, but be > TER';
  qoe(msg);error(msg);
end

%%%%---Compute Detrending Order --------- %%%%%%%
DTOrder = 0;
if(RemoveBaseline) DTOrder = 1; end
if(RemoveTrend)    DTOrder = 2; end

%%% --- Check that Percent Change is accompanied by DTOrder > 0 %%%
if(UsePercent & DTOrder == 0)
  msg = 'Percent Signal Change must be used with baseline or detrend options';
  qoe(msg);error(msg);  
end

%% -- Load the parfiles --- %%
Par = fmri_ldpar(ParFiles);
CondIdMap = par2condidmap(Par,NullCondId);
tmp = length(CondIdMap);
if(Nconds < 0)
  Nconds = tmp;
else
  if(Nconds ~= tmp)
    msg = sprintf('Nconds = %d but specified as %d\n',tmp,Nconds);
    qoe(msg);error(msg);
  end
end
Par = remap_par(Par,CondIdMap);

%% Check the number of conditions for each run %%
nCondPerRun = reshape1d(max(Par(:,2,:)));
fprintf(1,'Non-Null Conditions per Run: ');
fprintf(1,'%d ',nCondPerRun);
fprintf(1,'\n');

%%% ---- Load the fMRI Runs ------- %%%
Slice = fmri_ldbfile(InputFiles);
if(Synth)
  if(Seed < 0)
    Seed = sum(100*clock);
  end
  fprintf('Synthesizing input data with seed %d\n',Seed);
  randn('state',Seed);
  Slice = randn(size(Slice));
end
nRows = size(Slice,1);
nCols = size(Slice,2);
nTP   = size(Slice,3);
nV    = nRows*nCols;

%%%%-------- Rescale --------- %%%%%%%
if(RescaleTarget > 0)
  fprintf(1,'Rescaling to Mean = %g\n',RescaleTarget);
  for run = 1:nRuns
    runfile = deblank(InputFiles(run,:));
    stem = fmri_parsebfilename(runfile);
    meanvalfile = sprintf('%s.meanval',stem);
    l = dir(meanvalfile);
    if(size(l,1) == 0)
      msg = sprintf('Cannot open %s',meanvalfile);
      qoe(msg); error(msg);
    end
    meanval = load(meanvalfile);
    rescalefactor = RescaleTarget/meanval;
    fprintf('Run %2d: meanval = %g, rsfactor = %g\n',run,meanval,rescalefactor);
    Slice(:,:,:,run) = rescalefactor*Slice(:,:,:,run) ;
  end
end


% Compute the sub-sampling matrix
[Mss Rss]   = fmri_subsampmat(TR,nTP,TER);
Nts = nTP*Rss;

fprintf(1,'nRows = %d, nCols = %d, nTP = %d, nRuns = %d\n',...
        nRows,nCols,nTP,nRuns);
Nv = nRows*nCols;

%%%% --- Check for exclusions (including skips) ---%%%%
TPExclude   = zeros(nTP,nRuns);
ssTPExclude = zeros(Nts,nRuns);
if(nSkip > 0) tSkip = TR*[0:(nSkip-1)]';  %'
else          tSkip = [];
end

ntpxtot = 0;
for n = 1:nRuns,
  if(exist('TPExclFiles'))
    tpexcl = fmri_ldtpexcl(TPExclFiles(n,:));
    tpexcl = [tpexcl; tSkip;];
    tpexcl = unique(tpexcl);
  else
    tpexcl = tSkip;
  end

  fprintf('Run %2d Exclusions: ',n);
  fprintf('%2d ',tpexcl'); %'
  fprintf('\n');
  
  l = length(tpexcl);
  if(l ~= 0)
    TPExclude(floor(tpexcl/TR)+1,n)   = ones(l,1);
    ssTPExclude(floor(tpexcl/TER)+1,n) = ones(l,1);
  end
  ntpxtot = ntpxtot + l;
end

%%% -- Spatial Filtering ---- %%
if( ~exist('HanRadius') ) 
  HanFilter = [];
  HanRadius = 0;
  fSlice = Slice;
else
  if(HanRadius == 0)
    HanFilter = [];
    fSlice = Slice;
  elseif(HanRadius < 1)
    msg = sprintf('Error: HanRadius = %g, must be >= 1',HanRadius);
    qoe(msg);error(msg);
  else
    fprintf(1,'Using Spatial Filter, HanRad = %g\n',HanRadius);
    HanFilter = fmri_hankernel(HanRadius);
    disp(HanFilter);
    fSlice = fmri_spatfilter(Slice,HanFilter);
  end
end
%clear Slice;

%%%--- PreStim --- %%%%
nPreStim = floor(tPreStim/TER);
fprintf('nPreStim = %d\n',nPreStim);

%%%------ Slice Acqusition Delay Synch ---------%%%%%
if(~isempty(SynchFile))
  synch = fmri_ldsynch(SynchFile);
  ind = find(synch(:,1)==SliceNo);
  if(isempty(ind))
    msg = sprintf('Slice %d not found in SynchFile %s',SliceNo,SynchFile);
    qoe(msg);error(msg);
  end
  TSliceSynch = synch(ind,2);
else
  TSliceSynch = 0;
end

if(interleave)
  dt = TR/nslicesvol;
  nn = [ [0:2:nslicesvol-1]; [1:2:nslicesvol-1]];
  SliceSeq = reshape1d(nn'); %'
  fprintf(1,'SliceSeq: ');
  fprintf(1,'%2d ',SliceSeq);
  fprintf(1,'\n');

  ind = find(SliceSeq == SliceNo)-1;
  if(isempty(ind))
    msg = sprintf('Could not find SliceNo %d in SliceSeq',SliceNo);
    qoe(msg); error(msg);
  end
  TSliceSynch = dt*ind;
  fprintf('SliceNo=%d, TR=%g, nslicesvol=%d, dt=%g, ind=%d\n',...
        SliceNo,TR,nslicesvol,dt,ind);
end

fprintf('TSliceSynch = %g \n', TSliceSynch);

TimeOffset = TimeOffset - TSliceSynch;

%%%%-------- Generate Stim Conv Matricies --------- %%%%%%%
fprintf(1,'Creating Stim Conv Matricies \n');
if(NullCondId ~= 0)
  indNull = find(Par(:,2,:)== NullCondId);
  ind0    = find(Par(:,2,:)== 0);
  Par(:,ind0,:)    = NullCondId;
  Par(:,indNull,:) = 0;
end

Par(:,1,:) = Par(:,1,:) + TimeOffset;
[X nNNCond nPerCond] =  ...
   fmri_par2scm(Par,Nconds,Nts,TER,nHEst,tPreStim,ssTPExclude);
if(Rss ~= 1)
  X2 = zeros(nTP,nHEst*(Nconds-1),nRuns);
  for n = 1:nRuns,
    X2(:,:,n) = Mss * X(:,:,n);
  end
  X = X2;
  clear X2;
end

Xfir = X;
if(GammaFit == 1 ) 
  fprintf(1,'Using Gamma Fit: Delta = %g, Tau = %g\n',gfDelta,gfTau);
  X = fmri_scm2gcm(X,nNNCond,TR,tPreStim,gfDelta,gfTau);
  nHEst = 1;
end

Xorig = X;

y = fSlice;

Xdetrend1 = [];
Xdetrend2 = [];
if(DTOrder > 0)
  Xdetrend1 = [ones(nTP,1) zeros(nTP,nRuns-1)];
end
if(DTOrder > 1)
  Xdetrend2 = [[0:nTP-1]' zeros(nTP,nRuns-1)]; %'
end

%% Modify SCM to remove the global means %%%
Xdetrend = [];
X2 = X; 
clear X;
sumXtX = 0;
xx = [];
for n = 1:nRuns,
  x = X2(:,:,n);
  tmpxtx = x'*x; %'
  c = cond(tmpxtx);
  tr = trace(tmpxtx);
  eff = 1/trace(inv(tmpxtx));
  %fprintf(1,'Run %2d Condition XtX: %g, Eff = %g,\n',n,c,eff);
  sumXtX = sumXtX + X2(:,:,n)' * X2(:,:,n); %'

  if(DTOrder > 0) 
    Xdetrend = [Xdetrend1 Xdetrend2];
    tmp = repmat(~TPExclude(:,n), [1 DTOrder*nRuns]);
    Xdetrend = Xdetrend .* tmp;
    X(:,:,n)  = [X2(:,:,n) Xdetrend]; % add a column of ones
    Xdetrend1 = fmri_shiftcol(Xdetrend1,1,0) ;
    Xdetrend2 = fmri_shiftcol(Xdetrend2,1,0) ;
  end

  xx = [xx; X(:,:,n)];

end
if(DTOrder == 0) X = X2; end

condXtX = cond(sumXtX);
fprintf(1,' SumXtX Condition: %g, Eff = %g,\n',condXtX,1/trace(inv(sumXtX)));
if(condXtX > 10000)
  fprintf(1,' ERROR: paradigm is ill-conditioned.\n');
  fprintf(1,' Check your paradigm files for periodic stimuli\n');
  qoe(''); error('');
end

%fprintf(1,'Condition SumXtX: %g\n',cond(sumXtX));
fprintf(1,'Average Variance Reduction: %g\n',1/mean(diag(inv(sumXtX))));
fprintf(1,'Average StdErr   Reduction: %g\n',1/sqrt(mean(diag(inv(sumXtX)))));
%clear X2 Xdetrend Xdetrend1 Xdetrend2 tmp;

X0 = X;
fprintf(1,'Computing hCoVar\n');
ch          = fmri_hcovar(X);
fprintf(1,'Selective Sum\n');
ss          = fmri_selsum(X,y);
fprintf(1,'LMS Fit\n');
hest        = fmri_lmsfit(ch,ss);
%clear ss;
fprintf(1,'Residuals\n');
[eest sigest] = fmri_residual(y,X,hest);

DOF  = (nTP-DTOrder)*nRuns-nNNCond*nHEst - ntpxtot;
fprintf(1,'DOF = %d\n',DOF);

fprintf(1,'VoxVar\n');
voxvar = fmri_voxvar(eest,DOF,[],TPExclude);


fprintf(1,'Segmenting Brain and Air\n');
mnimg    = mean(y(:,:,:,1),3);
[ibrain iair thresh] = fmri_segment(mnimg);
if(~isempty(SegStem))
  segbrain = zeros(size(mnimg));
  segbrain(ibrain) = 1;
  fname = sprintf('%s_%03d.bshort',SegStem,SliceNo);
  fmri_svbfile(segbrain,fname);
end

%%------- Compute covariance matrix of residual error -------------%%
if(~isempty(ecovmtxstem))
  ecm = fmri_cvmstruct;
  ecm.d = TR;
  if(length(ibrain) == 0)
    fprintf(1,'INFO: no brain voxels found for slice %d',SliceNo);
    ecm.cvm = zeros(nTP*nRuns);
    ecm.n   = 0;
  else
    tic;
    fprintf(1,'INFO: computing residual error covariance matrix\n');

    %% Extract brain voxels %%
    tmp = reshape(eest,[Nv nTP*nRuns])'; %'
    
    if(BASegment) tmp = tmp(:,ibrain); end

    %% Temporal mean has already been removed from the res error %%

    %% Compute the covariance matrix %%
    ecm.n   =  size(tmp,2);
    ecm.cvm =  (tmp * tmp')/ecm.n; % '
    clear tmp;

    %% Save to bfile and .cvm %%
    slcstem = sprintf('%s_%03d',ecovmtxstem,SliceNo);
    fmri_svcvm(ecm,slcstem);
    fprintf(1,'INFO: finished computing error covariance matrix (%f)\n', toc);
  end

  %% Compute ideal ECM (noise spatial and temporally white) %%
  ecm = fmri_cvmstruct;
  ecm.d = TR;
  m = xx*inv(xx'*xx)*xx';
  ecm.cvm = eye(size(m)) - m;
  idealstem = sprintf('%s-ideal',ecovmtxstem);
  fmri_svcvm(ecm,idealstem);

  clear ecm, m;
end

%%- Compute spatially norm covariance matrix of residual error ----%%
if(~isempty(ecvmsnstem))
  ecm = fmri_cvmstruct;
  ecm.d = TR;
  if(length(ibrain) == 0)
    fprintf(1,'INFO: no brain voxels found for slice %d',SliceNo);
    ecm.cvm = zeros(nTP*nRuns);
    ecm.n   = 0;
  else
    tic;
    fprintf(1,'INFO: computing spat norm residual error covariance matrix\n');

    %% Extract brain voxels %%
    tmp = reshape(eest,[Nv nTP*nRuns])'; %'
    
    if(BASegment) tmp = tmp(:,ibrain); end

    %% Temporal mean has already been removed from the res error %%

    %% Compute per-voxel temporal stddev %%
    tmpstd1 = std(tmp);
    ind = find(tmpstd1==0);
    tmpstd1(ind) = 1;
    tmpstd2 = repmat(tmpstd1, [nTP*nRuns 1]);
   
    %% Spatially normalize residual error %%
    tmp = tmp./tmpstd2;
    clear tmpstd1 tmpstd2;

    %% Compute the covariance matrix %%
    ecm.n   =  size(tmp,2);
    ecm.cvm =  (tmp * tmp')/ecm.n; % '
    clear tmp;

    %% Save to bfile and .cvm %%
    slcstem = sprintf('%s_%03d',ecvmsnstem,SliceNo);
    fmri_svcvm(ecm,slcstem);
    fprintf(1,'INFO: finished computing error covariance matrix (%f)\n', toc);
  end

  clear ecm;
end

%%------- Compute covariance matrix of signal -------------%%
if(~isempty(scovmtxstem))
  scm = fmri_cvmstruct;
  scm.d = TR;
  if(length(ibrain) == 0)
    fprintf(1,'INFO: no brain voxels found for slice %d',SliceNo);
    scm.cvm = zeros(nTP*nRuns);
    scm.n   = 0;
  else
    tic;
    fprintf(1,'INFO: computing signal covariance matrix\n');

    %% Extract brain voxels %%
    tmp = reshape(sigest,[Nv nTP*nRuns])'; %'
    
    if(BASegment) tmp = tmp(:,ibrain); end

    %% Remove the temporal mean %%
    tmpmean = mean(tmp);
    tmp = tmp - repmat(tmpmean, [nTP*nRuns 1]);

    %% Compute the covariance matrix %%
    scm.n   =  size(tmp,2);
    scm.cvm =  (tmp * tmp')/scm.n; % '
    clear tmp;

    %% Save to bfile and .cvm %%
    slcstem = sprintf('%s_%03d',scovmtxstem,SliceNo);
    fmri_svcvm(scm,slcstem);
    fprintf(1,'INFO: finished computing signal covariance matrix (%f)\n', toc);
  end

  %% Compute ideal SCM (no signal, noise spatial and temporally white) %%
  iscm = fmri_cvmstruct;
  iscm.d = TR;
  iscm.cvm = xx*inv(xx'*xx)*xx';
  idealstem = sprintf('%s-ideal',scovmtxstem);
  fmri_svcvm(iscm,idealstem);

  clear scm, iscm;
end

if(NoiseType == 1 ) % Compensate for correlations in noise

  fprintf(1,'-------- Compensating for correlated noise -----\n');
  nRee = nNACorr;
  fprintf(1,'nRee = %d\n',nRee);

  fprintf(1,'Computing Error AutoCorrelation Function\n');
  if(BASegment) Ree = fmri_acorr(eest,nRee,repmat(ibrain,nRuns));
  else          Ree = fmri_acorr(eest,nRee);
  end

  avgRee = mean(Ree')';

  fprintf(1,'Fitting Average AutoCorrelation Function\n');
  [alpha0 rho0] = fmri_acorrfit(avgRee,nRee);
  alpha = repmat(alpha0, [1 nRuns]);
  rho   = repmat(rho0,   [1 nRuns]);

  fprintf(1,'run alpha rho\n');
  fprintf(1,'%3d %g    %g\n',[[1:nRuns]' alpha' rho']');

  ind = find(rho>=1);
  if(~isempty(ind))
    fprintf(1,'WARNING: some rhos are > 1 ... setting to .9\n')
    rho(ind) = .9;
  end

  fprintf(1,'Synthesizing AutoCorrelation Function\n');
  synthRee = fmri_acorrsynth(alpha,rho,nRee);
  fprintf(1,'Ree: Measured Synthetic\n');
  ReeErr = Ree-synthRee;
  fprintf(1,'Std of Ree Err %g\n',std(ReeErr))
  %[Ree synthRee]
  %synthRee = Ree;

  tt = TR*[-nRee:nRee]'; %'
  tmp = [tt Ree synthRee];
  fmt = repmat('%g ',[1 size(tmp,2)]);
  fmt = strcat(fmt,'\n');
  fid = fopen(ReeFile,'w');
  fprintf(fid,'# %3d %g    %g\n',[[1:nRuns]' alpha' rho']');
  if(BASegment)   fprintf(fid,'# nBrainVoxels %d\n',length(ibrain)); end
  fprintf(fid,fmt,tmp'); %'
  fclose(fid);

  fprintf(1,'Synthesizing Covariance Matrix\n');
  Ce = fmri_acorr2covmtx(synthRee,nTP);

  fprintf(1,'Recomputing hCoVar\n');
  ch          = fmri_hcovar(X,Ce);

  fprintf(1,'Checking hCoVar for Positive Definiteness\n');
  cheig = eig(ch);
  mineig = min(cheig);
  fprintf(1,'Eigenvalue range for hCovar: %g %g\n',mineig,max(cheig));
  if(mineig<0)
    fprintf(1,'WARNIING: hCoVar is not Positive Definite\n');
    fprintf(1,'   ... attempting to fix ...\n');
    ch = ch - eye(size(ch))*mineig*1.01;
    cheig = eig(ch);
    mineig = min(cheig);
    fprintf(1,'Eigenvalue range for hCovar: %g %g\n',mineig,max(cheig));
    if(mineig<0)
      msg = 'ERROR: hCoVar is not Positive Definite';
      qoe(msg); error(msg);
    end
  end
 
  fprintf(1,'Recomputing Selective Sum\n');
  ss          = fmri_selsum(X,y,Ce);

  fprintf(1,'Recomputing LMS Fit\n');
  hest        = fmri_lmsfit(ch,ss);
  clear ss;

  fprintf(1,'Recomputing Residuals\n');
  [eest sigest] = fmri_residual(y,X,hest);

  fprintf(1,'VoxVar\n');
  voxvar      = fmri_voxvar(eest,DOF,Ce);

  clear Ce;
end

evar   = std(reshape(eest, [nRows nCols nTP*nRuns]) ,[],3).^2;
sigvar = std(reshape(sigest, [nRows nCols nTP*nRuns]) ,[],3).^2;


if(~isempty(EResDir))
  for n = 1:nRuns
     fname = sprintf('%s/eres-r%03d_%03d.bfloat',EResDir,n,SliceNo);
     fmri_svbfile(eest(:,:,:,n),fname);
  end
end

if(~isempty(SignalDir))
  % Recompute sigest without trend and baseline %
  [tmp sigest] = fmri_residual(y,X,hest);
  clear tmp;
  for n = 1:nRuns
     fname = sprintf('%s/sigest-r%03d_%03d.bfloat',SignalDir,n,SliceNo);
     fmri_svbfile(sigest(:,:,:,n),fname);
  end
end

% Remove the information about the offset/trend from h %%
nTotEst = size(hest,3);
if(DTOrder > 0)
  n = nTotEst-DTOrder*nRuns;
  hoffset = hest(:,:,n+1:n+nRuns);
  if(DTOrder > 1)
    n2 = nTotEst-nRuns;
    htrend = hest(:,:,n2+1:n2+nRuns);
  end
  hest = hest(:,:,1:n);
  ch   =   ch(1:n,1:n);
  X = X(:,1:n,:);
  nTotEst = size(hest,3);
else
  hoffset = mean(squeeze(mean(y,3)),3);
end

SliceMean = mean(hoffset,3);
if(nRuns == 1) SliceStd = std(y,[],3);
else           SliceStd = sqrt(mean(squeeze(std(y,[],3)).^2,3));
end
fprintf(1,'Segmenting Threshold %g\n',thresh);
fprintf(1,'Brain Voxels: %d\n',length(ibrain));
fprintf(1,'Air   Voxels: %d\n',length(iair));
SegMask = zeros(nRows,nCols);
if(~isempty(ibrain)) SegMask(ibrain(:,1)) = 1;end
clear tmp;
tmp(:,:,1) = SegMask;
tmp(:,:,2) = SliceMean;
tmp(:,:,3) = SliceStd;
fmri_svbfile(SliceMean,segFile);
%clear SegMask SlcieStd tmp fSlice;

if(~isempty(ibrain)) 
SigPower   = mean(reshape1d(sigvar(ibrain(:,1))));
NoisePower = mean(reshape1d(evar(ibrain(:,1))));
SNR1 = SigPower/NoisePower;
SNR2 = mean(reshape1d(sigvar(ibrain(:,1))./evar(ibrain(:,1))));
%fprintf('SigPower = %g, NoisePower = %g, SNR1 = %g, SNR2 = %g\n',...
%        SigPower, NoisePower, SNR1, SNR2);
end

%clear eest;
%clear y;

%fmri_svbfile(ch,hcovFile);
clear tmp;

fprintf(1,'Saving dat File to \n');
fprintf(1,'   %s\n',datFile);
hd = fmri_hdrdatstruct;
hd.TR = TR;
hd.TER = TER;
hd.TimeWindow = TimeWindow;
hd.TPreStim = tPreStim;
hd.Nc = nNNCond + 1;
hd.Nh = nHEst;
hd.Nnnc = nNNCond;
hd.DOF  = DOF;
hd.Npercond  = nPerCond;
hd.Nruns = nRuns;
hd.Ntp = nTP;
hd.Nrows = nRows;
hd.Ncols = nCols;
hd.Nskip = nSkip;
hd.DTOrder = DTOrder;
hd.RescaleFactor = RescaleTarget;
hd.HanningRadius = HanRadius;
hd.BrainAirSeg   = BASegment;
hd.GammaFit      = GammaFit;
hd.gfDelta       = gfDelta;
hd.gfTau         = gfTau;
hd.NullCondId    = NullCondId;
hd.SumXtX        = sumXtX;
hd.nNoiseAC      = nNACorr;
hd.CondIdMap     = CondIdMap;
hd.hCovMtx       = ch;
fmri_svdat2(datFile,hd);

%fmri_svdat(datFile,TR,TimeWindow,tPreStim,nNNCond,DOF,nRuns,nTP,...
%        nRows,nCols,nSkip,DTOrder,RescaleTarget,HanRadius,nNACorr,BASegment,...
%        GammaFit,gfDelta,gfTau,NullCondId,sumXtX);
% clear sumXtX;

fprintf(1,'--- Saving Slice SelXAvg ------\n');
fprintf(1,'   %s\n',hAvgFile);
[ySA dof] = fmri_sxa2sa(voxvar,ch,hest,hd);
fmri_svbfile(ySA, hAvgFile); % save in selavg format %

%fprintf(1,'--- Saving dof File ------\n');
%fprintf(1,'   %s\n',dofFile);

%fid=fopen(deblank(dofFile),'w');
%if( fid == -1 )
%  msg = sprintf('Could not open dof file %s\n',dofFile);
%  qoe(msg);
%  error(msg);
%end
%fprintf(fid,'%5d %5d %5d\n',[ [0:nNNCond]; dof; dof-1]);
%fclose(fid);

if(UsePercent)
  fprintf(1,'Computing Percent Signal Change\n');
  ind0 = find(SliceMean==0);
  l0 = length(ind0);
  if(l0 ~= 0)
    indnz = find(SliceMean ~= 0);
    fprintf(1,'Found %d voxels with mean zero\n',length(ind0));
    bmin = min(abs(SliceMean(indnz)));
    SliceMean(ind0) = bmin;
    fprintf(1,'Resetting zero mean voxels with %g\n',bmin);
  end

  hest_orig = hest;
  voxvar_orig = voxvar;

  hest   = 100*(hest ./ repmat(SliceMean, [1 1 nTotEst]));
  voxvar = 100*100*(voxvar ./ (SliceMean.*SliceMean));

  fprintf(1,'--- Saving percent signal change to\n');
  fprintf(1,'      %s\n',PSCFile);
  [ySA dof] = fmri_sxa2sa(voxvar,ch,hest,hd);
  fmri_svbfile(ySA, PSCFile); % save in selavg format %

  fmri_svdat2(PSCDatFile,hd);

  %fid=fopen(deblank(PSCDOFFile),'w');
  %if( fid == -1 )
  %  msg = sprintf('Could not open dof file %s\n',PSCDOFFile);
    %qoe(msg);
  %  error(msg);
  %end
  %fprintf(fid,'%5d %5d %5d\n',[ [0:nNNCond]; dof; dof-1]);
  %fclose(fid);

end

fprintf(1,'-------- fmri_sxaslice: Done ----------\n');
