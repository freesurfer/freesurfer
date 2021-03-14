% fmri_sxavol.m
% Jan 26, 2000
% Modification of fmri_sxaslice to processes an entire volume
% 
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
% fmri_sxavol.m
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
fprintf(1,'  --- fmri_sxavol: Starting ------\n');
fprintf(1,'fmri_sxavol.m @FS_VERSION@\n');

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
VarNameList = strvcat(VarNameList,'ParFiles');

VarNameList = strvcat(VarNameList,'InputStems');
VarNameList = strvcat(VarNameList,'hAvgStem');
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

%%% Get the number of runs from the number of input stems %%%
nRuns = size(InputStems,1);

%% Compute the number of delay-points in the HDR to estimate
nHEst = floor(TimeWindow/TER);
if(nHEst < 1)
  msg = 'TimeWindow too small, must be > TER';
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
  end
end
Par = remap_par(Par,CondIdMap);

%% Check the number of conditions for each run %%
nCondPerRun = reshape1d(max(Par(:,2,:)));
fprintf(1,'Non-Null Conditions per Run: ');
fprintf(1,'%d ',nCondPerRun);
fprintf(1,'\n');

%% Check consistency of dimensions across runs %%
[NSlices nRows nCols nTP] = fmri_bvoldim(InputStems(1,:));
fprintf('Dim (%s): %d %d %d %d\n',InputStems(1,:),NSlices,nRows,nCols,nTP);
for Run = 2:nRuns
  [ns2 nr2 nc2 nt2] = fmri_bvoldim(InputStems(Run,:));
  fprintf('Dim (%s): %d %d %d %d\n',InputStems(Run,:),ns2,nr2,nc2,nt2);
  if(ns2 ~= NSlices | nRows ~= nr2 | nCols ~= nc2 | nt2 ~= nTP)
    msg = sprintf('Run %d has inconsistent dimensions\n',Run);
    qoe(msg);error(msg);
  end
end
nV    = nRows*nCols;

% Compute the sub-sampling matrix
[Mss Rss]   = fmri_subsampmat(TR,nTP,TER);
Nts = nTP*Rss;

fprintf(1,'nRows = %d, nCols = %d, nTP = %d, nRuns = %d\n',...
        nRows,nCols,nTP,nRuns);

%%%% --- Check for exclusions (including skips) ---%%%%
TPExclude   = zeros(nTP,nRuns);
ssTPExclude = zeros(Nts,nRuns);
if(nSkip > 0) tSkip = TR*[0:(nSkip-1)]';  %'
else          tSkip = [];
end

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

end

%%%--- PreStim --- %%%%
nPreStim = floor(tPreStim/TER);
fprintf('nPreStim = %d\n',nPreStim);

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
     %fmri_par2scm(Par,Nts,TER,nHEst,tPreStim,NullCondId,ssTPExclude);
  if(Rss ~= 1)
    X2 = zeros(nTP,nHEst*(Nconds-1),nRuns);
    for n = 1:nRuns,
      X2(:,:,n) = Mss * X(:,:,n);
    end
    X = X2;
    clear X2;
  end

if(GammaFit == 1 ) 
  fprintf(1,'Using Gamma Fit: Delta = %g, Tau = %g\n',gfDelta,gfTau);
  X = fmri_scm2gcm(X,nNNCond,TR,tPreStim,gfDelta,gfTau);
  nHEst = 1;
end

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
for n = 1:nRuns,
  c = cond(X2(:,:,n)'*X2(:,:,n)); %'
  fprintf(1,'Run %2d Condition XtX: %g\n',n,c);
  sumXtX = sumXtX + X2(:,:,n)' * X2(:,:,n); %'

  if(DTOrder > 0) 
    Xdetrend = [Xdetrend1 Xdetrend2];
    tmp = repmat(~TPExclude(:,n), [1 DTOrder*nRuns]);
    Xdetrend = Xdetrend .* tmp;
    X(:,:,n)  = [X2(:,:,n) Xdetrend]; % add a column of ones
    Xdetrend1 = fmri_shiftcol(Xdetrend1,1,0) ;
    Xdetrend2 = fmri_shiftcol(Xdetrend2,1,0) ;
  end
end
if(DTOrder == 0) X = X2; end

fprintf(1,'Condition SumXtX: %g\n',cond(sumXtX));
%clear X2 Xdetrend Xdetrend1 Xdetrend2 tmp;

%%%%%%%------------------------------------%%%%%%%%%%

Xorig = X;
for SliceNo = 0:NSlices-1
  fprintf('--- Processing Slice %2d ---\n',SliceNo);
  X = Xorig;

  InputFiles = [];
  for Run = 1:nRuns
    fname = sprintf('%s_%03d.bshort',deblank(InputStems(Run,:)),SliceNo);
    InputFiles = strvcat(InputFiles,fname);
  end

  %%% ---- Load the fMRI Runs ------- %%%
  if(~Synth)
    Slice = fmri_ldbfile(InputFiles);
  else
    if(Seed < 0) Seed = sum(100*clock);  end
    fprintf('Synthesizing input data with seed %d\n',Seed);
    randn('state',Seed);
    Slice = randn(nRuns,nRows,nCols,nTP);
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


  %%%%-------- Rescale --------- %%%%%%%
  if(RescaleTarget > 0)
    fprintf(1,'Rescaling to Mean = %g\n',RescaleTarget);
    [fSlice oldmean] = fmri_rescale(fSlice,RescaleTarget,TPExclude);
  else
    fprintf(1,'INFO: Skipping rescaling, RescaleTarget = %g\n',RescaleTarget);
  end
  y = fSlice;

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

  fprintf('TSliceSynch = %g\n',TSliceSynch);
  TimeOffset = TimeOffset - TSliceSynch;

  fprintf(1,'Computing hCoVar\n');
  ch          = fmri_hcovar(X);
  fprintf(1,'Selective Sum\n');
  ss          = fmri_selsum(X,y);
  fprintf(1,'LMS Fit\n');
  hest        = fmri_lmsfit(ch,ss);
  %clear ss;
  fprintf(1,'Residuals\n');
  [eest sigest] = fmri_residual(y,X,hest);

  DOF  = (nTP-DTOrder)*nRuns-nNNCond*nHEst;
  fprintf(1,'DOF = %d\n',DOF);

  fprintf(1,'VoxVar\n');
  voxvar = fmri_voxvar(eest,DOF,[],TPExclude);

  evar   = std(reshape(eest, [nRows nCols nTP*nRuns]) ,[],3).^2;
  sigvar = std(reshape(sigest, [nRows nCols nTP*nRuns]) ,[],3).^2;

  % Remove the information about the trend %%
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

  fprintf(1,'Segmenting Brain and Air\n');
  [ibrain iair thresh] = fmri_segment(hoffset);
  SliceMean = mean(hoffset,3);
  if(nRuns == 1) SliceStd = std(y,[],3);
  else           SliceStd = sqrt(mean(squeeze(std(y,[],3)).^2,3));
  end
  fprintf(1,'Segmenting Threshold %g\n',thresh);
  fprintf(1,'Brain Voxels: %d\n',length(ibrain));
  fprintf(1,'Air   Voxels: %d\n',length(iair));
  SegMask = zeros(nRows,nCols);
  SegMask(ibrain(:,1)) = 1;
  clear tmp;
  tmp(:,:,1) = SegMask;
  tmp(:,:,2) = SliceMean;
  tmp(:,:,3) = SliceStd;
  fmri_svbfile(tmp,segFile);
  %clear SegMask SlcieStd tmp fSlice;

  SigPower   = mean(reshape1d(sigvar(ibrain(:,1))));
  NoisePower = mean(reshape1d(evar(ibrain(:,1))));
  SNR1 = SigPower/NoisePower;
  SNR2 = mean(reshape1d(sigvar(ibrain(:,1))./evar(ibrain(:,1))));
  %fprintf('SigPower = %g, NoisePower = %g, SNR1 = %g, SNR2 = %g\n',...
  %        SigPower, NoisePower, SNR1, SNR2);

  if(NoiseType == 1) % Compensate for correlations in noise

    fprintf(1,'-------- Compensating for correlated noise -----\n');
    nRee = nNACorr;
    fprintf(1,'nRee = %d\n',nRee);

    fprintf(1,'Computing Error AutoCorrelation Function\n');
    if(BASegment) Ree = fmri_acorr(eest,nRee,ibrain);
    else          Ree = fmri_acorr(eest,nRee);
    end

    fprintf(1,'Fitting Average AutoCorrelation Function\n');
    [alpha rho] = fmri_acorrfit(Ree,nRee);
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
    eest        = fmri_residual(y,X,hest);

    fprintf(1,'VoxVar\n');
    voxvar      = fmri_voxvar(eest,DOF,Ce);

    clear Ce;
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

  fprintf(1,'--- Saving Slice SelXAvg ------\n');
  hAvgFile = sprintf('%s_%03d.bfloat',deblank(hAvgStem),SliceNo);
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
    PSCFile = sprintf('%s_%03d.bfloat',deblank(PSCStem),SliceNo);
    fprintf(1,'      %s\n',PSCFile);
    [ySA dof] = fmri_sxa2sa(voxvar,ch,hest,hd);
    fmri_svbfile(ySA, PSCFile); % save in selavg format %

    PSCDatFile = sprintf('%s.dat',deblank(PSCStem));
    fmri_svdat2(PSCDatFile,hd);

    %fid=fopen(deblank(PSCDOFFile),'w');
    %if( fid == -1 )
    %  msg = sprintf('Could not open dof file %s\n',PSCDOFFile);
    %  qoe(msg);
    %  error(msg);
    %end
    %fprintf(fid,'%5d %5d %5d\n',[ [0:nNNCond]; dof; dof-1]);
    %fclose(fid);

  end % if(UsePercent) %

end % Loop over slices %

fprintf(1,'-------- fmri_sxaslice: Done ----------\n');
