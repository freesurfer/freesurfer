function err = fast_selxavg3b(configfile)
% This is a new version of fast_selxavg3 with the only change being
% that it is written as a function that takes a config file
% argument. This can be compiled for use with the matlab run-time
% compiler. It (either compiled or uncompiled version) is run by
% selxavg3-sess instead of  fast_selxavg3.m
% To compile run
%   unalias matlab
%   set path = ( /usr/pubsw/common/matlab/V.V/bin/ $path )
%   mcc  -m -v -R -singleCompThread fast_selxavg3b.m
%   cp fast_selxavg3b $DEV/fsfast/bin/fast_selxavg3b.glnxa64
%


%
% fast_selxavg3.m
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

fprintf('starting fast_selxavg3b\n');

monly = 0;
err = 1;
fid = fopen(configfile,'r');
if(fid == -1) 
  fprintf('ERROR: cannot open %s\n',configfile);
  return; 
end
while(1)
  % scroll through any blank lines or comments %
  while(1)
    tline = fgetl(fid);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end
  % Read the key %  
  key = sscanf(tline,'%s',1);
  tlinesplit = splitstring(tline);
  switch(key)
   case 'sxa3pwd',  sxa3pwd = deblank(tlinesplit(2,:));
   case 'sxa3cmd',  sxa3cmd = deblank(tlinesplit(2,:));
   case 'okfile',   
    okfile  = deblank(tlinesplit(2,:));
    if(fast_fileexists(okfile)) delete(okfile); end
   case 'monly',  monly = sscanf(tline,'%*s %d',1);
   case 'perrun', perrun = sscanf(tline,'%*s %d',1);
   case 'jkrun',  jkrun = sscanf(tline,'%*s %d',1);
   case 'DoGLMFit', DoGLMFit = sscanf(tline,'%*s %d',1);
   case 'DoContrasts', DoContrasts = sscanf(tline,'%*s %d',1);
   case 'DoSynth', DoSynth = sscanf(tline,'%*s %d',1);
   case 'SynthNoiseAmp', SynthNoiseAmp = sscanf(tline,'%*s %f',1);
   case 'SynthNoiseAR1', SynthNoiseAR1 = sscanf(tline,'%*s %f',1);
   case 'SynthSignalAmp', SynthSignalAmp = sscanf(tline,'%*s %f',1);
   case 'SynthSeed', SynthSeed = sscanf(tline,'%*s %f',1);
   case 'sess', sess = sscanf(tline,'%*s %s',1);
   case 'ReduceToMask', ReduceToMask = sscanf(tline,'%*s %d',1);   
   case 'analysis', analysis = sscanf(tline,'%*s %s',1);
   case 'flacname', flacname = sscanf(tline,'%*s %s',1);
   case 'outtop', outtop = sscanf(tline,'%*s %s',1);
   case 'DoFWHM', DoFWHM = sscanf(tline,'%*s %d',1);
   case 'MatlabSaveRes', MatlabSaveRes = sscanf(tline,'%*s %d',1);
   case 'MatlabSaveYHat', MatlabSaveYHat = sscanf(tline,'%*s %d',1);
   case 'SaveResUnwhitened', SaveResUnwhitened = sscanf(tline,'%*s %d',1);
   case 'ConList', ConList = tlinesplit(2:end,:);
  end % switch
end % while (1)
fclose(fid);

sessname = basename(sess);
fprintf('\n');
fprintf('\n');
fprintf('#@# %s ###############################\n',sessname);
fprintf('%s\n',sess);
fprintf('-------------------------\n');
fprintf('fast_selxavg3b.m @FS_VERSION@\n');
fprintf('-------------------------\n');

if(isempty(outtop)) outtop = fast_dirname(sess); end
fprintf('outtop = %s\n',outtop);

SUBJECTS_DIR = getenv('SUBJECTS_DIR');
FSHOME = getenv('FREESURFER_HOME');

dof2 = 0; % in case there are no contrasts

ext = getenv('FSF_OUTPUT_FORMAT');
if(isempty(ext)) ext = 'nii'; end
fprintf('Extension format = %s\n',ext);

if(~isempty(analysis))
  flac0 = fast_ldanaflac(analysis);
else
  flac0 = fast_ldflac(flacname);
end
if(isempty(flac0))
  if(~monly) quit; end
  return; 
end
flac0.sxaversion = 'fast_selxavg3b.m @FS_VERSION@';

% remove non-mask when analyzing. This does not change the results
% at all, it just prevents the processing of voxels that are
% already zero. It reduces the memory load and speeds things
% up. For the surface, it creates identical results as when
% reduction is not done. For the volume, everything is the same
% except the fsnr and rho volumes. There are a few voxels at the
% edge of the volume that are not masked out when redction is not
% done. I have no idea why.
flac0.ReduceToMask = ReduceToMask; 

flac0.sess = sess;
flac0.nthrun = 1;
flac0 = flac_customize(flac0);
if(isempty(flac0)) 
  if(~monly) quit; end
  return; 
end
nruns = size(flac0.runlist,1);
fprintf('nruns = %d\n',nruns);
ncontrasts = length(flac0.con);
fprintf('autostimdur = %d\n',flac0.autostimdur);

outanadir0 = sprintf('%s/%s/%s/%s',outtop,sessname,flac0.fsd,flac0.name);

cputime0 = cputime;

%%%%%%%
if(perrun | jkrun) outer_runlist = [1:nruns];
else outer_runlist = 1;
end

for nthouter = outer_runlist

  outerrun = flac0.runlist(nthouter,:);
  if(perrun)      
    outanadir = sprintf('%s/pr%s',outanadir0,outerrun);
    nthrunlist = nthouter;
  elseif(jkrun)  
    outanadir = sprintf('%s/jk%s',outanadir0,outerrun);
    nthrunlist = setxor(nthouter,[1:nruns]);
  else            
    outanadir = outanadir0;
    nthrunlist = [1:nruns];
  end

  fprintf('\n\n');
  fprintf('outanadir = %s\n',outanadir);
  errmkd = mkdirp(outanadir);
  if(errmkd) 
    err = 1;
    return; 
  end

  LF = sprintf('%s/selxavg3.log',outanadir);
  fplf = fopen(LF,'w');
  fprintf(fplf,'cd %s\n',sxa3pwd);
  fprintf(fplf,'%s\n',sxa3cmd);
  
  xfile   = sprintf('%s/X.mat',outanadir);
  outresdir = sprintf('%s/res',outanadir);

  % Load and customize all the flacs
  clear runflac;
  allpar = [];
  for nthrun = nthrunlist
    flac = flac0;
    flac.nthrun = nthrun;
    flac = flac_customize(flac);
    if(isempty(flac)) 
      fprintf(fplf,'ERROR: loading flac\n');
      if(~monly) quit;  end
      return; 
    end
    runflac(nthrun).flac = flac;
    allpar = [allpar; runflac(nthrun).flac.par];
  end

  if(length(allpar)>0)
    % Make sure condition ids are contiguous starting at 1
    condidlist = unique(allpar(:,2));
    condidlist = setdiff(condidlist,0);
    fprintf('parfiles condition id list: ');
    fprintf('%2d ',condidlist);  fprintf('\n');
    fprintf(fplf,'parfiles condition id list: ');
    fprintf(fplf,'%2d ',condidlist);  fprintf(fplf,'\n');
    if(condidlist(1) ~= 1)
      fprintf('ERROR: condition ids must be contiguous starting at 1\n');
      fprintf(fplf,'ERROR: condition ids must be contiguous starting at 1\n');
      if(~monly) quit;  end
      return; 
    end
    if(length(find(diff(condidlist)~=1)))
      fprintf('ERROR: condition ids are not contiguous\n');
      fprintf(fplf,'ERROR: condition ids are not contiguous\n');
      if(~monly) quit;  end
      return; 
    end
    if(length(condidlist) ~= flac.ana.nconditions)
      fprintf('ERROR: found %d non-null conditions, expected %d\n',...
	      length(condidlist),flac.ana.nconditions);
      fprintf(fplf,'ERROR: found %d non-null conditions, expected %d\n',...
	      length(condidlist),flac.ana.nconditions);
      if(~monly) quit;  end
      return; 
    end
  end
  
  if(~isempty(flac0.mask))
    % Load the brain mask
    if(flac0.PerSession)
      % Native space
      if(~isempty(flac0.maskfspec))
	mask = MRIread(flac0.maskfspec);
	if(isempty(mask))
	  fprintf('ERROR: cannot load %s\n',flac0.maskfspec);
	  fprintf(fplf,'ERROR: cannot load %s\n',flac0.maskfspec);
	  return;
	end
      else
	tmp = MRIread(flac0.funcfspec,1);
	mask = tmp;
	mask.vol = ones(tmp.volsize);
	clear tmp;
      end
    else
      % Non-native space, load per-run mask
      mask = flac0.mri;
      mask.vol = 0;
      for nthrun = nthrunlist
	flac = runflac(nthrun).flac;
	runmask = MRIread(flac.maskfspec);
	if(isempty(runmask)) 
	  fprintf('ERROR: cannot load %s\n',flac.maskfspec);
	  fprintf(fplf,'ERROR: cannot load %s\n',flac.maskfspec);
	  return;
	end
	mask.vol = mask.vol + runmask.vol;
      end
      if(perrun) nmv = 1;
      else       nmv = nruns;
      end
      mask.vol = (mask.vol == nmv); % Take intersection
    end
  else
    fprintf('No mask used\n');
    hdr = MRIread(flac0.funcfspec,1);
    mask = hdr;
    mask.vol = ones(hdr.volsize);
  end
  
  % Save mask
  outmaskfile = sprintf('%s/mask.%s',outanadir,ext);
  MRIwrite(mask,outmaskfile);

  indmask = find(mask.vol);
  nmask = length(indmask);
  outmaskdatfile = sprintf('%s/nmask.dat',outanadir);
  fpmask = fopen(outmaskdatfile,'w');
  fprintf(fpmask,'%d\n',nmask);
  fclose(fpmask);
  
  indmaskout = find(~mask.vol);
  nslices = mask.volsize(3);
  nvox = prod(mask.volsize);
  fprintf('Found %d/%d (%4.1f) voxels in mask %d\n',nmask,nvox,100*nmask/nvox,flac0.ReduceToMask);
  fprintf(fplf,'Found %d/%d (%4.1f) voxels in mask %d\n',nmask,nvox,100*nmask/nvox,flac0.ReduceToMask);
  mri = mask; % save as template
  mri.vol = []; % blank
  if(flac0.ReduceToMask) nData = nmask;
  else                nData = nvox;
  end

  %---------------------------------------------%
  fprintf('Creating Design Matrix\n');
  fprintf(fplf,'Creating Design Matrix\n');
  Xt = [];
  Xn = [];
  tpindrun = [];
  DoMCFit = 1;
  mcAll = [];
  tic;
  for nthrun = nthrunlist
    flac = runflac(nthrun).flac;
    
    indtask = flac_taskregind(flac);
    Xtr = flac.X(:,indtask);
    nTaskRun = size(Xtr,2);
    
    indnuis = flac_nuisregind(flac);
    Xnr = flac.X(:,indnuis);
    nNuisRun = size(Xnr,2);
    
    Za = zeros(size(Xn,1),  size(Xnr,2));
    Zb = zeros(size(Xnr,1), size(Xn,2));
    
    Xt = [Xt; Xtr];
    Xn = [Xn Za; Zb Xnr];
    
    % Keep track of which tps belong to which runs
    tpindrun = [tpindrun; nthrun*ones(flac.ntp,1)];
    
    % Keep track of which regressors are the mean
    if(nthrun == nthrunlist(1)) reg0 = zeros(1,nTaskRun); end
    reg0run = zeros(1,nNuisRun);
    reg0run(1) = 1;
    reg0 = [reg0 reg0run];

    % Fit design to motion correction regressors
    if(isempty(flac.mc)) DoMCFit = 0; end
    if(DoMCFit) 
      if(size(flac.mc,1) ~= flac.ntp) DoMCFit = 0; 
      else mcAll = [mcAll; flac.mc]; 
      end
    end
    
    % This is for comparing cross-run FFX analysis vs this analysis 
    % in which the data and design matrices are concatenated.
    % Xrun = flac.X;
    %for nthcon = 1:ncontrasts
    %  if(~isempty(ConList))
    %	ind = strmatch(flac0.con(nthcon).name,ConList);
    %if(isempty(ind)) continue; end
    %  end
    %  C = flac.con(nthcon).C;
    %  M = C*inv(Xrun'*Xrun)*C';
    %  if(nthrun == 1) conffx(nthcon).Msum = 0; end
    %  conffx(nthcon).Msum = conffx(nthcon).Msum + M;
    %end

  end % loop over runs
  fprintf(' ... creation time = %6.3f sec\n',toc);
  fprintf(fplf,' ... creation time = %6.3f sec\n',toc);

  fprintf('DoMCFit = %d\n',DoMCFit);
  fprintf(fplf,'DoMCFit = %d\n',DoMCFit);
  
  % These are the true indices of the mean regressors
  ind0 = find(reg0);

  % Create the full design matrix
  X = [Xt Xn];
  nTask = size(Xt,2);
  nNuis = size(Xn,2);
  nX = size(X,2);
  ntptot = size(X,1);
  DOF = ntptot - nX;
  if(DOF <= 0)
    for fp = [1 fplf]
      fprintf(fp,'ERROR: DOF = %d <= 0\n',DOF);
      fprintf(fp,'You have too many conditions and/or too few time points\n');
      fprintf(fp,'Number of regressors = %d\n',nX);
      fprintf(fp,'Number of time points = %d\n',ntptot);
      fprintf(fp,'%s\n',outanadir);
    end
    return;
  end
    
  doffile = sprintf('%s/dof',outanadir);
  fp = fopen(doffile,'w');
  if(fp == -1)
    fprintf('ERROR: could not open %s\n',doffile);
    return;
  end
  fprintf(fp,'%d\n',DOF);
  fclose(fp);
  fprintf('ntptot = %d, nX = %d, DOF = %d\n',ntptot,nX,DOF);
  fprintf(fplf,'ntptot = %d, nX = %d, DOF = %d\n',ntptot,nX,DOF);

  % Check condition, normalize to distinguish from scaled
  Xsss = sqrt(sum(X.^2));
  Xn = X ./ repmat(Xsss,[ntptot 1]);
  XtX = Xn'*Xn;

  tmpxfile = sprintf('%s/Xtmp.mat',outanadir);
  fprintf('Saving X matrix to %s\n',tmpxfile);
  save(tmpxfile,'X','flac','Xsss','Xn','XtX','runflac');

  XCond = cond(XtX);
  fprintf('XCond = %g (normalized)\n',XCond);
  fprintf(fplf,'XCond = %g (normalized)\n',XCond);
  if(XCond > 1e6)
    fprintf('ERROR: design is ill-conditioned\n');
    fprintf(fplf,'ERROR: design is ill-conditioned\n');
    return;
  end

  B0 = inv(X'*X)*X';
  Ctask = [eye(nTask) zeros(nTask,nNuis)];
  nn = [1:ntptot]';
  
  if(flac0.fixacf & flac0.acfbins > 0 )
    fprintf('Computing compensation for resdual AR1 bias\n');
    R = eye(ntptot) - X*inv(X'*X)*X';
    [rfm.M rfm.rrho1 rfm.nrho1 rfm.nrho1hat] = fast_rfm2nrho1(R);
    fprintf('AR1 Correction M: %g %g\n',rfm.M(1),rfm.M(2));
    fprintf(fplf,'AR1 Correction M: %g %g\n',rfm.M(1),rfm.M(2));
    clear R;
  else
    rfm.M(1) = 0;
    rfm.M(2) = 1;
  end

  %---------------------------------------------%
  % The contrast matrices were originally computed assuming
  % only a single run's worth of nuisance regressors. Recompute.
  fprintf('Computing contrast matrices\n');
  flacC = flac0;
  for nthcon = 1:ncontrasts
    if(~isempty(ConList))
      ind = strmatch(flac0.con(nthcon).name,ConList);
      if(isempty(ind)) continue; end
    end
    indtask = flac_taskregind(flac0);
    C = flacC.con(nthcon).C;
    C = C(:,indtask);
    for nthrun = nthrunlist
      flac = runflac(nthrun).flac;
      Crun = flac.con(nthcon).C;
      indnuis = flac_nuisregind(flac);    
      Cnuis = Crun(:,indnuis);
      C = [C Cnuis];
    end
    [J K] = size(C);
    flacC.con(nthcon).C = C;
    M = C*inv(X'*X)*C';
    eff = 1/trace(M);
    vrf = 1/mean(diag(M));
    flacC.con(nthcon).M = M;
    flacC.con(nthcon).eff = eff;
    flacC.con(nthcon).vrf = vrf;

    % For FFX-vs-Concat Comparison (see above)
    %Mffx = conffx(nthcon).Msum/(flac0.nruns.^2);
    %vrfffx = 1/mean(diag(Mffx));
    %fprintf('%2d %-10s J=%2d  eff = %6.1f   vrf = %8.4f vrfffx = %8.4f r = %4.2f\n',...
    %nthcon,flac0.con(nthcon).name,J,eff,vrf,vrfffx,vrf/vrfffx);
    %fprintf('%2d %-10s J=%2d  eff = %6.1f   vrf = %8.4f\n');
  end

if(DoGLMFit)
  % Make the output dirs
  errmkd = mkdirp(outanadir);
  if(errmkd) 
    err = 1;
    return; 
  end
  errmkd = mkdirp(outresdir);      
  if(errmkd) 
    err = 1;
    return; 
  end

  % First pass thru the data to compute beta
  fprintf('OLS Beta Pass \n');
  fprintf(fplf,'OLS Beta Pass \n');
  tic;
  betamat0 = zeros(nX,nData);
  gmean = 0;
  yrun_randn = [];
  rawbeta = 0;
  rawrvar = 0;
  for nthrun = nthrunlist
    fprintf('  run %d    t=%4.1f\n',nthrun,toc);
    fprintf(fplf,'  run %d    t=%4.1f\n',nthrun,toc);
    flac = runflac(nthrun).flac;          
    indrun = find(tpindrun == nthrun);
    fprintf('reading data ...'); tic
    yrun = MRIread(flac.funcfspec);
    fprintf(' %g\n',toc);
    if(isempty(yrun))
      fprintf('ERROR: loading %s\n',funcfspec);
      fprintf(fplf,'ERROR: loading %s\n',funcfspec);
      return;
    end
    if(abs(yrun.tr/1000 - flac.TR) > .01)
      fprintf('\n\n');
      fprintf('ERROR: TR mismatch between analysis and data\n');
      fprintf('analysis TR = %g, data TR = %g\n',flac.TR,yrun.tr/1000);
      if(flac.OverrideTR == 0) return; end
      fprintf('BUT you have specified to continue anyway with TR = %g.\n',flac.TR);
      fprintf('\n\n');
    end
    if(yrun.volsize(1) ~= mask.volsize(1) | ...
       yrun.volsize(2) ~= mask.volsize(2) | ...
       yrun.volsize(3) ~= mask.volsize(3))
      fprintf('ERROR: dimension mismatch between mask and %dth run\n',nthrun);
      fprintf(fplf,'ERROR: dimension mismatch between mask and %dth run\n',nthrun);
      return;
    end
    yrun = fast_vol2mat(yrun);
    if(flac0.ReduceToMask) yrun = yrun(:,indmask); end
    if(~isempty(flac.TFmtx))
      % Temporal filter
      fprintf('Temporally filtering\n');
      fprintf(fplf,'Temporally filtering\n');
      yrun = flac.TFmtx * yrun;
    end
    % Compute mean and rvar of raw data for raw SFNR
    Xdt = fast_polytrendmtx(1,size(yrun,1),1,2);
    [rawbetarun rawrvarrun] = fast_glmfit(yrun,Xdt);
    rawbeta = rawbeta + rawbetarun;
    rawrvar = rawrvar + rawrvarrun;
    
    Brun = B0(:,indrun);

    betamat0 = betamat0 + Brun*yrun;

    fprintf('    Global Mean %8.2f\n',flac.globalmean);
    fprintf(fplf,'    Global Mean %8.2f\n',flac.globalmean);
    gmean = gmean + flac.globalmean;
    clear yrun;
    %pack; % not good with matlab 7.4
  end % loop over run
  gmean = gmean/length(nthrunlist);
  
  % Compute raw sfnr
  rawbeta = rawbeta/length(nthrunlist);
  rawrvar = rawrvar/length(nthrunlist);
  indtmp = find(rawrvar == 0);
  rawsfnr = rawbeta(1,:)./sqrt(rawrvar);
  rawsfnr(indtmp) = 0;
  if(flac0.ReduceToMask) 
    tmp = zeros(1,nvox);
    tmp(indmask) = rawsfnr;
    rawsfnr = tmp;
  end
  tmpmri = mri;
  tmpmri.vol = fast_mat2vol(rawsfnr,mri.volsize);
  fname = sprintf('%s/raw.fsnr.%s',outanadir,ext);
  MRIwrite(tmpmri,fname);
  % Save mean with-in mask raw fsnr
  rawfsnrmn = mean(tmpmri.vol(indmask));
  fname = sprintf('%s/raw.fsnr.dat',outanadir);
  fp = fopen(fname,'w');
  fprintf(fp,'%f\n',rawfsnrmn);
  fclose(fp);
  
  % Compute baseline
  betamn0 = mean(betamat0(ind0,:),1);
  % baseline0 = mri;
  % baseline0.vol = fast_mat2vol(betamn0,mri.volsize);
  % Compute Rescale Factor
  if(flac0.inorm ~= 0)
    if(flac0.ReduceToMask) gmean0 = mean(betamn0);
    else                gmean0 = mean(betamn0(indmask));
    end
    fprintf('Global In-Mask Mean = %g (%g)\n',gmean,gmean0);
    fprintf(fplf,'Global In-Mask Mean = %g (%g)\n',gmean,gmean0);
    %gmean = gmean0; % This will use old method
    RescaleFactor = flac0.inorm/gmean;
    fprintf('Rescale Target = %g\n',flac0.inorm);
    fprintf(fplf,'Rescale Target = %g\n',flac0.inorm);
  else
    RescaleFactor = 1;
  end
  fprintf('RescaleFactor = %g\n',RescaleFactor);
  fprintf(fplf,'RescaleFactor = %g\n',RescaleFactor);

  betamn0  = RescaleFactor*betamn0;
  betamat0 = RescaleFactor*betamat0;
  
  % Second pass thru the data to compute residual
  fprintf('OLS Residual Pass \n');
  fprintf(fplf,'OLS Residual Pass \n');
  tic;
  rsse = 0;
  rho1 = mri; 
  rho1.vol = zeros([mri.volsize length(nthrunlist)]);
  ErrCovMtx = 0;
  sstd = [];
  for nthrun = nthrunlist
    fprintf('  run %d    t=%4.1f\n',nthrun,toc);
    fprintf(fplf,'  run %d    t=%4.1f\n',nthrun,toc);
    flac = runflac(nthrun).flac;
    indrun = find(tpindrun == nthrun);
    fprintf('reading data ...'); tic
    yrun = MRIread(flac.funcfspec);
    fprintf(' %g\n',toc);
    yrun = fast_vol2mat(yrun);
    if(flac0.ReduceToMask) yrun = yrun(:,indmask); end
    if(~isempty(flac.TFmtx))
      % Temporal filter
      fprintf('Temporally filtering\n');
      fprintf(fplf,'Temporally filtering\n');
      yrun = flac.TFmtx * yrun;
    end
    yrun = RescaleFactor*yrun;
    Xrun = X(indrun,:);
    yhatrun = Xrun*betamat0;
    rrun = yrun - yhatrun;
    if(MatlabSaveYHat)
      outyhatdir = sprintf('%s/yhat',outanadir);
      errmkd = mkdirp(outyhatdir);	
      fname = sprintf('%s/yhat-%03d.%s',outyhatdir,nthrun,ext);
      fprintf('Saving yhat to %s\n',fname);
      rrunmri = mri;
      rrunmri.vol = fast_mat2vol(yhatrun,mri.volsize);
      MRIwrite(rrunmri,fname);
    end
    
    clear yhatrun;

    tmp =  sum(rrun.^2);
    indz  = find(tmp == 0); % keep zeros from screwing stuff up
    indnz = find(tmp ~= 0); % keep zeros from screwing stuff up
    
    yrunmn = mean(yrun,1);
    yrundm = yrun - repmat(yrunmn,[size(yrun,1) 1]);
    %sstdrun = std(yrundm(:,indnz),[],2); % spatial stddev at each TP
    sstdrun = std(rrun(:,indnz),[],2); % spatial stddev at each TP
    indtp0 = find(sstdrun < 10^-5); % set tps that are 0 to 1
    sstdrun(indtp0) = 1;
    sstd = [sstd; sstdrun];
    if(flac0.HeteroGCor) 
      rrun = rrun./repmat(sstdrun,[1 nData]); 
    end

    rsserun = sum(rrun.^2);
    rsserun(indz) = max(rsserun);
    rsse = rsse + rsserun;


    if(flac0.acfsvd > 0 & flac0.acfbins > 0)
      % For rho1/ar1 calculation, remove 1st two principle components
      % from residual. This is only for the calculation of rho1/ar1
      % so that non-stationary components do not mess up the 
      % calculation of the rho1 value
      nk = flac0.acfsvd;
      fprintf('Removing %d components prior to rho1 calc\n',nk);
      [uu ss vv] = fast_svd(rrun);
      rrunpca = uu(:,nk:end)*ss(nk:end,nk:end)*vv(:,nk:end)';
      rsserunpca = sum(rrunpca.^2);
      indz = find(rsserunpca == 0); % keep zeros from screwing stuff up
      rsserunpca(indz) = max(rsserunpca);
      rho1run = sum(rrunpca(1:end-1,:).*rrunpca(2:end,:))./rsserunpca;
      clear uu ss vv rrunpca;
    else
      rho1run = sum(rrun(1:end-1,:).*rrun(2:end,:))./rsserun;
    end
    if(flac0.ReduceToMask) 
      tmp = zeros(1,nvox);
      tmp(indmask) = rho1run;
      rho1.vol(:,:,:,nthrun) = fast_mat2vol(tmp,rho1.volsize); 
    else
      rho1.vol(:,:,:,nthrun) = fast_mat2vol(rho1run,rho1.volsize); 
    end

    fprintf('Saving rho1\n');
    fname = sprintf('%s/rho1.%s',outanadir,ext);
    MRIwrite(rho1,fname);
    
    if(SaveResUnwhitened)
      fprintf('Saving unwhitened residuals\n');
      fprintf(fplf,'Saving unwhitened residuals\n');
      fname = sprintf('%s/res-uw-%03d.%s',outresdir,nthrun,ext);
      rrunmri = mri;
      if(flac0.ReduceToMask) 
	tmp = zeros(size(rrun,1),nvox);
	tmp(:,indmask) = rrun;
	rrunmri.vol = fast_mat2vol(tmp,mri.volsize);
      else 
	rrunmri.vol = fast_mat2vol(rrun,mri.volsize);
      end
      MRIwrite(rrunmri,fname);
      fprintf('Computing ACF\n');
      fprintf(fplf,'Computing ACF\n');
      acfmat = fast_acorr(rrun);
      acfmat = acfmat(1:30,:);
      acf = mri;
      if(flac0.ReduceToMask) 
	tmp = zeros(size(acfmat,1),nvox);
	tmp(:,indmask) = acfmat;
	acf.vol = fast_mat2vol(tmp,acf.volsize);
      else
	acf.vol = fast_mat2vol(acfmat,acf.volsize);
      end
      fprintf('Saving ACF\n');
      fprintf(fplf,'Saving ACF\n');
      fname = sprintf('%s/acf-uw-%03d.%s',outresdir,nthrun,ext);      
      MRIwrite(acf,fname);
    end
    if(flac0.acfbins == 0)
      %fprintf('WARNING: unwhitened residuals are not intensity norm\n');
      if(MatlabSaveRes | DoFWHM)
	fname = sprintf('%s/res-%03d.%s',outresdir,nthrun,ext);
	rrunmri = mri;
	%rrunmri.vol = fast_mat2vol(yhatrun,mri.volsize);
	rrunmri.vol = fast_mat2vol(rrun,mri.volsize);
	MRIwrite(rrunmri,fname);
      end
    end

    if(flac.fsv3_whiten)
      % Compute Err Cov Mtx from unwhitened residuals
      fprintf('Computing ErrCovMtx\n');
      fprintf(fplf,'Computing ErrCovMtx\n');
      if(flac0.ReduceToMask) ErrCovMtxRun = rrun*rrun';
      else ErrCovMtxRun = rrun(:,indmask)*rrun(:,indmask)';
      end
      if(nthrun == 1) 
	ErrCovMtx = ErrCovMtxRun;
      else
	% This should fix case where all runs dont have 
	% same number of time points by making the 
	% final matrix the size of the min # of tps
	nErrCovMtx = size(ErrCovMtx,1);
	if(flac.ntp < nErrCovMtx)
	  ErrCovMtx = ErrCovMtx(1:flac.ntp,1:flac.ntp);
	end
	if(flac.ntp > nErrCovMtx)
	  ErrCovMtxRun = ErrCovMtxRun(1:nErrCovMtx,1:nErrCovMtx);
	end
	ErrCovMtx = ErrCovMtx + ErrCovMtxRun;
      end
    end
    
    clear yrun rrun;
    %pack;
  end % loop over run
  % Residual variance
  rvarmat0 = rsse/DOF;
  
  if(flac.fsv3_whiten)
    [fsv3W fsv3AutoCor] = fast_ecvm2wmtx(ErrCovMtx);
    fsv3AutoCor(end:ntptot) = 0; % pad
  end
  
  % Save AR1 mean
  rho1mn = mri;

  % This logic is needed for Octave
  if(size(rho1.vol,4)==1) rho1mn.vol= rho1.vol;
  else                    rho1mn.vol = mean(rho1.vol,4);
  end

  rho1mn.vol(indmaskout) = 0; % mask out
  rho1mnfile = sprintf('%s/rho1mn.%s',outanadir,ext);
  MRIwrite(rho1mn,rho1mnfile);

  % Apply spatial smoothing if desired
  if(flac.acffwhm > 0 & flac0.acfbins > 0)
    fprintf('Smoothing ACF\n');
    fprintf(fplf,'Smoothing ACF\n');
    rho1mnsmfile = sprintf('%s/rho1mn.sm.%s',outanadir,ext);
    opts = sprintf('--mask %s --i %s --o %s --fwhm %f --smooth-only',...
		   outmaskfile,rho1mnfile,rho1mnsmfile,flac.acffwhm);
    if(isempty(flac0.subject)) 
      cmd = sprintf('%s/bin/mri_fwhm %s',FSHOME,opts);
    else 
      cmd = sprintf('%s/bin/mris_fwhm %s --s %s --hemi %s --sd %s',...
		    FSHOME,opts,flac0.sourcesubject,flac0.hemi,SUBJECTS_DIR);
    end
    fprintf('%s\n',cmd);
    [err rescmd] = system(cmd);
    fprintf('%s\n',rescmd);
    if(err)
      fprintf('ERROR: %s\n',cmd);
      fprintf(fplf,'ERROR: %s\n',cmd);
      return;
    end
    % Reload smoothed acf
    rho1mn = MRIread(rho1mnsmfile);
  end
  
  % Apply bias correction
  nrho1mn = mri;
  nrho1mn.vol = rfm.M(1) + rfm.M(2)*rho1mn.vol;

  % Find places where the correction pushes AR1 > 0.90
  % Probably means there is something wrong with the scan
  ind = find(abs(nrho1mn.vol) > 0.90);
  fprintf('Found %d voxels with corrected AR1 > 0.90\n',length(ind));
  fprintf(fplf,'Found %d voxels with corrected AR1 > 0.90\n',length(ind));
  if(length(ind) > 0)
    % This is a bit of a hack to "fix" them - just rescale to
    % within +/-1
    tmp = nrho1mn.vol(ind);
    nrho1mn.vol(ind) = (.9 + .1*abs(tmp)/max(abs(tmp))).*sign(tmp);
  end
  
  % Save Correct AR1 mean
  fname = sprintf('%s/nrho1mn.%s',outanadir,ext);
  MRIwrite(nrho1mn,fname);

  % ---------------------------------------------------
  % Segment based on autocorrelation AR1
  acfseg = [];
  acfseg.vol = [];
  acfsegmn = [];
  nrho1segmn = [];
  if(flac0.acfbins > 0)
    fprintf('Whitening\n');
    fprintf(fplf,'Whitening\n');
    acfseg = mri;
    acfseg.vol = zeros(acfseg.volsize);
    [edge bincenter binmap] = fast_histeq(nrho1mn.vol(indmask), flac0.acfbins);
    acfseg.vol(indmask) = binmap;
    fname = sprintf('%s/acfseg.%s',outanadir,ext);
    MRIwrite(acfseg,fname);
    clear rvarmat0 betamat0;
    
    % Compute average ar1 in each seg and corresponding acf
    fprintf('Computing whitening matrices\n');
    fprintf(fplf,'Computing whitening matrices\n');
    tic;
    clear rho1segmn nalphasegmn acfsegmn S Sinv W;
    fprintf('Alloc Sinv: %d %d %d ... ',ntptot,ntptot,flac0.acfbins);
    fprintf(fplf,'Alloc Sinv: %d %d %d ... ',ntptot,ntptot,flac0.acfbins);
    Sinv = zeros(ntptot,ntptot,flac0.acfbins);
    fprintf(' done t = %g sec\n',toc);
    fprintf(fplf,' done t = %g sec\n',toc);
    %S    = zeros(ntptot,ntptot,flac0.acfbins);
    %W    = zeros(ntptot,ntptot,flac0.acfbins);
    fname = sprintf('%s/acfsegLUT.txt',outanadir);
    fp = fopen(fname,'w');
    for nthseg = 1:flac0.acfbins
      indseg = find(acfseg.vol==nthseg);
      nsegvox = length(indseg);
      if(nsegvox == 0) continue; end
      nrho1segmn(nthseg) = mean(nrho1mn.vol(indseg));
      fprintf(fp,'%2d \t AR1(%0.2f) \t %3d %3d %3d \t %7.4f\n',...
	      nthseg,nrho1segmn(nthseg),...
	      round(255*rand),round(255*rand),round(255*rand),...
	      nrho1segmn(nthseg));
      %nrho1segmn(nthseg) = 0; % No whitening
      if(~flac.fsv3_whiten)
	acfsegmn(:,nthseg) = nrho1segmn(nthseg).^(nn-1);
      else
	acfsegmn(:,nthseg) = fsv3AutoCor;
      end
      
      fprintf('  seg  %2d  %5d  nrho1 = %5.3f\n',...
	      nthseg,nsegvox,nrho1segmn(nthseg));
      fprintf(fplf,'  seg  %2d  %5d  nrho1 = %5.3f\n',...
	      nthseg,nsegvox,nrho1segmn(nthseg));
      for nthrun = nthrunlist
	indrun = find(tpindrun == nthrun);
	nnrun = 1:runflac(nthrun).flac.ntp;
	% Have to compute acf separately for each run in case ntp changes
	if(~flac.fsv3_whiten)
	  acfsegrun = nrho1segmn(nthseg).^(nnrun-1);
	else
	  acfsegrun = fsv3AutoCor(1:runflac(nthrun).flac.ntp);
	end
	Srun = toeplitz(acfsegrun);
	Sruninv = inv(Srun);
	Sinv(indrun,indrun,nthseg) = Sruninv;
	%Wrun = inv(chol(Srun)');
	%S(indrun,indrun,nthseg) = Srun;
	%W(indrun,indrun,nthseg) = Wrun;
      end % run
    end % if acfbins > 1
    fclose(fp);

    if(flac0.HeteroGCor) X = X./repmat(sstd,[1 nX]); end

    % First pass thru the data to compute beta
    fprintf('GLS Beta Pass \n');
    fprintf(fplf,'GLS Beta Pass \n');
    tic;
    betamat = zeros(nX,nData);
    for nthrun = nthrunlist
      fprintf('  run %d    t=%4.1f\n',nthrun,toc);
      fprintf(fplf,'  run %d    t=%4.1f\n',nthrun,toc);
      flac = runflac(nthrun).flac;
      indrun = find(tpindrun == nthrun);
      yrun = MRIread(flac.funcfspec);
      yrun = fast_vol2mat(yrun);
      if(flac0.ReduceToMask) yrun = yrun(:,indmask); end      
      if(flac0.HeteroGCor) yrun = yrun./repmat(sstd(indrun),[1 nData]); end
      if(~isempty(flac.TFmtx))
	% Temporal filter
	fprintf('Temporally filtering\n');
	fprintf(fplf,'Temporally filtering\n');
	yrun = flac.TFmtx * yrun;
      end
      yrun = RescaleFactor*yrun;  
      for nthseg = 0:flac0.acfbins
	%fprintf('     seg  %d    %g    ---------\n',nthseg,toc);
	if(flac0.ReduceToMask) indseg = find(acfseg.vol(indmask)==nthseg);
	else	            indseg = find(acfseg.vol==nthseg);
	end
	if(nthseg == 0)  B = B0;
	else   B = inv(X'*Sinv(:,:,nthseg)*X)*(X'*Sinv(:,:,nthseg));
	end
	Brun = B(:,indrun);
	betamat(:,indseg) = betamat(:,indseg) + Brun*yrun(:,indseg);
      end
      clear yrun;
      %pack;
    end
    clear Sinv;
    
    % Second pass thru the data to compute beta
    fprintf('GLS Residual Pass \n');
    fprintf(fplf,'GLS Residual Pass \n');
    errmkd = mkdirp(outresdir);
    if(errmkd) 
      err = 1;
      return; 
    end
    tic;
    rsse = 0;
    for nthrun = nthrunlist
      fprintf('  run %d    t=%4.1f\n',nthrun,toc);
      flac = runflac(nthrun).flac;
      indrun = find(tpindrun == nthrun);
      yrun = MRIread(flac.funcfspec);
      yrun = fast_vol2mat(yrun);
      if(flac0.ReduceToMask) yrun = yrun(:,indmask); end
      if(flac0.HeteroGCor) yrun = yrun./repmat(sstd(indrun),[1 nData]); end
      if(~isempty(flac.TFmtx))
	% Temporal filter
	fprintf('Temporally filtering\n');
	fprintf(fplf,'Temporally filtering\n');
	yrun = flac.TFmtx * yrun;
      end
      yrun = RescaleFactor*yrun;
      Xrun = X(indrun,:);
      yhatrun = Xrun*betamat;
      rrun = yrun - yhatrun;
      for nthseg = 1:flac0.acfbins % ok to skip 0
	if(flac0.ReduceToMask) indseg = find(acfseg.vol(indmask)==nthseg);
	else	            indseg = find(acfseg.vol==nthseg);
	end
	if(~flac.fsv3_whiten)
	  acfrun = nrho1segmn(nthseg).^([0:flac.ntp-1]');
	else
	  acfrun = fsv3AutoCor(1:flac.ntp);
	end
	Wseg = inv(chol(toeplitz(acfrun))');
	rrun(:,indseg) = Wseg*rrun(:,indseg);
      end
      rsserun = sum(rrun.^2);
      rsse = rsse + rsserun;
      
      if(MatlabSaveRes | DoFWHM)
	fname = sprintf('%s/res-%03d.%s',outresdir,nthrun,ext);
	rrunmri = mri;
	if(flac0.ReduceToMask) 
	  tmp = zeros(size(rrun,1),nvox);
	  tmp(:,indmask) = rrun;
	  rrunmri.vol = fast_mat2vol(tmp,mri.volsize);
	else 
	  rrunmri.vol = fast_mat2vol(rrun,mri.volsize);
	end
	MRIwrite(rrunmri,fname);
      end
     
      clear yrun yhatrun rrun rrunmri;
    end % run list
    rvarmat = rsse/DOF;
  else
    fprintf('Not Whitening\n');
    fprintf(fplf,'Not Whitening\n');
    rvarmat = rvarmat0;
    betamat = betamat0;
    W = [];
    clear rvarmat0 betamat0;
  end % acfbins > 0

 
  % Mask or unmask betas and rvars
  if(flac0.ReduceToMask)
    % Unmask
    tmp = zeros(size(betamat,1),nvox);
    tmp(:,indmask) = betamat;
    betamat = tmp;
    tmp = zeros(1,nvox);
    tmp(indmask) = rvarmat;
    rvarmat = tmp;
  else
    % Mask
    betamat(:,indmaskout) = 0;
    rvarmat(:,indmaskout) = 0;
  end
  
  if(DoMCFit) [betamc rvarmc] = fast_glmfit(mcAll,X);
  else betamc=[]; rvarmc=[];
  end
  
  save(xfile,'X','DOF','flac0','runflac','RescaleFactor',...
       'rfm','acfseg','nrho1segmn','acfsegmn','ErrCovMtx',...
       'yrun_randn','DoMCFit','mcAll','betamc','rvarmc','sstd');

  if(DoFWHM)
    fprintf('Concatenating residuals\n');
    fprintf(fplf,'Concatenating residuals\n');
    cmd = sprintf('%s/bin/mri_concat %s/res-???.%s --o %s/all.%s',...
		  FSHOME,outresdir,ext,outresdir,ext);    
    fprintf('%s\n',cmd);
    [err rescmd] = system(cmd);
    fprintf('%s\n',rescmd);
    fprintf(fplf,'%s\n',rescmd);
    if(err)
      printf('ERROR: %s\n',cmd);
      fprintf(fplf,'ERROR: %s\n',cmd);
      return;
    end
    fprintf('Computing FWHM\n');
    opts = sprintf('--mask %s --i %s/all.%s --sum %s/fwhm.sum --dat %s/fwhm.dat',...
		   outmaskfile,outresdir,ext,outanadir,outanadir);
    if(isempty(flac0.subject)) 
      cmd = sprintf('%s/bin/mri_fwhm %s',FSHOME,opts);
    else 
      cmd = sprintf('%s/bin/mris_fwhm %s --s %s --hemi %s --sd %s',...
		    FSHOME,opts,flac0.sourcesubject,flac0.hemi,SUBJECTS_DIR);
    end
    fprintf('%s\n',cmd);
    [err rescmd] = system(cmd);
    fprintf('%s\n',rescmd);
    fprintf(fplf,'%s\n',rescmd);
    if(err)
      fprintf('ERROR: %s\n',cmd);
      fprintf(fplf,'ERROR: %s\n',cmd);
      return;
    end
    if(MatlabSaveRes == 0)
      fprintf('Deleting residuals\n');
      cmd = sprintf('rm -f %s/res-???.%s %s/all.%s',outresdir,ext,outresdir,ext);    
      fprintf('%s\n',cmd);
      [err rescmd] = system(cmd);
      fprintf('%s\n',rescmd);
      fprintf(fplf,'%s\n',rescmd);
      if(err)
	fprintf('ERROR: %s\n',cmd);
	fprintf(fplf,'ERROR: %s\n',cmd);
	return;
      end
    end
  end
  
  % Save as ascii
  xascii = sprintf('%s/X.dat',outanadir);
  fmt = [repmat('%f ',[1 size(X,2)]) '\n'];
  fp = fopen(xascii,'w');
  fprintf(fp,fmt,X');
  fclose(fp);
  
  % Save baseline as both h-offset and meanfunc
  baseline = mri;
  baseline.vol = fast_mat2vol(mean(betamat(ind0,:),1),mri.volsize);
  fname = sprintf('%s/h-offset.%s',outanadir,ext);
  MRIwrite(baseline,fname);
  fname = sprintf('%s/meanfunc.%s',outanadir,ext);
  MRIwrite(baseline,fname);
  baselinemat = fast_vol2mat(baseline.vol);

  indz  = find(baseline.vol==0);
  indnz = find(baseline.vol~=0);
  fprintf('Found %d zero-valued voxels\n',length(indz));
  fprintf(fplf,'Found %d zero-valued voxels\n',length(indz));

  beta = mri;
  beta.vol = fast_mat2vol(betamat,beta.volsize);
  fname = sprintf('%s/beta.%s',outanadir,ext);
  MRIwrite(beta,fname);

  rvar = mri;
  rvar.vol = fast_mat2vol(rvarmat,rvar.volsize);
  fname = sprintf('%s/rvar.%s',outanadir,ext);
  MRIwrite(rvar,fname);

  rstd = mri;
  rstd.vol = sqrt(rvar.vol);
  fname = sprintf('%s/rstd.%s',outanadir,ext);
  MRIwrite(rstd,fname);

  fsnr = mri;
  fsnr.vol = zeros(fsnr.volsize);
  fsnr.vol(indnz) = baseline.vol(indnz)./rstd.vol(indnz);
  fname = sprintf('%s/fsnr.%s',outanadir,ext);
  MRIwrite(fsnr,fname);
  % Save mean with-in mask fsnr
  fsnrmn = mean(fsnr.vol(indnz));
  fname = sprintf('%s/fsnr.dat',outanadir);
  fp = fopen(fname,'w');
  fprintf(fp,'%f\n',fsnrmn);
  fclose(fp);
  
  fsnr2 = fsnr;
  fsnr2.vol = fsnr.vol.^2;
  fname = sprintf('%s/fsnr2.%s',outanadir,ext);
  MRIwrite(fsnr2,fname);
  
end % DoGLMFit

if(DoContrasts)
  fprintf('Computing contrasts\n');
  fprintf(fplf,'Computing contrasts\n');
  
  if(~DoGLMFit)
    fprintf('Loading previous GLM fit\n');
    fprintf(fplf,'Loading previous GLM fit\n');
    flac0tmp = flac0; % keep copy
    load(xfile);
    flac0 = flac0tmp; % So that it has all the contrasts
  
    fname = sprintf('%s/h-offset',outanadir);
    baseline = MRIread(fname);
    if(isempty(baseline)) return; end
    indz  = find(baseline.vol==0);
    indnz = find(baseline.vol~=0);
    baseline.vol(indz) = 1e9;
    baselinemat = fast_vol2mat(baseline.vol);

    fname = sprintf('%s/beta',outanadir);
    beta = MRIread(fname);
    if(isempty(beta)) return; end
    betamat = fast_vol2mat(beta.vol);
    clear beta.vol;

    fname = sprintf('%s/rvar',outanadir);
    rvar = MRIread(fname);
    if(isempty(rvar)) return; end
    rvarmat = fast_vol2mat(rvar.vol);
    clear rvar.vol;
    
    if(flac0.acfbins == 0) acfseg.vol = []; end
  end

  %---------------------------------------------------------------%
  fprintf('Starting contrasts\n');
  fprintf(fplf,'Starting contrasts\n');
  for nthcon = 1:ncontrasts
    if(~isempty(ConList))
      ind = strmatch(flac0.con(nthcon).name,ConList);
      if(isempty(ind)) continue; end
    end
    conname = flacC.con(nthcon).name;
    C = flacC.con(nthcon).C;
    [J K] = size(C);
    fprintf('%s J=%d -------------\n',conname,J);
    fprintf(fplf,'%s J=%d -------------\n',conname,J);
    if(J==1)
      [Fmat dof1 dof2 cesmat cesvarmat pccmat] = ...
	  fast_fratiow(betamat,X,rvarmat,C,acfsegmn,acfseg.vol(:));
    else
      [Fmat dof1 dof2 cesmat] = ...
	  fast_fratiow(betamat,X,rvarmat,C,acfsegmn,acfseg.vol(:));
    end
    pmat = FTest(dof1, dof2, Fmat);
    ind = find(rvarmat == 0); pmat(ind) = 1; % in case all tps are same value
    ind = find(pmat == 0);    pmat(ind) = eps(0); % for REALLY sig voxels
    fsigmat = -log10(pmat);

    % Contrast output
    outcondir = sprintf('%s/%s',outanadir,conname);
    if(exist(outcondir,'dir'))
      % Delete it if it exists
      fprintf('%s exists, deleting\n',outcondir);
      fprintf(fplf,'%s exists, deleting\n',outcondir);
      rmdir(outcondir,'s');
    end

    errmkd = mkdirp(outcondir);
    if(errmkd) 
      err = 1;
      return; 
    end

    fprintf('Saving efficiency %g\n',flacC.con(nthcon).eff);
    fname = sprintf('%s/efficiency.dat',outcondir);
    fid = fopen(fname,'w');
    fprintf(fid,'%g\n',flacC.con(nthcon).eff);
    fclose(fid);

    ces = mri;
    ces.vol = fast_mat2vol(cesmat,mri.volsize);
    fname = sprintf('%s/ces.%s',outcondir,ext);
    MRIwrite(ces,fname);
    
    tmp = zeros(J,nvox);
    tmp = 100*cesmat./repmat(abs(baselinemat),[J 1]);
    cespct = mri;
    cespct.vol = fast_mat2vol(tmp,cespct.volsize);
    fname = sprintf('%s/cespct.%s',outcondir,ext);
    MRIwrite(cespct,fname);
    
    fsig = mri;
    fsig.vol = fast_mat2vol(fsigmat,mri.volsize);
    fname = sprintf('%s/fsig.%s',outcondir,ext);
    MRIwrite(fsig,fname);

    Fvol = mri;
    Fvol.vol = fast_mat2vol(Fmat,mri.volsize);
    fname = sprintf('%s/F.%s',outcondir,ext);
    MRIwrite(Fvol,fname);

    if(J == 1)
      t = mri;
      t.vol = sqrt(Fvol.vol) .* sign(ces.vol);
      fname = sprintf('%s/t.%s',outcondir,ext);
      MRIwrite(t,fname);
      sig = mri;
      sig.vol = fsig.vol .* sign(ces.vol);
      fname = sprintf('%s/sig.%s',outcondir,ext);
      MRIwrite(sig,fname);

      cesvar = mri;
      cesvar.vol = fast_mat2vol(cesvarmat,mri.volsize);
      fname = sprintf('%s/cesvar.%s',outcondir,ext);
      MRIwrite(cesvar,fname);
      cesvarpct = mri;
      cesvarpct.vol = (100.^2)*cesvar.vol./(baseline.vol.^2);
      fname = sprintf('%s/cesvarpct.%s',outcondir,ext);
      MRIwrite(cesvarpct,fname);
    
      if(0)
	cesstd = mri;
	cesstd.vol = fast_mat2vol(sqrt(cesvarmat),mri.volsize);
	fname = sprintf('%s/cesstd.%s',outcondir,ext);
	MRIwrite(cesstd,fname);
	cesstdpct = mri;
	cesstdpct.vol = (100.^2)*cesstd.vol./(baseline.vol.^2);
	fname = sprintf('%s/cesstdpct.%s',outcondir,ext);
	MRIwrite(cesstdpct,fname);
      end
    
      cnr = mri;
      cnr.vol = zeros(cnr.volsize);
      cnr.vol(indnz) = ces.vol(indnz) ./ sqrt(rvar.vol(indnz));
      fname = sprintf('%s/cnr.%s',outcondir,ext);
      MRIwrite(cnr,fname);
    
      pcc = mri;
      pcc.vol = fast_mat2vol(pccmat,pcc.volsize);
      fname = sprintf('%s/pcc.%s',outcondir,ext);
      MRIwrite(pcc,fname);
    end

    zscoremat = fast_p2z(pmat/2); % Div by 2 to make it one-sided
    if(J==1) zscoremat = zscoremat.*sign(cesmat); end
    zscore = mri;
    zscore.vol = fast_mat2vol(zscoremat,zscore.volsize);
    fname = sprintf('%s/z.%s',outcondir,ext);
    MRIwrite(zscore,fname);

    if(flac.IsRetinotopy | flac.IsABBlocked)
      if(strcmp(conname,'eccen') | strcmp(conname,'polar') | strcmp(conname,'fund'))
	cesreal = ces;	cesreal.vol = ces.vol(:,:,:,1);
	cesimag = ces;	cesimag.vol = ces.vol(:,:,:,2);
	if(~isempty(flac0.subject) & strcmp(flac0.hemi,'rh') & ...
	  strcmp(conname,'polar'))
	  % Make Left Hor Meridian have angle 0
	  cesreal.vol = -cesreal.vol;
	end
	mag = ces; mag.vol = sqrt(cesreal.vol.^2 + cesimag.vol.^2);
	phz = ces;
	phz.vol = atan2(cesimag.vol,cesreal.vol);
	if(strcmp(conname,'eccen'))
	  % Force eccen angle to be 0-2pi
	  ind = find(phz.vol < 0);
	  phz.vol(ind) = phz.vol(ind) + 2*pi;
	end
	fname = sprintf('%s/imag.%s',outcondir,ext);
	MRIwrite(cesimag,fname);
	fname = sprintf('%s/real.%s',outcondir,ext);
	MRIwrite(cesreal,fname);
	fname = sprintf('%s/angle.%s',outcondir,ext);
	MRIwrite(phz,fname);
	fname = sprintf('%s/mag.%s',outcondir,ext);
	MRIwrite(mag,fname);
	% For compatibility with tksurfer color wheel
	cesimag.vol = sin(phz.vol).*fsig.vol;
	fname = sprintf('%s/cwmap-imag.%s',outcondir,ext);
	MRIwrite(cesimag,fname);
	cesreal.vol = cos(phz.vol).*fsig.vol;
	fname = sprintf('%s/cwmap-real.%s',outcondir,ext);
	MRIwrite(cesreal,fname);
	% Mask angle by fsig
	indz = find(abs(fsig.vol)<2);
	phz.vol(indz) = 0;
	fname = sprintf('%s/angle.masked.%s',outcondir,ext);
	MRIwrite(phz,fname);
      end
    end
    
    if(J > 1)
      % Compute CES amplitude as the sqrt of sum of the squares
      if(J < 20) % 20 to prevent it from going crazy
	fprintf('Computing CES Magnitude\n');
	cesmag = mri;
	cesmag.vol = sqrt(sum(ces.vol.^2,4));
	fname = sprintf('%s/cesmag.%s',outcondir,ext);
	MRIwrite(cesmag,fname);
	
	cesmagpct = mri;
	cesmagpct.vol = 100*cesmag.vol./abs(baseline.vol);
	fname = sprintf('%s/cesmagpct.%s',outcondir,ext);
	MRIwrite(cesmagpct,fname);
      end
    
      tsigmatall = [];
      cesvarmatall = [];
      for nthj = 1:J
	Cj = C(nthj,:);
	[Fmat dof1 dof2 cesmat cesvarmat] = ...
	    fast_fratiow(betamat,X,rvarmat,Cj,acfsegmn,acfseg.vol(:));
	pmat = FTest(dof1, dof2, Fmat);
	ind = find(pmat == 0); pmat(ind) = 1;
	tsigmat = -log10(pmat) .* sign(cesmat);
	tsigmatall(nthj,:) = tsigmat;
	cesvarmatall(nthj,:) = cesvarmat;
      end
      tsigall = mri;
      tsigall.vol = fast_mat2vol(tsigmatall,mri.volsize);
      fname = sprintf('%s/sig.%s',outcondir,ext);
      MRIwrite(tsigall,fname);

      [sigmin rmin]  = max(abs(tsigmatall));
      indmin = sub2ind(size(tsigmatall),rmin,1:nvox);
      % log10(J) is bonf cor 
      sigmin = sign(tsigmatall(indmin)).*(sigmin - log10(J)); 
      tminsig = mri;
      tminsig.vol = fast_mat2vol(sigmin,mri.volsize);
      fname = sprintf('%s/minsig.%s',outcondir,ext);
      MRIwrite(tminsig,fname);

      rminmri = mri;
      rminmri.vol = fast_mat2vol(rmin,mri.volsize);
      fname = sprintf('%s/iminsig.%s',outcondir,ext);
      MRIwrite(rminmri,fname);
      
      cesvarall = mri;
      cesvarall.vol = fast_mat2vol(cesvarmatall,mri.volsize);
      fname = sprintf('%s/cesvar.%s',outcondir,ext);
      MRIwrite(cesvarall,fname);

      tmp = zeros(J,nvox);
      tmp = (100.^2)*cesvarmatall./repmat(baselinemat.^2,[J 1]);
      cesvarpct = mri;
      cesvarpct.vol = fast_mat2vol(tmp,cesvarpct.volsize);
      fname = sprintf('%s/cesvarpct.%s',outcondir,ext);
      MRIwrite(cesvarpct,fname);

    end % if(J>1)

    if(DoMCFit)
      fprintf('Testing Motion Correction Parameters\n');
      fprintf(fplf,'Testing Motion Correction Parameters\n');
      [F Fsig ces] = fast_fratio(betamc,X,rvarmc,C);
      Fsig = -log10(Fsig);
      minsig = max(Fsig);
      if(J==1) Fsig = sign(ces).*Fsig; end
      fprintf('%6.2f ',Fsig)
      fprintf('\n');
      % ------- sig ---------
      fname = sprintf('%s/mc.sig.dat',outcondir);
      fpMC = fopen(fname,'w');
      fprintf(fpMC,'%6.2f ',Fsig);
      fprintf(fpMC,'\n');
      fclose(fpMC);
      % ------- minsig ---------
      fname = sprintf('%s/mc.minsig.dat',outcondir);
      fpMC = fopen(fname,'w');
      fprintf(fpMC,'%6.2f\n',minsig);
      fclose(fpMC);
      % ------- ces ---------
      fname = sprintf('%s/mc.ces.dat',outcondir);
      fpMC = fopen(fname,'w');
      for nthj = 1:J
	fprintf(fpMC,'%10.6f ',ces(nthj,:));
	fprintf(fpMC,'\n');
      end
      fclose(fpMC);
    end

  end % contrast list

  if(flac.IsRetinotopy)    
    fname = sprintf('%s/eccen/fsig',outanadir);
    eccen = MRIread(fname);
    fname = sprintf('%s/polar/fsig',outanadir);
    polar = MRIread(fname);
    conj = eccen;
    conjvol = min(eccen.vol,polar.vol);
    fsdir = sprintf('%s/fieldsign',outanadir);
    err = mkdirp(fsdir);
    fname = sprintf('%s/fsig.%s',fsdir,ext);
    MRIwrite(conj,fname);
  end
    
end % DoContrasts

%------------------------------------------------------%
% Check to make sure that each task ev as the same number
% of regressors before saving the h.dat file
evtaskind = flac_evtaskind(flac0);
% number of regressors in each ev
nregperev = zeros(length(evtaskind),1);
for n = 1:length(evtaskind)
  nregperev(n) = flac0.ev(evtaskind(n)).nreg;
end
nregperevunique = length(unique(nregperev));

if(~isempty(analysis) & DoGLMFit & nregperevunique==1 & ...
   strcmp(flac.designtype,'event-related'))
  % Construct selxavg-style h.dat strucutre for backwards compat
  SumXtX = Ctask*X'*X*Ctask';
  NTaskAvgs = nTask;
  eres_std = sqrt(rvarmat);
  evtaskind = flac_evtaskind(flac0);
  evtask1 = flac0.ev(evtaskind(1));
  hd = fmri_hdrdatstruct;
  hd.TR = flac0.TR;
  if(~isempty(evtask1.psdwin)) 
    hd.TER = evtask1.psdwin(3);
    hd.TimeWindow = evtask1.psdwin(2)-evtask1.psdwin(1);
    hd.TPreStim   = -evtask1.psdwin(1);
  else
    hd.TER = flac0.TR;
    hd.TimeWindow = 0;
    hd.TPreStim   = 0;
  end
  hd.Nc = length(evtaskind)+1;
  Nc = hd.Nc;
  hd.Nnnc = hd.Nc - 1;
  hd.Nh = nTask/hd.Nnnc;
  Navgs_per_cond = hd.Nh;
  hd.DOF = dof2;
  hd.Npercond = 0; % N presentations per cond, who cares?
  hd.Nruns = length(nthrunlist);
  hd.Ntp = ntptot;
  hd.Nrows = mri.volsize(1);
  hd.Ncols = mri.volsize(2);
  hd.Nskip = 0; % Who cares?
  tmp = flac_evindex(flac0,'Poly');
  if(~isempty(tmp))
    hd.DTOrder = flac0.ev(tmp).params(1) + 1;
  else
    hd.DTOrder = 1;
  end
  hd.RescaleFactor = RescaleFactor; 
  hd.HanningRadius = 0.0;
  hd.BrainAirSeg = 0;
  hd.NullCondId    = 0;
  hd.SumXtX        = SumXtX;
  hd.nNoiseAC      = 0;
  hd.CondIdMap     = [0:hd.Nc-1];
  hCovMtx          = Ctask*inv(X'*X)*Ctask';
  hd.hCovMtx       = hCovMtx;
  hd.WhitenFlag    = 0;
  hd.runlist       = nthrunlist;
  hd.funcstem      = flac0.funcstem;
  hd.parname       = flac0.parfile;
  hd.GammaFit = 0; %?
  hd.gfDelta  = [];%?
  hd.gfTau    = [];%?
  %if(s.spmhrf > -1) % Hack
  %  hd.GammaFit = s.spmhrf + 1;
  %  hd.gfDelta = -1*ones(1,s.spmhrf + 1);
  %  hd.gfTau   = -1*ones(1,s.spmhrf + 1);
  %end;
  fname = sprintf('%s/h.dat',outanadir);
  fprintf('Saving h.dat to %s\n',fname);
  fmri_svdat3(fname,hd);
  % hd2 = fmri_lddat3('sxa.dat');

  % -------- Convert to selavg format -------------- %
  hhattmp = betamat(1:NTaskAvgs,:);
  hhattmp = [zeros(Navgs_per_cond,nvox); hhattmp]; % Add zero for cond 0
  hhattmp2 = reshape(hhattmp,[Navgs_per_cond Nc nvox]);
  hstd = sqrt( (diag(hCovMtx).*diag(SumXtX)) * eres_std.^2);
  hstdtmp = hstd(1:NTaskAvgs,:); % Remove offset and baseline
  hstdtmp = [repmat(eres_std, [Navgs_per_cond 1]); hstdtmp]; % Add 0 for cond 0
  hstdtmp2 = reshape(hstdtmp,[Navgs_per_cond Nc nvox]);
  
  %--- Merge Averages and StdDevs ---%
  tmp = zeros(Navgs_per_cond,2,Nc,nvox);
  tmp(:,1,:,:) = hhattmp2;
  tmp(:,2,:,:) = hstdtmp2;
  tmp = reshape(tmp,[Navgs_per_cond*2*Nc mri.volsize ]);
  tmp = permute(tmp,[2 3 4 1]);
  hsxa = mri;
  hsxa.vol = tmp;
  fname = sprintf('%s/h.%s',outanadir,ext);
  MRIwrite(hsxa,fname);
end % sxa format

end % outer run loop

if(exist('okfile','var'))  fmri_touch(okfile); end

fprintf('cputime-min %f\n',(cputime-cputime0)/60);
fprintf(fplf,'cputime-min %f\n',(cputime-cputime0)/60);

fprintf('fast_selxavg3 done for %s\n',sess);
fprintf(fplf,'fast_selxavg3 done for %s\n',sess);
fclose(fplf);
err = 0;

return;
