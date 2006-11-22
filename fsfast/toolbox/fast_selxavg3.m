% fast_selxavg3.m
% $Id: fast_selxavg3.m,v 1.11 2006/11/22 00:15:01 greve Exp $

% Synth - noise and signal
% MultiVar t
% Why does analysis polyfit have too hi an order?
% Inorming residual when not whitening?
% Force choice autostim/noautostim 
% Force choice on whitening or not
% Force choice of mask

% JK?
% Per-run?
% Save LUT for ACF Seg

% analysis
% sessdir

monly       = 1;
DoGLMFit    = 1;
DoContrasts = 1;
DoSynth     = 0;
SynthSeed   = -1;

analysis = '';
flacname = '';
outtop = '';

%sess  = 'tl20000621';
%flacname = 'flac/edp.flac';
%analysis = 'edp2';
%analysis = 'main3-fir';
% outtop = '/space/greve/1/users/greve/workmem-ana';

sess  = 'dng';
flacname = 'flac/rest.flac';

%sess = '/autofs/space/annecy_014/users/kdevaney/subj27/color_orientation_20060403';
%analysis = 'subj27_color';
%outtop = '/space/greve/1/users/greve/kd';

sessname = basename(sess);
%outtop = dirname(sess);
if(isempty(outtop)) outtop = '.'; end

ext = getenv('FSF_OUTPUT_FORMAT');
if(isempty(ext)) ext = 'bhdr'; end
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

flac0.sess = sess;
flac0.nthrun = 1;
flac0 = flac_customize(flac0);
if(isempty(flac0)) 
  if(~monly) quit; end
  return; 
end
outanadir = sprintf('%s/%s/%s/%s-new',outtop,sessname,flac0.fsd,flac0.name);
xfile = sprintf('%s/X.mat',outanadir);
outresdir = sprintf('%s/res',outanadir);

nruns = size(flac0.runlist,1);
%nruns = 2;
fprintf('nruns = %d\n',nruns);
ncontrasts = length(flac0.con);

if(DoSynth)
  if(SynthSeed < 0) SynthSeed = sum(100*clock); end
  fprintf('SynthSeed = %10d\n',SynthSeed);
end
yrun_randn = [];


% Load the brain mask
mask = MRIread(flac0.maskfspec);
if(isempty(mask))
  fprintf('ERROR: cannot load %s\n',flac0.maskfspec);
  %return;
end
if(0)
fname = sprintf('%s/bold/005/f.bhdr',sess);
mri = MRIread(fname);
mask = mri;
mask.vol = ones(mri.volsize);
end

indmask = find(mask.vol);
nmask = length(indmask);
nslices = mask.volsize(3);
nvox = prod(mask.volsize);
fprintf('Found %d/%d (%4.1f) voxels in mask\n',nmask,nvox,100*nmask/nvox);
mri = mask; % save as template

% Create a volume with vox val = the slice 
svol = zeros(mri.volsize);
for s = 1:nslices,  svol(:,:,s) = s; end

%---------------------------------------------%
fprintf('Creating Design Matrix\n');
Xt = [];
Xn = [];
tpindrun = [];
tic;
for nthrun = 1:nruns
  flac = flac0;
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  if(isempty(flac)) 
    if(~monly) quit;  end
    return; 
  end

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
  if(nthrun == 1) reg0 = zeros(1,nTaskRun); end
  reg0run = zeros(1,nNuisRun);
  reg0run(1) = 1;
  reg0 = [reg0 reg0run];
end
% These are the actual indices of the mean regressors
ind0 = find(reg0);

% Create the full design matrix
X = [Xt Xn];
nTask = size(Xt,2);
nNuis = size(Xn,2);
nX = size(X,2);
ntptot = size(X,1);
DOF = ntptot - nX;
B0 = inv(X'*X)*X';
Ctask = [eye(nTask) zeros(nTask,nNuis)];
nn = [1:ntptot]';
R = eye(ntptot) - X*inv(X'*X)*X';

fprintf('Computing compensation for resdual AR1 bias\n');
[rfm.M rfm.rrho1 rfm.nrho1 rfm.nrho1hat] = fast_rfm2nrho1(R);
fprintf('AR1 Correction M: %g %g\n',rfm.M(1),rfm.M(2));

%---------------------------------------------%
% The contrast matrices were originally computed assuming
% only a single run's worth of nuisance regressors. Recompute.
flacC = flac0;
for nthcon = 1:ncontrasts
  indtask = flac_taskregind(flac0);
  C = flacC.con(nthcon).C;
  C = C(:,indtask);
  for nthrun = 1:nruns
    flac = flac_customize(flac0);
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
  fprintf('%2d %-10s   eff = %6.1f   vrf = %6.1f\n',...
	  nthcon,flac0.con(nthcon).name,eff,vrf);
end

if(DoGLMFit)
  % Make the output dirs
  err = mkdirp(outanadir);
  if(err) return; end
  err = mkdirp(outresdir);      
  if(err) return; end

  % First pass thru the data to compute beta
  fprintf('OLS Beta Pass \n');
  tic;
  betamat0 = zeros(nX,nvox);
  for nthrun = 1:nruns
    fprintf('  run %d    t=%4.1f\n',nthrun,toc);
    flac = flac0;
    flac.nthrun = nthrun;
    flac = flac_customize(flac);
    if(~DoSynth)
      yrun = MRIread(flac.funcfspec);
      yrun = fast_vol2mat(yrun);
    else
      yrun_randn(:,nthrun) = randn('state'); % save state
      yrun = MRIread(flac.funcfspec,1);
      yrun = randn(yrun.nframes,nvox);
    end
    if(isempty(yrun))
      fprintf('ERROR: loading %s\n',funcfspec);
      return;
    end
    indrun = find(tpindrun == nthrun);
    Brun = B0(:,indrun);
    betamat0 = betamat0 + Brun*yrun;
    clear yrun;
    pack;
  end
  
  % Compute baseline
  betamn0 = mean(betamat0(ind0,:));
  % baseline0 = mri;
  % baseline0.vol = fast_mat2vol(betamn0,mri.volsize);
  
  % Compute Rescale Factor
  if(~isempty(flac0.inorm))
    gmean = mean(betamn0(indmask));
    RescaleFactor = flac0.inorm/gmean;
  else
    RescaleFactor = 1;
  end
  if(DoSynth) RescaleFactor = 1; end
  fprintf('RescaleFactor = %g\n',RescaleFactor);
  
  % Second pass thru the data to compute residual
  fprintf('OLS Residual Pass \n');
  tic;
  rsse = 0;
  rho1 = mri;
  rho2 = mri; % not really rho2
  for nthrun = 1:nruns
    fprintf('  run %d    t=%4.1f\n',nthrun,toc);
    flac = flac0;
    flac.nthrun = nthrun;
    flac = flac_customize(flac);
    if(~DoSynth)
      yrun = MRIread(flac.funcfspec);
      yrun = fast_vol2mat(yrun);
    else
      randn('state',yrun_randn(:,nthrun))
      yrun = MRIread(flac.funcfspec,1);
      yrun = randn(yrun.nframes,nvox);
    end
    indrun = find(tpindrun == nthrun);
    Xrun = X(indrun,:);
    yhatrun = Xrun*betamat0;
    rrun = yrun - yhatrun;
    rsserun = sum(rrun.^2);
    rsse = rsse + rsserun;
    rho1run = sum(rrun(1:end-1,:).*rrun(2:end,:))./rsserun;
    rho1.vol(:,:,:,nthrun) = fast_mat2vol(rho1run,rho1.volsize);
    %rho2run = sum(rrun(1:end-2,:).*rrun(3:end,:))./rsserun;
    %rho2.vol(:,:,:,nthrun) = fast_mat2vol(rho2run,rho2.volsize);
    if(flac0.acfbins == 0)
      fprintf('WARNING: unwhitened residuals are not intensity norm\n');
      fname = sprintf('%s/res-%03d.%s',outresdir,nthrun,ext);
      rrunmri = mri;
      rrunmri.vol = fast_mat2vol(rrun,mri.volsize);
      MRIwrite(rrunmri,fname);
    end
    clear yrun rrun;
    pack;
  end
  % Residual variance
  rvarmat0 = rsse/DOF;
  
  % Save input mask
  fname = sprintf('%s/mask.%s',outanadir,ext);
  MRIwrite(mask,fname);
  
  % Save AR1 maps
  fname = sprintf('%s/rho1.%s',outanadir,ext);
  MRIwrite(rho1,fname);
  rho1mn = mri;
  rho1mn.vol = mean(rho1.vol,4);
  fname = sprintf('%s/rho1mn.%s',outanadir,ext);
  MRIwrite(rho1mn,fname);

  % Apply bias correction
  nrho1mn = mri;
  nrho1mn.vol = rfm.M(1) + rfm.M(2)*rho1mn.vol;
  clear rho1 rho1mn;

  % Save AR2 maps
  %fname = sprintf('%s/rho2.%s',outanadir,ext);
  %MRIwrite(rho2,fname);
  %rho2mn = mri;
  %rho2mn.vol = mean(rho2.vol,4);
  %fname = sprintf('%s/rho2mn.%s',outanadir,ext);
  %MRIwrite(rho2mn,fname);
  %nrho = rho2mn.vol./rho1mn.vol;
  %nalpha = rho1mn.vol./nrho;
  %indtmp = find(nrho > 1);
  %nrho(indtmp) = 0;
  %nalpha(indtmp) = 0;
  
  % ---------------------------------------------------
  % Segment based on autocorrelation AR1
  acfseg = [];
  acfseg.vol = [];
  acfsegmn = [];
  nrho1segmn = [];
  if(flac0.acfbins > 0)
    fprintf('Whitening\n');

    acfseg = mri;
    acfseg.vol = zeros(acfseg.volsize);
    if(1)
      a = betamn0(indmask);
      a = a/std(a);
      b = rvarmat0(indmask);
      b = b/std(b);
      c = nrho1mn.vol(indmask)';
      ykm = [a; b; c];
      kmeans0 = rand(3,flac0.acfbins);
      [kmeans kmap] = fast_kmeans(ykm,flac0.acfbins,kmeans0);
      acfseg.vol(indmask) = kmap;
    else
      [edge bincenter binmap] = fast_histeq(nrho1mn.vol(indmask), flac0.acfbins);
      %[edge bincenter binmap] = fast_histeq(betamn0(indmask), flac0.acfbins);
      acfseg.vol(indmask) = binmap;
    end
    fname = sprintf('%s/acfseg.%s',outanadir,ext);
    MRIwrite(acfseg,fname);

    clear rvarmat0 betamat0;
    
    % Compute average ar1 in each seg and corresponding acf
    fprintf('Computing whitening matrices\n');
    tic;
    clear rho1segmn nrho2segmn nalphasegmn acfsegmn S Sinv W;
    for nthseg = 1:flac0.acfbins
      indseg = find(acfseg.vol==nthseg);
      nsegvox = length(indseg);
      nrho1segmn(nthseg) = mean(nrho1mn.vol(indseg));
      %nrho1segmn(nthseg) = 0; % No whitening
      fprintf('  seg  %2d  %5d  nrho1 = %5.3f (t=%4.1f)\n',....
	      nthseg,nsegvox,nrho1segmn(nthseg),toc);
      acfsegmn(:,nthseg) = nrho1segmn(nthseg).^(nn-1);
      S(:,:,nthseg) = toeplitz(acfsegmn(:,nthseg));
      Sinv(:,:,nthseg) = inv(S(:,:,nthseg));
      W(:,:,nthseg) = inv(chol(S(:,:,nthseg))');
    end

    % First pass thru the data to compute beta
    fprintf('GLS Beta Pass \n');
    tic;
    betamat = zeros(nX,nvox);
    for nthrun = 1:nruns
      fprintf('  run %d    t=%4.1f\n',nthrun,toc);
      flac = flac0;
      flac.nthrun = nthrun;
      flac = flac_customize(flac);
      if(~DoSynth)
	yrun = MRIread(flac.funcfspec);
	yrun = fast_vol2mat(yrun);
      else
	randn('state',yrun_randn(:,nthrun))
	yrun = MRIread(flac.funcfspec,1);
	yrun = randn(yrun.nframes,nvox);
      end
      yrun = RescaleFactor*yrun;  
      indrun = find(tpindrun == nthrun);
      for nthseg = 0:  flac0.acfbins
	%fprintf('     seg  %d    %g    ---------\n',nthseg,toc);
	indseg = find(acfseg.vol==nthseg);
	if(nthseg == 0)  B = B0;
	else      B = inv(X'*Sinv(:,:,nthseg)*X)*X'*Sinv(:,:,nthseg);
	end
	Brun = B(:,indrun);
	betamat(:,indseg) = betamat(:,indseg) + Brun*yrun(:,indseg);
      end
      clear yrun;
      pack;
    end
    
    % Second pass thru the data to compute beta
    fprintf('GLS Residual Pass \n');
    err = mkdirp(outresdir);
    if(err) return; end
    tic;
    rsse = 0;
    for nthrun = 1:nruns
      fprintf('  run %d    t=%4.1f\n',nthrun,toc);
      flac = flac0;
      flac.nthrun = nthrun;
      flac = flac_customize(flac);
      if(~DoSynth)
	yrun = MRIread(flac.funcfspec);
	yrun = fast_vol2mat(yrun);
      else
	randn('state',yrun_randn(:,nthrun))
	yrun = MRIread(flac.funcfspec,1);
	yrun = randn(yrun.nframes,nvox);
      end
      yrun = RescaleFactor*yrun;
      indrun = find(tpindrun == nthrun);
      Xrun = X(indrun,:);
      yhatrun = Xrun*betamat;
      rrun = yrun - yhatrun;
      for nthseg = 1:flac0.acfbins % ok to skip 0
	indseg = find(acfseg.vol==nthseg);
	acfrun = nrho1segmn(nthseg).^([0:flac.ntp-1]');
	Wseg = inv(chol(toeplitz(acfrun))');
	rrun(:,indseg) = Wseg*rrun(:,indseg);
      end
      rsserun = sum(rrun.^2);
      rsse = rsse + rsserun;
      
      fname = sprintf('%s/res-%03d.%s',outresdir,nthrun,ext);
      rrunmri = mri;
      rrunmri.vol = fast_mat2vol(rrun,mri.volsize);
      MRIwrite(rrunmri,fname);
      
      clear yrun rrun rrunmri;
      pack;
    end
    rvarmat = rsse/DOF;
  else
    fprintf('Not Whitening\n');
    rvarmat = rvarmat0*(RescaleFactor.^2);
    betamat = betamat0*RescaleFactor;
    clear rvarmat0 betamat0;
  end % acfbins > 0

  save(xfile,'X','flac0','RescaleFactor',...
       'rfm','acfseg','acfsegmn','nrho1segmn',...
       'DoSynth','SynthSeed','yrun_randn');
  
  baseline = mri;
  baseline.vol = fast_mat2vol(mean(betamat(ind0,:)),mri.volsize);
  fname = sprintf('%s/h-offset.%s',outanadir,ext);
  MRIwrite(baseline,fname);

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
end % DoGLMFit

if(DoContrasts)
  
  if(~DoGLMFit)
    load(xfile);
  
    fname = sprintf('%s/h-offset',outanadir);
    baseline = MRIread(fname);
    if(isempty(baseline)) return; end

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
  fprintf('Computing contrasts\n');
  for nthcon = 1:ncontrasts
    C = flacC.con(nthcon).C;
    [J K] = size(C);
    if(J==1)
      [Fmat dof1 dof2 cesmat cesvarmat] = ...
	  fast_fratiow(betamat,X,rvarmat,C,acfsegmn,acfseg.vol(:));
    else
      [Fmat dof1 dof2 cesmat ] = ...
	  fast_fratiow(betamat,X,rvarmat,C,acfsegmn,acfseg.vol(:));
    end
    pmat = FTest(dof1, dof2, Fmat);
    ind = find(pmat == 0); pmat(ind) = 1;
    fsigmat = -log10(pmat);
    
    outcondir = sprintf('%s/%s',outanadir,flacC.con(nthcon).name);
    err = mkdirp(outcondir);
    if(err) return; end

    ces = mri;
    ces.vol = fast_mat2vol(cesmat,mri.volsize);
    fname = sprintf('%s/ces.%s',outcondir,ext);
    MRIwrite(ces,fname);
    
    fsig = mri;
    fsig.vol = fast_mat2vol(fsigmat,mri.volsize);
    fname = sprintf('%s/fsig.%s',outcondir,ext);
    MRIwrite(fsig,fname);

    if(J == 1)
      t = mri;
      t.vol = sqrt(fsig.vol) .* sign(ces.vol);
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
    end
  
  end
end % DoContrasts

%------------------------------------------------------%
if(~isempty(analysis) & ~DoGLMFit)

% Construct selxavg-style h.dat strucutre for backwards compat
SumXtX = Ctask*X'*X*Ctask';
NTaskAvgs = nTask;
eres_std = sqrt(rvarmat);
evtaskind = flac_evtaskind(flac0);
evtask1 = flac0.ev(evtaskind(1));
hd = fmri_hdrdatstruct;
hd.TR = flac0.TR;
hd.TER = evtask1.psdwin(3);
hd.TimeWindow = evtask1.psdwin(2)-evtask1.psdwin(1);
hd.TPreStim   = -evtask1.psdwin(1);
hd.Nc = length(evtaskind)+1;
Nc = hd.Nc;
hd.Nnnc = hd.Nc - 1;
hd.Nh = nTask/hd.Nnnc;
Navgs_per_cond = hd.Nh;
hd.DOF = dof2;
hd.Npercond = 0; % N presentations per cond, who cares?
hd.Nruns = nruns;
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
hd.runlist       = [1:nruns];
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

return;
end
