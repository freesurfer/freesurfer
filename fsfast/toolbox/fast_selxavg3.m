% fast_selxavg3.m
% $Id: fast_selxavg3.m,v 1.2 2006/11/16 05:57:33 greve Exp $

% analysis
% sessdir

monly = 1;
analysis = 'edp2';
sess  = 'tl20000621';

flac0 = fast_ldanaflac(analysis);
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

% Load the brain mask
mask = MRIread(flac0.maskfspec);
if(isempty(mask))
  fprintf('ERROR: cannot load %s\n',flac0.maskfspec);
  return;
end
indmask = find(mask.vol);
nmask = length(indmask);
nvox = prod(mask.volsize);
fprintf('Found %d/%d (%g) voxels in mask\n',nmask,nvox,100*nmask/nvox);
mri = mask; % save as template

nruns = size(flac0.runlist,1);
fprintf('nruns = %d\n',nruns);

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

% Load in the raw data
tic;
y = [];
for nthrun = 1:nruns
  fprintf('run %d    %g\n',nthrun,toc);
  flac = flac0;
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  yrun = MRIread(flac.funcfspec);
  if(isempty(yrun))
    fprintf('ERROR: loading %s\n',funcfspec);
    return;
  end
  y = [y; fast_vol2mat(yrun.vol)];
  clear yrun;
end

fprintf('Performing OLS GLM Fit  '); tic;
[beta0 rvar0 vdof r0] = fast_glmfitw(y,X);
fprintf(' ... done %g\n',toc);

baseline = mri;
baseline.vol = fast_mat2vol(mean(beta0(ind0,:)),mri.volsize);

% Compute the ar1 for each run separately
ar1 = mri;
%ar2 = mri; % Not really ar2, just the 2nd lag
for nthrun = 1:nruns
  indrun = find(tpindrun == nthrun);
  rrun = r0(indrun,:);
  sserun = sum(rrun.^2);
  ar1run = sum(rrun(1:end-1,:).*rrun(2:end,:))./sserun;
  ar1.vol(:,:,:,nthrun) = fast_mat2vol(ar1run,ar1.volsize);
  %ar2run = sum(rrun(1:end-2,:).*rrun(3:end,:))./sserun;
  %ar2.vol(:,:,:,nthrun) = fast_mat2vol(ar2run,ar2.volsize);
end
% Compute ar1 averaged across runs
ar1mn = mri;
ar1mn.vol = mean(ar1.vol,4);
%ar2mn = mri;
%ar2mn.vol = mean(ar2.vol,4);

% Segment based on autocorrelation AR1
nacfsegs = 20;
[edge bincenter binmap] = fast_histeq(ar1mn.vol(indmask),nacfsegs);
acfseg = mri;
acfseg.vol = zeros(acfseg.volsize);
acfseg.vol(indmask) = binmap;

% Compute average ar1 in each seg and corresponding acf
nn = [1:ntptot]';
clear ar1segmn acfsegmn;
for nthseg = 1:nacfsegs
  indseg = find(acfseg.vol==nthseg);
  ar1segmn(nthseg) = mean(ar1mn.vol(indseg));
  acfsegmn(:,nthseg) = ar1segmn(nthseg).^(nn-1);
end
    
fprintf('Performing GLS GLM Fit  '); tic;
[beta rvar vdof r] = fast_glmfitw(y,X,acfsegmn,acfseg.vol(:));
fprintf(' ... done %g\n',toc);

baseline2 = mri;
baseline2.vol = fast_mat2vol(mean(beta(ind0,:)),mri.volsize);

% The contrast matrices were originally computed assuming
% only a single run's worth of nuisance regressors. Recompute.
% And compute the contrasts.
ncontrasts = length(flac0.con);
for nthcon = 1:ncontrasts
  indtask = flac_taskregind(flac0);
  C = flac0.con(nthcon).C;
  Ctask = C(:,indtask);
  [J K] = size(C);
  ZC = zeros(J,nNuis);
  C = [Ctask ZC];
  flac0.con(nthcon).C = C;
  eff = 1/trace(C*inv(X'*X)*C');
  vrf = 1/mean(diag(C*inv(X'*X)*C'));
  fprintf('%2d %-10s   eff = %6.1f   vrf = %6.1f\n',...
	  nthcon,flac0.con(nthcon).name,eff,vrf);
  [F dof1 dof2 ces cescvm] = fast_fratiow(beta,X,rvar,C,acfsegmn,acfseg.vol(:));
  p = FTest(dof1, dof2, F);
  ind = find(p == 0); p(ind) = 1;
  sig(:,:,:,nthcon) = -log10(fast_mat2vol(p,mri.volsize));
end

% Save contrasts
% Convert beta to sxa, save


