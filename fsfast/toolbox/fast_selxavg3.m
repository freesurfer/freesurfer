% fast_selxavg3.m
% $Id: fast_selxavg3.m,v 1.5 2006/11/17 06:28:45 greve Exp $

% Save ACF Seg Means
% Break-out contrasts?
% Save X.mat
% Allow turning off whitening

% Save final flac.
% JK?
% Per-run?
% Save LUT for ACF Seg

% analysis
% sessdir
monly = 1;

sess  = 'tl20000621';

flacname = 'flac/edpnomc.flac';
analysis = '';
%analysis = 'edp2';
%analysis = 'main3-fir';

%sess = '/autofs/space/annecy_014/users/kdevaney/subj27/color_orientation_20060403';
%analysis = 'subj27_color';
%outtop = '/space/greve/1/users/greve/kd';

sessname = basename(sess);
outtop = dirname(sess);

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

outanadir = sprintf('%s/%s/%s/%s-new',outtop,sessname,flac0.fsd,analysis);
mkdirp(outanadir);

% Load the brain mask
mask = MRIread(flac0.maskfspec);
if(isempty(mask))
  fprintf('ERROR: cannot load %s\n',flac0.maskfspec);
  %return;
end
%fname = sprintf('%s/bold/005/f.bhdr',sess);
%mri = MRIread(fname);
%mask = mri;
%mask.vol = ones(144,194,26);
indmask = find(mask.vol);
nmask = length(indmask);
nslices = mask.volsize(3);
nvox = prod(mask.volsize);
fprintf('Found %d/%d (%4.1f) voxels in mask\n',nmask,nvox,100*nmask/nvox);
mri = mask; % save as template

% Create a volume with vox val = the slice 
svol = zeros(mri.volsize);
for s = 1:nslices,  svol(:,:,s) = s; end

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
DOF = ntptot - nX;
B0 = inv(X'*X)*X';
Ctask = [eye(nTask) zeros(nTask,nNuis)];

% First pass thru the data to compute beta
fprintf('OLS Beta Pass \n');
tic;
beta0 = 0;
for nthrun = 1:nruns
  fprintf('  run %d    %g\n',nthrun,toc);
  flac = flac0;
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  yrun = MRIread(flac.funcfspec);
  yrun = fast_vol2mat(yrun);
  if(isempty(yrun))
    fprintf('ERROR: loading %s\n',funcfspec);
    return;
  end
  indrun = find(tpindrun == nthrun);
  Brun = B0(:,indrun);
  beta0 = beta0 + Brun*yrun;
  clear yrun;
  pack;
end

% Second pass thru the data to compute residual
fprintf('OLS Residual Pass \n');
tic;
rsse = 0;
ar1 = mri;
for nthrun = 1:nruns
  fprintf('  run %d    %g\n',nthrun,toc);
  flac = flac0;
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  yrun = MRIread(flac.funcfspec);
  indrun = find(tpindrun == nthrun);
  Xrun = X(indrun,:);
  yhatrun = Xrun*beta0;
  rrun = fast_vol2mat(yrun) - yhatrun;
  rsserun = sum(rrun.^2);
  rsse = rsse + rsserun;
  ar1run = sum(rrun(1:end-1,:).*rrun(2:end,:))./rsserun;
  ar1.vol(:,:,:,nthrun) = fast_mat2vol(ar1run,ar1.volsize);
  clear yrun rrun;
  pack;
end
rvar0 = rsse/DOF;

betamn0 = mean(beta0(ind0,:));
if(~isempty(flac0.inorm))
  gmean = mean(betamn0(indmask));
  RescaleFactor = flac0.inorm/gmean;
else
  RescaleFactor = 1;
end
fprintf('RescaleFactor = %g\n',RescaleFactor);
baseline = mri;
baseline.vol = fast_mat2vol(betamn0,mri.volsize);

fname = sprintf('%s/ar1.%s',outanadir,ext);
MRIwrite(ar1,fname);
ar1mn = mri;
ar1mn.vol = mean(ar1.vol,4);
fname = sprintf('%s/ar1mn.%s',outanadir,ext);
MRIwrite(ar1mn,fname);

% Segment based on autocorrelation AR1
nacfsegs = 3;
[edge bincenter binmap] = fast_histeq(ar1mn.vol(indmask),nacfsegs);
acfseg = mri;
acfseg.vol = zeros(acfseg.volsize);
acfseg.vol(indmask) = binmap;
fname = sprintf('%s/acfseg.%s',outanadir,ext);
MRIwrite(acfseg,fname);

% Compute average ar1 in each seg and corresponding acf
fprintf('Computing whitening matrices\n');
tic;
nn = [1:ntptot]';
clear ar1segmn acfsegmn S Sinv W;
for nthseg = 1:nacfsegs
  indseg = find(acfseg.vol==nthseg);
  ar1segmn(nthseg) = mean(ar1mn.vol(indseg));
  %ar1segmn(nthseg) = 0; % No whitening
  fprintf('  seg  %2d  ar1 = %g  (%g)\n',nthseg,ar1segmn(nthseg),toc);
  acfsegmn(:,nthseg) = ar1segmn(nthseg).^(nn-1);
  S(:,:,nthseg) = toeplitz(acfsegmn(:,nthseg));
  Sinv(:,:,nthseg) = inv(S(:,:,nthseg));
  W(:,:,nthseg) = inv(chol(S(:,:,nthseg))');
end

% First pass thru the data to compute beta
fprintf('GLS Beta Pass \n');
tic;
betamat = zeros(nX,nvox);
for nthrun = 1:nruns
  fprintf('  run %d    %g\n',nthrun,toc);
  flac = flac0;
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  yrun = MRIread(flac.funcfspec);
  yrun = RescaleFactor*fast_vol2mat(yrun);  
  indrun = find(tpindrun == nthrun);
  for nthseg = 0:nacfsegs
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
outresdir = sprintf('%s/res',outanadir);
mkdirp(outresdir);
tic;
rsse = 0;
for nthrun = 1:nruns
  fprintf('  run %d    %g\n',nthrun,toc);
  flac = flac0;
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  yrun = MRIread(flac.funcfspec);
  yrun = RescaleFactor*fast_vol2mat(yrun);
  indrun = find(tpindrun == nthrun);
  Xrun = X(indrun,:);
  yhatrun = Xrun*betamat;
  rrun = yrun - yhatrun;
  for nthseg = 1:nacfsegs % ok to skip 0
    indseg = find(acfseg.vol==nthseg);
    acfrun = ar1segmn(nthseg).^([0:flac.ntp-1]');
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

fname = sprintf('%s/mask.%s',outanadir,ext);
MRIwrite(mask,fname);

baseline2 = mri;
baseline2.vol = fast_mat2vol(mean(betamat(ind0,:)),mri.volsize);
fname = sprintf('%s/h-offset.%s',outanadir,ext);
MRIwrite(baseline2,fname);

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

% The contrast matrices were originally computed assuming
% only a single run's worth of nuisance regressors. Recompute.
% And compute the contrasts. Need to handle nuisance contrasts.
fprintf('Computing contrasts  ');
ncontrasts = length(flac0.con);
for nthcon = 1:ncontrasts
  indtask = flac_taskregind(flac0);
  C = flac0.con(nthcon).C;
  Ctaskcon = C(:,indtask);
  [J K] = size(C);
  ZC = zeros(J,nNuis);
  C = [Ctaskcon ZC];
  flac0.con(nthcon).C = C;
  eff = 1/trace(C*inv(X'*X)*C');
  vrf = 1/mean(diag(C*inv(X'*X)*C'));
  fprintf('%2d %-10s   eff = %6.1f   vrf = %6.1f\n',...
	  nthcon,flac0.con(nthcon).name,eff,vrf);
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

  outcondir = sprintf('%s/%s',outanadir,flac0.con(nthcon).name);
  mkdirp(outcondir);

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

if(isempty(analysis)) return; end

% Construct selxavg-style h.dat strucutre for backwards compat
SumXtX = Ctask*X'*X*Ctask';
Nc = hd.Nc;
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

fprintf('Performing GLS GLM Fit  '); 
%[beta rvar vdof] = fast_glmfitw(y,X,acfsegmn,acfseg.vol(:));
% Do it slice by slice to avoid memory problems
betamat = zeros(nX,nvox);
rvarmat = zeros(1,nvox);
rmat    = zeros(ntptot,nvox);
for s = 1:nslices
  indslice = find(svol == s);
  yslice = y(:,indslice);
  segslice = acfseg.vol(indslice);
  [bs rvs vdof rs] = fast_glmfitw(yslice,X,acfsegmn,segslice);
  betamat(:,indslice) = bs;
  rvarmat(:,indslice) = rvs;
  rmat(:,indslice)    = rs;
  clear bs rvs rs yslice;
  pack;
end
fprintf(' ... done %g\n',toc);
clear y;



% Compute smoothess here? Or save r and compute smoothess later?
r = mri;
r.vol = fast_mat2vol(rmat,r.volsize);
clear rmat;
fname = sprintf('%s/res.%s',outanadir,ext);
MRIwrite(r,fname);
clear r;

