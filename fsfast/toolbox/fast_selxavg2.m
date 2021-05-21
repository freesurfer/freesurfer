function r = fast_selxavg2(varargin)
% r = fast_selxavg2(varargin)
%
% For compatibility with version 1:
%  DOF ignores tpexcl
%  baseline is that of the first run
%
% Incompatibilities:
%   1. Error bars on FIR plots wont be exactly the same when using
%      skip or exluding time points 
%
% Still need to test with:
%   FIR
%   spm hrf
%   gamma with different exp
%   nyq
%   no ext reg
%   time offset?
%   TER
%
% Still need to implement
%   save eres, signal
%   slice-timing correction?
%   vox-wise pct signal change (?)
%
% Update
%   downstream components to look for nii
%   saving betas, h.dat. Update downstream.
%
% Add
%   functional connectivity
%   multiple external reg
%   res fwhm
%   new parfile format
%   compute all contrasts?
%
% Redo
%   inorm
%   whitening
%   fourier
%   change "ces" to "gamma"
%   


%
% fast_selxavg2.m
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



tic;
version = 'fast_selxavg2.m @FS_VERSION@';
fprintf(1,'%s\n',version);
r = 1;
outfmt = 'nii';

%% Print useage if there are no arguments %%
if(nargin == 0)
  print_usage;
  return;
end

%% Make sure $FREESURFER_HOME/matlab is in the path
if(exist('MRIread','file') ~= 2)
  fprintf(['WARNING: it does not look like $FREESURFER_HOME/matlab is' ...
	   ' in your matlab path, so I am going to add it.']);
  FSH = getenv('FREESURFER_HOME');
  fshmatlab = sprintf('%s/matlab',FSH);
  fprintf('Adding %s\n',fshmatlab);
  path(path,fshmatlab);
  clear FSH fshmatlab;
  % Should add it here
  return;
end

%% Parse the arguments %%
s = parse_args(varargin);
if(isempty(s)) return; end
s = check_params(s);
if(isempty(s)) return; end

% Directory of the output
outvolpath = fast_dirname(deblank(s.hvol));
fsd = fast_dirname(deblank(outvolpath));

sxa_print_struct(s,1);

TR  = s.TR;
TER = s.TER;
TW  = s.TotWin;
TPS = s.PreStimWin;
RmBaseline = s.RmBaseline;
RmTrend    = s.RmTrend;
QTrendFit  = s.QTrendFit;
RescaleTarget = s.RescaleTarget;
GammaFit = s.GammaFit;
gfDelta = s.gfDelta;
gfTau   = s.gfTau;
gfAlpha = s.gfAlpha;
nskip = s.nSkip;
firstslice = s.firstslice;
nslices = s.nslices;
TimeOffset = s.TimeOffset;
HanRadius = s.HanRad;
AcqOrder = s.AcqOrder;
SynthSeed = s.SynthSeed;

parfilelist = s.parlist;

instemlist = s.invollist;

tpxlist = s.tpxlist;

hstem     = s.hvol;
eresdir   = s.eresdir;
sigestdir = s.sigestdir;
pctstem   = s.pscvol;
fomnibusstem = s.fomnibusvol;
pomnibusstem = s.pomnibusvol;

%-------------------------------------------------%
if(isempty(s.maskid))
  s.maskid = sprintf('%s/masks/brain',fsd);
  fprintf('INFO: no mask specified, using %s\n',s.maskid);
end      
if(~isempty(s.maskid))
  fprintf('INFO: loading mask %s\n',s.maskid);
  tmp = MRIread(s.maskid);
  if(isempty(tmp))
    fprintf('ERROR: could not load %s\n',s.maskid);
    return;
  end
  mask = tmp.vol;
  clear tmp;
  indmask = find(mask);
  nmasktot = length(indmask);
  fprintf('INFO: mask has %d points\n',nmasktot);
else
  mask = [];
end


%-------------------------------------------------%
nruns = size(parfilelist,1);
Nfir = round(TW/TER);

if(SynthSeed < 0) SynthSeed = sum(100*clock); end
fprintf('SynthSeed = %10d\n',SynthSeed);
randn('state',SynthSeed); 
rand('state',SynthSeed); 

%-- Determine Number of Conditions across all runs----------%
% Get list of unqiue condition IDs for all runs %
condlist = [];
for run = 1:nruns
  par = fmri_ldpar(deblank(parfilelist(run,:)));
  if(isempty(par))
    fprintf('ERROR: reading parfile %s\n',deblank(parfilelist(run,:)));
    return;
  end
  condid = par(:,2);
  clrun = unique(par(:,2));
  condlist = unique([clrun; condlist]);
end

% Count number of presentations per condition
% Remove -1 and 0 %
ind = find(condlist ~= -1 & condlist ~= 0);
condlist = condlist(ind);

Nnnc = length(condlist); % excludes null %
Nc = Nnnc + 1;

fprintf(1,'Conditions Found (%d): ',Nnnc);
fprintf(1,'%2d ',condlist);
fprintf(1,'\n');

if(max(abs(diff(condlist))) > 1)
  fprintf('ERROR: the condition identifiers as found the the paradigm files\n');
  fprintf('       do not appear to be contiguous.\n');
  return;
end

% Check for holes in the list of condition numbers %
if(~isempty(find(diff(condlist)~=1)))
  fprintf(2,'fast_selxavg2: missing conditions\n');
  return;
end

% Count the number per condition %
parall = [];
Npercond= zeros(Nc,1);
for run = 1:nruns
  par = fmri_ldpar(deblank(parfilelist(run,:)));
  npc = [0];
  fprintf(1,'Run %2d: ',run);
  for c = condlist'  %'
    nc = length(find(par(:,2)==c)); 
    npc = [npc nc];
    fprintf(1,'%3d ',nc);
  end
  fprintf(1,'\n');
  Npercond = Npercond + npc'; %'
  parall = [parall; par];
end

SubSampRate = round(TR/TER);

% Get basic info from the first run %
instem = deblank(instemlist(1,:));
y0 = MRIread(instem,1);
if(isempty(y0))
  fprintf('ERROR: could not load %s\n',instem);
  return;
end
nrows   = y0.volsize(1);
ncols   = y0.volsize(2);
nslices = y0.volsize(3);
ntrs    = y0.nframes;
nvox = nrows*ncols*nslices;

%------------------------------------------------
% Build the design matrix across all runs 
fprintf('Building design matrix\n');
X = [];
Xfirall = [];
indTPExcl = [];
runframes = zeros(nruns,2);
n1 = 1;
for run = 1:nruns
  fprintf(1,'     Run %d/%d, %g \n',run,nruns,toc);

  instem = deblank(instemlist(run,:));
  tmp = MRIread(instem,1);
  if(isempty(tmp))
    fprintf('ERROR: could not load %s\n',instem);
    return;
  end
  ntrs = tmp.nframes;
  clear tmp;
  n2 = n1 + ntrs - 1;
  runframes(run,1) = n1; % start 
  runframes(run,2) = n2; % stop

  % Time Point Exclusion %
  if(~isempty(tpxlist))
    TPExcludeFile = deblank(tpxlist(run,:));
    if(~strcmp(TPExcludeFile,'noexcl')) 
      indTPExclRun = fast_ldtpexcl(TPExcludeFile,TR,ntrs,nskip);
      indTPExcl = [indTPExcl; (indTPExclRun+(n1-1))];
      if(s.debug)
	fprintf(1,'       Excluding %d Points: ',ntpx);
	fprintf(1,'%d ',indTPExcl);
	fprintf(1,'\n');
      end
    end
  end
  
  Xdrift = [];
  if(s.PFOrder < 0)
    % Create Baseline/Trend Components of Convolution Matrix %
    Xbaseline = []; Xtrend = []; Xqtrend  = [];
    if(RmBaseline) Xbaseline = fast_baselinemtx(run,ntrs,nruns); end
    if(RmTrend)    Xtrend    = fast_trendmtx(run,ntrs,nruns); end
    if(QTrendFit)  Xqtrend   = fast_quadtrendmtx(run,ntrs,nruns); end
    Xdrift = [Xbaseline Xtrend Xqtrend];
  else
    Xdrift  = fast_polytrendmtx(run,ntrs,nruns,s.PFOrder);
  end
  
  if(s.nyqreg)
    Xnyq = ones(ntrs,1);
    Xnyq(2:2:end) = -1;
    Xnyq = Xnyq - mean(Xnyq); %Make sure it's demeaned
  else
    Xnyq = [];
  end
  
  if(~isempty(s.extreglist))
    extregstem = deblank(s.extreglist(run,:));
    extreg = fmri_ldbvolume(extregstem);
    if(isempty(extreg))
      fprintf('ERROR: could not load %s\n',extregstem);
      return;
    end
    if(size(extreg,3)~=1) extreg = squeeze(extreg)'; %'
    else                  extreg = squeeze(extreg);
    end
    if(s.nextreg < 0) s.nextreg = size(extreg,2); end
    if(s.nextreg > size(extreg,2))
      fprintf('ERROR: %s does not have enough regressors\n',extregstem);
      return;
    end
    % Remove mean of External Regressor %
    extreg = extreg(:,1:s.nextreg);
    extreg = extreg - repmat(mean(extreg), [ntrs 1]);
    extreg = extreg./repmat(std(extreg), [ntrs 1]);
    if(s.extregorthog)
      extreg = ( eye(ntrs) - Xpar*inv(Xpar'*Xpar)*Xpar') * extreg;
    end
    z = zeros(size(extreg));
    extregrun = [repmat(z,[1 (run-1)]) extreg repmat(z,[1 (nruns-run)])];
  else
    extregrun = [];
  end

  % Load paradigm for this run %
  par = fmri_ldpar(deblank(parfilelist(run,:)));
  if(s.autostimdur)
    par = parboxcar(par,TER,[],ntrs*TR);
    if(isempty(par))
      fprintf('ERROR: applying auto stimulus duration\n');
      return;
    end
  end
  if(~isempty(s.stimdur))
    par = parboxcar(par,TER,s.stimdur);
    if(isempty(par))
      fprintf('ERROR: applying stimulus duration\n');
      return;
    end
  end
      
  % Adjust for Time Offset %
  par(:,1) = par(:,1) + TimeOffset;

  % Convert paradigm to FIR stimulus convolution matrix %
  Xfir = fmri_par2scm(par,Nc,SubSampRate*ntrs,TER,Nfir,TPS);

  % For Sub-TR Estimation %
  if(TR ~= TER)
    Xfirtmp = Xfir;
    nn = [1:SubSampRate:size(Xfirtmp,1)];
    Xfir = Xfirtmp(nn,:);
  end
  Xfirall = [Xfirall; Xfir];
  
  % Tranform for Fitting to Gamma Function or SPM HRF %
  if(GammaFit > 0)
    Xpar = fmri_scm2gcm(Xfir,Nnnc,TER,TPS,gfDelta,gfTau,gfAlpha);
    Navgs_per_cond = length(gfDelta);
  elseif(s.spmhrf > -1)
    tspmhrf = TER*[0:Nfir-1]'-TPS;
    hspmhrf = fast_spmhrf(tspmhrf);
    Aspmhrf = hspmhrf;
    dhspmhrf = hspmhrf;
    for nderiv = 1:s.spmhrf
      % Divide by TER for gradient.
      % Multiply by 2.6 to bring 1st deriv to amp of 1
      dhspmhrf = 2.6*gradient(dhspmhrf)/TER;
      Aspmhrf = [Aspmhrf dhspmhrf];
    end
    A = [];
    for c = 1:Nnnc
      A = fast_blockdiag2(A,Aspmhrf);
    end
    Xpar = Xfir*A;
    Navgs_per_cond = s.spmhrf+1;
  else
    Xpar = Xfir;
    Navgs_per_cond = Nfir;
  end
  Ntask = size(Xpar,2);
  
  % Create final Convolution Matrix for ith run %
  Xi = [Xpar Xdrift extregrun Xnyq];
  X = [X; Xi];
  n1 = n2+1;
end % runs
fprintf(1,'done creating design matrix t=%g\n',toc);
ntrstot = size(X,1);

%------------------------------------------------------
% Add time point exclude regressors
ntpx = length(indTPExcl);
if(ntpx > 0)
  Xexcl = zeros(ntrstot,ntpx);
  indtmp = sub2ind(size(Xexcl),indTPExcl,[1:ntpx]');
  Xexcl(indtmp) = 1;
  X = [X Xexcl];
end
fprintf('Excluding a  total of %d points out of %d\n',ntpx,ntrstot);

Ntot = size(X,2);
Nnuis = Ntot - Ntask;

fspec = sprintf('%s/X.mat',outvolpath);
Xfinal=X;
save(fspec,'Xfinal');

%---------------------------------------------------------
% Test the design matrix condition
c = cond(X);
fprintf('Design matrix conditioned is %g\n',c);
if(c > 10000)
  fprintf('ERROR: matrix is ill-conditioned cond=%g\n',c);
  return;
end
ntrstot = size(X,1);

% Get regressor indices for task
TaskInd = [1:Ntask];
% Get regressor indices for baseline
BaselineRegInd = (Ntask) + [0:nruns-1]*(s.PFOrder+1) + 1;

%------------------------------------------------
% Load the data
if(s.SynthSeed == 0)
  fprintf('Loading data\n');
  y = zeros(ntrstot,nvox);
  for run = 1:nruns
    fprintf(1,'     Run %d/%d, t=%g \n',run,nruns,toc);
    
    instem = deblank(instemlist(run,:));
    tmp = MRIread(instem);
    if(isempty(tmp))
      fprintf('ERROR: could not load %s\n',instem);
      return;
    end
    yrun = fast_vol2mat(tmp.vol);
    clear tmp;
    
    if(RescaleTarget > 0)
      MeanValFile = sprintf('%s.meanval',instem);
      [RescaleFactor MeanVal]=fast_rescalefactor(MeanValFile, RescaleTarget);
      fprintf(1,'       Rescaling Global Mean %g,%g,%g\n',...
	      MeanVal,RescaleTarget,RescaleFactor);
      yrun = RescaleFactor * yrun;
    else
      RescaleFactor = 1;
    end
    
    n1 = runframes(run,1);
    n2 = runframes(run,2);
    y(n1:n2,:) = yrun;
    
  end % runs
  fprintf(1,'done loading data t=%g\n',toc);
else
  fprintf('Synthesizing data seed=%d\n',s.SynthSeed);
  y = randn(ntrstot,nvox);
  if(s.AutoWhiten)
    fprintf(1,'Synthesizing with AR1=.3 t=%g\n',toc);
    racf = .3.^[0:ntrstot-1]';
    Far1 = chol(toeplitz(racf));
    y = Far1*y;
  end
  RescaleFactor = 1;
end

%--------------------------------------------------
% Fit the data
fprintf(1,'Fitting GLM t=%g\n',toc);
[beta, rvar, vdof, res] = fast_glmfitw(y,X);
DOF = vdof+ntpx;
rvar = rvar*(vdof/DOF); % compatible with version1
fprintf(1,'Done fitting t=%g\n',toc);

%--------------------------------------------------
% Compute mean baseline image across all runs
b0 = mean(beta(BaselineRegInd,:));

%--------------------------------------------------
% AR1 computation
fprintf(1,'Computing AR1 t=%g ... \n',toc);
sse = sum(res.^2);
indz = find(sse == 0);
sse(indz) = 10^10;
ar1 = sum(res(1:end-1,:) .* res(2:end,:)) ./ sse;
tmp = y0;
tmp.vol = fast_mat2vol(ar1,y0.volsize);
fspec = sprintf('%s/rar1.nii',outvolpath);
MRIwrite(tmp,fspec);

%--------------------------------------------------
% Segment mean image
b0mn = mean(b0(indmask));
fprintf('Global mean = %g\n',b0mn);
fprintf('Computing segs t=%g\n',toc);
nsegs = 20;
xseg = fast_histeq(b0(indmask),nsegs);
hseg = hist(b0(indmask),xseg);
dxseg = xseg(2)-xseg(1);
segmap = y0;
segmap.vol = zeros(y0.volsize);
nperseg = zeros(nsegs,1);
for nthseg = 1:nsegs
  yminseg = xseg(nthseg);
  ymaxseg = xseg(nthseg+1);
  if(nthseg == 1)     yminseg = -inf; end
  if(nthseg == nsegs) ymaxseg = +inf; end
  indseg = find(yminseg < b0(indmask) & b0(indmask) <= ymaxseg);
  nperseg(nthseg) = length(indseg);
  segmap.vol(indmask(indseg)) = nthseg;
  fprintf('seg %2d  %3d  %8.3f %8.3f\n',nthseg,nperseg(nthseg),yminseg,ymaxseg);
end
fspec = sprintf('%s/acfseg.nii',outvolpath);
MRIwrite(segmap,fspec);

%----------------------------------------------------------
fprintf('Computing NACF   (%6.1f)\n',toc);
R = eye(ntrstot) - X*inv(X'*X)*X';% Residual forming matrix 
nn = [1:ntrstot]';
clear ar1seg nacfseg;
for nthseg = 1:nsegs
  indseg = find(segmap.vol==nthseg);
  rar1seg = mean(ar1(indseg));
  racfkjw = fast_yacf_kjw([1 rar1seg]',R);
  ar1seg(nthseg) = racfkjw(2);
  nacfseg(:,nthseg) = ar1seg(nthseg).^(nn-1);
  fprintf('  nthseg = %2d, rar1 = %6.3f, nar1 = %6.3f\n',...
	  nthseg,rar1seg,ar1seg(nthseg));
end
fspec = sprintf('%s/acfseg.lut',outvolpath);
fp = fopen(fspec,'w');
fprintf(fp,'# Color LUT for Autocorrelation Function (ACF) Segmentation\n');
fprintf(fp,'# Index  Name    R G B   Nvoxels\n');
for nthseg = 0:nsegs
  if(nthseg == 0)
    fprintf(fp,' 0    NonBrain           0   0   0    0\n');
  else
    fprintf(fp,'%2d   Bin%02d.AR%02d   %3d %3d %3d    %5d\n',...
	    nthseg,nthseg,round(100*ar1seg(nthseg)),...
	    round(255*rand),round(255*rand),round(255*rand),...
	    nperseg(nthseg));
  end
end
fclose(fp);

%--------------------------------------------------
% ReFit the data with whitening
fprintf(1,'Refitting GLM with whitening t=%g\n',toc);
[beta, rvar, vdof, res] = fast_glmfitw(y,X,nacfseg,segmap.vol);
DOF = vdof+ntpx;
rvar = rvar*(vdof/DOF); % compatible with version1
fprintf(1,'Done fitting t=%g\n',toc);

%--------------------------------------------------
% Omnibus Contrast
fprintf(1,'Computing omnibus contrast t=%g\n',toc);
Comni = [eye(Ntask) zeros(Ntask,Nnuis)];
[F dof1 dof2 ces] = fast_fratiow(beta,X,rvar,Comni,nacfseg,segmap.vol);
fprintf(1,'Done omnibus t=%g\n',toc);
dof2 = DOF; % compatible with version 1
qsum = sum(ces,1);
p = FTest(dof1, dof2, F);
sig = -sign(qsum).*log10(p);
tmp = y0;
tmp.vol = fast_mat2vol(sig,y0.volsize);
fspec = sprintf('%s/omnibus/fsig.%s',outvolpath,outfmt);
MRIwrite(tmp,fspec);

%--------------------------------------------------
% Save betas
tmp = y0;
tmp.vol = fast_mat2vol(beta,y0.volsize);
fspec = sprintf('%s/beta.nii',outvolpath);
MRIwrite(tmp,fspec);

% Save rvar
tmp = y0;
tmp.vol = fast_mat2vol(rvar,y0.volsize);
fspec = sprintf('%s/rvar.nii',outvolpath);
MRIwrite(tmp,fspec);

% Save rstd
tmp = y0;
tmp.vol = fast_mat2vol(sqrt(rvar),y0.volsize);
fspec = sprintf('%s/rstd.nii',outvolpath);
MRIwrite(tmp,fspec);

% Compute and save mean/baseline image across all runs
%b0 = mean(beta(BaselineRegInd,:));
b0 = beta(BaselineRegInd(1),:); % first run (compat with version 1)
tmp = y0;
tmp.vol = fast_mat2vol(b0,y0.volsize);
fspec = sprintf('%s/h-offset.%s',outvolpath,outfmt);
MRIwrite(tmp,fspec);
fspec = sprintf('%s/baseline.nii',outvolpath);
MRIwrite(tmp,fspec);


% ---------------------------------------------------------
% Save the .dat file %
fname = sprintf('%s.dat',hstem);
SumXtXTmp  = Comni*X'*X*Comni';
hCovMtxTmp = Comni*inv(X'*X)*Comni';
hd = fmri_hdrdatstruct;
hd.TR  = TR;
hd.TER = TER;
hd.TimeWindow = TW;
hd.TPreStim = TPS;
hd.Nc = Nc;
hd.Nh = Navgs_per_cond;
hd.Nnnc = Nnnc;
hd.DOF= DOF;
hd.Npercond= Npercond;
hd.Nruns = nruns;
hd.Ntp = ntrstot;
hd.Nrows = nrows;
hd.Ncols = ncols;
hd.Nskip = nskip;
if(s.PFOrder < 0)
  hd.DTOrder = RmBaseline+RmTrend+QTrendFit;
else
  hd.DTOrder = s.PFOrder + 1;
end
hd.RescaleFactor = RescaleFactor;
hd.HanningRadius = 0.0;
hd.BrainAirSeg = 0;
hd.GammaFit = GammaFit ;
hd.gfDelta  = gfDelta;
hd.gfTau    = gfTau;
if(s.spmhrf > -1) % Hack
  hd.GammaFit = s.spmhrf + 1;
  hd.gfDelta = -1*ones(1,s.spmhrf + 1);
  hd.gfTau   = -1*ones(1,s.spmhrf + 1);
end;
%hd.gfAlpha      = gfAlpha; % This is not saved yet
hd.NullCondId    = 0;
hd.SumXtX        = SumXtXTmp;
hd.nNoiseAC      = 0;
hd.CondIdMap     = [0:Nc-1];
hd.hCovMtx       = hCovMtxTmp;
hd.WhitenFlag    = s.WhitenFlag;
hd.runlist       = getrunlist(s.invollist);
hd.funcstem      = basename(deblank(s.invollist(1,:)));
hd.parname       = s.parname;
if(~isempty(s.extreglist))
  hd.extregstem  = basename(deblank(s.extreglist(1,:)));
  hd.nextreg  = s.nextreg;
  hd.extortho = s.extregorthog;
end
fmri_svdat3(fname,hd);

hsxa = fast_beta2sxa(beta,rvar,Nc-1,Navgs_per_cond,X);
ntmp = Navgs_per_cond * 2 * Nc;
hsxa = reshape(hsxa,[ntmp nvox]);
tmp = y0;
tmp.vol = fast_mat2vol(hsxa,y0.volsize);
fspec = sprintf('%s/h.%s',outvolpath,outfmt);
MRIwrite(tmp,fspec);

%---------------------------------------------------------------
fprintf(1,'fast_selxavg2 done %g\n',toc);
r = 0;

return;
%---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%-----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%


%------------- Print Usage ---------------------%
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_selxavg2\n');
  fprintf(1,'     -i   invol ... \n');
  fprintf(1,'     -p   parfile ... \n');
  fprintf(1,'     -parname  parfile \n');
  fprintf(1,'     -extreg   extregfile \n');
  fprintf(1,'     -nextreg  number of external regressors to use\n');
  fprintf(1,'     -parname  parfile \n');
  fprintf(1,'     -tpx tpxfile ... \n');
  fprintf(1,'     -whtmtx whitening matrix file \n');
  fprintf(1,'     -o   hdrstem \n');
  fprintf(1,'     -psc pscstem \n');
  fprintf(1,'     -fomnibus stem \n');
  fprintf(1,'     -pomnibus stem \n');
  fprintf(1,'     -TR   TR\n');
  fprintf(1,'     -TER  TER\n');
  fprintf(1,'     -timewindow totwin  \n');
  fprintf(1,'     -prewindow  prewin  \n');
  fprintf(1,'     -nobaseline  \n');
  fprintf(1,'     -detrend  \n');
  fprintf(1,'     -qtrendfit  \n');
  fprintf(1,'     -rescale  target \n');
  fprintf(1,'     -nskip  n \n');
  fprintf(1,'     -hanrad radius \n');
  fprintf(1,'     -fwhm   width \n');
  fprintf(1,'     -ipr    inplaneres \n');
  fprintf(1,'     -gammafit delta tau \n');
  fprintf(1,'     -timeoffset t \n');
  fprintf(1,'     -acqorder  <linear or interleaved> \n');
  fprintf(1,'     -firstslice sliceno : 0 \n');
  fprintf(1,'     -nslices    nslices : auto \n');
  fprintf(1,'     -eresdir    dir \n');
  fprintf(1,'     -sigestdir  dir \n');
  fprintf(1,'     -synth      seed \n');
  fprintf(1,'     -cfg        file \n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = sxa_struct
  s.invollist      = '';
  s.parlist        = '';
  s.autostimdur    = 0; % extract stim durations from par
  s.stimdur        = []; % stimulus duration
  s.nruns          = 0;
  s.parname        = '';
  s.extregfile     = '';
  s.extreglist     = '';
  s.nextreg        = -1;
  s.extregorthog   =  0;
  s.tpxlist        = '';
  s.WhtnMtxFile     = '';
  s.AutoWhiten     = 0;
  s.NoAutoWhiten   = 0;
  s.SecondPass     = 0;
  s.TauMaxWhiten   = 0;
  s.NMaxWhiten     = 0;
  s.PctWhiten      = 0;
  s.LPFFlag        = 0;
  s.HPF            = [];
  s.WhitenFlag     = 0;
  s.maskid         = []; % for whitening only
  s.hvol           = '';
  s.betavol        = '';
  s.fomnibusvol    = '';
  s.pomnibusvol    = '';
  s.ErrCovMtxStem  = '';
  s.SaveErrCovMtx  = 0;
  s.pscvol   = '';
  s.TR    = '';
  s.TER    = '';
  s.TotWin      = '';
  s.PreStimWin  = 0;
  s.PostStimWin = '';
  s.SegBrainAir = 1;
  s.RmBaseline = 1;
  s.RmTrend    = 0;
  s.QTrendFit  = 0;
  s.PFOrder    = -1;
  s.RescaleTarget = 0;
  s.nSkip  = 0;
  s.FWHM = 0;
  s.InPlaneRes = 0;
  s.HanRad = 0;
  s.gfDelta = [];
  s.gfTau = [];
  s.gfAlpha = 2;
  s.spmhrf = -1;
  s.TimeOffset = 0;
  s.AcqOrder = '';
  s.SynthSeed = 0;
  s.cfgfile = '';
  s.verbose = 0;
  s.firstslice = 0;
  s.nslices    = -1;
  s.eresdir    = '';
  s.sigestdir  = '';
  s.acfdir  = '';
  s.snrdir  = '';
  s.debug = 0;
  s.loginput = 0;
  s.funcstem = '';
  s.nyqreg = 0; % nyquist regressor
return;

%--------------------------------------------------%
% Parse the arguments from the config file %
function argscfg = parse_cfg(args)
  argscfg = args;
  cfgfile = '';
  nargs = length(args);
  narg = 1;
  while(narg <= nargs)
    flag = deblank(args{narg});
    narg = narg + 1;
    if(strcmp(flag,'-cfg'))
      arg1check(flag,narg,nargs);
      cfgfile = args{narg};
      break;
    end
  end

  if(~isempty(cfgfile))
    fid = fopen(cfgfile,'r');
    if(fid == -1)
      fprintf(2,'ERROR: cannot open %s\n',cfgfile);
      argscfg = []; return;
    end
    [s n] = fscanf(fid,'%s',1);
    while(n ~= 0)
      nargs = nargs + 1;;
      argscfg{nargs} = s;
      [s n] = fscanf(fid,'%s',1);
    end
  end

return

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = sxa_struct;
  inputargs = parse_cfg(varargin{1});
  ninputargs = length(inputargs);

  narg = 1;
  while(narg <= ninputargs)

    flag = deblank(inputargs{narg});
    narg = narg + 1;
    %fprintf(1,'Argument: %s\n',flag);
    if(~isstr(flag))
      flag
      fprintf(1,'ERROR: All Arguments must be a string\n');
      error;
    end

    switch(flag)

      case '-i',
        arg1check(flag,narg,ninputargs);
        s.invollist = strvcat(s.invollist,inputargs{narg});
        narg = narg + 1;

      case '-p',
        arg1check(flag,narg,ninputargs);
        s.parlist = strvcat(s.parlist,inputargs{narg});
        narg = narg + 1;

      case '-extreg',
        arg1check(flag,narg,ninputargs);
        s.extregfile = inputargs{narg};
        narg = narg + 1;

      case '-extregorthog',
        s.extregorthog = 1;

      case '-nextreg',
        arg1check(flag,narg,ninputargs);
        s.nextreg = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-parname'},
        arg1check(flag,narg,ninputargs);
        s.parname = inputargs{narg};
        narg = narg + 1;

      case '-stimdur',
        arg1check(flag,narg,ninputargs);
        s.stimdur = sscanf(inputargs{narg},'%g',1);
        narg = narg + 1;

      case {'-autostimdur'},
        s.autostimdur = 1;

      case {'-tpx','-tpexclfile'}
        arg1check(flag,narg,ninputargs);
        s.tpxlist = strvcat(s.tpxlist,inputargs{narg});
        narg = narg + 1;

	% Arg is max ACF delay
      case {'-autowhiten'} 
        arg1check(flag,narg,ninputargs);
        s.TauMaxWhiten = sscanf(inputargs{narg},'%f',1);
        if(s.TauMaxWhiten > 0) s.AutoWhiten = 1; end
        narg = narg + 1;

      case {'-mask'},
        arg1check(flag,narg,ninputargs);
        s.maskid = inputargs{narg};
        narg = narg + 1;

      case {'-o','-h'},
        arg1check(flag,narg,ninputargs);
        s.hvol = inputargs{narg};
        narg = narg + 1;

      case {'-beta'},
        arg1check(flag,narg,ninputargs);
        s.betavol = inputargs{narg};
        narg = narg + 1;

      case {'-fomnibus'},
        arg1check(flag,narg,ninputargs);
        s.fomnibusvol = inputargs{narg};
        narg = narg + 1;

      case {'-pomnibus'},
        arg1check(flag,narg,ninputargs);
        s.pomnibusvol = inputargs{narg};
        narg = narg + 1;

      case {'-ecovmtx'},
        arg1check(flag,narg,ninputargs);
        s.ErrCovMtxStem = inputargs{narg};
        narg = narg + 1;

      case {'-svecovmtx','-sverrcovmtx','-svecvm','-svacf'},
        s.SaveErrCovMtx = 1;

      case {'-psc','-percent'},
        arg1check(flag,narg,ninputargs);
        s.pscvol = inputargs{narg};
        narg = narg + 1;

      case {'-TR'}
        arg1check(flag,narg,ninputargs);
        s.TR = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-TER'}
        arg1check(flag,narg,ninputargs);
        s.TER = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-timewindow','-totwin','-tw'}
        arg1check(flag,narg,ninputargs);
        s.TotWin = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-prewindow','-prewin','-prestim'}
        arg1check(flag,narg,ninputargs);
        s.PreStimWin = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-postwindow','-postwin','-poststim'}
        arg1check(flag,narg,ninputargs);
        s.PostStimWin = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-timeoffset'}
        arg1check(flag,narg,ninputargs);
        s.TimeOffset = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-basegment'}
        s.SegBrainAir = 1;
  
      case {'-nobasegment'}
        s.SegBrainAir = 0;
  
      case {'-nobaseline'}
        s.RmBaseline = 0;
  
      case {'-baseline'}
        s.RmBaseline = 1;
  
      case {'-detrend'}
        s.RmTrend = 1;
  
      case {'-qtrendfit'}
        s.QTrendFit = 1;
  
      case {'-lpf'}
        s.LPFFlag = 1;
  
      case {'-hpf'}
        arg2check(flag,narg,ninputargs);
        s.HPF = sscanf(inputargs{narg},'%f %f',1);
        narg = narg + 1;

      case {'-polyfit'}
        arg1check(flag,narg,ninputargs);
        s.PFOrder = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-rescale'}
        arg1check(flag,narg,ninputargs);
        s.RescaleTarget = sscanf(inputargs{narg},'%f',1);
	fprintf('RescaleTarget = %g\n',s.RescaleTarget);
        narg = narg + 1;

      case {'-nskip'}
        arg1check(flag,narg,ninputargs);
        s.nSkip = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-hanrad'}
        arg1check(flag,narg,ninputargs);
        s.HanRad = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-fwhm'}
        arg1check(flag,narg,ninputargs);
        s.FWHM = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-ipr'}
        arg1check(flag,narg,ninputargs);
        s.InPlaneRes = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-gammafit'}
        arg2check(flag,narg,ninputargs);
        gfDelta = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
        gfTau   = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
        s.gfDelta = [s.gfDelta gfDelta];
        s.gfTau   = [s.gfTau   gfTau];

      case {'-gammaexp'}
        arg1check(flag,narg,ninputargs);
        s.gfAlpha = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case '-spmhrf', % Argument is number of derivatives
        arg1check(flag,narg,ninputargs);
        s.spmhrf = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-acqorder',
        arg1check(flag,narg,ninputargs);
        s.AcqOrder = inputargs{narg};
        narg = narg + 1;

      case {'-firstslice', '-fs'}
        arg1check(flag,narg,ninputargs);
        s.firstslice = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nslices', '-ns'}
        arg1check(flag,narg,ninputargs);
        s.nslices = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-eresdir',
        arg1check(flag,narg,ninputargs);
        s.eresdir = inputargs{narg};
        narg = narg + 1;

      case '-acfdir',
        arg1check(flag,narg,ninputargs);
        s.acfdir = inputargs{narg};
        narg = narg + 1;

      case '-snrdir',
        arg1check(flag,narg,ninputargs);
        s.snrdir = inputargs{narg};
        narg = narg + 1;

      case {'-sigestdir','-signaldir'}
        arg1check(flag,narg,ninputargs);
        s.sigestdir = inputargs{narg};
        narg = narg + 1;

      case '-cfg',
        % This is actually handled by parse_cfg
        arg1check(flag,narg,ninputargs);
        narg = narg + 1;

      case '-synth', 
        arg1check(flag,narg,ninputargs);
        s.SynthSeed = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-nyqreg',
        s.nyqreg = 1;

      case '-verbose',
        s.verbose = 1;

      case '-log',
        s.loginput = 1;

      case {'-whtnmtx'}
        arg1check(flag,narg,ninputargs);
        s.WhtnMtxFile = strvcat(s.WhtnMtxFile,inputargs{narg});
        narg = narg + 1;
      case {'-whiten'} 
        s.WhitenFlag = 1; 
        % Requires that -whtnmtx be specified for each run
        % the whtn matrix will be stored in matlab4 format
        % in the variable named W.
      case {'-noautowhiten'} % To ease recursive calls
        s.NoAutoWhiten = 1;
        s.AutoWhiten = 0;
        s.TauMaxWhiten = 0;
        s.SecondPass = 1;

      % ignore these guys %
     case {'-monly', '-nullcondid','-umask','-sveres',...
	   '-svsignal','-svsnr','-sxamver'},
        arg1check(flag,narg,ninputargs);
        narg = narg + 1;

      case {'-debug','-echo'}, % ignore
        s.debug = 1;

      otherwise
        fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
        s = [];
        return;

    end % --- switch(flag) ----- %

  end % while(narg <= ninputargs)

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check that there is at least one more argument %%
function arg1check(flag,nflag,nmax)
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    error;
  end
return;
%--------------------------------------------------%
%% Check that there are at least two more arguments %%
function arg2check(flag,nflag,nmax)
  if(nflag > nmax-1 ) 
    fprintf(1,'ERROR: Flag %s needs two arguments',flag);
    error;
  end
return;


%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');

  s.nruns = size(s.invollist,1);
  npars = size(s.parlist,1);
  ntpxs = size(s.tpxlist,1);

  if(s.nruns < 1) 
    fprintf(2,'ERROR: No input volumes specified\n');
    s=[]; return;
  end

  if(s.nslices < 0)
    instem = deblank(s.invollist(1,:));
    [s.nslices nrows ncols ntrs] = fmri_bvoldim(instem);
    if(s.nslices == 0) 
      fprintf(2,'ERROR: Volume %s does not exist\n',instem);
      s=[]; return;
    end      
  end

  if(npars ~= 0 & ~isempty(s.parname) ) 
    fprintf(2,'ERROR: Cannot specify both -p and -parname\n');
    s=[]; return;
  end

  if(npars == 0 & isempty(s.parname) ) 
    fprintf(2,'ERROR: No paradigm specified\n');
    s=[]; return;
  end

  if( ~isempty(s.parname) ) 
    for n = 1:s.nruns
      involpath = fast_dirname(deblank(s.invollist(n,:)));
      par = sprintf('%s/%s',involpath,s.parname);
      s.parlist = strvcat(s.parlist,par);
    end
    npars = size(s.parlist,1);
  end

  if(npars ~= s.nruns)
    fprintf(2,'ERROR: Number of input volumes (%d) and paradigms (%d) differ\n',...
                  s.nruns,npars);
    s=[]; return;
  end

  if(ntpxs ~= 0 & ntpxs ~= s.nruns)
    fprintf(2,'ERROR: Number of input volumes (%d) and tpexcl files (%d) differ\n',...
                  s.nruns,ntpxs);
    s=[]; return;
  end

  if(~isempty(s.extregfile) ) 
    for n = 1:s.nruns
      involpath = fast_dirname(deblank(s.invollist(n,:)));
      extregtmp = sprintf('%s/%s',involpath,s.extregfile);
      s.extreglist = strvcat(s.extreglist,extregtmp);
    end
  end

  if(size(s.hvol,1) ~= 1)
    fprintf(2,'ERROR: No output volume specified\n');
    s = []; return;
  end

  if(length(s.TR) == 0)
    fprintf(2,'ERROR: No TR specified\n');
    s = []; return;
  end

  if(length(s.TotWin) == 0)
    fprintf(2,'ERROR: No Time Window specified \n');
    s = []; return;
    %fprintf(2,'INFO: No Time Window specified ...\n');
    %fprintf(2,' Setting to 20 sec\n');
    %s.TotWin = 20;
  end

  if(length(s.AcqOrder) > 0)
    if(~strcmpi(s.AcqOrder,'Linear') & ~strcmpi(s.AcqOrder,'Interleaved'))
     fprintf(2,'ERROR: Acquisition Order %s unknown (Linear or Interleaved)\n',...
              s.AcqOrder);
      s = []; return;
    end
  end

  if(length(s.TER) == 0) s.TER = s.TR; end

  if(s.firstslice < 0) 
    fprintf('ERROR: firstslice (%d) < 0',s.firstslice);
    s = []; return;
  end

  s.GammaFit = length(s.gfDelta);

  if(s.SaveErrCovMtx)
    s.ErrCovMtxStem = sprintf('%s-ecvm',s.hvol);
  end

  if(s.FWHM > 0 & s.HanRad > 0)
    fprintf('ERROR: Cannot specify both -hanrad and -fwhm\n');
    s = []; return;
  end

  if(s.FWHM > 0 & isempty(s.InPlaneRes))
    fprintf('ERROR: Need -ipr with -fwhm\n');
    s = []; return;
  end

  if(s.FWHM > 0 )
    s.HanRad = pi*s.FWHM/(2*s.InPlaneRes*acos(.5));
  end

  fprintf('AutoStimDur: %d\n',s.autostimdur);
  fprintf('StimDur: %g\n',s.stimdur);
  
  if(s.autostimdur & ~isempty(s.stimdur))
    fprintf('ERROR: cannot specify both autostimdur and stimdur\n');
    s = []; return;
  end
    
return;

%--------------------------------------------------%
%% Print data structure
function s = sxa_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'Number of Runs: %d\n',s.nruns);

  fprintf(fid,'Input Volume List\n');
  for n = 1:size(s.invollist,1),
    fprintf(fid,'  %d  %s\n',n,s.invollist(n,:));    
  end

  fprintf(fid,'Input Pardigm File List\n');
  for n = 1:size(s.parlist,1),
    fprintf(fid,'  %d  %s\n',n,s.parlist(n,:));    
  end

  if(~isempty(s.tpxlist))
    fprintf(fid,'TP Exclude File List\n');
    for n = 1:size(s.tpxlist,1),
      fprintf(fid,'  %d  %s\n',n,s.tpxlist(n,:));    
    end
  end

  fprintf(fid,'Output Volume  %s\n',s.hvol);
  if(~isempty(s.betavol))
     fprintf(fid,'Beta Volume  %s\n',s.betavol);
  end
  if(~isempty(s.fomnibusvol))
    fprintf(fid,'F Omnibus Volume  %s\n',s.fomnibusvol);
  end
  if(~isempty(s.pomnibusvol))
    fprintf(fid,'Sig Omnibus Volume  %s\n',s.pomnibusvol);
  end
  fprintf(fid,'TR    %f\n',s.TR);
  fprintf(fid,'TER   %f\n',s.TER);
  fprintf(fid,'Total   Window  %g\n',s.TotWin);
  fprintf(fid,'PreStim Window  %g\n',s.PreStimWin);
  fprintf(fid,'Remove Baseline %d\n',s.RmBaseline);
  fprintf(fid,'Remove Trend    %d\n',s.RmTrend);
  fprintf(fid,'Remove QTrend   %d\n',s.QTrendFit);
  fprintf(fid,'Rescale Target  %g\n',s.RescaleTarget);
  fprintf(fid,'nSkip           %d\n',s.nSkip);
  fprintf(fid,'InPlane Res     %g\n',s.InPlaneRes);
  fprintf(fid,'FWHM            %g\n',s.FWHM);
  fprintf(fid,'Hanning Radius  %g\n',s.HanRad);
  fprintf(fid,'Time Offset     %g\n',s.TimeOffset);
  if(~isempty(s.AcqOrder))
    fprintf(fid,'Acquistion Order %s\n',s.AcqOrder);
  end

  fprintf(fid,'GammaFit        %d\n',s.GammaFit);
  for n = 1:s.GammaFit
    fprintf(fid,'%d  %g  %g\n',n,s.gfDelta,s.gfTau);
  end
  fprintf('GammaFit Alpha: %g\n',s.gfAlpha);
  fprintf('SPM HRF: %g\n',s.spmhrf);

  fprintf(fid,'Seg Brain/Air   %d\n',s.SegBrainAir);
  fprintf(fid,'SynthSeed       %d\n',s.SynthSeed);

  if(~isempty(s.ErrCovMtxStem))
    fprintf(fid,'ErrCovMtx Stem   %s\n',s.ErrCovMtxStem);
  end

  if(~isempty(s.WhtnMtxFile))
    fprintf(fid,'WhtnMtx File   %s\n',s.WhtnMtxFile);
  end

  if(~isempty(s.extregfile))
    fprintf(fid,'ExtReg File   %s\n',s.extregfile);
    fprintf(fid,'NExtReg       %d\n',s.nextreg);
    fprintf(fid,'ExtRegOrthog  %d\n',s.extregorthog);
  end

  fprintf(fid,'firstslice   %d\n',s.firstslice);
  fprintf(fid,'nslices      %d\n',s.nslices);
  fprintf(fid,'nyqreg       %d\n',s.nyqreg);

return;
%--------------------------------------------------%
function runlist = getrunlist(invollist)
  nruns = size(invollist,1);
  runlist = [];
  for run = 1:nruns
    invol = deblank(invollist(run,:));
    tmp = fast_dirname(invol);
    runid = basename(tmp);
    runno = sscanf(runid,'%d');
    runlist = [runlist runno];
  end
return;

