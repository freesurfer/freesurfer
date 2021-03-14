function ok = flacrfx(flac,contrast)
% ok = flacrfx(flac,<contrast>)
%
% Random effects averaging of the first-level analysis. Uses
% weighted least squares. Algorithm is correct even when 
% contrast matrix has multiple rows.
%
% flac is a customized flac.
% contrast is name of a contrast, or empty for all.
%
% Saves results to flac/rfx.
%
%
%


%
% flacrfx.m
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

ok = 0;
if(nargin < 1 | nargin > 2)
  fprintf('ok = flacffx(flac,<contrast>)\n');
  return;
end

flac.nthrun = 1;
flac = flac_customize(flac);
if(isempty(flac)) return; end

nruns = size(flac.runlist,1);
if(nruns == 1)
  fprintf('ERROR: cannot RFx with only 1 run\n');
  return;
end

ncon = length(flac.con);
if(~exist('contrast','var')) contrast = []; end
if(isempty(contrast))
  nthconlist = [1:ncon];
else
  nthconlist = flac_conindex(contrast,flac);
  if(isempty(nthconlist)) 
    fprintf('ERROR: contrast %s does not exist\n',contrast);
    return; 
  end
end

flarfxdir = sprintf('%s/%s/%s/rfx',flac.sess,flac.fsd, ...
		    flac.name);
mkdirpcmd = sprintf('mkdir -p %s',flarfxdir);
unix(mkdirpcmd);

matfile = sprintf('%s/%s/%s/%s/flac.mat',flac.sess,flac.fsd, ...
		    flac.name,flac.runlist(1,:));
flac.mat = load(matfile);
nseg = size(flac.mat.nacfseg,2);
Nv = prod(flac.mri.volsize);

% Load beta and rvar for all runs at the start
betaruns = [];
rvarruns = [];
for nthrun = 1:nruns
  flac.nthrun = nthrun;
  flac = flac_customize(flac);

  betarun = MRIread(flac.betafspec);
  if(isempty(betarun)) return; end
  betarun.vol = fast_vol2mat(betarun.vol);
  betaruns(:,:,nthrun) = betarun.vol; %

  rvarrun = MRIread(flac.rvarfspec);
  rvarrun.vol = fast_vol2mat(rvarrun.vol);
  rvarruns(1,:,nthrun) = rvarrun.vol; %
end
betamn = mean(betaruns,3);

outfspec = sprintf('%s/beta.mgh',flarfxdir);
mri = flac.mri;
mri.vol = fast_mat2vol(betamn,mri.volsize);
MRIwrite(mri,outfspec);

% Now go through each contrast -----------------------
for nthcon = nthconlist
  
  C = flac.con(nthcon).C;
  J = size(C,1);
  fprintf('  contrast %d %s, J=%d\n',nthcon,flac.con(nthcon).name,J);
  condir = sprintf('%s/%s',flarfxdir,flac.con(nthcon).name);
  mkdirpcmd = sprintf('mkdir -p %s',condir);
  unix(mkdirpcmd);

  beta_rfx = zeros(J,Nv);
  rvar_rfx = zeros(1,Nv);
  F_rfx    = zeros(1,Nv);
  Fsig_rfx = zeros(1,Nv);
  for nthseg = 0:nseg
    indseg = find(flac.mat.acfseg.vol == nthseg);
    
    y_rfx = [];
    M_rfx = [];
    X_rfx = [];
    for nthrun = 1:nruns
      flac.nthrun = nthrun;
      flac = flac_customize(flac);
      
      gamrun = C*betaruns(:,indseg,nthrun);
      y_rfx = [y_rfx; gamrun];

      rvarrunsegmn = mean(rvarruns(1,indseg,nthrun));
      
      % Compute the covariance for gamma
      Xrun = flac.X;
      if(nthseg > 0)
	nacf = flac.mat.nacfseg(:,nthseg);
	Srun = toeplitz(nacf);
	Mrun = rvarrunsegmn*inv(C*Xrun'*inv(Srun)*Xrun*C');
      else
	Mrun = rvarrunsegmn*inv(C*Xrun'*Xrun*C');
      end
      M_rfx = fast_blockdiag2(M_rfx,Mrun);
      
      X_rfx = [X_rfx; eye(J)];
    end % run

    % Construct the GLM for the RFx analysis
    W_rfx  = inv(chol(M_rfx)');
    y_rfxw = W_rfx*y_rfx;
    X_rfxw = W_rfx*X_rfx;
    C_rfx  = eye(size(X_rfxw,2));

    [betatmp rvartmp vdof] = fast_glmfitw(y_rfxw,X_rfxw);
    [Ftmp dof1 dof2 ces] = fast_fratiow(betatmp,X_rfxw,rvartmp,C_rfx);
    Fsigtmp = FTest(dof1, dof2, Ftmp);

    beta_rfx(:,indseg) = betatmp;
    rvar_rfx(:,indseg) = rvartmp;
    F_rfx(:,indseg)    = Ftmp;
    Fsig_rfx(:,indseg) = Fsigtmp;

  end % seg
  
  outfspec = sprintf('%s/gam.mgh',condir);
  mri = flac.mri;
  mri.vol = fast_mat2vol(beta_rfx,mri.volsize);
  MRIwrite(mri,outfspec);
  
  outfspec = sprintf('%s/gamvar.mgh',condir);
  mri = flac.mri;
  mri.vol = fast_mat2vol(rvar_rfx,mri.volsize);
  MRIwrite(mri,outfspec);
  
  outfspec = sprintf('%s/f.mgh',condir);
  mri = flac.mri;
  mri.vol = fast_mat2vol(F_rfx,mri.volsize);
  MRIwrite(mri,outfspec);
  
  outfspec = sprintf('%s/fsig.mgh',condir);
  mri = flac.mri;
  mri.vol = fast_mat2vol(Fsig_rfx,mri.volsize);
  indnz = find(mri.vol ~= 0);
  mri.vol(indnz) = -log10(mri.vol(indnz));
  MRIwrite(mri,outfspec);
  
end % contrast

ok = 1;
return;

