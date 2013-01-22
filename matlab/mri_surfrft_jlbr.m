function r = mri_surfrft_jlbr(yfile,glmdir,vwthresh,sgn,subject,hemi)
% r = mri_surfrft_jlbr(yfile,glmdir,vwthresh,<sgn>,<subject>,<hemi>)
%
% Performs GLM analysis and RFT correction on the output of mri_glmfit
% when the analysis has been performed on the surface. The full
% analysis is recomputed.
%
% This is based on Jorge Louis Bernal-Rusiel's RFT program (search
% for RFT_FDR in www.nitrc.org). Both tools use Keith Worsley's
% SurfStat program www.math.mcgill.ca/keith/surfstat.
%
% yfile -- input data. This is the file passed to mri_glmfit with
%   the --y option.
% glmdir -- GLM directory (argument to the --glmdir option). 
% vwthresh -- voxel-wise p-value threshold (ie, the cluster-forming threshold)
% sgn -- sign of the contrast (+1 or -1, default is +1)
% subject - subject whose surface the analysis is done on. This
%   will be read from the GLMDIR by default.
% hemi - hemisphere the analysis is done on. This  will be read
%   from the GLMDIR by default.
%
% The analysis will be run for each contrast in the GLM dir. The
% output will be sig.cw.pos.rft.mgh (for cluster-corrected signficance
% maps with pos contrast) and sig.vw.pos.rft.mgh (for voxel-corrected
% signficance maps with pos contrast). 
%
% Limitations: does not work with per-voxel regressors, weighted-least
% squares, fixed effects, or multi-variate contrasts.
%
% $Id: mri_surfrft_jlbr.m,v 1.2.2.1 2013/01/22 20:59:09 nicks Exp $

%
% Original Author: Jorge Louis Bernal-Rusiel and Douglas Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/01/22 20:59:09 $
%    $Revision: 1.2.2.1 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

r = 1;
if(nargin < 3 | nargin > 6)
  fprintf('r = mri_surfrft_jlbr(yfile,glmdir,vwthresh,<sgn>,<subject>,<hemi>)\n');
  return;
end

r = 0;

SUBJECTS_DIR = getenv('SUBJECTS_DIR');
surfname = 'white';

if(~exist('sgn','var')) sgn = +1; end
if(sgn == +1) sgnstring = 'pos'; end
if(sgn == -1) sgnstring = 'neg'; end

if(~exist('subject','var'))
  % Get subject and hemi from 'surface' file. New with FS 5.1
  fname = sprintf('%s/surface',glmdir);
  fp = fopen(fname,'r');
  subject = fscanf(fp,'%s',1);
  hemi = fscanf(fp,'%s',1);
  fclose(fp);
  %[subject hemi] = textread(fname,'%s %s');
  %subject = char(subject);
  %hemi = char(hemi);
end

% Load the surface
spath = sprintf('%s/%s/surf/%s.%s',SUBJECTS_DIR,subject,hemi,surfname);
surf = SurfStatReadSurf(spath,'b',2);

% Load raw data
ymri = MRIread(yfile);
y = fast_vol2mat(ymri);
mri = ymri;

% Load mask
fname = sprintf('%s/mask',glmdir);
maskmri = MRIread(fname);
mask = logical(fast_vol2mat(maskmri));

% Load design matrix
Xfile = sprintf('%s/Xg.dat',glmdir);
X = load(Xfile);

% Solve the linear model
tX = term(X,'X');
slm = SurfStatLinMod(y, tX, surf,1,.01,0.1);
slm.k = 1;

% pairs of connected vtxnos (does not excl mask so nans are
% present)
if(0)
edg=SurfStatEdg(surf); 
ntp = size(X,1);
nX  = size(X,2);
R = eye(ntp) - X*inv(X'*X)*X';
yr = R*y;
ysse = sum(yr.^2); % residual
yrn = yr./repmat(sqrt(ysse),[ntp 1]); % voxel-norm resid
d = yrn(:,edg(:,1))-yrn(:,edg(:,2)); % norm resid diff bet pairs
resl = sum(d.^2); % resel size
end

% Get a list of contrasts
flist = dir(glmdir);
conlist = '';
for n = 1:length(flist)
  if(~flist(n).isdir) continue; end
  cdat = sprintf('%s/%s/C.dat',glmdir,flist(n).name);
  if(~fast_fileexists(cdat)) continue; end
  conlist = strvcat(conlist,flist(n).name);
end

% Get p-values for each contrast, save output
ncon = size(conlist,1);
for nthcon = 1:ncon
  conname = deblank(conlist(nthcon,:));
  cdat = sprintf('%s/%s/C.dat',glmdir,conname);  
  C = sgn*load(cdat);
  slmC = SurfStatT(slm, C);

  [pval peak clus] = SurfStatP(slmC, mask, vwthresh);

  mri.vol = fast_mat2vol(-log10(pval.C)*sgn,mri.volsize);
  fname = sprintf('%s/%s/sig.cw.%s.mgh',glmdir,conname,sgnstring);  
  MRIwrite(mri,fname);
  mri.vol = fast_mat2vol(-log10(pval.P)*sgn,mri.volsize);
  fname = sprintf('%s/%s/sig.vw.%s.mgh',glmdir,conname,sgnstring);  
  MRIwrite(mri,fname);

  fname = sprintf('%s/%s/sig.cw.%s.dat',glmdir,conname,sgnstring);
  fp = fopen(fname,'w');
  fprintf(fp,'# vwthresh %g\n',vwthresh);
  nclusters = length(clus.P);
  for nthcluster = 1:nclusters
    fprintf(fp,'%2d %7.6f %5d\n',nthcluster,clus.P(nthcluster),clus.nverts(nthcluster));
  end
  fclose(fp);
  
  % This gives same t as KJW
  %[beta rvar] = fast_glmfit(y,X);
  %[F, Fsig, ces, cesvar] = fast_fratio(beta,X,rvar,C);
  %t = sqrt(F).*sign(ces);
end


return;
