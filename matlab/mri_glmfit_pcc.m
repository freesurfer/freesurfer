function r = mri_glmfit_pcc(glmdir)
% r = mri_glmfit_pcc(glmdir)
%
% Computes the partial correlation coefficient (PCC) for each contrast
% in the output dir of mri_glmfit. It also computes a z-score map. The
% z is "two-sided" meaning that if it is converted back to a p, it
% must be done with a two-sized computation to get the original
% p-value. This matches how selxavg3 creates z-scores (and I think FSL
% too).  The output will be two files called pcc.ext and z.ext where
% ext is the file format extension found in the glmdir folder.
%

%
% Original Author: Douglas Greve
%
% Copyright © 2021
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
if(nargin < 1 | nargin > 1)
  fprintf('r = mri_glmfit_pcc(glmdir)');
  return;
end

r = 0;

% Load design matrix
Xfile = sprintf('%s/Xg.dat',glmdir);
X = load(Xfile);

fname = sprintf('%s/beta',glmdir);
beta = MRIread(fname);
betamat = fast_vol2mat(beta);
[fspec fstem fmt] = MRIfspec(fname);

fname = sprintf('%s/rvar',glmdir);
rvar = MRIread(fname);
rvarmat = fast_vol2mat(rvar);

fname = sprintf('%s/mask',glmdir);
mask = MRIread(fname);
indmask = find(mask.vol);
nmask = length(indmask);

% Get a list of contrasts
flist = dir(glmdir);
conlist = '';
for n = 1:length(flist)
  if(~flist(n).isdir) continue; end
  cdat = sprintf('%s/%s/C.dat',glmdir,flist(n).name);
  if(~fast_fileexists(cdat)) continue; end
  conlist = strvcat(conlist,flist(n).name);
end

% Compute partial correlation coefficient (PCC) and Z
ncon = size(conlist,1);
for nthcon = 1:ncon
  conname = deblank(conlist(nthcon,:));
  cdat = sprintf('%s/%s/C.dat',glmdir,conname);  
  C = load(cdat);
  rhomat = fast_glm_pcc(betamat(:,indmask),X,C,rvarmat(indmask));
  rho = beta;
  rho.vol = zeros(size(rho.vol(:,:,:,1)));
  rho.vol(indmask) = rhomat;
  fname = sprintf('%s/%s/pcc.%s',glmdir,conname,fmt);  
  MRIwrite(rho,fname);

  fname = sprintf('%s/%s/sig',glmdir,conname);  
  sig = MRIread(fname);
  p = sign(sig.vol(indmask)).*(10.^(-abs(sig.vol(indmask))));
  z = sig;
  z.vol = zeros(size(z.vol));
  % div 2 converts to a one-sided, but it keeps the sign
  % So, if p=.02, it becomes .01, the z is 2.326
  % If the z is converted back to a p, it must be done
  % with a two-sized computation to get the original p-value
  z.vol(indmask) = fast_p2z(p/2); 
  fname = sprintf('%s/%s/z.%s',glmdir,conname,fmt);  
  MRIwrite(z,fname);
end

return;
