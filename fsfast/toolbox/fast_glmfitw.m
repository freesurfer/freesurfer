function [beta, rvar, vdof, r] = fast_glmfitw(y,X,nacf,nacfmap)
% [beta, rvar, vdof, r] = fast_glmfitw(y,X,<nacf>,<nacfmap>)
%
% Fits a GLM, possibly including prewhitening
%
% nacf is the noise autocorrelation function. Can have multiple
%   columns. If so, need nacfmap.
% nacfmap is a map of which column of nacf is to be applied
%   to which voxels in y. Voxels that have nacfmap=0 will
%   be treated as white.
%
% This gives the same answer as fast_glmfit when no whiting.
%
% See also: fast_fratiow, FTest, fast_glmfit, fast_fratio.
%
%


%
% fast_glmfitw.m
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

if(nargin < 2 | nargin > 4)
  fprintf('[beta, rvar, vdof, r] = fast_glmfitw(y,X,<nacf>,<nacfmap>)\n');
  return;
end

[nf nv] = size(y);
if(size(X,1) ~= nf)
  fprintf('ERROR: X and y have different number of frames\n');
  return;
end

if(~exist('nacf','var'))      nacf = []; end
if(~exist('nacfmap','var'))   nacfmap = []; end

if(~isempty(nacf))
  if(size(nacf,2) > 1 & isempty(nacfmap))
    fprintf('ERROR: must supply nacfmap with multiple nacf\n');
    return;
  end
  if(size(nacf,1) ~= nf)
    fprintf('ERROR: nframes in nacf must equal that in y and X\n');
    return;
  end
end

nbeta = size(X,2);
vdof = nf-nbeta;

if(isempty(nacf))
  beta = (inv(X'*X)*X')*y;
  r = y - X*beta;
  rvar = sum(r.^2)/vdof;
  return;
end

usematrix = 0;
%if(usematrix) fprintf('Using matrix filtering\n');
%else          fprintf('Using fft filtering\n');
%end

% Only gets here if nacf is non-empty
if(isempty(nacfmap)) nacfmap = ones(nv,1); end
nbins = size(nacf,2);
beta = zeros(nbeta,nv);
if(nargout == 4) r = zeros(nf,nv); end
rvar = zeros(1,nv);
rsse = zeros(1,nv);
X_fft = fft(X,2*nf);
for nthbin = 0:nbins
  indbin = find(nacfmap==nthbin);
  nbin = length(indbin);
  %fprintf(' nthbin = %d, nbin = %d \n',nthbin,nbin);
  if(isempty(indbin)) continue; end
  if(nthbin ~= 0)
    nacfbin = nacf(:,nthbin);
    if(usematrix)
      % Matrix formulation
      W = chol(inv(toeplitz(nacfbin)));
      ybin = W*y(:,indbin);
      Xbin = W*X;
    else
      % FFT formulation
      ybin = y(:,indbin);
      nacfbin_fft = fft(nacfbin,2*nf);
      ybin_fft = fft(ybin,2*nf)./repmat(conj(nacfbin_fft),[1 nbin]);
      ybin = real(ifft(ybin_fft));
      ybin = ybin(1:nf,:);
      clear ybin_fft;
      Xbin_fft = X_fft./repmat(conj(nacfbin_fft),[1 nbeta]);
      Xbin = real(ifft(Xbin_fft));
      Xbin = Xbin(1:nf,:);
      clear Xbin_fft;
    end
    beta(:,indbin) = (inv(Xbin'*Xbin)*Xbin')*ybin;
    rbin = ybin - Xbin*beta(:,indbin);
  else
    Xbin = X;
    beta(:,indbin) = (inv(Xbin'*Xbin)*Xbin')*y(:,indbin);
    rbin = y(:,indbin) - Xbin*beta(:,indbin);
  end
  if(nargout == 4) r(:,indbin) = rbin; end
  rsse(indbin) = sum(rbin.^2);
  clear rbin;
end

rvar = rsse/vdof;

return;  
  
