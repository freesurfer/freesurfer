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
% $Id: fast_glmfitw.m,v 1.2 2004/11/14 22:33:26 greve Exp $

if(nargin < 2 | nargin > 4)
  fprintf('[beta, rvar, vdof, r] = fast_glmfitw(y,X,<nacf>,<nacfmap>)\n');
  return;
end

[nf nv] = size(y);
if(size(X,1) ~= nf)
  fprintf('ERROR: X and y have different number of frames\n');
  return;
end

if(~exist('nacf','var'))    nacf = []; end
if(~exist('nacfmap','var')) nacfmap = []; end

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

% Only gets here if nacf is non-empty
if(isempty(nacfmap)) nacfmap = ones(nv,1); end
nbins = size(nacf,2);
beta = zeros(nbeta,nv);
r    = zeros(nf,nv);
rvar = zeros(1,nv);
for nthbin = 0:nbins
  indbin = find(nacfmap==nthbin);
  if(isempty(indbin)) continue; end
  if(nthbin ~= 0)
    nacfbin = nacf(:,nthbin);
    %W = inv(chol(toeplitz(nacfbin))');
    W = chol(inv(toeplitz(nacfbin)));
    ybin = W*y(:,indbin);
    Xbin = W*X;
  else
    ybin = y(:,indbin);
    Xbin = X;
  end
  beta(:,indbin) = (inv(Xbin'*Xbin)*Xbin')*ybin;
  r(:,indbin) = ybin - Xbin*beta(:,indbin);
end

rvar = sum(r.^2)/vdof;

return;  
  
