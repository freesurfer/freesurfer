function [F, dof1, dof2, ces, cescvm] = fast_fratiow(beta,X,rvar,C,nacf,nacfmap)
% [F dof1 dof2 ces cescvm] = fast_fratiow(beta,X,rvar,C,<nacf>,<nacfmap>)
%
% beta - GLM regression coefficients from fast_glmfitw
% rvar - residual error variance from GLM from fast_glmfitw
% X - design matrix used with fast_glmfit.
% nacf is the noise autocorrelation function used with fast_glmfit.
%  Can have multiple columns. If so, need nacfmap.
% nacfmap is a map of which column of nacf is to be applied
%   to which voxels in y as used with fast_glmfitw. Voxels that 
%   have nacfmap=0 will be treated as white.
% 
% F - F-ratio. The Fsig is not computed here. Use:
%   Fsig = FTest(dof1, dof2, F, dof2max);
% ces - contrast effect size = C*beta.
% cescvm - contrast effect size covariance matrix = rvar * inv(C*inv(X'*W*X)*C');
%
% See also: fast_glmfitw, FTest, fast_glmfit, fast_fratio.
%
% $Id: fast_fratiow.m,v 1.4 2005/06/19 05:22:53 greve Exp $

if(nargin < 4 | nargin > 6)
  fprintf('[F dof1 dof2 ces cescvm] = fast_fratiow(beta,X,rvar,C,<nacf>,<nacfmap>)\n');
  return;
end

[nf nbeta] = size(X);
nv = size(beta,2);

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

dof1 = size(C,1); % J %
dof2 = nf-nbeta;

% Handle case where some voxels are zero
indz  = find(rvar == 0);
indnz = find(rvar ~= 0);
if(~isempty(indz))
  beta = beta(:,indnz);
  rvar = rvar(:,indnz);
end
nnz = length(indnz);

% Contast Effect Size
ces = C*beta;

if(isempty(nacf))
  % Covariance matrix of contrast effect size
  cescvmr = inv(C*inv(X'*X)*C');
  if(dof1 ~= 1) F = (sum(ces .* (cescvmr*ces))./rvar)/dof1;
  else          F = ((ces.^2)./rvar)*(cescvmr/dof1);
  end
  if(nargout == 5) 
    invcescvmr = inv(cescvmr);
    cescvm = invcescvmr(:) * rvar; % outer product
  end
else
  if(isempty(nacfmap)) nacfmap = ones(nnz,1); 
  else                 nacfmap = nacfmap(indnz); 
  end
  nbins = size(nacf,2);
  F = zeros(1,nnz);
  for nthbin = 0:nbins
    indbin = find(nacfmap==nthbin);
    if(isempty(indbin)) continue; end
    if(nthbin ~= 0)
      nacfbin = nacf(:,nthbin);
      %W = inv(chol(toeplitz(nacfbin))');
      W = chol(inv(toeplitz(nacfbin)));
      Xbin = W*X;
    else
      Xbin = X;
    end
    cescvmr  = inv(C*inv(Xbin'*Xbin)*C');
    cesbin  = ces(:,indbin);
    rvarbin = rvar(indbin);
    if(dof1 ~= 1) Fbin = (sum(cesbin .* (cescvmr*cesbin))./rvarbin)/dof1;
    else          Fbin = ((cesbin.^2)./rvarbin)*(cescvmr/dof1);
    end
    F(indbin) = Fbin;
    if(nargout == 5) 
      invcescvmr = inv(cescvmr);
      cescvm(:,indbin) = invcescvmr(:) * rvarbin; % outer product
    end
  end
end

% If there were voxels with rvar==0, fix them
if(~isempty(indz))
  F0 = zeros(1,nv);
  F0(indnz) = F;
  F = F0;
  ces0 = zeros(dof1,nv);
  ces0(:,indnz) = ces;
  ces = ces0;
end

return;



