function [F, dof1, dof2, ces, cescvm, pcc] = fast_fratiow(beta,X,rvar,C,nacf,nacfmap)
% [F dof1 dof2 ces cescvm pcc] = fast_fratiow(beta,X,rvar,C,<nacf>,<nacfmap>);
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
% pcc - partial correlation coefficient (for univariate contrasts)
%
% See also: fast_glmfitw, FTest, fast_glmfit, fast_fratio.
%
%


%
% fast_fratiow.m
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

pcc = [];

if(nargin < 4 | nargin > 6)
  fprintf('[F dof1 dof2 ces cescvm] = fast_fratiow(beta,X,rvar,C,<nacf>,<nacfmap>)\n');
  return;
end

[nf nbeta] = size(X);
nvox = size(beta,2);

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

dof1 = size(C,1); % J %
dof2 = nf-nbeta;

% Handle case where some voxels are zero
indz  = find(rvar == 0);
indnz = find(rvar ~= 0);
nnz = length(indnz);
if(~isempty(indz))
  beta = beta(:,indnz);
  rvar = rvar(:,indnz);
end

% Contast Effect Size
ces = C*beta;
usematrix = 0;

%if(usematrix) fprintf('Using matrix filtering\n');
%else          fprintf('Using fft filtering\n');
%end
X_fft = fft(X,2*nf);
if(isempty(nacf))
  % Covariance matrix of contrast effect size
  cescvmr = inv(C*inv(X'*X)*C');
  if(dof1 ~= 1) F = (sum(ces .* (cescvmr*ces))./rvar)/dof1;
  else          
    F = ((ces.^2)./rvar)*(cescvmr/dof1);
    if(nargout >= 6) pcc = fast_glm_pcc(beta,X,C,rvar); end
  end
  if(nargout >= 5) 
    invcescvmr = inv(cescvmr);
    cescvm = zeros(dof1.^2,nnz);
    cescvm = invcescvmr(:) * rvar; % outer product
  end
else
  if(isempty(nacfmap)) nacfmap = ones(nnz,1); 
  else                 nacfmap = nacfmap(indnz); 
  end
  nbins = size(nacf,2);
  F = zeros(1,nnz);
  if(nargout == 6) pcc = zeros(1,nnz); end
  cescvm = zeros(dof1,nnz); % bug was here (nvox instead of nnz)
  for nthbin = 0:nbins
    indbin = find(nacfmap==nthbin);
    nbin = length(indbin);
    if(isempty(indbin)) continue; end
    if(nthbin ~= 0)
      nacfbin = nacf(:,nthbin);
      if(usematrix)
	% Matrix formulation
	W = chol(inv(toeplitz(nacfbin)));
	Xbin = W*X;
      else
	% FFT formulation
	nacfbin_fft = fft(nacfbin,2*nf);
	Xbin_fft = X_fft./repmat(conj(nacfbin_fft),[1 nbeta]);
	Xbin = real(ifft(Xbin_fft));
	Xbin = Xbin(1:nf,:);
	clear Xbin_fft;
      end
    else
      Xbin = X;
    end
    cescvmr  = inv(C*inv(Xbin'*Xbin)*C');
    cesbin  = ces(:,indbin);
    rvarbin = rvar(indbin);
    if(dof1 ~= 1) Fbin = (sum(cesbin .* (cescvmr*cesbin))./rvarbin)/dof1;
    else          
      Fbin = ((cesbin.^2)./rvarbin)*(cescvmr/dof1);
      if(nargout == 6) pccbin = fast_glm_pcc(beta(:,indbin),Xbin,C,rvarbin); end
    end
    F(indbin) = Fbin;
    if(nargout >= 5) 
      invcescvmr = inv(cescvmr);
      cescvm(:,indbin) = invcescvmr(:) * rvarbin; % outer product
    end
    if(nargout >= 6) pcc(indbin) = pccbin; end
  end
end

% If there were voxels with rvar==0, fix them
if(~isempty(indz))
  F0 = zeros(1,nvox);
  F0(indnz) = F;
  F = F0;
  ces0 = zeros(dof1,nvox);
  ces0(:,indnz) = ces;
  ces = ces0;
  if(nargout >= 5) 
    cescvm0 = zeros(dof1,nvox);
    cescvm0(:,indnz) = cescvm;
    cescvm = cescvm0;
  end
  if(nargout >= 6) 
    pcc0 = zeros(1,nvox);
    pcc0(indnz) = pcc;
    pcc = pcc0;
  end
end

return;



