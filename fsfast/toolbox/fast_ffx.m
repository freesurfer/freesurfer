function [p, gamma, gammavar, F] = fast_ffx(means,varmeans,dofs,X,C);
% [p gamma gammavar F] = fast_ffx(means,varmeans,dofs,X,C);
%
% Group-level fixed effects analysis
%   means    is nsamples-by-nvoxels
%   varmeans is nsamples-by-nvoxels
%   dof = sum(dofs) [if dofs is a scalar, then dof = nsamples*dofs]
%   X is nsamples-by-nX (if not spec or empty, then X = ones(nsamples,1))
%   C is J-by-nX (if not spec or empty, then C = eye(nX))
% 
%   gamma    = inv(X'*X)*X'*means;  (J-by-nvox)
%   gammavar = inv(X'*X)*X'*diag(varmeans)*X*inv(X'*X); (J^2-by-nvox)
%   F = gamma'*inv(gammavar)*gamma/J; (1-by-nvox)
%   p = FTest(J, dof, F); (1-by-nvox)
%
% When X is the ones vector and C=1, then gives the same result as
% fast_ffx_osgm.
%

%
% fast_ffx.m
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
p = []; 
gamma = [];
gammavar = [];
F = [];

if(nargin < 3 | nargin > 6)
  fprintf('[p gamma gammavar F] = fast_ffx(means,varmeans,dofs,X,C);\n');
  return;
end

[nframes nvox] = size(means);

% If X does not exist or is empty, set to all 1s (osgm)
if(~exist('X','var')) X = []; end
if(isempty(X)) X = ones(nframes,1); end
if(size(X,1) ~= nframes)
  fprintf('ERROR: dimension mismatch between X and data\n')';
  return;
end
XtX = X'*X;
condXtX = cond(XtX);
if(condXtX > 1000)
  fprintf('ERROR: XtX condition is %g\n',condXtX);
  return;
end
nX = size(X,2);

% Contrast Matrix
if(~exist('C','var')) C = []; end
if(isempty(C)) C = eye(nX); end
if(size(C,2) ~= nX)
  fprintf('ERROR: dimension mismatch between X and C\n')';
  return;
end
J = size(C,1);

% If dofs is a scalar, mult by nframes to get total dof, otherwise
% just sum the dofs (make sure the frames are the same)
ndofs = length(dofs);
if(ndofs == 1) dof = dofs*nframes;
else 
  if(ndofs ~= nframes)
    fprintf('ERROR: dimension mismatch between dofs and frames\n');
    return;
  end
  dof = sum(dofs); 
end

% Finally, estimate
T = (inv(XtX)*X');
beta = T*means;
gamma = C*beta;
CT = C*T;
CTt = CT';
F = zeros(1,nvox);
gammavar = zeros(J*J,nvox);
for v = 1:nvox
  % Probably a more efficient way to do this
  S = diag(varmeans(:,v)); % CovMtx across samples
  vgammavar = CT*S*CTt;
  F(v) = (gamma(:,v)' * inv(vgammavar) * gamma(:,v))/J;
  gammavar(:,v) = vgammavar(:);
end

p = FTest(J, dof, F);


return;

