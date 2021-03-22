function [pvs, u, acf] = fast_sim_acfsvd(nf, nv, X, nmax)
% [pvs u acf] = fast_sim_acfsvd(nf, nv, X, nmax)
%
% Computes the SVD of the autocorrelation function (ACF)
% of white, gaussian noise with nf points and nv exemplars.
% If X is included, X is used to detrend the white noise.
% If nmax is included, the svd is computed over the first
% nmax ACF delays. The first delay is always exluded 
% since it is always 1.
%
% pvs - percent variance spanned by each eigenvectors
% u - eigenvectors.
% acf - without delay 1


%
% fast_sim_acfsvd.m
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

pvs = [];

if(nargin < 2 | nargin > 4)
  fprintf('[pvs u]= fast_sim_acfsvd(nf, nv, X, nmax)\n');
  return;
end

if(nf > nv)
  fprintf('ERROR: nf > nv\n');
  return;
end

if(nargin < 4) nmax = nf; end

if(nargin > 2)
  if(nf ~= size(X,1))
    fprintf('ERROR: nf/X dimension mismatch\n');
    return;
  end
  R = eye(nf) - X*inv(X'*X)*X';
  y = R*randn(nf,nv);
else
  y = randn(nf,nv);
end


acf = fast_acorr(y);
acf = acf(2:nmax,:);
Macf = acf*acf'; %'
[u s v] = svd(Macf);
ds = diag(s);
pvs = 100*ds/sum(ds);


return;


