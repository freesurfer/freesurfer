function [mrgl, neigused, mrglcond] = fast_cvm_pvergl(m,pvemin)
% [mrgl, neigused, mrglcond] = fast_cvm_pvergl(m,pvemin)
%
% Regularization by setting a minimum on the eigenvalues equal
% to the the eigenvalue of the nth eigenvector needed to
% account for pvemin percent of the variance.
%


%
% fast_cvm_pvergl.m
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

mrgl = [];
if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: [mrgl neigused mrglcond]= fast_cvm_pvergl(m,pvemin)';
  qoe(msg);error(msg);
end

if(nargin == 1)
  pvemin = 100;
end

% Decompose %
[u s v] = svd(m);

% Vector of Eigenvalues %
meig = diag(s);
neig = length(meig);

% Cumulative Percent Variance Explained
pve = 100*cumsum(meig)/sum(meig);

% Pickout enough eigenvalues to explain the min desired variance %
neigused = length(find(pve < pvemin));

% rglularize %
mineig = meig(neigused);
meig2 = meig;
meig2(neigused+1:neig) = mineig;

s2 = diag(meig2);

mrgl = u*s2*v'; %'

mrglcond = cond(mrgl);

return;
