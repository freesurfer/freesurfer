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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.2 $
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
