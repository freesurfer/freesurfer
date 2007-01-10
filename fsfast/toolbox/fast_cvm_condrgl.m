function [mrgl, neigused, pve] = fast_cvm_condrgl(m,condmin)
% [mrgl, neigused] = fast_cvm_condrgl(m,condmin)
%
% Regularization by specifying a minimum condition number for
% the resulting regularized matrix.
%


%
% fast_cvm_condrgl.m
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
  msg = 'USAGE: [mrgl, neigused] = fast_cvm_condrgl(m,condmin)';
  qoe(msg);error(msg);
end

if(nargin == 1)
  condmin = 0;
end

if(condmin < 1) 
  msg = 'Cannot specify a condition less than 1';
  qoe(msg);error(msg);
end

% Decompose %
[u s v] = svd(m);

% Vector of Eigenvalues (in descending order) %
meig = diag(s); 
neig = length(meig);

% Only look at the non-zero eigenvalues %
inz = find(meig ~= 0);

% Condition after keeping each successive eigenvector
condeig = meig(1)./meig(inz);

% Pickout enough eigenvalues to achive min desired condition %
neigused = length(find(condeig < condmin));

% regularize %
mineig = meig(neigused);
meig2 = meig;
meig2(neigused+1:neig) = mineig;

s2 = diag(meig2);

mrgl = u*s2*v'; %'

% Percent variance explained by neigused %
pve = 100*sum(meig(1:neigused))/sum(meig);

return;
