function [im, neigused] = invillcond(m,pvemin)
% [im neigused]= invillcond(m,pvemin)


%
% invillcond.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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

im = [];
if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: [im neigused]= invillcond(m,pvemin)';
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

% Pickout enough eigenvalues to explain the min desired %
neigused = length(find(pve < pvemin));

% regularize %
mineig = meig(neigused);
meig2 = meig;
meig2(neigused+1:neig) = mineig;

invs = diag(meig2.^(-1));

im = u*invs*v'; %'


return;
