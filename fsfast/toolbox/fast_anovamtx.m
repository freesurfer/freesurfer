function M = fast_anovamtx(nlevels,factorlist)
% M = fast_anovamtx(nlevels,<factorlist>)
% 
% nlevels is a vector of the number of levels in each 
% factor. nfactors = length(nlevels);
%
% factorlist are the factors to test for an interaction
% across. To test the main effect of a factor, list
% only that factor in the factorlist.
%
% NOTE: I don't think the main effect is working properly.
%
% The order of estimates must be last factor fastest.
% See fast_glm2anova_mtx.
%
%


%
% fast_anovamtx.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.3 $
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

M = [];

if(nargin < 1 | nargin > 2)
  fprintf('M = fast_anovamtx(nlevels,<factorlist>)\n');
  return;
end

nfactors = length(nlevels);
if(exist('factorlist') ~= 1) factorlist = []; end
if(isempty(factorlist)) factorlist = [1:nfactors]; end

% Check that the factor list is ok
if(max(factorlist) > nfactors)
  fprintf('ERROR: bad factor list (max)\n');
  return;
end
if(min(factorlist) < 1)
  fprintf('ERROR: bad factor list (min)\n');
  return;
end

M = 1;
for f = 1:nfactors
  inlist = ~isempty(find(f==factorlist));
  if(inlist) Mf = omnibusmtx(nlevels(f));
  else       Mf = sumvect(nlevels(f));
  end
  M = kron(M,Mf);
end

return;

%-----------------------------------------------%
function O = omnibusmtx(levels)
O = [eye(levels-1) -ones(levels-1,1)];
return;

%-----------------------------------------------%
function s = sumvect(levels)
s = ones(1,levels);
return;






