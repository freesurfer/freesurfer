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






