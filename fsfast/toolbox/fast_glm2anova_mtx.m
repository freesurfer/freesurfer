function M = fast_glm2anova_mtx(cflmap,nlevels,ncond,vr,nesttot)
% M = fast_glm2anova_mtx(cflmap,nlevels,<ncond>,<vr>,<nesttot>);
%
% Creates a matrix that will convert a glm parameter vector (beta)
% into a vector of population means suitable for applying
% an ANOVA contrast matrix to. See fast_anovamtx.m.
%
% Each row of the cflmap lists the condition number and the 
% location in the factor tree to which that condition belongs.
% The number of columns will be the number of factors + 1.
% The number of rows must be the product of the number of
% levels (this is the number of population means in the ANOVA).
%
% Eg:  5 2 3 1 would indicate that there are 3 factors and
% that condition 5 is the 2nd branch of factor 1, the 3rd
% branch of factor 2, and the 1st branch of factor 3.
%
% nlevels is the list of the number of levels of each factor.
%
% ncond is an optional argument of the number of conditions.
% If unspecified, the number of conditions is the maximum
% condition number found in the cfl map.
%
%


%
% fast_glm2anova_mtx.m
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

if(nargin < 2 | nargin > 5)
  fprintf('M = fast_glm2anova_mtx(cflmap,nlevels,<ncond>,<vr>,<nesttot>)\n');
  return;
end

[ok, nfactors, npopmeans, popcid] = checkcflmap(cflmap,nlevels);
if(~ok) return; end

if(exist('ncond')~=1) ncond = []; end
if(isempty(ncond))    ncond = max(cflmap(:,1)); end
if(max(cflmap(:,1)) > ncond)
  fprintf('ERROR: cfl cond number exceeds number of conditions\n');
  return;
end


if(exist('vr')~=1) vr = []; end
if(isempty(vr))    vr = 1; end
nestcond = ncond*length(vr);

if(exist('nesttot')~=1) nesttot = []; end
if(isempty(nesttot)) nesttot = nestcond; end

if(nesttot < nestcond)
  fprintf('ERROR: nesttot\n');
  return;
end

Du = zeros(npopmeans,ncond);

for n = 1:npopmeans
  c = popcid(n);
  Du(n,c) = 1;
end

M = Du;

return;

%----------------------------------------------------------------%
function [ok, nfactors, npopmeans, popcid] = checkcflmap(cflmap,nlevels)

ok = 0;
nfactors = length(nlevels);
npopmeans = prod(nlevels);
popcid = [];

if(size(cflmap,2) ~= nfactors+1)
  fprintf('ERROR: mismatch in the number of factors\n');
  return;
end
if(size(cflmap,1) ~= npopmeans)
  fprintf('ERROR: mismatch in the number of population means\n');
  return;
end
if(~isempty(find(reshape1d(cflmap)<1)))
  fprintf('ERROR: cflmap has a component less than one\n');
  return;
end

% check that the factor levels do not exceed maximum %
for f = 1:nfactors
  m = max(cflmap(:,f+1));
  if(m > nlevels(f))
    fprintf('ERROR: factor %d level exceeds limit\n',f);
    return;
  end
end

% check that there are no condition replications %
tmp = unique(cflmap(:,1));
if(length(tmp) ~= length(cflmap(:,1)))
  fprintf('ERROR: cflmap has condition replications\n');
  return;
end

% check that each branch is represented %
popcid = zeros(npopmeans,1);
bv = ones(nfactors,1); % branch vector starts as all ones
for n = 1:npopmeans
  %fprintf('%2d '); fprintf('%d ',bv); fprintf('\n');

  % check that the current branch vector exists and is not replicated
  cfltmp = cflmap;
  for f = 1:nfactors
    ind = find(cfltmp(:,f+1) == bv(f));
    if(isempty(ind))
      fprintf('ERROR: cannot find the following branch vector (f=%d)\n',f);
      fprintf('%d ',bv); fprintf('\n');
      return; 
    end
    if(f==nfactors)
      if(length(ind) > 1)
	fprintf('ERROR: following branch vector is replicated\n');
	fprintf('%d ',bv); fprintf('\n');
	return;
      end
      popcid(n) = cfltmp(ind,1);
      %fprintf('%2d ',n); fprintf('%d ',bv); fprintf('  %d\n',cfltmp(ind,1));
    end
    cfltmp = cfltmp(ind,:); 
  end

  % Increment the current branch vector %
  f = nfactors;
  while(1)
    bv(f) = bv(f) + 1;
    if(bv(f) <= nlevels(f)) break; end
    bv(f) = 1;
    f = f - 1;
    if(f==0) break; end
  end
  
end

ok = 1;

return;
