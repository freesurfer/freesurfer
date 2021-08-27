function [indbest, cpvs] = fast_tspvsrank(dt,ds,nbest)
% [indbest cpvs] = fast_tspvsrank(dt,ds,<nbest>)
%
% Ranks each column in the source data (ds) according to the percent
% variance in the target data (dt) that is spanned by that column.
% After the best column is chosen, both the source and target data are
% orthogonalized with respect to that column so that the second best
% is measured with respect to the variance unaccounted for by the
% first, etc.
%
% If nbest is set, the only the nbest are ranked instead of ranking
% all columns of ds.
%
% indbest - columns ranked from best to worst
% cpvs - cumulative percent variance spanned by columns. If dt
%   and ds are independent, then cpvs will increment by 1/dim(dt)
%   for each column of ds.
%
% Eg, if indbest(1) = 5 and cpvs(1) = 10.2, then the 5th column of
% ds will span 10.2% of the variance in dt.
%
%


%
% fast_tspvsrank.m
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

indbest = [];
cpvs = [];

if(nargin < 2 | nargin > 3)
  fprintf('[indbest, cpvs] = fast_tspvsrank(dt,ds,nbest)\n');
  return;
end

if(size(dt,1) ~= size(ds,1))
  fprintf('ERROR: dt and ds have a different number of rows\n');
  return;
end

ndt = size(dt,2);
nds = size(ds,2);
if(~exist('nbest','var')) nbest = nds; end

% This yields the same result but is faster
if(size(dt,1) < size(dt,2))
  % Number of rows < cols, reduce to rows-by-rows
  [u s] = svd(dt);
  diags = diag(s)';
  dt = u .* repmat(diags, [size(u,2) 1]);
end

% Rescale dt so that var is 100
%dt = 10*dt/sqrt(mean(dt(:).^2));
dt = 10*dt/sqrt(mean(var(dt)));

indbest = [];
for nthbest = 1:nbest
  %fprintf('    nthbest = %d (%g)\n',nthbest,toc);

  % Go through each col of source data
  for nthds = 1:nds 
    c = ds(:,nthds); % candidate column
    if(nthbest > 1 & ~isempty(find(indbest==nthds))) 
      % Skip this column if it is already one of the best
      rtvarmn(nthds) = 10^10;
      continue; 
    end
    % Orthog dt wrt c
    rt = dt - c*((inv(c'*c)*c')*dt); 
    % Compute residual variance
    %rtvarmn(nthds) = mean(rt(:).^2);
    rtvarmn(nthds) = mean(var(rt));
  end

  % Find source column that reduce var the most
  [rtvarmnmin imin] = min(rtvarmn);
  indbest(nthbest) = imin;
  cpvs(nthbest) = 100-rtvarmnmin;

  if(nthbest < nbest)
    % Now orthog ds and dt wrt best
    c = ds(:,imin);
    ds = ds - c*((inv(c'*c)*c')*ds);
    dt = dt - c*((inv(c'*c)*c')*dt);
  end

end

return;




