function [mreg, mdim] = fast_svdregpct(m,pct)
% [mreg, mdim] = fast_svdregpct(m,pct)
%
% Regularizes a matrix by choosing enough eigenvectors to span pct
% percent of the variance.
%
%


%
% fast_svdregpct.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

mreg = [];
if(nargin ~= 2)
  fprintf('[mreg, mdim] = fast_svdregpct(m,pct)\n');
  return;
end

[u s v] = svd(m);
ds = diag(s);
pvs = 100*ds/sum(ds);
cpvs = cumsum(pvs);

mdim = min(find(cpvs > pct));
ds2 = ds;
ds2(mdim:end) = ds2(mdim);

mreg = u * diag(ds2) * v';

return;

















