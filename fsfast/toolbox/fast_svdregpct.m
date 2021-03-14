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

















