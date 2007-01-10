function pthresh = fast_fdrthresh(p,fdr)
% pthresh = fast_fdrthresh(p,fdr)
%
% p = list of p values between -1 and 1
% fdr = false discovery rate, between 0 and 1
%
% Based on Tom's FDR.m from 
%   http://www.sph.umich.edu/~nichols/FDR/FDR.m
% The threshold returned from this function is based on an 
% assumption of "independence or positive dependence",
% which should be "reasonable for imaging data".
%
%
%


%
% fast_fdrthresh.m
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

pthresh = [];

if(nargin ~= 2)
  fprintf('pthresh = fast_fdrthresh(p,fdr)\n');
  return;
end

p = sort(abs(p(:)));
Nv = length(p(:));
nn = [1:Nv]';

cVID = 1; % Not sure what this is for

imax = max( find(p <= fdr*nn/Nv ) );
if(~isempty(imax))
  %fprintf('imax = %d\n',imax);
  pthresh = p(imax);
else
  pthresh = min(p);
end

return;
