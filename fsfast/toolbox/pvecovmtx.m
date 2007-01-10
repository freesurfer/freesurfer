function [pve, s] = pvecovmtx(m)
% [pve s] = pvecovmtx(m)
%
% Percent Variance Explained.


%
% pvecovmtx.m
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

if(nargin ~= 1) 
  msg = 'USAGE: [pve s] = pvecovmtx(m)';
  qoe(msg); error(msg);
end

s = svd(m);

pve = 100*cumsum(s)/sum(s);

return;
