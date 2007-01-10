function [cnd, mineig, S] = fast_acfcond(acf,taper)
% [cnd mineig S] = fast_acfcond(acf,<taper>)


%
% fast_acfcond.m
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

cnd = [];
mineig = [];
S = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('[cnd mineig S] = fast_acfcond(acf,<taper>)\n');
  return;
end

if(exist('taper') ~= 1) taper = []; end
if(~isempty(taper))
  acf = acf .* taper;
end

S = toeplitz(acf);
mineig = min(eig(S));
cnd = cond(S);

return;







