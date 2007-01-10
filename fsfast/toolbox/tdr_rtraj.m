function r = tdr_rtraj(nrsamples)
% r = tdr_rtraj(nrsamples)
%
% recon space 'trajectory' (returned as a row vector) 
% 
% Starts at -(nrsamples/2-1) (eg, -63 for 128 or -31 for 64)
% Ends at    +nrsamples/2    (eg, +64 for 128 or +32 for 64)
% Increments by 1 each sample.
% Passes through 0 at nrsamples/2.
%
% See also tdr_uniform_phtraj.m
%
%
%


%
% tdr_rtraj.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:35 $
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

r = [];
if(nargin ~= 1)
  fprintf('r = tdr_rtraj(nrsamples)\n');
  return;
end

% Starts at -63, passes thru 0 at nc/2, ends at +64
rstart = -(nrsamples/2-1);
rend   = +nrsamples/2;
r = [rstart:rend];

return;
