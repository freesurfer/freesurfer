function Xbaseline = fast_baselinemtx(run,ntrs,nruns)
% Xbaseline = fast_baselinemtx(run,ntrs,nruns)


%
% fast_baselinemtx.m
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

if(nargin ~= 3)
  msg = 'USAGE: Xbaseline = fast_baselinemtx(run,ntrs,nruns)';
  qoe(msg);error(msg);
end

v = ones(ntrs,1);

Xbaseline        = zeros(ntrs,nruns);
Xbaseline(:,run) = v;

return;
