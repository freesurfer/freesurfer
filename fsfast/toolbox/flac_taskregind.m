function taskregind = flac_taskregind(flac)
% taskregind = flac_taskregind(flac);
%
% Returns the column indices of the task-related regressors.
%
%


%
% flac_taskregind.m
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

taskregind = [];
if(nargin ~= 1)
  fprintf('taskregind = flac_taskregind(flac)\n');
  return;
end

nev = length(flac.ev);

for nthev = 1:nev
  if(strcmp(flac.ev(nthev).type,'task'))
    evregind = flac_evregind(flac,nthev);
    taskregind = [taskregind evregind];
  end
end

return;







