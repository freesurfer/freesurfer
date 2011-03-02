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
%    $Date: 2011/03/02 00:04:05 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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







