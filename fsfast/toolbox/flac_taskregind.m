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







