function taskregind = flac_taskregind(flac)
% taskregind = flac_taskregind(flac);
%
% Returns the column indices of the task-related regressors.
%
% $Id: flac_taskregind.m,v 1.1 2004/10/17 18:31:37 greve Exp $

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







