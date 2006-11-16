function evtaskind = flac_evtaskind(flac)
% evtaskind = flac_evtaskind(flac)
%
% Returns the indices of the task EVs.
%
% $Id: flac_evtaskind.m,v 1.1 2006/11/16 06:23:30 greve Exp $

evtaskind = [];
if(nargin ~= 1)
  fprintf('evtaskind = flac_evtaskind(flac)\n');
  return;
end

nev = length(flac.ev);
for nthev = 1:nev
  if(strcmp(flac.ev(nthev).type,'task'))
    evtaskind = [evtaskind nthev];
  end
end


return;

