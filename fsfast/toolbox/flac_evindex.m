function index = flac_evindex(flac,evname)
% index = flac_evindex(flac,evname) 
%
% Retuns the index of an EV from its name. Returns empty if evname
% not found.
%
% $Id: flac_evindex.m,v 1.1 2004/10/18 05:52:24 greve Exp $

index = [];

if(nargin ~= 2)
  fprintf('index = flac_evindex(flac,evname)\n');
  return;
end

nevs = length(flac.ev);
for nthev = 1:nevs
  if(strcmp(flac.ev(nthev).name,evname))
    index = nthev;
    return;
  end
end

return;


















