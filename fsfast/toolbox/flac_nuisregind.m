function nuisregind = flac_nuisregind(flac)
% nuisregind = flac_nuisregind(flac);
%
% Returns the column indices of the nuissance regressors.
%
% $Id: flac_nuisregind.m,v 1.1 2004/11/02 05:15:59 greve Exp $

nuisregind = [];
if(nargin ~= 1)
  fprintf('nuisregind = flac_nuisregind(flac)\n');
  return;
end

nev = length(flac.ev);

for nthev = 1:nev
  if(strcmp(flac.ev(nthev).type,'nuis'))
    evregind = flac_evregind(flac,nthev);
    nuisregind = [nuisregind evregind];
  end
end

return;







