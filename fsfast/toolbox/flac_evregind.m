function evregind = flac_evregind(flac,nthev)
% evregind = flac_evregind(flac,nthev);
%
% Returns the columns of the regressors in the full design matrix for
% the given EV.
%
% $Id: flac_evregind.m,v 1.2 2004/10/25 04:48:58 greve Exp $

evregind = [];
if(nargin ~= 2)
  fprintf('evregind = flac_evregind(flac,nthev);\n');
  return;
end

nev = length(flac.ev);
if(nthev > nev) 
  fprintf('ERROR: requested nthEV=%d > nEVs=%d\n',nthev,nev);
  return;
end

nregcum = 1;
for mthev = 1:nthev-1
  nregcum = nregcum + flac.ev(mthev).nreg;
end

evregind = [nregcum:nregcum+flac.ev(nthev).nreg-1];

return;

