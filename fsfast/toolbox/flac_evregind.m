function evregind = flac_evregind(flac,nthev)
% evregind = flac_evregind(flac,nthev);
%
% Returns the columns of the regressors in the full design matrix for
% the given EV.
%
%


%
% flac_evregind.m
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

