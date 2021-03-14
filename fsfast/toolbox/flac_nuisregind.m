function nuisregind = flac_nuisregind(flac)
% nuisregind = flac_nuisregind(flac);
%
% Returns the column indices of the nuissance regressors.
%
%


%
% flac_nuisregind.m
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







