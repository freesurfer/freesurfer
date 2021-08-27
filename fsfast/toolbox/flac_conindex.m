function ind = flac_conindex(conname,flac)
% ind = flac_conindex(conname,flac)
% Returns the index of the given contrast name in the flac
%


%
% flac_conindex.m
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

ind = [];
if(nargin ~= 2)
  fprintf('ind = flac_conindex(conname,flac)\n');
  return;
end

ncon = length(flac.con);
for nthcon = 1:ncon
  if(strcmp(flac.con(nthcon).name,conname))
    ind = nthcon;
    return;
  end
end

% Will return empty ind if gets here.

return;
















