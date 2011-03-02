function index = flac_evindex(flac,evname)
% index = flac_evindex(flac,evname) 
%
% Retuns the index of an EV from its name. Returns empty if evname
% not found.
%
%


%
% flac_evindex.m
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


















