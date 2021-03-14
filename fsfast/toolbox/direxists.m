function r = direxists(dirpath)
% r = direxists(dirpath)
%
% returns 1 if dirpath exists, 0 else
%
%


%
% direxists.m
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

r = -1;
if(nargin ~= 1)
  fprintf('r = direxists(dirpath)\n');
  return;
end

d = dir(dirpath);
if(size(d,1) == 0) r = 0;
else               r = 1;
end

return;





