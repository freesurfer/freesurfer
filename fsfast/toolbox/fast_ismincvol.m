function isminc = fast_ismincvol(volid)
% isminc = fast_ismincvol(volid)
%
% Returns 1 if volume is in MINC format.
% 
%


%
% fast_ismincvol.m
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

isminc = 0;

if(nargin ~= 1)
  msg = 'USAGE: r = fast_ismincvol(volid)'
  qoe(msg); error(msg);
end

volid = deblank(volid);
len = length(volid);

% Cannot be MINC if less than 5 characters %
if(len < 5) isminc = 0; return; end

last4 = volid(len-3:len);
if(strcmp(last4,'.mnc')) isminc = 1; return; end

if(len < 8) isminc = 0; return; end

last7 = volid(len-6:len);
if(strcmp(last7,'.mnc.gz')) isminc = 1; return; end

isminc = 0;

return;
