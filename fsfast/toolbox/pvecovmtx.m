function [pve, s] = pvecovmtx(m)
% [pve s] = pvecovmtx(m)
%
% Percent Variance Explained.


%
% pvecovmtx.m
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

if(nargin ~= 1) 
  msg = 'USAGE: [pve s] = pvecovmtx(m)';
  qoe(msg); error(msg);
end

s = svd(m);

pve = 100*cumsum(s)/sum(s);

return;
