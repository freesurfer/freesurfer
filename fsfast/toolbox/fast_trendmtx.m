function Xtrend = fast_trendmtx(run,ntrs,nruns)
% Xbasline = fast_trendmtx(run,ntrs,nruns)


%
% fast_trendmtx.m
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

if(nargin ~= 3)
  msg = 'USAGE: Xbasline = fast_trendmtx(run,ntrs,nruns)';
  qoe(msg);error(msg);
end

v = [0:ntrs-1]'; %'
v = v - mean(v);
v = v./sqrt(sum(v.^2));

Xtrend        = zeros(ntrs,nruns);
Xtrend(:,run) = v;

return;
