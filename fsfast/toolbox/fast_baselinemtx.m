function Xbaseline = fast_baselinemtx(run,ntrs,nruns)
% Xbaseline = fast_baselinemtx(run,ntrs,nruns)


%
% fast_baselinemtx.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:03 $
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

if(nargin ~= 3)
  msg = 'USAGE: Xbaseline = fast_baselinemtx(run,ntrs,nruns)';
  qoe(msg);error(msg);
end

v = ones(ntrs,1);

Xbaseline        = zeros(ntrs,nruns);
Xbaseline(:,run) = v;

return;
