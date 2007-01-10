function parbox = fmri_boxcar(par,boxduration,TR)
%
% parbox = fmri_boxcar(par,boxduration,TR)
%


%
% fmri_boxcar.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(nargin ~= 3)
  msg = 'USAGE: parbox = fmri_boxcar(par,boxduration,TR)';
  qoe(msg); error(msg);
end

Nruns = size(par,3);
Nbox  = boxduration/TR - 1;
boxcar = ones(Nbox,1);

for run = 1:Nruns
  p = par(:,:,run);
  stimid = unique(p(:,2)); 
  nstim  = length(stimid);
  tmin   = min(p(:,1));
  tmax   = max(p(:,1));
  t      = [tmin:TR:tmax+boxduration]';

  for stim = 1:nstim
    istim = find(p(:,2)==stimid(stim));
    tstim = p(istim,1);
    itstim = (tstim-tmin)/2 + 1;
    s = zeros(size(t));
    s(itstim) = 1;
    sb = conv(s,boxcar);
    t2 = tmin + TR*[0:length(sb)-1]';
keyboard
  end

end

return;
