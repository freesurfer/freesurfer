function [nstd, hr, hi] = fmri_sfaunpack(sfaslice)
% [nstd hr hi] = fmri_sfaunpack(sfaslice)
%
% sfaslice: nrows ncols nplanes
%
%


%
% fmri_sfaunpack.m
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
  msg = '[nstd hr hi] = fmri_sfaunpack(sfaslice)';
  qoe(msg); error(msg);
end

[nrows ncols nplanes] = size(sfaslice);

nstd = sfaslice(:,:,1);

nfreq = (nplanes-1)/2;
indreal = [2:2+(nfreq-1)];
indimag = indreal + nfreq;
hr = sfaslice(:,:,indreal);
hi = sfaslice(:,:,indimag);


return;
