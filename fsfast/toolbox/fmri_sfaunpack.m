function [nstd, hr, hi] = fmri_sfaunpack(sfaslice)
% [nstd hr hi] = fmri_sfaunpack(sfaslice)
%
% sfaslice: nrows ncols nplanes
%
% $Id: fmri_sfaunpack.m,v 1.1 2003/03/04 20:47:40 greve Exp $

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
