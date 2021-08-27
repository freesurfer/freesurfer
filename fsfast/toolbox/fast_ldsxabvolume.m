function [havg, eresvar, sxadat] = fast_ldsxabvolume(sxastem)
%
% [havg eresvar sxadat] = fast_ldsxabvolume(sxabfile)
%
% This function reads in the selxavg values from the given bvolume
% assuming that the data are stored in selavg format.
%
% havg - (nslices,nrows,ncols,Nhtot)
% eresvar - (nslices,nrows,ncols) - residual error variance
% sxadat - info from the .dat file
%
%
%
% See also: fast_svsxabfile()


%
% fast_ldsxabvolume.m
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

havg = [];
eresvar = [];
sxadat = [];

if(nargin ~= 1) 
  msg = 'USAGE: [havg eresvar sxadat] = fast_ldsxabvolume(sxastem)';
  qoe(msg); error(msg);
end

datfile = sprintf('%s.dat',sxastem);
sxadat = fmri_lddat3(datfile);
if(isempty(sxadat)) return; end

[nslices nrows ncols nt endian bext sxadat] = fmri_bvoldim(sxastem);
eresvar = zeros(nslices,nrows,ncols);
havg = zeros(nslices,nrows,ncols,sxadat.Nh*sxadat.Nnnc);

for slice = 0:nslices-1
  fname = sprintf('%s_%03d.bfloat',sxastem,slice);
  [havg(slice+1,:,:,:) eresvar(slice+1,:,:)] = fast_ldsxabfile(fname);
end

return;
