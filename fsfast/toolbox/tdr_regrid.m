function [kdatgrid, rgmap] = tdr_regrid(kx,ky,kdat,kgridsize)
% [kdatgrid rgmap] = tdr_regrid(kx,ky,kdat,kgridsize)
% 
% Regrids using nearest neighbor
%
%


%
% tdr_regrid.m
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


kdatgrid = [];

if(nargin ~= 4)
  fprintf('[kdatgrid rgmap] = tdr_regrid(kx,ky,kdat,kgridsize)\n');
  return;
end

xkgridsize = kgridsize(2);
ykgridsize = kgridsize(1);

kr2max = max((kx.^2 + ky.^2));
kxmax = max(kx);
kxmin = min(kx);
kymax = max(ky);
kymin = min(ky);


gx0 = round(xkgridsize/2) + 1;
gy0 = round(ykgridsize/2) + 1;

tmp = [1:xkgridsize];
tmp = tmp - tmp(gx0);
tmp = tmp/max(abs(tmp));
kxg = kxmax * tmp;

tmp = [1:ykgridsize];
tmp = tmp - tmp(gy0);
tmp = tmp/max(abs(tmp));
kyg = kymax * tmp;

kdatgrid = zeros(kgridsize);
rgmap    = zeros(kgridsize);

for gx = 1:xkgridsize
  for gy = 1:ykgridsize
    kr2 = kxg(gx).^2 + kyg(gy).^2;
    if(kr2 > kr2max) continue; end
    dx = kx - kxg(gx);
    dy = ky - kyg(gy);
    d2 = dx.^2 + dy.^2;
    [blah nnbr] = min(d2);
    kdatgrid(gy,gx) = kdat(nnbr);
    rgmap(gy,gx) = nnbr;
  end
end

return;





