function [iv, im, tszmos] = mosind2volind(im, szvol, tszmos)
% [iv tszmos] = mosind2volind(im, szvol, tszmos)


%
% mosind2volind.m
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

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: [iv tszmos] = mosind2volind(im, szvol, tszmos)';
  error(msg);
end

szvol = szvol(1:3);
Nvr = szvol(1);
Nvc = szvol(2);
Nvs = szvol(3);

% Size of Mosaic measured in Tiles %
if(nargin == 2) tszmos = []; end
tszmos = defmossize(Nvs, tszmos);
Ntr = tszmos(1);
Ntc = tszmos(2);

% Size of Mosaic measured in Elements %
Nmr = Ntr*Nvr;
Nmc = Ntc*Nvc;
szmos = [Nmr Nmc];

[rm cm] = ind2sub(szmos,im);
[rv cv sv im] = mossub2volsub(rm,cm,szvol,tszmos);

iv = sub2ind(szvol, rv, cv, sv);

return;
