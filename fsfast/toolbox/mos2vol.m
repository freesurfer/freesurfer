function [vol, tszmos] = mos2vol(mos, szvol, tszmos)
%
% [vol tszmos] = mos2vol(mos, szvol, <tszmos>)
%
% Given a mosaic  (rows, cols, planes), produces a 
% volume (rows, cols, slices, planes). If the mosaic has been padded
% with blank images, those will be excluded from the volume.
%
% szvol - size of the volume (rows, cols, slices, <planes>).
%
% tszmos - size (rows, cols) of the mosaic measured in tiles (optional).
% If tszmos is not specified, a default one will be computed using
% the function defmossize.
%
% See also: mos2vol vol2mos mosind2volind mossub2volsub 
%           volind2mosind volsub2mossub defmossize
%
%


%
% mos2vol.m
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
  msg = 'USAGE: [vol tszmos] = mos2vol(mos, szvol, <tszmos>)';
  error(msg);
end

Np  = size(mos,3);

Nvr = szvol(1);
Nvc = szvol(2);
Nvs = szvol(3);
szvol = [Nvr Nvc Nvs];
Nv = prod(szvol(1:3));

if(nargin == 2) tszmos = []; end
tszmos = defmossize(Nvs, tszmos);
Ntr = tszmos(1);
Ntc = tszmos(2);
Nmr = Ntr*Nvr;
Nmc = Ntc*Nvc;
szmos = [Nmr Nmc];
Nm  = prod(szmos);

im = [1:Nm];
[iv im] = mosind2volind(im,szvol,tszmos);

tmp = find(iv <= Nv);
iv = iv(tmp);
im = im(tmp);
mos = reshape(mos, [Nm Np]);
vol = zeros(Nv,Np);
vol(iv,:) = mos(im,:);
vol = reshape(vol, [Nvr Nvc Nvs Np]);


return;
