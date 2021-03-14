function [mos, tszmos] = vol2mos(vol, tszmos)
% [mos tszmos] = vol2mos(vol, tszmos)
%
% Given a volume (rows, cols, slices, planes), produces a mosaic
% (rows, cols, planes). vol can be an mristruct.
%
% tszmos - size (rows, cols) of the mosaic measured in tiles (optional).
% If tszmos is not specified, a default one will be computed using
% the function defmossize.  For example, if there are 16 slices, each
% 64X64, then the mosaic will be 4 tiles by 4 tiles which is 256 rows 
% by 256 columns.
%
% If the total number of tiles exceeds the number of volume slices,
% the mosaic is padded with blank images.
%
% See also: mos2vol vol2mos mosind2volind mossub2volsub 
%           volind2mosind volsub2mossub defmossize
%
%


%
% vol2mos.m
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: [mos tszmos] = vol2mos(vol, <tszmos>)';
  error(msg);
end

if(isfield(vol,'vol'))  vol = vol.vol(:,:,:,1); end

% Get the dimensions of the volume %
Nvr = size(vol,1);
Nvc = size(vol,2);
Nvs = size(vol,3);
szvol = [Nvr Nvc Nvs];
Nv = prod(szvol(1:3)); % Total number of pixels in the volume

% Get number of planes 
Np = size(vol,4);

% Compute the dimensions of the mosaic %
% First compute the mosaic size in tiles
if(nargin == 1) tszmos = []; end
tszmos = defmossize(Nvs, tszmos);
Ntr = tszmos(1);
Ntc = tszmos(2);
% Now compute the mosaic size in pixels
Nmr = Ntr*Nvr;
Nmc = Ntc*Nvc;
szmos = [Nmr Nmc];

% Total number of pixels in the volume. Note this may not
% be the same as Nv if the mosaic has been padded with blank
% images.
Nm  = prod(szmos);

% List of all volume indices
iv = [1:Nv];

% List of all mosiac indices that correspond, index-by-index
% to those in iv.
im = volind2mosind(iv,szvol,tszmos);

% Reshape the volume into 2D: volume elements X planes
% This makes it easy to map the volume elements to mosaic elements
vol = reshape(vol, [Nv Np]);

% Create a 2D mosaic: volume elements X planes
mos = zeros(Nm, Np);
mos(im,:) = vol(iv,:);
clear vol;

% Reshape the mosaic back to 3D: rows, cols, planes
mos = reshape(mos, [Nmr Nmc Np]);

return;
