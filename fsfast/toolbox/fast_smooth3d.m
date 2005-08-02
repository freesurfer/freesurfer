function [volsm, G] = fast_smooth3d(vol,cfwhm,rfwhm,sfwhm)
% [volsm, G] = fast_smooth3d(vol,cfwhm,rfwhm,sfwhm)
% 
% 3D gaussian smoother.
%
% cfwhm - fwhm for cols
% rfwhm - fwhm for rows
% sfwhm - fwhm for slice
%
% Note: does not attempt to handle wrap-around
%
% $Id: fast_smooth3d.m,v 1.3 2005/08/02 00:13:47 greve Exp $

if(nargin ~= 4)
  fprintf('[volsm, G] = fast_smooth3d(vol,cfwhm,rfwhm,sfwhm)\n');
  return;
end

if(isfield(vol,'vol')) vol = vol.vol; end

rstd = rfwhm/sqrt(log(256.0));
cstd = cfwhm/sqrt(log(256.0));
sstd = sfwhm/sqrt(log(256.0));

[nrows ncols nslices nframes] = size(vol);

% Do the rows:
if(rstd > 0)
  vol = reshape(vol,[nrows ncols*nslices*nframes]);
  vol = fast_smooth1d(vol,rstd);
  vol = reshape(vol,[nrows ncols nslices nframes]);
end

% Do the cols:
if(cstd > 0)
  vol = permute(vol,[2 1 3 4]);
  vol = reshape(vol,[ncols nrows*nslices*nframes]);
  vol = fast_smooth1d(vol,cstd);
  vol = reshape(vol,[ncols nrows nslices nframes]);
  vol = permute(vol,[2 1 3 4]);
end

% Do the slices
if(sstd > 0)
  vol = permute(vol,[3 2 1 4]);
  vol = reshape(vol,[nslices ncols*nrows*nframes]);
  vol = fast_smooth1d(vol,sstd);
  vol = reshape(vol,[nslices ncols nrows nframes]);
  vol = permute(vol,[3 2 1 4]);
end

volsm = vol;

return


volsz = [nrows ncols nslices];
volsz2 = 2*volsz;

y = ([1:nrows]' - (nrows/2 + 1))/rstd;
yy = repmat(y,[1 ncols nslices]);
x = ([1:ncols]' - (ncols/2 + 1))/cstd;
xx = repmat(x,[1 nrows nslices]);
xx = permute(xx, [2 1 3]);
z = ([1:nslices]' - (nslices/2 + 1))/sstd;
zz = repmat(z,[1 nrows ncols]);
zz = permute(zz, [2 3 1]);

G = exp( -(xx.^2)/2 + -(yy.^2)/2 + -(zz.^2)/2 );
G = G/max(G(:));

kG    = fftn(G,volsz);

% Compute scaling factor. Gotta be a better way.
% Assures that if you put in a volume of ones,
% you get out a volume of ones.
v1 = ones(nrows,ncols,nslices);
kv1  = fftn(v1(:,:,:),volsz);
v1sm = fftshift(ifftn(kv1.*kG));
v1sm = abs(v1sm+1e5) - 1e5;
s = 1/mean(v1sm(:));

volsm = zeros(size(vol));
for frame = 1:nframes
  kvol  = fftn(vol(:,:,:,frame));
  volsm(:,:,:,frame) = fftshift(ifftn(kvol.*kG));
end

if( max(abs(imag(vol(:)))) == 0 )
  % Only take the magnitude if input is totally real
  volsmmax = max(abs(volsm(:)));
  % Shift up to avoid tasking abs of neg number
  volsm = volsm + 10*volsmmax;
  volsm = abs(volsm);
  volsm = volsm - 10*volsmmax;
end

volsm = s*volsm;
