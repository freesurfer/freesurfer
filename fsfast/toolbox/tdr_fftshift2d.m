function vols = tdr_fftshift2d(vol)
% vols = tdr_fftshift2d(vol)
%
% Performs 2d fft shift on each slice/frame
%
% $Id: tdr_fftshift2d.m,v 1.1 2005/03/19 00:26:05 greve Exp $

if(nargin ~= 1)
  fprintf('vols = tdr_fftshift2d(vol)\n');
  return;
end

[nr nc ns nf] = size(vol);
vols = zeros(size(vol));

for f = 1:nf
  for s = 1:ns
    vols(:,:,s,f) = fftshift(vol(:,:,s,f));
  end
end


return;


