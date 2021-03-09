function rimg1d = rimgfunc(kimg1d,rimgsize,kimgsize,coilprof)
% rimg1d = rimgfunc(kimg1d,rimgsize,kimgsize,<coilprof>)
% backprojection model


%
% rimgfunc.m
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

if(~exist('coilprof','var')) coilprof = []; end

Nk = prod(kimgsize);
Faccel = rimgsize(1)/kimgsize(1);

if(isempty(coilprof))
  kimg = reshape(kimg1d,kimgsize);
  if(Faccel ~= 1)
    kimgtmp = zeros([kimgsize(1)*Faccel kimgsize(2)]);
    kimgtmp(1:Faccel:end,:) = kimg;
    kimg = kimgtmp;
  end
  rimg1d = ifft2(fftshift(kimg)); 
  rimg1d = rimg1d(:);
else
  ncoils = size(coilprof,3);
  if(Faccel ~= 1)
    kimg1dtmp = zeros(Faccel*Nk*ncoils,1);
    kimg1dtmp(1:Faccel:end) = kimg1d;
    kimg1d = kimg1dtmp;
  end
  dk = Nk*Faccel;
  rimg = 0;
  nk1 = 1;
  for nthcoil = 1:ncoils
    nk2 = nk1 + dk - 1;
    kimgcoil = reshape(kimg1d(nk1:nk2),[kimgsize(1)*Faccel kimgsize(2)]);
    rimgcoil = ifft2(fftshift(kimgcoil)) .* conj(coilprof(:,:,nthcoil)); 
    rimg = rimg + rimgcoil;
    nk1 = nk2 + 1;
  end
  rimg1d = rimg(:);
end

%rimg1d = ((Faccel.^2)*Nk)*rimg1d;
rimg1d = ncoils*Nk*rimg1d;

return;
