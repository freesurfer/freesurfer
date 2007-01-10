function kimg1d = kimgfunc(rimg1d,rimgsize,kimgsize,coilprof)
% kimg1d = kimgfunc(rimg1d,rimgsize,kimgsize,<coilprof>)
% Forward model


%
% kimgfunc.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(~exist('coilprof','var')) coilprof = []; end

rimg = reshape(rimg1d,rimgsize);

OverSampRate = kimgsize(2)/rimgsize(2);
if(OverSampRate > 1)
  rimgzpad = zeros(rimgsize(1),kimgsize(2));
  Nzpadtot = kimgsize(2) - rimgsize(2);
  Nzpadside = Nzpadtot/2;
  Nc1 = Nzpadside+1;
  Nc2 = Nc1 + rimgsize(2);
end

if(isempty(coilprof))
  if(OverSampRate > 1)
    rimgzpad(:,Nc1:Nc2) = rimg;
    kimg1d = fft2(rimgzpad);
  else
    kimg1d = fftshift(fft2(rimg));
  end
  kimg1d = kimg1d(:);
else
  kimg1d = [];
  ncoils = size(coilprof,3);
  for nthcoil = 1:ncoils
    rimgcoil = rimg .* coilprof(:,:,nthcoil);
    %rimgcoil = rimg;
    if(OverSampRate > 1)
      rimgzpad(:,Nc1:Nc2) = rimgcoil;
      kimgcoil = fftshift(fft2(rimgzpad));
    else
      kimgcoil = fftshift(fft2(rimgcoil));
    end
    kimg1d = [kimg1d; kimgcoil(:)];
  end
  Faccel = rimgsize(1)/kimgsize(1);
  if(Faccel ~= 1)
    kimg1d = kimg1d(1:Faccel:end);
  end
end

return;
