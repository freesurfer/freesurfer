function vols = tdr_fftshift2d(vol)
% vols = tdr_fftshift2d(vol)
%
% Performs 2d fft shift on each slice/frame
%
%


%
% tdr_fftshift2d.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:07 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

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


