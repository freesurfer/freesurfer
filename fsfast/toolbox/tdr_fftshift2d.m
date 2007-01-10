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
%    $Date: 2007/01/10 22:02:35 $
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


