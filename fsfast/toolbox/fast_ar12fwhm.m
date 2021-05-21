function fwhm = fast_ar12fwhm(ar1,d)
% fwhm = fast_ar12fwhm(ar1,d)
%
% Converts an AR1 to FWHM. 
% d is the voxel size. FWHM will be in units of d.
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

fwhm = [];
if(nargin ~= 1 & nargin ~= 2)
  fprintf('ar1 = fast_ar12fwhm(ar1,d)\n');
  return;
end

fwhm = sqrt(log(256.0))*d./sqrt(-4*log(ar1));

return
