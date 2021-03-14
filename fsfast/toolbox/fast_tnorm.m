function ftnorm = fast_tnorm(f)
% imgtnorm = fast_tnorm(f)
% 
% Temporal Normalization f: N1xN2x...Nt
%


%
% fast_tnorm.m
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

if(nargin ~= 1)
  msg = 'imgtnorm = fast_tnorm(f)';
  qoe(msg);error(msg);
end

szf = size(f);
ndimf = length(szf);

aszf = szf(1:ndimf-1);
nvox = prod(aszf);
nt   = szf(ndimf);

f = reshape(f, [nvox nt]);
fmn = mean(f,2);
fstd = std(f,[],2);

ftnorm = (f - repmat(fmn,[1 nt])) ./ repmat(fstd,[1 nt]) ;

ftnorm = reshape(ftnorm,szf);

return;
