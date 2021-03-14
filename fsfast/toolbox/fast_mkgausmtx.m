function G = fast_mkgausmtx(nstddev,len)
% G = fast_mkgausmtx(nstddev,len)
% Create a len-by-len guassian filter matrix 


%
% fast_mkgausmtx.m
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

G = [];

if(nargin ~= 2)
  fprintf('G = fast_mkgausmtx(nstddev,len)\n');
  return;
end

for n = 1:len
  g = fast_gaussian(n,nstddev,len);
  % g = g/sqrt(sum(g.^2));
  g = g/sum(abs(g));
  G = [G; g];
end

return;
