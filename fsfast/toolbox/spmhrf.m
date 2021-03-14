function h = spmhrf(t,a1,b1,a2,b2,c)
%
% h = spmhrf(t,a1,b1,a2,b2,c)
%


%
% spmhrf.m
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

h = [];

if(nargin ~= 1 & nargin ~= 6)
  fprintf('USAGE: h = spmhrf(t,a1,b1,a2,b2,c)\n');
  return;
end

if(nargin == 1)
  a1 = 6;
  b1 = 0.9;
  a2 = 12;
  b2 = 0.9;
  c  = .35;
end

d1 = a1*b1;
d2 = a2*b2;

h =    (t/d1).^a1  .* exp(-(t-d1)/b1) - ...
    c*((t/d2).^a2) .* exp(-(t-d2)/b2);

return;



