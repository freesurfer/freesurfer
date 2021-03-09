function [pvs, u, ds] = fast_pvs(y)
% [pvs, u, s] = fast_pvs(y)
% 
% Computes percent variance spanned by each of the eigenvectors
%


%
% fast_pvs.m
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

pvs = [];
u = [];
ds = [];

if(nargin ~= 1)
  fprintf('[pvs, u, s] = fast_pvs(y)\n');
  return;
end

My = y*y'; %'
[u s v] = svd(My);
ds = diag(s);
pvs = 100*ds/sum(ds);

return;
