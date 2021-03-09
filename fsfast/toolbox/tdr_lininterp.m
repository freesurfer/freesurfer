function ysamp = tdr_lininterp(y,tsamp,dt,toffset)
% tdr_lininterp(y,tsamp,<dt>,<toffset>)
%
% 1D linear interpolation.
%
% y is uniformly sampled starting at t=toffset with
% increment dt. If dt is not specified, dt=1. If
% toffset is not specified, toffset = 0;
%
% ysamp is y sampled at tsamp.
%
%


%
% tdr_lininterp.m
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

ysamp = [];

if(nargin < 2 | nargin > 4)
  fprintf('tdr_lininterp(y,tsamp,<dt>,<toffset>)\n');
  return;
end
  
if(exist('dt') ~= 1) dt = 1; end
if(exist('toffset') ~= 1) toffset = 0; end

y = reshape1d(y);

sztsamp = size(tsamp);
tsamp = reshape1d(tsamp) - toffset;

ny = length(y);

isamp = tsamp/dt;
isamp1 = floor(isamp) + 1;
isamp2 = isamp1 + 1;
tsamp1 = dt*(isamp1-1);

n = length(find(isamp1 < 1));
if(n ~= 0)
  fprintf('ERROR: time samples less than zero\n');
  return;
end

n = length(find(isamp1 > ny));
if(n ~= 0)
  fprintf('ERROR: time samples beyond end\n');
keyboard
  
  return;
end

m = (y(isamp2) - y(isamp1))/dt;
ysamp = m .* (tsamp - tsamp1) + y(isamp1);

ysamp = reshape(ysamp,sztsamp);

return;

  
  
  
  
  
