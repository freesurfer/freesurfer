function [order, X] = fast_polyorder(ntp,TR,fCutOffHz)
% [order X] = fast_polyorder(ntp,TR,fCutOffHz)
%
% Computes the order of polynomial regressors needed to achieve a
% certain high-pass cutoff frequency given the number of times points
% and TR. The method was derived empirically.
%
% X = fast_polytrendmtx(1,ntp,1,order);%  
%

%
% fast_polyorder.m
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

order = [];
if(nargin ~= 3)
  fprintf('[order X] = fast_polyorder(ntp,TR,fCutOffHz)\n')';
  return;
end

b1 = [.3654 .3405]';
bexp = b1/(TR*ntp);
% To compute expected cutfoff frequency = bexp(1)+order*bexp(2);
order = round((fCutOffHz-bexp(1))/bexp(2));
if(nargout ==2)
  X = fast_polytrendmtx(1,ntp,1,order);
end

return;

