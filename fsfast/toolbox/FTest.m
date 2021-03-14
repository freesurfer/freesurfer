function p = FTest(dof1, dof2, F, dof2max)
%
% p = FTest(dof1, dof2, F, <dof2max>)
%
% Computes p-value given F-value. p and F can be vectors.
% dof1 = dof of numerator (Number of Rows in RM)
% dof2 = dof of denominator (DOF)
%
% Ref: Numerical Rec in C, pg 229.
%
%
%


%
% FTest.m
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



if(nargin ~= 3 & nargin ~= 4)
  msg = 'Usage: p = FTest(dof1, dof2, F, <dof2max>)';
  qoe(msg);error(msg);
  return;
end

if(length(dof1) > 1)
  error('dof1 must be a scalar');
end
if(length(dof2) > 1)
  error('dof2 must be a scalar');
end

if(exist('dof2max') ~= 1) dof2max = []; end
if(isempty(dof2max))      dof2max = dof2; end
dof2 = min(dof2,dof2max);

z = dof2./(dof2 + dof1 * F);
p = betainc(z, dof2/2, dof1/2);

return;
