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
% $Id: FTest.m,v 1.1 2003/03/04 20:47:40 greve Exp $
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

if(nargin == 4) dof2 = min(dof2,dof2max); end

z = dof2./(dof2 + dof1 * F);
p = betainc(z, dof2/2, dof1/2);

return;
