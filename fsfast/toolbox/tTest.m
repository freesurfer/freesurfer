function p = tTest(dof, t, dofApprx)
% 
% pSig = tTestB(dof, t, <dofApprx>)
%
% If dofApprx is not set, tTest() computes the
% significance level of t given dof using the
% formula:
%    z = dof./(dof + t .^ 2);
%    pSig = betainc(z,dof/2,0.5);
% When dofApprx is set, an approximation is used
% when dof exceeds dofApprx:
%    p = erfc(abs(t)/sqrt(2.0));
% This can speed the function considerably when
% the dof is large.
%
% Ref: Numerical Rec in C, pg 229.
%
%
% $Id: tTest.m,v 1.1 2003/03/04 20:47:41 greve Exp $


if(nargin < 2)
  error('Usage: p = tTest(dof, t, <dofApprx>)');
  return;
end

if(length(dof) > 1)
  error('dof must be a scalar');
end

% If dofApprx is unspecified, do not Approx %
if(nargin == 2) dofApprx = -1; end

dofhalf = dof/2;
if(dof < dofApprx | dofApprx < 0)
  z = dof./(dof + t .^ 2);
  p = betainc(z,dof/2,0.5);
else
  p = erfc(abs(t)/sqrt(2.0));
end

return;
