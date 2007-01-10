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
%


%
% tTest.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:35 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%


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
