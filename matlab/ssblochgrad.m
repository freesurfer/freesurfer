function [dSdT2s, dSdT1, dSdPD, dSdTE, dSdFA, dSdTR] = ssblochgrad(tr,te,fa,t1,t2s,pd)
% [dSdT2s dSdT1 dSdPD dSdTE dSdFA dsdTR] = ssblochgrad(tr,te,fa,t1,t2s,<pd>)
%
% Gradients of the steady-state Bloch equation wrt MR and acquisition paramaters.
%
% tr = repetition time
% te = echo time
% fa = flip angle (radians)
% t1 = T1
% t2s = T2 star
% pd = proton density (default = 1)
%
% Time units don't matter as long as they are consitent
%
% if fa=[], set to ernst angle:
%   fa = acos(exp(-TR/T1))
% 
% See also ssbloch
%  
% $Id: ssblochgrad.m,v 1.1 2007/12/23 22:17:21 greve Exp $  

%
% ssblochgrad.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/12/23 22:17:21 $
%    $Revision: 1.1 $
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

s = [];
if(nargin < 5 | nargin > 6)
  fprintf('[dSdT2s dSdT1 dSdPD dSdTE dSdFA dsdTR] = ssblochgrad(tr,te,fa,t1,t2s,<pd>)\n');
  return;
end

if(~exist('pd','var')) pd = 1; end
if(isempty(fa)) fa = acos(exp(-tr./t1)); end % ernst

etrt1 = exp(-tr./t1);
N = pd .* sin(fa) .* (1-etrt1) .* exp(-te./t2s) ; % numerator
D = 1./(1-cos(fa).*etrt1);                        % denominator
s = N.*D;

% -------------------------------------------------------
% dSdT2s
dSdT2s = s.*(-te).*(-1./(t2s.^2));

% dSdT1
dNdT1 = -(pd .* sin(fa) .* exp(-te./t2s) .* (tr.*etrt1) ./ (t1.^2));
dDdT1 = cos(fa) .* (tr.*etrt1) .* (1./(t1.^2)) .* (D.^2);
dSdT1 = (dNdT1.*D + dDdT1.*N);

% dSdPD
dSdPD = s./pd;

% -------------------------------------------------------
% dSdTE
dSdTE = s .* (-1./t2s);

% dSdFA
dNdFA = pd .* cos(fa) .* (1-etrt1) .* exp(-te./t2s);
dDdFA = -etrt1 .* sin(fa) .* (D.^2);
dSdFA  = (dNdFA.*D + dDdFA.*N);

% sSdTR
dNdTR = pd .* sin(fa) .* exp(-te./t2s) .* etrt1 ./ t1;
dDdTR = -(D.^2) .* cos(fa) .* etrt1 ./ t1;
dSdTR = (dNdTR.*D + dDdTR.*N);

% -------------------------------------------------------


return;






