function [s,fa] = ssbloch(tr,te,fa,t1,t2s,pd)
% [s fa] = ssbloch(tr,te,fa,t1,t2s,<pd>)
%
% Steady-state Bloch Equation
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
% From: Wansapura, et al, J MAG RESE IMG 9:531 538 (1999)
%  At 3T, 
%  Gray:  T1 = 1331ms, T2* = 42-52 ms
%  White: T1 =  832ms, T2* = 45-48 ms
%  CSF:   T1 = 4163ms                  (Chen Proc SMRI, 2001)
%  Caud:  T1 = 1271ms                  (Chen Proc SMRI, 2001)
%   


%
% ssbloch.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/02/16 19:46:14 $
%    $Revision: 1.4 $
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


% This is something from Gary G. Mildly related.  Using IR, collect
% data with TR >> expected_T1, e.g. 5s, and with TIs of 100 200 500
% 1000 3000 ms.  Fit to
% 
% S(n) = So*(1-2*(1+eps)*exp(-TI(n)/T1) for every voxel,
% 
% where So is proton density, and eps is an error term that accounts
% for inaccuracies in flip angle.  I have a program that gives me the
% three images (T1, So, eps).  eps is usually < 0.05, but depends on
% coil and Bo.


s = [];
if(nargin < 5 | nargin > 6)
  fprintf('[s fa]= ssbloch(tr,te,fa,t1,t2s,<p>)\n');
  return;
end

if(~exist('pd','var')) pd = 1; end
if(isempty(fa)) fa = acos(exp(-tr./t1)); end

etrt1 = exp(-tr./t1);
s = pd .* sin(fa) .* (1-etrt1) .* exp(-te./t2s) ./ (1-cos(fa).*etrt1);

return;






