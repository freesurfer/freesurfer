function dof = fmri_dof(X,DTOrder)
%
% Compute degrees of freedom.
%
% dof = fmri_dof(X,DTOrder)
%
% X - stim conv matricies: nTP x nTotEst x nRuns, where nTotEst
%     is (nHEst + nPreStim)*nNonNullCond.
% DTOrder - number of trends removed.
%
% The dof is the number of observations (nRuns*nTP) minus
% the number of parameters fit (HOrder*nNNCond + nRuns*PPOrder).
%
%
%


%
% fmri_dof.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

if(nargin ~= 2)
  msg = 'USAGE: dof = fmri_dof(X,DTOrder)';
  qoe(msg);
  error(msg);
end

[nTP nTotEst nRuns] = size(X);

dof = nRuns*nTP - nRuns*DTOrder - nTotEst;

return;
