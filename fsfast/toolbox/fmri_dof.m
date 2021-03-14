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

if(nargin ~= 2)
  msg = 'USAGE: dof = fmri_dof(X,DTOrder)';
  qoe(msg);
  error(msg);
end

[nTP nTotEst nRuns] = size(X);

dof = nRuns*nTP - nRuns*DTOrder - nTotEst;

return;
