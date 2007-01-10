function Xexptrend = fast_exptrendmtx(run,ntrs,nruns,ntau)
% Xexptrend = fast_exptrendmtx(run,ntrs,nruns,ntau)
%
% Exponential trend, orthogonalized wrt mean and linear trend
% 
% x =  (1-exp(-nthtr/ntau)) [then orthoginalized]
%
% ntau is the decay constant measured in TRs.


%
% fast_exptrendmtx.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

if(nargin ~= 4)
  msg = 'USAGE: Xexptrend = fast_exptrendmtx(run,ntrs,nruns,ntau)';
  qoe(msg);error(msg);
end

m = ones(ntrs,1);
t = [0:ntrs-1]'; %'

et = 1 - exp(-t/ntau); 

% Orthoganalize %
f = [m t];
fnot = eye(ntrs) - f*inv(f'*f)*f';
et = fnot*et;

% Normalize magnitude %
et = et./sqrt(sum(et.^2));

Xexptrend        = zeros(ntrs,nruns);
Xexptrend(:,run) = et;

return;
