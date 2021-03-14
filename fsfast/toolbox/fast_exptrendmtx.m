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
