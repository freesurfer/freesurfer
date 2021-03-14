function [sig, nkeep] = fmri_acorrsig(Ravg,Rstd,pthresh)
%
% [sig, nkeep] = fmri_acorrsig(Ravg,Rstd,pthresh)
%
% Compute the significance of an autocorrelation function
% at each delay given the avg and std at each point. Also
% returns the number of contiguous components after delay
% zero which have signficances below pthresh.
%
%


%
% fmri_acorrsig.m
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

[nR nRuns] = size(Ravg);

ir = [((nR+1)/2 +2):nR];

for r = 1:nRuns,
  tRee = Ravg(ir,r) ./ Rstd(ir,r);
  sig(:,r) = erfc(abs(tRee)/sqrt(2.0));
  nkeep(r) = max(find(sig(:,r)<pthresh));
  if(isempty(nkeep(r))) nkeep(r)=1; end
end

return;

