function [sig, nkeep] = fmri_acorrsig(Ravg,Rstd,pthresh)
%
% [sig, nkeep] = fmri_acorrsig(Ravg,Rstd,pthresh)
%
% Compute the significance of an autocorrelation function
% at each delay given the avg and std at each point. Also
% returns the number of contiguous components after delay
% zero which have signficances below pthresh.
%
% $Id: fmri_acorrsig.m,v 1.1 2003/03/04 20:47:39 greve Exp $

[nR nRuns] = size(Ravg);

ir = [((nR+1)/2 +2):nR];

for r = 1:nRuns,
  tRee = Ravg(ir,r) ./ Rstd(ir,r);
  sig(:,r) = erfc(abs(tRee)/sqrt(2.0));
  nkeep(r) = max(find(sig(:,r)<pthresh));
  if(isempty(nkeep(r))) nkeep(r)=1; end
end

return;

