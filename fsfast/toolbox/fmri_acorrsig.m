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

[nR nRuns] = size(Ravg);

ir = [((nR+1)/2 +2):nR];

for r = 1:nRuns,
  tRee = Ravg(ir,r) ./ Rstd(ir,r);
  sig(:,r) = erfc(abs(tRee)/sqrt(2.0));
  nkeep(r) = max(find(sig(:,r)<pthresh));
  if(isempty(nkeep(r))) nkeep(r)=1; end
end

return;

