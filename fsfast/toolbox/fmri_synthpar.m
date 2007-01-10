function Par = fmri_synthpar(nStimPerCond, TR, nRuns)
%
% Par = fmri_synthpar(nStimPerCond, TR, nRuns)
%
%
%
%


%
% fmri_synthpar.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: Par = fmri_synthpar(nStimPerCond, TR, <nRuns>)';
  qoe(msg);error(msg);
end

if(nargin == 2) nRuns = 1; end

nCond = length(nStimPerCond);
nTP   = sum(nStimPerCond);

for r = 1:nRuns
  StimSeq = StimSessionSeq(nStimPerCond, 1);
  Par(:,:,r) = [TR*[0:(nTP-1)]' StimSeq ];
end

return


