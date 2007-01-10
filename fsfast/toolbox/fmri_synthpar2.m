function par = fmri_synthpar2(nStimPerCond, TR, Rss, nRuns)
%
% par = fmri_synthpar2(nStimPerCond, TR, Rss, nRuns)
%
%


%
% fmri_synthpar2.m
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

if(nargin < 2 | nargin > 4)
  msg = 'USAGE: Par = fmri_synthpar2(nStimPerCond, TR, <Rss>, <nRuns>)';
  qoe(msg);error(msg);
end

if(nargin == 2 ) Rss = 1; end

if(nargin == 2 | nargin == 3) nRuns = 1; end

Ntp   = sum(nStimPerCond);
nStimPerCond(1) = nStimPerCond(1) + (Rss-1)*Ntp;
Ntp   = sum(nStimPerCond);

nCond = length(nStimPerCond);
TS = TR/Rss;

for r = 1:nRuns
  StimSeq = StimSessionSeq(nStimPerCond, 1);
  par(:,:,r) = [TS*[0:(Ntp-1)]' StimSeq ];
end

return


