function par = fmri_synthpar2(nStimPerCond, TR, Rss, nRuns)
%
% par = fmri_synthpar2(nStimPerCond, TR, Rss, nRuns)
%
% $Id: fmri_synthpar2.m,v 1.1 2003/03/04 20:47:40 greve Exp $

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


