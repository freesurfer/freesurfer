function Par = fmri_synthpar(nStimPerCond, TR, nRuns)
%
% Par = fmri_synthpar(nStimPerCond, TR, nRuns)
%
%
%
% $Id: fmri_synthpar.m,v 1.1 2003/03/04 20:47:40 greve Exp $

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


