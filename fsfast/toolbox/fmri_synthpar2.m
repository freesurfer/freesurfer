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
%    $Date: 2011/03/02 00:04:06 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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


