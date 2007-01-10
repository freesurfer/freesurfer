function [t1, hAvg, hStd, hplt] = hview(hsa,Cond,nHEst,TR,tPre,CtrlCond)
%
% [t, hAvg, hStd, hplt] = hview(hsa,Cond,nHEst,TR,tPre,<CtrlCond>)
%
%
%


%
% hview.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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


if(nargin ~= 6) CtrlCond = -1; end

nPreStim = floor(tPre/TR);
nHEst = nHEst + nPreStim;

nCondTot   = length(hsa)/(2*nHEst); % includes fixation
nCondView = length(Cond);

ih = [1:nHEst];
ind = [];
for m = 1:nCondView,
  n = Cond(m)+1;
  ind = [ind (ih+2*nHEst*(n-1))];
end

hAvg = reshape(hsa(ind),[nHEst nCondView]);
hStd = reshape(hsa(ind+nHEst),[nHEst nCondView]);

if(CtrlCond > -1) 
  hAvgCtrl = squeeze(hsa(2*CtrlCond*nHEst + ih));
  hAvg = hAvg - repmat(hAvgCtrl,[nCondView 1]);
end

t1 = TR*[0:nHEst-1]' - tPre;
t = reshape(repmat(t1,[nCondView 1]),nHEst,nCondView);

hplt = errorbar(t,hAvg,hStd,'o-');
co = strvcat('blue','green','red','cyan','mag','yell','blk');

tit = ' ';
for m = 1:nCondView,
  n = Cond(m);
  tit = cat(2,tit,'Cond', int2str(n),':  ',deblank(co(m,:)),', ');
end
hold;
plot(t1,zeros(nHEst,1),'k-.');
hold;

title(tit);
xlabel('Time (sec)');
set(gca,'XLim',[min(t1)-TR max(t1)+TR]);

if(nargout < 4) clear hplt; end
if(nargout < 3) clear hStd; end
if(nargout < 2) clear hAvg; end
if(nargout < 1) clear t; end

return;
