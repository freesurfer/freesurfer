function [t1, hAvg, hStd, hplt] = hview(hsa,Cond,nHEst,TR,tPre,CtrlCond)
%
% [t, hAvg, hStd, hplt] = hview(hsa,Cond,nHEst,TR,tPre,<CtrlCond>)
%
%
% $Id: hview.m,v 1.1 2003/03/04 20:47:41 greve Exp $


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
