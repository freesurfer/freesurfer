function [Xfir, Nc] = fast_par2Xfir(par,ntrs,TR,TER,TPreStim,TimeWindow,W)
% [Xfir Nc] = fast_par2Xfir(par,ntrs,TR,TER,TPreStim,TimeWindow,W)

Xfir = [];
Nc = [];

if(nargin ~= 7)
fprintf('[Xfir Nc] = fast_par2Xfir(par,ntrs,TR,TER,TPreStim,TimeWindow,W)\n');
return;
end

% Get the number of conditions %
Nc = fast_par2nconds(par);

Xfir = [];
for c = 1: Nc
  indc = find(par(:,2)==c);
  tPres = par(indc,1);
  if(~isempty(W)) Wc = W(indc);
  else Wc = [];
  end
  Xfirc = fast_sched2Xfir(tPres,ntrs,TR,TER,TPreStim,TimeWindow,Wc);
  Xfir = [Xfir Xfirc];
end

return;
