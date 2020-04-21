function [rvt trvt] = fast_rvt(t,resp)
% [rvt trvt] = fast_rvt(t,resp)
%
% Respiratory volume per unit time (RVT) as defined in Birn, et al, NI
% 31 (2006) 1536-1548. Note that the units are arbitrary. Note that
% this is different than the "respiratory volume" measure from
% Change, et al, NI 44 (2009) 857-869.
%

rvt = [];
if(nargin ~= 2)
  fprintf('[rvt trvt] = fast_rvt(t,resp)\n')
  return;
end

indrespP = peakfinder(+resp);
indrespN = peakfinder(-resp);

nP = length(indrespP);
nN = length(indrespN);
nUse = min(nP,nN)-1;

clear rvt trvt;
nth = 1;
for nthP = 1:nUse
  tnP = t(indrespP(nthP));
  dt = t(indrespN) - tnP;
  ind = find(dt > 0);
  if(isempty(ind)) break; end
  nthN = ind(1);
  tnN = t(indrespN(nthN));
  rP = resp(indrespP(nthP));
  rN = resp(indrespN(nthN));
  rvt(nth) = (rP-rN)/(tnN-tnP);
  trvt(nth) = (tnP+tnN)/2;
  nth = nth + 1;
end

return;
