function [X M] = fast_retroicor(t,indpeak,M)
% [X M] = fast_retroicor(t,indpeak,<M>)
%
% Based on Glover, Li, and Ress, Image-Based Method for Retrospective
% Correction of Physiological Motion Effects in fMRI: RETROICOR. MRM,
% 2000, 44:162-167.
%
% $Id: fast_retroicor.m,v 1.1 2008/09/15 21:00:04 greve Exp $

X = [];

if(nargin < 2 | nargin > 3)
  fprintf('X = fast_retroicor(t,indpeak,<M>)\n');
  return;
end

tpeak = t(indpeak);
npeaks = length(indpeak);
nt = length(t);

if(nargin == 2)
  % Compute M
  dt0 = t(2)-t(1);
  idelta = indpeak(2:end)-indpeak(1:end-1);
  M = floor(mean(idelta)/4);
end

X = zeros(nt,2*M);

for nthpeak = 1:npeaks-1
  t1 = tpeak(nthpeak);
  t2 = tpeak(nthpeak+1);
  dt = t2-t1;
  k1 = indpeak(nthpeak);
  k2 = indpeak(nthpeak+1)-1;
  tk = t(k1:k2)-t(k1);
  phik = 2*pi*tk/dt;
  for h = 1:M
    c = 1 + 2*(h-1);
    X(k1:k2,c) = cos(h*phik);
    s = 2 + 2*(h-1);
    X(k1:k2,s) = sin(h*phik);
  end
end

return;



