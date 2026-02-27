function [X M] = fast_retroicor(t,indpeak,M,peakamp,troughamp)
% [X M] = fast_retroicor(t,indpeak,<M>,<peakamp>,<troughamp>)
%
% Based on Glover, Li, and Ress, Image-Based Method for Retrospective
% Correction of Physiological Motion Effects in fMRI: RETROICOR. MRM,
% 2000, 44:162-167.
%

X = [];

if(nargin < 2 | nargin > 5)
  fprintf('[X M] = fast_retroicor(t,indpeak,<M>,<peakamp>,<troughamp>)');
  return;
end

tpeak = t(indpeak);
npeaks = length(indpeak);
nt = length(t);

if(~exist('M','var')) M = []; end
if(isempty(M))
  % Compute M
  dt0 = t(2)-t(1);
  idelta = indpeak(2:end)-indpeak(1:end-1);
  M = floor(mean(idelta)/4);
end
if(~exist('peakamp','var')) peakamp = []; end
if(isempty(peakamp)) peakamp = ones(npeaks,1); end
if(~exist('troughamp','var')) troughamp = []; end
if(isempty(troughamp)) troughamp = -ones(npeaks,1); end

baseline = (peakamp + troughamp)/2;
amp = (peakamp - troughamp)/2;

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
    X(k1:k2,c) = amp(nthpeak)*cos(h*phik) + baseline(nthpeak);
    s = 2 + 2*(h-1);
    X(k1:k2,s) = amp(nthpeak)*sin(h*phik) + baseline(nthpeak);
  end
  if(nthpeak == 1)
    % Time before the 1st peak
    k1 = 1;
    k2 = indpeak(1);
    tk = t(k1:k2)-t(k2);
    phik = 2*pi*tk/dt;
    for h = 1:M
      c = 1 + 2*(h-1);
      X(k1:k2,c) = amp(nthpeak)*cos(h*phik) + baseline(nthpeak);
      s = 2 + 2*(h-1);
      X(k1:k2,s) = amp(nthpeak)*sin(h*phik) + baseline(nthpeak);
    end
  end
  if(nthpeak == npeaks-1)
    % Time after the last peak
    k1 = indpeak(npeaks);
    k2 = nt;
    tk = t(k1:k2)-t(k1);
    phik = 2*pi*tk/dt;
    for h = 1:M
      c = 1 + 2*(h-1);
      X(k1:k2,c) = amp(nthpeak)*cos(h*phik) + baseline(nthpeak);
      s = 2 + 2*(h-1);
      X(k1:k2,s) = amp(nthpeak)*sin(h*phik) + baseline(nthpeak);
    end
  end
  

end

return;



