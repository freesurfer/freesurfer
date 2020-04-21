function [u, cpvs, pvs, y] = fast_aareg(nf,TR,fmax,fdelta)
% [u, cpvs, pvs, y] = fast_aareg(nf,TR,fmax,<fdelta>)
% anti-aliasing reg

u = [];
pvs = [];
cpvs = [];
y = [];
if(nargin ~= 3 & nargin ~= 4)
  fprintf('[u, cpvs, pvs, y] = fast_aareg(nf,TR,fmax,<fdelta>)\n');
  return;
end

% nf = 100;
% TR = 2;
% fmax = 60/60; % 90 beats per min
% fdelta = .05;

t = TR*[0:nf-1]';
if(nargin == 3) fdelta = []; end
if(isempty(fdelta))
  [fftaxis fdelta] = fast_fftaxis(nf,TR);
end

fTR = 1/2; % sample rate
fnyq = fTR/2; % nyquist frequency
f = [fnyq:fdelta:fmax];
ph = 2*pi*t*f;
y = [cos(ph) sin(ph)];
y = y - repmat(mean(y),[nf 1]); % remove mean
[u s v] = svd(y);
ds2 = diag(s.^2);
pvs = 100*ds2/sum(ds2);
cpvs = cumsum(pvs);

%plot(cpvs,'+-');
return;









