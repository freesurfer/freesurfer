function [im, tszmos] = volind2mosind(iv, szvol, tszmos)
% [im tszmos] = volind2mosind(iv, szvol, tszmos)

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: [im tszmos] = volind2mosind(iv, szvol, <tszmos>)';
  error(msg);
end

szvol = szvol(1:3);
Nvr = szvol(1);
Nvc = szvol(2);
Nvs = szvol(3);

% Size of Mosaic measured in Tiles %
if(nargin == 2) tszmos = []; end
tszmos = defmossize(Nvs, tszmos);
Ntr = tszmos(1);
Ntc = tszmos(2);

[rv cv sv] = ind2sub(szvol,iv);
[rm cm] = volsub2mossub(rv,cv,sv,szvol,tszmos);

Nmr = Ntr*Nvr;
Nmc = Ntc*Nvc;
szmos = [Nmr Nmc];

im = sub2ind(szmos, rm, cm);

return;
