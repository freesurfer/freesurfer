function [kvec, gvec, K, icenter] = kspacevector2(nReadout,tDwell,Tru,Tft,Trd,Tds,Tdelay)
% [kvec, gvec] = kspacevector2(nReadout,tDwell,Tru,Tft,Trd,Tds,<Tdelay>)
%
% Tru = total ramp up time from 0 gradient.
% Tds = time after start of ramp up until sampline starts
% Tdelay = delay waveform
%
% $Id: kspacevector2.m,v 1.2 2003/11/06 19:44:39 greve Exp $


kvec = [];
gvec = [];

if(nargin ~= 6 & nargin ~= 7)
  fprintf('[kvec, gvec] = kspacevector2(nReadout,tDwell,Tru,Tft,Trd,Tds,Tdelay)\n');
  return;
end

if(exist('Tdelay') ~= 1) Tdelay = []; end
if(isempty(Tdelay)) Tdelay = 0; end

% Example values
% nReadout = 128;
% tDwell   = 3.2; % usec
% Tru      = 140;
% Tft      = 190; % usec
% Trd      = 140;
% Tds      =  30;

Ttrap = Tru + Tft + Trd;
Ntrap = round(Ttrap/tDwell);
t = tDwell * [0:Ntrap-1];
gvec0 = trapezoid(t,Tru,Tft,Trd,Tdelay);
kvec0 = cumtrapz(gvec0);
kvec0mid = (kvec0(1)+kvec0(end))/2;
kvec0 = kvec0 - kvec0mid;
kvec0mid = (kvec0(1)+kvec0(end))/2;

gvec00 = trapezoid(t,Tru,Tft,Trd,-Tdelay);
kvec00 = cumtrapz(gvec00);
kvec00mid = (kvec00(1)+kvec00(end))/2;
kvec00 = kvec00 - kvec00mid;
kvec00mid = (kvec00(1)+kvec00(end))/2;

ind0A = max(find(kvec0 < 0));
ind0B = min(find(kvec0 > 0));
kvec0A  = kvec0(ind0A);
kvec0Ad = kvec0mid - kvec0A;
kvec0B  = kvec0(ind0B);
kvec0Bd = - kvec0mid;
icenter0 = ind0A + kvec0Ad/(kvec0B-kvec0A)-1;
%fprintf('icenter0 = %g\n',icenter0)

%indadc = [1:nReadout] + round((Tds+Tdelay)/tDwell);
indadc = [1:nReadout] + round(Tds/tDwell);
gvec = gvec0(indadc);
kvec = kvec0(indadc);

kvec = pi*kvec/abs(kvec(1));
kvecmid = (kvec(1)+kvec(end))/2;
tvec = t(indadc);
[m kcenter] = min(abs(kvec-0));

indA = max(find(kvec < 0));
indB = min(find(kvec > 0));
kvecA  = kvec(indA);
kvecAd = kvecmid - kvecA;
kvecB  = kvec(indB);
kvecBd = kvecB - kvecmid;

icenter = indA + kvecAd/(kvecB-kvecA)-1;
%fprintf('icenter = %g\n',icenter)

rvec = [0:nReadout-1];
%rvec = rvec - rvec(kcenter);
rvec = rvec - icenter;

phmat = (kvec' * rvec); % outer product
K = exp(-i * phmat);

return;

%tadc = tDwell * [0:nReadout-1] + Tds;
%gvec = trapezoid(tadc,Tru,Tft,Trd,Tdelay);

% Integrate gvec to get kvec, rescale
kcenter = nReadout/2 + 1;
kvec = cumsum(gvec);
kvec = kvec - kvec(kcenter);
kvec = pi*kvec/abs(kvec(1));

return;






