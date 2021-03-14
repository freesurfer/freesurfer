function y = tdr_intpertrap(t,Tru,Tft,Trd)
% y = tdr_intpertrap(t,Tru,Tft,Trd)
%
% Exact integral of a periodic trapezoid waveform at time t.
%   Tru - ramp up time of trapezoid
%   Tft - flat time of trapezoid
%   Tdu - ramp down time of trapezoid
%
% The trapezoid starts at zero, linearly ramps up for a time Tru,
% stays flat for Tft, then linearly ramps back to zero over the
% next Trd, then ramps negatively for Tru, then stays flat for Tft,
% then ramps back to zero over the next Tru. The period is
% Tru+Tft+Trd. The integral is shifted and scaled so that it equals -1
% at t=0, peaks at +1 at t=Tru+Tft+Trd=Tperiod/2. The zero-crossing
% is a little hard to compute when Tru != Trd, othersize the
% zero-crossing will be located midway between the peak and trough.
%
%
%


%
% tdr_intpertrap.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

y = [];
if(nargin ~= 4)
  fprintf('y = tdr_intpertrap(t,Tru,Tft,Trd)\n');
  return;
end

% Period of the waveform
Tper = 2*(Tru+Tft+Trd);

% It's periodic, so force 0 <= t < Tper
t = rem(t,Tper); 

% Max of the waveform prior to shifting and scaling
ymax = Tru/2 + Tft + Trd/2;

y = zeros(size(t));

% Ramp up %
ind = find(t >= 0 & t < Tru);
if(~isempty(ind))
  dt = t(ind);
  y(ind) = (dt.^2)/(2*Tru);
end

% Flat top %
ind = find(t >= Tru & t < Tru+Tft);
if(~isempty(ind))
  dt = t(ind) - Tru;
  y(ind) = dt + Tru/2;
end

% Ramp down %
ind = find(t >= Tru+Tft & t <= Tper/2);
if(~isempty(ind))
  dt = t(ind) - (Tru+Tft);
  y(ind) = -(dt.^2)/(2*Trd) + dt + Tru/2 + Tft;
end

% Have to shift and rescale prior to recursive call %
y = (y-ymax/2)/(ymax/2);

% For times in the secon half of the period %
ind = find(t > Tper/2);
if(~isempty(ind))
  dt = Tper - t(ind);
  y(ind) = tdr_intpertrap(dt,Tru,Tft,Trd);
end


return;
