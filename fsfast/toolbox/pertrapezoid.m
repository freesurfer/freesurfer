function y = pertrapezoid(t,Tru,Tft,Trd,Td)
% y = pertrapezoid(t,Tru,Tft,Trd,<Td>)
% 
% Computes the value of the periodic trapezoid waveform at 
% time t. t can be a vector.
%
% Tru = Ramp Up Duration
% Tft = Flattop Duration
% Trd = Ramp Down Duration
%
% Td is an optional delay.
% 
% See also trapezoid.m
% 
%


%
% pertrapezoid.m
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
if(nargin ~= 4 & nargin ~= 5)
  fprintf('y = pertrapezoid(t,Tru,Tft,Trd,<Td>)\n');
  return;
end

if(exist('Td') ~= 1) Td = []; end
if(isempty(Td)) Td = 0; end

tFTStart = Tru; % Start of Flattop
tRDStart = Tru + Tft; % Start of Ramp Down
tEnd     = Tru + Tft + Trd; % End of Ramp Down

% Subtract delay
t = t - Td;

% Make periodic %
tCycle = 2*tEnd;
t = rem(t,tCycle);

y = zeros(size(t));

% Positive part of the cycle
indpos = find(t < tEnd);
tpos = t(indpos);
y(indpos) = trapezoid(tpos,Tru,Tft,Trd,0);

% Negative part of the cycle
indneg = find(t >= tEnd);
tneg = t(indneg) - tEnd;
y(indneg) = -trapezoid(tneg,Tru,Tft,Trd,0);

return;












