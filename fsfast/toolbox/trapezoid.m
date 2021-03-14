function y = trapezoid(t,Tru,Tft,Trd,Td)
% y = trapezoid(t,Tru,Tft,Trd,<Td>)
% 
% Computes the value of the trapezoid waveform at time t
% t can be a vector.
%
% Tru = Ramp Up Duration
% Tft = Flattop Duration
% Trd = Ramp Down Duration
%
% Td is an optional delay.
% 
% For times before Td and after the end, the result is 0.
% The amplitude at the flat top is 1.
%
%


%
% trapezoid.m
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
  fprintf('y = trapezoid(t,Tru,Tft,Trd,<Td>)\n');
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
%tCycle = tEnd;
%t = rem(t,tCycle);

y = zeros(size(t));

%ind = find(t > 0 & t < tFTStart);
ind = find(t <= tFTStart);
if(~isempty(ind)) y(ind) = t(ind)/Tru; end

ind = find(t > tFTStart & t < tRDStart);
if(~isempty(ind)) y(ind) = 1; end

%ind = find(t > tRDStart & t < tEnd);
ind = find(t >= tRDStart );
if(~isempty(ind)) y(ind) = -(t(ind)-tEnd)/Trd; end


return;












