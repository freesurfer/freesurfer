function [phx, phy, epipar] = tdr_epi_phtraj(arg1)
% [phx, phy, epipar] = tdr_epi_phtraj(arg1)
% 
% Assumes phase space trajectory starts in upper right (phx>0 and
% phy>0), then goes to upper left on the first line (phx>0 and
% phy<0), then drops down one line, then goes to the right, etc.
%

phx = [];
phy = [];
epipar = [];
if(nargin ~= 1)
  fprintf('[phx phy epipar] = tdr_epi_phtraj(measasc)\n');
  fprintf('[phx phy epipar] = tdr_epi_phtraj(epipar)\n');
  return;
end

% The 1st arg is either the path to a meas.asc file or an epipar struct
if(isstr(arg1))
  measasc = arg1;
  epipar.tDwell = tdr_measasc(measasc,'sRXSPEC.alDwellTime[0]'); % nsec
  if(isempty(epipar.tDwell)) return; end
  epipar.tDwell = epipar.tDwell/1000; % convert to usec
  epipar.tRampUp   = tdr_measasc(measasc,'alRegridRampupTime[0]');  % us
  epipar.tFlat     = tdr_measasc(measasc,'alRegridFlattopTime[0]'); % us
  epipar.tRampDown = tdr_measasc(measasc,'alRegridRampdownTime[0]'); % us
  epipar.tDelSamp  = tdr_measasc(measasc,'alRegridDelaySamplesTime[0]'); % us
  epipar.eshospacing = tdr_measasc(measasc,'sFastImaging.lEchoSpacing'); % us
  epipar.nkcols    = tdr_measasc(measasc,'iNoOfFourierColumns'); 
  epipar.nkrows    = tdr_measasc(measasc,'iNoOfFourierLines'); 
else
  epipar = arg1;
end


% This is the phase trajectory for a column. It starts at -pi and
% goes to pi-delta, passing thru 0 between (nkcols/2) and
% (nkcols/2)+1. 
phx0 = kspacevector2(epipar.nkcols,epipar.tDwell,epipar.tRampUp,...
		     epipar.tFlat,epipar.tRampDown,epipar.tDelSamp,0);

phx = repmat(phx0,[epipar.nkrows 1]);
% Assume first phx line starts pos and goes neg
phx(1:2:end,:) = fliplr(phx(1:2:end,:)); 

tmp = [0:epipar.nkrows-1]'/epipar.nkrows;
tmp = tmp - tmp(round(epipar.nkrows/2));
% Assume first phy col starts pos and goes neg
phy0 = -2*pi*tmp;
phy = repmat(phy0, [1 epipar.nkcols]);

phx = reshape1d(phx');
phy = reshape1d(phy');


return
