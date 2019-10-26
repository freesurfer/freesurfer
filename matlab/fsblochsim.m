function blochstruct = fsblochsim(blochstruct)
% blochstruct = fsblochsim(blochstruct)
% Perform a simulation of the Bloch equations given a pulse
% sequence. To get the structure, see, eg, mpr2pseq(). 

% yRead = irepisynth(tEvent,FlipEvent,IsReadOut)
% tEvent - time of excitation/readout, same units as T1
% FlipEvent - flip angle at Event, in deg
% IsReadOut - 1 if event is a readout
% T1 - long relaxation, same units as tEvent
% M0 - optional
% sigma - noise stddev optional


s = blochstruct; % copy into s for easy handling

FlipEventRad = s.FlipEvent*pi/180;
tEvent = s.tEvent;
nEvents = length(tEvent);
nT1 = length(s.T1);

if(~isfield(s,'Mz')) s.Mz = []; end

nthRead = 1;
M0 = 1;
Mz = M0; 
for n = 1:nEvents
  if(n>1) dt = tEvent(n)-tEvent(n-1);
  else    dt = 0; 
  end
  % Let Mz evolve since the last event
  Mz = M0 - (M0-Mz).*exp(-dt./s.T1);
  theta = FlipEventRad(n);
  if(theta ~= 0) 
    % Get xy magnetization
    Mxy = Mz*sin(theta);
    % Now apply pulse 
    Mz = Mz*cos(theta); 
  end
  if(s.IsReadOut(n))
    s.yRead(nthRead,:) = Mxy; % *Mxy*exp(-dt/s.T2s) 
    nthRead = nthRead + 1;
  end
  s.Mz(end+1,:) = Mz;
end

iRead = find(s.IsReadOut);
s.tRead = tEvent(iRead);

blochstruct = s;

return;

