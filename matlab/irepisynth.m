function irepistruct = irepisynth(irepistruct)
% irepistruct = irepisynth(irepistruct)
% plot(s.tRead,s.yRead)

% yRead = irepisynth(tEvent,FlipEvent,IsReadOut)
% tEvent - time of excitation/readout, same units as T1
% FlipEvent - flip angle at Event, in deg
% IsReadOut - 1 if event is a readout
% T1 - long relaxation, same units as tEvent
% M0 - optional
% sigma - noise stddev optional

s = irepistruct; % copy into s for easy handling

indSlice = find(s.EventSliceNo == s.sliceno | s.EventSliceNo < 0);

FlipEventRad = s.eff*s.FlipEvent(indSlice)*pi/180;
tEvSlice = s.tEvent(indSlice);
nEvents = length(tEvSlice);
nT1 = length(s.T1);

%nRead = length(iRead);

%s.yRead = zeros(nRead,nT1);
s.yMz = zeros(nEvents,nT1);

nthRead = 1;
M0 = 1; % M0 just controls the scale
Mz = M0; 
for n = 1:nEvents
  if(n>1) 
    dt = tEvSlice(n)-tEvSlice(n-1);
    ex = exp(-dt./s.T1);
    %ex = (.9*ex + .1*exp(-dt./(s.T1/10))); % biexponential
    Mz = M0 - (M0-Mz).*ex;
  end
  theta = FlipEventRad(n);
  Mz = Mz*cos(theta);
  if(s.IsReadOut(indSlice(n)))
    s.yRead(nthRead,:) = abs(Mz * sin(theta));
    nthRead = nthRead + 1;
  end
  s.yMz(n,:) = Mz;
end

iRead = find(s.IsReadOut(indSlice));
s.tRead = tEvSlice(iRead);

if(s.sigma > 0)
  % Add rician noise, pretty slow to use laguerreL() to cache
  snr = s.yRead/s.sigma;
  s.yRead = interp1(s.lagsnr,s.lag,snr)*s.sigma*sqrt(pi/2);
  %yy = laguerreL(0.5,-0.5*(snr.^2))*s.sigma*sqrt(pi/2);
  %plot(yy,s.yRead,'.',yy,yy)
end

irepistruct = s;

return;

