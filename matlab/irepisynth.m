function irepistruct = irepisynth(irepistruct,M0)
% irepistruct = irepisynth(irepistruct,M0)
% plot(s.tRead,s.yRead)

% yRead = irepisynth(tEvent,FlipEvent,IsReadOut)
% tEvent - time of excitation/readout, same units as T1
% FlipEvent - flip angle at Event, in deg
% IsReadOut - 1 if event is a readout
% T1 - long relaxation, same units as tEvent
% M0 - optional
% sigma - noise stddev optional

if(nargin == 1)
  M0 = 1; % M0 just controls the scale, but snr important below
end
  
s = irepistruct; % copy into s for easy handling

% s.EventSliceNo = -1 for inversions which apply to all slices
indSlice = find(s.EventSliceNo == s.sliceno | s.EventSliceNo < 0);

FlipEventRad = s.eff*s.FlipEvent(indSlice)*pi/180;
tEvSlice = s.tEvent(indSlice);
nEvents = length(tEvSlice);
nT1 = length(s.T1);

s.yMz = zeros(nEvents,nT1);

nthRead = 1;
Mz = M0; 
for n = 1:nEvents
  theta = FlipEventRad(n);
  IsInv = abs(abs(theta)-pi)<.01;
  dt = 0;
  ex0 = 1;
  biex = 1;
  f = 0;
  if(n>1) 
    dt = tEvSlice(n)-tEvSlice(n-1);
    ex = exp(-dt./s.T1);
    ex0 = ex;
    biex = 0;
    if(length(s.biexp)>0) 
      f = s.biexp(1);
      T1biexp = s.biexp(2);
      biex = exp(-dt./T1biexp); 
      ex = ex + f*biex;
    end
    Mz = M0 - (M0-Mz).*ex;
  end
  Mz = Mz*cos(theta);
  if(0 & ~IsInv & n > 0 & n < 34)
  fprintf('%3d %7.1f %3.0f %6.1f %6.4f %6.4f %7.3f\n',...
	  n,tEvSlice(n),theta*180/pi,dt,(1-f)*ex0,f*biex,Mz);
  end
  if(s.IsReadOut(indSlice(n)))
    s.yRead(nthRead,:) = abs(Mz * sin(theta));
    nthRead = nthRead + 1;
  end
  s.yMz(n,:) = Mz;
end

iRead = find(s.IsReadOut(indSlice));
s.tRead = tEvSlice(iRead);

if(0)
  % Barral, MRM, 2010, Robust Methodology ofr invivo T1 Mapping
  % Assumes theta1=180
  % Does not seem to help, makes things a little worse
  indSlice = find(s.EventSliceNo == s.sliceno & s.IsReadOut );
  TSI = s.TI(indSlice);
  TSI2 = repmat(TSI,[1 length(s.T1)]);
  TBI2 = s.TBI * ones(size(TSI2));
  T1b = repmat(s.T1,[s.ntp 1]);
  N = 1+exp(-TBI2./T1b)-2*exp(-TSI2./T1b);
  D = 1+cos(s.ROFlip*pi/180)*exp(-TBI2./T1b);
  s.yRead = abs(sin(s.ROFlip*pi/180).*N./D);
end

if(s.sigma > 0)
  % Add rician noise, pretty slow to use laguerreL() to cache
  snr = s.yRead/s.sigma;
  s.yRead = interp1(s.lagsnr,s.lag,snr)*s.sigma*sqrt(pi/2);
  %yy = laguerreL(0.5,-0.5*(snr.^2))*s.sigma*sqrt(pi/2);
  %plot(yy,s.yRead,'.',yy,yy)
end

irepistruct = s;

return;

