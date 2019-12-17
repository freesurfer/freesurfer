function s = mpr2pseq(TIlist,alpha,Nalpha,TE,AlphaSpacing,RO,t0,s)
% s = mpr2pseq(TIlist,alpha,Nalpha,TE,AlphaSpacing,RO,t0,s)
% Convert MPRAGE parameters into a pulse sequence for a single
%   inversion pulse. These can be concatenated to create a more
%   complete pulse sequence.
% TIlist - list of inversion times (for a typical MPRAGE there is
%   only one inversion)
% alpha - readout flip angle in degrees
% Nalpha - number of readouts
% TE - echo time
% AlphaSpacing - time between readout pulses (sometimes called
%   "echo spacing", but this is confusing)
% RO - flag to indicate whether a readout actually happened. There
%  might not be a readout during the "dummy" scans when trying to 
%  equilibrium
% t0 - temporal offset used when concatenating sequences 
% s - structure used to store the pulse sequence
% 
% This can be converted into readouts with something like
% s = mpr2pseq(TIlist,alpha,Nalpha,TE,AlphaSpacing,RO,1);
% s.T1 = 1400;
% s = blochsimdng(s);
%
% TI=1250;
% alpha = 7;
% Nalpha = 176;
% TE = 1.69;
% AlphaSpacing = 9.8;
% TR = 2530;
% t0 = 0;
% s = mpr2pseq(TI,alpha,Nalpha,TE,AlphaSpacing,0,t0);
% s = mpr2pseq(TI,alpha,Nalpha,TE,AlphaSpacing,1,t0+TR,s);
% [s.tEvent/1000 s.FlipEvent s.IsReadOut]

if(~exist('s','var')) s = []; end
if(isempty(s))
  s.tEvent = [];
  s.FlipEvent = [];
  s.IsReadOut = [];
  nthEv = 0;
else
  nthEv = length(s.tEvent);
end
t = t0;
nthEv = nthEv + 1;

ReadOutTime = Nalpha*AlphaSpacing;

% Inversion
s.tEvent(nthEv,1) = t;
s.FlipEvent(nthEv,1) = 180;
s.IsReadOut(nthEv,1) = 0;
nthEv = nthEv + 1;

% Alphas
nthTI = 0;
for TI = TIlist
  nthTI = nthTI + 1;
  t = t0 + TI - ReadOutTime/2.0;
  for nthAlpha = 1:Nalpha
    s.tEvent(nthEv,1) = t;
    s.FlipEvent(nthEv,1) = alpha(nthTI);
    s.IsReadOut(nthEv,1) = 0;
    nthEv = nthEv + 1;  
    if(RO)
      s.tEvent(nthEv,1) = t+TE;
      s.FlipEvent(nthEv,1) = 0;
      s.IsReadOut(nthEv,1) = nthTI;
      nthEv = nthEv + 1;  
    end
    t = t + AlphaSpacing;
  end
end

return
