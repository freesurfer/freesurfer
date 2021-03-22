function wpar = fast_wparnonlin(par,Tpc,tau,attnmax)
%
% wpar = fast_wparnonlin(par,Tpc,tau,attnmax)
%
% Computes a nonlinear weighting for each presentation based 
% on how quickly one stimulus follows another. The nonlinear
% function is based on a gamma function (delta,tau)
%
%  wpar(n) = (1 - attnmax * exp(-dt/tau) );
%
% where attnmax (0,1) is the maximum amount of attenuation, and
% dt is the time between the end of the (n-1)th stimulus and
% the begining of the nth stimulus. If dt=0 (ie, the stimuli 
% abut each other, then the attenuation is max (ie, 
% par = (1-attnmax)). As dt gets larger, the attenuation drops 
% (ie, the weight goes closer to 1). The rate at which the
% attenuation drops is controlled by tau.


%
% fast_wparnonlin.m
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

wpar = [];

if(nargin ~= 4)
  fprintf('USAGE: wpar = fast_wparnonlin(par,Tpc,tau,attnmax)\n');
  return;
end

tPres  = par(:,1);
StimId = par(:,2);
Nstim = length(tPres);

wpar = zeros(Nstim,1);
wpar(1) = 1;
tEndPrev = tPres(1) + Tpc(StimId(1));
for n = 2:Nstim
  dt = tPres(n) - tEndPrev;
  wpar(n) = (1 - attnmax * exp(-dt/tau) );
  tEndPrev = tPres(n) + Tpc(StimId(n));
  % fprintf('%2d %g %g\n',n,dt,wpar(n));
end

return;
