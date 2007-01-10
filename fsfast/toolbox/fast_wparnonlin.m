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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
