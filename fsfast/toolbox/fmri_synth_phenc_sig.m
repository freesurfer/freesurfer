function signal = fmri_synth_phenc_sig(fund, TR, Ntp, phase, delay, tau)
%
% signal = fmri_synth_phenc_sig(fund, TR, Ntp, phase, delay, tau)
%
% phase in radians
%
%


%
% fmri_synth_phenc_sig.m
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

if(nargin ~= 6)
  msg = 'USAGE: signal = fmri_synth_phenc_sig(fund, TR, Ntp, phase, delay, tau)';
  qoe(msg); error(msg);
end

phase = reshape1d(phase)'; %'
delay = reshape1d(delay)'; %'

if(isempty(tau)) tau = -1; end

if(prod(size(tau)) ~= 1)
  msg = 'tau must be scalar';
  qoe(msg); error(msg);
end

Nv1 = length(phase);
Nv2 = length(delay);
Nv3 = length(tau);

if( (Nv1 ~= 1 & Nv2 ~= 1 & Nv1 ~= Nv2) | ...
    (Nv1 ~= 1 & Nv3 ~= 1 & Nv1 ~= Nv3) | ...
    (Nv2 ~= 1 & Nv3 ~= 1 & Nv2 ~= Nv3) )
  size(phase)
  size(delay)
  size(tau)
  msg = 'phase, delay, and/or tau have inconsitent dimensions'
  qoe(msg); error(msg);
end

Nv = max([Nv1 Nv2 Nv3]);
%fprintf('Nv = %d\n',Nv);

Trun    = TR*Ntp;
Ncycles = Trun*fund;
Tcycle  = 1/fund;
t   = TR*[0:Ntp-1]'; %'

%% Pure sinusoid %%
if(tau < 0) 
  total_phase = phase + delay*(2*pi)/Tcycle;
  arg = 2*pi*t*fund - total_phase;  % use '-' because its a delay
  s = cos(arg);
  signal = repmat(s, [1 Nv]);
  return;
end

%% If it gets this far, use gamma function %%
t2d = repmat(t, [1 Nv]);
phase_delay = phase * Tcycle/(2*pi);
time_to_peak = 2*tau + delay;
time_to_peak = delay;
total_delay =  phase_delay + time_to_peak;
totdel2d = repmat(total_delay, [Ntp 1]);

signal = 0;
for n = 1:Ncycles;
  cycle_delay = (n-1) * Tcycle;
  dt  = (t2d - (cycle_delay + totdel2d))/tau;

  ind = find(dt < 0);
  dt(ind) = 0;

  if(length(dt) > 0)
    s = (dt.^2) .* exp(-dt);
    signal = signal + s;
  end

end

return;
