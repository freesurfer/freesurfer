function [G, Ip] = computeUpsampled(I, parms)

Iu = upsample(I',parms.dt);
Iu = conv([0 0 Iu(1:end-2)], [1 2 3 2 1]*parms.dt/9, 'same')';
parms.dt = 1;
parms.time_series = 1:parms.T2;
[G, Ip] = compute_BG_from_insulin(Iu, parms);
