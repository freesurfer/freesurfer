function pd = mp2rage_to_pd(inv1, inv2, t1, TR, TI1, TI2, AlphaSpacing, FA1, FA2, Nalpha)
% pd = mp2rage_to_pd(inv1, inv2, t1, TR, TI1, TI2, AlphaSpacing, FA1, FA2, Nalpha)
% Compute PD from mp2rage INV1, INV2, and T1 outputs. 
% Temporal units must be the same across t1, TR, TI1, TI2, and AlphaSpacing
% inv1 - the first inversion output
% inv2 - the second inversion output
% t1 - quantitative t1 map. 
% TR - time between first inversions (eg, 5000ms)
% TI1, TI2 - inversion times
% AlphaSpacing - time between alpha (readout) pulses (eg, 7ms)
% FA1, FA2 - flip angle of inversion readouts in degrees
% Nalpha - number of alpha pulses for each inversion
% Simulates inv1 and inv2 with pd=1, then fits pd.

TE = 0;

TA = TI1-(AlphaSpacing*Nalpha/2);
TB = TI2-(AlphaSpacing*Nalpha/2)-(TA+AlphaSpacing*Nalpha);
TC = TR-(2*AlphaSpacing*Nalpha+TA+TB);
FA1rad = FA1*pi/180;
FA2rad = FA2*pi/180;

[x1 x2] = ssmp2rage(t1,TR,AlphaSpacing,FA1rad,FA2rad,Nalpha,TA,TB,TC,TI1,TI2,1);
% abs() because scanner does abs()
x1 = 100*abs(x1);
x2 = 100*abs(x2);

pd = (inv1.*x1+inv2.*x2)./(x1.^2 + x2.^2 + eps);
ind = find(isnan(pd));
pd(ind) = 0;

return

