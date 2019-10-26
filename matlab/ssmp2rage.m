function [y1 y2 mzss] = ssmp2rage(T1,MTR,TR,a1,a2,n,TA,TB,TC,TI1,TI2,inveff)
% [y1 y2 mzss] = mp2ragess(T1,MTR,TR,a1,a2,n,TA,TB,TC,TI1,TI2,inveff)
%   MTR = MRRAGE TR (eg, 5000ms).  
%   TR = time between alpha pulses (eg, % 8ms) 
%
% 4/5/2017 Based on Marques NI 2010. In some cases could get it to
% agree very well with Bloch sims. mzss was always positive whereas
% the actual simulation was negative because of the inversion. It
% appears that using n/2 instead of n/2-1 matches better (see
% below). If number of alphas (n) is odd, then need to interpolate to
% get exactly at the TI. mzss, y1 and y2 were accurate under the
% conditions of Fig 3. I was able to reproduce Fig 3 well, except that
% I had to use TI2=2200, not 3200 (2200 was found to be optimal
% indicated in the text above the fig).  The TRs in Fig3 disagree with
% the optimal, but they don't make much difference. One thing that is
% different is that in the upper left hand corner, the curve does not
% asymptot at -0.5 for high T1 values; instead it curves back around
% (nonunique). It may be that they just did not include that part in
% the plot or that there is still something wrong with my
% implementation. I think this effect is evident in Fig 6. Might be
% present in the blochsim too. The blochsim and theoretical equations
% had differences when the T1 was very low or the flip angle was high
% (>20ish). At low T1s (eg, <400), there are numeric instabilities,
% eg, E1=exp(-TC/T1) can become very close to 0 (eg, 10e-62), so
% changed code to use vpa(); this makes the range better, but it is
% still fairly unstable at low values (eg, T1=80).

E1 = exp(-TR./T1);

% These values can become very small so use variable precision
% It slows things down
%EA = vpa(exp(-TA./T1));
%EB = vpa(exp(-TB./T1));
%EC = vpa(exp(-TC./T1));
EA = (exp(-TA./T1));
EB = (exp(-TB./T1));
EC = (exp(-TC./T1));

ca1 = cos(a1);
ca2 = cos(a2);

ca1n = (ca1.*E1).^n;
ca2n = (ca2.*E1).^n;

f1 = (1-E1).*(1-ca1n)./(1-ca1.*E1);
f2 = (1-E1).*(1-ca2n)./(1-ca2.*E1);

g1 = (1-EA).*ca1n+f1;
g2 =         ca2n+f2;

mzd = 1 + inveff.*((cos(a1).*cos(a2)).^n).*exp(-MTR./T1);

mzss = (( (g1.*EB + (1-EB)).*ca2n + f2).*EC+(1-EC))./mzd;

ca1n2 = (ca1.*E1).^(n/2); % paper calls for n/2-1, but n/2 works better
ca2n2 = (ca2.*E1).^(n/2);
ca2n2b = (ca2.*E1).^(-n/2);

h1 = (-inveff.*mzss.*EA + (1-EA));
y1 = sin(a1).*(h1.*ca1n2 + (1-E1).*(1-ca1n2)./(1-ca1.*E1));

h2 = (mzss-(1-EC))./(EC.*ca2n2);
y2 = sin(a2).*(h2 - (1-E1).*(ca2n2b-1)./(1-ca2.*E1));


return