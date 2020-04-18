function pw = lme_realizedPower(stats,C,alpha)                                                                    
% pw = lme_realizedPower(stats,C,alpha) 
%
% Power associated with a contrast C for the actual realization of the 
% linear mixed-effects model (Depends on the Statistics toolbox).
%
% Input
% stats: Structure obtained from any of lme_fit_EM, lme_fit_FS or
% lme_fit_NR functions.
% C: Contrast matrix.
% alpha: Significance level of the test. Default 0.05.
%
% Output
% pw: Power of the test.
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
if nargin < 2
    error('Too few inputs');
elseif nargin < 3
    alpha = 0.05;
end;
X = stats.X;
Zcols = stats.Zcols;
ni = stats.ni;
m = length(ni);
q = length(Zcols);
n = sum(ni);
Z = zeros(n,m*q);
posi = 1;
for i=1:m
    posf = posi+ni(i)-1;
    Z(posi:posf,(i-1)*q+1:i*q) = X(posi:posf,Zcols);
    posi = posf+1;
end;
Bhat = stats.Bhat;
CovBhat = stats.CovBhat;
c = rank(C);
ve = n-rank([X,Z]);
nc = Bhat'*C'*(C*CovBhat*C')^-1*C*Bhat;
cv = finv(1-alpha,c,ve);
pw = 1-ncfcdf(cv,c,ve,nc);
 
