function fstats = lme_F(stats,C)
% fstats = lme_F(stats,C) 
% 
% Inference for the fixed-effects in the linear mixed-effects model(Depends 
% on the Statistics toolbox).
%
% Input
% stats: Structure obtained from any of the model fitting functions: 
% lme_fit_EM, lme_fit_FS and lme_fit_NR.
% C: Contrast matrix.
%
% Output
% fstats.F: F-Statistic.
% fstats.pval: P-value of the F-Statistic.
% fstats.sgn: Sign of the contrast.
% fstats.df: Degrees of freedom of the F-Statistic.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:11 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:11 $
%    $Revision: 1.1.2.2 $
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
if nargin < 2
    error('Too few inputs'); 
end;
X = stats.X;
if size(C,2) ~= size(X,2)
     error(['The number of colums in C must be equal to the number of ' ...
                                         'colums in the design matrix X']); 
end;  
ni = stats.ni;
Zcols = stats.Zcols;
q = length(Zcols);
nth = q*(q+1)/2+1;
W = stats.W;
V = stats.SIGMA;
CBhat = stats.CovBhat;
L = chol(stats.Dhat);
phi = sqrt(stats.phisqhat);
%Computation of Pis,Qijs and the expected information matrix EI.
[EI,Pth,Qthth] = lme_EI(X,Zcols,W,CBhat,V,L,phi,ni);
invEI = EI\eye(nth);
%Estimation of the bias in the covariance matrix and computation of the 
%F-statistic and the degrees of freedom of the test.
Bias = 0;
OM = C'*(C*CBhat*C')^-1*C;
A1 = 0; A2 = 0;
Term1 = OM*CBhat;
Term2 = CBhat*OM*CBhat;
for k=1:nth
    Pk = squeeze(Pth(k,:,:));
    for j=1:nth
        Qkj = squeeze(Qthth(k,j,:,:));       
        Pj = squeeze(Pth(j,:,:));
        Bias = Bias + invEI(k,j)*(Qkj-Pk*CBhat*Pj);
        A1 = A1+invEI(k,j)*trace(Term1*Pk*CBhat)*trace(Term1*Pj*CBhat);
        A2 = A2+invEI(k,j)*trace(Term1*Pk*Term2*Pj*CBhat);
    end;
end;
szC = size(C,1);
Bdf = (A1+6*A2)/(2*szC);
g = ((szC+1)*A1-(szC+4)*A2)/((szC+2)*A2);
d = 3*szC+2*(1-g);
c1 = g/d;
c2 = (szC-g)/d;
c3 = (szC+2-g)/d;
EF = (1-A2/szC)^-1;
VF = (2/szC)*(1+c1*Bdf)/((1-c2*Bdf)^2*(1-c3*Bdf));
ro = VF/(2*EF^2);
m = 4+(szC+2)/(szC*ro-1);
l = m/(EF*(m-2));
Bias = CBhat*Bias*CBhat;
CovBhat = CBhat + 2*Bias;
Bhat = stats.Bhat; 
F = l*Bhat'*C'*(C*CovBhat*C')^-1*C*Bhat/szC;
if F<0
    F = 0;
end;
fstats.F = F;
fstats.pval = 1-fcdf(F,szC,m);
fstats.sgn = sign(C*Bhat);
fstats.df = [szC m];
