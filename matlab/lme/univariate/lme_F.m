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
% Original Author: Jorge Luis Bernal Rusiel 
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
p = size(X,2);
ni = stats.ni;
Zcols = stats.Zcols;
Z = X(:,Zcols);
q = length(Zcols);
nth = q*(q+1)/2+1;
W = stats.W;
V = stats.SIGMA;
CBhat = stats.CovBhat;
L = chol(stats.Dhat);
phi = sqrt(stats.phisqhat);
% Computation of Rijs
Rthth = zeros(nth,nth,p,p);
jk = 0;
for k = 1:q
    for j = 1:k
        jk = jk+1;
        
        uv = 0;
        for v = 1:q
            for u = 1:v
                uv = uv+1;
                
                posi = 1; SumR = 0;
                for i = 1:length(ni)
                    posf = posi+ni(i)-1;
                    
                    Xi = X(posi:posf,:); Zi = Z(posi:posf,:); Wi = W(posi:posf,1:ni(i));
                    Ekj = zeros(q,q); Ekj(k,j) = 1; Euv = zeros(q,q); Euv(u,v) = 1;
                    Ai = Zi*Ekj*Euv*Zi';
                    Ri = Xi'*Wi*(Ai+Ai')*Wi*Xi;
                    SumR = SumR+Ri;
                    
                    posi = posf+1;
                end
                Rthth(jk,uv,:,:) = SumR;
            end
        end
    end
end
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
        Qkj = squeeze(Qthth(k,j,:,:)); Rkj = squeeze(Rthth(k,j,:,:));  
        Pj = squeeze(Pth(j,:,:));
        Bias = Bias + invEI(k,j)*(Qkj-Pk*CBhat*Pj-0.25*Rkj);
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
fstats.pval = max([1-fcdf(F,szC,m), 1e-30]);
fstats.sgn = sign(C*Bhat);
fstats.df = [szC m];
