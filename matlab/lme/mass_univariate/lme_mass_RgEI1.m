function [EI,Pth,Qthth] = lme_mass_RgEI1(X,Zcols,W,CBhat,L,phi,ni,G,GDa,GDb)
% [EI,Pth,Qthth] = lme_mass_RgEI1(X,Zcols,W,CBhat,L,phi,ni,G,GDa,GDb)
%
% Expected information matrix of the restricted log-likelihood of a whole
% region. This is less computationally efficient than lme_mass_RgEI 
% but more numerically stable.
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% W: Inverses of the estimated temporal covariance matrices for each 
% subject stacked in W.
% CBhat: Asymptotic covariance matrix of the fixed effects.
% L: Cholesky factor of the covariance matrix of the random effects (D).
% phi: Within-subject standard deviation of the errors.
% ni: Vector whose entries are the number of repeated measures for each
% subject in the study (ordered according to X).
% G: Spatial covariance matrix.
% GDa: Derivative of the spatial covariance matrix for the first spatial 
% parameter.
% GDb: Derivative of the spatial covariance matrix for the second spatial 
% parameter (empty for spatial models with a single parameter).
%
% Output
% EI: Expected information matrix.
% Pth,Qthth: Matrices that are useful for inferences on the fixed effects.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:10 $
% Original Author: Jorge Luis Bernal Rusiel
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:10 $
%    $Revision: 1.1.2.2 $
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
m = length(ni);
n = sum(ni);
nv = size(G,1);
q = length(Zcols);
nth = q*(q+1)/2+1;
p = size(X,2);
Pth = zeros(nth,p,p);
Qthth = zeros(nth,nth,p,p);
Der = zeros(nth,n,max(ni));
%Computation of the first order derivatives of the temporal covariance matrix
jk = 0;
for k=1:q
    for j=1:k
        jk = jk + 1;
        posi = 1;
        for i=1:m
            posf = posi+ni(i)-1;
            Zi = X(posi:posf, Zcols);
            Zki = Zi(:,k);        
            Mjki = Zki*L(j,:)*Zi';
            Mjki = Mjki + Mjki';          
            Der(jk,posi:posf,1:ni(i)) = Mjki;
            posi = posf+1;
        end;
    end;
end;
posi = 1; 
for i=1:m
    posf = posi+ni(i)-1;
    Der(nth,posi:posf,1:ni(i)) = 2*phi*eye(ni(i));
    posi = posf+1;
end;
%Computation of Pis,Qijs and the expected information matrix EI.
for j=1:nth
    posi = 1; Pj = 0;
    Bj = squeeze(Der(j,:,:));
    for i=1:m
        posf = posi+ni(i)-1;
        Wi = W(posi:posf,1:ni(i));
        Pj = Pj - X(posi:posf,:)'*Wi*Bj(posi:posf,1:ni(i))*Wi*X(posi:posf,:);
        posi = posf+1;
    end;
    Pth(j,:,:) = Pj;
end;
EI = zeros(nth+2,nth+2);
%Expected information among Lijs (including phi)
for k=1:nth
    Bk = squeeze(Der(k,:,:));
    Pk = squeeze(Pth(k,:,:));
    for j=1:k                                 
        Bj = squeeze(Der(j,:,:));
        Pj = squeeze(Pth(j,:,:));
        posi = 1; Qkj = 0;
        traceBkj = 0; 
        for i=1:m
            posf = posi+ni(i)-1;
            Wi = W(posi:posf,1:ni(i));
            Bkji = Wi*Bk(posi:posf,1:ni(i))*Wi*Bj(posi:posf,1:ni(i));
            traceBkj = traceBkj + trace(Bkji); 
            Bkji = Bkji*Wi;
            Qkji = X(posi:posf,:)'*Bkji*X(posi:posf,:);
            Qkj = Qkj + Qkji;
            posi = posf+1;
        end;
        Qthth(k,j,:,:) = Qkj;
        Qthth(j,k,:,:) = Qkj;
        EI(k,j) = nv*(traceBkj - trace(CBhat*(2*Qkj-Pk*CBhat*Pj)));
        EI(j,k) = EI(k,j);
    end;
end;
%Expected information between Lijs (including phi) and a (first spatial parameter)
Mauxa = G\GDa;
trMauxa = trace(Mauxa);
for k=1:nth
    Bk = squeeze(Der(k,:,:));
    Pk = squeeze(Pth(k,:,:));
    posi = 1; 
    traceBk = 0;
    for i=1:m
        posf = posi+ni(i)-1;
        Wi = W(posi:posf,1:ni(i));
        Bki = Wi*Bk(posi:posf,1:ni(i));
        traceBk = traceBk + trace(Bki); 
        posi = posf+1;
    end;
    EI(k,nth+1) = trMauxa*(traceBk + trace(CBhat*Pk));
    EI(nth+1,k) = EI(k,nth+1);
end
%Expected information between a and a
EI(nth+1,nth+1) = (n-p)*trace(Mauxa*Mauxa);
if ~isempty(GDb)
    %Expected information between Lijs (including phi) and b (second spatial parameter)
    Mauxb = G\GDb;
    trMauxb = trace(Mauxb);
    for k=1:nth
        Bk = squeeze(Der(k,:,:));
        Pk = squeeze(Pth(k,:,:));
        posi = 1;
        traceBk = 0;
        for i=1:m
            posf = posi+ni(i)-1;
            Wi = W(posi:posf,1:ni(i));
            Bki = Wi*Bk(posi:posf,1:ni(i));
            traceBk = traceBk + trace(Bki);
            posi = posf+1;
        end;
        EI(k,nth+2) = trMauxb*(traceBk + trace(CBhat*Pk));
        EI(nth+2,k) = EI(k,nth+2);
    end
    %Expected information between b and b
    EI(nth+2,nth+2) = (n-p)*trace(Mauxb*Mauxb);
    %Expected information between a and b
    EI(nth+1,nth+2) = (n-p)*trace(Mauxa*Mauxb);
    EI(nth+2,nth+1) = EI(nth+1,nth+2);
else
    EI = EI(1:nth+1,1:nth+1);
end;
EI = 0.5*EI;


