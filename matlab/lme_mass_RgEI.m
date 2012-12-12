function [EI,Pth,Qthth] = lme_mass_RgEI(X,Zcols,W,CBhat,L,phi,ni,GD,invG)
% [EI,Pth,Qthth] = lme_mass_RgEI(X,Zcols,W,CBhat,L,phi,ni,GD,invG)
%
% Expected information matrix of the restricted log-likelihood of a whole
% region. 
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
% GD: Derivative of the spatial covariance matrix.
% invG: Inverse of the spatial covariance matrix.
%
% Output
% EI: Expected information matrix.
% Pth,Qthth: Matrices that are useful for inferences on the fixed effects.
%
% $Revision: 1.2 $  $Date: 2012/12/12 22:58:12 $
% Original Author: Jorge Luis Bernal Rusiel
% CVS Revision Info:
%    $Author: vinke $
%    $Date: 2012/12/12 22:58:12 $
%    $Revision: 1.2 $
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
m = length(ni);
n = sum(ni);
nv = size(invG,1);
q = length(Zcols);
nth = q*(q+1)/2+1;
tnth = nth+1;
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
EI = zeros(tnth,tnth);
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
%Expected information among Lijs (including phi) and a (the spatial parameter)
Maux = GD*invG;
trMaux = trace(Maux);
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
    EI(k,tnth) = trMaux*(traceBk + trace(CBhat*Pk));
    EI(tnth,k) = EI(k,tnth);
end
%Expected information between a and a
EI(tnth,tnth) = (n-p)*trace(Maux*Maux);
EI = 0.5*EI;
