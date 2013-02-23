function [EI,Pth,Qthth] = lme_EI(X,Zcols,W,CBhat,SIGMA,L,phi,ni)
% [EI,Pth,Qthth] = lme_EI(X,Zcols,W,CBhat,SIGMA,L,phi,ni)
%
% Expected information matrix for the restricted log-likelihood. 
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% W: Inverses of the estimated marginal covariance matrices for each 
% subject stacked in W.
% CBhat: Asymptotic covariance matrix of the fixed effects.
% SIGMA: Estimated marginal covariance matrices for each subject 
% stacked in SIGMA. 
% L: Cholesky factor of the covariance matrix of the random effects (D).
% phi: Within-subject standard deviation of the errors.
% ni: Vector whose entries are the number of repeated measures for each
% subject in the study (ordered according to X).
%
% Output
% EI: Expected information matrix.
% Pth,Qthth: Matrices that are useful for inferences on the fixed effects.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:11 $
% Original Author: Jorge Luis Bernal Rusiel
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:11 $
%    $Revision: 1.1.2.2 $
% Reference: Kenward MG and Roger JH, 1997. Small sample inference for fixed
% effects from restricted maximum likelihood. Biometrics,Vol. 53, No.3.
%
m = length(ni);
n = sum(ni);
q = length(Zcols);
nth = q*(q+1)/2+1;
p = size(X,2);
Pth = zeros(nth,p,p);
Qthth = zeros(nth,nth,p,p);
Der = zeros(nth,n,size(W,2));
%Computation of the first order derivatives of W
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
            Wi = W(posi:posf,1:ni(i));
            Der(jk,posi:posf,1:ni(i)) = -Wi*Mjki*Wi;
            posi = posf+1;
        end;
    end;
end;
posi = 1;
for i=1:m
    posf = posi+ni(i)-1;
    Wi = W(posi:posf,1:ni(i));
    Der(nth,posi:posf,1:ni(i)) = -2*phi*Wi*Wi;
    posi = posf+1;
end;
%Computation of Pis,Qijs and the expected information matrix EI.
for j=1:nth
    posi = 1; Pj = 0;
    Bj = squeeze(Der(j,:,:));
    for i=1:m
        posf = posi+ni(i)-1;
        Pj = Pj + X(posi:posf,:)'*Bj(posi:posf,1:ni(i))*X(posi:posf,:);
        posi = posf+1;
    end;
    Pth(j,:,:) = Pj;
end;
EI = zeros(nth,nth);
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
            Vi = SIGMA(posi:posf,1:ni(i));
            Bkji = Bk(posi:posf,1:ni(i))*Vi*Bj(posi:posf,1:ni(i));
            traceBkj = traceBkj + trace(Bkji*Vi);
            Qkji = X(posi:posf,:)'*Bkji*X(posi:posf,:);
            Qkj = Qkj + Qkji;
            posi = posf+1;
        end;
        Qthth(k,j,:,:) = Qkj;
        Qthth(j,k,:,:) = Qkj;
        EI(k,j) = 0.5*(traceBkj - trace(CBhat*(2*Qkj-Pk*CBhat*Pj)));
        EI(j,k) = EI(k,j);
    end;
end;
end
