function gr = lme_Gradient(X,Zcols,W,invH,L,phi,re,ni)
% gr = lme_Gradient(X,Zcols,W,invH,L,phi,re,ni)
% 
% Gradient vector for the restricted log-likelihood.
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% W: Inverses of the estimated marginal covariance matrices for each 
% subject stacked in W.
% invH: Asymptotic covariance matrix of the fixed effects.
% L: Cholesky factor of the covariance matrix of the random effects (D).
% phi: Within-subject standard deviation of the errors.
% re: Residuals;
% ni: Vector whose entries are the number of repeated measures for each
% subject in the study (ordered according to X).
%
% Output
% gr: Gradient vector.
%
% Original Author: Jorge Luis Bernal Rusiel 
%
m = length(ni);
q = length(Zcols);
Z = X(:,Zcols);
nth = q*(q+1)/2+1;
%Derivatives for L
gr = zeros(nth,1);
jk = 0;
for k=1:q
    for j=1:k
        jk = jk + 1;
        posi = 1;
        a1 = 0; a2 = 0; M1 = 0; 
        for i=1:m
            posf = posi+ni(i)-1;
            Zi = Z(posi:posf,:);
            Wi = W(posi:posf,1:ni(i));
            Xi = X(posi:posf,:);
            rei = re(posi:posf);
            Mjki = Zi(:,k)*L(j,:)*Zi';
            Mjki = Mjki + Mjki';
            Maux = Mjki*Wi;
            a1 = a1 + trace(-Maux);
            a2 = a2 + rei'*Wi*Maux*rei;
            M1 = M1 + Xi'*Wi*Maux*Xi;
            posi = posf+1;
        end;
        a3 = trace(-invH*M1);
        gr(jk) = (a1+a2-a3)/2;
    end;
end;
%Derivative for phi
posi = 1;
a1 = 0; a2 =0; M1 = 0;
for i=1:m
    posf = posi+ni(i)-1;
    Wi = W(posi:posf,1:ni(i));
    Xi = X(posi:posf,:);
    rei = re(posi:posf);
    a1 = a1 + trace(-Wi);
    Maux = Wi*Wi;
    a2 = a2 + rei'*Maux*rei;
    M1 = M1 + Xi'*Maux*Xi;
    posi = posf+1;
end;
a3 = trace(-invH*M1);
gr(nth) = phi*(a1+a2-a3);
