function He = lme_Hessian(X,Zcols,W,CBhat,L,phi,re,ni)
% He = lme_Hessian(X,Zcols,W,CBhat,L,phi,re,ni)
% 
% Hessian matrix of the restricted log-likelihood.
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% W: Inverses of the estimated marginal covariance matrices for each 
% subject stacked in W.
% CBhat: Asymptotic covariance matrix of the fixed effects.
% L: Cholesky factor of the covariance matrix of the random effects (D).
% phi: Within-subject standard deviation of the errors.
% re: Residuals;
% ni: Vector whose entries are the number of repeated measures for each
% subject in the study (ordered according to X).
%
% Output
% He:Hessian matrix.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:11 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:11 $
%    $Revision: 1.1.2.2 $
%
m = length(ni);
q = length(Zcols);
phisq = phi*phi;
Z = X(:,Zcols);
nth = q*(q+1)/2+1;
He = zeros(nth,nth);
jk = 0;
for k=1:q
    for j=1:k
        jk = jk + 1;
        nextp = j; pr = jk;
        for r=k:q
            for p=nextp:r
                posi = 1;
                a1 = 0; a2 = 0; M1 = 0; M2 = 0; M3 = 0;
                for i=1:m
                    posf = posi+ni(i)-1;
                    Zi = Z(posi:posf,:);
                    Zki = Z(posi:posf,k);
                    Zri = Z(posi:posf,r);
                    Wi = W(posi:posf,1:ni(i));
                    Xi = X(posi:posf,:);
                    rei = re(posi:posf);
                    Mjki = Zki*L(j,:)*Zi';
                    Mjki = Mjki + Mjki';
                    Mpri = Zri*L(p,:)*Zi';
                    Mpri = Mpri + Mpri';
                    Zkri = Zki*Zri';
                    Zkri = Zkri + Zkri';
                    M1 = M1 + Xi'*Wi*Mpri*Wi*Xi;
                    M2 = M2 + Xi'*Wi*Mjki*Wi*Xi;
                    Maux1 = Mpri*Wi*Mjki;
                    Maux2 = Mjki*Wi*Mpri;
                    if p == j
                        a1 = a1 + trace((Maux2-Zkri)*Wi);
                        Maux3 = Wi*(Maux1-Zkri+Maux2)*Wi;
                        a2 = a2 + rei'*Maux3*rei;
                        M3 = M3 + Xi'*Maux3*Xi;
                    else
                        a1 = a1 + trace(Maux2*Wi);
                        Maux3 = Wi*(Maux1+Maux2)*Wi;
                        a2 = a2 + rei'*Maux3*rei;
                        M3 = M3 + Xi'*Maux3*Xi;
                    end;
                    posi = posf+1;
                end;
                a3 = trace(CBhat*(M3-M1*CBhat*M2));
                He(jk,pr) = (a1-a2-a3)/2;
                He(pr,jk) = He(jk,pr);
                pr = pr + 1;
            end;
            if j == k
                nextp = 1;
            end;
        end;
    end;
end;

jk = 0;
for k=1:q
    for j=1:k
        jk = jk+1;
        posi = 1;
        a1 = 0; a2 = 0; M1 = 0; M2 = 0; M3 = 0;
        for i=1:m
            posf = posi+ni(i)-1;
            Zi = Z(posi:posf,:);
            Zki = Z(posi:posf,k);
            Wi = W(posi:posf,1:ni(i));
            Xi = X(posi:posf,:);
            rei = re(posi:posf);
            Mjki = Zki*L(j,:)*Zi';
            Mjki = Mjki + Mjki';
            a1 = a1 + trace(Mjki*Wi*Wi);
            Maux1 = Wi*(Wi*Mjki+Mjki*Wi)*Wi;
            a2 = a2 + rei'*Maux1*rei;
            M1 = M1 + Xi'*Wi*Wi*Xi;
            M2 = M2 + Xi'*Wi*Mjki*Wi*Xi;
            M3 = M3 + Xi'*Maux1*Xi;
            posi = posf+1;
        end;
        a3 = trace(CBhat*(M3-M1*CBhat*M2));
        He(jk,nth) = phi*(a1-a2-a3);
        He(nth,jk) = He(jk,nth);
    end;  
end;
posi = 1;
a1 = 0; a11 = 0; a2 = 0; a22 = 0;
M1 = 0; M2 = 0;
for i=1:m
    posf = posi+ni(i)-1;
    Wi = W(posi:posf,1:ni(i));
    Xi = X(posi:posf,:);
    rei = re(posi:posf);
    a1 = a1 + trace(Wi*Wi);
    a11 = a11 + trace(Wi);
    a2 = a2 + rei'*Wi*Wi*Wi*rei;
    a22 = a22 + rei'*Wi*Wi*rei;
    M1 = M1 + Xi'*Wi*Wi*Xi;
    M2 = M2 + Xi'*Wi*Wi*Wi*Xi;  
    posi = posf+1;
end;
a1 = 2*phisq*a1-a11;
a2 = 4*phisq*a2-a22;
a3 = 2*phisq*trace(CBhat*(2*M2-M1*CBhat*M1));
He(nth,nth) = a1-a2-a3;
