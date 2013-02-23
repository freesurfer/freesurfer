function gr = lme_mass_RgGradient(X,Zcols,W,invH,L,phi,re,ni,invG,GDa,GDb)
% gr = lme_mass_RgGradient(X,Zcols,W,invH,L,phi,re,ni,invG,GDa,GDb)
% 
% Gradient vector of the restricted log-likelihood for a whole region.
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% W: Inverses of the estimated temporal covariance matrices for each 
% subject stacked in W.
% invH: Asymptotic covariance matrix of the fixed effects.
% L: Cholesky factor of the covariance matrix of the random effects (D).
% phi: Within-subject standard deviation of the errors.
% re: Residuals;
% ni: Vector whose entries are the number of repeated measures for each
% subject in the study (ordered according to X).
% invG: Inverse of the spatial covariance matrix.
% GDa: Derivative of the spatial covariance matrix for the first spatial 
% parameter.
% GDb: Derivative of the spatial covariance matrix for the second spatial 
% parameter (empty for spatial models with a single parameter).
%
% Output
% gr: Gradient vector.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:10 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:10 $
%    $Revision: 1.1.2.2 $
%
m = length(ni);
q = length(Zcols);
Z = X(:,Zcols);
nth = q*(q+1)/2+1;
nv = size(invG,1);
%Log-likelihoood derivative for L
gr = zeros(nth+1,1);
jk = 0;
for k=1:q
    for j=1:k
        jk = jk + 1;
        posi = 1; posir = 1;
        a1 = 0; a2 = 0; M1 = 0; 
        for i=1:m
            posf = posi+ni(i)-1;
            posfr = posir+ni(i)*nv-1;
            Zi = Z(posi:posf,:);
            Wi = W(posi:posf,1:ni(i));
            Xi = X(posi:posf,:);
            rei = re(posir:posfr);
            Mjki = Zi(:,k)*L(j,:)*Zi';
            Mjki = Mjki + Mjki';
            Maux = Mjki*Wi;
            a1 = a1 - trace(Maux);
            a2 = a2 + rei'*kron(invG,Wi*Maux)*rei;
            M1 = M1 + Xi'*Wi*Maux*Xi;
            posi = posf+1;
            posir = posfr+1;
        end;
        a3 = -trace(invH*M1);
        gr(jk) = (nv*a1+a2-nv*a3)/2;
    end;
end;
%Log-likelihoood derivative for phi
posi = 1; posir = 1;
a1 = 0; a2 =0; M1 = 0;
for i=1:m
    posf = posi+ni(i)-1;
    posfr = posir+ni(i)*nv-1;
    Wi = W(posi:posf,1:ni(i));
    Xi = X(posi:posf,:);
    rei = re(posir:posfr);
    a1 = a1 - trace(Wi);
    Maux = Wi*Wi;
    a2 = a2 + rei'*kron(invG,Maux)*rei;
    M1 = M1 + Xi'*Maux*Xi;
    posi = posf+1;
    posir = posfr+1;
end;
a3 = -trace(invH*M1);
gr(nth) = phi*(nv*a1+a2-nv*a3);
%Log-likelihoood derivative for a
posi = 1; posir = 1;
a2 =0; 
Maux1 = GDa*invG;
Maux2 = invG*Maux1;
for i=1:m
    posf = posi+ni(i)-1;
    posfr = posir+ni(i)*nv-1;
    Wi = W(posi:posf,1:ni(i));
    rei = re(posir:posfr);
    a2 = a2 + rei'*kron(Maux2,Wi)*rei;
    posi = posf+1;
    posir = posfr+1;
end;
aux = -trace(Maux1);
a1 = sum(ni)*aux;
a3 = size(X,2)*aux;
gr(nth+1) = (a1+a2-a3)/2;
if ~isempty(GDb)
    %Log-likelihoood derivative for b
    posi = 1; posir = 1;
    a2 =0;
    Maux1 = GDb*invG;
    Maux2 = invG*Maux1;
    for i=1:m
        posf = posi+ni(i)-1;
        posfr = posir+ni(i)*nv-1;
        Wi = W(posi:posf,1:ni(i));
        rei = re(posir:posfr);
        a2 = a2 + rei'*kron(Maux2,Wi)*rei;
        posi = posf+1;
        posir = posfr+1;
    end;
    aux = -trace(Maux1);
    a1 = sum(ni)*aux;
    a3 = size(X,2)*aux;
    gr(nth+2) = (a1+a2-a3)/2;
else
    gr = gr(1:nth+1);
end;



