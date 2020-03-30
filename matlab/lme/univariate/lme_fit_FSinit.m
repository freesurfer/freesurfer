function [D0,phisq0,st] = lme_fit_FSinit(X,Zcols,y,ni,e)
% [phisq0,D0,st] = lme_fit_FSinit(X,Zcols,y,ni,e)
%
% Starting values for linear mixed-effects estimation.This function is 
% intended to be used to provide starting values for the lme_fit_NR  
% function by mean of some initial iterations of the Fisher scoring algorithm.
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% y: Ordered data vector (according to X).
% ni: Vector whose entries are the number of repeated measures for each
% subject (ordered according to X).
% e: Convergence epsilon. Default 10^-3;
%
% Output
% phisq0: Estimated within-subject variability.
% D0: Estimated random effects covariance matrix.
% st: Termination state (1 for convergence and 0 otherwise).
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%   
if nargin < 4 
    error('Too few inputs');   
elseif nargin < 5
    e = 10^-3;
end;
nit = 50;
st = 1;
m = length(ni);
p = size(X,2);
q = length(Zcols);
ind = [false(q*q,1);true];
for k=1:q
    for j=1:k
        ind((k-1)*q+j) = true;
    end;
end;
n = sum(ni);
Z = X(:,Zcols);
W = zeros(n,max(ni));
SIGMA = W;

%Starting values
[D,phisq] = lme_fit_init(X,Zcols,y,ni);
L = chol(D);
phi = sqrt(phisq);
theta = [vec(L);phi];

%% Iterations
tf = true;
it = 0;
display('Starting initial Fisher scoring iterations');
while tf
    it = it+1;
    %Computation of W = SIGMA^-1 and H.
    posi = 1; H = 0; Term = 0;
    scInvD = D\eye(q)*phisq;
    for i=1:m
        posf = posi+ni(i)-1;
        Zi = Z(posi:posf,:);
        Wi = (eye(ni(i))-Zi/(Zi'*Zi+scInvD)*Zi')/phisq;
        W(posi:posf,1:ni(i)) = Wi;
        SIGMA(posi:posf,1:ni(i)) = Zi*D*Zi'+ eye(ni(i))*phisq;
        Xi = X(posi:posf,:);
        Ti = Xi'*Wi;
        H = H + Ti*Xi;
        Term = Term + Ti*y(posi:posf);
        posi = posf+1;
    end;
    invH = H\eye(p);
   %Estimation
    Bhat = invH*Term;
    r = y-X*Bhat;
    posi = 1; lreml = 0; 
    for i=1:m
        posf = posi+ni(i)-1;
        Wi = W(posi:posf,1:ni(i));
        ri = r(posi:posf);
        lreml = lreml + log(det(Wi))-ri'*Wi*ri;
        posi = posf+1;
    end;   
    gr = lme_Gradient(X,Zcols,W,invH,L,phi,r,ni);
    EI = lme_EI(X,Zcols,W,invH,SIGMA,L,phi,ni);
    theta(ind) = theta(ind) + EI\gr;
    L = reshape(theta(1:end-1),q,q);
    D = L'*L;
    phi = theta(end);
    phisq = phi*phi;
    
    %Restricted log-likelihood
    lreml = 0.5*(lreml - log(det(H)));
    display(['Likelihood at FS iteration ' num2str(it) ' : ' num2str(lreml)]);
    eps = norm(gr);
    display(['Gradient norm: ' num2str(eps)]);
       
    %Termination
    if (it==nit) || (eps<e)
        tf = false;
        phisq0 = phisq;
        D0 = D;
        if it == nit
            st = 0;
            display(['Initial FS does not converge after ' num2str(nit)...
                                                        ' iterations!!!']);
        end;
    end;
end


