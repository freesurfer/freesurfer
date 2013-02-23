function [stats,st] = lme_FSfit(X,Zcols,y,ni,e)
% [stats,st] = lme_FSfit(X,Zcols,y,ni,e)
%
% Linear mixed-effects estimation by the Fisher scoring algorithm.
%
% Input
% X: Ordered (according to time for each subject) design matrix (nmxp, nm 
% total # of maps, p # of fixed effects). 
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% y: Ordered data vector (according to X).
% ni: Vector whose entries are the number of repeated measures for each
% subject (ordered according to X).
% e: Convergence epsilon (gradient's norm). Default 10^-1;
%
% Output
% stats.Bhat: Estimated vector of the population regresion parameters.
% stats.CovBhat: Estimated covariance matrix of the population regresion 
% parameters.
% stats.phisqhat: Estimated within-subject variability.
% stats.Dhat = Estimated random effects covariance matrix.
% stats.Zcols: Same as Zcols in the input.
% stats.re: Residuals;
% stats.ni: Same as ni in the input.
% stats.invEI: Inverse of the expected information matrix of the restricted
% log-likelihood.
% stats.Pth, stats.Qthth: Matrices that are useful for inferences on the 
% fixed effects.
% stats.lreml: Values of the resticted maximum log-likelihood across the
% optimization process.
% st: Termination state (1 for convergence and 0 otherwise).
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
if nargin < 4 
    error('Too few inputs');   
elseif nargin < 5
    e = 10^-1;
end;
try
    Z = X(:,Zcols);
catch em
    error(['The colums of X specify in Zcols are not correct: ' em.message]);
end
nit = 20;
st = 1;
m = length(ni);
p = size(X,2);
q = length(Zcols);
nth = q*(q+1)/2+1;
ind = [false(q*q,1);true];
for k=1:q
    for j=1:k
        ind((k-1)*q+j) = true;
    end;
end;
n = sum(ni);
W = zeros(n,max(ni));
SIGMA = W;
lr = zeros(1,nit);

%Starting values
[D,phisq] = lme_fit_init(X,Zcols,y,ni);
L = chol(D);
phi = sqrt(phisq);
theta = [vec(L);phi];

%% Iterations
tf = true;
it = 0;
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
    lreml = 0.5*(lreml - log(det(H)));
    gr = lme_Gradient(X,Zcols,W,invH,L,phi,r,ni);
    [EI,Pth,Qthth] = lme_EI(X,Zcols,W,invH,SIGMA,L,phi,ni);
    invEI = EI\eye(nth);
    theta(ind) = theta(ind) + invEI*gr;
    eps = norm(gr);
    lr(it) = lreml;

    %Termination
    if (it==nit) || (eps<e)
        tf = false;
        stats = struct('Bhat',Bhat,'CovBhat',invH,'phisqhat',phisq,'Dhat',D,...
        'Zcols',Zcols,'invEI',invEI,'Pth',Pth,'Qthth',Qthth,'lreml',lr(1:it));
        if it == nit
            st = 0;
        end;
    else
        L = reshape(theta(1:end-1),q,q);
        D = L'*L;
        phi = theta(end);
        phisq = phi*phi;
    end;
end


