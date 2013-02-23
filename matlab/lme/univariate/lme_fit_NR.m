function [stats,st] = lme_fit_NR(X,Zcols,y,ni,e)
% [stats,st] = lme_fit_NR(X,Zcols,y,ni,e)
%
% Linear mixed-effects estimation by Newton-Raphson optimization.
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% y: Ordered data vector (according to X).
% ni: Vector whose entries are the number of repeated measures for each
% subject (ordered according to X).
% e: Convergence epsilon (gradient's norm). Default 10^-3;
%
% Output
% stats.Bhat: Estimated vector of the population regresion parameters.
% stats.CovBhat: Estimated covariance matrix of the population regresion 
% parameters.
% stats.bihat: Estimated subject especific random effects. 
% stats.Covbihat: Estimated covariance of the subject especific random 
% coefficients.
% stats.phisqhat: Estimated within-subject variability.
% stats.SIGMA: Estimated marginal covariance matrices for each subject 
% stacked in SIGMA. 
% stats.W: Inverses of the estimated marginal covariance matrices for each 
% subject stacked in W.
% stats.Dhat = Estimated random effects covariance matrix.
% stats.X: Design matrix.
% stats.Zcols: Same as Zcols in the input.
% stats.re: Residuals;
% stats.ni: Same as ni in the input.
% st: Termination state (1 for convergence and 0 otherwise).
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
tic;
if nargin < 4 
    error('Too few inputs');   
elseif nargin < 5
    e = 10^-3;
end;
try
    Z = X(:,Zcols);
catch Me
    error(['The colums of X specify in Zcols are not correct: ' Me.message]);
end
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
nimax = max(ni);
W = zeros(n,nimax);

%Starting values
%[D,phisq] = lme_fit_EMinit(X,Zcols,y,ni,0.1);
[D,phisq] = lme_fit_FSinit(X,Zcols,y,ni,1);
L = chol(D);
phi = sqrt(phisq);
theta = [vec(L);phi];

%% Iterations
tf = true;
it = 0;
display('Starting Newton-Raphson iterations');
while tf 
    it = it+1;
    %Computation of W = SIGMA^-1 and H.
    posi = 1; H = 0; Term = 0;
    invL = L\eye(q);
    scInvD = invL*invL'*phisq;
    for i=1:m
        posf = posi+ni(i)-1;
        Zi = Z(posi:posf,:);
        Wi = (eye(ni(i))-Zi/(Zi'*Zi+scInvD)*Zi')/phisq;
        W(posi:posf,1:ni(i)) = Wi;
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
    He = lme_Hessian(X,Zcols,W,invH,L,phi,r,ni);
    theta(ind) = theta(ind) - He\gr;
    
    %Restricted log-likelihood
    lreml = 0.5*(lreml - log(det(H)));
    display(['Likelihood at NR iteration ' num2str(it) ' : ' num2str(lreml)]);
    eps = norm(gr);
    display(['Gradient norm: ' num2str(eps)]);
    
    %Termination
    if (it==nit) || (eps<e)
        tf = false; 
        D = L'*L;
        SIGMA = zeros(n,nimax);
        Cbihat = zeros(q*m,q);
        bihat = zeros(q,m);
        posi = 1; 
        for i=1:m
            posf = posi+ni(i)-1;
            Zi = Z(posi:posf,:);
            SIGMA(posi:posf,1:ni(i)) = Zi*D*Zi'+ eye(ni(i))*phisq;
            Wi = W(posi:posf,1:ni(i));
            Xi = X(posi:posf,:);
            Pi = Wi - Wi*Xi*invH*Xi'*Wi;
            ri = r(posi:posf);
            bihat(:,i) = D*Zi'*Wi*ri;
            Cbihat((i-1)*q+1:i*q,:) = D-D*Zi'*Pi*Zi*D;
            posi = posf+1;
        end;
        stats = struct('Bhat',Bhat,'CovBhat',invH,'bihat',bihat,...
             'Covbihat',Cbihat,'phisqhat',phisq,'SIGMA',SIGMA,'W',W,...
             'Dhat',D,'X',X,'Zcols',Zcols,'re',r,'ni',ni,'lreml',lreml);
        if it == nit
            st = 0;
            display(['Algorithm does not converge after ' num2str(nit)...
                                                        ' iterations!!!']);
        end;
    else
        L = reshape(theta(1:end-1),q,q);
        phi = theta(end);
        phisq = phi*phi;
    end;
end
et = toc;
display(['Total elapsed time is ' num2str(et) ' seconds']);



