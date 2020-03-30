function [D0,phisq0,st] = lme_fit_EMinit(X,Zcols,y,ni,e)
% [phisq0,D0,st] = lme_fit_EMinit(X,Zcols,y,ni,e)
%
% Starting values for linear mixed-effects estimation.This function is 
% intended to be used to provide starting values for the lme_fit_NR and 
% lme_fit_FS functions by mean of some initial iterations of the 
% expectation maximization algorithm.
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% y: Ordered data vector (according to X).
% ni: Vector whose entries are the number of repeated measures for each
% subject (ordered according to X).
% e: Convergence epsilon. Default 10^-6;
%
% Output
%
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
    e = 10^-6;
end;
nit = 500;
st = 1;
m = length(ni);
p = size(X,2);
q = length(Zcols);
n = sum(ni);
nimax = max(ni);
W = zeros(n,nimax);

%Starting values
[D,phisq] = lme_fit_init(X,Zcols,y,ni);

%% Iterations
lreml0 = -10^10; 
tf = true;
it = 0;
display('Starting Expectation Maximization iterations');
while tf 
    it = it+1;
    %Computation of W = SIGMA^-1 and H.
    posi = 1; H = 0; Term = 0;
    scInvD = D\eye(q)*phisq;
    for i=1:m
        posf = posi+ni(i)-1;
        Zi = X(posi:posf,Zcols);
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
    Dhat = zeros(q,q);
    phisqhat = 0;
    posi = 1; lreml = 0; 
    for i=1:m
        posf = posi+ni(i)-1;
        Zi = X(posi:posf,Zcols);
        Wi = W(posi:posf,1:ni(i));
        Xi = X(posi:posf,:);
        ri = r(posi:posf);
        bihat = D*Zi'*Wi*ri;
        Pi = Wi - Wi*Xi*invH*Xi'*Wi;
        phisqhat = phisqhat + (ri-Zi*bihat)'*(ri-Zi*bihat)+phisq*...
                                             trace(eye(ni(i))-phisq*Pi);
        Dhat = Dhat+bihat*bihat'+D-D*Zi'*Pi*Zi*D;
        lreml = lreml + log(det(Wi))-ri'*Wi*ri;
        posi = posf+1;
    end;
    phisqhat = phisqhat/n;
    Dhat = Dhat/m;
    %Restricted log-likelihood
    lreml = 0.5*(lreml - log(det(H)));
    eps = abs(lreml-lreml0);
    display(['Likelihood at EM iteration ' num2str(it) ' : ' num2str(lreml)]);
    display(['Epsilon: ' num2str(eps)]);
   
    %Termination
    if (it==nit) || (eps<e)
        tf = false;
        phisq0 = phisqhat; 
        D0 = Dhat;
        if it == nit
            st = 0;
            display(['Algorithm does not converge after ' num2str(nit)...
                                                        ' iterations!!!']);
        end;
    else
        lreml0 = lreml;
        D = Dhat;
        phisq = phisqhat;
    end;   
end

