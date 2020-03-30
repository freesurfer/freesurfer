function [stats,st] = lme_fit_EM(X,Zcols,y,ni,e)
% [stats,st] = lme_fit_EM(X,Zcols,y,ni,e)
%
% Linear mixed-effects estimation by expectation maximization.
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
% stats.ni: Same as ni in the input.;
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
tic;   
if nargin < 4 
    error('Too few inputs');   
elseif nargin < 5
    e = 10^-6;
end;
try
    Z = X(:,Zcols);
catch Me
    error(['The colums of X specify in Zcols are not correct: ' Me.message]);
end
nit = 1000;
st = 1;
m = length(ni);
p = size(X,2);
q = length(Zcols);
n = sum(ni);
nimax = max(ni);
W = zeros(n,nimax);
bihat = zeros(q,m);
theta3 = zeros(q*q+1,3);

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
    Dhat = zeros(q,q);
    phisqhat = 0;
    posi = 1; lreml = 0; 
    for i=1:m
        posf = posi+ni(i)-1;
        Zi = Z(posi:posf,:);
        Wi = W(posi:posf,1:ni(i));
        Xi = X(posi:posf,:);
        ri = r(posi:posf);
        bihat(:,i) = D*Zi'*Wi*ri;
        Pi = Wi - Wi*Xi*invH*Xi'*Wi;
        phisqhat = phisqhat + (ri-Zi*bihat(:,i))'*(ri-Zi*bihat(:,i))+phisq*...
                                             trace(eye(ni(i))-phisq*Pi);
        Dhat = Dhat+bihat(:,i)*bihat(:,i)'+D-D*Zi'*Pi*Zi*D;
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
    %Aitken acceleration
    theta3(:,1) = theta3(:,2);
    theta3(:,2) = theta3(:,3);
    theta3(:,3) = [vec(Dhat);phisqhat];
    if mod(it,10)==0 
        lhat = mean((theta3(:,3)-theta3(:,2))./(theta3(:,2)-theta3(:,1)));
        if (lhat > 0) && (lhat < 1)
           theta3(:,3) = theta3(:,2) + (theta3(:,3)-theta3(:,2))/(1-lhat);
           Dhat = reshape(theta3(1:end-1,3),q,q);
           phisqhat = theta3(end,3);
        end;
    end; 
    
    %Termination
    if (it==nit) || (eps<e)
        tf = false;
        SIGMA = zeros(n,nimax);
        Cbihat = zeros(q*m,q);
        posi = 1; 
        for i=1:m
            posf = posi+ni(i)-1;
            Zi = Z(posi:posf,:);
            SIGMA(posi:posf,1:ni(i)) = Zi*D*Zi'+ eye(ni(i))*phisq;
            Wi = W(posi:posf,1:ni(i));
            Xi = X(posi:posf,:);
            Pi = Wi - Wi*Xi*invH*Xi'*Wi;
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
        lreml0 = lreml;
        D = Dhat;
        phisq = phisqhat;
    end;   
end
toc;
