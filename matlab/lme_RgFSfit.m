function [stats,st,a] = lme_RgFSfit(X,Zcols,Y,ni,Theta0,Dist,model,e)
% [stats,st,a] = lme_RgFSfit(X,Zcols,Y,ni,Theta0,Dist,model,e)
%
% Linear mixed-effects estimation by the Fisher scoring algorithm for a 
% whole region. This function is intended to only be called from other 
% functions to perform region-wise mass-univariate analyses (in opposition 
% to lme_mass_RgFSfit).
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% Y: Ordered data Matrix (n x nv, n=#total number of scans and nv=# of 
% vertices/voxels).
% ni: Vector whose entries are the number of repeated measures for each
% subject (ordered according to X) .
% Theta0: Matrix whose colums are initial estimators of the covariance 
% components at each location.
% Dist: Square symmetric matrix with distances among all locations inside
% the region.
% model: Spatial model for the covariance matrix of the region. It can be
% 'exp' or 'gauss'. Default 'exp'.
% e: Convergence epsilon (gradient's norm). Default 10^-1;
%
% Output
% stats: Structure array containing statistiscs for every voxel/vertex 
% inside the region (see lme_FSfit for more details on these statistics).
% st: Array containing the termination state for each location.
% (1 for convergence and 0 otherwise).
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
if nargin < 6
    error('Too few inputs');
elseif nargin < 8
    e = 10^-1;
    if nargin < 7
        model = 'exp';
    end;
end;
try
    Z = X(:,Zcols);
catch em
    error(['The colums of X specify in Zcols are not correct: ' em.message]);
end
m = length(ni);
n = sum(ni);
p = size(X,2);
q = length(Zcols);
nth = q*(q+1)/2+1;
ind = [false(q*q,1);true(2,1)];
for k=1:q
    for j=1:k
        ind((k-1)*q+j) = true;
    end;
end;
W = zeros(n,max(ni));
%Starting values
mTheta0 = mean(Theta0,2);
D = zeros(q*q,1);
D(ind(1:end-2)) = mTheta0(1:end-1);
D = reshape(D,q,q);
D = D + triu(D,1)';
phisq = mTheta0(end);
L = chol(D);
phi = sqrt(phisq);
nv = size(Dist,1);
if strcmpi(model,'gauss')
    Dist = Dist.^2;
    inita = 0.3;
else
    inita = 0.003;
end;
a = inita;
theta = [vec(L);phi;a];

%% Iterations
nit = 20;
lr = zeros(1,nit);
Y = reshapeM(Y,m,ni,nv);
r = Y;
tf = true;
it = 0;
while tf 
    it = it+1;
    %Computation of W = SIGMA^-1 and H.
    posi = 1; H = 0;
    invL = L\eye(q);
    scInvD = invL*invL'*phisq;
    for i=1:m
        posf = posi+ni(i)-1;
        Zi = Z(posi:posf,:);
        Wi = (eye(ni(i))-Zi/(Zi'*Zi+scInvD)*Zi')/phisq;
        W(posi:posf,1:ni(i)) = Wi;
        Xi = X(posi:posf,:);
        H = H + Xi'*Wi*Xi;
        posi = posf+1;
    end;
    invH = H\eye(p);
    %Estimation
    posi = 1; posiY = 1; Bhat = 0;
    for i=1:m
        posf = posi+ni(i)-1;
        posfY = posiY+ni(i)*nv-1;
        Xi = X(posi:posf,:);
        Wi = W(posi:posf,1:ni(i));
        Bhat = Bhat + kron(eye(nv),invH*Xi'*Wi)*Y(posiY:posfY);
        posi = posf+1;
        posiY = posfY+1;
    end;
    %spatial derivative
    [G,GD] = sptDer(Dist,a,model);
    invG = G\eye(nv);
    %log-likelihood and residuals
    posi = 1; posiY = 1;
    lreml = 0; term = 0;
    for i=1:m
        posf = posi+ni(i)-1;
        posfY = posiY+ni(i)*nv-1;
        Xi = X(posi:posf,:);
        Wi = W(posi:posf,1:ni(i));
        ri = Y(posiY:posfY)-kron(eye(nv),Xi)*Bhat;
        r(posiY:posfY) = ri;
        term = term + ri'*kron(invG,Wi)*ri;
        lreml = lreml + log(det(Wi));
        posi = posf+1;
        posiY = posfY+1;
    end;
    lreml = 0.5*((n-p)*log(det(invG)) + nv*(lreml-log(det(H))) - term);
    %Fisher scoring
    gr = lme_mass_RgGradient(X,Zcols,W,invH,L,phi,r,ni,GD,invG);
    [EI,Pth,Qthth] = lme_mass_RgEI(X,Zcols,W,invH,L,phi,ni,GD,invG);
    invEI = EI\eye(nth+1);
    theta(ind) = theta(ind) + invEI*gr;
    eps = norm(gr);
    lr(it) = lreml;
 
    %Termination
    if (it==nit) || (eps<e)
        tf = false;
        if it == nit
            st = zeros(nv,1);
        else
            st = ones(nv,1);
        end;
        D = L'*L;
        invEI = invEI(1:nth,1:nth);
        lr = lr(1:it);
        stats(1:nv) = struct('Bhat',zeros(p,1),'CovBhat',invH,...
            'phisqhat',phisq,'Dhat',D,'Zcols',Zcols,'invEI',invEI,'Pth',...
            Pth,'Qthth',Qthth,'lreml',lr);
        for i=1:nv
            stats(i).Bhat = Bhat((i-1)*p+1:i*p);
        end;
    else
        L = reshape(theta(1:end-2),q,q);
        phi = theta(end-1);
        phisq = phi*phi;
        a = theta(end);
        if a < 0
            a = inita;
        end;
    end;
end
end



%% Auxiliar functions

function RM = reshapeM(Y,m,ni,nv)
posi = 1; posiM = 1;
RM = zeros(sum(ni)*nv,1);
for i=1:m
    posf = posi+ni(i)-1;
    posfM = posiM+ni(i)*nv-1;
    RM(posiM:posfM) = reshape(Y(posi:posf,:),ni(i)*nv,1);
    posi = posf+1;
    posiM = posfM+1;
end;
end


function [G,GD] = sptDer(Dist,a,model)
if strcmpi(model,'gauss')
    G = exp(-a^2.*Dist);
    GD = -2*a.*Dist.*G;
elseif strcmpi(model,'exp')
    G = exp(-a.*Dist);
    GD = -Dist.*G;
else
    error('Valid values for the model parameter are ''gauss'' or ''exp''');
end;
end


