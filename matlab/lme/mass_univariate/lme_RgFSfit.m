function [stats,st,a,b] = lme_RgFSfit(X,Zcols,Y,ni,Theta0,Dist,model,e)
% [stats,st,a,b] = lme_RgFSfit(X,Zcols,Y,ni,Theta0,Dist,model,e)
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
% 'exp' or 'sph' or 'gauss'. Default 'exp'.
% e: Convergence epsilon (gradient's norm). Default 1;
%
% Output
% stats: Structure array containing statistiscs for every voxel/vertex 
% inside the region (see lme_FSfit for more details on these statistics).
% st: Array containing the termination state for each location.
% (1 for convergence and 0 otherwise).
% a: Estimate of the first parameter of the spatial correlation matrix.
% b: Estimate of the second parameter of the spatial correlation matrix (
% empty for spatial models with a single parameter).
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%   
if nargin < 6
    error('Too few inputs');
elseif nargin < 8
    e = 0.1;
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
qsq = q*q;
nth = (qsq+q)/2+1;
ind = [false(qsq,1);true(2,1)];
for k=1:q
    for j=1:k
        ind((k-1)*q+j) = true;
    end;
end;
W = zeros(n,max(ni));
SIGMA = W;
%Starting values
mTheta0 = mean(Theta0,2);
D = zeros(qsq,1);
D(ind(1:end-2)) = mTheta0(1:end-1);
D = reshape(D,q,q);
D = D + triu(D,1)';
phisq = mTheta0(end);
L = chol(D);
phi = sqrt(phisq);
nv = size(Dist,1);
b = [];
if strcmpi(model,'gauss')
    if nv>2
        ind = [ind;true];
        b = 0.01;
    else
        model = 'exp';
        Dist = Dist.*Dist;
    end
end
a = 0.05;
theta = [vec(L);phi;a;b];


%% Iterations
nit = 35;
lr = zeros(1,nit);
Y = reshapeM(Y,m,ni,nv);
r = Y;
tf = true;
it = 0;
while tf 
    it = it+1;
    %Computation of W = SIGMA^-1 and H.
    posi = 1; H = 0;
    scInvD = inv(D)*phisq;
    for i=1:m
        posf = posi+ni(i)-1;
        Zi = Z(posi:posf,:);
        Wi = (eye(ni(i))-Zi/(Zi'*Zi+scInvD)*Zi')/phisq;
        SIGMA(posi:posf,1:ni(i)) = Zi*D*Zi'+ eye(ni(i))*phisq;
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
    [G,GDa,GDb] = sptDer(Dist,a,b,model);
    %log-likelihood and residuals
    posi = 1; posiY = 1;
    lreml = 0; term = 0;
    for i=1:m
        posf = posi+ni(i)-1;
        posfY = posiY+ni(i)*nv-1;
        Xi = X(posi:posf,:);
        Wi = W(posi:posf,1:ni(i));
        SIGMAi = SIGMA(posi:posf,1:ni(i));
        ri = Y(posiY:posfY)-kron(eye(nv),Xi)*Bhat;
        r(posiY:posfY) = ri;
        term = term + ri'*(kron(G,SIGMAi)\ri);
        lreml = lreml + log(det(Wi));
        posi = posf+1;
        posiY = posfY+1;
    end;   
    lreml = 0.5*((n-p)*log(1/det(G)) + nv*(lreml-log(det(H))) - term);
    %Fisher scoring
    if (cond(G) < 1e+10)
        invG = G\eye(nv);
        gr = lme_mass_RgGradient(X,Zcols,W,invH,L,phi,r,ni,invG,GDa,GDb);
        [EI,Pth,Qthth] = lme_mass_RgEI(X,Zcols,W,invH,L,phi,ni,invG,GDa,GDb);
    else
        gr = lme_mass_RgGradient1(X,Zcols,SIGMA,W,invH,L,phi,r,ni,G,GDa,GDb);
        [EI,Pth,Qthth] = lme_mass_RgEI1(X,Zcols,W,invH,L,phi,ni,G,GDa,GDb);
    end
    if (cond(EI) < 1e+10)
        invEI = EI\eye(size(EI,1));
    else
        [Vtemp,Dtemp] = eig(EI);
        invEI = Vtemp*diag(1./max(diag(Dtemp),1e-5))*Vtemp';
    end
    theta(ind) = theta(ind) + invEI*gr;
    gnorm = norm(gr);
    lr(it) = lreml;
    
    %Termination
    if (it==nit) || (gnorm<e) 
        tf = false;
        if gnorm < e
            st = ones(nv,1);
        else
            st = zeros(nv,1);
        end;
        invEI = invEI(1:nth,1:nth);
        stats(1:nv) = struct('Bhat',zeros(p,1),'CovBhat',invH,...
            'phisqhat',phisq,'Dhat',D,'Zcols',Zcols,'invEI',invEI,'Pth',...
            Pth,'Qthth',Qthth,'lreml',lr(1:it));
        for i=1:nv
            stats(i).Bhat = Bhat((i-1)*p+1:i*p);
        end;
    else
        L = reshape(theta(1:qsq),q,q);
        D = L'*L;
        phi = theta(qsq+1);
        phisq = phi*phi;
        a = theta(qsq+2);
        if strcmpi(model,'gauss')
            b = theta(end);
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


function [G,GDa,GDb] = sptDer(Dist,a,b,model)
if strcmpi(model,'gauss')
    Distsq = Dist.*Dist;
    G = exp(-b*b*Dist-a*a*Distsq);
    try
        chol(G);
    catch
        [Vtemp,Dtemp] = eig(G);
        G = Vtemp*diag(max(diag(Dtemp),1e-5))*Vtemp';
    end;
    GDa = -2*a*Distsq.*G;
    GDb = -2*b*Dist.*G;
elseif strcmpi(model,'exp')
    G = exp(-a*a*Dist);
    try
        chol(G);
    catch
        [Vtemp,Dtemp] = eig(G);
        G = Vtemp*diag(max(diag(Dtemp),1e-5))*Vtemp';
    end;
    GDa = -2*a*Dist.*G;
    GDb = [];
elseif strcmpi(model,'sph')
    C = a*a*Dist;
    G = C; Ind = C<=1; 
    C3 = C(Ind).^3;
    G(Ind) = 1-0.5*(3*C(Ind)-C3);
    G(C>1) = 0;
    try
        chol(G);
    catch
        [Vtemp,Dtemp] = eig(G);
        G = Vtemp*diag(max(diag(Dtemp),1e-5))*Vtemp';
    end;
    GDa = G;
    GDa(Ind) = -3*(a*Dist(Ind)-C3./a);  
    GDb = [];
else
    error('Valid spatial parametric models are ''exp'',''gauss'', or ''sph''');
end;
end

