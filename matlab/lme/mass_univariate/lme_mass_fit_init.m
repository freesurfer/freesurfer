function [Theta0,Re0] = lme_mass_fit_init(X,Zcols,Y,ni,maskvtx,prs)
% [Theta0,Re0] = lme_mass_fit_init(X,Zcols,Y,ni,maskvtx,prs) 
% 
% Starting values at each location for the linear mixed-effects iterative 
% estimation. These starting values are based on the ordinary least squares
% estimators (OLS).
%
% Input
% X: Ordered (according to time for each subject) design matrix (nmxp, nm 
% total # of maps, p # of fixed effects). 
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% Y: Data matrix (nmxnv, nv #vertices) whos colums represent the ordered 
% data vector (according to X) for each voxel/vertex.
% ni: Vector whose entries are the number of repeated measures for each
% subject (ordered according to X).
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
% prs: Number of workers for parallel computing. Default numcores;
%
% Output
% Theta0: Matrix whose colums are estimators of the covariance components at.
% each location.
% Re0: Matrix of residual errors at each location (from OLS).
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
tic;
if nargin < 4
    error('Too few inputs');
elseif nargin < 6
    prs = feature('numcores');
    if nargin < 5
        maskvtx = [];
    end;
end;
nv0 = size(Y,2);
if isempty(maskvtx)
    maskvtx = 1:nv0;
end;
Y = Y(:,maskvtx);
m = length(ni);
p = size(X,2);
q = length(Zcols);
n = sum(ni);
nv = size(Y,2);
ind = false(q*q,1);
for k=1:q
    for j=1:k
        ind((k-1)*q+j) = true;
    end;
end;
nth = q*(q+1)/2+1;
Theta1 = zeros(nth,nv);
phisqd = (n-(m-1)*q-p);
posi = 1;
t = zeros(m*q,q); 
t1 = 0; 
for i=1:m
    posf = posi+ni(i)-1;
    Zi = X(posi:posf,Zcols);
    t2 = pinv(Zi'*Zi);
    t((i-1)*q+1:i*q,:) = t2;
    t1 = t1 + t2;
    posi = posf+1;
end;
if license('test','distrib_computing_toolbox')
    if verLessThan('matlab','8.2.0.29')
        if (prs==1) || (matlabpool('size') ~= prs)
            if (matlabpool('size') > 0)
                matlabpool close;
            end;
            if (prs>1)
                matlabpool(prs);
            end;
        end;
    else
        pc = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(gcp('nocreate'))
            if (prs>1)
                parpool(prs);
            end;
        elseif (pc.NumWorkers  ~= prs) || (prs==1)
            if  ~isempty(gcp('nocreate'))
                delete(gcp('nocreate'))
            end
            if (prs>1)
                parpool(prs);
            end;
        end;
    end;
else
    display(' ');
    display('Warning: Parallel Computing Toolbox missing, things will be real slow ...');
end;   
display(' ');
display('Computing initial values ...');
sinvX = pinv(X);
Bhat = sinvX*Y;
parfor j=1:nv
    %Estimation
    b = Bhat(:,j);
    y = Y(:,j);
    t2 = 0; t3 = 0; t4 = 0; 
    posi = 1;
    for i=1:m
        posf = posi+ni(i)-1;
        Zi = X(posi:posf,Zcols);
        Xi = X(posi:posf,:);
        yi = y(posi:posf);
        ri = yi-Xi*b;
        bihat = t((i-1)*q+1:i*q,:)*Zi'*ri;
        t2 = t2 + yi'*yi-bihat'*Zi'*ri;
        t3 = t3 + Xi'*yi;
        t4 = t4 + bihat*bihat';
        posi = posf+1;
    end;
    phisq = (t2-b'*t3)/phisqd;
    if phisq <= 0
        phisq = 1;
    end;
    D = (t4-phisq*t1)/m;
    try
        chol(D);
    catch
        % Handling non-positive definite initial D;
        [EV,ED] = eig(D);
        ED(ED<0) = 10^-4;
        D = EV*ED*EV^-1;
    end
    vecD = vec(D);
    Theta1(:,j) = [vecD(ind);phisq];
end
Theta0 = zeros(nth,nv0);
Theta0(:,maskvtx) = Theta1;
Re0 = zeros(n,nv0);
Re0(:,maskvtx) = Y-X*Bhat;
display('Done');
et = toc;
display(['Elapsed time is ' num2str(et/60) ' minutes.']);
