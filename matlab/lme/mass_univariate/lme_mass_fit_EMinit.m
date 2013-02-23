function [Theta0,Re0] = lme_mass_fit_EMinit(X,Zcols,Y,ni,maskvtx,nit,prs)
% [Theta0,Re0] = lme_mass_fit_EMinit(X,Zcols,Y,ni,maskvtx,nit,prs)
% 
% Starting values at each location for the linear mixed-effects iterative 
% estimation. These starting values are based on the ordinary least squares
% estimators plus some fast expectation maximization iterations. This
% function should give better estimators than the lme_mass_fit_init
% function but at the cost of extra computational time.
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
% nit: Number of Expectation Maximization iterations. Default 5.
% prs: Number of workers for parallel computing. Default 8;
%
% Output
% Theta0: Matrix whose colums are estimators of the covariance components at.
% each location.
% Re0: Matrix of residual errors at each location.
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
tic;
if nargin < 4
    error('Too few inputs');
elseif nargin < 7
    prs = 8;
    if nargin < 6
        nit = 5;
        if nargin < 5
            maskvtx = [];
        end;
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
nimax = max(ni);
nv = size(Y,2);
nth = q*(q+1)/2+1;
nDth = q*q;
ind = false(nDth,1);
for k=1:q
    for j=1:k
        ind((k-1)*q+j) = true;
    end;
end;
[Theta1,Re1] = lme_mass_fit_init(X,Zcols,Y,ni,[],prs);
display(' ');
display('Starting Expectation Maximization iterations ...');
fname = parfor_progress('init',nv);
parfor j=1:nv
    % Iterations
    W = zeros(n,nimax);
    ths = Theta1(:,j);
    D = zeros(nDth,1);
    D(ind) = ths(1:end-1);
    D = reshape(D,q,q);
    D = D + triu(D,1)';
    phisq = ths(end);
    for k=1:nit 
        %Computation of W = SIGMA^-1 and H.
        posi = 1; H = 0; Term = 0;
        y = Y(:,j);
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
        Re1(:,j) = r;
        Dhat = zeros(q,q);
        phisqhat = 0;
        posi = 1; 
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
            posi = posf+1;
        end;
        D = Dhat/m;
        phisq = phisqhat/n;
    end
    vecD = vec(D);
    Theta1(:,j) = [vecD(ind);phisq];
    if mod(j,2000) == 0
        parfor_progress(fname,2000);
    end;
end
parfor_progress(fname,0);
Theta0 = zeros(nth,nv0);
Theta0(:,maskvtx) = Theta1;
Re0 = zeros(n,nv0);
Re0(:,maskvtx) = Re1;
if (matlabpool('size') > 0)
    matlabpool close;
end;
et = toc;
display(['Total elapsed time is ' num2str(et/60) ' minutes.']);

