function [stats,st] = SStat_mass_CoxPH(X,d,t,Y,maskvtx,prs,e)
% [stats,st] = SStat_mass_CoxPH(X,d,t,Y,maskvtx,prs,e)
%
% Vertex/voxel-wise parameter estimation for the Cox proportional hazards 
% model. Survival time ties are handled using Efron's method. 
%
% Input
% X: Design Matrix with the time-independent covariates. (mxp, m # of
% subjects, p # of covariates). 
% d: Logical vector (mx1) indicating censorship status (1 if the subject got 
% the failure event or 0 otherwise).
% t: Vector (mx1) whose entries are the survival and censored times (ordered 
% according to X).
% Y: Data matrix (mxnv, nv #vertices) whos colums represent the ordered 
% data vector (according to X) for each voxel/vertex.
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
% prs: Number of workers for parallel computing. Default 8;
% e: Convergence epsilon (gradient's norm). Default 10^-2;
%
% Output
% stats: Structure array containing statistiscs for every voxel/vertex (see
% SStat_CoxPH for more details on these statistics).
% st: Array containing the termination state for each voxel/vertex
% (1 for convergence and 0 otherwise).
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Kleinbaum, D.G., Klein, M., 2005. Survival analysis. A self-
% learning approach, second edition. New York: Springer..
%   
tic;
%Validation of inputs
if nargin < 4
    error('Too few inputs');
elseif nargin < 7
    e = 0.01;
    if nargin < 6
        prs = 8;
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
nv = size(Y,2);
[m,p] = size(X);
if (length(d)~=m) || (length(t)~=m) || (size(Y,1)~=m)
    error(['The design matrix X, censorship status vector d, time'...
        ' vector t and data matrix Y must all have the same number of rows.']);
end;
%Sort the data by time. If there is a tie between a failure time and a
%censored time then the failure time goes first.
st_ix = find(d==1);
t1 = t(st_ix);
[t1,t1_ix] = sort(t1);
X1 = X(st_ix(t1_ix),:);
Y1 = Y(st_ix(t1_ix),:);
cs_ix = find(d==0);
if ~isempty(cs_ix)
    t2 = t(cs_ix);
    [t2,t2_ix] = sort(t2);
    X2 = X(cs_ix(t2_ix),:);
    Y2 = Y(cs_ix(t2_ix),:);
    count1 = 1; count2 = 1; i = 0;
    while (count1 <= length(t1)) && (count2 <= length(t2))
        i = i + 1;
        if t1(count1) <= t2(count2)
            X(i,:) = X1(count1,:);
            Y(i,:) = Y1(count1,:);
            d(i) = 1;
            t(i) = t1(count1);
            count1 = count1 + 1;
        else 
            X(i,:) = X2(count2,:);
            Y(i,:) = Y2(count2,:);
            d(i) = 0;
            t(i) = t2(count2);
            count2 = count2 + 1;
        end;
    end;
    if (count1 > length(t1))
        X(i+1:end,:) = X2(count2:end,:);
        Y(i+1:end,:) = Y2(count2:end,:);
        d(i+1:end) = 0;
        t(i+1:end) = t2(count2:end);
    else
        X(i+1:end,:) = X1(count1:end,:);
        Y(i+1:end,:) = Y1(count1:end,:);
        d(i+1:end) = 1;
        t(i+1:end) = t1(count1:end);
    end;
else
    X = X1;
    Y = Y1;
    t = t1;
end;
%indices of unique failure times in ft_ix (last index when ties happen)
st_ix = find(d==1);
[ft,ft_ix] = unique(t(st_ix),'last');
ft_ix = st_ix(ft_ix);
nft = length(ft);
%handling ties in failure times
nties = zeros(1,nft);
for j=1:nft
    i = 1;
    while (ft_ix(j)-i>0) && (ft(j)==t(ft_ix(j)-i))
        i = i + 1;
    end;
    nties(j) = i-1;
end;
%Initialization
p = p+1;
stats1(1:nv) = struct('Bhat',zeros(p,1),'zscore',zeros(p,1),'pval',zeros(p,1),...
    'CovBhat',zeros(p,p),'llh',0);
st1 = false(nv,1);
if (prs==1) || (matlabpool('size') ~= prs)
    if (matlabpool('size') > 0)
        matlabpool close;
    end;
    if (prs>1)
        matlabpool(prs);
    end;
end;
display(' ');
display('Starting model fitting at each location ...');
fname = parfor_progress('init',nv);
parfor j=1:nv
    XD =  [X Y(:,j)];
    Bhat = zeros(p,1);
    % Iterations
    nit = 20;
    gnorm = e+1;
    it = 1;
    warning off all
    while (gnorm>e) && (it<=nit)
        gr = SStat_Gradient(XD,Bhat,ft_ix,nties);
        He = SStat_Hessian(XD,Bhat,ft_ix,nties);
        if (cond(He) < 1e+10)
            invHe = He\eye(p);
        else
            [Vtemp,Dtemp] = eig(He);
            invHe = Vtemp*diag(1./max(diag(Dtemp),1e-5))*Vtemp';
        end
        Bhat = Bhat - invHe*gr;
        gnorm = norm(gr);
        it = it+1;
    end;
    llh = SStat_Likelihood(XD,Bhat,ft_ix,nties);
    z_sc = Bhat./sqrt(diag(-invHe));
    pval = 2*(1-normcdf(abs(z_sc),0,1));
    stats1(j) = struct('Bhat',Bhat,'zscore',z_sc,'pval',pval,'CovBhat',-invHe,'llh',llh);
    if (gnorm<=e)
        st1(j) = true;
    end;
    if mod(j,1000) == 0
        parfor_progress(fname,1000);
    end;
end;
parfor_progress(fname,0);
if (matlabpool('size') > 0)
    matlabpool close;
end;
stats(1:nv0) = struct('Bhat',[],'zscore',[],'pval',[],'CovBhat',[],'llh',[]);
st = false(nv0,1);
for i=1:nv
    stats(maskvtx(i)) = stats1(i);
    st(maskvtx(i)) = st1(i);
end;
%Summary
display('  ');
display('Summary:');
display(['Algorithm did not converge at ' num2str(100*sum(st1==false)/nv)...
    ' percent of the total number of locations.']);
et = toc;
display(['Total elapsed time is ' num2str(et) ' seconds']);
end









%% Likelihood, Gradient and Hessian

function llk = SStat_Likelihood(X,Bhat,ft_ix,nties)
% 
% Log-likelihood value.
%
% Input
% X: Ordered design Matrix (according to survival time).
% Bhat: Estimated vector of the population regression parameters.
% ft_ix: Failure time indices in X (last index if any tie).
% nties: Number of ties for each survival time.
%
% Output
% llk: Log-likelihood value.
%
llk = 0;
nft = length(ft_ix);
for j=1:nft
    term = 0;
    lpr = X(ft_ix(j)-nties(j):ft_ix(j),:)*Bhat;
    aux1 = sum(exp(lpr));
    aux2 = sum(exp(X(ft_ix(j)-nties(j):end,:)*Bhat));
    for l=0:nties(j)
        term = term + log(aux2-l*aux1/(nties(j)+1));
    end;
    llk = llk + sum(lpr) - term;
end;
end


function gr = SStat_Gradient(X,Bhat,ft_ix,nties)
% 
% Gradient vector for the log-likelihood.
%
% Input
% X: Ordered design Matrix (according to survival time).
% Bhat: Estimated vector of the population regression parameters.
% ft_ix: Failure time indices in X (last index if any tie).
% nties: Number of ties for each survival time.
%
% Output
% gr: Gradient vector.
%
p = size(X,2);
gr = zeros(p,1);
nft = length(ft_ix);
for j=1:nft
    term1 = sum(X(ft_ix(j)-nties(j):ft_ix(j),:),1);
    term2 = exp(X(ft_ix(j)-nties(j):end,:)*Bhat);
    term3 = exp(X(ft_ix(j)-nties(j):ft_ix(j),:)*Bhat);
    term4 = term2'*X(ft_ix(j)-nties(j):end,:);
    term5 = term3'*X(ft_ix(j)-nties(j):ft_ix(j),:);
    term = 0;  
    for l=0:nties(j)
        term = term + (term4-l*term5/(nties(j)+1))/(sum(term2)-l*sum(term3)/(nties(j)+1));
    end;
    gr = gr + (term1 - term)';
end;
end


function He = SStat_Hessian(X,Bhat,ft_ix,nties)
% 
% Hessian matrix for the log-likelihood.
%
% Input
% X: Ordered design Matrix (according to survival time).
% Bhat: Estimated vector of the population regression parameters.
% ft_ix: Failure time indices in X (last index if any tie).
% nties: Number of ties for each survival time.
%
% Output
% He: Hessian matrix.
%
[m,p] = size(X);
He = zeros(p,p);
nft = length(ft_ix);
for j=1:nft
    term1 = 0; 
    for i=ft_ix(j)-nties(j):m
        term1 = term1 + exp(X(i,:)*Bhat)*X(i,:)'*X(i,:);
    end;
    term2 = 0; 
    for i=ft_ix(j)-nties(j):ft_ix(j)
        term2 = term2 + exp(X(i,:)*Bhat)*X(i,:)'*X(i,:);
    end;
    term3 = sum(exp(X(ft_ix(j)-nties(j):end,:)*Bhat));
    term4 = sum(exp(X(ft_ix(j)-nties(j):ft_ix(j),:)*Bhat));
    term5 = exp(X(ft_ix(j)-nties(j):m,:)*Bhat)'*X(ft_ix(j)-nties(j):m,:);
    term6 = exp(X(ft_ix(j)-nties(j):ft_ix(j),:)*Bhat)'*X(ft_ix(j)-nties(j):ft_ix(j),:);
    term = 0;  
    for l=0:nties(j)
        Z = term5 - l*term6/(nties(j)+1);
        phi = term3 - l*term4/(nties(j)+1);
        term = term + (term1-l*term2/(nties(j)+1))/phi - Z'*Z/(phi*phi);                   
    end;
    He = He - term;
end;
end


