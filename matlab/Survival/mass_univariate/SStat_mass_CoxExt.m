function [stats,st] = SStat_mass_CoxExt(sID_ext,X_ext,d_ext,t_ext,Y,maskvtx,prs,e)
% [stats,st] = SStat_mass_CoxExt(sID_ext,X_ext,d_ext,t_ext,Y,maskvtx,prs,e)
%
% Vertex/voxel-wise parameter estimation for the extended Cox model. This 
% function uses as input the output of SStat_X_ext but with X_ext having 
% new time-dependent columns added by the user (eg. the product of a column
% of X_ext with t_ext. Survival time ties are handled using Efron's method. 
%
% Input
% sID_ext: Extended subjects' IDs (nx1).
% X_ext: Extended design matrix (nxp).
% d_ext: Extended censorship status vector (nx1).
% t_ext: Extended survival time vector (nx1).
% Y: Time-dependent data matrix (nxnv, nv #vertices) whos colums represent 
% the ordered data vector (according to X_ext) for each voxel/vertex.
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
if nargin < 5
    error('Too few inputs');
elseif nargin < 8
    e = 0.01;
    if nargin < 7
        prs = 8;
        if nargin < 6
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
[n,p] = size(X_ext);
if (length(sID_ext)~=n) || (length(d_ext)~=n) || (length(t_ext)~=n) || (size(Y,1)~=n)
        error(['The design matrix X_ext, censorship status vector d_ext, time'...
        ' vector t_ext, subject ID vector sID_ext and data matrix Y must all'...
        ' have the same number of rows.']); 
end;
%indices of unique failure times in ft_ix (last index when ties happen)
st_ix = find(d_ext==1);
[~,ft_ix] = unique(t_ext(st_ix),'last');
ft_ix = st_ix(ft_ix);
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
    XD =  [X_ext Y(:,j)];
    Bhat = zeros(p,1);
    % Iterations
    nit = 20;
    gnorm = e+1;
    it = 1;
    warning off all
    while (gnorm>e) && (it<=nit)
        gr = SStat_Gradient(XD,t_ext,Bhat,ft_ix);
        He = SStat_Hessian(XD,t_ext,Bhat,ft_ix);
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
    llh = SStat_Likelihood(XD,t_ext,Bhat,ft_ix);
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

function llk = SStat_Likelihood(X_ext,t_ext,Bhat,ft_ix)
% 
% Log-likelihood value.
%
% Input
% X_ext: Extended design matrix.
% t_ext: Extended survival time vector.
% Bhat: Estimated vector of the population regression parameters.
% ft_ix: Failure time indices in X_ext (last index if any tie).
%
% Output
% llk: Log-likelihood value.
%
llk = 0;
nft = length(ft_ix);
for j=1:nft
    term = 0;
    ties = ft_ix(t_ext(ft_ix(j))==t_ext(ft_ix));
    nties = length(ties)-1;
    lpr = X_ext(ties,:)*Bhat;
    aux1 = sum(exp(lpr));
    aux2 = sum(exp(X_ext(t_ext(ft_ix(j))==t_ext,:)*Bhat));
    for l=0:nties
        term = term + log(aux2-l*aux1/(nties+1));
    end;
    llk = llk + sum(lpr) - term;
end;
end


function gr = SStat_Gradient(X_ext,t_ext,Bhat,ft_ix)
% 
% Gradient vector for the log-likelihood.
%
% Input
% X_ext: Extended design matrix.
% t_ext: Extended survival time vector.
% Bhat: Estimated vector of the population regression parameters.
% ft_ix: Failure time indices in X_ext (last index if any tie).
%
% Output
% gr: Gradient vector.
%
p = size(X_ext,2);
gr = zeros(p,1);
nft = length(ft_ix);
for j=1:nft
    ties = ft_ix(t_ext(ft_ix(j))==t_ext(ft_ix));
    nties = length(ties)-1;
    riskset = t_ext(ft_ix(j))==t_ext;
    term1 = sum(X_ext(ties,:),1);
    term2 = exp(X_ext(riskset,:)*Bhat);
    term3 = exp(X_ext(ties,:)*Bhat);
    term4 = term2'*X_ext(riskset,:);
    term5 = term3'*X_ext(ties,:);
    term = 0;  
    for l=0:nties
        term = term + (term4-l*term5/(nties+1))/(sum(term2)-l*sum(term3)/(nties+1));
    end;
    gr = gr + (term1 - term)';
end;
end


function He = SStat_Hessian(X_ext,t_ext,Bhat,ft_ix)
% 
% Hessian matrix for the log-likelihood.
%
% Input
% X_ext: Extended design matrix.
% t_ext: Extended survival time vector.
% Bhat: Estimated vector of the population regression parameters.
% ft_ix: Failure time indices in X_ext (last index if any tie).
%
% Output
% He: Hessian matrix.
%
p = size(X_ext,2);
He = zeros(p,p);
nft = length(ft_ix);
for j=1:nft
    ties = ft_ix(t_ext(ft_ix(j))==t_ext(ft_ix));
    nties = length(ties)-1;
    rsk_ix = find(t_ext(ft_ix(j))==t_ext);  
    m = length(rsk_ix);
    term1 = 0;
    for i=1:m
        term1 = term1 + exp(X_ext(rsk_ix(i),:)*Bhat)*X_ext(rsk_ix(i),:)'*X_ext(rsk_ix(i),:);
    end;
    term2 = 0; 
    for i=1:nties
        term2 = term2 + exp(X_ext(ties(i),:)*Bhat)*X_ext(ties(i),:)'*X_ext(ties(i),:);
    end;
    term3 = sum(exp(X_ext(rsk_ix,:)*Bhat));
    term4 = sum(exp(X_ext(ties,:)*Bhat));
    term5 = exp(X_ext(rsk_ix,:)*Bhat)'*X_ext(rsk_ix,:);
    term6 = exp(X_ext(ties,:)*Bhat)'*X_ext(ties,:);
    term = 0;  
    for l=0:nties
        Z = term5 - l*term6/(nties+1);
        phi = term3 - l*term4/(nties+1);
        term = term + (term1-l*term2/(nties+1))/phi - Z'*Z/(phi*phi);                   
    end;
    He = He - term;
end;
end


