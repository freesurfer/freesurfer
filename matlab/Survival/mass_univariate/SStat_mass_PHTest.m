function [pval,rho] = SStat_mass_PHTest(X,d,t,Y,stats_PH,maskvtx)
% [pval,rho] = SStat_mass_PHTest(X,d,t,Y,stats_PH,maskvtx)
%
% Vertex/voxel-wise Schoenfeld residuals test for the proportional hazards 
% assumption.
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
% stats_PH: Structure array containing statistics obtained from SStat_mass_CoxPH.
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
%
% Output
% pval: Matrix of p-values (pxnv) indicating the probability of the PH  
% assumption for each covariate and at each location.
% rho: Matrix of correlations (pxnv) between Schoenfeld residuals and ranked   
% failure times at each location.
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Kleinbaum, D.G., Klein, M., 2005. Survival analysis. A self-
% learning approach, second edition. New York: Springer..
%
%Validation of inputs
if nargin < 5
    error('Too few inputs');
elseif nargin < 6
    maskvtx = [];
end;
nv0 = size(Y,2);
if (length(stats_PH)~=nv0)
    error(['The structure array stats_PH must have the same length as the'...
        ' number of colums of Y.']);
end;
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
%handling ties in failure times by substracting a very small random number
rand('state',sum(100*clock));
for j=1:nft
    i = 1;
    while (ft_ix(j)-i>0) && (ft(j)==t(ft_ix(j)-i))
        i = i + 1;
    end;
    nties = i-1;
    tt = t(ft_ix(j)-nties:ft_ix(j)) - 1e-5*(1 + rand(nties+1,1));
    tX = X(ft_ix(j)-nties:ft_ix(j),:);
    tY = Y(ft_ix(j)-nties:ft_ix(j),:);
    [stt,stt_ix] = sort(tt);
    t(ft_ix(j)-nties:ft_ix(j)) = stt;
    X(ft_ix(j)-nties:ft_ix(j),:) = tX(stt_ix,:);
    Y(ft_ix(j)-nties:ft_ix(j),:) = tY(stt_ix,:);
end;
nft = length(st_ix);
p = p+1;
pval = ones(p,nv);
rho = zeros(p,nv);
r = zeros(nft,p);
display('Checking the proportional hazards assumption at each location ...');
for j=1:nv
    if ~isempty(stats_PH(j).Bhat)
        if length(stats_PH(j).Bhat)==p
            XD =  [X Y(:,j)];
            for i=1:nft
                lprv = exp(XD(st_ix(i):end,:)*stats_PH(j).Bhat);
                r(i,:) = XD(st_ix(i),:) - (lprv'*XD(st_ix(i):end,:))./sum(lprv);
            end;
            [rho(:,j),pval(:,j)] = corr(r,[1:nft]');
        else
            error([' The number of elements of stats(' num2str(j) ').Bhat'...
                ' is different from ' num2str(p)]);
        end;       
    end;
end

