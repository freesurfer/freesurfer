function [stats,st] = lme_mass_fit(X,Xcols,Xrows,Zcols,Y,ni,prs,e)
% [stats,st] = lme_mass_fit(X,Xcols,Xrows,Zcols,Y,ni,prs,e)
%
% Location-wise linear mixed-effects estimation. Allows to have different 
% models across locations (this is useful when there are either missing  
% values at some locations or the number of fixed effects and/or the number
% of random effects varies across locations). The model at each location is
% always a subset of the "maximal model" given by X.
%
% Input
% X: Ordered (according to time for each subject) design matrix (nmxp, nm 
% total # of maps, p # of fixed effects) for a "maximal linear model". It 
% is the more complex model under consideration across locations.
% Xcols: Matrix (nvxp, nv #vertices) whos rows contain ones or zeros 
% indicating which colums of X are/are not going to be considered at each
% voxel/vertex respectively. Set Xcols to the empty matrix ([]) to consider
% all colums of X at each voxel/vertex.
% Xrows: Matrix (nmxnv) whos colums contain ones or zeros indicating which 
% rows of X are/are not going to be considered at each voxel/vertex
% respectively. Set Xrows to the empty matrix ([]) to consider all rows 
% of X at each voxel/vertex.
% Zcols: Matrix (nvxp) whos rows contain ones or zeros indicating which 
% colums of X are/are not going to be considered as random effects at each
% voxel/vertex. Each row of Zcols can only contain ones at those entries 
% for which the corresponding entries in Xcols is also one (the random 
% effects' design matrix must be a subset of X). If Zcols is a vector then
% the same random effects are considered across locations. 
% Y: Data matrix (nmxnv) whos colums represent the ordered data vector
% (according to X) for each voxel/vertex.
% ni: Vector whose entries are the number of repeated measures for each
% subject (ordered according to X) .
% e: Convergence epsilon (gradient's norm). Default 10^-1;
% prs: Number of workers for parallel computing. Default 8;
%
% Output
% stats: Structure array containing statistiscs for every voxel/vertex (see
% lme_massfit_FS for more details on these statistics).
% st: Array containing the termination state for voxel/vertex
% (1 for convergence and 0 otherwise).
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

% tic;
%Validation of inputs
if nargin < 6
    error('Too few inputs');
elseif nargin < 8
    e = 10^-1;   
    if nargin < 7
       prs = 8; 
    end;
end;
[nm,nv] = size(Y);
if (nm ~= size(X,1)) 
    error(['The design matrix X and the data matrix Y must have the same'...
        ' number of rows.']);
end;
p = size(X,2);
if isempty(Xcols)
    Xcols = ones(nv,p,'int8');
else
    if (p ~= size(Xcols,2))
        error(['The second input matrix Xcols must have the same # of'...
            ' columns as the design matrix X']);
    end;
    if (nv ~= size(Xcols,1))
        error(['The second input matrix Xcols must have the same # of '...
            'rows as the # of columns of the data matrix Y']);
    end;
    Xcols = int8(Xcols);
end;
if (p ~= size(Zcols,2))
    error(['The fourth input matrix Zcols must have the same # of colums'...
        ' as the design matrix X']);
end;
if size(Zcols,1) == 1
    rfcols = find(Zcols == 1);
    Zcols = zeros(nv,p,'int8');
    Zcols(:,rfcols) = ones(nv,length(rfcols),'int8');
else
    if (nv ~= size(Zcols,1))
        error(['The fourth input matrix Zcols must have the same # of '...
            'rows as the # of columns of the data matrix Y']);
    end;
    Zcols = int8(Zcols);    
end;
if ~isempty(Xrows)
    if (size(Xrows,1) ~= nm)
        error(['The number of rows of the third input matrix Xrows must'...
            ' be the same as the number of rows of the design matrix X.']);
    end;
    if (size(Xrows,2) ~= nv)
        error(['The number of columns of the third input matrix Xrows must'...
            ' be the same as the number of colums of the data matrix Y.']);
    end;
    Xrows = int8(Xrows);
end;
if (sum(ni) ~= nm) 
    error(['The total number of measurements, indicated by sum(ni), must'...
       'be the same as the number of rows of the design matrix X']);
end;
for i=1:nv
    cols = Xcols(i,Zcols(i,:) == 1);
    if cols ~= ones(1,length(cols))
        error(['Inconsistence at location ' num2str(i) ' Each row of '... 
            'Zcols must has ones for a subset of the one-elements in '...
            'the corresponding row of Xcols.']);
    end;
end;
%Initialization
stats(1:nv) = struct('Bhat',[],'CovBhat',[],'phisqhat',[],'Dhat',[],...
             'Zcols',[],'invEI',[],'Pth',[],'Qthth',[],'lreml',-10^10);
st = false(nv,1);
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
display(' ');
fn = parfor_progress('init',nv);
if ~isempty(Xrows)
    sID = zeros(nm,1);
    m = length(ni);
    count = 0;
    for i=1:m
        sID(count+1:count+ni(i)) = i;
        count = count+ni(i);
    end;
    parfor i=1:nv
        %Build the model for this specific location
        rows = find(Xrows(:,i) == 1);
        XM = X(rows,Xcols(i,:) == 1);
        sIDM = sID(rows);
        mM = length(unique(sIDM));
        niM = ones(mM,1,'uint8');
        count = 1;
        for j=2:length(rows)
            if sIDM(j-1)==sIDM(j)
                niM(count) = niM(count) + 1;
            else
                count = count + 1;
            end;
        end;
        rfcols = find(Zcols(i,:) == 1);
        warning off all
        %Estimation
        try
            [stats(i), st(i)] = lme_FSfit(XM,rfcols,Y(rows,i),niM,e);
        catch em
            display(['Location ' num2str(i) ': ' em.message]);
        end;
        displayConvergence(st(i),stats(i).lreml,i);
        stats(i).lreml = stats(i).lreml(end);
        if mod(i,20) == 0
            parfor_progress(fn,20);
        end;
    end;
else
    parfor i=1:nv
        %Build the model for this specific location
        XM = X(:,Xcols(i,:) == 1);
        rfcols = find(Zcols(i,:) == 1);
        warning off all
        %Estimation
        try
            [stats(i), st(i)] = lme_FSfit(XM,rfcols,Y(:,i),ni,e);
        catch em
            display(['Location ' num2str(i) ': ' em.message]);
        end;
        displayConvergence(st(i),stats(i).lreml,i);
        stats(i).lreml = stats(i).lreml(end);
        if mod(i,20) == 0
            parfor_progress(fn,20);
        end;
    end;
end;
parfor_progress(fn,0);
if (matlabpool('size') > 0)
    matlabpool close;
end;
%et = toc;
%display(['Elapsed time is ' num2str(et/60) ' minutes.']);
end





function  displayConvergence(st,lreml,i)
    if st && (lreml(end) >= lreml(1))
        display(['Location ' num2str(i) ': Convergence at iteration ' ...
        num2str(length(lreml)) '. Initial and final likelihoods: '...
         num2str(lreml(1)) ', ' num2str(lreml(end)) '.']);
    elseif st && (lreml(end) < lreml(1))
        display(['Location ' num2str(i) ': Convergence to a saddle point '...
        'at iteration ' num2str(length(lreml)) '. Initial and final'... 
        ' likelihoods: ' num2str(lreml(1)) ', ' num2str(lreml(end)) '.']);      
    else
        display(['Location ' num2str(i) ': Algorithm did not converge.'...
           ' Initial and final likelihoods: ' num2str(lreml(1)) ', ' ...
             num2str(lreml(end)) '.']);
    end;
end

