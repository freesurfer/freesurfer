function [stats,st] = lme_mass_fit_Rgw(X,Zcols,Y,ni,Th0,Rgs,Surf,fname,...
                                                          Dtype,sptm,prs,e)                                                              
% [stats,st] = lme_mass_fit_Rgw(X,Zcols,Y,ni,Th0,Rgs,Surf,fname,Dtype,sptm,prs,e)
%
% Region-wise linear mixed-effects estimation. This function estimates the
% parameters at each location by pooling covariance components over all 
% locations inside homogeneous regions of the errors. To achive that goal 
% the spatial covariance inside each region is modeled using Exponential or 
% Gaussian spatial correlation models (the exponential model seems to work
% best in neuroimaging).
%
% Input
% X: Ordered design matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% Y: Ordered data matrix (n x nv, n=#total number of scans and nv=# of 
% vertices/voxels).
% ni: Vector whose entries are the number of repeated measures for each
% subject (ordered according to X).
% Th0: Matrix whose colums are initial estimators of the covariance 
% components at each location.
% Rgs: 1 x nv segmentation vector containing a region number assigned to
% each location.
% Surf: Surface which represent a coordinate system where the analysis is 
% to be done. It is a structure with Surf.tri = t x 3 matrix of triangle 
% indices, 1-based, t=#triangles and Surf.coord = 3 x nv matrix of 
% coordinates, nv=#vertices.
% fname: File name to save outputs. Default [] (no output is saved to any
% file).
% Dtype: Type of distances to be computed among the surface nodes. It is
% 'euc' for Euclidean or 'surf' for geodesic distances along the surface. 
% Default 'euc'.
% sptm: Spatial model for the covariance matrix inside the regions. It can 
% be 'exp' or 'gauss'. Default 'exp'.
% prs: Number of workers for parallel computing. Default 8;
% e: Convergence epsilon (gradient's norm). Default 10^-1;
%
% Output
% stats: Structure array containing statistiscs for every voxel/vertex 
% (see lme_FSfit for more details on these statistics).
% st: Array containing the termination state for each location.
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
tic;
%Validation of inputs
if nargin < 7
    error('Too few inputs');
elseif nargin < 12
        e = 10^-1;
        if nargin < 11
            prs = 8;
            if nargin < 10
                sptm = 'exp';
                if nargin < 9
                    Dtype = 'euc';
                    if nargin < 8
                        fname = [];
                    end;
                end;
            end;
        end;
end;
if (size(Y,1) ~= size(X,1)) 
    error('X and Y must have the same number of rows.');
end;
n = sum(ni);
if (n ~= size(X,1)) 
    error(['The total number of measurements, indicated by sum(ni), must'...
       'be the same as the number of rows of the design matrix X']);
end;
nv0 = size(Y,2);
if (nv0 ~= size(Th0,2)) 
    error('Y and Theta0 must have the same number of colums.');
end;
%to make slightly different two vertices with same coord in FreeSurfer
%meshes
[~,loc] = unique(Surf.coord','rows');
pos = find(~ismember(1:nv0,loc));
if ~isempty(pos)
    Surf.coord(3,pos) = Surf.coord(3,pos) + 0.001;
end;
maskvtx = find(Rgs ~= 0);
nv = length(maskvtx);
Rgs = Rgs(maskvtx);
Y = Y(:,maskvtx);
Th0 = Th0(:,maskvtx);
Rgnums = unique(Rgs);
nRg = length(Rgnums);
Rglth = zeros(1,nRg);
for i=1:nRg
    Rglth(i) = sum(Rgs == Rgnums(i));
end;
%Sort regions by length in descend order
[~,ix] =  sort(Rglth,'descend');
aux = Rgs;
for i=1:nRg
    Rgs(aux == Rgnums(ix(i))) = i;
end;
%Balance region size across workers
[bRgind,prsnRg] = balanceRgind(prs,nRg);
%Initialization
statstruct = struct('Bhat',[],'CovBhat',[],'phisqhat',[],'Dhat',[],...
    'Zcols',[],'invEI',[],'Pth',[],'Qthth',[],'lreml',-10^10);
if (prs==1) || (matlabpool('size') ~= prs)
    if (matlabpool('size') > 0)
        matlabpool close;
    end;
    if (prs>1)
        matlabpool(prs);
    end;
end;
prsnumloc = zeros(prs,1);
parfor j=1:prs
    for i=1:prsnRg(j)
        nRgvtx = sum(Rgs == bRgind(j,i));
        prsnumloc(j) = prsnumloc(j) + nRgvtx;
    end;
end;
nrf = size(Th0,1);
for j=1:prs
    Rgstats{j}(1:prsnumloc(j)) = statstruct;
    Rgsts{j} = false(prsnumloc(j),1);
    Yj{j} = zeros(n,prsnumloc(j));
    Th0j{j} = zeros(nrf,prsnumloc(j));
end;
for j=1:prs
    posi = 0;
    for i=1:prsnRg(j)
        RgVtxInd = find(Rgs == bRgind(j,i));
        nRgvtx = length(RgVtxInd);
        posf = posi+nRgvtx;
        Yj{j}(:,posi+1:posf) = Y(:,RgVtxInd); 
        Th0j{j}(:,posi+1:posf) = Th0(:,RgVtxInd); 
        posi = posf;
    end;
end;
clear Y Th0 
%Estimation
display(' ');
display('Starting model fitting at each region ...');
display(' ');
sRgst = 0;
fn = parfor_progress('init',nv);
parfor j=1:prs
    warning off all
    progress = 0;
    posi = 0;
    for i=1:prsnRg(j)
        RgVtxInd = find(Rgs == bRgind(j,i));
        nRgvtx = length(RgVtxInd);
        posf = posi+nRgvtx;
        try
            if nRgvtx > 1
                Dist = lme_mass_RgDist(Surf,maskvtx,RgVtxInd,Dtype);
                [Rgstats{j}(posi+1:posf),Rgsts{j}(posi+1:posf)] = lme_RgFSfit(X,...
                    Zcols,Yj{j}(:,posi+1:posf),ni,Th0j{j}(:,posi+1:posf),Dist,sptm);
            else
                [Rgstats{j}(posf),Rgsts{j}(posf)] = lme_FSfit(X,...
                    Zcols,Yj{j}(:,posf),ni,e);
            end
        catch em
            display(['Region ' num2str(bRgind(j,i)) ': ' em.message]);
        end;
        if convergence(Rgsts{j}(posf),Rgstats{j}(posf).lreml,bRgind(j,i));
            sRgst = sRgst + 1;
        end;
        progress = progress + nRgvtx;
        if mod(bRgind(j,i),10) == 0
            parfor_progress(fn,progress);
            progress = 0;
        end;
        posi = posf;
    end;
end;
parfor_progress(fn,0);
if (matlabpool('size') > 0)
    matlabpool close;
end;
stats(1:nv0) = statstruct;
st = false(nv0,1);
for j=1:prs
    posi = 0;
    for i=1:prsnRg(j)
        RgVtxInd = find(Rgs == bRgind(j,i));
        nRgvtx = length(RgVtxInd);
        posf = posi+nRgvtx;
        Rgstatsji = Rgstats{j}(posi+1:posf);
        Rgstsji = Rgsts{j}(posi+1:posf);
        for k=1:nRgvtx           
            stats(maskvtx(RgVtxInd(k))) = Rgstatsji(k);
            stats(maskvtx(RgVtxInd(k))).lreml = stats(maskvtx(RgVtxInd(k))).lreml(end);
            st(maskvtx(RgVtxInd(k))) = Rgstsji(k);
        end;
       posi = posf; 
    end;
end;
%Summary
display('  ');
display('Summary:');
display(['Algorithm did not converge at ' num2str((nRg-sRgst)*100/nRg)...
    ' percent of the total number of regions.']);
et = toc;
display(['Elapsed time is ' num2str(et/60) ' minutes.']);
if ~isempty(fname)
    display('  ');
    display('Saving output ...');
    save(fname, 'stats','st','-v7.3');
end;
end






%% AUXILIAR FUNCTIONS

function tf = convergence(Rgst,lreml,i)
tf = false;
if Rgst  && (lreml(end) >= lreml(1))
    display(['Region ' num2str(i) ': Convergence at iteration ' ...
        num2str(length(lreml)) '. Initial and final likelihoods'...
        ': ' num2str(lreml(1)) ', ' num2str(lreml(end)) '.']);
    tf = true;
elseif Rgst && (lreml(end) < lreml(1))
    display(['Region ' num2str(i) ': Convergence to a saddle point '...
        'at iteration ' num2str(length(lreml)) '. Initial and '...
        'final likelihoods: ' num2str(lreml(1)) ', ' ...
        num2str(lreml(end)) '.']);
else
    display(['Region ' num2str(i) ': Algorithm did not converge.'...
        ' Initial and final likelihoods: ' num2str(lreml(1))...
        ', ' num2str(lreml(end)) '.']);
end;
end




function [bRgind,prsnRg] = balanceRgind(prs,nRg)
bRgind = zeros(prs,ceil(nRg/prs));
prsnRg = zeros(prs,1);
reverse = false;
i = 1;
while i <= nRg   
    if reverse
        j = prs;
        while (i <= nRg) && (j >=1);
            prsnRg(j) = prsnRg(j) + 1;
            bRgind(j,prsnRg(j)) = i;
            j = j - 1;
            i = i + 1;
        end;
        reverse = false;
    else
        j = 1;
        while (i <= nRg) && (j <= prs)
            prsnRg(j) = prsnRg(j) + 1;
            bRgind(j,prsnRg(j)) = i;
            j = j + 1;
            i = i + 1;
        end;
        reverse = true;
    end;
end;
end

