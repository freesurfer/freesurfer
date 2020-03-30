function [stats,st] = lme_mass_fit_vw(X,Zcols,Y,ni,maskvtx,fname,prs,e,Xrows)
% [stats,st] = lme_mass_fit_vw(X,Zcols,Y,ni,maskvtx,fname,prs,e,Xrows)
%
% Vertex/voxel-wise linear mixed-effects estimation. 
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
% fname: File name to save outputs. Default [] (no output is saved to any
% file).
% prs: Number of workers for parallel computing. Default numcores;
% e: Convergence epsilon (gradient's norm). Default 10^-1;
% Xrows: Optional matrix (nmxnv) whos colums contain ones or zeros indicating 
% which rows of X are/are not going to be considered at each voxel/vertex
% respectively (this is useful when there are missing data values at some 
% locations. Default, all rows of X are considered at each voxel/vertex.
%
% Output
% stats: Structure array containing statistiscs for every voxel/vertex (see
% lme_FSfit for more details on these statistics).
% st: Array containing the termination state for each voxel/vertex
% (1 for convergence and 0 otherwise).
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
tic;
%Validation of inputs
if nargin < 4
    error('Too few inputs');
elseif nargin < 9
    Xrows = [];
    if nargin < 8
        e = 10^-1;
        if nargin < 7
            prs = feature('numcores');
            if nargin < 6
                fname = [];
                if nargin < 5
                    maskvtx = [];
                end;
            end;
        end;
    end;
end;
nv0 = size(Y,2);
if isempty(maskvtx)
    maskvtx = 1:nv0;
end;
Y = Y(:,maskvtx);
p = size(X,2);
nv = size(Y,2);
rfcols = Zcols;
Zcols = zeros(1,p,'int8');
Zcols(:,rfcols) = ones(1,length(rfcols),'int8');
if ~isempty(Xrows)
    Xrows = int8(Xrows(:,maskvtx));
end;
[stats1,st1] = lme_mass_fit(X,[],Xrows,Zcols,Y,ni,prs,e);
stats(1:nv0) = struct('Bhat',[],'CovBhat',[],'phisqhat',[],'Dhat',[],...
    'Zcols',[],'invEI',[],'Pth',[],'Qthth',[],'lreml',-10^10);
st = false(nv0,1);
for i=1:nv
    stats(maskvtx(i)) = stats1(i);
    st(maskvtx(i)) = st1(i);
end;

%Summary
display('  ');
display('Summary:');
display(['Algorithm did not converge at ' num2str(100*(1-sum(st==1)/nv))...
    ' percent of the total number of locations.']);
et = toc;
display(['Total elapsed time is ' num2str(et/60) ' minutes.']);
if ~isempty(fname)
    display('  ');
    display('Saving output ...');
    save(fname, 'stats','st','-v7.3');
end;

