function [Rfx,nrfx,Bhat] = lme_mass_rfx(stats,X,Zcols,Y,ni,maskvtx)
% [Rfx,nrfx,Bhat] = lme_mass_rfx(stats,X,Zcols,Y,ni,maskvtx)
%
% Returns the subject-specific random effects estimates at each vertex.
%
% Input
% stats: Structure array containing statistics for every voxel/vertex 
% (generated with either lme_mass_fit_Rgw or lme_mass_fit_vw).
% X: Ordered design matrix (according to time for each subject).
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
% Zcols: Vector with the indices of the colums of X that are considered as
% random effects.
% Y: Ordered data matrix (n x nv, n=total number of scans and nv=number of 
% vertices/voxels).
% ni: Vector whose m entries are the number of repeated measures for each
% subject (ordered according to X).
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
%
% Output
% Rfx: Estimated subject-especific random effects matrix (m x nrfx*nv). The columns 
% of this matrix are grouped by vertex. For example if there are two random effects
% in the model then the first two columns contain the subject-specific random effect 
% coefficients for the first vertex, then the next two columns contain the 
% subject-specific random effect coefficients for the second vertex and so on...
% nrfx: Number of random effects (length(Zcols)).
% Bhat: Population-level regression coefficients in stats stacked in one matrix.
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
nv0 = size(Y,2);
if isempty(maskvtx)
    maskvtx = 1:nv0;
end;
m = length(ni);
nrfx = length(Zcols);
Rfx = zeros(m,nv0*nrfx);
Bhat = zeros(nrfx,nv0);
for j=1:length(maskvtx)
      %Computation of W = SIGMA^-1 
    D = stats(maskvtx(j)).Dhat;
    phisq = stats(maskvtx(j)).phisqhat;
    r = Y(:,maskvtx(j))-X*stats(maskvtx(j)).Bhat;
    posi = 1; 
    scInvD = D\eye(nrfx)*phisq;
    ijcol = (maskvtx(j)-1)*nrfx + 1; fjcol = ijcol + nrfx - 1;
    for i=1:m
        posf = posi+ni(i)-1;
        Zi = X(posi:posf,Zcols);
        Wi = (eye(ni(i))-Zi/(Zi'*Zi+scInvD)*Zi')/phisq;
        ri = r(posi:posf);
        Rfx(i,ijcol:fjcol) = (D*Zi'*Wi*ri)';
        posi = posf+1;      
    end;  
    Bhat(:,maskvtx(j)) = stats(maskvtx(j)).Bhat(Zcols);
end;

