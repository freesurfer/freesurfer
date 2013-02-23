function pval = lme_mass_LR(statsfull,statsred,q)
% pval = lme_mass_LR(statsfull,statsred,q) 
%
% Likelihood ratio test for the random effects. It can be used to test if a
% model with q+1 random effects is significantly better than a model with q
% random effects.
%
% Input
% statsfull: Structure array containing statistiscs for every voxel/vertex
% (see lme_FSfit for more details on these statistics) for the full 
% model (the one with q+1 random effects).
% statsred: Structure array containing statistiscs for every voxel/vertex
% for the reduced model (the one with q random effects).
% q: Number of random effects in the reduced model.
%
% Output
% pval: P-value of the test at each voxel/vertex (based on a 50:50 mixture  
% of chi-squared distributions with q and q+1 degrees of freedom). 
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

if nargin < 3 
    error('Too few inputs');   
end;

nv = length(statsfull);
pval = zeros(1,nv);
for i=1:nv
    lrstats = lme_LR(statsfull(i).lreml,statsred(i).lreml,q);
    pval(i) = lrstats.pval;
end;
