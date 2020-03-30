function pw = lme_plannedPower(Zi,ZiCol,Dhat,phisqhat,effsz,dr,sz1,sz2,alpha)                                                                   
% pw = lme_plannedPower(Zi,ZiCol,Dhat,phisqhat,effsz,dr,sz1,sz2,alpha)                                                                  
%
% Power calculation for a planned balanced design of a linear mixed-effects
% model for a prospective comparison of one random effect between two groups
% (eg.intercept or regression slope). Depends on the Statistics toolbox.
%
% Input
% Zi: Common random effects design matrix for the subjects (a balanced
% design is expected in advance).
% ZiCol: Colum of Zi corresponding to the effect to be compared between the 
% two groups.
% Dhat: Estimated random effects covariance matrix (eg. from previous studies).
% phisqhat: Estimated intra-subject variability (eg. from previous studies).
% effsz: Effect size (can be the size of a detected effect in a previous
% study).
% dr: Expected drop out rate (proportion of subjects who drop out before the
% completion of the study). Must be greater than zero and lesser than one.
% sz1: Sample size of the first group.
% sz2: Sample size of the second group.
% alpha: Significance level of the test. Default 0.05.
%
% Output
% pw: Power of the test.
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
if nargin < 8
    error('Too few inputs');
elseif nargin < 9
    alpha = 0.05;
end;
if (dr <= 0) || (dr >= 1)
    error('Input "dr" must be greater than zero and lesser than one');
end;
CovBihat = phisqhat*(Zi'*Zi)^-1+Dhat;
phisq = CovBihat(ZiCol,ZiCol);
n = 0.5*(sz1+sz2)*(1-dr);
zpw = sqrt(2*sz1*sz2*n*effsz^2/phisq)/(sz1+sz2)-norminv(1-alpha/2,0,1);
pw = normcdf(zpw,0,1);


