function [sz1,sz2] = lme_plannedSampleSize(Zi,ZiCol,Dhat,phisqhat,effsz,...
                                                         dr,pw,alpha,gr_pr)
% [sz1,sz2] = lme_plannedSampleSize(Zi,ZiCol,Dhat,phisqhat,effsz,dr,pw,alpha,gr_pr)                                                                                                                
%
% Sample size calculation for a planned balanced design of a linear
% mixed-effects model for a prospective comparison of one random effect 
% between two groups (eg.intercept or regression slope). Depends on the
% Statistics toolbox.
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
% pw: Power of the test. Default 0.8.
% alpha: Significance level of the test. Default 0.05.
% gr_pr: This parameter determines the group proportions in the final sample
% size. It can be used to requiere a different sample size for each group: 
% total_sz*gr_pr and total_sz*(1-gr_pr). Default 0.5 (equal size).
% 
%
% Output
% sz1: Requiered sample size for the first group.
% sz2: Requiered sample size for the second group.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:11 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:11 $
%    $Revision: 1.1.2.2 $
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
if nargin < 6
    error('Too few inputs');
elseif nargin < 9
    gr_pr = 0.5;
    if nargin < 8
        alpha = 0.05;
        if nargin < 7
            pw = 0.8;
        end;
    end
end;
if (dr <= 0) || (dr >= 1)
    error('Input "dr" must be greater than zero and lesser than one');
end;
if (gr_pr <= 0) || (gr_pr >= 1)
    error('Input "p" must be greater than zero and lesser than one');
end;
if (pw <= 0) || (pw > 1)
    error('Input "pw" must be greater than zero and lesser than or equal one');
else
    CovBihat = phisqhat*(Zi'*Zi)^-1+Dhat;
    phisq = CovBihat(ZiCol,ZiCol);
    n = (norminv(1-alpha/2,0,1)+norminv(pw,0,1))^2*phisq/...
        (2*gr_pr*(1-gr_pr)*effsz^2);
    sz1 = 2*n*gr_pr/(1-dr);
    sz2 = 2*n*(1-gr_pr)/(1-dr);
end;

