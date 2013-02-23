function [D0,phisq0] = lme_fit_init(X,Zcols,y,ni)
% [D0,phisq0] = lme_fit_init(X,Zcols,y,ni) 
% 
% Starting values for linear mixed-effects iterative estimation.
% These starting values are based on the ordinary least squares estimators.
%
% Input
% X: Ordered design Matrix (according to time for each subject).
% Zcols: Vector with the indices of the colums of X that will be considered
% as random effects.
% y: Ordered data vector (according to X).
% ni: Vector whose entries are the number of repeated measures for each
% subject (ordered according to X).
%
% Output
% phisq0: Estimated within-subject variability.
% D0: Estimated random effects covariance matrix.
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
m = length(ni);
p = size(X,2);
q = length(Zcols);
n = sum(ni);
Z = X(:,Zcols);
Bhat = pinv(X)*y;
t1 = 0; t2 = 0; t3 = 0; t4 = 0;
posi = 1;
for i=1:m
    posf = posi+ni(i)-1;
    Zi = Z(posi:posf,:);
    Xi = X(posi:posf,:);
    yi = y(posi:posf);
    ri = yi-Xi*Bhat;
    t = pinv(Zi'*Zi);
    t1 = t1 + t;
    bihat = t*Zi'*ri;
    t2 = t2 + yi'*yi-bihat'*Zi'*ri;
    t3 = t3 + Xi'*yi;
    t4 = t4 + bihat*bihat';
    posi = posf+1;
end;
phisq0 = (t2-Bhat'*t3)/(n-(m-1)*q-p);
if phisq0 <= 0
    phisq0 = 1;
end;
D0 = (t4-phisq0*t1)/m;
try
    chol(D0);
catch
% Handling non-positive definite initial D;
    [EV,ED] = eig(D0);
    ED(ED<0) = 10^-4;
    D0 = EV*ED*EV^-1;
end
