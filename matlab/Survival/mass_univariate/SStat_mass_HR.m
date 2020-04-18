function [hr,pval] = SStat_mass_HR(x1,x2,y1,y2,stats)
% [hr,pval] = SStat_mass_HR(x1,x2,y1,y2,stats)
%
% Vertex/voxel-wise hazard ratio estimates for Cox models.
%
% Input
%
% x1: Row vector with the covariate values for the first group (1xp). 
% x2: Row vector with the covariate values for the second group (1xp) .
% y1: Row vector with the data for the first group (1xnv, nv=#vertices);
% y2: Row vector with the data for the second group (1xnv);
% stats: Structure array containing statistics obtained with any of 
% SStat_mass_CoxPH,SStat_mass_CoxStratPH and SStat_mass_CoxExt.
%
% Output
% hr: Hazard ratio values.
% pval: P-values for the hr values.
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Kleinbaum, D.G., Klein, M., 2005. Survival analysis. A self-
% learning approach, second edition. New York: Springer..
%   
if nargin < 5
    error('Too few inputs');
end;
if (length(x1)~=length(x2))
    error('Vectors x1 and x2 must have the same length.');
end;
nv = length(stats);
hr = zeros(1,nv);
pval = ones(1,nv);
for j=1:nv
    if ~isempty(stats(j).Bhat)
        if length(stats(j).Bhat)==length([x1 y1(j)])
            [hr(j),pval(j)] = SStat_HR([x1 y1(j)],[x2 y2(j)],stats(j));
        else
           error(['stats(' num2str(j) ').Bhat must have the same length as'...
            ' vectors x1 and x2']);
        end;
    end;
end;
