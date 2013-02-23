function pthresh = lme_mass_FDR(p,fdr)
% pthresh = lme_mass_FDR(p,fdr)
% This function has slightly modified the freesurfer's fast_fdrthresh
% function (jbernal modification 2010)
%
% p = list of p values between -1 and 1
% fdr = false discovery rate, between 0 and 1
%
% Based on Tom's FDR.m from 
%   http://www.sph.umich.edu/~nichols/FDR/FDR.m
% The threshold returned from this function is based on an 
% assumption of "independence or positive dependence",
% which should be "reasonable for imaging data".
%
% $Id: lme_mass_FDR.m,v 1.1.2.2 2013/02/23 21:08:10 nicks Exp $
%

if(nargin ~= 2)
  fprintf('pthresh = lme_mass_FDR(p,fdr)\n');
  return;
end
% pthresh = [];
p = sort(abs(p(:)));
Nv = length(p(:));
nn = [1:Nv]';
imax = max(find(p <= fdr*nn/Nv));
if(~isempty(imax))
  %fprintf('imax = %d\n',imax);
  pthresh = p(imax);
else
    %This is just to not return the min(p) in this case
    pthresh = min(p)/10;
end
return;