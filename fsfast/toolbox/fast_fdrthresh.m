function pthresh = fast_fdrthresh(p,fdr)
% pthresh = fast_fdrthresh(p,fdr)
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
% $Id: fast_fdrthresh.m,v 1.1 2004/10/30 00:36:45 greve Exp $
%

pthresh = [];

if(nargin ~= 2)
  fprintf('pthresh = fast_fdrthresh(p,fdr)\n');
  return;
end

p = sort(abs(p(:)));
Nv = length(p(:));
nn = [1:Nv]';

cVID = 1; % Not sure what this is for

imax = max( find(p <= fdr*nn/Nv ) );
if(~isempty(imax))
  %fprintf('imax = %d\n',imax);
  pthresh = p(imax);
else
  pthresh = min(p);
end

return;