function [dim, pvs, pvsw] = fast_estdimsvd(s,pvsthresh)
% [dim, pvs, pvsw] = fast_estdimsvd(s,pvsthresh)
%
% Estimates the dimension of a data set from the
% s matrix of the SVD (y = u*s*v') based on where
% the Percent Variance Spanned (PVS) (ie, eigenspectrum)
% crosses that of white noise.
%

if(nargin ~= 1 & nargin ~= 2)
  fprintf('dim = fast_estdimsvd(s,<pvsthresh>)\n');
  return;
end

if(nargin ~= 2) pvsthresh = 0; end

nf = size(s,1);
ds = diag(s);
pvs = 100*ds/sum(ds);

% Simulate white noise process %
w = randn(nf,10*nf);
Mw = w*w'; 
[uw sw blah] = svd(Mw);
dsw = diag(sw);
pvsw = 100*dsw/sum(dsw);

% This is the difference in the eigen spectra
d = 100*(pvs-pvsw)./pvsw;

dim = max(find(d > pvsthresh));

%nn = 2:20;
%plot(nn,pvsw(nn),nn,pvs(nn));
%keyboard

return;

















