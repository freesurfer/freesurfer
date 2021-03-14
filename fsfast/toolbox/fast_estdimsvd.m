function [dim, pvs, pvsw] = fast_estdimsvd(s,pvsthresh)
% [dim, pvs, pvsw] = fast_estdimsvd(s,pvsthresh)
%
% Estimates the dimension of a data set from the
% s matrix of the SVD (y = u*s*v') based on where
% the Percent Variance Spanned (PVS) (ie, eigenspectrum)
% crosses that of white noise.
%


%
% fast_estdimsvd.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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

















