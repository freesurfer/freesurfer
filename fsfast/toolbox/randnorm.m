function [y, H] = randnorm(m, n, Sigma)
%
% [y H] = randnorm(m, n, Sigma)
%
% Creates an m-by-n matrix y whose columns are
% distributed N(0,Sigma).  H is the cholesky 
% decomposition of Sigma. y = H * randn(m,n).
%
%


%
% randnorm.m
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

if(nargin ~= 3) 
  msg = 'USAGE: [y H] = randnorm(m, n, Sigma)';
  qoe(msg); error(msg);
end

q = size(Sigma,1) ; 
if( q > m)
  tmp = Sigma(1:m,1:m);
  Sigma = tmp;
  clear tmp;
elseif(q < m)
  tmp = zeros(m,m);
  tmp(1:q,1:q) = Sigma;
  Sigma = tmp;
  clear tmp;
end

H = chol(Sigma);

y = H * randn(m,n);

return ;
