function [y, H] = randnorm(m, n, Sigma)
%
% [y H] = randnorm(m, n, Sigma)
%
% Creates an m-by-n matrix y whose columns are
% distributed N(0,Sigma).  H is the cholesky 
% decomposition of Sigma. y = H * randn(m,n).
%
% $Id: randnorm.m,v 1.1 2003/03/04 20:47:41 greve Exp $

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
