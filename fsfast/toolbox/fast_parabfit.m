function [par, yhat] = fast_parabfit(x,y)
% [par yhat] = fast_parabfit(x,y)
%
% Fit a parabola.
%
%  y = y0 + m*(x-x0).^2
%
%  Fits for parameters par = [y0 m x0]
%
% $Id: fast_parabfit.m,v 1.1 2003/10/01 05:45:05 greve Exp $

par = [];
yhat = [];

if(nargin ~= 2)
  fprintf('[par yhat] = fast_parabfit(x,y)\n');
  return;
end

% First, fist a linear model
x = reshape1d(x);
nx = length(x);
X = [ones(nx,1) x x.^2];
beta = inv(X'*X)*X'*y;
yhat = X*beta;

% Compute parameters from linear coefficients
m = beta(3);
x0 = -beta(2)/(2*m);
y0 = beta(1)-m*x0^2;
par = [y0 m x0];

return;






