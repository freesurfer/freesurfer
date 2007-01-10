function [par, yhat] = fast_parabfit(x,y)
% [par yhat] = fast_parabfit(x,y)
%
% Fit a parabola.
%
%  y = y0 + m*(x-x0).^2
%
%  Fits for parameters par = [y0 m x0]
%
%


%
% fast_parabfit.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

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






