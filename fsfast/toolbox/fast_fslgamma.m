function h = fast_fslgamma(t,meanlag,stddev)
% h = fast_fslgamma(t,<meanlag>,<stddev>)
%
% FSL's Gamma Function
%
% meanlag = 6 sec by default
% stddev  = 3 sec by default
% No normalization is done (yet)
%
% $Id: fast_fslgamma.m,v 1.1 2007/06/19 18:49:38 greve Exp $

h = [];
if(nargin < 1 | nargin > 3)
  fprintf('h = fast_fslgamma(t,meanlag,stddev)\n');
  return;
end

if(~exist('meanlag')) meanlag = 6; end
if(~exist('stddev'))  stddev = 3; end

a = (meanlag/stddev)^2;
b = meanlag/(stddev^2);
h = pdf_gamma(t,a,b);

return;



