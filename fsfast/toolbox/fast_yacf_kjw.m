function [yacf, M] = fast_yacf_kjw(racf,R)
% yacf = fast_yacf_kjw(racf,R)
% 
% Remove bias from residual autocorrelation function using Keith
% Worsley's fmristats method. The bias is induced by projecting out
% the task components to form the residual.
%
% racf is the autocor function of the residual. This
%  can be nf X nv. See fast_acorr.m
% R is the residual forming matrix
% yacf is an approximation of the autocor function
%  of the original noise.
% M is the correction matrix.
%
% Notes:
%  1. As the number of frames increases, the computation is 
%     more intense and the correction matrix is more singular.
%     There is some happy medium where the accuracy of the 
%     ACF is best, but it depends upon the design among other
%     things.
%
% See also: fast_acorr, fast_kjw_mtx.
%
% $Id: fast_yacf_kjw.m,v 1.4 2004/04/06 16:22:59 greve Exp $

if(nargin ~= 2)
  fprintf('yacf = fast_yacf_kjw(racf,R)\n');
  return;
end

[nf nv] = size(racf);
M = fast_kjw_mtx(R,nf);

yacf = inv(M)*racf;
yacf = yacf./repmat(yacf(1,:),[nf 1]);

return;
