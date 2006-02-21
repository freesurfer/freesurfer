function [p,p2] = fast_trf_sup_surf(t,dof,fwhm,area)
% [p p2] = fast_trf_sup_surf(t,dof,fwhm,area)
%
% Based on Moo Chung, et al, Cortical thickness analysis in autism
% with heat kernel smoothing. NeuroImage, 25, 1256-1265. 2005.
%
% $Id: fast_trf_sup_surf.m,v 1.1 2006/02/21 01:30:52 greve Exp $

p = [];
if(nargin ~= 4)
  fprintf('[p p2] = fast_trf_sup_surf(t,dof,fwhm,area)\n');
  return;
end

p0 = tTest(dof,t)/2;

gstd = fwhm/sqrt(log(256));
lambda = 1/(2*gstd*gstd);

tmp = t.*(1 + (t.^2)/dof).^(-(dof-1)/2);

f = (lambda*gamma((dof+1)/2)) /...
    ( ((2*pi)^(3/2)) * sqrt(dof/2) * gamma(dof/2) );
p2 = f*tmp;

p = 2*p0 + area*p2/2;
%p = p0;

return;















