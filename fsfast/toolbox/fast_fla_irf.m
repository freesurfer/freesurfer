function [irf, irftaxis] = fast_fla_irf(fla,beta)
% [irf irftaxis] = fast_fla_irf(fla,beta)
%
% Returns the Impulse Response Function (IRF) for the current
% fx (ie, nthfx) and the time axis for the IRF. Uses beta
% to weight the components of the assumption matrix. The
% nthfx must be an ERM.
%
% $Id: fast_fla_irf.m,v 1.1 2003/08/21 03:43:03 greve Exp $

irf = [];
irftaxis = [];

if(nargin ~= 2)
  fprintf('[irf irftaxis] = fast_fla_irf(fla,beta)\n');
  return;
end

if(fla.nthfx > length(fla.fxlist))
  fprintf('fla_irf: nthfx=%d too big\n',fla.nthfx);
  return;
end

if(~fast_fxcfg('iserm',fla))
  fprintf('fla_irf: nthfx=%d is not an ERM\n',fla.nthfx);
  return;
end

fx = fla.fxlist(fla.nthfx).fx;
A = fast_fxcfg('irfmatrix',fla);

if(strcmp(fx.fxtype,'fixed'))
  ind = fx.regind{1};
else
  ind = fx.regind{fla.nthrun};
end

irf = A*beta(ind,:);
irftaxis = fast_fxcfg('irftaxis',fla);

return;