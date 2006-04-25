function cost = fast_raygauscost(params,y)
% cost = fast_raygauscost(p,y)
% params = [alpha Rmu Gmu Gstd];
%
% Cost function for Rayleigh-Gaussian Mixture Model
%
% [p pR pG] = fast_raygaus(y,alpha,Rmu,Gmu,Gstd);
% cost = -sum(log10(p));
%
% optparams = fminsearch('fast_raygauscost',params,[],y);
%
% See also fast_raygaus, fast_raygauscost.
%
% $Id: fast_raygauscost.m,v 1.1 2006/04/25 01:31:40 greve Exp $

alpha = params(1);
Rmu   = params(2);
Gmu   = params(3);
Gstd   = params(4);

if(alpha < 0 | alpha > 1 | Rmu < 0 | Gmu < 0 | Gstd < 0)
  cost = 10e10;
  return;
end

p = fast_raygaus(y,alpha,Rmu,Gmu,Gstd);
ind = find(p == 0);
if(~isempty(ind))
  cost = 10e10;
  return;
end

cost = -sum(log10(p));

return;










