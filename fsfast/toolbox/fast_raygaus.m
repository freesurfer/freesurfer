function [p,pR,pG] = fast_raygaus(y,alpha,Rmu,Gmu,Gstd)
% [p pR pG] = raygaus(y,alpha,Rmu,Gmu,Gstd)
% [p pR pG] = raygaus(y,params)
%
% Computes the probability of getting y from
% a Rayleigh-Gaussian mixture.
%
% alpha - fraction in rayleigh
% Rmu - mean of rayleigh
% Gmu - mean of gaussian
% Gstd - std of gaussian
%
% params = [alpha Rmu Gmu Gstd];
%
% p = alpha*pR + (1-alpha)*pG
%
% See also: pdf_rayleigh, gaussian, fast_raygausinit, 
%   fast_raygauscost.
%
% $Id: fast_raygaus.m,v 1.1 2006/04/25 01:31:40 greve Exp $

p = [];
if(nargin ~= 2 & nargin ~= 5)
  fprintf('[p pR pG] = fast_raygaus(y,alpha,Rmu,Gmu,Gstd)\n');
  fprintf('[p pR pG] = fast_raygaus(y,params)\n');
  return;
end

if(nargin == 2)
  params = alpha;
  alpha = params(1);
  Rmu   = params(2);
  Gmu   = params(3);
  Gstd   = params(4);
end

pR = pdf_rayleigh(y,Rmu);
pG = gaussian(y,Gmu,Gstd);

p = alpha*pR + (1-alpha)*pG;

return;
