function [pdf, pdfvar] = pdf_rayleigh(x,mean)
% [pdf, pdfvar] = pdf_rayleigh(x,mean)
%
% Rayleigh Probability distribution/density function sampled at
% x. If mean is not specified, it defaults to 1. The variance
% is dependent on the mean and is (4/pi-1)*mean.^2, which is
% returned in pdfvar. 
%
% Rayleigh dist noise can be constructed from complex white noise:
%   y = abs(randn(10000,1) +i*randn(10000,1));
%   The var will be the var of the white noise times (2-pi/2)
%
% $Id: pdf_rayleigh.m,v 1.2 2004/09/30 19:49:22 greve Exp $

pdf = [];
if(nargin < 1 | nargin > 2)
  fprintf('[pdf pdfvar] = pdf_rayleigh(x,<mean>)\n');
  return;
end

if(~exist('mean','var')) mean = []; end
if(isempty(mean)) mean = 1; end

pdfvar = (4/pi-1)*(mean.^2);

alpha = mean/sqrt(pi/2);

pdf = x .* exp(-(x.^2)./(2*(alpha.^2))) ./ (alpha.^2);

return;