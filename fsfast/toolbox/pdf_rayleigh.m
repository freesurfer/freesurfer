function [pdf, pdfvar] = pdf_rayleigh(x,mean,log10flag)
% [pdf, pdfvar] = pdf_rayleigh(x,<mean>,<log10flag>)
%
% Rayleigh Probability distribution/density function sampled at
% x. If mean is not specified, it defaults to 1. The variance
% is dependent on the mean and is (4/pi-1)*mean.^2, which is
% returned in pdfvar. 
%
% If log10flag is not 0, then -log10 of the pdf is returned, ie,
%  pdf_rayleigh(x) = 10^-(pdf_rayleigh(x,[],1))
%
% Rayleigh dist noise can be constructed from complex white noise:
%   y = abs(randn(10000,1) +i*randn(10000,1));
%   The var will be the var of the white noise times (2-pi/2)
%
%


%
% pdf_rayleigh.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

pdf = [];
if(nargin < 1 | nargin > 3)
  fprintf('[pdf pdfvar] = pdf_rayleigh(x,<mean>,<log10flag>)\n');
  return;
end

if(~exist('mean','var')) mean = []; end
if(isempty(mean)) mean = 1; end

if(~exist('log10flag','var')) log10flag = 0; end
if(isempty(log10flag)) log10flag=0; end

pdfvar = (4/pi-1)*(mean.^2);
alpha = mean/sqrt(pi/2);

if(~log10flag)
  pdf = x .* exp(-(x.^2)./(2*(alpha.^2))) ./ (alpha.^2);
else
  pdf = log10(x) + -(x.^2).*(log10(exp(1))./(2*(alpha.^2))) - log10(alpha.^2);
  pdf = -pdf; % Actually -log10
end

return;
