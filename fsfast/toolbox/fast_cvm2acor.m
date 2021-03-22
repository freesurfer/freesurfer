function acor = fast_cvm2acor(cvm,biased) 
%
% acor = fast_cvm2acor(cvm) 
%
% Computes the autocorrelation function from the 
% covariance matrix


%
% fast_cvm2acor.m
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

acor = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('USAGE: acor = fast_cvm2acor(cvm,<biased>)\n');
  return;
end

if(exist('biased') ~= 1) biased = []; end
if(isempty(biased)) biased = 0; end

ntrs = size(cvm,1);

acor = zeros(ntrs,1);
for n = 1:ntrs, 
  acor(n) = mean(diag(cvm,n-1));
end

acor = acor/acor(1);

if(biased)
  acor = acor.*([ntrs:-1:1]');
  acor = acor/max(acor);
end


return;
