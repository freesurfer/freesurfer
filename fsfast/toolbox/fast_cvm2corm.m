function cormtx = fast_cvm2corm(cvm)
% cormtx = fast_cvm2corm(cvm)
% Converts a covariance matrix into a correlation coefficient
% matrix. It computes the correlation coefficient for compoent
% i,j as cvm(i,j) / sqrt(cvm(i,i) * cvm(j,j))


%
% fast_cvm2corm.m
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

cormtx = [];

if(nargin ~= 1)
  msg = 'USAGE: cormtx = fast_cvm2corm(cvm)';
  fprintf('%s\n',msg);
  return;
end

d = diag(cvm); % n X 1
n = length(d);

m = sqrt(repmat(d, [1 n]) .* repmat(d', [n 1])); %'

cormtx = cvm ./ m;

return;
