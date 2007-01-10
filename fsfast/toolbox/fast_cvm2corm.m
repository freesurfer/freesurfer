function cormtx = fast_cvm2corm(cvm)
% cormtx = fast_cvm2corm(cvm)
% Converts a covariance matrix into a correlation coefficient
% matrix. It computes the correlation coefficient for compoent
% i,j as cvm(i,j) / sqrt(cvm(i,i) * cvm(j,j))


%
% fast_cvm2corm.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
