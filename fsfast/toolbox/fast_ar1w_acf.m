function acf = fast_ar1w_acf(alpha, rho, nf)
% acf = fast_ar1w_acf(alpha, rho, nf)
%
% Autocorrelation function of AR1+White Noise
%   acf(0) = 1
%   acf(n) = (1-alpha)*rho^n


%
% fast_ar1w_acf.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.3 $
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

acf = [];

if(nargin ~= 3)
  fprintf('acf = fast_ar1w_acf(alpha, rho, nf)\n');
  return;quit;
end

nv = size(alpha,2);

nn = repmat([0:nf-1]',[1 nv]);
rho = repmat(rho,[nf 1]);
acf = (rho.^nn) .* repmat((1-alpha),[nf 1]);
acf(1,:) = 1;

return;
