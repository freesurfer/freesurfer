function y = gaussian(x,mean,stddev)
%
% Gaussian Distribution Function
%
% y = gaussian(x,mean,stddev)
%
% y = exp( -((x-mean).^2)/(2*(stddev.^2))) / (stddev * sqrt(2*pi));
% y = y/sum(y);


%
% gaussian.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.4 $
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

if(nargin ~= 3) 
  msg = 'USAGE: y = gaussian(x,mean,stddev)'
  error(msg);
end

y = exp( -((x-mean).^2)./(2*(stddev.^2))) ./ (stddev * sqrt(2*pi));

return
