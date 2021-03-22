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

if(nargin ~= 3) 
  msg = 'USAGE: y = gaussian(x,mean,stddev)'
  error(msg);
end

y = exp( -((x-mean).^2)./(2*(stddev.^2))) ./ (stddev * sqrt(2*pi));

return
