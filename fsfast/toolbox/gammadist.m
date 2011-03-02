function y = gammadist(x,avg,stddev)
% Gamma Distribution
%
% y = gammadist(x,avg,stddev)
%
% This probably does not work.


%
% gammadist.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:07 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
  msg = 'USAGE: y = gammadist(x,avg,stddev)'
  error(msg);
end

z = (x-avg)/stddev;
y = (z.^2) .* exp(-z);
ind = find(z<0);
y(ind) = 0;
y = y/sum(y);

return
