function [slope, offset] = linfit(y,x)
%
% [slope offset] = linfit(y,x)
%
% Solves the equation:
%   y = slope * x + offset
%
%


%
% linfit.m
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

if(nargin ~= 2)
  msg = 'USAGE: [slope offset] = linfit(y,x)';
  qoe(msg);error(msg);
end 

if(length(y) ~= length(x))
  msg = 'Length of x and y must be the same';
  qoe(msg);error(msg);
end

x = reshape1d(x);
y = reshape1d(y);

x2 = [x ones(size(x))];

m = pinv(x2)*y;

slope  = m(1);
offset = m(2);

return;
