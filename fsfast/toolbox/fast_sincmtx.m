function M = fast_sincmtx(len,delta,rownorm,hw)
% M = fast_sincmtx(len,delta,<rownorm>,<hw>)
%
% Create a sinc interpolation matrix. Given a vector y, yhat = M*y
% will be the sinc interpolation of y shifted backwards by delta,
% where delta is a fration of the distance between y. If rownorm is
% set to non-zero, then each row of the matrix M is normalized so that
% the mean is 1. hw is the width of the hann window. The wider the
% window, the less the effect of the window. Setting hw=0 is like
% setting hw = infinity.
%

% fast_sincmtx.m
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

M = [];

if( nargin < 2 | nargin > 4 )
  fprintf('M = fast_sincmtx(len,delta,<rownorm>)\n');
  return;
end
if(exist('rownorm') ~= 1) rownorm = []; end
if(isempty(rownorm)) rownorm = 0; end

if(exist('hw') ~= 1) hw = []; end
if(isempty(hw))      hw = 0; end

M = zeros(len);
i = 1:len;
for n = 1:len
  x = i - n - delta;
  if(hw > 0) r = sinc(x) .* fast_hann(x,hw);
  else       r = sinc(x);
  end
  if(rownorm) r = r/sum(r); end
  M(n,:) = r;
end

return;

%-----------------------------------------------------------%
function y = sinc(x)

y      = ones(size(x));
inz    = find(x ~= 0);
y(inz) = sin(pi*x(inz)) ./ (pi*x(inz));

return;

