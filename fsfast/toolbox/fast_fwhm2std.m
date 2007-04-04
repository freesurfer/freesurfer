function gstd = fast_fwhm2std(fwhm)
% gstd = fast_fwhm2std(fwhm)
% gstd = fwhm/sqrt(log(256.0))
%
% $Id: fast_fwhm2std.m,v 1.1 2007/04/04 01:59:36 greve Exp $

gstd = [];

if(nargin ~= 1)
  fprintf('gstd = fast_fwhm2std(fwhm)\n');
  fprintf('gstd = fwhm/sqrt(log(256.0))\n');
  return;
end

  
gstd = fwhm/sqrt(log(256.0))

return;






