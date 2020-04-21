function gstd = fast_fwhm2std(fwhm)
% gstd = fast_fwhm2std(fwhm)
% gstd = fwhm/sqrt(log(256.0))
%

gstd = [];

if(nargin ~= 1)
  fprintf('gstd = fast_fwhm2std(fwhm)\n');
  fprintf('gstd = fwhm/sqrt(log(256.0))\n');
  return;
end

  
gstd = fwhm/sqrt(log(256.0))

return;






