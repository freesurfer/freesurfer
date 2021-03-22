function pdfx = pdf_gamma(x,a,b)
% pdfx = pdf_gamma(x,a,b)
% Gamma distribution evaluated at x for parameters a and b
%
% pdfx = (b.^2) .* (x.^(a-1)) .* exp(-b.*x) ./ gamma(a);
% mean = a/b
% var  = a/(b^2)
% mode = (a-1)/b --> this the value of x at the peak
%
% When creating a hemodynamic response (eg, as in FSL), then
%   a = (tm/sigma)^2
%   b = tm/(sigma^2)
% where sigma = "Stddev" and tm = "Mean lag"
%
%


%
% pdf_gamma.m
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

pdfx = [];

if(nargin ~= 3)
  fprintf('pdfx = pdf_gamma(x,a,b)\n');
  return;
end

nx = length(x);
pdfx = zeros(nx,1);
ind = find(x>0);
pdfx(ind) = (b.^2) .* (x(ind).^(a-1)) .* exp(-b.*x(ind)) ./ gamma(a);

return;



