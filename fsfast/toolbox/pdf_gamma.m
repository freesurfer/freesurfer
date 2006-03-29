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
% $Id: pdf_gamma.m,v 1.2 2006/03/29 23:22:38 greve Exp $

pdfx = [];

if(nargin ~= 3)
  fprintf('pdfx = pdf_gamma(x,a,b)\n');
  return;
end

pdfx = (b.^2) .* (x.^(a-1)) .* exp(-b.*x) ./ gamma(a);

return;



