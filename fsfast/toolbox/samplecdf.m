function r = samplecdf(rdim,cdf,xcdf)
% r = samplecdf(rdim,cdf,xpdf);
%
% Samples random numbers from the given cdf. xcdf 
% gives the abscissa for each point on the cdf.
%
% $Id: samplecdf.m,v 1.1 2004/02/12 18:43:31 greve Exp $

r = [];

if(nargin ~= 3)
  fprintf('r = samplecdf(rdim,cdf,xpdf)\n');
  return;
end

nr = prod(rdim);
r = zeros(nr,1);

for n = 1:nr
  u = rand;
  [m i] = min(abs(cdf-u));
  r(n) = xcdf(i);
end

r = reshape(r,[rdim 1]);


return;