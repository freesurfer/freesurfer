function r = samplecdf(rdim,cdf,xcdf)
% r = samplecdf(rdim,cdf,xpdf);
%
% Samples random numbers from the given cdf. xcdf 
% gives the abscissa for each point on the cdf.
%
%


%
% samplecdf.m
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
