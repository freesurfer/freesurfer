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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
