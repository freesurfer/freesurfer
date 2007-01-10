function [g,x] = fast_gaussian(nmean,nstddev,len,twosided)
% [g,x] = fast_gaussian(nmean,nstddev,len,<twosided>)


%
% fast_gaussian.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

g = [];
x = [];

if(nargin ~= 3 & nargin ~= 4)
  fprintf('[g,x] = fast_gaussian(nmean,nstddev,len,<twosided>)\n');
  return;
end

if(exist('twosided') ~= 1) twosided = 0; end

nvar = nstddev.^2;
f = 1/sqrt(2*pi*nvar);

if(0)
  x = ([1:len] - nmean);
  g = f * exp ( -(x.^2)/(2*nvar) );
else
  x = ([1:len] - nmean)/nstddev;
  g = f * exp ( -(x.^2)/2 );
end

if(twosided)
  grev = fliplr(g(2:end));
  g = [grev g];
  x = [-fliplr(x(2:end)) x];
end


return
