function [g,x] = fast_gaussian(nmean,nstddev,len,twosided)
% [g,x] = fast_gaussian(nmean,nstddev,len,<twosided>)


%
% fast_gaussian.m
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
