function [px, xbin, xdelta] = samp2pdf(x, b)
% 
% [px xbin xdelta] = samp2pdf(x,b)
%
% Computes the probability distribution function from a sample.
% x is the 1-D vector of samples. b is either the number of bins
% that the pdf will be divided into or a 3-component vector
% [minbin deltabin binmax].
%
% xbin is the vector of
% bin centers, and px is the vector of probabilities that a
% sample will fall into the corresponding bin.  If the number
% of bins is specified, the bins are
% uniformly distributed between the min and max of the sample.
%
% Author: Douglas N. Greve, MGH-NMR
% Date:   12/30/98


%
% samp2pdf.m
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

if( nargin ~= 2 )
  error('USAGE: [px xbin xdelta] = samp2pdf(x,b)');
end
if( isempty(x) )
  error('x is empty matrix');
end

if(length(b) == 1) 
  nbins = b;
  if( nbins < 1 )
    error('nbins must be greater than zero');
  end
  xmax = max(x);
  xmin = min(x);
  xdelta = (xmax-xmin)/nbins;
elseif(length(b) == 3)
  xmin   = b(1);
  xdelta = b(2);
  xmax   = b(3);
else
  error('Second argument must be either the number of bins or range');
end

if( length(x) == 1 )
  warning('only one element in sample vector');
  xbin = x;
  px = 1;
  return;
end


xbin = [xmin : xdelta : xmax];

xHist = hist(x,xbin);

px = xHist/(xdelta*length(x));

return;
