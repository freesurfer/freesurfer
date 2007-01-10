function G = fast_mkgausmtx(nstddev,len)
% G = fast_mkgausmtx(nstddev,len)
% Create a len-by-len guassian filter matrix 


%
% fast_mkgausmtx.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

G = [];

if(nargin ~= 2)
  fprintf('G = fast_mkgausmtx(nstddev,len)\n');
  return;
end

for n = 1:len
  g = fast_gaussian(n,nstddev,len);
  % g = g/sqrt(sum(g.^2));
  g = g/sum(abs(g));
  G = [G; g];
end

return;
