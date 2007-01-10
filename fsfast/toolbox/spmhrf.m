function h = spmhrf(t,a1,b1,a2,b2,c)
%
% h = spmhrf(t,a1,b1,a2,b2,c)
%


%
% spmhrf.m
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

h = [];

if(nargin ~= 1 & nargin ~= 6)
  fprintf('USAGE: h = spmhrf(t,a1,b1,a2,b2,c)\n');
  return;
end

if(nargin == 1)
  a1 = 6;
  b1 = 0.9;
  a2 = 12;
  b2 = 0.9;
  c  = .35;
end

d1 = a1*b1;
d2 = a2*b2;

h =    (t/d1).^a1  .* exp(-(t-d1)/b1) - ...
    c*((t/d2).^a2) .* exp(-(t-d2)/b2);

return;



