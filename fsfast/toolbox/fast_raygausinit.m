function params = fast_raygausinit(y)
% params = fast_raygausinit(y)
% params = [alpha Rmu Gmu Gstd];
%
%


%
% fast_raygausinit.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
%    $Revision: 1.4 $
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

if(nargin ~= 1)
  fprintf('params = fast_raygausinit(y)\n');
  return;
end

if(0)
%md = mode(y);
ysort = sort(y);
ny = length(y);
md = ysort(round(ny/2));
indray = find(y < md);
indgaus = find(y >= md);
end

ymn = mean(y);
indray  = find(y < ymn);
indgaus = find(y > ymn);

alpha = .25;
Rmu = mean(y(indray));
Gmu = mean(y(indgaus));
Gstd = std(y(indgaus));

params = [alpha Rmu Gmu Gstd];

return;
