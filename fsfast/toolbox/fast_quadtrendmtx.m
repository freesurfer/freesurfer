function Xqtrend = fast_quadtrendmtx(run,ntrs,nruns)
% Xqtrend = fast_quadtrendmtx(run,ntrs,nruns)
%
% Quadratic trend - centered half-way into the run; mean zero.
% This vector will be orthogonal to those produced by
% fast_baselinemtx fast_trendmtx 


%
% fast_quadtrendmtx.m
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

if(nargin ~= 3)
  msg = 'USAGE: Xqtrend = fast_quadtrendmtx(run,ntrs,nruns)';
  qoe(msg);error(msg);
end

t = [0:ntrs-1]'; %'
t = t - mean(t);
v = t.^2;
v = v - mean(v);
v = v./sqrt(sum(v.^2));

Xqtrend        = zeros(ntrs,nruns);
Xqtrend(:,run) = v;

return;
