function [hbar, htxt] = barsegstats(segnames,segstats)
%  [hbar htxt] = barsegstats(segnames,segstats)
%
% [segnames segindex segstats] = load_segstats('aseg.stats','bert');
% [hbar htxt] = barsegstats(segnames,segstats(:,3))
% 
% 


%
% barsegstats.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/02/15 23:51:41 $
%    $Revision: 1.3 $
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


mn = min(segstats(:));
mn = 0;
mx = max(segstats(:));
d = mx-mn;
a = mn+.02*d;

nseg = size(segnames,1);
x = 1:nseg;
hbar = bar(segstats,'r');
for n = 1:nseg
  htxt(n) = text(n+.25,a,segnames(n,:));
  set(htxt(n),'rotation',90);
  set(htxt(n),'fontsize',4)
end


return;
