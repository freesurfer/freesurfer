function hplot = plot_detrend5(t,y)
% hplot = plot_detrend5(t,y)
%
% This function plots y2 vs t, where y2 is the
% detrended version of y. y is detrended with
% a 5th order polynomial. This can be used in
% conjunction with yakview by passing plot_detrend5
% with the -rf flag.


%
% plot_detrend5.m
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

ntrs = length(t);
Xdrift  = fast_polytrendmtx(1,ntrs,1,5);
Tdrift = Xdrift*inv(Xdrift'*Xdrift)*Xdrift';
Edrift = eye(ntrs) - Tdrift;
y2 = Edrift*y;
%hplot = plot(t,y-mean(y),t,y2);
hplot = plot(t,y2);

return;
