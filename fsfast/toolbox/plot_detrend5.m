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

ntrs = length(t);
Xdrift  = fast_polytrendmtx(1,ntrs,1,5);
Tdrift = Xdrift*inv(Xdrift'*Xdrift)*Xdrift';
Edrift = eye(ntrs) - Tdrift;
y2 = Edrift*y;
%hplot = plot(t,y-mean(y),t,y2);
hplot = plot(t,y2);

return;
