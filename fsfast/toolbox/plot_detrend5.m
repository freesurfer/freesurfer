function hplot = plot_detrend5(t,y)
% hplot = plot_detrend5(t,y)
%
% This function plots y2 vs t, where y2 is the
% detrended version of y. y is detrended with
% a 5th order polynomial. This can be used in
% conjunction with yakview by passing plot_detrend5
% with the -rf flag.

ntrs = length(t);
Xdrift  = fast_polytrendmtx(1,ntrs,1,5);
Tdrift = Xdrift*inv(Xdrift'*Xdrift)*Xdrift';
Edrift = eye(ntrs) - Tdrift;
y2 = Edrift*y;
%hplot = plot(t,y-mean(y),t,y2);
hplot = plot(t,y2);

return;
