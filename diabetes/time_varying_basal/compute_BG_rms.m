function error_ret = compute_BG_rms(G_target, G_timeseries, error_type)
% error_ret = compute_BG_rms(G_target, G_timeseries, error_type)
% error_type: 0 - RMS
% error_type: 1 - variable
% error_type: 2 - L1

if (nargin < 3)
   error_type = 0 ;
end

switch error_type,
case 0,
   error_ret = rms(G_timeseries - G_target) ;
case 1,
   ind = find(G_timeseries > 140 | G_timeseries < 80) ;
   error_ret = norm(G_timeseries(ind)-G_target) ;
case 2,
     error_ret = sum(abs(G_timeseries-G_target))/length(G_timeseries);
end





