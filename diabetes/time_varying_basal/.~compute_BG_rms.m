function rms_error = compute_BG_rms(G_target, G_timeseries, error_type)
% rms_error = compute_BG_rms(G_target, G_timeseries, error_type)
% error_type: 0 - RMS
% error_type: 1 - variable

if (nargin < 3)
   error_type = 0 ;
end

if (error_type == 0)
   rms_error = rms(G_timeseries - G_target) ;
end

if (error_type == 1)
   ind = find(G_timeseries > 140 | G_timeseries < 80) ;
   rms_error = norm(G_timeseries(ind)-G_target) ;
end


