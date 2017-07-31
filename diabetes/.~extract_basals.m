function basal_time_series = extract_basals(S, sind, time_series)
% basal_time_series = extract_basals(S, sind, time_series)

define_medtronic_columns

ntps = length(time_series) ;

basal_time_series = zeros(ntps,1);

basal = 0 ;
for t=1:ntps
    if (strcmp(char(S(sind).rows(t,raw_type_ind)), 'BasalProfileStart') == 1)
       str = char(S(1).rows(2,raw_val_ind));
       basal = sscanf(str, 'PATTERN_NAME=standard, PROFILE_INDEX=0, RATE=%f') ;
    end
    if (t > 1 && time_series(t,1) == time_series(t-1,1))
        basal_time_series(t,1) = 0 ;   % basal already entered for this t
     else
	if (isnan(time_series(t)))
	  basal_time_series(t,1) = 0 ;
	else
	  basal_time_series(t,1) = basal ;
	end
     end
     basal_time_series(t,1) = basal ;
end 


temp_basal_amount_vec = cell2mat(S(sind).rows(:,temp_basal_amount_ind)) ;
temp_basal_duration_vec = cell2mat(S(sind).rows(:,temp_basal_duration_ind)) ;
for t=1:ntps
    if (isnan(temp_basal_amount_vec(t)) == 0)
       temp_basal = temp_basal_amount_vec(t) ;
       temp_basal_duration = temp_basal_duration_vec(t) ;
       disp(sprintf('temp basal time index %d', t));
       disp(sprintf('new temp basal %2.3f U/hr detected at time %d for %d min', temp_basal,time_series(t),temp_basal_duration)) ;
%        if (strcmp(char(S(sind).rows(t,temp_basal_type_ind)), 'Percent') == 1)
%       str = char(S(1).rows(2,raw_val_ind));
%       basal = sscanf(str, 'PATTERN_NAME=standard, PROFILE_INDEX=0, RATE=%f') ;
    end
%    basal_time_series(t,1) = basal ;
end 

