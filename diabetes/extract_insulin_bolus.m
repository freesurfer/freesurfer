function insulin_time_series = extract_insulin_bolus(S, sind, time_series)
% insulin_time_series = extract_insulin(S, sind, time_series)
% convert insulin bolus inputs to a nx2 array
% col 1 is the insulin given (0 if none)
% col 2 is the time since the last insulins
insulin_input_ind = 13 ;

insulin_vec = cell2mat(S(sind).rows(:,insulin_input_ind)) ;
ntps = length(time_series) ;

insulin_time_series = zeros(ntps,2);

last_insulin_time = -1 ;
for t=1:ntps
    if (isnan(insulin_vec(t)))	
         insulin_time_series(t,1) = 0 ;
       	  if (last_insulin_time > 0)
	     insulin_time_series(t,2) = time_series(t) - last_insulin_time ;
	  else
	     insulin_time_series(t,2) = -1 ;
          end
    else  % a real insulin values
    	  insulin_time_series(t,1) = insulin_vec(t) ;
	  insulin_time_series(t,2) = 0 ;   % time since last insulins 
	  last_insulin_time = time_series(t) ;	  
    end
end 
