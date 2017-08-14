function carb_time_series = extract_carbs(S, sind, time_series)
% carb_time_series = extract_carbs(S, sind, time_series)
% convert carb inputs to a nx2 array
% col 1 is the carbs given (0 if none)
% col 2 is the time since the last carbs
carb_input_ind = 24 ;

carb_vec = cell2mat(S(sind).rows(:,carb_input_ind)) ;
ntps = length(time_series) ;

carb_time_series = zeros(ntps,2);

last_carb_time = -1 ;
for t=1:ntps
    if (isnan(carb_vec(t)) | (carb_vec(t) == 0))
         carb_time_series(t,1) = 0 ;
       	  if (last_carb_time > 0)
	     carb_time_series(t,2) = time_series(t) - last_carb_time ;
	  else
	     carb_time_series(t,2) = -1 ;
          end
    else  % a real carb values
    	  carb_time_series(t,1) = carb_vec(t) ;
	  carb_time_series(t,2) = 0 ;   % time since last carbs 
	  last_carb_time = time_series(t) ;	  
	  for (t1=max(t-1,1):-1:max(t-20,1))
	      if (time_series(t1) == time_series(t))
	      	 carb_time_series(t1,2) = 0 ;
	      end
          end
    end
end 
