function active_insulin_time_series = extract_active_insulin(S, sind,time_series, bolus_vec, basal_vec,iob_table)
% active_insulin_time_series = extract_active_insulin(S, sind,time_series, bolus_vec, basal_vec,iob_table)


ntps = length(time_series) ;

active_insulin_time_series = zeros(ntps,1);


% track the start time of the bolus and the amount for caculating
% current active
max_boluses = 1000;
active_boluses = zeros(max_boluses, 2) ;
num_active_boluses = 0 ;

for t=1:ntps
    if (isnan(time_series(t)))
            continue ;
    end
    % zero out boluses that are now over
    for bolus_ind = 1:num_active_boluses
      deltat = round(time_series(t) - active_boluses(bolus_ind,2) + 1) ;
      bolus_over = (deltat > length(iob_table)) ;
      if (bolus_over == 0)
      	 bolus_over = (iob_table(deltat) == 0) ;
      end
      if (bolus_over)
       	  active_boluses(bolus_ind,1) = 0 ;
       	  active_boluses(bolus_ind,2) = 0 ;
      end
    end
    insulin = basal_vec(t)/60 + bolus_vec(t) ;
    if (insulin > 0)
       for bolus_ind = 1:num_active_boluses
       	   if (active_boluses(bolus_ind,1) == 0)   % found an empty slot
	      current_ind = bolus_ind ;
	      break ;
	   end
       end    
       if (num_active_boluses == 0) 
         current_ind = 1 ;
       end
       if (current_ind >= num_active_boluses)
         num_active_boluses = num_active_boluses + 1 ;
         current_ind = num_active_boluses ;
       end
       active_boluses(current_ind,1) = insulin ;
       active_boluses(current_ind,2) = time_series(t) ;
    end

    % now add up all the active boluses  to compute the amount of active insulin
    active_insulin = 0 ;
    for bolus_ind = 1:num_active_boluses
        deltat = round(time_series(t) - active_boluses(bolus_ind,2) + 1) ;
	bolus = active_boluses(bolus_ind,1) ;
	if (bolus > 0)  % otherwise deltat might be out of bounds
	    	active_insulin = active_insulin + bolus * iob_table(deltat)/100 ;
        end
    end
    active_insulin_time_series(t,1) = active_insulin ;
end 


