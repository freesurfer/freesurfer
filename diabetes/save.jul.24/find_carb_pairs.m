min_time_for_metabolized_carbs = 190 ;
max_time_for_metabolized_carbs = 4*60; ;
ind = 1 ;
clear('M', 'F', 'dt') ;
carb_time = -1 ;
ind = 1 ;
for t=1:ntps

    if (t > 1 & ((time_series(t) - time_series(t-1))>30))
       carb_time = -1 ;   % some period of no measurements, don't use for training
    end
    if ((carb_vec(t,1) > 0) & (isnan(sensor_vec(t)) == 0))
        carbs = carb_vec(t,1) ;
        carb_time = t ;
    end

    end_time = 0 ;
    if (carb_vec(t,2) > min_time_for_metabolized_carbs & carb_time > 0)
       if (t == ntps || carb_vec(t,2) >= max_time_for_metabolized_carbs)
       	  end_time = t;
       elseif  ((carb_vec(t+1,2) < carb_vec(t,2)) & (isnan(sensor_vec(t)) ==0)) % next time is new carb infusion
       	  end_time = t;
       elseif (isnan(bg_vec(t)) == 0)
       	  end_time = t ;
	end
    end

    % if there was a meter measurement in the middle of the time
    %  period then don't use this pair as the subject may have taken
    % some action based on it (like giving carbs) that isn't recorded
    if ((carb_vec(t,2) > 5) & (carb_vec(t,2) < min_time_for_metabolized_carbs) & (carb_time > 0) & (bg_vec(t) > 0))
%       disp(sprintf('bg_vec(%d) = %2.0f', t, bg_vec(t)));
       carb_time = -1 ;
    end

    if (end_time)
        % compute total insulin absobed in that time
        F(ind) = carbs ;

	% use meter reading if available
	meter_reading = -1 ;
	for t1=max(1,carb_time-20):min(ntps,carb_time+20)
	    if ((abs(time_series(t1) - time_series(carb_time))<5) & (isnan(bg_vec(t1)) == 0))
	       meter_reading = bg_vec(t1) ;
	       break ;
	    end
	end
	if (meter_reading > 0)
	   M(ind,1) = meter_reading ;
	   measured(ind,1) = 1 ;
	else
	   M(ind,1) = sensor_vec(carb_time) ;
	   measured(ind,1) = 0 ;
        end
	meter_reading = -1 ;
	for t1=max(1,carb_time-20):min(ntps,carb_time+20)
	    if ((abs(time_series(t1) - time_series(carb_time))<5) & (isnan(bg_vec(t1)) == 0))
	       meter_reading = bg_vec(t1) ;
	       break ;
	    end
	end
	if (meter_reading > 0)
	   measured(ind,2) = 1 ;
           M(ind,2) = meter_reading ;
	else
	   measured(ind,2) = 0 ;
	   M(ind,2) = sensor_vec(t) ;
	end
        insulin = 0 ;
        for t1=max(1,carb_time-20):min(ntps,t+10)
            if (time_series(t1) < time_series(carb_time))  % still looking for start of this time
                continue ;
            end
            if (time_series(t1) > time_series(t))
                break ;
            end
            deltat = round(time_series(t) - time_series(t1) + 1) ;
            if (deltat > length(iob_table))
                pct_active = 0 ;
            else
                pct_active = iob_table(deltat)/100 ;
            end
            insulin = insulin + bolus_vec(t1,1)*(1-pct_active) ;
            if (bolus_vec(t1,1) > 0)
%                disp(sprintf('adding insulin vec at t=%2.0f: %2.2f = %%2.2f', t1,bolus_vec(t1,1), insulin));
            end
        end

	start_times(ind) = carb_time ;
	end_times(ind) = t ;	        
        D(ind) = insulin ;
        dt(ind) = (time_series(t) - time_series(carb_time) )/60 ;
%        disp(sprintf('%d: BG of %2.0f-->%2.0f at t=%2.0f (%2.0f), %2.1f hours since %d g carb intake at %2.0f (%2.0f), insulin=%2.2f', ind,M(ind,1), M(ind,2), t,time_series(t),dt(ind),F(ind),carb_time,time_series(carb_time),D(ind))) ;
        ind = ind + 1 ;
        carb_time = -1 ;
    end
end



