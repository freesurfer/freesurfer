dbcleind = 1 ;
carb_time = -1 ;
for t=1:ntps
    if (carb_vec(t,1) > 0)
        carbs = carb_vec(t,1) ;
        carb_time = t ;
        M(ind,1) = sensor_vec(t) ;
    end
    if (carb_vec(t,2) > 120 & carb_time > 0)
        % compute total insulin absobed in that time
        F(ind) = carbs ;
        M(ind,2) = sensor_vec(t) ;
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
                disp(sprintf('adding insulin vec at t=%d: %d = %d', t1,bolus_vec(t1,1), insulin));
            end
        end
        
        D(ind) = insulin ;
        dt(ind) = time_series(t) - time_series(carb_time) ;
        disp(sprintf('BG of %f-->%f at t=%d (%d), %d minutes since %d carb intake, insulin=%d', M(ind,1), M(ind,2), t,time_series(t),dt(ind),F(ind),D(ind))) ;
        ind = ind + 1 ;
        carb_time = -1 ;
    end
end



