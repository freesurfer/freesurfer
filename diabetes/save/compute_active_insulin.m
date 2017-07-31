function active_insulin_time_series = extract_active_insulin(S, sind,time_series, bolus_vec, basal_vec,iob_table)
% active_insulin_time_series = extract_active_insulin(S, sind,time_series, bolus_vec, basal_vec,iob_table)


ntps = length(time_series) ;

active_insulin_time_series = zeros(ntps,1);


% track the start time of the bolus and the amount for caculating
% current active
max_boluses = 1000;
active_boluses = zeros(max_boluses, 2) ;
num_active_boluses = 0 ;
dt = 5 ;   % model basal as small boluses given in this time interval

basal_vec = basal_vec * dt / 60 ;  % turn it into amount to give every dt minutes

last_basal_time = -(60*3) ; % assume a 3 hour active life for the basal
start_index = 4 ;
for t=start_index:ntps
    % integrate the basal given in this time period
     if (t > 1)
        for t1=last_basal_time+dt:dt:(time_series(t)-dt)
            [next_index, active_boluses, num_active_boluses] = add_bolus(active_boluses, basal_vec(t-1), t1, num_active_boluses,iob_table) ;
            last_basal_time = t1 ;
        end
    end
    if (time_series(t)+.01 >= last_basal_time+dt)
        [next_index, active_boluses, num_active_boluses] = add_bolus(active_boluses, basal_vec(t), time_series(t), num_active_boluses,iob_table) ;
        last_basal_time = time_series(t) ;
    end
    if (bolus_vec(t) > 0)
        [next_index, active_boluses, num_active_boluses] = add_bolus(active_boluses, bolus_vec(t), time_series(t), num_active_boluses,iob_table) ;
    end
    
    % now add up all the active boluses  to compute the amount of active insulin
    active_insulin = 0 ;
    for bolus_ind = 1:num_active_boluses
        deltat = round(time_series(t) - active_boluses(bolus_ind,2) + 1) ;
        if (deltat > length(iob_table))  % expired bolus
            continue ;
        end
        bolus = active_boluses(bolus_ind,1) ;
        if (bolus > 0)  % otherwise deltat might be out of bounds
            active_insulin = active_insulin + bolus * iob_table(deltat)/100 ;
        end
    end
    active_insulin_time_series(t,1) = active_insulin ;
end

for t=1:start_index-1
    active_insulin_time_series(t) = active_insulin_time_series(start_index) ;
end


