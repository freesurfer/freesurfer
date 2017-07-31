function [next_index, active_boluses, max_boluses] = add_bolus(active_boluses,bolus_amount,bolus_time,max_boluses,iob_table)
% [next_index, active_boluses, max_boluses] = add_bolus(active_boluses,bolus-amount,bolus_time,max_boluses,iob_table)

next_index = -1 ;
if (bolus_amount*60/5 > 1)
    disp('big bolus');
end

for bolus_index=1:max_boluses
    deltat = round(bolus_time-active_boluses(bolus_index,2)+1) ;
    if (deltat > length(iob_table)) 
       	pct_active = 0 ;
    else
	pct_active = iob_table(deltat) ;
    end
    if (pct_active == 0)
       next_index = bolus_index ;
       break ;
    end
end

if (next_index < 0)
   max_boluses = max_boluses + 1 ;
   next_index = max_boluses  ;
end
if (bolus_time >0 & bolus_time<11620)
   disp(sprintf('bolus added at index %2.2d: amount = %f, time = %3.3d', next_index, bolus_amount*60/5,bolus_time));
end
active_boluses(next_index,1) = bolus_amount ;
active_boluses(next_index,2) = bolus_time ;

