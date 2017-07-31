function [insulin_active, insulin_absorbed] = compute_total_insulin_timecourse(insulin_schedule, time,K_fluid_to_blood, K_blood_to_fluid, K_blood_to_cells, dt)
  % insulin = compute_total_insulin_timecourse(insulin_schedule,time,K_fluid_to_blood, K_blood_to_fluid, K_blood_to_cells, dt)
% time is in minutes

if (nargin < 6)
   dt = 1 ;
end

insulin_active = zeros(size(insulin_schedule));
insulin_absorbed = zeros(size(insulin_schedule));
total_time = time(end)-time(1);

insulin_fluid = 0 ;
insulin_blood = 0 ;
insulin_used = 0 ;

% exponential model
i = 1 ;
for t=time(1):dt:time(end)
    if (t+dt > time(i))
      insulin_fluid = insulin_fluid +  insulin_schedule(i) ;
      insulin_active(i) = insulin_blood ;
      insulin_absorbed(i) = insulin_used ;
      i = i + 1 ;
      if (i > length(time))
	break ;
      end
    end
    delta_fluid = -insulin_fluid*K_fluid_to_blood*dt ;
    delta_blood_to_fluid = -insulin_blood*K_blood_to_fluid*dt ;
    delta_blood_to_cells = -insulin_blood*K_blood_to_cells*dt ;
    insulin_used = insulin_used - delta_blood_to_cells ;
    insulin_blood = insulin_blood + delta_blood_to_cells + delta_blood_to_fluid - delta_fluid;
    insulin_fluid = insulin_fluid + delta_fluid - delta_blood_to_fluid ;
end



