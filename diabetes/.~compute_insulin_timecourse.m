function insulin = compute_insulin_timecourse(bolus, total_time,K_fluid_to_blood, K_blood_to_fluid, K_blood_to_cells, dt)
% insulin = compute_insulin_timecourse(bolus, total_time,K_fluid_to_blood, K_blood_to_fluid, K_blood_to_cells, dt)


if (nargin < 6)
   dt = 1 ;
end


insulin = zeros(total_time/dt,1) ;


time = 1:dt:total_time ;

insulin_fluid = bolus ;
insulin_blood = 0 ;
insulin_used = 0 ;

% exponential model
for t=1:dt:total_time
    delta_fluid = -insulin_fluid*K_fluid_to_blood*dt ;
    delta_blood_to_fluid = -insulin_blood*K_blood_to_fluid*dt ;
    delta_blood_to_cells = -insulin_blood*K_blood_to_cells*dt ;
    insulin_used = insulin_used - delta_blood_to_cells ;
    insulin_blood = insulin_blood + delta_blood_to_cells + delta_blood_to_fluid - delta_fluid;
    insulin_fluid = insulin_fluid + delta_fluid - delta_blood_to_fluid ;
    insulin(t) = bolus - insulin_used ;
end



