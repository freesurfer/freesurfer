function [insulin_active, insulin_absorbed] = compute_insulin_absorbed(insulin_given, time,k_fluid_to_blood, k_blood_to_fluid, k_blood_to_cells, dt)
  %  [insulin_active, insulin_absorbed] = compute_insulin_absorbed(insulin_given, time,k_fluid_to_blood, k_blood_to_fluid, k_blood_to_cells, dt)


ntp = length(time) ;
if (nargin < 6)
   dt = 1 ;
end
it =  compute_total_insulin_timecourse(insulin_given, k_fluid_to_blood, k_blood_to_fluid, k_blood_to_cells, dt) ;
for i=2:ntp
   active_insulin(i) = bolus*it(round((i-1)*dt*60)) ;
   insulin_absorbed(i) = bolus - active_insulin(i) ;
end

