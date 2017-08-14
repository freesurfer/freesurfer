BG = [143.5 147.5 155.5 166 216 222.5 240 248.5 264.5 294.5 309.5 327.5 332 330.5 326 325 326 326.5 324 310.5 312.5 306.5 296 293.5 280.5 273.5];

min_per_hour = 60 ;
dt  = 10;
time = [-10:dt:240];
time = time ./ min_per_hour ;
dt = dt ./ min_per_hour ;

ntp = length(time) ;
carbs = 50;
bolus = 1.65;

time_to_metabolize = -100/min_per_hour ;
interval=find(time>=time_to_metabolize) ;

% given at 9:01 AM.
% food eaten at 8:56


estimateIOBTable;


active_insulin(1) = 0 ;
active_insulin(2) = bolus ;

for i=2:ntp
	t = round((i-2)*dt*min_per_hour) ;
  active_insulin(i) = bolus*iob_table4(t+1)/100;
  A(i,1) = active_insulin(i-1) ;
  A(i,2) = active_insulin(i) ;

  insulin_absorbed(i) = bolus - active_insulin(i) ;
  M(i,1) = BG(i-1) ;
  M(i,2) = BG(i) ;

end
carbs_eaten = ones(size(active_insulin))*carbs;
doses = zeros(size(active_insulin));
doses(2) = bolus ;
carbs_eaten(2) = carbs ;

SI_best = 50 ; % insulin sensitivity, BG/insulin
e_best = -.2 ;
C_best = 10 ; % carb ratio, carbs/insulin
SC_best = SI_best / C_best ;
k_s2b_best = .005 ;
k_insulin_best = .025 ;
it =  compute_insulin_timecourse(1, 60*(time(end)-time(1)), k_insulin_best, 0 , k_insulin_best) ;
for i=2:ntp
   active_insulin(i) = bolus*it(round((i-1)*dt*60)) ;
   insulin_absorbed(i) = bolus - active_insulin(i) ;
end

carbs_metabolized = predict_carb_uptake(carbs, .0, k_s2b_best, time) ;


rms_best = compute_bg_error(BG(1),BG(interval), insulin_absorbed(interval), doses(interval), carbs_metabolized(interval), time(interval), SI_best, e_best, SC_best); 

M0 = BG(1) ;

disp(sprintf('**** new optimum %2.3f found at de=%2.2f, k=%2.3f, SC=%2.1f, SI=%2.0f, C=%2.1f *****', rms_best,e_best,k_s2b_best,SC_best,SI_best, C_best));

for k_stomach_to_blood=.006:.001:.1
  carbs_metabolized = predict_carb_uptake(carbs, .0, k_stomach_to_blood, time) ;
  disp(sprintf('searching k = %2.3f', k_stomach_to_blood)) ;
  for k_insulin=.005:0.001:.04
     it =  compute_insulin_timecourse(1, 60*(time(end)-time(1)), k_insulin, 0, k_insulin) ;
     for i=2:ntp
	t = round((i-2)*dt*min_per_hour) ;
        active_insulin(i) = bolus*it(round((i-1)*dt*60)) ;
        insulin_absorbed(i) = bolus - active_insulin(i) ;
      end

  for de=-0.5:.025:0.5
    for C=5:.1:30
      for SI=10:1:200
	     SC = SI / C ;
             rms = compute_bg_error(BG(1), BG(interval), insulin_absorbed(interval), doses(interval), carbs_metabolized(interval), time(interval), SI, de, SC); 
           if (rms < rms_best)
	     rms_best = rms ; 
             k_insulin_best = k_insulin ;
             it =  compute_insulin_timecourse(1, 60*(time(end)-time(1)), k_insulin_best, 0 , k_insulin_best) ;
             for i=2:ntp
               active_insulin(i) = bolus*it(round((i-1)*dt*60)) ;
               insulin_absorbed(i) = bolus - active_insulin(i) ;
             end
	     e_best = de ;
	     C_best = SI/SC ;
	     SI_best = SI ;
	     SC_best = SC ;
             k_s2b_best = k_stomach_to_blood ;
             disp(sprintf('**** new optimum %2.3f found at de=%2.2f, k_carbs=%2.3f, k_insulin=%2.4f, SC=%2.1f, SI=%2.0f, C=%2.1f *****', rms,e_best,k_s2b_best,k_insulin_best,SC,SI, C_best));
             [rms, BG_p] = compute_bg_error(BG(1),BG(interval), insulin_absorbed(interval), doses(interval), carbs_metabolized(interval), time(interval), SI_best, e_best, SC_best);
             figure(1) ;
             plot(time*60, BG, 'r') ;
             title(sprintf('k=%2.3f, RMS=%2.1f', k_s2b_best, rms_best)) ;
             hold on ;
             plot(time(interval)*60, BG_p, 'g') ;
             hold off;
             drawnow ;
             figure(2) ;
             plotyy(time, active_insulin, time, carbs_metabolized) ;
             drawnow ;
           end
      end
   end
  end
  it =  compute_insulin_timecourse(1, 60*(time(end)-time(1)), k_insulin_best, 0 , k_insulin_best) ;
  for i=2:ntp
     active_insulin(i) = bolus*it(round((i-1)*dt*60)) ;
     insulin_absorbed(i) = bolus - active_insulin(i) ;
  end
  carbs_metabolized = predict_carb_uptake(carbs, .0, k_s2b_best, time) ;
  [rms, BG_p] = compute_bg_error(BG(1),BG(interval), insulin_absorbed(interval), doses(interval), carbs_metabolized(interval), time(interval), SI_best, e_best, SC_best);
  figure(1) ;
  plot(time*60, BG, 'r') ;
  title(sprintf('k=%2.3f, RMS=%2.1f', k_s2b_best, rms_best)) ;
  hold on ;
  plot(time(interval)*60, BG_p, 'g') ;

  legend('actual BG', 'model BG')
  xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold')
  ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold')
  set(gca, 'fontsize', 16, 'fontweight', 'bold')

  hold off;
  figure(2) ;
  [ax, h1, h2] = plotyy(time, active_insulin, time, carbs_metabolized) ;
  legend('insulin', 'carbs')
  xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold')
  set(ax(1), 'fontsize', 16, 'fontweight', 'bold')
  set(ax(2), 'fontsize', 16, 'fontweight', 'bold')
  set(ax(1), 'ylabel', 'active insulin (units)') ;
  set(ax(2), 'ylabel', 'carbs metabolized (grams)') ;
  drawnow ;
end
end


