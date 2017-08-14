BG = [143.5 147.5 155.5 166 216 222.5 240 248.5 264.5 294.5 309.5 327.5 332 330.5 326 325 326 326.5 324 310.5 312.5 306.5 296 293.5 280.5 273.5];

carbs = 50;
bolus = 1.65;
min_per_hour = 60 ;
dt  = 10;
time = [-10:dt:240];
time = time ./ min_per_hour ;
dt = dt ./ min_per_hour ;
basal1 = .95*dt ;
basal2 = .75*dt ;



basal_time = -10:dt:time(end) ;
basal = ones(size(basal_time))*basal1 ;
[foo, basal_change_index] = min(abs(basal_time*60-30));
[foo, bolus_index] = min(abs(basal_time*60));
basal(basal_change_index:end) = basal(basal_change_index:end) *basal2/basal1;
total_insulin = basal ;
total_insulin(bolus_index) = total_insulin(bolus_index)+bolus ;
bolus_schedule = zeros(size(time)) ;
bolus_schedule(2) = bolus ;
[foo,time_start_index] = min(abs(basal_time-time(1)));
[foo,time_end_index] = min(abs(basal_time-time(end)));

pump_basal = basal ;

time_to_metabolize = -100/min_per_hour ;

% given at 9:01 AM.
% food eaten at 8:56


estimateIOBTable;

SI_best = 50 ; % insulin sensitivity, BG/insulin
b1_best = basal1 ;
b2_best = basal2 ;
C_best = 10 ; % carb ratio, carbs/insulin
SC_best = SI_best / C_best ;
k_s2b_best = .005 ;
k_insulin_best = .025 ;
basal = ones(size(basal_time))*b1_best ;
basal(basal_change_index:end) = basal(basal_change_index:end) *b2_best/b1_best;
[bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin_best,.001,k_insulin_best);
basal_absorbed = babs(time_start_index:time_end_index) ;
[pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin_best,.001,k_insulin_best);
pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);
basal_difference = pump_basal_absorbed - basal_absorbed  ;
[active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin_best,.001, k_insulin_best);

carbs_metabolized = predict_carb_uptake(carbs, .0, k_s2b_best, time) ;


rms_best = compute_bg_error(BG(1),BG, bolus_absorbed, carbs_metabolized, time, SI_best, basal_difference, SC_best); 

M0 = BG(1) ;

disp(sprintf('**** new optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.3f, k_insulin=%2.4f, SC=%2.1f, SI=%2.0f, C=%2.1f *****', rms_best,b1_best/dt, b2_best/dt,k_s2b_best,k_insulin_best,SC_best,SI_best, C_best));

rms_best = sqrt((length(time) * 400^2)/length(time)) ;

k_s2b_min = .005 ;
k_s2b_max = .05 ;
k_s2b_step = 0.001 ;

k_s2b_num = 1+round((k_s2b_max-k_s2b_min)/k_s2b_step) ;
k_s2b = k_s2b_min:k_s2b_step:k_s2b_max ;
k_s2b_table = zeros(length(k_s2b), length(time)) ;

ci = 1 ;
for k=k_s2b
  k_s2b_table(ci,:) = predict_carb_uptake(carbs, .0, k, time) ;
  ci = ci+1 ;
end


C_range = 5:.5:20 ;
SI_range = 10:1:100 ;
Clen = length(C_range) ;
SIlen = length(SI_range) ;
rms_table = zeros(Clen, SIlen);

pool = parpool(20) ;

for b1=.65*dt:.025*dt:1.1*dt
  for b2=.5*dt:.025*dt:1*dt
    disp(sprintf('searching basal1 = %2.3f, basal2 = %2.3f', b1/dt, b2/dt)) ;
    basal = ones(size(basal_time))*b1 ;
    basal(basal_change_index:end) = basal(basal_change_index:end) *b2/b1;

    for k_insulin=.01:0.001:.03
      [active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin,.001,k_insulin);
      if (active_insulin(end) > 0.1*bolus)
        continue
      end
      [bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin,.001,k_insulin);
      basal_absorbed = babs(time_start_index:time_end_index) ;
      [pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin,.001,k_insulin);
      pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);
      basal_difference = pump_basal_absorbed - basal_absorbed  ;
      ci = 1 ;
      for k_stomach_to_blood=k_s2b
%       carbs_metabolized = predict_carb_uptake(carbs, .0, k_stomach_to_blood, time) ;
        carbs_metabolized = k_s2b_table(ci,:) ;
        ci = ci + 1 ;
        if (carbs_metabolized(end) < .9*carbs)
          continue
        end
        parfor CRi=1:Clen
          C = C_range(CRi) ;
          for SIi=1:SIlen
	    SI = SI_range(SIi) ;
	    SC = SI / C ;
            rms = compute_bg_error(BG(1), BG, bolus_absorbed, carbs_metabolized, time, SI, basal_difference, SC); 
            rms_table(CRi,SIi) = rms ;
          end
        end

	[test_min_rms,index] = min(rms_table(:)) ;
        if (test_min_rms < rms_best)
	     rms_best = test_min_rms ; 
             CRi = mod(index,size(rms_table,1)) ;
             SIi = ceil(index/size(rms_table,1));
             C = C_range(CRi) ;
             SI = SI_range(SIi) ;
             k_insulin_best = k_insulin ;
             k_s2b_best = k_stomach_to_blood ;
             b1_best = b1 ;
             b2_best = b2 ;
             C_best = C ;
	     SI_best = SI ;
	     SC_best = SI_best/C_best ;

             [pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin_best,.001,k_insulin_best);
             pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);
             basal = ones(size(basal_time))*b1_best ;
             basal(basal_change_index:end) = basal(basal_change_index:end) *b2_best/b1_best;
             [bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin_best,.001, k_insulin_best);
             basal_absorbed = babs(time_start_index:time_end_index) ;
             basal_difference = pump_basal_absorbed - basal_absorbed  ;

             [active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin_best,.001,k_insulin_best);
             disp(sprintf('**** new optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.3f, k_insulin=%2.4f, SC=%2.1f, SI=%2.1f, C=%2.2f *****', rms_best, b1_best/dt, b2_best/dt, k_s2b_best, k_insulin_best, SC_best,SI_best, C_best));
            [rms, BG_p] = compute_bg_error(BG(1),BG, bolus_absorbed,  carbs_metabolized, time, SI_best, basal_difference, SC_best);
             figure(1) ;
             plot(time*60, BG, 'r') ;
             ch = get(gca, 'children') ;
             title(sprintf('k=%2.3f, RMS=%2.1f', k_s2b_best, rms_best)) ;
             hold on ;
             plot(time*60, BG_p, 'g') ;
             hold off;
             drawnow ;
             figure(2) ;
             plotyy(60*time, active_insulin, time, carbs_metabolized) ;
             drawnow ;
           end
        end
     end
  end
  basal = ones(size(basal_time))*b1_best ;
  basal(basal_change_index:end) = basal(basal_change_index:end) *b2_best/b1_best;
  [bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin_best,.001, k_insulin_best);
  basal_absorbed = babs(time_start_index:time_end_index) ;
  [pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin_best,.001,k_insulin_best);
  pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);
  basal_difference = pump_basal_absorbed - basal_absorbed  ;
  [active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin_best,.001,k_insulin_best);
  carbs_metabolized = predict_carb_uptake(carbs, .0, k_s2b_best, time) ;
  [rms, BG_p] = compute_bg_error(BG(1),BG, bolus_absorbed, carbs_metabolized, time, SI_best, basal_difference, SC_best);
  figure(1) ;
  clf ;
  plot(time*60, BG, 'r') ;
  title(sprintf('B1=%2.3f, B2=%2.3f, k=%2.3f, RMS=%2.1f', b1_best/dt, b2_best/dt, k_s2b_best, rms_best)) ;
  hold on ;
  plot(time*60, BG_p, 'g') ;

  legend('actual BG', 'model BG')
  xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold')
  ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold')
  set(gca, 'fontsize', 16, 'fontweight', 'bold')

  hold off;
  figure(2) ;
  clf ;
  [ax, h1, h2] = plotyy(time*60, active_insulin, time*60, carbs_metabolized) ;
  legend('insulin', 'carbs')
  xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold')
  set(ax(1), 'fontsize', 16, 'fontweight', 'bold')
  set(ax(2), 'fontsize', 16, 'fontweight', 'bold')
  axes(ax(1)) ; ylabel('active insulin (units)') ;
  axes(ax(2)) ; ylabel('carbs metabolized') ;
  drawnow ;
end

delete(pool) ;
