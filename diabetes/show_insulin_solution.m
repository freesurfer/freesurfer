function  show_insulin_solution(parms)


carbs = parms.carbs ;
carb_delay = parms.carb_delay ;
BG = parms.BG ;
k_insulin_best = parms.k_insulin ;
k_carbs_best = parms.k_carbs ;
b1_best = parms.b1 ;
b2_best = parms.b2 ;
C_best = parms.C ;
SI_best = parms.SI ;
SC_best = SI_best/C_best ;
bci_best = parms.bci ;
basal_time = parms.basal_time;
time_start_index = parms.time_start_index ;
time_end_index = parms.time_end_index ;
time = parms.time ;
pump_basal = parms.pump_basal ;
b1_start_ind = parms.b1_start_ind ;
bolus_schedule = parms.bolus_schedule ;
[CR, CRi]=min(abs(parms.C_range-C_best));
dt = parms.dt ;
zero_ind = parms.zero_ind ;
bolus = parms.bolus ;
carbs_metabolized = predict_carb_uptake(carbs, .0, k_carbs_best, time,carb_delay) ;

[pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin_best,.001,k_insulin_best);
pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);

basal = pump_basal ;
basal(b1_start_ind:bci_best) = b1_best ;
basal(bci_best:end) = b2_best;
	    
[bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin_best,.001, k_insulin_best);
basal_absorbed = babs(time_start_index:time_end_index) ;
basal_difference = pump_basal_absorbed - basal_absorbed  ;
	    
[active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin_best,.001,k_insulin_best);
[rms, BG_p] = compute_bg_error(BG(1),BG, bolus_absorbed,  carbs_metabolized, time, SI_best, basal_difference, SC_best);
disp(sprintf('**** new optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.4f, k_insulin=%2.4f, SC=%2.1f, SI=%2.1f, C=%2.2f, BCI=%2.1f *****', rms, b1_best/dt, b2_best/dt, k_carbs_best, k_insulin_best, SC_best,SI_best, C_best,(bci_best-zero_ind)*dt));
figure(1) ;
plot(time*60, BG, 'r') ;
title(sprintf('B1=%2.3f, B2=%2.3f, k=%2.3f, RMS=%2.1f', b1_best/dt, b2_best/dt, k_carbs_best, rms)) ;
hold on ;
plot(time*60, BG_p, 'g') ;
set(gca, 'fontsize', 16, 'fontweight', 'bold','xlim', [time(1)*60-10 time(end)*60+10]);
ch = get(gca, 'children') ;
set(ch, 'linewidth', 6) ;
hold off;
	    
figure(2) ;
clf ;
[ax, h1, h2] = plotyy(60*time, bolus_absorbed, 60*time, carbs_metabolized) ;
set(ax, 'fontsize', 16, 'fontweight', 'bold','xlim', [time(1)*60-10 time(end)*60+10]);
legend('insulin', 'carbs')
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold')
set(ax(1), 'fontsize', 16, 'fontweight', 'bold')
set(ax(2), 'fontsize', 16, 'fontweight', 'bold')
axes(ax(1)) ; ylabel('insulin absorbed (units)') ;
axes(ax(2)) ; ylabel('carbs metabolized (g)') ;
set(get(ax(1), 'children'), 'linewidth', 6) ;
set(get(ax(2), 'children'), 'linewidth', 6) ;
set(ax(1), 'ylim', [0 bolus]) ;

drawnow ;

