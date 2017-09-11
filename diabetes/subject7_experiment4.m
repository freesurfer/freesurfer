%

BG0 = [166.5 163.5  159.5 168 166.5 164 159.5 156.5 151 147 147 154 ...
       176 191.5 210.5 222.5 228.5 235.5 260.5 247.5 237.5 229.5 ...
       216.5 213 197];

min_per_hour = 60 ;
carb_delay = 95/min_per_hour;
protein = 22 ;
protein_to_carb_ratio = 0.1 ;
carbs = 50 + + protein * protein_to_carb_ratio ;
dt  = 10;
time = [0:dt:240];
time = time ./ min_per_hour ;
dt = dt ./ min_per_hour ;

if (dt == 10/60)
  BG = BG0;
else

  BG = zeros(size(time)) ;
  BG(1:2:end) = BG0 ;
  for i=2:2:length(BG)-1
 	BG(i) = (BG(i-1) + BG(i+1))/2;
  end
end

% subject specific stuff
basal1 = .6*dt ;
basal2 = .6*dt ;

basal_time = -10:dt:time(end) ;  % in hours
[foo, zero_ind] = min(abs(basal_time*60-0));  % time 0 in spreadsheet

% time0 = 8:29
% basals: 
% 0-2:        0.375
% 2-5:30:     0.35
% 5:30-7:30:  0.5
% 7:30-end    0.6

zero_time = 8.5 ; % 8:30 or so
basal = zeros(size(basal_time));
[foo, basal2_index] = min(abs(zero_time+basal_time-2));
[foo, basal3_index] = min(abs(zero_time+basal_time-5.5));
[foo, basal4_index] = min(abs(zero_time+basal_time-7.5));
basal(1:basal2_index-1) = .375 ;
basal(basal2_index:basal3_index-1) = .35 ;
basal(basal3_index:basal4_index-1) = .5 ;
basal(basal4_index:end) = .6 ;
basal = basal * dt ;  % divide it up into 10min intervals
% bolus times
% 9:19   0.45
% 9:55    0.05
% 10:00  5.4

bolus = 5.4;
[foo, bolus1_index] = min(abs(zero_time+time-9.33));
[foo, bolus2_index] = min(abs(zero_time+time-9.92));
[foo, bolus3_index] = min(abs(zero_time+time-10));
bolus_schedule = zeros(size(time)) ;
bolus_schedule(bolus1_index) = bolus_schedule(bolus1_index)+0.45 ;
bolus_schedule(bolus2_index) = bolus_schedule(bolus2_index)+0.05 ;
bolus_schedule(bolus3_index) = bolus_schedule(bolus3_index)+5.4 ;

[foo,time_start_index] = min(abs(basal_time-time(1)));
[foo,time_end_index] = min(abs(basal_time-time(end)));

pump_basal = basal ;

% given at 8:52 AM. (time 0)
% food eaten at 8:56


estimateIOBTable;

SI_best = 100 ; % insulin sensitivity, BG/insulin
b1_best = basal1 ;
b2_best = basal2 ;
C_best = 20 ; % carb ratio, carbs/insulin
SC_best = SI_best / C_best ;
k_carbs_best = .005 ;
bci_best = 1 ;
k_insulin_best = .025 ;
basal = pump_basal ;
basal(zero_ind:bci_best-1) = b1_best ;
basal(bci_best:end) = b2_best;
[bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin_best,.001,k_insulin_best);
basal_absorbed = babs(time_start_index:time_end_index) ;
[pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin_best,.001,k_insulin_best);
pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);
basal_difference = pump_basal_absorbed - basal_absorbed  ;
[active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin_best,.001, k_insulin_best);

carbs_metabolized = predict_carb_uptake(carbs, .0, k_carbs_best, time,carb_delay) ;


rms_best = compute_bg_error(BG(1),BG, bolus_absorbed, carbs_metabolized, time, SI_best, basal_difference, SC_best); 

M0 = BG(1) ;

disp(sprintf('**** new optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.3f, k_insulin=%2.4f, SC=%2.1f, SI=%2.0f, C=%2.1f *****', rms_best,b1_best/dt, b2_best/dt,k_carbs_best,k_insulin_best,SC_best,SI_best, C_best));

rms_best = sqrt((length(time) * 400^2)/length(time)) ;

k_carb_range = 0.0025:0.0025:0.04;
k_carbs_table = zeros(length(k_carb_range), length(time)) ;
ci = 1 ;
for k=k_carb_range
  k_carbs_table(ci,:) = predict_carb_uptake(carbs, .0, k, time,carb_delay) ;
  ci = ci+1 ;
end

start_ind = 1 ;
for k=k_carb_range
  if (k_carbs_table(start_ind,end) > .5*carbs)
    break ;
  end
  start_ind = start_ind + 1 ;
end
%k_carbs_table = k_carbs_table(start_ind:end,:);
%k_carb_range = k_carb_range(start_ind:end) ;

basal_search_range = 0.25;
b1_range=(basal1-basal_search_range*dt):.025*dt:(basal1+3*basal_search_range*dt);
b2_range=(basal2-basal_search_range*dt):.025*dt:(basal2+basal_search_range*dt);
bci_range = (zero_ind+(round(-1.0/dt)):round(0.5/dt):(zero_ind+(round(3.5/dt))));
% bci_range = basal_change_index ;
C_range = 5:1:15 ;
SI_range = 40:2:150 ;
SC_range = 5:.5:30;
Clen = length(C_range) ;
SIlen = length(SI_range) ;
rms_table = zeros(Clen, SIlen);

k_insulin_range=.005:0.0025:.04;
start_ind = 1 ;
for k_insulin=k_insulin_range
   [active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin,.001,k_insulin);
   if (bolus_absorbed(end) >= 0.25*bolus)
    break ;
  end
   start_ind = start_ind + 1 ;
end
% k_insulin_range = k_insulin_range(start_ind:end) 


pool = parpool(12) ;

for bci=bci_range
  for b1=b1_range
    disp(sprintf('searching BCI=%2.1f (%d), basal1 = %2.3f', (bci-zero_ind)*dt, bci,b1/dt)) ;
    for b2=b2_range
      if (b2>b1)
%	continue ;
      end
      
      basal = pump_basal ;
      basal(zero_ind:bci-1) = b1 ;
      basal(bci:end) = b2;
      
      for k_insulin=k_insulin_range
	[active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin,.001,k_insulin);
	[bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin,.001,k_insulin);
	basal_absorbed = babs(time_start_index:time_end_index) ;
	[pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin,.001,k_insulin);
	pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);
	basal_difference = pump_basal_absorbed - basal_absorbed  ;
	ci = 1 ;
	for k_carbs=k_carb_range
	  %       carbs_metabolized = predict_carb_uptake(carbs, .0, k_carbs, time,carb_delay) ;
	  carbs_metabolized = k_carbs_table(ci,:) ;
	  ci = ci + 1 ;
	  parfor CRi=1:Clen
	    C = C_range(CRi) ;
	    for SIi=1:SIlen
	      SI = SI_range(SIi) ;
	      SC = SI / C ;
	      if (SC >= SC_range(1) & SC < SC_range(end))
		rms = compute_bg_error(BG(1), BG, bolus_absorbed, carbs_metabolized, time, SI, basal_difference, SC); 
	      else 
		rms = sqrt((length(time) * 400^2)/length(time)) ;
              end
	      rms_table(CRi,SIi) = rms ;
	    end
	  end
	  
	  [test_min_rms,index] = min(rms_table(:)) ;
	  if (test_min_rms < rms_best)
	    rms_best = test_min_rms ; 
	    CRi = mod(index,size(rms_table,1)) ;
	    if (CRi == 0)
	      CRi = size(rms_table,1) ;
	    end

	    SIi = ceil(index/size(rms_table,1));
	    C = C_range(CRi) ;
	    SI = SI_range(SIi) ;
	    k_insulin_best = k_insulin ;
	    k_carbs_best = k_carbs ;
	    b1_best = b1 ;
	    b2_best = b2 ;
	    C_best = C ;
	    SI_best = SI ;
	    SC_best = SI_best/C_best ;
	    bci_best = bci ;
	    
	    [pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin_best,.001,k_insulin_best);
	    pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);
	    basal = pump_basal ;
	    basal(zero_ind:bci_best-1) = b1_best ;
	    basal(bci_best:end) = b2_best;
	    [bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin_best,.001, k_insulin_best);
	    basal_absorbed = babs(time_start_index:time_end_index) ;
	    basal_difference = pump_basal_absorbed - basal_absorbed  ;
	    
	    [active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin_best,.001,k_insulin_best);
	    disp(sprintf('**** new optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.3f, k_insulin=%2.4f, SC=%2.1f, SI=%2.1f, C=%2.2f, BCI=%2.1f *****', rms_best, b1_best/dt, b2_best/dt, k_carbs_best, k_insulin_best, SC_best,SI_best, C_best,(bci_best-zero_ind)*dt));
	    [rms, BG_p] = compute_bg_error(BG(1),BG, bolus_absorbed,  carbs_metabolized, time, SI_best, basal_difference, SC_best);
	    figure(1) ;
	    plot(time*60, BG, 'r') ;
	    title(sprintf('B1=%2.3f, B2=%2.3f, k=%2.3f, RMS=%2.1f', b1_best/dt, b2_best/dt, k_carbs_best, rms_best)) ;
	    hold on ;
	    plot(time*60, BG_p, 'g') ;
	    set(gca, 'fontsize', 16, 'fontweight', 'bold','xlim', [time(1)*60-10 time(end)*60+10]);
	    ch = get(gca, 'children') ;
	    set(ch, 'linewidth', 6) ;
	    hold off;
	    drawnow ;
	    
	    figure(2) ;
	    clf ;
	    [ax, h1, h2] = plotyy(60*time, bolus_absorbed, 60*time, carbs_metabolized) ;
	    set(gca, 'fontsize', 16, 'fontweight', 'bold','xlim', [time(1)*60-10 time(end)*60+10]);
	    legend('insulin', 'carbs')
	    xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold')
	    set(ax(1), 'fontsize', 16, 'fontweight', 'bold')
	    set(ax(2), 'fontsize', 16, 'fontweight', 'bold')
	    axes(ax(1)) ; ylabel('insulin absorbed (units)') ;
	    axes(ax(2)) ; ylabel('carbs metabolized (g)') ;
	    set(get(ax(1), 'children'), 'linewidth', 6) ;
	    set(get(ax(2), 'children'), 'linewidth', 6) ;
	    set(ax(1), 'ylim', [0 bolus]) ;
%	    drawnow ;
	  end
	end
      end
    end
    
    basal = pump_basal ;
    basal(zero_ind:bci_best-1) = b1_best ;
    basal(bci_best:end) = b2_best;
    [bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin_best,.001, k_insulin_best);
    basal_absorbed = babs(time_start_index:time_end_index) ;
    [pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin_best,.001,k_insulin_best);
    pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);
    basal_difference = pump_basal_absorbed - basal_absorbed  ;
    [active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin_best,.001,k_insulin_best);
    carbs_metabolized = predict_carb_uptake(carbs, .0, k_carbs_best, time,carb_delay) ;
    [rms, BG_p] = compute_bg_error(BG(1),BG, bolus_absorbed, carbs_metabolized, time, SI_best, basal_difference, SC_best);
    figure(1) ;
    clf ;
    plot(time*60, BG, 'r') ;
    set(gca, 'fontsize', 16, 'fontweight', 'bold','xlim', [time(1)*60-10 time(end)*60+10]);
    title(sprintf('B1=%2.3f, B2=%2.3f, k=%2.3f, RMS=%2.1f', b1_best/dt, b2_best/dt, k_carbs_best, rms_best)) ;
    hold on ;
    plot(time*60, BG_p, 'g') ;
    ch = get(gca, 'children') ;
    set(ch, 'linewidth', 6) ;
    
    legend('measured BG', 'model BG')
    xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold')
    ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold')
    set(gca, 'fontsize', 16, 'fontweight', 'bold','xlim', [time(1)*60-10 time(end)*60+10]);
		    
    
    hold off;
    figure(2) ;
    clf ;
    [ax, h1, h2] = plotyy(time*60, bolus_absorbed, time*60, carbs_metabolized) ;
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
  end
end


delete(pool) ;
disp(sprintf('**** global optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.3f, k_insulin=%2.4f, SC=%2.1f, SI=%2.1f, C=%2.2f, BCI=%2.1f hours, carb delay %d min *****', rms_best, b1_best/dt, b2_best/dt, k_carbs_best, k_insulin_best, SC_best,SI_best, C_best,dt*(bci_best-zero_ind),carb_delay*60));
