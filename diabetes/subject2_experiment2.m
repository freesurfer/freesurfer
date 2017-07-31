% **** global optimum 6.750 found at B1=0.625, B2=0.625
% k_carbs=0.006, k_insulin=0.0210, SC=4.8, SI=19.0, C=4.00, BCI=0.0
% hours, carb delay 0 min ***** 

% best rms for delay=5: **** global optimum 6.144 found at
% B1=0.675, B2=0.675 k_carbs=0.008, k_insulin=0.0220, SC=6.2,
% SI=31.0, C=5.00, BCI=0.0 hours, carb delay 5 min *****

% *** global optimum 5.756 found at B1=0.700, B2=0.700
% k_carbs=0.023, k_insulin=0.0340, SC=3.7, SI=26.0, C=7.00, BCI=0.0
% hours, carb delay 10 min ***** 

% **** global optimum 5.653 found at B1=0.700, B2=0.700
% k_carbs=0.050, k_insulin=0.0390, SC=2.2, SI=18.0, C=8.00, BCI=0.0
% hours, carb delay 15 min ***** 

% time 0 is 9:38, earlier bolus of 1.2 given at 7:10
first_bolus_time = -(2*60 + (38-10))/60;
first_bolus = 1.2 ;

BG0_1 = [116.5 111.5 110 118 125 128 137.5 128 115 97.15 87.05 94 94 ...
       90.8 90.6 84.3 78.7 81.15 72.25 68.2 72 67.9 56.9 52.45 50.6 ...
       47.9];

BG0 = [66 57.5 59.45 64.7 74.85 89.75 85.85 90.65 87.45 100.5 104.5 ...
       107.5 116.5 130.5 137 127 122.5 118.5 108.5 105.5 102 92.55 ...
       87.9 84.3 81.05 76.85];

min_per_hour = 60 ;
carb_delay = 15/min_per_hour
carbs = 50;
protein = 22 ;
protein_to_carb_ratio = 0.1 ;
carbs = carbs + protein * protein_to_carb_ratio ;
dt  = 10;
time = [-10:dt:240];
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
bolus = 7.7 ;
basal1 = .7*dt ;
basal2 = .7*dt ;


pct_parm_change = 0.5 ;

basal_time = -10:dt:time(end) ;
basal = ones(size(basal_time))*basal1 ;
[foo, basal_change_index] = min(abs(basal_time*60-10));
[foo, bolus_index] = min(abs(basal_time*60-0));  % bolus at time 0
[foo, first_bolus_index] = min(abs(basal_time-first_bolus_time));  % bolus at 7:10
basal(first_bolus_index) = basal(first_bolus_index) + first_bolus ;
basal(basal_change_index:end) = basal(basal_change_index:end) *basal2/basal1;
total_insulin = basal ;
total_insulin(bolus_index) = total_insulin(bolus_index)+bolus ;
bolus_schedule = zeros(size(time)) ;
[foo, bolus_index] = min(abs(time*60-0));
bolus_schedule(bolus_index) = bolus ;
[foo,time_start_index] = min(abs(basal_time-time(1)));
[foo,time_end_index] = min(abs(basal_time-time(end)));

pump_basal = basal ;

% given at 8:52 AM. (time 0)
% food eaten at 8:56


estimateIOBTable;

SI_best = 50 ; % insulin sensitivity, BG/insulin
b1_best = basal1 ;
b2_best = basal2 ;
C_best = 10 ; % carb ratio, carbs/insulin
SC_best = SI_best / C_best ;
k_carbs_best = .005 ;
bci_best = basal_change_index ;
k_insulin_best = .025 ;
basal = ones(size(basal_time))*b1_best ;
basal(basal_change_index:end) = basal(basal_change_index:end) *b2_best/b1_best;
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

k_carb_range = 0.005:0.001:0.05;
k_carbs_table = zeros(length(k_carb_range), length(time)) ;

ci = 1 ;
for k=k_carb_range
  k_carbs_table(ci,:) = predict_carb_uptake(carbs, .0, k, time,carb_delay) ;
  ci = ci+1 ;
end

b1_range = .6*dt:.025*dt:.8*dt ;
b2_range = .5*dt:.025*dt:.85*dt ;
bci_range = basal_change_index;
bci_range = (basal_change_index-(round(2.5/dt)):round(0.5/dt):(basal_change_index+(round(1.5/dt))));
k_insulin_range=.005:0.001:.04;
C_range = 2:1:25 ;
SI_range = 40:1:90 ;
SC_range = 1.5:.5:15;
Clen = length(C_range) ;
SIlen = length(SI_range) ;
rms_table = zeros(Clen, SIlen);

pool = parpool(10) ;

for bci=bci_range
  for b1=b1_range
    disp(sprintf('searching BCI=%2.1f (%d), basal1 = %2.3f', dt*(bci-basal_change_index), bci,b1/dt)) ;
    for b2=b2_range
      basal = ones(size(basal_time))*b1 ;
      %    basal(basal_change_index:end) = basal(basal_change_index:end) *b2/b1;
      basal(bci:end) = basal(bci:end) *b2/b1;
      
      for k_insulin=k_insulin_range
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
	for k_carbs=k_carb_range
	  %       carbs_metabolized = predict_carb_uptake(carbs, .0, k_carbs, time,carb_delay) ;
	  carbs_metabolized = k_carbs_table(ci,:) ;
	  ci = ci + 1 ;
	  if (carbs_metabolized(end) < .9*carbs)
	    continue
	  end
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
	    basal = ones(size(basal_time))*b1_best ;
	    basal(bci_best:end) = basal(bci_best:end) *b2_best/b1_best;
	    [bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin_best,.001, k_insulin_best);
	    basal_absorbed = babs(time_start_index:time_end_index) ;
	    basal_difference = pump_basal_absorbed - basal_absorbed  ;
	    
	    [active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin_best,.001,k_insulin_best);
	    disp(sprintf('**** new optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.3f, k_insulin=%2.4f, SC=%2.1f, SI=%2.1f, C=%2.2f, BCI=%2.1f *****', rms_best, b1_best/dt, b2_best/dt, k_carbs_best, k_insulin_best, SC_best,SI_best, C_best,dt*(bci_best-basal_change_index)));
	    [rms, BG_p] = compute_bg_error(BG(1),BG, bolus_absorbed,  carbs_metabolized, time, SI_best, basal_difference, SC_best);
	    figure(1) ;
	    plot(time*60, BG, 'r') ;
	    title(sprintf('B1=%2.3f, B2=%2.3f, k=%2.3f, RMS=%2.1f', b1_best/dt, b2_best/dt, k_carbs_best, rms_best)) ;
	    hold on ;
	    plot(time*60, BG_p, 'g') ;
	    ch = get(gca, 'children') ;
	    set(ch, 'linewidth', 6) ;
	    hold off;
	    drawnow ;
	    
	    figure(2) ;
	    clf ;
	    [ax, h1, h2] = plotyy(60*time, bolus_absorbed, 60*time, carbs_metabolized) ;
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
    
    basal = ones(size(basal_time))*b1_best ;
    basal(bci_best:end) = basal(bci_best:end) *b2_best/b1_best;
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
    title(sprintf('B1=%2.3f, B2=%2.3f, k=%2.3f, RMS=%2.1f', b1_best/dt, b2_best/dt, k_carbs_best, rms_best)) ;
    hold on ;
    plot(time*60, BG_p, 'g') ;
    ch = get(gca, 'children') ;
    set(ch, 'linewidth', 6) ;
    
    legend('measured BG', 'model BG')
    xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold')
    ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold')
    set(gca, 'fontsize', 16, 'fontweight', 'bold')
    
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
disp(sprintf('**** global optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.3f, k_insulin=%2.4f, SC=%2.1f, SI=%2.1f, C=%2.2f, BCI=%2.1f hours, carb delay %d min *****', rms_best, b1_best/dt, b2_best/dt, k_carbs_best, k_insulin_best, SC_best,SI_best, C_best,dt*(bci_best-basal_change_index),carb_delay*60));
