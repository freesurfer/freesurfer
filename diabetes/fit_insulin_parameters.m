b1_range=(basal1_min*dt):.025*dt:(basal1_max*dt);
b2_range=(basal2_min*dt):.025*dt:(basal2_max*dt);
ind = find(b1_range>=0) ;
b1_range = b1_range(ind) ;
ind = find(b2_range>=0) ;
b2_range = b2_range(ind) ;


if (dt == 10/60)
  BG = BG0;
else

  BG = zeros(size(time)) ;
  BG(1:2:end) = BG0 ;
  for i=2:2:length(BG)-1
 	BG(i) = (BG(i-1) + BG(i+1))/2;
  end
end


pump_basal = basal ;

estimateIOBTable;

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

disp(sprintf('**** new optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.4f, k_insulin=%2.4f, SC=%2.1f, SI=%2.0f, C=%2.1f *****', rms_best,b1_best/dt, b2_best/dt,k_carbs_best,k_insulin_best,SC_best,SI_best, C_best));

rms_best = sqrt((length(time) * 400^2)/length(time)) ;

k_carbs_range = 0.015:0.001:0.04;
k_carbs_range = 0.0025:0.0001:0.02;
k_carbs_range = 0.0025:0.001:0.02;
k_carbs_table = zeros(length(k_carbs_range), length(time)) ;
ci = 1 ;
for k=k_carbs_range
  k_carbs_table(ci,:) = predict_carb_uptake(carbs, .0, k, time,carb_delay) ;
  ci = ci+1 ;
end

start_ind = 1 ;
for k=k_carbs_range
  if (k_carbs_table(start_ind,end) > pct_carbs_absorbed_at_end*carbs)
    break ;
  end
  start_ind = start_ind + 1 ;
end
k_carbs_table = k_carbs_table(start_ind:end,:);
k_carbs_range = k_carbs_range(start_ind:end) ;

Clen = length(C_range) ;
SIlen = length(SI_range) ;
rms_table = zeros(Clen, SIlen);

k_insulin_range=.01:0.0025:.15;
start_ind = 1 ;
for k_insulin=k_insulin_range
   [active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin,.001,k_insulin);
   if (bolus_absorbed(end) >= pct_insulin_absorbed_at_end*bolus)
     break ;
  end
   start_ind = start_ind + 1 ;
end
k_insulin_range = k_insulin_range(start_ind:end) ;


% these are static (they don't change with the optimum)
parms.BG = BG ;	    
parms.zero_time = zero_time ;
parms.carbs = carbs ;
parms.carb_delay = carb_delay ;
parms.dt = dt ;
parms.time_start_index = time_start_index ;
parms.time_end_index = time_end_index ;
parms.pump_basal = pump_basal ;	    
parms.b1_start_ind = b1_start_ind ;
parms.bolus_schedule = bolus_schedule ;
parms.basal_time = basal_time;
parms.time = time ;
parms.zero_ind = zero_ind ;
parms.bolus = bolus ;
parms.C_range = C_range ;
parms.k_carbs_table = k_carbs_table ;
	    

disp(sprintf('searching b1 %2.3f->%2.3f, b2 %2.3f->%2.3f, bci %2.1f->%2.1f',...
	     b1_range(1)/dt, b1_range(end)/dt, b2_range(1)/dt, ...
	     b2_range(end)/dt, (bci_range(1)-zero_ind)*dt, ...
	     (bci_range(end)-zero_ind)*dt));
disp(sprintf('\tk_insulin %2.4f->%2.4f, k_carbs %2.4f->%2.4f, C %2.0f->%2.0f, SI %2.0f->%2.0f, SC %2.1f->%2.1f\n',...
	     k_insulin_range(1), k_insulin_range(end),...
	     k_carbs_range(1), k_carbs_range(end),...
	     C_range(1), C_range(end),...
	     SI_range(1), SI_range(end),...
	     SC_range(1), SC_range(end)));
	     
pool = parpool([10 20], 'IdleTimeout', 120) ;

for bci=bci_range
  for b1=b1_range
    disp(sprintf('searching BCI=%2.1f (%d), basal1 = %2.3f', (bci-zero_ind)*dt, bci,b1/dt)) ;
    for b2=b2_range
      
      basal = pump_basal ;
      basal(b1_start_ind:bci) = b1 ;
      basal(bci:end) = b2;
      
      for k_insulin=k_insulin_range
	[active_insulin, bolus_absorbed] = compute_total_insulin_timecourse(bolus_schedule, 60*time, k_insulin,.001,k_insulin);
	[bact, babs] = compute_total_insulin_timecourse(basal, 60*basal_time, k_insulin,.001,k_insulin);
	basal_absorbed = babs(time_start_index:time_end_index) ;
	[pump_basal_active_long, pump_basal_absorbed_long] = compute_total_insulin_timecourse(pump_basal, 60*basal_time, k_insulin,.001,k_insulin);
	pump_basal_absorbed = pump_basal_absorbed_long(time_start_index:time_end_index);
	basal_difference = pump_basal_absorbed - basal_absorbed  ;
	ci = 1 ;
	for k_carbs=k_carbs_range
	  carbs_metabolized = k_carbs_table(ci,:) ;
	  ci = ci + 1 ;
	  parfor CRi=1:Clen
	    C = C_range(CRi) ;
	    for SIi=1:SIlen
	      SI = SI_range(SIi) ;
	      SC = SI / C ;
	      if (SC >= SC_range(1) & SC <= SC_range(length(SC_range)))
		rms = compute_bg_error(BG(1), BG, bolus_absorbed, ...
				       carbs_metabolized, time, SI, basal_difference, SC); 
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

	    parms.carbs_metabolized = carbs_metabolized ;
	    parms.k_insulin_best = k_insulin ;
	    parms.b1 = b1 ;
	    parms.b2 = b2 ;
	    parms.C = C ;
	    parms.SI = SI ;
	    parms.SC = SI_best/C_best ;
	    parms.bci = bci ;
	    parms.basal = basal ;
	    parms.k_insulin = k_insulin_best ;
	    parms.k_carbs = k_carbs_best ;

	    best_parms = parms ;
	    show_insulin_solution(parms) ;
	  
	  end
	end
      end
    end

    show_insulin_solution(best_parms) ;
  end
end


show_insulin_solution(best_parms) ;
delete(pool) ;
disp(sprintf('**** global optimum %2.3f found at B1=%2.3f, B2=%2.3f k_carbs=%2.4f, k_insulin=%2.4f, SC=%2.1f, SI=%2.1f, C=%2.2f, BCI=%2.1f hours, carb delay %d min *****', rms_best, b1_best/dt, b2_best/dt, k_carbs_best, k_insulin_best, SC_best,SI_best, C_best,dt*(bci_best-zero_ind),carb_delay*60));
