parms.time_series=1:5*60 ;

parms.G_target = 120 ;   % desired blood glucose (Gp) level
parms.carb_sensitivity = 7 ;
parms.insulin_sensitivity = 180 ;  % glucose/unit insulin


parms.insulin_delay = 15 ;


parms.k1 = .021 ;   % rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
parms.k2 = .001 ;   % rate at which insulin moves from plasma to interstitial fluid (units/min)
parms.k3 = .021 ;   % rate at which insulin moves from fluid to plasma (units/min)
parms.k4 = .05 ;    % rate at which carbs are metabolized from stomach to blood (grams/min)
parms.k5 = parms.insulin_sensitivity ;   % insulin sensitivity (fit from medronic IOB table- amount BG is lowered by insulin)
parms.k6 = .5 ;    % rate at which liver drips glucose into blood plasma (glucose/min)
parms.low_limit = 70 ;
estimparms.high_limit = 180 ;
parms.basal = parms.k6/parms.insulin_sensitivity;   % basal insulin
parms.insulin = parms.basal*ones(size(parms.time_series)) ;


% compute derivatives of change in basals
parms.carb_delay = 0 ;
parms.carb_grams = 0 ;
parms.G0 = parms.G_target ;
baseline = simulate_timecourse(parms) ;
del = .1 ;
endt = max(parms.time_series) ;
for t=1:endt-1
    delta = parms.insulin(t)*del ;
    parms.insulin(t) = parms.insulin(t) - delta ;
    parms.insulin(t+1) = parms.insulin(t+1) + delta ;
    outputs = simulate_timecourse(parms) ;
    derivatives(t,:) = (outputs.Gp_t - baseline.Gp_t)/del ;
    parms.insulin(t) = parms.insulin(t) + delta ;
    parms.insulin(t+1) = parms.insulin(t+1) - delta ;
end


parms.carb_delay = 30 ;
parms.carb_grams = 50 ;
parms.G0 = 120 ;
bolus = parms.carb_grams/outputs.carb_ratio;                    % insulin to counteract carbs
bolus = bolus + (parms.G0-parms.G_target)/parms.insulin_sensitivity ; % insulin to correct for current BG level

parms.insulin(parms.insulin_delay) = parms.insulin(parms.insulin_delay) + bolus ;

baseline = simulate_timecourse(parms) ;
baseline_insulin = parms.insulin ;

plot(parms.time_series, baseline.Gp_t) ; hold on ;
ln = line([0 max(parms.time_series)],[parms.G_target parms.G_target]);
set(ln, 'linestyle', '-.') ; hold off ;
disp(sprintf('mean BG = %2.1f (min = %2.1f, max = %2.1f)', mean(baseline.Gp_t),min(baseline.Gp_t),max(baseline.Gp_t))) ;

%parms.insulin(parms.insulin_delay) = parms.basal ;
%parms.insulin(1:60) = parms.insulin(1:60) + bolus/60 ;

min_rms = compute_BG_rms(parms.G_target, baseline.Gp_t) ;
prev_rms = min_rms ;
dt = .000001 ;
total_insulin = sum(parms.insulin) ;
disp(sprintf('iter %4.4d: rms %2.1f',0,min_rms));
for iter=1:1000
    outputs = simulate_timecourse(parms) ;
    start_rms = compute_BG_rms(parms.G_target, outputs.Gp_t) ;
    [mx,mxt]= max(outputs.Gp_t) ;
    [mn,mnt]= min(outputs.Gp_t) ;
    ind = find(outputs.Gp_t > parms.high_limit) ;
    if (mn > parms.low_limit+5 & length(ind) > 0)
    if (0)
       end_high=ind(end) ;
       for t=end_high:-1:2
         delta = parms.insulin(t)*.1 ;
	 parms.insulin(t) = parms.insulin(t) - delta ;
	 parms.intuslin(t-1) = parms.insulin(t-1)+delta ;
       end
    else
      for t=1:endt-1
        dI = 0 ;
        for tn=t+1:endt
	    error = (parms.G_target - outputs.Gp_t(t)) ;
            dI = dI - error * derivatives(t, tn) ;
        end
        dI = dI/(endt-t) ;
        dI_t(t) = dI ;
        parms.insulin(t) = parms.insulin(t) + dI * dt ;
      end
    end
    ind = find(parms.insulin < 0) ;
    parms.insulin(ind) = zeros(size(ind)) ;
    parms.insulin = parms.insulin * total_insulin / sum(parms.insulin) ;
    outputs = simulate_timecourse(parms) ;
    end_rms = compute_BG_rms(parms.G_target, outputs.Gp_t) ;
    if (mod(iter,10) ==0)
        pct_change = 100*(prev_rms - end_rms) / (prev_rms) ;
        disp(sprintf('iter %4.4d: rms %2.1f --> %2.1f (%2.1f%%)',iter,prev_rms, end_rms,pct_change));
	prev_rms = end_rms ;
    end
end

plot(baseline.Gp_t, 'r') ;
hold on ;
plotyy(parms.time_series, outputs.Gp_t, parms.time_series,parms.insulin) ; 
ln=line([0 endt], [parms.G_target parms.G_target]);
set(ln, 'linestyle', '-.') ; hold off ;
ln=line([0 endt], [70 70]);
set(ln, 'linestyle', '-.', 'color', 'r') ; hold off ;
hold off; 


if 0
	    if (outputs.Gp_t(t) < parms.low_limit)
	       error = error* 10 ;
	    elseif (outputs.Gp_t(t) < parms.G_target)
               error = 0 ;
	    elseif (outputs.Gp_t(t) > parms.high_limit)
	       error = error* 2 ;
            end
end
