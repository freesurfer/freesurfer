error_type = 1 ;
parms.time_series=1:5*60 ;

parms.G_target = 100 ;   % desired blood glucose (Gp) level
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
parms.high_limit = 180 ;
parms.basal = parms.k6/parms.insulin_sensitivity;   % basal insulin
parms.insulin = parms.basal*ones(size(parms.time_series)) ;


% compute derivatives of change in basals
parms.carb_delay = 0 ;
parms.carb_grams = 0 ;
parms.G0 = parms.G_target ;
del_insulin = parms.basal/2 ;

baseline = simulate_timecourse(parms) ;
endt = max(parms.time_series) ;
for t=1:0
    for t2=1:endt
        parms.insulin(t) = parms.insulin(t) - del_insulin ;
	parms.insulin(t2) = parms.insulin(t2) + del_insulin ;
	outputs = simulate_timecourse(parms) ;
	derivatives(t,t2,:) = (outputs.Gp_t - baseline.Gp_t) ;
        parms.insulin(t) = parms.insulin(t) + del_insulin ;
	parms.insulin(t2) = parms.insulin(t2) - del_insulin ;
   end
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

min_rms = compute_BG_rms(parms.G_target, baseline.Gp_t,error_type) ;
prev_rms = min_rms ;
total_insulin = sum(parms.insulin) ;
disp(sprintf('iter %4.4d: rms %2.1f',0,min_rms));
for iter=1:20000
    outputs = simulate_timecourse(parms) ;
    start_rms = compute_BG_rms(parms.G_target, outputs.Gp_t,error_type) ;
    best_rms = start_rms ;
    best_t1 = -1 ;
    best_t2 = -1 ;
    for t1=1:endt
        for t2=t1+1:endt
            Gp_t = outputs.Gp_t + squeeze(derivatives(t1,t2,:))' ;
            Gp_t_r = outputs.Gp_t - squeeze(derivatives(t1,t2,:))' ;
	    rms = compute_BG_rms(parms.G_target, Gp_t,error_type) ;
            del_rms = start_rms - rms ;
            if (rms < best_rms & (min(Gp_t) > parms.low_limit))
                best_t1 = t1 ;
                best_t2 = t2 ;
                best_rms = rms ;
            elseif (start_rms + del_rms < best_rms & (min(Gp_t_r) > parms.low_limit))  % transfer in other direction is best
                best_t1 = t2 ;
                best_t2 = t1 ;
                best_rms = start_rms+del_rms ;
            end
        end
    end
    if (best_t1 > 0)
        parms.insulin(best_t1) = parms.insulin(best_t1) - del_insulin ;
        parms.insulin(best_t2) = parms.insulin(best_t2) + del_insulin ;
    end
    
    outputs = simulate_timecourse(parms) ;
    end_rms = compute_BG_rms(parms.G_target, outputs.Gp_t,error_type) ;
    if (mod(iter,10) ==0)
        pct_change = 100*(prev_rms - end_rms) / (prev_rms) ;
        disp(sprintf('iter %4.4d: rms %2.2f --> %2.2f (%2.1f%%), mean BG = %2.1f',iter,prev_rms, end_rms,pct_change,mean(outputs.Gp_t)));
        prev_rms = end_rms ;
	if (pct_change < .1)
	   break ;
	end
    end
end

h1 = plot(baseline.Gp_t, 'r') ;
hold on ;
h2 = plotyy(parms.time_series, outputs.Gp_t, parms.time_series,parms.insulin) ; 
ln=line([0 endt], [parms.G_target parms.G_target]);
set(ln, 'linestyle', '-.','linewidth',6) ; 
ln=line([0 endt], [parms.low_limit parms.low_limit]);
set(ln, 'linestyle', '-.', 'color', 'r','linewidth',6) ; 
hold off; 
set(gca, 'fontsize', 16, 'fontweight', 'bold') ;
set(h2(2), 'fontsize', 16, 'fontweight', 'bold') ;
set(get(gca,'children'), 'linewidth', 6) 
set(get(h2(1),'children'), 'linewidth', 6) 
set(get(h2(2),'children'), 'linewidth', 6) 
ylabel('blood glucose level (mg/dl)', 'fontsize', 16, 'fontweight', 'bold');
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
set(h2(2), 'fontsize', 16, 'fontweight', 'bold') ;
set(get(h2(2),'ylabel'), 'string', 'Insulin (U)', 'fontsize', 16, 'fontweight', 'bold')
set(gca, 'ytick', [100:100:300],'ylim',[0 325])
legend('normal BG response', 'optimal BG response', 'normal BG level', 'low BG level', 'optimal insulin')    
disp(sprintf('mean BG = %2.1f (min = %2.1f, max = %2.1f)', mean(outputs.Gp_t),min(outputs.Gp_t),max(outputs.Gp_t))) ;

if 0

last_dI_t = zeros(size(last_dI_t)) ;
for iter=1:220
    outputs = simulate_timecourse(parms) ;
    start_rms = compute_BG_rms(parms.G_target, outputs.Gp_t,error_type) ;
    for t=1:20
        dI = 0 ;
        for tn=t+1:endt
	    error = -(parms.G_target - outputs.Gp_t(tn)) ;
            dI = dI + error * derivatives(t, tn) ;
%	    disp(sprintf('del(%2.0f) = %2.1f * %2.5f = %2.3f', tn,error, derivatives(t,tn), error * derivatives(t, tn)));
        end
        dI = dI/(endt-t) ;
        dI_t(t) = dI ;
        parms.insulin(t) = parms.insulin(t) + dI * dt ;
    end
    dI_t = dI_t + momentum*last_dI_t ;
%    parms.insulin = parms.insulin + dI_t * dt ;
    last_dI_t = dI_t ;
    ind = find(parms.insulin < 0) ;
    parms.insulin(ind) = zeros(size(ind)) ;
    parms.insulin = parms.insulin * total_insulin / sum(parms.insulin) ;
    outputs = simulate_timecourse(parms) ;
    end_rms = compute_BG_rms(parms.G_target, outputs.Gp_t,error_type) ;
    if (mod(iter,10) ==0)
        pct_change = 100*(prev_rms - end_rms) / (prev_rms) ;
        disp(sprintf('iter %4.4d: rms %2.1f --> %2.1f (%2.6f%%)',iter,prev_rms, end_rms,pct_change));
	prev_rms = end_rms ;
    end
end
disp(sprintf('mean BG = %2.1f (min = %2.1f, max = %2.1f)', mean(outputs.Gp_t),min(outputs.Gp_t),max(outputs.Gp_t))) ;

end
