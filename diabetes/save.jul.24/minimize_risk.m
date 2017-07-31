

parms.G0 = 200 ;
parms.time_series = 1:600;
parms.G_target = 100;   % desired blood glucose (Gp) level
parms.low_limit = 70;
parms.high_limit = 180;
parms.insulin_sensitivity = 180 ;  % glucose/unit insulin
parms.k1 = .021 ;   % (old k4) rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
parms.k2 = .001 ;   % (old k3) rate at which insulin moves from plasma to interstitial fluid (units/min)
parms.k3 = .021 ;   % (old k5) rate at which insulin moves from fluid to plasma (units/min)
parms.k4 = .021 ;    % (old k1) rate at which carbs are metabolized from stomach to blood (grams/min)
parms.k5 = parms.insulin_sensitivity ;   % insulin sensitivity (fit from medronic IOB table- amount BG is lowered by insulin)
parms.k6 = .5 ;    % (old k2) rate at which liver drips glucose into blood plasma (glucose/min)
parms.k7 = .05;    % (old k6) relative rate at which muscle glucose is refilled (1/min)
parms.T1 = 300 ;
parms.T2 = 600 ;

parms.carb_sensitivity = 7 ;
parms.basal = parms.k6/parms.insulin_sensitivity;   % basal insulin 
parms.carb_grams = 80 ;
parms.carb_delay = 31 ;
parms.insulin_delay = 1 ;
parms.max_BG = 700;
carb_ratio = parms.insulin_sensitivity/parms.carb_sensitivity ;

parms.nsamples = 1000 ;

hyper_parms.G0_std = .0 ;   % if the 95% of readings are within 10% of BG this is about right
hyper_parms.insulin_sensitivity_std = .0 ;
hyper_parms.carb_sensitivity_std = .0 ;
hyper_parms.carb_std = .0 ;
hyper_parms.basal_std = .0 ;

hyper_parms.G0_std = .05 ;   % if the 95% of readings are within 10% of BG this is about right
hyper_parms.insulin_sensitivity_std = .05 ;
hyper_parms.carb_sensitivity_std = .05 ;
hyper_parms.carb_std = .05 ;
hyper_parms.basal_std = .05 ;
hyper_parms.carb_delay_std = 0 ;
hyper_parms.carb_grams_std = 0 ;


erate = 10 ;
Glow = 50 ;
sigma = 20;
escale = 1000;

G = 1:parms.max_BG;
loss_vec = ((G-parms.G_target).^2./sigma.^2)';
eg = escale*exp(-(G-Glow)/erate)' ;
loss_vec = (loss_vec+eg) ;
parms.risk_matrix = repmat(loss_vec, [1 parms.T1]) ;


[BG_init,insulin_init] = compute_BG_and_insulin(parms) ;
parms.insulin = insulin_init ;
BG_prob = estimate_BG_probability(parms, hyper_parms);

parms.nsamples = 1000 ;
fminopts = optimset('fminsearch') ;
fminopts = optimset(fminopts, 'Display', 'iter','MaxIter', 1000);
insulin_opt = fminsearch(@compute_insulin_schedule_loss, insulin_init(1:parms.T1), fminopts,parms, hyper_parms) ;


parms.insulin = insulin_opt ;
BGopt_prob = estimate_BG_probability(parms, hyper_parms);
BG_opt = simulate_timecourse(parms) ;

%%%%% plot the risk minimizing uncertainty
figure('name', 'risk minimizing', 'position', [786 604 1379  689]) ;
h2 = imagesc(BGopt_prob);
set(gca, 'ylim', [0 350], 'xlim', [0 300], 'fontsize', 16, 'fontweight', 'bold');
axis xy
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold') ;
colorbar;
title('uncertainty curves for risk-minimizing schedule', 'fontsize', 16, 'fontweight', 'bold') ;
print -dtiff risk_minimizing_uncertainty.tif

%%%%% plot the deteriministic uncertainty
chigh = max(BGopt_prob(:)) ;
figure('name', 'L1 minimizing', 'position', [2174   601 1248  691]) ;
h1 = imagesc(BG_prob, [0 chigh]);
set(gca, 'ylim', [0 350], 'xlim', [0 300], 'fontsize', 16, 'fontweight', 'bold');
axis xy
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold') ;
colorbar;
title('uncertainty curves for deterministic schedule', 'fontsize', 16, 'fontweight', 'bold') ;
print -dtiff deterministic_L1_uncertainty.tif

%%%% plot the timecourses
figure('name', 'BG for risk and L1-minimizing schedules') ;
plot(BG_opt.Gp_t, 'g')
hold on;
plot(BG_init, 'r')
hold off;
legend('risk minimizing', 'deterministic');
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold') ;
set(gca, 'xlim', [0 300], 'fontsize', 16, 'fontweight', 'bold');
ch = get(gca, 'children') ;
set(ch, 'linewidth', 4);
print -dtiff BG_timecourses.tif

figure('name', 'insulin schedules');
plot(insulin_opt - insulin_init, 'g') ;
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('insulin (A/U)', 'fontsize', 16, 'fontweight', 'bold') ;
set(gca, 'xlim', [0 300], 'fontsize', 16, 'fontweight', 'bold');
ch = get(gca, 'children') ;
set(ch, 'linewidth', 4);
title('risk minimizing - deterministic insulin schedule', 'fontsize', 16, 'fontweight', 'bold');

if 0
carb_ratio = parms.insulin_sensitivity/parms.carb_sensitivity ;
parms.insulin = parms.basal*ones(size(parms.time_series)) ;
bolus = parms.carb_grams/carb_ratio;                    % insulin to counteract carbs
bolus = bolus + (parms.G0-parms.G_target)/parms.insulin_sensitivity ; % insulin to correct for current BG level
parms.insulin(parms.insulin_delay) = parms.insulin(parms.insulin_delay) + bolus ;

end


if 0
k4 = .021;   % rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
k3 = .001;   % rate at which insulin moves from plasma to interstitial fluid (units/min)
k5 = .021;   % rate at which insulin moves from fluid to plasma (units/min)
k1 = .05;    % rate at which carbs are metabolized from stomach to blood (grams/min)
k1 = .02;    % rate at which carbs are metabolized from stomach to blood (grams/min)
k6 = .05;    % relative rate at which muscle glucose is refilled (1/min)
SC = 7;      % Carb sensitivity
k2 = .5 ;    % rate at which liver drips glucose into blood plasma (glucose/min)
end
