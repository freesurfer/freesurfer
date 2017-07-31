parms.dt = 5 ;
parms.G0 = 100 ;
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
parms.T1 = 5*60 ;
parms.T2 = 10*60 ;
parms.time_series = 1:parms.dt:parms.T2;

parms.carb_sensitivity = 7 ;
parms.basal = parms.k6/parms.insulin_sensitivity;   % basal insulin 
parms.carb_grams = 80 ;
parms.carb_delay = 61 ;
parms.insulin_delay = 1 ;
parms.max_BG = 700;
carb_ratio = parms.insulin_sensitivity/parms.carb_sensitivity ;


parms.nsamples = 1000 ;

hyper_parms.G0_std = .0 ;   % if the 95% of readings are within 10% of BG this is about right
hyper_parms.insulin_sensitivity_std = .0 ;
hyper_parms.carb_sensitivity_std = .0 ;
hyper_parms.carb_std = .0 ;

hyper_parms.G0_std = .05 ;   % if the 95% of readings are within 10% of BG this is about right
hyper_parms.insulin_sensitivity_std = .05 ;
hyper_parms.carb_sensitivity_std = .05 ;
hyper_parms.carb_std = .05 ;
hyper_parms.basal_std = .05 ;
hyper_parms.basal_std = .01 ;    % reduce basal uncertainty
hyper_parms.carb_delay_std = 0 ;
hyper_parms.carb_grams_std = 0 ;

if 0
%make them small for testing
hyper_parms.G0_std = .025 ;
hyper_parms.insulin_sensitivity_std = .025 ;
hyper_parms.carb_sensitivity_std = .025 ;
hyper_parms.carb_std = .0 ;
hyper_parms.basal_std = .01 ;
end

erate = 5 ;
Glow = 50 ;
sigma = 20;
escale = 1000;
quadratic_loss_cost = 10 ;

G = 1:parms.max_BG;
test_sigma = 1000;
loss_vec = ((G-parms.G_target).^2./sigma.^2)';
ind = find(loss_vec>140);
loss_vec(ind) = loss_vec(ind)*5 ;
eg = escale*exp(-(G-Glow)/erate)' ;
loss_vec = (quadratic_loss_cost*loss_vec+eg) ;
% no_test_matrix  = repmat(no_test_prob, [parms.max_BG 1]);
parms.risk_matrix = repmat(loss_vec, [1 parms.T1]) ;

%parms.risk_matrix = repmat(loss_vec, [1 parms.T1/parms.dt]) ;
%no_test_prob = exp(-([1:parms.dt:parms.T1].^2) ./ (2*test_sigma.^2));
%parms.risk_matrix = repmat(loss_vec, [1 parms.T1/parms.dt]) .* no_test_matrix;


parms.nsamples = 1000 ;
parms.rand = randn(parms.nsamples, 6) ;

[BG_init,insulin_init] = compute_BG_and_insulin(parms) ;
parms.insulin = insulin_init ;
BG_prob = estimate_BG_probability(parms, hyper_parms);
loss_init = compute_insulin_schedule_loss(insulin_init, parms, hyper_parms)


insulin_opt = minimize_loss_functional(insulin_init, parms, hyper_parms,10,10);
save('Qopt.mat') ;
insulin_opt = minimize_loss_functional(insulin_opt, parms, hyper_parms,10,1);
save('Qopt2.mat') ;
insulin_opt2 = minimize_loss_functional(insulin_opt, parms, hyper_parms,10,.1);
insulin_opt = insulin_opt2 ;
save('Qopt2.mat') ;
insulin_opt3 = minimize_loss_functional(insulin_opt2, parms, hyper_parms,10,.1);
insulin_opt = insulin_opt3 ;
save('Qopt3.mat') ;
insulin_opt4 = minimize_loss_functional(insulin_opt3, parms, hyper_parms,10,.1);
save('Qopt4.mat') ;
insulin_opt5 = minimize_loss_functional(insulin_opt4, parms, hyper_parms,10,.1);
insulin_opt = insulin_opt5;
save('Qopt5.mat') ;
insulin_opt6 = minimize_loss_functional(insulin_opt5, parms, hyper_parms,10,.1);

insulin_opt = insulin_opt6;
save('Qopt6.mat') ;
insulin_opt7 = minimize_loss_functional(insulin_opt6, parms, hyper_parms,10,.1);
save('Qopt7.mat') ;
insulin_opt8 = minimize_loss_functional(insulin_opt7, parms, hyper_parms,10,1); 
insulin_opt = insulin_opt8 ;
save('Qopt8.mat') ;
insulin_opt9 = minimize_loss_functional(insulin_opt8, parms, hyper_parms,40,1); 
insulin_opt = insulin_opt9 ;
insulin_opt = minimize_loss_functional(insulin_opt9, parms, hyper_parms,40,.1); 
save('Qopt9.mat') ;

parms.insulin = insulin_opt ;
BGopt_prob = estimate_BG_probability(parms, hyper_parms);
clear('BG_opt') ;
%BG_opt = compute_BG_from_insulin(insulin_opt, parms) ;

plot_uncertainty_and_timecourses

if 0
figure('name', 'insulin schedules');
plot(insulin_opt - insulin_init, 'g') ;
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('insulin (A/U)', 'fontsize', 16, 'fontweight', 'bold') ;
set(gca, 'xlim', [0 300/parms.dt], 'fontsize', 16, 'fontweight', 'bold');
ch = get(gca, 'children') ;
set(ch, 'linewidth', 4);
title('risk minimizing - deterministic insulin schedule', 'fontsize', 16, 'fontweight', 'bold');
end

if 0
carb_ratio = parms.insulin_sensitivity/parms.carb_sensitivity ;
parms.insulin = parms.basal*ones(size(parms.time_series))' ;
bolus = sum(insulin_init) - sum(parms.insulin);
%bolus = parms.carb_grams/carb_ratio;                    % insulin to counteract carbs
%bolus = bolus + (parms.G0-parms.G_target)/parms.insulin_sensitivity ; % insulin to correct for current BG level
insulin_index = round((parms.carb_delay-parms.insulin_delay)/parms.dt);
parms.insulin(insulin_index) = parms.insulin(insulin_index) + bolus ;
insulin_fixed = parms.insulin ;
BG_fixed = computeUpsampled(insulin_fixed, parms);

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

if 0

fminopts = optimset('fminsearch') ;
fminopts = optimset(fminopts, 'Display', 'iter','MaxIter', 1000);
insulin_opt = fminsearch(@compute_insulin_schedule_loss, insulin_init, fminopts,parms, hyper_parms) ;

total_insulin = sum(insulin_init(1:parms.T2/parms.dt)) ;
fconopts = optimset('fmincon') ;
fconopts = optimset(fconopts, 'Display', 'iter','MaxIter', 1000, 'MaxFunEvals',20000,'Algorithm','interior-point','tolcon',1);
fconopts = optimset(fconopts, 'Display', 'iter','MaxIter', 1000, 'MaxFunEvals',20000,'Algorithm','sqp','tolcon',1,'tolx',1e-15);
fconopts = optimset(fconopts, 'Display', 'iter','MaxIter', 1000, 'MaxFunEvals',20000,'Algorithm','active-set','tolcon',1,'tolx',1e-15);
[insulin_opt,loss_opt,eflag] = fmincon(@compute_insulin_schedule_loss, insulin_init, [], [], ones(1,length(insulin_init)),total_insulin, zeros(size(insulin_init)),[],[],fconopts,parms, hyper_parms) ;

% remove total insulin constraint
[insulin_opt,fminval,eflag] = fmincon(@compute_insulin_schedule_loss,insulin_init, [], [], [], [], zeros(size(insulin_init)),[],[],fconopts,parms, hyper_parms) ;

% remove non-negativeity constraint
[insulin_opt,fminval,eflag] = fmincon(@compute_insulin_schedule_loss, insulin_init, [], [], ones(1,length(insulin_init)),total_insulin, [],[],[],fconopts,parms, hyper_parms) ;

end

