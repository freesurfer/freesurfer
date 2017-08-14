T1 = 300;
parms.G0 = 200 ;
T2 = 600;
parms.G_target = 100;   % desired blood glucose (Gp) level
low_limit = 70;
high_limit = 180;
parms.insulin_sensitivity = 180 ;  % glucose/unit insulin
parms.k1 = .021 ;   % (old k4) rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
parms.k2 = .001 ;   % (old k3) rate at which insulin moves from plasma to interstitial fluid (units/min)
parms.k3 = .021 ;   % (old k5) rate at which insulin moves from fluid to plasma (units/min)
parms.k4 = .021 ;    % (old k1) rate at which carbs are metabolized from stomach to blood (grams/min)
parms.k5 = parms.insulin_sensitivity ;   % insulin sensitivity (fit from medronic IOB table- amount BG is lowered by insulin)

parms.k6 = .5 ;    % (old k2) rate at which liver drips glucose into blood plasma (glucose/min)
parms.k7 = .05;    % (old k6) relative rate at which muscle glucose is refilled (1/min)
parms.low_limit = low_limit ;
parms.high_limit = high_limit ;
parms.carb_sensitivity = 7 ;
parms.time_series = 1:T2;
parms.basal = parms.k6/parms.insulin_sensitivity;   % basal insulin 
parms.carb_grams = 80 ;
parms.carb_delay = 16 ;
parms.insulin_delay = 1 ;
carb_ratio = parms.insulin_sensitivity/parms.carb_sensitivity ;

parameter_error = .1 ;
%parms.insulin_sensitivity = parms.insulin_sensitivity * (1+parameter_error);
%carb_ratio = carb_ratio * (1+parameter_error) ;
parms.basal = parms.basal * (1-parameter_error) ;

parms.insulin = parms.basal*ones(size(parms.time_series)) ;
bolus = parms.carb_grams/carb_ratio;                    % insulin to counteract carbs
bolus = bolus + (parms.G0-parms.G_target)/parms.insulin_sensitivity ; % insulin to correct for current BG level
parms.insulin(parms.insulin_delay) = parms.insulin(parms.insulin_delay) + bolus ;

true_parms = parms ;
true_timecourse = simulate_timecourse(parms) ;
tt = true_timecourse ;

sigma = .05 ;  % if the 95% of readings are within 10% of BG this is about right
measured_timecourse = true_timecourse;
measured_timecourse.Gp_t = true_timecourse.Gp_t + randn(1,size(true_timecourse.Gp_t,2)).*sigma.*true_timecourse.Gp_t;

times = [120:10:T1];
for i=1:length(times)
    t = times(i) ;
    dt(i) = t ;
    A(i,1) = true_timecourse.Ip_t(1)+bolus ;
    A(i,2) = true_timecourse.Ip_t(t) ;
    M(i,1) = true_timecourse.Gp_t(1) ;	
    M(i,2) = true_timecourse.Gp_t(t) ;	
    F(i) = parms.carb_grams ;
end

bg = predict_bg(true_timecourse.Gp_t(1), parms.insulin_sensitivity, A(1,1)-A(1,2), 0,times(1),parms.carb_grams,parms.carb_sensitivity);

[Ibest, Sbest, ebest, Cbest] = calibrate_pump(A, M, F, dt) ;

nsamples = 500 ;
for n=1:nsamples
    measured_timecourse.Gp_t = true_timecourse.Gp_t + randn(1,size(true_timecourse.Gp_t,2)).*sigma.*true_timecourse.Gp_t;
    input_timecourse(n,:) = measured_timecourse.Gp_t ;
    for i=1:length(times)
        t = times(i) ;
    	dt(i) = t ;
        A(i,1) = measured_timecourse.Ip_t(1)+bolus ;
        A(i,2) = measured_timecourse.Ip_t(t) ;
        M(i,1) = measured_timecourse.Gp_t(1) ;	
        M(i,2) = measured_timecourse.Gp_t(t) ;	
        F(i) = parms.carb_grams ;
    end
    [Ibest, Sbest, ebest, Cbest] = calibrate_pump(A, M, F, dt) ;

    noisy_parms(n) = parms ;
    noisy_parms(n).insulin_sensitivity = Sbest ;
    noisy_parms(n).k5 = Sbest ;
    noisy_parms(n).carb_sensitivity = Ibest ;
    noisy_parms(n).basal = parms.basal + ebest ;
    noisy_carb_ratio(n) = noisy_parms(n).insulin_sensitivity/noisy_parms(n).carb_sensitivity ;
    noisy_bolus(n) = noisy_parms(n).carb_grams/noisy_carb_ratio(n);                    % insulin to counteract carbs
    noisy_bolus(n) = noisy_bolus(n) + (noisy_parms(n).G0-noisy_parms(n).G_target)/noisy_parms(n).insulin_sensitivity ; % insulin to correct for current BG level

    parms.insulin(parms.insulin_delay) = parms.basal + noisy_bolus(n) ;
    noisy_timecourse(n) = simulate_timecourse(parms) ;
    nS(n) = Sbest ; nI(n) = Ibest ; ne(n) = ebest ; nC(n) = Cbest ;
    ncr(n,:) = noisy_carb_ratio(n) ;
    nt(n,:) = noisy_timecourse(n).Gp_t ;
end

% convert to percentages
pnc = abs(nC - carb_ratio) / carb_ratio ;

mean_timecourse = mean(nt) ;
std_timecourse = std(nt) ;
confidence95 = 2*std_timecourse ;
plot(mean_timecourse, 'g')  ;
hold on ;
plot(mean_timecourse+confidence95, 'r--')  ;
%plot(tps, mt(tps), 'r--')  ;
plot(mean_timecourse-confidence95, 'r--')  ;
ch = get(gca,'children') ;
set(ch(1:2), 'linewidth', 2);
set(ch(3), 'linewidth', 5);
set(gca, 'fontsize', 16, 'fontweight', 'bold') ;
hold off;
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('blood glucose (mg/dL)', 'fontsize', 16, 'fontweight', 'bold') ;
