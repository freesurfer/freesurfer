function [BG_t] = estimate_BG_probability(true_parms, hyper_parms)
% function [BG_t] = estimate_BG_probability(true_parms, hyper_parms)

BG_t = zeros(true_parms.max_BG, true_parms.T1/true_parms.dt) ;

nsamples = true_parms.nsamples ;


for n=1:nsamples
    parms = true_parms ;
    parms.basal = parms.basal*(1+hyper_parms.basal_std*true_parms.rand(n,1)) ;
    parms.insulin_sensitivity = parms.insulin_sensitivity * (1+hyper_parms.insulin_sensitivity_std*true_parms.rand(n,2));	
    parms.carb_sensitivity = parms.carb_sensitivity * (1+hyper_parms.carb_sensitivity_std*true_parms.rand(n,3));	
    parms.carb_delay(1) = parms.carb_delay(1) * (1+hyper_parms.carb_delay_std*true_parms.rand(n,4));	
    parms.carb_grams(1) = parms.carb_grams(1) * (1+hyper_parms.carb_grams_std*true_parms.rand(n,5));	
    parms.G0 = parms.G0 * (1+hyper_parms.G0_std*true_parms.rand(n,6));	
    carb_ratio = parms.insulin_sensitivity/parms.carb_sensitivity ;

%    timecourse.Gp_t = compute_BG_from_insulin(parms.insulin, parms) ;
    timecourse.Gp_t = compute_Upsampled(parms.insulin, parms) ;
    for t=1:true_parms.T1/true_parms.dt
    	ind = min(parms.max_BG, max(1,round(timecourse.Gp_t(t)))) ;
    	BG_t(ind,t) = BG_t(ind,t) + 1 ;
    end
end
BG_t = BG_t ./ nsamples ;






if 0 



true_parms = parms ;
parameter_error = .1 ;
%parms.insulin_sensitivity = parms.insulin_sensitivity * (1+parameter_error);
%carb_ratio = carb_ratio * (1+parameter_error) ;
parms.basal = parms.basal * (1-parameter_error) ;

parms.insulin = parms.basal*ones(size(parms.time_series)) ;
bolus = parms.carb_grams/carb_ratio;                    % insulin to counteract carbs
bolus = bolus + (parms.G0-parms.G_target)/parms.insulin_sensitivity ; % insulin to correct for current BG level
parms.insulin(parms.insulin_delay) = parms.insulin(parms.insulin_delay) + bolus ;

true_timecourse.Gp_t = simulate_timecourse(parms.insulin, parms) ;
tt = true_timecourse ;


measured_timecourse = true_timecourse;
measured_timecourse.Gp_t = true_timecourse.Gp_t + randn(1,size(true_timecourse.Gp_t,2)).*sigma.*true_timecourse.Gp_t;

times = [120:10:T2];
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




end
