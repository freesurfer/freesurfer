% Ib - Insulin in the blood
% If - Insulin in interstilial fluid
% BG - blood glucose level
% C - Carb ratio (number of carbs accounted for my one unit of insulin)
% Is - Insulin sensitivity (amount BG drops per unit of insulin)
% Fs - Food sensitivity - amount BG rises per gram of carb
% carbs - time course of carb intake


dt = 1 ;
time=1:dt:4*60;

insulin_sensitivity = 150 ;
carb_sensitivity = 7 ;
carb_ratio = insulin_sensitivity / carb_sensitivity ;

insulin_track = zeros(size(time)) ;
carbs_track = zeros(size(time)) ;
blood_glucose_track = zeros(size(time)) ;

BG(1) = 120 ;

carbs_stomach(1) = 40 ;
Kcs = 3 / 60 ;
Ki =  2 /(60) ;
carbs = 40 ;
insulin = carbs / carb_ratio ;
bolus = insulin ;
estimateIOBTable ;
[Ki, best_it] = fit_insulin_timecourse(iob_table) ;

for t=time
    delta_carb = carbs * Kcs * dt ;
    carbs = carbs - delta_carb ;
    carb_track(t) = carbs ;

    delta_insulin = insulin * Ki * dt ;
    if (t > length(iob_change_table))
       delta_insulin = 0 ;
    else
        delta_insulin = bolus*iob_change_table(t)/100 ;
    end
    BG = BG + (carb_sensitivity*delta_carb) - (insulin_sensitivity*delta_insulin) ;
    BG_track(t) = BG ;
    insulin = insulin - delta_insulin ;
    insulin_track(t) = insulin ;
end

subplot(211)
plotyy(time, carb_track, time, insulin_track) ;
subplot(212)
plot(time, BG_track) ;
