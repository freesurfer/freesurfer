parms.insulin = parms.basal*ones(size(parms.time_series)) ;
bolus = parms.carb_grams/carb_ratio;                    % insulin to counteract carbs
bolus = bolus + (parms.G0-parms.G_target)/parms.insulin_sensitivity ; % insulin to correct for current BG level
parms.insulin(parms.insulin_delay) = parms.insulin(parms.insulin_delay) + bolus ;
