function outputs = simulate_timecourse(parms)

% outputs = simulate_timecourse(parms)
% parms has:
% parms.k1 = .021 ;   % rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
% parms.k1 = .03 ;   % rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
% parms.k2 = .001 ;   % rate at which insulin moves from plasma to interstitial fluid (units/min)
% parms.k3 = .021 ;   % rate at which insulin moves from fluid to plasma (units/min)
% parms.k4 = .05 ;    % rate at which carbs are metabolized from stomach to blood (grams/min)
% parms.k5 = insulin_sensitivity ;   % insulin sensitivity (fit from medronic IOB table- amount BG is lowered by insulin)
% parms.k6 = .5 ;    % rate at which liver drips glucose into blood plasma (glucose/min)
% parms.carb_delay = 20 ;
% parms.insulin_delay = 1 ;
% parms.carb_grams = 50 ;
% parms.time_series=1:5*60 ;
% parms.G_target = 120 ;   % desired blood glucose (Gp) level
% parms.carb_sensitivity = 7 ;
% parms. insulin_sensitivity = 180 ;  % glucose/unit insulin
% parms.G0 = starting BG level

outputs.carb_ratio = parms.insulin_sensitivity/parms.carb_sensitivity ;

Ip = 0 ;   % insulin in plasma
If = 0 ;   % insulin in interstitial fluid

% first estimate steady state basal insulin levels in plasma and in
% fluid assuming current basal has been constant for a while
for  t=1:parms.T1
    % deltas are all positive
    delta_Ip_to_cells = parms.k1*Ip ;
    delta_If_to_blood = parms.k3*If ;
    delta_Ip_to_fluid = parms.k2*Ip ;

    Ip = Ip -delta_Ip_to_cells - delta_Ip_to_fluid + delta_If_to_blood ;
    If = If + delta_Ip_to_fluid - delta_If_to_blood + parms.basal ;
end

Gred = 0 ;
Ginc = 0 ;
if (parms.carb_delay == 0)
   Gs = parms.carb_grams*parms.carb_sensitivity;   % glucose in   stomach
else
   Gs = 0 ;
end

Iu = 0 ;                    % total insulin used by cells
Gp = parms.G0 ;

for  t=1:parms.T1
    if (t == parms.carb_delay & Gs == 0)
      Gs =  parms.carb_grams*parms.carb_sensitivity;
    end
    Ib = parms.insulin(t) ;
    outputs.Gs_t (t) = Gs ;
    outputs.Gp_t (t) = Gp ;
    outputs.Ip_t (t) = Ip ;
    outputs.If_t (t) = If ;
    outputs.Iu_t(t) = Iu ;

    % deltas are all positive
    delta_Ip_to_cells = parms.k1*Ip ;
    delta_glucose_to_blood = parms.k4*Gs ;
    delta_If_to_blood = parms.k3*If ;
    delta_Ip_to_fluid = parms.k2*Ip ;

    Ip = Ip -delta_Ip_to_cells - delta_Ip_to_fluid + delta_If_to_blood ;
    If = If + delta_Ip_to_fluid - delta_If_to_blood + Ib ;
    Gs = Gs - delta_glucose_to_blood ;
    Iu = Iu + delta_Ip_to_cells ;
    Gred = Gred + delta_Ip_to_cells*parms.k5 ;
    Ginc = Ginc + delta_glucose_to_blood ;
    Gp = Gp + delta_glucose_to_blood + (-delta_Ip_to_cells*parms.k5) + parms.k6 ;
end

