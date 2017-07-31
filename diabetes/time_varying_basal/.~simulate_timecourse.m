time_series=1:5*60 ;

G_target = 120 ;   % desired blood glucose (Gp) level
carb_sensitivity = 7 ;
insulin_sensitivity = 180 ;  % glucose/unit insulin
carb_ratio = insulin_sensitivity/carb_sensitivity ;


carb_delay = 20 ;
carb_grams = 50 ;
steady_state_basal_in_plasma = 0.13214 ;
steady_state_basal_in_fluid = 0.13844 ;

%uncomment these to disable basal
%steady_state_basal_in_plasma = 0 ;
%steady_state_basal_in_fluid = 0 ;
%parms.k6=0 ;


parms.k1 = .021 ;   % rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
parms.k1 = .03 ;   % rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
parms.k2 = .001 ;   % rate at which insulin moves from plasma to interstitial fluid (units/min)
parms.k3 = .021 ;   % rate at which insulin moves from fluid to plasma (units/min)
parms.k4 = .05 ;    % rate at which carbs are metabolized from stomach to blood (grams/min)
parms.k5 = insulin_sensitivity ;   % insulin sensitivity (fit from medronic IOB table- amount BG is lowered by insulin)
parms.k6 = .5 ;    % rate at which liver drips glucose into blood plasma (glucose/min)

Ib = .001852*parms.k6*insulin_sensitivity/60;   % basal insulin
Ib_orig = Ib ;

Gs = carb_grams*carb_sensitivity;   % glucose in stomach
bolus = carb_grams/carb_ratio;                   % insulin in interstitial fluid

Gred = 0 ;
Ginc = 0 ;
%Ip = 0 ;                    % insulin in plasma
Gs = 0 ;

Ip = steady_state_basal_in_plasma ;
If = bolus+steady_state_basal_in_fluid ;

Gp = G_target ;                  % glucose in bloodparms.k6
Iu = 0 ;                    % total insulin used by cells

for  t=time_series
    if (t == carb_delay & Gs == 0)
      Gs =  carb_grams*carb_sensitivity;
    end
    if (t > 0 & t <= 20)
       Ib = 4*Ib_orig ;
    elseif (t > 150 & t <= 230)
       Ib = 0 ;
    else
	Ib = Ib_orig ;
    end
    Gs_t (t) = Gs ;
    Gp_t (t) = Gp ;
    Ip_t (t) = Ip ;
    If_t (t) = If ;
    Iu_t(t) = Iu ;

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

plot(time_series, Gp_t) ; hold on ;
ln = line([0 max(time_series)],[G_target G_target]);
set(ln, 'linestyle', '-.') ; hold off ;
disp(sprintf('mean BG = %2.1f', mean(Gp_t))) ;
