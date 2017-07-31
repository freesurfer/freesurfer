function loss = compute_insulin_schedule_loss(insulin, parms, hyper_parms)
% function loss = compute_insulin_schedule_loss(insulin, parms, hyper_parms)
parms.insulin = insulin ;
BG_t = estimate_BG_probability(parms, hyper_parms);

loss = compute_loss(BG_t, parms.G_target, parms.risk_matrix) ;
