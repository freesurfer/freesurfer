function loss = compute_insulin_schedule_loss(parms, hyper_parms)
% function loss = compute_insulin_schedule_loss(parms, hyper_parms)
BG_t = estimate_BG_probability(parms, hyper_parms);
