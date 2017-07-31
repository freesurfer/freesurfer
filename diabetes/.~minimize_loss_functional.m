function [insulin_opt] =  minimize_loss_functional(insulin_init,parms,hyper_parms)
%function [insulin_opt] =  minimize_loss_functional(insulin_init,parms, hyper_parms)

min_del_insulin = insulin_init(1) ;

max_t = parms.T1/parms.dt ;

insulin_opt = insulin_init;
loss_init = compute_insulin_schedule_loss(insulin_init, parms,hyper_parms);
loss_opt = loss_init ;

for step=1:100
    insulin_test = insulin_opt ;
    for t=1:max_t 
    	if (insulin_opt(t) > 0)
	   min_del = min(insulin_opt(t), min_del_insulin) ;
	   max_del = insulin_opt(t) ;
	   if (t < max_t)
	      insulin_test(t+1) = insulin_opt(t+1) + min_del ;
	      insulin_test(t) = insulin_opt(t) - min_del ;
	      loss_test = compute_insulin_schedule_loss(insulin_init, parms, hyper_parms)
	      if (loss_test < loss_opt)
	      	 insulin_opt = insulin_test ;
		 loss_opt = loss_test ;
	      end
	   end
	   if (t > 1)
	      insulin_test(t-1) = insulin_opt(t-1) + min_del ;
	      insulin_test(t) = insulin_opt(t) - min_del ;
	      loss_test = compute_insulin_schedule_loss(insulin_init, parms, hyper_parms)
	      if (loss_test < loss_opt)
	      	 insulin_opt = insulin_test ;
		 loss_opt = loss_test ;
	      end
	   end
	end
    end
    disp(sprintf('step %3.3d: loss %2.4f', step, loss_opt)) ;
end
