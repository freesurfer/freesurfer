function [insulin_opt] =  minimize_loss_functional(insulin_init,parms,hyper_parms,max_steps,del_scale)
%function [insulin_opt] =  minimize_loss_functional(insulin_init,parms, hyper_parms,max_steps,del_scale)

if (nargin < 4)
   max_steps = 1000;
end
if (nargin < 5)
   del_scale = 1 ;
end


min_del_insulin = insulin_init(1)*del_scale ;

max_t = parms.T1/parms.dt ;

insulin_opt = insulin_init;
loss_init = compute_insulin_schedule_loss(insulin_init, parms,hyper_parms);
loss_opt = loss_init ;
disp(sprintf('step %3.3d: loss %2.4f', 0, loss_opt)) ;

total_insulin = sum(insulin_opt) ;
for step=1:max_steps
    nchanged = 0 ;
    insulin_test = insulin_opt ;
    last_loss = loss_opt ;
    for t=1:max_t 
    	if (insulin_opt(t) > 0)
	   min_del = min(insulin_opt(t), min_del_insulin) ;
	   max_del = insulin_opt(t)/10 ;
	   if (t < max_t)
	      insulin_test(t+1) = insulin_opt(t+1) + min_del ;
	      insulin_test(t) = insulin_opt(t) - min_del ;
	      loss_test = compute_insulin_schedule_loss(insulin_test, parms, hyper_parms);
	      if (loss_test < loss_opt)
	      	 insulin_opt = insulin_test ;
		 nchanged = nchanged + 1 ;
		 loss_opt = loss_test ;
	      else
		insulin_test(t) = insulin_opt(t) ;
		insulin_test(t+1) = insulin_opt(t+1) ;
	      end
	   end
	   if (t > 1 & insulin_opt(t) > 0)
	      min_del = min(insulin_opt(t), min_del_insulin) ;
	      insulin_test(t-1) = insulin_opt(t-1) + min_del ;
	      insulin_test(t) = insulin_opt(t) - min_del ;
	      loss_test = compute_insulin_schedule_loss(insulin_test, parms, hyper_parms);
	      if (loss_test < loss_opt)
	      	 insulin_opt = insulin_test ;
		 loss_opt = loss_test ;
		 nchanged = nchanged + 1 ;
	      else
		insulin_test(t) = insulin_opt(t) ;
		insulin_test(t-1) = insulin_opt(t-1) ;
	      end
	   end
	end
    end
    pct_change = 100*(last_loss - loss_opt)/last_loss ;
    disp(sprintf('step %3.3d: loss %2.4f (%2.1f%%), changed %d', step, loss_opt, pct_change, nchanged)) ;
    if (nchanged == 0)
       break ;
    end
end
