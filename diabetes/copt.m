
clear('loss') ;
di = .02 ;

parms.nsamples = 128 ;
colors = ['r','b','g','c','m','y','k','ko'];
parms.rand = randn(2^15, 6) ;
for n=1:8
    delta_insulin = di:di:insulin_init(2) ;
    for i = 1:length(delta_insulin)
        del = delta_insulin(i) ;
	insulin_tmp = insulin_init;
	insulin_tmp(2) = insulin_init(2) - del ;
    	insulin_tmp(3) = insulin_init(3) + del ;
    	loss(n,i) = compute_insulin_schedule_loss(insulin_tmp, parms,hyper_parms);
     end
     parms.nsamples = parms.nsamples*2 ;
     hold on;
     plot(delta_insulin, loss(n,:),colors(n))
     hold off

end

