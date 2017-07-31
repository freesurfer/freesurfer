function [best_k, best_it] = fit_insulin_timecourse(iob_table) 


min_rms = 10000000 ;
best_k = 0 ;
for k=.0001:.0001:1
    it =  compute_insulin_timecourse(1, 3*60, k);
    test_rms = rms(it-iob_table(1:length(it))/100) ;
    if (test_rms < min_rms)
       best_it = it ;
       min_rms = test_rms ;
       best_k = k ;
%       disp(sprintf('new best %2.4f found, min rms %2.4f', k,test_rms)) ;
    end
 end
