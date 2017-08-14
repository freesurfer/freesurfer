function [best_k1, best_k2, best_k3, best_it] = fit_insulin_timecourse(iob_table) 


min_rms = 10000000 ;
best_k1 = 0 ;
best_k2 = 0 ;
best_k3 = 0 ;
for k1=.001:.001:1
  if (mod(k1, .001) == 0)
     disp(sprintf('k1 = %f', k1)) ;
     drawnow ;
  end
  for k2=.001:.001:1
    for k3=.001:.001:1
      it =  compute_insulin_timecourse(1, 3*60, k1, k2, k3);
      test_rms = rms(it'-iob_table(1:length(it))/100) ;
      if (test_rms < min_rms)
         best_it = it ;
         min_rms = test_rms ;
         best_k1 = k1 ;
         best_k2 = k2 ;
         best_k3 = k3 ;
         disp(sprintf('new best (%2.4f, %2.4f, %2.4f) found, min rms %2.4f', best_k1,best_k2,best_k3,test_rms)) ;
      end
    end
  end
 end

best_k1 = 0.021 ;
best_k2 = 0.001 ;
best_k3 = 0.021 ;
