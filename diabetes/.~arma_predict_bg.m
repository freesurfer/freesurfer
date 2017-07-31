
clear all ;
read_medtronic_data

for subject=1:nsubjects

clear('Smat', 'St', 'Bt', 'Svec','B') ;
sensor_glucose = cell2mat(S(subject).rows(:,sensor_glucose_ind)) ;
bg_glucose = cell2mat(S(subject).rows(:,bg_reading_ind)) ;

tps = size(sensor_glucose,1);

nparms = 10;     % # of timepoints in training


max_train = 0 ; 
for t=nparms+1:tps
    if (isnan(bg_glucose(t)) == 0)	
       p = 1 ;
       for dt=0:100   % search backwards for real CGM readings
       	   if (dt >= t)
	      break ;
	   end
       	   cgm = sensor_glucose(t-dt) ;
	   if (isnan(cgm) == 0)
	      p = p + 1 ;
	      if (p > nparms)
	      	 break ;
 	      end
           end
        end	
	if (p < nparms)
	   continue ;
	end
       max_train = max_train + 1 ;
    end
end

% count the number of possible training events
Smat = zeros(max_train, nparms) ;
ntrain = 1 ;
for t=nparms+1:tps
    if (isnan(bg_glucose(t)) == 0)
       p = 1 ;
       for dt=0:100   % search backwards for real CGM readings
       	   if (dt >= t)
	      break ;
	   end
       	   cgm = sensor_glucose(t-dt) ;
	   if (isnan(cgm) == 0)
	      Svec(p,1) = cgm ;
	      p = p + 1 ;
	      if (p > nparms)
	      	 break ;
 	      end
           end
        end	
	if (p < nparms)
	   continue ;
	end
       B(ntrain,1) = bg_glucose(t) ;
       Smat(ntrain,:) = Svec';
       ntrain = ntrain + 1 ;
   end
end
Smat_saved = Smat ;
%[Smat_filled,Smat,replace] = fill_holes(Smat_saved) ;
%B = B(find(replace == 0));
Smat_filled = Smat ;

max_train = size(Smat,1) ;
disp(sprintf('%d training blood-glucose samples found after data cleaning', max_train)) ;

for row=1:max_train
    Svec = Smat_filled(row,:) ;
    St = Smat_filled([1:row-1 row+1:end],:) ;
    Bt = B([1:row-1 row+1:end]) ;
    p = pinv(St)*Bt ;
    predicted(row) = Svec*p ;
    err_predicted(row) = predicted(row) - B(row) ;
    err_cgm(row) = Svec(1) - B(row) ;
end

disp(sprintf('subject %d: mean CGM err = %2.0f, mean predicted err = %2.0f\n', subject,mean(abs(err_cgm)),mean(abs(err_predicted))));
end
