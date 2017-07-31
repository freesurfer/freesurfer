
clear all ;
if (exist('S') == 0)
   S = read_medtronic_data('allCSV.mat')  ;
end
nsubjects = size(S,2) ;
define_medtronic_columns

for subject=1:nsubjects

clear('Smat', 'St', 'Bt', 'Svec','B') ;
sensor_glucose = cell2mat(S(subject).rows(:,sensor_glucose_ind)) ;
bg_glucose = cell2mat(S(subject).rows(:,bg_reading_ind)) ;

tps = size(sensor_glucose,1);

nparms = 15 ;     % # of timepoints in training


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

ndeltas = 5 ;
St = zeros(max_train-1, nparms+ndeltas) ;
for row=1:max_train
    Svec = Smat_filled(row,:) ;
    St(1:max_train-1,1:nparms) = Smat_filled([1:row-1 row+1:end],:) ;
    for n=nparms+1:nparms+ndeltas
    	Svec(n) = Svec(n-nparms)-Svec(n-nparms+1) ;
	for r2=1:max_train-1
	    St(r2,n) = St(r2,n-nparms)-St(r2,n-nparms+1) ;
	end
    end
    Bt = B([1:row-1 row+1:end]) ;
    p = pinv(St)*Bt ;
    predicted(row) = Svec*p ;
    err_predicted(row) = predicted(row) - B(row) ;
    err_cgm(row) = Svec(1) - B(row) ;
end

disp(sprintf('subject %d: mean CGM err = %2.0f, mean predicted err = %2.0f\n', subject,mean(abs(err_cgm)),mean(abs(err_predicted))));
end
