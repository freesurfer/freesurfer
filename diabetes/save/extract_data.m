read_medtronic_data

subject = 2 ;
sensor_glucose = cell2mat(S(subject).rows(:,sensor_glucose_ind)) ;
bg_glucose = cell2mat(S(subject).rows(:,bg_reading_ind)) ;

tps = size(sensor_glucose,1);

nparms = 2;     % # of timepoints in training


max_train = 0 ; 
for t=nparms+1:tps
    if (isnan(bg_glucose(t)) == 0)
       max_train = max_train + 1 ;
    end
end

% count the number of possible training events
Smat = zeros(max_train, nparms) ;
ntrain = 1 ;
for t=nparms+1:tps
    if (isnan(bg_glucose(t)) == 0)
       keyboard
       for p=1:nparms
	   Svec(p,1) = sensor_glucose(t-(p-1));
        end
       B(ntrain,1) = bg_glucose(t) ;
       Smat(ntrain,:) = Svec';
       ntrain = ntrain + 1 ;
   end
end
Smat_saved = Smat ;
[Smat_filled,Smat,replace] = fill_holes(Smat_saved) ;
B = B(find(replace == 0));

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
