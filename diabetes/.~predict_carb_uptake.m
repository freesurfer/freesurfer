function carbs_timecourse = predict_carb_uptake(carbs_t0, k_blood_to_stomach, k_stomach_to_blood, time,carb_time)
% carb_timecourse = predict_carb_uptake(carbs_t0, k_blood_to_stomach, k_stomach_to_blood, time,carb_time)

  

carbs_timecourse = zeros(size(time)) ;

time = time * 60 ;   % change it to minutes
t0 = min(time) ;
t1 = max(time) ;
dt = .05 ;
i = 1 ;
ic = 2 ;
stomach = carbs_t0 ;
blood = 0 ;
carbs_timecourse(ic) = 0 ;
for t=t0:dt:t1
   carbs_highres(i) = blood ;
   if (t >= 0)
      dblood =   dt *  k_blood_to_stomach * blood ;
      dstomach = dt *k_stomach_to_blood * stomach ;
      stomach = stomach + dblood - dstomach ;
      blood = blood + dstomach - dblood ;
   end
    
   if ( ( (t-dt) < time(ic)) & (t >= time(ic)))
     carbs_timecourse(ic) = blood ;
     ic = ic + 1 ;
     if (ic > length(time))
       break ;
   end
  end
  i = i+1;
end

if (nargin > 4)
  [min_time, min_ind] = min(abs(time-carb_time));
  [ztime, zind] = min(abs(time-0));
  min_ind = min_ind - zind ;
  carb_timecourse_new = zeros(size(carbs_timecourse)) ;
  carb_timecourse_new(min_ind+1:end) = carbs_timecourse(1:[end-min_ind]);
  carbs_timecourse = carb_timecourse_new ;
end
